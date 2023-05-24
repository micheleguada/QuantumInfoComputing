import tensornetwork as tn
import numpy as np

from tqdm import tqdm as tqdm
from tqdm import trange
from time import time
from datetime import datetime
import csv
import os


class KMeansMPS():
    """ 
    Class to be used to perform a clusterization of a dataset with MPS based KMeans approach.
    The difference with the classical KMeans algorithm is that the centroids are represented by a Matrix
    Product State (MPS) with a fixed bond dimension and the data, through a feature map, are represented by
    a separable MPS.
    
    Parameters:
    - mps_size [int]     : the number of tensors that compose the centroids; should be equal to the number 
                           of features in the data;
    - phys_dim [int]     : the dimension of the physical indexes;
    - bond_dim [int]     : the dimension of the auxiliary indexes;
    - canonicalize [bool]: if initialize centroids in canonical form.
    
    Attributes:
    - centroids [list]   : the list of the centroids; centroids are objects of the class 'Centroid';
    - max_iter  [int]    : maximum number of KMeans iterations;
    - kmeans_pp [bool]   : if to use the initialization procedure based on the classical kmeans_++
                           (see 'init_centroids' for details)
    - num_pts [int]      : the number of points to be used in the initialization procedure (used only
                           if kmeans_pp = True);
    - training_loss [list]: stores the total loss at every KMeans iteration;
    - debug_ [bool]      : simple debug flag;
    - verbose_ [bool]    : simple verbosity flag;
    
    Methods:
    - import_data         : it imports the data contained in 'dataset' and also apply the selected feature map.
                            There are two possible feature_map: 'trig' for trigonometric feature map,
                            and 'lin', linear feature map.
    - _trig_feature_map   : it applies the trigonometric feature map to the passed data
    - _linear_feature_map : it applies a linear feature map to the passed data
    - Total_Loss          : it computes the total loss by summing the losses of the singular clusters
    - init_centroids      : it initializes the centroids using a simplyfied kmeans_++
    - train               : it performs the kmeans iterations to minimize the cost function
    - assign_labels       : it assign labels to unseen data based on the distance from the centroids
    
    """
    def __init__(self,
                 n_classes,
                 mps_size,
                 phys_dim,
                 bond_dim,
                 canonicalize = True,
                 ):

        #initialization
        self.mps_size = mps_size
        self.bond_dim = bond_dim
        self.phys_dim = phys_dim
        self.n_classes = n_classes
        
        # creating the centroids as random initialized MPS
        phys_list = [ phys_dim for i in range(mps_size) ]
        bond_list = [ bond_dim for i in range(mps_size-1) ]
        self.centroids = [ Centroid.random(phys_list,bond_list,dtype=np.float64,canonicalize=canonicalize) for i in range(n_classes) ]
        
        # other parameters
        self.max_iter  = 10
        self.kmeans_pp = True
        self.num_pts   = 20
        
        self.training_loss = []
        
        # debugging
        self.debug_    = False
        self.verbose_  = True
        
    
    def import_data(self, dataset, feature_map="trig"):
        """ 
        This function takes care of importing the data into the model and apply the requested feature map.
        The function also computes and store the min and max values of training data columns to allow 
        also their use during prediction (see 'assign_labels' below).
        
        dataset    : the input data in a 2D array-like shape;
        feature_map: it is a string variable indicating which feature map has to be applied;
        """
        
        if feature_map == "trig":
            self.mins  = np.amin(dataset, axis=0)
            dataset = dataset - self.mins
            
            self.maxs  = np.amax(dataset, axis=0)
            dataset = dataset/self.maxs

            # parametrize input data into separable MPSs with trigonometric feature map
            self.data = self._trig_feature_map(dataset)
            
        elif feature_map == "lin":
            self.mins  = np.amin(dataset, axis=0)
            dataset = dataset - self.mins
            
            self.maxs  = np.amax(dataset, axis=0)
            dataset = dataset/self.maxs
            
            # linear feature map
            self.data = self._linear_feature_map(dataset)
            
        else:
            print("Invalid feature map. Data not imported!")
            return
        
        # storing the name of the applied feature map
        self.applied_feature_map = feature_map
        
        # initializing labels
        self.labels = np.zeros(self.data.shape[0]).astype(np.int8)
        
        return
    
    
    def _trig_feature_map(self, dataset):
        """trigonometric feature map for the input data """
        
        temp = []
        for row in dataset:
            
            coss = np.cos(0.5*np.pi*row)
            sins = np.sin(0.5*np.pi*row)
            temp.append( np.stack((coss, sins), axis=-1).reshape(len(row),2) )
        
        return np.array(temp)
    
    
    def _linear_feature_map(self, dataset):
        """linear feature map """
        
        temp = []
        for row in dataset:
            up     = row.astype(np.float64)
            bottom = 1. - up
            temp.append( np.stack((up, bottom), axis=-1).reshape(len(row),2) )
        
        return np.array(temp)

    
    def Total_Loss(self):
        """ It computes the total loss by summing up the losses of all the centroids. """
        
        Loss = 0.
        for cc in range(self.n_classes):
            Loss += self.centroids[cc].final_losses[-1].astype(np.float64)
                
        return Loss[0]
    
    
    def init_centroids(self):
        """ 
        Initialization of the centroids by with a procedure similar to the classic 'kmeans_++' procedure:
        - the first centroid is initialized with a data point selected randomly;
        - then, any other centroid is initialized iteratively:
            * for a batch of points, all the minimum distances between the point and the already computed centroids
              computed and stored in 'dists';
            * the new centroid is initialized as the point with the greatest minimum distance.
            
        The only parameter of this procedure is:
        self.num_pts: the size of the batch of points to use
        
        This procedure is used only if self.kmeans_pp is True.
        
        """
        
        # initialization of first centroid with a random point
        self.centroids[0].sweep(self.data[np.random.randint(0,self.data.shape[0],1)])
        
        for jj in range(1, self.n_classes):
            batch = np.random.randint(0,self.data.shape[0],self.num_pts)
            dists = np.zeros(self.num_pts)
            for ii in range(self.num_pts):    
                # taking the minimum distance between points and every already initialized centroid
                dists[ii] = np.min( [self.centroids[kk].distance(self.data[batch[ii]]) for kk in range(jj)] )
                
            # selecting the point with the maximum distance from the closest centroid
            point_next_id = batch[ np.argmax(dists) ]
            self.centroids[jj].sweep(self.data[[point_next_id]])
        
        return
    
    
    def train(self):
        """ This function defines the KMeans training loop. 
        It has three main steps (the same as for a classical kmeans):
        - 1) assignment of data to centroids according to their distance
        - 2) recompute the centroids with the new points
        - 3) check exit condition, that in this case is whether the new labels are equal to the old ones or not.
        """
        
        labels_equals = False
        index = 0
        
        if self.kmeans_pp:
            self.init_centroids()
        
        # Loop 
        while not labels_equals and index < self.max_iter:
            old_labels = np.copy(self.labels)
            
            if self.verbose_:
                print("Iteration #: ", index)
                
            # Step 1): assignment of data to centroids
            if self.verbose_:
                print("1) Assigning data to centroids...")
                
            for ii in range(len(self.data)):
                distances = np.zeros(self.n_classes)
                for cc in range(self.n_classes):
                    distances[cc] = self.centroids[cc].distance(self.data[ii])  #shape data[ii] : (mps_size, 2)
            
                ## assigning label according to minimum distance
                self.labels[ii] = np.argmin(distances)
            
            if self.debug_:
                print("DEBUG -> labels: ", self.labels)
            
            # Step 2): calculate new centroids
            if self.verbose_:
                print("2) Calculating new centroids...")
            for cc in range(self.n_classes):
                
                if self.verbose_: print("  Calculating centroid #: ", cc)
                assigned_ids = (self.labels == cc)
                if np.sum(assigned_ids) > 0:
                    if self.verbose_: print("    # assigned data: ", np.sum(assigned_ids))
                    self.centroids[cc].sweep(self.data[assigned_ids,:])
                
            # exit condition
            labels_equals = np.array_equal(old_labels, self.labels)
            index += 1
            
            #saving loss
            self.training_loss.append(self.Total_Loss())
        
        return
    
    
    def assign_labels(self, test_data, apply_feature_map=True):
        """ 
        This function is used to assign labels to previously unseen data (aka test_data).
        It also (optionally) takes care of applying the same rescaling and feature map to the test_data as done
        for the training data.
        """
               
        if apply_feature_map:
            
            if self.applied_feature_map == "trig":
                # rescaling data the same as the training set
                test_data = (test_data -self.mins) / self.maxs
                
                featured = self._trig_feature_map(test_data)
            elif self.applied_feature_map == "lin":
                # rescaling data the same as the training set
                test_data = (test_data -self.mins) / self.maxs
                
                featured = self._linear_feature_map(test_data)
            else:
                print("Invalid feature map reference. Skipping test_data labeling...")
                return [-1]
         
        labels = []
        for row in featured:
            dists = np.zeros(self.n_classes)
            for cc in range(self.n_classes):
                dists[cc] = self.centroids[cc].distance(row)
            labels.append(np.argmin(dists))
        
        return labels
    

#-------------------------------------------------------------------------------------------------#


class Centroid(tn.FiniteMPS):
    """ 
    Class that contain the definition of a centroid as a MPS for quantum inspired K-Means clustering. 
    
    Methods:
    - _init_AB     : it initializes the cache of tensor contraction with data
    - _update_AB   : it updates the cache after a minimization step
    - _minimize    : it performs the optimization of a tensor of the Centroid MPS
    - distance     : it computes the distance from the centroids of the given data 'Xn'
    - sweep        : it sweep back and forth iteratively in order to optimize all the 
                     tensors of the chain
    - LossFunction : it computes the loss of the centroid by summing up the distances of
                     the data assigned to it.
    """    
    
    def __init__(self,
                tensors,
                center_position = None,
                canonicalize = True,
                backend = None):

        super().__init__(tensors=tensors,
                        center_position=center_position,
                        canonicalize=canonicalize,
                        backend=backend,
                        )
        
        # diagnostics
        self.initial_losses = []
        self.final_losses   = []
        self.Ndata          = []
        
        
    def _init_AB(self, data):
        """
        Initialization of the cache of tensor contraction that is useful during training sweeps.
        """

        if self.center_position != 0:
            print("Right-canonicalizing the MPS...")
            self.position(0) 
            
        ## cumulants initialization [list of matrices]
        self.AB = [np.ones((data.shape[0],1))] + [0]*(self.__len__() -1) + [np.ones((1,data.shape[0]))]
        
        ## MPS is in right-canonical form, so we start contractions from the right (#j = b_l+1, #i=sigma_l+1, #k = b_l+2)
        for idx in range(self.__len__()-1, self.center_position, -1): # first: N-1; last: 1
            self.AB[idx] = self.backend.einsum("kn,jik,ni -> jn", self.AB[idx+1], self.tensors[idx], data[:,idx,:])
            
        return


    def _update_AB(self, data, going_right):
        """ 
        Update self.AB after the optimization of a tensor of the MPS.
        - going_right [bool] : if the direction of the last sweep is 'right'.
        """
        kk = self.center_position
        
        if going_right: #j = b_l ,  #i=sigma_l-1,  #k = b_l-1
            self.AB[kk  ] = self.backend.einsum("nk,kij,ni -> nj", self.AB[kk-1], self.tensors[kk-1], data[:,kk-1,:])         
            
        else:           #j = b_l+1, #i=sigma_l+1,  #k = b_l+2
            self.AB[kk+1] = self.backend.einsum("kn,jik,ni -> jn", self.AB[kk+2], self.tensors[kk+1], data[:,kk+1,:])  

        return
    
    
    def distance(self, Xn):
        """
        Calculate distance between the centroid and data instance. The distance is calculated 
        as: 1 - the overlap resulting from the full contraction of the MPS with the datum. 
        
        Xn : row containing vectors of feature map for one data instance.
        """       
        
        PP = self.backend.einsum("i,jik -> jk"  , Xn[-1,:], self.tensors[-1])
        for idx in range(self.__len__()-2, -1, -1):
            # auxiliary index contraction
            BB = self.backend.einsum("lij,jk -> lik", self.tensors[idx], PP)
            
            # physical index contraction
            PP = self.backend.einsum("i,jik -> jk"  , Xn[idx,:], BB)

        overlap = PP.reshape(1)
        
        return 1. - overlap
    
    
    def LossFunction(self, data):
        """ It computes the Loss of the centroid as the sum of all the distances of the centroid from its points. """
        
        Loss = 0.
        for ii in range(data.shape[0]):
            Loss += self.distance(data[ii])
                
        return Loss
    
    
    def _minimize(self, data):
        """ 
        It performs the optimization of the tensor based on the assigned data. Note that the updated 
        tensor is the current orthogonality center of the MPS.
        """
        kk = self.center_position
        
        #j = b_l   #k = b_l+1    #i = sigma_l   #n=n        
        ### Note: data shape [#data row, feature col, mapped feature]
        self.tensors[kk]  = self.backend.einsum("nj,ni,kn -> jik", self.AB[kk], data[:,kk,:], self.AB[kk+1])
        
        return


    def sweep(self, data):
        """ 
        This function performs the sweeping back and forth in order to recompute the centroid.
        """
        
        assert self.center_position == 0
        
        self.Ndata.append(data.shape[0])
        self.initial_losses.append(self.LossFunction(data))
        
        # initializing cache for contractions
        self._init_AB(data)
        
        # going right
        for idx in range(self.__len__()-1): 
            # minimize tensor and shift position onward
            self._minimize(data)
            ZZ = self.position( idx+1 )

            # update the cache
            self._update_AB(data, going_right=True)
            
        # going left
        for idx in range(self.__len__()-1, 0, -1):
            # minimize tensor and shift position backward
            self._minimize(data)
            ZZ = self.position( idx-1 )
            
            # update the cache
            self._update_AB(data, going_right=False) 
            
        self.final_losses.append(self.LossFunction(data))
        
        return
