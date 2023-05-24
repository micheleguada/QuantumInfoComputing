import tensornetwork as tn
import numpy as np

from tqdm import tqdm as tqdm
from time import time
from datetime import datetime
import csv
import os


class GenerativeMPS(tn.FiniteMPS):
    """
    Class that allows to use a Matrix Product State as a generative model. It is built to inherit 
    from the tensornetwork class: 'FiniteMPS', as this allows to use some of its function and attributes.
    In particular:
    - random          : it initializes the MPS with random numbers in canonical form with the desired physical
                        and bond dimensions.
    - position        : it shifts the orthogonality center position using QR decomposition
    - backend.svd     : it performs a SVD decomposition and can optionally truncate the results (see 
                         'rebuild_bond' for details).
    - center_position : it stores the orthogonality center of the MPS
    
    Main attributes:
    - merged_tensor  : is variable that stores the merged tensor during training
    - data           : contains the data to be used for learning
    - nbatch         : number of batches in which splitting the data for training
    - left_site      : is the left site to be merged with the one to right during training
    - max_bond_dim   : maximum bond dimension allowed.
    - max_error      : maximum allowed error for truncation
    - descent_steps  : number of steps to do during gradient descent.
    - lr             : learning rate factor
    - use_Adam       : if to use Adam Optimizer during gradient descent
    - debug_         : simple debug flag
    - log_file       : the name of the log file. Today's date will be added to the name automatically.
    
    Main methods:
    - import_data     : import the data into the class and calculate the batch size using the 
                        passed argument 'nbatch' (default is 1)
    - train           : it starts the training loop by calling the appropriate functions in sequence.
                        it also create a log file in .csv format where training statistics are saved.
    - training_Loss   : calculate the Negative Log Likelihood (NLL) on the training data
    - test_Loss       : calculate NLL on the passed argument 'test_data'.
    - generate_sample : generates a new sample from the learned distribution
    - save_MPS_model  : saves the tensors of the MPS into a .npy file.
    
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
        
        # initializations
        self.merged_tensor = None
        self.data          = None
        self.nbatch        = None
        self.left_site     = self.center_position    # left one of the two sites that will be merged during training
        
        # hyper-parameters
        self.max_bond_dim  = 600
        self.max_error     = 1e-3
        self.descent_steps = 20
        self.lr            = 1e-4    #learning rate
        self.use_Adam      = True
        
        # diagnostics
        self.train_loss    = []
        self.test_loss     = []
        self.epochs_time   = []
        self.right_time    = []
        self.left_time     = []
        self.log_file      = "MPS_Logs_"

        # debugging utils
        self.debug_ = False # debugging flag to use inside the class


    def import_data(self, dataset, nbatch=1):
        """ 
        Imports the data into the model and calculate the batchsize to be used in training. 
        - dataset   : the input data, in a 2D array_like shape
        - nbatch    : number of batches
        Note: dataset must contain binarized pixels
        """

        self.data   = dataset
        self.nbatch = nbatch
        self.batchsize = self.data.shape[0]//self.nbatch

        return


    def merge_bond(self):
        """
        Merge two adiacent nodes to enable the training of the MPS with adaptive bond dimension.
        """
        kk = self.left_site
        
        self.merged_tensor = self.backend.einsum("ijk,klm->ijlm", self.tensors[kk], self.tensors[kk+1])

        return


    def rebuild_bond(self, going_right):
        """
        Restore the MPS after training the merged tensor.
        - going_right [bool] : if the direction of sweeping is 'right'
        """
        assert self.merged_tensor is not None
        kk = self.left_site
        
        # performing Singular Value Decomposition on merged_tensor with truncation
        U, s, V, cutted = self.backend.svd(self.merged_tensor,
                                           pivot_axis = 2,      #defines where to split the edges ('legs') of the merged tensor
                                           max_singular_values  = self.max_bond_dim,
                                           max_truncation_error = self.max_error,
                                           relative = True,     #if True, use max_truncation_error relative to the highest singular value
                                           )      

        if self.debug_: print("Removed singular values: ", len(cutted), " | S_rank: ", len(s))        

        if going_right: # needed to easily mantain canonicalization of the MPS
            V  = self.backend.einsum("i,ijk -> ijk", s, V)
            V /= self.backend.norm(V)
        else:
            U  = self.backend.einsum("kji,i -> kji", U, s)
            U /= self.backend.norm(U)

        # returning the splitted tensors to the MPS
        self.tensors[kk]   = U
        self.tensors[kk+1] = V

        self.merged_tensor = None
        
        # it is possible to move the center manually (instead of calling self.position(kk+/-1)) because
        # the procedure above mantain the canonicalization of the MPS
        ### Note: if not going_right: left_site = center_position -1  
        self.left_site       = kk -(1-2*int(going_right)) 
        self.center_position = self.center_position -(1-2*int(going_right)) 
        
        return


    def init_cumulants(self):
        """
        Initialize the cache (here called 'cumulants') for tensor contractions. 
        """       
        assert self.merged_tensor is None
        assert self.data is not None
        
        if self.center_position != 0:
            print("Right-canonicalizing the MPS...")
            self.position(0) 
            self.left_site = self.center_position
            
        # cumulants initialization [list of matrices]
        self.cumulants = [np.ones((self.data.shape[0],1))] + [0]*(self.__len__() -2) + [np.ones((1,self.data.shape[0]))]
            
        # MPS is in right-canonical form, so we start contractions from the right 
        ### Note: tensors are contracted with data by passing pixels values as 'self.data[:,idx]'
        for idx in range(self.__len__() -1, self.center_position+1, -1):   #right cumulants            
            self.cumulants[idx-1] = self.backend.einsum("jik,ki-> ji", self.tensors[idx][:, self.data[:,idx], :],
                                                                       self.cumulants[idx])
            
        return


    def update_cumulants(self, going_right):
        """
        Update self.cumulants after the training of the two splitted tensors. 
        Note that, thanks to mixed canonical form, it is only needed to update the cumulants 
        for the current site.
        - going_right [bool] : if the direction of the last sweep is 'right'
        """
        assert self.merged_tensor is None
        kk = self.left_site
        
        if going_right:
            self.cumulants[kk  ] = self.backend.einsum("ij,jik-> ik", self.cumulants[kk-1],
                                                                      self.tensors[kk-1][:, self.data[:,kk-1], :])
        else:
            self.cumulants[kk+1] = self.backend.einsum("jik,ki-> ji", self.tensors[kk+2][:, self.data[:,kk+2], :],
                                                                      self.cumulants[kk+2])
        return


    def _gradient_cumulants(self):
        """ 
        This function computes and returns the gradient of the NLL with respect to the merged tensor
        using a batch of the imported data.
        """
        
        assert self.merged_tensor is not None

        # getting a random batch from the training dataset
        indexes = np.random.choice(self.data.shape[0], size=self.batchsize, replace=False)
        batch = self.data[indexes,:]

        kk = self.left_site        
        
        # phi = d(Psi)/d(Aĸ,ĸ+1)
        phi = self.backend.einsum("ij,ki->ijk", self.cumulants[kk][indexes,:], 
                                                self.cumulants[kk+1][:,indexes])
        
        # psi = full contraction of the MPS with the merged_tensor and data in the batch
        ### note: merged tensor is contracted with the corresponding pixels in the batch by passing their values as indexes.
        ### this can be done because it is assumed that input data are binarized pixels.
        psi = self.backend.einsum("ij,jik,ki->i", self.cumulants[kk][indexes,:], 
                                                  self.merged_tensor[:,batch[:,kk],batch[:,kk+1],:], 
                                                  self.cumulants[kk+1][:,indexes])
        
        if (psi == 0).any():
            print("Gradient Descent: psi has some zero values. Skipping gradient computation at bond:", kk)
            return 0.
        
        psi_inv = 1./ psi
            
        # Computing gradient
        gradient = np.zeros( (self.bond_dimensions[kk],
                              self.physical_dimensions[kk],
                              self.physical_dimensions[kk+1],
                              self.bond_dimensions[kk+2]) )
        if self.debug_:
            print("gradient shape: ", gradient.shape, 
                  "\npsi_inv shape : ", psi_inv.shape, 
                  "\nphi shape     : ", phi.shape)

        for ii,jj in [(ii,jj) for ii in range(self.physical_dimensions[kk]) for jj in range(self.physical_dimensions[kk+1])]:
             mask = (batch[:,kk]==ii) * (batch[:,kk+1]==jj)                                                                 
             gradient[:,ii,jj,:] = self.backend.einsum("ijk,i->jk", phi[mask,:,:], 
                                                                    psi_inv[mask] )
        
        gradient = 2.*(self.merged_tensor - gradient/self.batchsize)
        
        return gradient


    def _bondtrain(self, going_right):
        """Training on the current bond
        - going_right: if the direction of the sweeping is 'right'
        """
        assert self.merged_tensor is None

        # merging two adiacent tensors
        self.merge_bond()
        
        # moments for Adam optimization procedure (used only if self.use_Adam = True)
        mt = 0.
        vt = 0.
        
        for tt in range(1, self.descent_steps+1):
            if self.debug_: print("descent_step: ", tt)
            
            # computing batch gradient
            grad = self._gradient_cumulants()
            
            # if True, the update is recomputed using Adam procedure
            if self.use_Adam:
                grad, mt, vt = self._AdamOptimizer(grad, mt, vt, tt)
                
            # updating and normalizing the merged_tensor
            self.merged_tensor -= self.lr * grad
            self.merged_tensor /= self.backend.norm(self.merged_tensor)
            
        # restoring the two tensors of the MPS from the merged tensor
        self.rebuild_bond(going_right)
        
        # updating the cumulants cache using the new computed tensors
        self.update_cumulants(going_right)

        return
    
    
    def _AdamOptimizer(self, grad, mt, vt, step):
        """ Adam optimization procedure.
        Implementation based on the original paper available at:
            https://arxiv.org/pdf/1412.6980.pdf
        """
        eps = 1e-8
        beta1 = 0.9
        beta2 = 0.999
        
        # computing 1st and 2nd moments
        mm = beta1*mt + (1.- beta1) * grad
        vv = beta2*vt + (1.- beta2) * grad**2
        
        # bias correction
        m_hat = mm / (1.- beta1**step)
        v_hat = vv / (1.- beta2**step)
        
        # compute the update
        update = m_hat / (self.backend.sqrt(v_hat) + eps)
        
        return update, mm, vv


    def train(self, epochs, test_data=None):
        """
        Training over several epochs. Before entering the loop, it ensures that the MPS is in right-canonical
        form and also creates the log file if it does not exit yet.
        - epochs [int] : number of epochs to be done. An epoch corresponds to a complete sweep from left to 
                         right and back to left.
        - test_data    : data to be used to compute the test NLL during training.
        """

        if self.center_position != 0 or self.left_site != 0:
            print("Right-canonicalizing the MPS and redoing cumulants initialization...")
            self.position(0) 
            self.left_site = self.center_position
            self.init_cumulants()
            self.merged_tensor = None

        # creating or opening the log file 
        fullname = self.log_file + datetime.now().strftime("%d-%m-%Y") + ".csv"
        start_train_time = datetime.now().strftime("%H.%M")
        
        if os.path.exists(fullname):
            log = open(fullname, mode="a", newline='') #open log file
            log_writer = csv.writer(log, delimiter=',')
        else:
            log = open(fullname, mode="w", newline='') #create log file
            log_writer = csv.writer(log, delimiter=',')
            log_writer.writerow( ["#Epoch", "EpochTime[s]", "RightTime[s]", "LeftTime[s]", "TrainingLoss",   # writing header
                                  "TestLoss", "HighestBondDim", "MeanBondDim", "MaxBondDimAllowed", 
                                  "LearningRate", "MaxErrorAllowed", "DescentSteps", "StartTrainTime", 
                                  "DataSize", "BatchSize"] )
            log.flush()
            
        # creating the progress bar for the training loop
        with tqdm(total=epochs, mininterval=1.) as pbar:
            pbar.set_description("epochs")
            
            for ep in range(epochs):  # training loop
                start_ep = time()
                
                with tqdm(total=2*self.__len__()-4, mininterval=1.) as bar_bond:
                    start_right = time()
                    
                    # sweeping from left to right
                    for bond in range(0, self.__len__()-2):
                        bar_bond.set_description("bond_right")
                        self._bondtrain(going_right=True)
                        bar_bond.update(1)
                    
                    start_left = time()
                    
                    # adjusting the center position due to the inverted direction of sweeping
                    self.center_position += 1

                    # sweeping from right to left
                    for bond in range(self.__len__() -2, 0, -1):
                        bar_bond.set_description("bond_left")
                        self._bondtrain(going_right=False)
                        bar_bond.update(1)
                    
                pbar.update(1)
                
                end_ep = time()   
                
                # readjusting center position
                self.position(0)
                
                # store diagnostics
                self.epochs_time.append(end_ep - start_ep)
                self.right_time.append(start_left - start_right)
                self.left_time.append(end_ep - start_left)
                self.train_loss.append(self.training_Loss())
                self.test_loss.append(self.test_Loss(test_data))
                
                # writing in log file
                log_writer.writerow( [ep, (end_ep-start_ep), (start_left-start_right), (end_ep-start_left), 
                                      self.train_loss[-1], self.test_loss[-1], np.max(self.bond_dimensions), 
                                      np.mean(self.bond_dimensions), self.max_bond_dim, self.lr, 
                                      self.max_error, self.descent_steps, start_train_time, 
                                      self.data.shape[0], self.batchsize] )
                log.flush()
                            
        log.close() #close log file                

        return


    def training_Loss(self):
        """Compute the average Negative Log-Likelihood (NLL) on the training set"""
        
        Psi_cumulants_2 = (self.give_Psi_cumulants())**2
        if (Psi_cumulants_2 == 0.).any():
            print("training_Loss: Psi_cumulants has zero entries. Skipping Training Loss computation and returning: -1.")
            trainLoss = -1  
        else:
            trainLoss = -np.log( Psi_cumulants_2 ).mean()   

        return trainLoss


    def test_Loss(self, test_data):
        """Compute the average Negative Log-Likelihood (NLL) on 'test_data'"""
        
        if test_data is None:
            return -10
        
        Psi_test_2 = (self.give_Psi(test_data))**2
        if (Psi_test_2 == 0.).any():
            print("test_Loss: Psi_test has zero entries. Skipping Test Loss computation and returning: -1.")
            testLoss = -1  
        else:
            testLoss = -np.log( Psi_test_2 ).mean()   

        return testLoss


    def give_Psi(self, states):
        """
        Compute the full contraction of the MPS with the data contained in 'states' variable.
        - states  : data in the same format as training data.
        """
        assert self.merged_tensor is None

        # check mps is right-canonical, otherwise shift center to 0
        if self.center_position != 0:
            self.position(0)
        
        psi = self.tensors[-1][:, states[:,0], 0]
        for idx in range(self.__len__()-1, -1, -1):
            psi = self.backend.einsum("kil,li->ki", self.tensors[idx][:, states[:,idx], :], psi)
            
        return psi


    def give_Psi_cumulants(self):
        """Compute the full contraction of the MPS with the training data stored in 'self.data'."""
        assert self.merged_tensor is None
        assert self.data is not None
        
        kk = self.left_site
        
        psi_cumuls = self.backend.einsum("ij,jik,kil,li-> i", self.cumulants[kk],
                                                              self.tensors[kk  ][:, self.data[:,kk  ], :],
                                                              self.tensors[kk+1][:, self.data[:,kk+1], :],
                                                              self.cumulants[kk+1])
        return psi_cumuls


    def generate_sample(self):
        """ It generate a new sample from the learned probability distribution. """

        if self.center_position != 0:
            print("Right-canonicalizing the MPS...")
            self.position(0)

        pixels = np.empty((self.__len__(),), dtype=np.int8)

        prev_psi = np.ones((1,))
        for idx in range(self.__len__()):
            if self.debug_: print(prev_psi.shape, " ^ ", self.tensors[idx][:,1,:].shape)
            psi = np.dot(prev_psi, self.tensors[idx][:,1,:])

            # conditional probability P(K|K-1) = P(K)/P(K-1)
            prob = ( np.linalg.norm(psi)/np.linalg.norm(prev_psi) )**2 
            if prob > np.random.rand():
                pixels[idx] = 1
                prev_psi = psi
            else:
                pixels[idx] = 0
                prev_psi = np.dot(prev_psi, self.tensors[idx][:,0,:])

        return pixels

        
    def save_MPS_model(self, filename, add_datetime=False):
        """
        Save the trained tensors of the MPS model in '.npy' file.
        - filename [str]      : prefix of the file name
        - add_datetime [bool] : if True, date and time will be attached to 'filename'. Default is False.
        
        The file will be named: [filename](_datetime).npy
        Example: MPSmodel_20-4-2021_23.31.npy
        """
        
        if add_datetime:
            filename = filename + datetime.now().strftime("_%d-%m-%Y_%H.%M")
        
        np.save(filename, self.tensors, allow_pickle=True)
        
        return
    
        
