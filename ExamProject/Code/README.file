## MPS based unsupervised learning - Quantum Information and Computing project ##
year    : 2020/2021
Authors : Alessandro Lambertini,
          Michele Guadagnini


This folder contains the implementations of two methods for unsupervised learning by mean of Matrix Product States (MPS):
- 'GenerativeMPS' contained in the file: 'GenerativeMPS.py'
- 'K-MeansMPS' contained in the file: 'KMeansMPS.py'
Both the classes are built on top of the 'tensornetwork' library, described in:
    [https://arxiv.org/pdf/1905.01330.pdf]
To install the library it is only needed to type:
    pip install tensornetwork

We provide also two ready to use applications:
- A jupyter notebook for the k-means algorithm in the file: KMeansMPS_applications.ipynb
- A script to be run on terminal for the generative model in the file: Gen_MPS_train_from_terminal.py

######################################################################################################################################

GenerativeMPS:

    Class that allows to use a Matrix Product State as a generative model. It is built to inherit
    from the tensornetwork class: 'FiniteMPS', as this allows to use some of its function and attributes.
    In particular:
    - random          : it initializes the MPS with random numbers in canonical form with the desired physical
                        and bond dimensions.
    - position        : it shifts the orthogonality center position using QR decomposition
    - backend.svd     : it performs a SVD decomposition and can optionally truncate the results.
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

#----------------------------------------------------------------------------------------------------#

Below we report a simple usage example of the class:

* initialization of the MPS:
    Npix = data.shape[1]   # number of pixels of the image
    phys_dim = 2
    bond_dim = 2

    phys_list = [phys_dim]*Npix     # list of physical dimensions
    bond_list = [bond_dim]*(Npix-1) # list of bond dimensions

    mps = GenerativeMPS.random(phys_list, bond_list, dtype=np.float64, canonicalize=True)


* import data and train the model:
    datasize = 1000
    nbatches = 10
    mps.import_data(data[0:datasize,:], nbatch=nbatches)

    # initialization of the cumulants cache
    mps.init_cumulants()

    # training loop
    testsize = 200
    epochs = 20
    mps.train(epochs=epochs, test_data=data[datasize:datasize+testsize,:])


* generate a set of 9 samples from the trained model:
    Nimages = 9
    size    = 20   #number of pixels on each axis (must be the same used in training)
    images = [mps.generate_sample().reshape(size,size) for i in range(Nimages)]


* save and load a trained model:
    # saving
    mps.save_MPS_model("MPS_model", add_datetime=True)

    #loading
    tensors_list = numpy.load("model_name", allow_pickle=True)
    mps = GenerativeMPS(tensors_list, canonicalize=True)


######################################################################################################################################


K-MeansMPS:

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

#----------------------------------------------------------------------------------------------------#

Below we report a simple usage example of the class:

    # KMeansMPS model initialization
    kmeans_mps = KMeansMPS(n_classes = Ncs,
                           mps_size  = dataset.shape[1],
                           phys_dim  = 2,
                           bond_dim  = 1,
                           canonicalize = True)

    # parameters
    kmeans_mps.max_iter = 20
    kmeans_mps.num_pts  = 50
    kmeans_mps.kmeans_pp= True
    kmeans_mps.verbose_ = False

    #importing data
    kmeans_mps.import_data(dataset, feature_map="trig")

    #training and get assigned labels
    kmeans_mps.train()
    assigned_labels = kmeans_mps.labels

    #predict labels on unseen data
    kmeans_mps.assign_labels(test_data)
