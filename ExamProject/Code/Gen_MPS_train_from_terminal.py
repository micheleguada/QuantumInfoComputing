import numpy as np
import tensornetwork as tn
import gzip
import matplotlib.pyplot as plt
import pandas as pd
from skimage.transform import rescale

from GenerativeMPS import GenerativeMPS   #our model

# Import data MNIST
f = open('Data/train-images.idx3-ubyte','rb')

image_size = 28
num_images = 2000

# reading
f.read(16)
buf = f.read(image_size * image_size * num_images)
data = np.frombuffer(buf, dtype=np.uint8).astype(np.float32)

# from array to image 28x28
data = data.reshape(num_images, image_size, image_size)

# rescaling to 20 x 20
rescaled_size = 20
data_rescaled = np.zeros((num_images, rescaled_size, rescaled_size))
for ii in range(num_images):
    data_rescaled[ii,:] = rescale(data[ii,:,:], rescaled_size/image_size)

# back to array
data_rescaled = data_rescaled.reshape(num_images, rescaled_size**2)

# binarization of images
data_rescaled = np.floor(data_rescaled/255. + 0.5).astype(np.int8)

data = data_rescaled
data.shape


Npix = data.shape[1]   # number of pixels of the image
phys_dim = 2
bond_dim = 2

phys_list = [phys_dim]*Npix
bond_list = [bond_dim]*(Npix-1)

#print("phys: ", phys_list, "length: ", len(phys_list))
#print("bond: ", bond_list, "length: ", len(bond_list))

mps = GenerativeMPS.random(phys_list, bond_list, dtype=np.float64, canonicalize=True) #create MPS network

mps.import_data(data[0:1000], nbatch=10)
mps.init_cumulants()
mps.position(0)

mps.descent_steps = 20
mps.lr = 0.0001
mps.max_bond_dim = 400
mps.max_error    = 0.001
#mps.use_Adam     = False

#mps.debug_ = True

mps.train(epochs=5, test_data=data[1001:1201])

filename = "MPSmodel_epoch20_400"
np.save(filename, mps.tensors, allow_pickle=True)
