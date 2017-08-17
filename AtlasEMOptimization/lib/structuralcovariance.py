import numpy as np
import sys
from pyminc.volumes.factory import *
from skimage.morphology import watershed
from skimage.feature import peak_local_max


def get_data(files, label_file, mask_file):
    m = volumeFromFile(mask_file, labels=True)
    l = volumeFromFile(label_file, labels=True)
    num_files = len(files)
    num_voxels = m.data[m.data > 0.5].shape[0]
    
    # Get jacobian determinants data
    data_matrix = np.empty([num_files, num_voxels], dtype=np.float32)
    data_matrix_size = sys.getsizeof(data_matrix)/(1024**3)
    sys.stdout.write("Allocated a matrix of size {size:.2f} GB\n"
          .format(size=data_matrix_size))
    sys.stdout.write("Loading file: ")
    for i in range(num_files):
        sys.stdout.write("{i} ".format(i=i))
        sys.stdout.flush()
        f = files[i]
        v = volumeFromFile(f, dtype='float')
        data_matrix[i, :] = v.data[m.data > 0.5]
        v.closeVolume()
    sys.stdout.write("done.\n")
    sys.stdout.flush()
    
    # Return and close
    out = {"data_matrix"    : data_matrix,
           "label_data"     : l[m.data > 0.5],
           "mask_data"      : m.data}   
    l.closeVolume()
    m.closeVolume()
    
    return(out)
    

def compute_correlations(data_object):
    existing_labels = np.unique(data_object['label_data'])
    existing_labels = np.delete(existing_labels, np.where(existing_labels==0)[0])
    
    data_matrix = data_object["data_matrix"]
    label_data = data_object["label_data"]
    mask_data = data_object["mask_data"]
    
    num_regions = len(existing_labels)
    num_voxels = mask_data[mask_data > 0.5].shape[0]
    
    correlation_matrix = np.empty([num_regions, num_voxels], dtype=np.float32)
    correlation_matrix_size = sys.getsizeof(correlation_matrix)/(1024**3)
    sys.stdout.write("Allocated a matrix of size {size:.2f} GB\n"
          .format(size=correlation_matrix_size))
    
    n = float(data_matrix.shape[0])
    sys.stdout.write("Working on label: ")
    for i in range(num_regions):
        l = existing_labels[i]
        seed_values = data_matrix[:,label_data==l].mean(axis=1)
        s = sum(seed_values)
        s2 = sum(seed_values**2)
        sq = s**2
        sys.stdout.write("{i} ".format(i=i))
        sys.stdout.flush()
        for j in range(num_voxels):
            if (j % 10000)==0:
                sys.stdout.write(".")
                sys.stdout.flush()
            t = sum(data_matrix[:, j])
            t2 = sum(data_matrix[:, j]**2)
            tq = t**2
            st = sum(seed_values * data_matrix[:, j])
            row = ((n*st)-(s*t))/(np.sqrt((n*s2)-(sq))*np.sqrt((n*t2)-(tq)))
        correlation_matrix[i, :] = row
            
    sys.stdout.write("done.\n")
    sys.stdout.flush()
    return(correlation_matrix)