import numpy as np
from scipy.ndimage.measurements import label as cclabel

def relabel(data_object, correlation_matrix):
    
    labels = np.argmax(correlation_matrix, axis=0) + 1
    unique_labels = np.unique(labels)
    
    mask_data = data_object["mask_data"]
    new_labels = np.copy(mask_data)
    new_labels[mask_data > 0.5] = labels
    
    for l in unique_labels:
        max_label_index = np.argmax(correlation_matrix[l-1,:])
        
        labeled_array, num_features = cclabel(new_labels==l)
        hold_cclabel = labeled_array[mask_data > 0.5][max_label_index] 
        labeled_array[labeled_array != hold_cclabel] = 0
        new_labels[(new_labels==l) & (labeled_array==0)] = 0
    
    data_object["label_data"] = new_labels[mask_data > 0.5]
    return(data_object)

def writelabels(mask_file, outfile, labels):
    maskvol = volumeFromFile(mask_file)
    outvol = volumeLikeFile(mask_file, outfile)
    outvol.data[maskvol.data > 0.5] = labels
    outvol.writeFile()
    outvol.closeVolume()
    maskvol.closeVolume()