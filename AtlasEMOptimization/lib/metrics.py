import numpy as np

def overlap(new_labels, old_labels):
    ol = np.equal(new_labels, old_labels)
    return(float(len(np.where(ol)[0]))/len(ol))