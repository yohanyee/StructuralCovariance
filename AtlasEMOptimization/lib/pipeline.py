from structuralcovariance import get_data, compute_correlations
from segmentation import relabel, writelabels
from metrics import overlap
import os

def EMOptimize(files, label_file, mask_file, tol=1e-6, max_passes=10, outdir="/tmp"):
    
    # Initialize
    data_object = get_data(files, label_file, mask_file)
    correlation_matrix = compute_correlations(data_object)
    
    # Save starting point
    outfile = os.path.join(outdir, "AEMO-labels-starting.mnc")
    writelabels(mask_file, outfile, data_object["label_data"])
    
    
    # Pass zero and save
    data_object = relabel(data_object, correlation_matrix)
    outfile = os.path.join(outdir, "AEMO-labels-pass_test.mnc")
    writelabels(mask_file, outfile, data_object["label_data"])
    
    # EM algorithm
    olap = 0
    olap_diff = 1
    p = 1
    
    while ((olap_diff > tol) & (p <= max_passes)):
        # Expectation 
        print("Pass {p}: expectation step".format(p=p))
        correlation_matrix = compute_correlations(data_object)
        
        # Maximization
        print("Pass {p}: maximization step".format(p=p))
        old_labels = np.copy(data_object["label_data"])
        data_object = relabel(data_object, correlation_matrix)
        
        # Save
        print("Pass {p}: saving labels".format(p=p))
        outfile = os.path.join(outdir, 
                               "AEMO-labels-pass_{p:04d}.mnc".format(p=p))
        writelabels(mask_file, outfile, data_object["label_data"])
        
        # Compute metric
        print("Pass {p}: computing overlap".format(p=p))
        olap_new = overlap(data_object["label_data"], old_labels)
        olap_diff = abs(olap_new - olap)
        olap = olap_new
        print("Overlap is {ol}".format(ol=olap))
        
        convergence_datafile = os.path.join(outdir, "convergence_data.txt")
        f = open(convergence_datafile, 'a')
        f.write('Pass {p:04d}: overlap is {ol}'.format(p=p, ol=olap))  # python will convert \n to os.linesep
        f.close()
        
        # Iterate
        p += 1
        
    print("Optimization complete after {p} passes.".format(p=p))
        
        
    