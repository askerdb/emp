import pickle as pkl
import nimfa
import time
import numpy as np
import umap
import umap.plot
import pickle as pkl
def row_filter_intensity(X, bin_names, threshold = 1/1000):
    colsums = np.array(X.sum(axis = 0)).flatten()
    for i in range(X.shape[1]):
        X[:, i] = X[:, i]/colsums[i]
    rowsums = np.array(X.sum(axis = 1)).flatten()
    rowkeep = rowsums > threshold
    X = X[rowkeep, :]
    bin_names = [x for (x, v) in zip(bin_names, rowkeep) if v]
    return((X, bin_names))

def row_filter_intensity_csc(X, bin_names, threshold = 1/1000):
    X = X.tocsc()
    colsums = np.array(X.sum(axis = 0)).flatten()
    for i in range(X.shape[1]):
        X[:, i] = X[:, i]/colsums[i]
    rowsums = np.array(X.sum(axis = 1)).flatten()
    rowkeep = rowsums > threshold
    X = X[rowkeep, :]
    bin_names = [x for (x, v) in zip(bin_names, rowkeep) if v]
    return((X, bin_names))


X, bins, col_names = pkl.load(open("1907_EMPv2_INN_GNPS.mgf.pkl", 'rb'))
start = time.time()
X_f, bins_f = row_filter_intensity_csc(X, bins, threshold = 1/1)
print('It took {0:0.1f} seconds to do csc filter'.format(time.time() - start))
print(X_f.shape)
#start = time.time()
#nmf_model = nimfa.Nmf(X_f, rank=2)
#nmf_fit = nmf_model()
print('It took {0:0.1f} seconds to do NMF processes'.format(time.time() - start))

start = time.time()
mapper = umap.UMAP(metric='cosine', random_state=42, low_memory=True).fit(X_f)
umap.plot.points(mapper)
print('It took {0:0.1f} seconds to do NMF processes'.format(time.time() - start))
pkl.dump(mapper, open("1907_EMPv2_INN_GNPS.mgf.pkl", "wb"))
