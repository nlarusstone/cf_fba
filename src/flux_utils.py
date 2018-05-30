import numpy as np
import pandas as pd
from functools import partial
from sklearn.preprocessing import maxabs_scale, minmax_scale, normalize, scale, robust_scale
from sklearn.model_selection import train_test_split
np.random.seed(42)

def scale_data(data, scale_type, in_place=True):
    if not in_place:
        scaled_data = data.copy()
    else:
        scaled_data = data
    if 'norm' in scale_type:
        scale_func = scale
    elif 'robust' in scale_type:
        scale_func = robust_scale
    elif 'maxabs' in scale_type:
        scale_func = maxabs_scale
    elif 'negone' in scale_type:
        scale_func = partial(minmax_scale, feature_range=(-1, 1))
    else:
        scale_func = minmax_scale
    if len(data.shape) < 3:
        print len(data.shape)
        scaled_data = scale_func(scaled_data, copy=(not in_place))
    else:
        if 'flux' in scale_type:
            for i in range(data.shape[2]):
                    scaled_data[:, :, i] = scale_func(scaled_data[:, :, i])
        elif 'exp' in scale_type:
            for i in range(data.shape[1]):
                    scaled_data[:, i, :] = scale_func(scaled_data[:, i, :])
        else:
            raise('No scale direction')
    
    if not in_place:
        return scaled_data

def resample(data, n_rows, n_experiments, n_rxns):
    resamp_data = np.empty((n_rows, n_experiments, n_rxns))
    for j in range(n_rows):
        inds = np.random.choice(range(data.shape[0]), size=n_experiments, replace=True)
        for i in range(n_experiments):
            ind = inds[i]
            flxs = data[ind, i, :]
            resamp_data[j][i][:] = flxs
        if j % (n_rows / 10) == 0:
            print j
    print 'done'
    return resamp_data

def biased_resample(sorted_samp_data, n_rows=2000):
    #sorted_samp_data = samp_data_scaled.copy()
    #for i in range(n_experiments):
    #    exp = samp_data_scaled[:, i, :]
    #    sorted_samp_data[:, i, :] = samp_data_scaled[exp[:,btol_col].argsort(), i, :]
    np.random.seed(42)
    resamp_data = np.empty((n_rows, n_experiments, max_sz[0]))
    for j in range(n_rows):
        #exps = []
        btol_val = 0
        for i in range(n_experiments):
            inds = (sorted_samp_data[:, i, btol_col] >= btol_val).nonzero()[0]
            if inds.any():
                ind = np.random.choice(inds[:5], size=1, replace=True)[0]
            else:
                ind = np.argmax(sorted_samp_data[:, i, btol_col])
            btol_val = sorted_samp_data[ind, i, btol_col]
            #exps.append(btol_val)
            flxs = sorted_samp_data[ind, i, :]
            resamp_data[j][i][:] = flxs
        #print scipy.stats.pearsonr(exps, df['AVG.1'])
        if j % (n_rows / 10) == 0:
            print j
    print 'done'
    return resamp_data

def gen_train_test(data, y=None):
    #train_ind = np.random.choice(data.shape[0], size=int(0.9 * data.shape[0]), replace=False)
    #test_ind = list(set(range(data.shape[0])) - set(train_ind))
    #y = np.array(range(41) * data.shape[0])
    if not y is None:
        return train_test_split(data, y, random_state=42)
    else:
        return train_test_split(data, random_state=42)
    #return data[train_ind], y[train_ind], data[test_ind], y[test_ind]

def read_data(dir_name, resamp=True, scale=None, n_rows=50000):
    df = pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    n_experiments = df.shape[0]
    flux_dfs = []
    btol_col = None
    max_sz = (None, None)
    for i in range(n_experiments):
        d = pd.read_csv('{0}/fluxes_{1}'.format(dir_name, i), index_col=0)
        new_btol_col = d.columns.get_loc('DM_btol_c')
        if btol_col:
            assert(new_btol_col == btol_col)
        else:
            btol_col = new_btol_col
        if d.shape[1] > max_sz[0]:
            max_sz = (d.shape[1], d.columns)
        flux_dfs.append(d)
    cols = max_sz[1]
    fluxes = []
    for d in flux_dfs:
        vals = d.reindex(columns=max_sz[1]).fillna(0)
        assert(vals.columns.get_loc('DM_btol_c') == btol_col)
        fluxes.append(vals.values)
    if resamp:
        samp_data = np.stack(fluxes, axis=1)
        samp_data = resample(samp_data, n_rows, n_experiments, max_sz[0])
    else:
        samp_data = np.vstack(fluxes)
        y = np.array([[i] * fluxes[0].shape[0] for i in range(n_experiments)]).reshape(samp_data.shape[0])

    if scale:
        scale_data(samp_data, scale, in_place=True)

    if not resamp:
        X_train, X_test, y_train, y_test = gen_train_test(samp_data, y)
    else:
        X_train, X_test = gen_train_test(samp_data)
        y_train, y_test = None, None
    print 'Min: {0}, Max: {1}'.format(np.min(X_train), np.max(X_train))
    #np.savez_compressed('../data/fluxes_resampled', train=X_train, test=X_test)
    return X_train, y_train, X_test, y_test, btol_col, cols
