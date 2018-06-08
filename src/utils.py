from functools import partial
import Bio.SeqUtils as su
import cobra
import cobra.test
import numpy as np
from sklearn.preprocessing import maxabs_scale, minmax_scale, normalize, scale, robust_scale
from sklearn.model_selection import train_test_split

# For the Karim model, adds necessary reactions for the pathway
def add_but(model):
    alc_dehydr = cobra.Reaction(id='ALCDBUT', name='Alcohol dehydrogenase (butanal)', subsystem='c')
    model.add_reaction(alc_dehydr)
    alc_dehydr.add_metabolites({'btcoa_c': -1, 'h_c': -1, 'nadh_c': -1, 'nad_c': 1, 'btal_c': 1})

    butanol = cobra.Metabolite(id='btol_c', name='1-butanol', compartment='c', charge=0, formula='C4H9OH')
    but_synth = cobra.Reaction(id='BUTSYN', name='Butanol synthesis', subsystem='c')
    model.add_reaction(but_synth)
    but_synth.add_metabolites({'btal_c': -1, 'h_c': -1, 'nadh_c': -1, 'nad_c': 1, butanol: 1})
    return model

# amt in g, vol in mL, mw in g/mol
def calc_conc(amt, vol, mw=None, seq=None, seq_type=None):
    # seq can be DNA or protein or an amino acid
    if seq:
        mw = su.molecular_weight(seq, seq_type)
    elif not mw:
        raise Exception('Need a molecular weight for non-DNA')
    conc = (amt * 1000) / (vol * mw)
    # returns Molar concentrations
    return conc

def conc_dilution(start_conc, vol_add, tot_vol):
    return start_conc * (vol_add / tot_vol)

# Scales flux data in a number of different ways
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

# Generates stacked flux datasets by resampling from an input flux distribution
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

# Resamples from flux samples in a biased manner to maximize correlation with the experimental dataset
def biased_resample(sorted_samp_data, n_rows=2000):
    np.random.seed(42)
    resamp_data = np.empty((n_rows, n_experiments, max_sz[0]))
    for j in range(n_rows):
        btol_val = 0
        for i in range(n_experiments):
            inds = (sorted_samp_data[:, i, btol_col] >= btol_val).nonzero()[0]
            if inds.any():
                ind = np.random.choice(inds[:5], size=1, replace=True)[0]
            else:
                ind = np.argmax(sorted_samp_data[:, i, btol_col])
            btol_val = sorted_samp_data[ind, i, btol_col]
            flxs = sorted_samp_data[ind, i, :]
            resamp_data[j][i][:] = flxs
        if j % (n_rows / 10) == 0:
            print j
    print 'done'
    return resamp_data

# Generates train and test datasets
def gen_train_test(data, y=None):
    if not y is None:
        return train_test_split(data, y, random_state=42)
    else:
        return train_test_split(data, random_state=42)

# Adds noise to a biased dataset
def add_noise(biased_resamp_data, enc, gen):
    noise_arr = np.logspace(start=-5, stop=1, num=10)
    corrs = []
    for noise in noise_arr:
        noisy_data = biased_resamp_data.copy()
        for i in range(n_experiments):
            s = noisy_data[:, i, :].shape
            noisy_data[:, i, :] += np.random.normal(scale=noise, size=s)
        corr = pred(noisy_data, enc, gen)
        corrs.append(corr)
    return zip(noise_arr, corrs) + [(0, pred(biased_resamp_data, enc, gen))]
