import Bio.PDB.Polypeptide as pdb
import numpy as np
import pandas as pd
import datasets
import utils

def get_conc(aa_f='../data/aa_mix.csv', nrg_f='../data/energy_mix.csv', cfps_start='../data/cfps_start.csv', seq_f='../genes/rfp.txt', cfps_final=None,
                rxn_amt=5, batch_size=50):
    n_batch = batch_size / rxn_amt
    if cfps_final:
        cfps_conc = pd.read_csv(cfps_final, index_col='compound')
    else:
        aa_mix = pd.read_csv('../data/aa_mix.csv', index_col='AA')
        nrg_mix = pd.read_csv('../data/energy_mix.csv', index_col='compound')
        with open('../genes/rfp.txt', 'r') as f:
            seq = f.read()
        cfps_conc = pd.read_csv('../data/cfps_start.csv', index_col='compound')
        dna_conc = utils.calc_conc(0.000750, 0.00496, seq=seq, seq_type='DNA')

        cfps_conc['amt'] = cfps_conc['amt'] / n_batch
        aa_mix['start_conc'] = aa_mix.apply(lambda row: utils.calc_conc(row['weight_add'], 1, 
                                                                  seq=pdb.three_to_one(row.name.upper()), seq_type='protein'), axis=1)
        aa_mix['conc_add'] = utils.conc_dilution(aa_mix['start_conc'], aa_mix['vol_add'], aa_mix['vol_add'].sum())
        pi_conc = utils.calc_conc(0.15, 5, mw=611.77)
        nrg_mix['start_conc'] = nrg_mix.apply(lambda row: utils.calc_conc(row['amt'], row['fill'], mw=row['mw']), axis=1)
        nrg_mix['conc_add'] = utils.conc_dilution(nrg_mix['start_conc'], nrg_mix['vol_add'], nrg_mix['vol_add'].sum())

        for cmpnd, vals in nrg_mix.iterrows():
            cfps_conc.loc[cmpnd] = [vals['conc_add'], (5.0 / n_batch)]
        for aa, vals in aa_mix.iterrows():
            aa = aa.lower()
            cfps_conc.loc[aa] = [vals['conc_add'], (10.0 / n_batch)]
        cfps_conc.loc['GENE'] = [dna_conc, (4.96 / n_batch)]
        cfps_conc['final_conc'] = utils.conc_dilution(cfps_conc['start_conc'], cfps_conc['amt'], rxn_amt)
    return cfps_conc

def get_exp_data(froot):
    # '../data/17_5_18_T7_mRFP_NLS.CSV'
    if froot == 'karim':
        df = pd.read_csv('../data/{0}_data.CSV'.format(froot))
        conds_norm = df
        conds_norm['OUT'] = conds_norm['AVG.1'] / conds_norm['AVG.1'].max()
        conds_norm.drop(columns=['AVG', 'STD', 'AVG.1', 'STD.1', 'Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)
    else:
        df = pd.read_csv('../data/{0}_data.CSV'.format(froot), skiprows=6)
        print df.shape
        gain_diff = df.shape[0] / 5
        times = df["Unnamed: 1"]
        df.drop('Unnamed: 1', inplace=True, axis=1)
        # Bad data
        if froot == 'manual':
            df.drop('E09', inplace=True, axis=1)
        # Remove negative control
        df.drop(df.columns[-2:], inplace=True, axis=1)
        
        gain2 = gain_diff
        outs = df[gain2:gain2+gain_diff - 1].mean(axis=0)
        # '../data/17_5_18_exp_conditions.csv'
        conds = pd.read_csv('../data/{0}_exp_conditions.csv'.format(froot))
        if froot == 'manual':
            conds.drop(conds.shape[0] - 1, inplace=True)
            conds_full = conds.reindex(np.repeat(conds.index.values, 2)).reset_index()
            conds_full = conds_full.drop(32).reset_index(drop=True)
        else:
            conds_full = conds.reindex(np.repeat(conds.index.values, 2)).reset_index()
        conds_full['OUT'] = outs.reset_index(drop=True)
        conds_avg = conds_full.groupby('index').mean()
        conds_norm = conds_avg
        conds_norm['OUT'] = conds_norm['OUT'] / conds_norm['OUT'].max()
        if 'level_0' in conds_norm.columns:
            conds_norm.drop('level_0', axis=1, inplace=True)
    conds_norm.to_csv('../data/{0}_EXPERIMENT.csv'.format(froot))
    return conds_norm

def get_test_data(froot, txtl, resamp, latent_dim, layer_szs, use_corr=True, n_epochs=200, batch_size=256, scale='flux_zero'):
    from keras.models import load_model
    flat = not resamp
    if not resamp:
        use_corr = False
    if froot == 'karim':
        assert(not txtl)
    fname = 'epochs={0}_batch={1}_dimension={2}_corr={3}_scale={4}_froot={5}_txtl={6}_nlayers={7}_resamp={8}{9}.h5'.format(
        n_epochs, batch_size, latent_dim, use_corr, scale, froot, txtl, len(layer_szs), resamp, '_lastlayer={0}'.format(layer_szs[-1]) if layer_szs[-1] != 256 else '')
    print('Load models {0}'.format('../models/encoder_' + fname))
    encoder = load_model('../models/encoder_{0}'.format(fname))
    generator = load_model('../models/generator_{0}'.format(fname))

    data_f = '../data/{0}{1}_{2}_fluxes'.format(froot, '_txtl' if txtl else '', 'stacked' if resamp else 'flat')
    X_train, y_train, X_test, y_test, obj_col, cols, y_vals_d = datasets.get_dataset(data_f)
    print 'Read in data from {0}'.format(data_f)
    y_vals = np.array(y_vals_d)
    X_test = np.array(X_test)
    y_test = np.array(y_test)

    test_enc = encoder.predict(X_test)
    print 'Encoded data'
    test_dec = generator.predict(test_enc)
    print 'Decoded data'

    return encoder, generator, X_test, y_test, obj_col, cols, y_vals_d, test_enc, test_dec

def read_flux_data(dir_name, flux_root, n_experiments, obj_name, resamp=True, scale=None, n_rows=50000):
    flux_dfs = []
    obj_col = None
    max_sz = (None, None)
    for i in range(n_experiments):
        d = pd.read_csv('{0}/{1}_{2}'.format(dir_name, flux_root, i), index_col=0)
        new_obj_col = d.columns.get_loc(obj_name)
        if obj_col:
            assert(new_obj_col == obj_col)
        else:
            obj_col = new_obj_col
        if d.shape[1] > max_sz[0]:
            max_sz = (d.shape[1], d.columns)
        flux_dfs.append(d)
    cols = max_sz[1]
    fluxes = []
    for d in flux_dfs:
        vals = d.reindex(columns=max_sz[1]).fillna(0)
        assert(vals.columns.get_loc(obj_name) == obj_col)
        fluxes.append(vals.values)
    if resamp:
        samp_data = np.stack(fluxes, axis=1)
        samp_data = utils.resample(samp_data, n_rows, n_experiments, max_sz[0])
    else:
        samp_data = np.vstack(fluxes)
        y = np.array([[i] * fluxes[0].shape[0] for i in range(n_experiments)]).reshape(samp_data.shape[0])

    if scale:
        utils.scale_data(samp_data, scale, in_place=True)

    if not resamp:
        X_train, X_test, y_train, y_test = utils.gen_train_test(samp_data, y)
    else:
        X_train, X_test = utils.gen_train_test(samp_data)
        y_train, y_test = None, None
    print 'Min: {0}, Max: {1}'.format(np.min(X_train), np.max(X_train))
    #np.savez_compressed('../data/fluxes_resampled', train=X_train, test=X_test)
    return X_train, y_train, X_test, y_test, obj_col, cols
