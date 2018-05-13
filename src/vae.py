import keras
from keras.models import Sequential, Model
from keras.layers import Input, Dense, Activation, Lambda
from keras.callbacks import EarlyStopping
from keras import backend as K
from keras import metrics
import numpy as np
from matplotlib import cm
from collections import Counter
import tensorflow as tf
import scipy.stats
import pandas as pd
from tensorflow.contrib.metrics import streaming_pearson_correlation
from sklearn.preprocessing import maxabs_scale, minmax_scale, normalize, scale
np.random.seed(42)

parser = argparse.ArgumentParser(description='VAE on metabolic fluxes')
parser.add_argument('latent_dimension', metavar='d', type=int, help='latent dimensionality for vae',
                    default=2)
parser.add_argument('batch_size', metavar='b', type=int, help='batch_size',
                    default=256)
parser.add_argument('n_epochs', metavar='n', type=int, help='Number of epochs',
                    default=100)

def gen_train_test(data):
    train_ind = np.random.choice(data.shape[0], size=int(0.9 * data.shape[0]), replace=False)
    test_ind = list(set(range(data.shape[0])) - set(train_ind))
    min_val = np.min(data)
    max_val = np.max(data)
    # mn = np.mean(flat_data)
    #std = np.std(flat_data)
    #norm = lambda x: (x - mn) / std
    #scale = lambda x: (x + abs(min_val)) / (abs(min_val) + max_val)
    #data_scaled = scale(data)
    y = np.array(range(41) * data.shape[0])
    return data[train_ind], y[train_ind], data[test_ind], y[test_ind]

#'../data/flux_samps_2k/fluxes_{0}'
def read_data(dir_name, resamp=True, n_rows=50000):
    flux_dfs = []
    btol_col = None
    for i in range(n_experiments):
        d = pd.read_csv('../data/flux_samps_2k/fluxes_{0}'.format(i), index_col=0)
        new_btol_col = d.columns.get_loc('DM_btol_c')
        if btol_col:
            assert(new_btol_col == btol_col)
        else:
            btol_col = new_btol_col
        flux_dfs.append(d)
    max_sz = (None, None)
    for d in flux_dfs:
        if d.shape[1] > max_sz[0]:
            max_sz = (d.shape[1], d.columns)
    fluxes = []
    for d in flux_dfs:
        vals = d.reindex(columns=max_sz[1]).fillna(0)
        assert(vals.columns.get_loc('DM_btol_c') == btol_col)
        # TODO try different scaling technique?
        scaled_vals = minmax_scale(vals.values)
        fluxes.append(scaled_vals)
    samp_data = np.stack(fluxes, axis=1)
    if resamp:
        resamp_data = np.empty((n_rows, n_experiments, max_sz[0]))
        for j in range(n_rows):
            #exps = []
            inds = np.random.choice(range(samp_data.shape[0]), size=n_experiments, replace=True)
            for i in range(n_experiments):
                ind = inds[i]
                flxs = samp_data[ind, i, :]
                resamp_data[j][i][:] = flxs
                #exps.append(flxs)
            #resamp_data_l.append(exps)
            if j % (n_rows / 10) == 0:
                print j
        samp_data = resamp_data
    X_train, y_train, X_test, y_test = gen_train_test(samp_data)
    #np.savez_compressed('../data/fluxes_resampled', train=X_train, test=X_test)
    return X_train, y_train, X_test, y_test

def sampling(args):
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape=(K.shape(z_mean)[0], n_experiments, latent_dim), mean=0.,
                              stddev=epsilon_std)
    return z_mean + K.exp(z_log_var / 2) * epsilon

def pearson_corr(y_pred, y_true):
    cent_pred = y_pred - K.mean(y_pred)
    cent_tr = y_true - K.mean(y_true)

    std_pred = K.std(y_pred)
    std_tr = K.std(y_true)

    return K.mean(cent_pred*cent_tr)/(std_pred*std_tr)

def build_vae(X_shape, n_experiments, targets, output_ind, latent_dim=2, batch_size=100):
    encoded_dim1 = 1024
    encoded_sz = 256
    # Encoder network
    x = Input(shape=(n_experiments, X_shape,))
    #h = Dense(encoded_dim1, activation='relu')(x)
    h = Dense(encoded_sz, activation='relu')(x)#h)
    z_mean = Dense(latent_dim)(h)
    z_log_var = Dense(latent_dim)(h)
    
    # Sample points from latent space
    z = Lambda(sampling, output_shape=(n_experiments,latent_dim,))([z_mean, z_log_var])
    
    # Decoder network
    decoder_h = Dense(encoded_sz, activation='relu')
    #decoder_h2 = Dense(encoded_dim1, activation='relu')
    decoder_mean = Dense(X_shape, activation='sigmoid')
    h_decoded = decoder_h(z)
    #h_decoded2 = decoder_h2(h_decoded)
    x_decoded_mean = decoder_mean(h_decoded)#2)

    # end-to-end autoencoder
    vae = Model(x, x_decoded_mean)
    #output_flux = x_decoded_mean[:, :, output_ind]
    #experiment_loss = scipy.stats.spearmanr(targets, output_flux)
    #experiment_loss = pearson_corr(output_flux, targets)#streaming_pearson_correlation(output_flux, targets)
    xent_loss = X_shape * metrics.mean_squared_error(x, x_decoded_mean)
    kl_loss = - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)
    vae_loss = K.mean(xent_loss + kl_loss)# + experiment_loss
    #vae_loss = K.sum(vae_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer='rmsprop')
    vae.summary()

    # encoder, from inputs to latent space
    encoder = Model(x, z_mean)

    # generator, from latent space to reconstructed inputs
    decoder_input = Input(shape=(n_experiments, latent_dim,))
    _h_decoded = decoder_h(decoder_input)
    #_h_decoded2 = decoder_h2(_h_decoded)
    _x_decoded_mean = decoder_mean(_h_decoded)#2)
    generator = Model(decoder_input, _x_decoded_mean)
    return vae, encoder, generator


if __name__ == '__main__':
    args = parser.parse_args()

    df = pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)
    n_experiments = df.shape[0]

    dat = np.load('../data/fluxes_resampled.npz')
    X_train, X_test = dat['train'], dat['test']
    latent_dim = args.latent_dimension
    epsilon_std = 1.0

    batch_size = args.batch_size
    X_shape, n_experiments = X_train.shape[2], df.shape[0]
    targets = tf.convert_to_tensor(df['AVG.1'].values, dtype=tf.float32)
    output_ind = btol_col
    vae, encoder, generator = build_vae(X_shape, n_experiments, targets, output_ind, batch_size)
    es = EarlyStopping(patience=10)
    vae.fit(X_train, shuffle=True, epochs=n_epochs, batch_size=batch_size,
            validation_data=(X_test, None), callbacks=[es])
    vae.save('vae_epochs={0}_batch={1}_dimension={2}.h5'.format(n_epochs, batch_size, latent_dim))
    encoder.save('encoder_epochs={0}_batch={1}_dimension={2}.h5'.format(n_epochs, batch_size, latent_dim))
    generator.save('generator_epochs={0}_batch={1}_dimension={2}.h5'.format(n_epochs, batch_size, latent_dim))

    x_test_encoded = encoder.predict(X_test, batch_size=batch_size)
    x_test_gen = generator.predict(x_test_encoded, batch_size=batch_size)

    corrs = []
    for i in range(x_test_gen.shape[0]):
        corr = scipy.stats.pearsonr(x_test_gen[i, :, btol_col], df['AVG.1'])
        corrs.append(corr)
    print 'Mean correlation: ', np.mean(corrs, axis=0)
