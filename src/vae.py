import keras
from keras.models import Sequential, Model
from keras.layers import Input, Dense, Activation, Lambda
from keras.callbacks import EarlyStopping, Callback
from keras import backend as K
from keras import metrics
import numpy as np
from matplotlib import cm
from collections import Counter
import tensorflow as tf
import scipy.stats
import pandas as pd
from tensorflow.contrib.metrics import streaming_pearson_correlation
from functools import partial
from sklearn.preprocessing import maxabs_scale, minmax_scale, normalize, scale, robust_scale
import argparse
import pickle
import flux_utils as futils
import create_dataset as dataset
np.random.seed(42)

parser = argparse.ArgumentParser(description='VAE on metabolic fluxes')
parser.add_argument('-d', '--dim', metavar='d', type=int, help='latent dimensionality for vae',
                    default=2)
parser.add_argument('-b', '--batch', metavar='b', type=int, help='batch_size',
                    default=256)
parser.add_argument('-n', '--epochs', metavar='n', type=int, help='Number of epochs',
                    default=100)
parser.add_argument('--corr', dest='corr', help='Correlation loss', default=True, action='store_true')
parser.add_argument('--no-corr', dest='corr', help='Correlation loss', action='store_false')
parser.add_argument('-s', '--scale', metavar='s', type=str, default='flux_zero', help='type of data scaling')
parser.add_argument('-l','--layers', nargs='+', help='Layer sizes', type=int, default=[1024])


parser.add_argument('-f', '--froot', type=str, default='hand')
parser.add_argument('-r', '--resamp', dest='resamp', help='Resample', default=False, action='store_true')
parser.add_argument('--no-resamp', dest='resamp', help='Dont resample', action='store_false')
parser.add_argument('-t', '--txtl', dest='txtl', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', action='store_false')

class LossHistory(Callback):
    def __init__(self):
        self.recon_losses = []
        self.kl_losses = []
        self.corr_losses = []
        
    def on_epoch_end(self, epoch, logs={}):
        y_true = self.validation_data[0]
        y_pred = self.model.predict(self.validation_data[0])
        cor_loss = -1 * corr_loss(y_vals, y_pred[:, :, obj_col.value], be=np)
        xent_loss = y_true.shape[-1] * np.mean(np.square(y_true - y_pred), axis=-1)
        inputs = [K.learning_phase()] + self.model.inputs
        zvar = K.function(inputs=inputs, outputs=[self.model.get_layer('z_log_var').output])
        zmn = K.function(inputs=inputs, outputs=[self.model.get_layer('z_mean').output])
        z_log_var = zvar([0, self.validation_data[0]])[0]
        z_mean = zmn([0, self.validation_data[0]])[0]
        kl_loss = - 0.5 * np.sum(1 + z_log_var - np.square(z_mean) - np.exp(z_log_var), axis=-1)
        print "Reconstruction loss: {0}".format(np.mean(xent_loss))
        print "KL loss: {0}".format(np.mean(kl_loss))
        print "Corr loss: {0}".format(cor_loss)
        self.recon_losses.append(np.mean(xent_loss))
        self.kl_losses.append(np.mean(kl_loss))
        self.corr_losses.append(cor_loss)

def sampling(args):
    epsilon_std = 1.0
    z_mean, z_log_var = args
    if flat:
        epsilon = K.random_normal(shape=(K.shape(z_mean)[0], K.shape(z_mean)[1]), mean=0.,
                                  stddev=epsilon_std)
    else:
        epsilon = K.random_normal(shape=(K.shape(z_mean)[0], K.shape(z_mean)[1], K.shape(z_mean)[2]), mean=0.,
                                  stddev=epsilon_std)
    return z_mean + K.exp(z_log_var / 2) * epsilon

def corr_loss(y_pred, y_true, be=K):
    cent_pred = y_pred - be.mean(y_pred)
    cent_tr = y_true - be.mean(y_true)

    std_pred = be.std(y_pred)
    std_tr = be.std(y_true)

    return be.mean(cent_pred*cent_tr)/(std_pred*std_tr)

def calc_kl_loss(z_log_var, z_mean):
    return - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)

def build_vae(X_shape, n_experiments, targets, output_ind, layer_sizes=[256], latent_dim=2, batch_size=256, use_corr=True, scale_type=None, flat=False):
    # Encoder network
    if flat:
        x = Input(shape=(X_shape,))
    else:
        x = Input(shape=(n_experiments, X_shape,))
    h = Dense(layer_sizes[0], activation='relu')(x)
    for layer_sz in layer_sizes[1:]:
        h = Dense(layer_sz, activation='relu')(h)
    z_mean = Dense(latent_dim, name='z_mean')(h)
    z_log_var = Dense(latent_dim, name='z_log_var')(h)
    
    # Sample points from latent space
    if flat:
        z = Lambda(sampling, output_shape=(latent_dim,))([z_mean, z_log_var])
    else:
        z = Lambda(sampling, output_shape=(n_experiments,latent_dim,))([z_mean, z_log_var])
    
    # Decoder network
    decoder_h = Dense(layer_sizes[-1], activation='relu')
    decoder_layers = []
    for layer_sz in layer_sizes[::-1][1:]:
        decoder_layers.append(Dense(layer_sz, activation='relu'))
    if not scale_type or 'zero' in scale_type:
        act = 'sigmoid'
    elif 'negone' in scale_type or 'maxabs' in scale_type:
        act = 'tanh'
    else:
        act = 'linear'
    print act
    decoder_mean = Dense(X_shape, activation=act)
    h_decoded = decoder_h(z)
    for decoder_layer in decoder_layers:
        h_decoded = decoder_layer(h_decoded)
    x_decoded_mean = decoder_mean(h_decoded)

    # end-to-end autoencoder
    vae = Model(x, x_decoded_mean)
    if use_corr:
        output_flux = x_decoded_mean[:, :, output_ind]
        experiment_loss = -1 * corr_loss(targets, output_flux)#streaming_pearson_correlation(output_flux, targets)
    xent_loss = x.shape[-1].value * metrics.mean_squared_error(x, x_decoded_mean)
    kl_loss_val = calc_kl_loss(z_log_var, z_mean)
    vae_loss = K.mean(xent_loss + kl_loss_val)
    if use_corr:
        vae_loss += experiment_loss
    #vae_loss = K.sum(vae_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer='rmsprop')
    vae.summary()

    # encoder, from inputs to latent space
    encoder = Model(x, z_mean)

    # generator, from latent space to reconstructed inputs
    if flat:
        decoder_input = Input(shape=(latent_dim,))
    else:
        decoder_input = Input(shape=(n_experiments, latent_dim,))
    _h_decoded = decoder_h(decoder_input)
    for decoder_layer in decoder_layers:
        _h_decoded = decoder_layer(_h_decoded)
    _x_decoded_mean = decoder_mean(_h_decoded)
    generator = Model(decoder_input, _x_decoded_mean)
    return vae, encoder, generator


if __name__ == '__main__':
    args = parser.parse_args()
    print args
    batch_size = args.batch
    n_epochs = args.epochs
    use_corr = args.corr
    flat = not args.resamp
    if flat:
        use_corr = False
    scale = args.scale
    latent_dim = args.dim
    layer_szs = args.layers

    #df = pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    #df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)
    #n_experiments = df.shape[0]

    #dat = np.load('../data/fluxes_resampled.npz')
    #X_train, X_test = dat['train'], dat['test']
    #X_train, y_train, X_test, y_test, btol_col, cols = futils.read_data('../data/flux_samps_2k', scale=scale)

    if args.froot == 'karim' and args.txtl:
        quit()
    fname = '../data/{0}{1}_{2}_fluxes'.format(args.froot, '_txtl' if args.txtl else '', 'stacked' if args.resamp else 'flat')
    X_train, y_train, X_test, y_test, obj_col, cols, y_vals_d = dataset.get_dataset(fname)
    y_vals = np.array(y_vals_d)

    X_shape, n_experiments = X_train.shape[-1], y_vals.shape[0]
    targets = tf.convert_to_tensor(y_vals, dtype=tf.float32)
    vae, encoder, generator = build_vae(X_shape, n_experiments, targets, obj_col.value, layer_szs, latent_dim, batch_size, use_corr, scale, flat)
    es = EarlyStopping(patience=10)
    cbs = [es]
    if use_corr:
        lh = LossHistory()
        cbs.append(lh)
    vae.fit(X_train, shuffle='batch', epochs=n_epochs, batch_size=batch_size,
            validation_data=(X_test, None), callbacks=cbs)
    fname = 'epochs={0}_batch={1}_dimension={2}_corr={3}_scale={4}_froot={5}_txtl={6}_nlayers={7}_resamp={8}.h5'.format(n_epochs, batch_size, latent_dim, use_corr, scale, args.froot, args.txtl, len(layer_szs), args.resamp)
    vae.save('../models/vae_{0}'.format(fname))
    encoder.save('../models/encoder_{0}'.format(fname))
    generator.save('../models/generator_{0}'.format(fname))
    if use_corr:
        with open('../models/losses_{0}'.format(fname), 'w') as f:
            pickle.dump(file=f, obj={'recon_losses': lh.recon_losses, 'kl_losses': lh.kl_losses, 'corr_losses': lh.corr_losses})

    x_test_encoded = encoder.predict(X_test, batch_size=batch_size)
    x_test_gen = generator.predict(x_test_encoded, batch_size=batch_size)

    corrs = []
    for i in range(x_test_gen.shape[0]):
        corr = scipy.stats.pearsonr(x_test_gen[i, :, obj_col.value],)
        corrs.append(corr)
    print 'Mean correlation: ', np.mean(corrs, axis=0)
