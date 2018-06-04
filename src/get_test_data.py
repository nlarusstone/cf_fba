from keras.models import load_model
import create_dataset as dataset
import numpy as np

#n_epochs = 200
#batch_size = 256
#latent_dim = 10
#use_corr = True
#scale = 'flux_zero'
#froot = 'hand'
#txtl = True
#resamp = True
#flat = not resamp
#layer_szs = [1024, 512, 256]

def get_test_data(froot, txtl, resamp, latent_dim, layer_szs, use_corr=True, n_epochs=200, batch_size=256, scale='flux_zero'):
    flat = not resamp
    if not resamp:
        use_corr = False
    if froot == 'karim':
        assert(not txtl)
    fname = 'epochs={0}_batch={1}_dimension={2}_corr={3}_scale={4}_froot={5}_txtl={6}_nlayers={7}_resamp={8}{9}.h5'.format(
        n_epochs, batch_size, latent_dim, use_corr, scale, froot, txtl, len(layer_szs), resamp, '_lastlayer=1024' if layer_szs[-1] == 1024 else '')
    print('Load models {0}'.format('../models/encoder_' + fname))
    encoder = load_model('../models/encoder_{0}'.format(fname))
    generator = load_model('../models/generator_{0}'.format(fname))

    data_f = '../data/{0}{1}_{2}_fluxes'.format(froot, '_txtl' if txtl else '', 'stacked' if resamp else 'flat')
    X_train, y_train, X_test, y_test, obj_col, cols, y_vals_d = dataset.get_dataset(data_f)
    print 'Read in data from {0}'.format(data_f)
    y_vals = np.array(y_vals_d)
    X_test = np.array(X_test)
    y_test = np.array(y_test)

    test_enc = encoder.predict(X_test)
    print 'Encoded data'
    test_dec = generator.predict(test_enc)
    print 'Decoded data'

    return encoder, generator, X_test, y_test, obj_col, cols, y_vals_d, test_enc, test_dec
