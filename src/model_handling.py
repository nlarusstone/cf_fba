import pickle
import itertools

def compare_losses(category, opts):
    cats = ['epochs', 'batch', 'dimension', 'corr', 'scale']
    dvs = [100, 256, 2, True, 'flux_zero']
    idx = cats.index(category)
    del cats[idx]
    del dvs[idx]
    defdict = dict(zip(cats, dvs))
    for opt in opts:
        curdict = defdict.copy()
        curdict.update({category: opt})
        with open('../models/losses_epochs={epochs}_batch={batch}_dimension={dimension}_corr={corr}_scale={scale}.h5'.format(
            **curdict
        ), 'r') as f:
            losses = pickle.load(f)
            print opt, [ls[-1] for k, ls in losses.items()]

#scale_type = ['norm', 'robust', 'maxabs', 'negone', 'zero']
#dirs = ['flux', 'exp']
#compare_losses('scale', ["{0}_{1}".format(d, st) for d, st in itertools.product(dirs, scale_type)])
#compare_losses('dimension', [2, 3, 5, 7, 10])
