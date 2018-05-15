import pandas as pd
from cobra.flux_analysis import sample
import multiprocessing
import cobra
import utils

if __name__ == '__main__':
    model = cobra.io.read_sbml_model('../models/ecoli_cf_base.sbml')

    df = pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)

    for i, row in df.iterrows():
	print i
	model_i = utils.add_reagents_to_model(model, row)
	samples_i = sample(model_i, 2000, processes=multiprocessing.cpu_count() - 1)
	samples_i.to_csv('../data/fluxes_{0}'.format(i))
