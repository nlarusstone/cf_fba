import pandas as pd
import Bio.PDB.Polypeptide as pdb
import Bio.SeqUtils as su

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

if __name__ == '__main__':
    aa_mix = pd.read_csv('../data/aa_mix.csv', index_col='AA')
    nrg_mix = pd.read_csv('../data/energy_mix.csv', index_col='compound')
    with open('../genes/rfp.txt', 'r') as f:
	seq = f.read()
    cfps_conc = pd.read_csv('../data/cfps_start.csv', index_col='compound')

    aa_mix['start_conc'] = aa_mix.apply(lambda row: calc_conc(row['weight_add'], 1, 
							      seq=pdb.three_to_one(row.name.upper()), seq_type='protein'), axis=1)
    aa_mix['conc_add'] = conc_dilution(aa_mix['start_conc'], aa_mix['vol_add'], aa_mix['vol_add'].sum())
    nrg_mix['start_conc'] = nrg_mix.apply(lambda row: calc_conc(row['amt'], row['fill'], mw=row['mw']), axis=1)
    nrg_mix['conc_add'] = conc_dilution(nrg_mix['start_conc'], nrg_mix['vol_add'], nrg_mix['vol_add'].sum())
    dna_conc = calc_conc(0.000750, 0.00496, seq=seq, seq_type='DNA')

    for cmpnd, vals in nrg_mix.iterrows():
	cfps_conc.loc[cmpnd] = [vals['conc_add'], 5]
    for aa, vals in aa_mix.iterrows():
	cfps_conc.loc[aa] = [vals['conc_add'], 10]
    cfps_conc.loc['GENE'] = [dna_conc, 4.96]
    cfps_conc['final_conc'] = conc_dilution(cfps_conc['start_conc'], cfps_conc['amt'], 50)
    return cfps_conc
