import json
import os

def convert_rxn(julia_rxn, debug_dict):
    rxn_name = julia_rxn['reaction_name']
    debug_entry = debug_dict[rxn_name]
    metabs = {react["symbol"]: int(-1 * float(react["stoichiometry"])) for react in julia_rxn["list_of_reactants"]}
    metabs.update({product["symbol"]: int(float(product["stoichiometry"])) for product in julia_rxn["list_of_products"]})
    cobra_rxn = {
        'id': rxn_name,
        'name': rxn_name,
        'lower_bound': debug_entry['lower_bound'],
        'upper_bound': debug_entry['upper_bound'],
        'metabolites': metabs
    }
    return cobra_rxn

def convert_metab(julia_metab, atom_dict):
    sym = julia_metab['species_symbol']
    cobra_metab = {
        'compartment': julia_metab['species_compartment'],
        'id': sym,
        'name': sym,
        'formula': atom_dict[sym] if sym in atom_dict else None
    }
    return cobra_metab

def parse_atom_file(f_atom):
    atom_dict = {}
    with open(f_atom, 'r') as f:
        for line in f.readlines():
            l = line.strip().split(',')
            if len(l) < 7:
                continue
            metab_name = l[0]
            metab_formula = 'C{0}H{1}N{2}O{3}P{4}S{5}'.format(*[l[i] for i in range(1, 7)])
            #C,H,N,O,P,S
            atom_dict[metab_name] = metab_formula
    return atom_dict

def parse_bounds(f_bounds, rxn_dict, idx_map):
    with open(f_bounds, 'r') as f_bounds:
        if f_bounds == 'DataDictionary.jl':
            bounds = False
            for i, l in enumerate(f_bounds.readlines()):
                line = l.strip().split() 
                if '];' in line:
                    bounds = False
                if bounds:
                    lb, ub = float(line[0]), float(line[1])
                    idx = int(line[6]) - 1
                    rxn_name = idx_map[i]
                    rxn_dict[rxn_name]['lower_bound'] = lb 
                    rxn_dict[rxn_name]['upper_bound'] = ub
                if 'default_bounds_array' in line:
                    bounds = True
        else:
            for idx, l in enumerate(f_bounds.readlines()):
                line = l.strip().split('\t') 
                lb, ub = map(float, line)
                rxn_name = idx_map[idx]
                rxn_dict[rxn_name]['lower_bound'] = lb 
                rxn_dict[rxn_name]['upper_bound'] = ub
    return rxn_dict
            
def parse_debug_file(f_debug, f_bounds):
    rxn_dict = {}
    idx_map = {}
    with open(f_debug, 'r') as f_debug:
        for line in f_debug.readlines():
            try:
                idx, rxn = line.strip().split(' ', 1)
            except ValueError:
                break
            rxn_name, rxn_str = rxn.split('::')
            # Julia is 1 indexed
            idx_map[int(idx) - 1] = rxn_name
            rxn_dict[rxn_name] = {'formula': rxn_str}
    rxn_dict = parse_bounds(f_bounds, rxn_dict, idx_map)
    return rxn_dict

def julia_to_cobra_json(folder, fout=None):
    f_atom = os.path.join(folder, 'Atom.txt')
    atom_dict = parse_atom_file(f_atom)
    f_debug = os.path.join(folder, 'Debug.txt')
    f_bounds = os.path.join(folder, 'Bounds.txt') #DataDictionary.jl#Bounds.txt
    rxn_dict = parse_debug_file(f_debug, f_bounds)
    f_model = os.path.join(folder, 'Reactions.json')
    with open(f_model, 'r') as f:
        julia_model = json.load(f)

    cobra_model = {'id': 'cell_free_varner'}
    cobra_metabs, cobra_rxns = [], []
    cmpts = set()
    for julia_metab in julia_model["list_of_species"]:
        cobra_metab = convert_metab(julia_metab, atom_dict)
        cmpts.add(cobra_metab['compartment'].lower())
        cobra_metabs.append(cobra_metab)
    for julia_rxn in julia_model["list_of_reactions"]:
        cobra_rxn = convert_rxn(julia_rxn, rxn_dict)
        cobra_rxns.append(cobra_rxn)
    # NOTE: we put everything in cytoplasm right now
    cmpts = set(['cytoplasm'])
    cobra_model['metabolites'] = cobra_metabs
    cobra_model['reactions'] = cobra_rxns
    cobra_model['genes'] = []
    cobra_cmpts = {}
    for cmpt in cmpts:
        cobra_cmpts[cmpt[0]] = cmpt.capitalize()
    cobra_model['compartments'] = cobra_cmpts
    if fout:
        with open(fout, 'w') as outfile:
            json.dump(cobra_model, outfile)
    return cobra_model

julia_to_cobra_json('../3rd_Party_Code/Sequence-Specific-FBA-CFPS-Publication-Code/Staging/', '../models/varner_rfp.json')
