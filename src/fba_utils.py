import re
import cobra

def get_aa_metab(model, aa, cmpt='c'):
    return model.metabolites.query('{0}__._{1}'.format(aa, cmpt))

def replace_metab(model, metab):
    new_id = re.sub(r'_.$', '_c', metab.id)
    try:
        cyt = model.metabolites.get_by_id(new_id)
    except:
        cyt = metab
        cyt.id = new_id
        cyt.compartment = 'c'
    return cyt

# Changes objective function of model
# Works in place
def change_obj(model, metab=None, rxn=None):
    if rxn is not None:
        try:
            model.objective = model.reactions.get_by_id(rxn)
        except:
            raise Exception('Objective reaction does not exist')
    elif metab is not None:
        try:
            met = model.metabolites.get_by_id(metab)
        except:
            raise Exception('Metabolite does not exist')
        if (not 'DM_{0}'.format(met.id) in model.reactions):
            model.add_boundary(metabolite=met, type='demand')
            model.objective = model.reactions.get_by_id('DM_{0}'.format(met.id))
    else:
        print 'Need to specify either a metabolite or reaction as an objective'
        raise Exception
