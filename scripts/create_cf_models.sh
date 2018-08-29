cd ../src
#python convert.py --no-txtl
#python convert.py --txtl -o PROTEIN_export_RFP
#python convert.py -c '../data/karim_concs.csv' -s 'karim' --no-txtl -o btol_c
#model_l=('ColiPruned.xml' 'ColiCore.xml' 'ColiGS.xml' 'ColiPrunedComp_withoutConsRel.xml')
model_l=('ColiPruned.xml' 'ColiGS.xml' 'ColiPrunedComp_withoutConsRel.xml')
for model_f in "${model_l[@]}"; do
    python convert.py --no-txtl -m ../bio_models/${model_f} -o Ec_biomass_iAF1260_core_59p81M
    #python convert.py --txtl -m '../bio_models/${model_f}'
    python convert.py --no-txtl -m ../bio_models/${model_f} -c '../data/karim_concs.csv' -s 'karim' -o btol_c
done
cd -
