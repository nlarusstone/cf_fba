scale=('norm' 'robust' 'maxabs' 'negone' 'zero')
direction=('exp')
latents=(3 5 7 10 20 50)
# 'exp')
#for sc in "${scale[@]}"; do
#	for dirc in "${direction[@]}"; do
#		python vae.py -d 2 -b 256 -n 100 -s ${dirc}_${sc}
##	done
#done
#python vae.py -d 2 -b 256 -n 100
for latent in "${latents[@]}"; do
    python vae.py -d ${latent} -b 256 -n 100 -s 'flux_zero'
done
python vae.py -d 2 -b 128 -n 100 -s 'flux_zero'
python vae.py -d 2 -b 256 -n 100 -s 'flux_zero' --no-corr
python vae.py -d 2 -b 256 -n 500 -s 'flux_zero' 
