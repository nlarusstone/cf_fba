#scale=('norm' 'robust' 'maxabs' 'negone' 'zero')
#direction=('exp')
#latents=(3 5 7 10 20 50)
froots=('echo' 'hand' 'karim')
latents=(2 10)
#froots=('echo')
#latents=(2)
layers1=(1024 1024 1024)
#layers1=(1024 512 256)
#layers2=(1024 256)
typeset -A layers
layers=([l1]="${layers1[@]}")
# [l2]="${layers2[@]}")
# 'exp')
#for sc in "${scale[@]}"; do
#	for dirc in "${direction[@]}"; do
#		python vae.py -d 2 -b 256 -n 100 -s ${dirc}_${sc}
##	done
#done
#python vae.py -d 2 -b 256 -n 100
for froot in "${froots[@]}"; do
    for latent in "${latents[@]}"; do
        for layer in "${!layers[@]}"; do
            python vae.py -d ${latent} -b 256 -n 200 --layers ${layers[$layer]} --no-resamp --no-txtl -f ${froot} &
            python vae.py -d ${latent} -b 256 -n 200 --layers ${layers[$layer]} --resamp --no-txtl -f ${froot} &
            python vae.py -d ${latent} -b 256 -n 200 --layers ${layers[$layer]} --no-resamp --txtl -f ${froot} &
            python vae.py -d ${latent} -b 256 -n 200 --layers ${layers[$layer]} --resamp --txtl -f ${froot} &
        done
    done
done
#python vae.py -d 2 -b 128 -n 100 -s 'flux_zero'
#python vae.py -d 2 -b 256 -n 100 -s 'flux_zero' --no-corr
#python vae.py -d 2 -b 256 -n 500 -s 'flux_zero' 
