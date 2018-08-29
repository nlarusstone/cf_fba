cd ../src
froots=('echo' 'manual' 'karim')
latents=(10)
layers1=(1024 1024 1024)
layers2=(1024 512 256)
typeset -A layers
layers=([l1]="${layers1[@]}" [l2]="${layers2[@]}")
for froot in "${froots[@]}"; do
    for latent in "${latents[@]}"; do
        for layer in "${!layers[@]}"; do
            python threshold.py -d ${latent} --layers ${layers[$layer]} --no-txtl -f ${froot} &
            python threshold.py -d ${latent} --layers ${layers[$layer]} --txtl -f ${froot} &
        done
    done
done
#python threshold.py -d 2 --layers 1024 1024 1024
cd -
