samps=(2000)
# 200 10)
for samp in "${samps[@]}"; do
    python flux_sample.py -s ${samp} -d 'karim' -f 'karim' -r 25 -a 0 -b 25 --no-txtl
done
for samp in "${samps[@]}"; do
    python flux_sample.py -s ${samp} -f 'echo' -r 2 -a 0.5 --no-txtl
    python flux_sample.py -s ${samp} -f 'echo' -r 2 -a 0.5 --txtl
done
for samp in "${samps[@]}"; do
    python flux_sample.py -s ${samp} -f 'manual' -r 5 -a 1 --no-txtl
    python flux_sample.py -s ${samp} -f 'manual' -r 5 -a 1 --txtl
done
