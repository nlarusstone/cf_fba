cd ../src
# IO: cf_io
# Create base models: convert.py
#python convert.py --no-txtl &&
# Flux utils: fba_utils
# Flux sample (creates experimental condition models and samples): flux_sample.py
#python flux_sample.py --no-txtl -m iJO1366 -d nls -f manual -r 5 -a 1 -b 50 -s 2000 && 
# Create dataset: datasets.py
#python datasets.py -r  --no-txtl -f manual -d nls &&
# Run vae: vae.py
#python vae.py -d 2 -b 256 -n 20 --layers 1024 1024 1024 --resamp --no-txtl -f manual
# Reduce model: threshold.py
python threshold.py -d 2 -n 20 --corr -l 1024 1024 1024 -f manual -s nls -m iJO1366 --no-txtl --resamp
cd -
