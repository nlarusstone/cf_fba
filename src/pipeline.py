# IO: cf_io
# Create base models: convert.py
python convert.py --txtl -o PROTEIN_export_RFP
# Flux utils: fba_utils
# Flux sample (creates experimental condition models and samples): flux_sample.py
python flux_sample.py --txtl -m iJO1366 -d nls -f manual -r 5 -a 1 -b 50 -s 2000
# Create dataset: datasets.py
python datasets.py -r  --txtl -f manual -d nls
# Run vae: vae.py
python vae.py -d 2 -b 256 -n 200 --layers 1024 1024 1024 --resamp --txtl -f manual
# Test correlation: correlation.py
