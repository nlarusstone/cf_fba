cd ../src
froots=('echo' 'manual' 'karim')
#for froot in "${froots[@]}"; do
#    python create_dataset.py -f ${froot} --no-resamp --no-txtl &
#    python create_dataset.py -f ${froot} --no-resamp --txtl &
#done
for froot in "${froots[@]}"; do
    python create_dataset.py -f ${froot} --resamp --no-txtl &
    python create_dataset.py -f ${froot} --resamp --txtl &
done
cd -
