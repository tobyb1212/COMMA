#!/bin/bash

#cd /COMMA/data/meta_analysis/fastq/rna_linna_kuosmanen
##########
# get /COMMA/data/meta_analysis/fastq/rna_linna_kuosmanen/download_unique_samples.sh
# get *R*

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate interstellar

cd /COMMA/data/meta_analysis/fastq/rna_linna_kuosmanen/

for encrypted_file in *.c4gh; do
    decrypted_file="${encrypted_file%.c4gh}" 
    crypt4gh decrypt --sk "$CRYPT4GH_KEY" < "$encrypted_file" > "$decrypted_file"
    echo "Decrypted $encrypted_file"
done

cd /COMMA/data/meta_analysis/fastq/spatial_linna_kuosmanen/

for encrypted_file in *.c4gh; do
    decrypted_file="${encrypted_file%.c4gh}" 
    crypt4gh decrypt --sk "$CRYPT4GH_KEY" < "$encrypted_file" > "$decrypted_file"
    echo "Decrypted $encrypted_file"
done
