# /bin/bash

inputGWAS=$1
outputPrefix=$2
eur_w_ld_chr='~/eur_w_ld_chr'

ldsc.py --h2 $inputGWAS --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $outputPrefix

