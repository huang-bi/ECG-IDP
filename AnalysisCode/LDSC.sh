# /bin/bash

inputGWAS1=$1
inputGWAS2=$2
outputPrefix=$3
eur_w_ld_chr='~/eur_w_ld_chr'

ldsc.py --rg $inputGWAS1,$inputGWAS2 --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $outputPrefix
