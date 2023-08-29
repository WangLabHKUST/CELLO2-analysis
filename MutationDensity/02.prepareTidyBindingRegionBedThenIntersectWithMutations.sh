#! /bin/bash
bed=$1 #e.g. U87.MYC.SRX129069.05.bed

###prepare the regions near the given peaks, then intersect with mutations


pre=${bed%.*}
pre=${pre/U87./}
echo -e $pre"\t"$bed
for sdi in {123,42,5,666,7890}; do
  echo $sdi
  for ii in {1..9} $(seq 10 10 90) $(seq 100 100 900) $(seq 1000 1000 10000) $(seq 20000 20000 100000) $(seq 500000 500000 2000000) {5000000,10000000}; do
    i=$(printf "%.0f" $ii)
    echo ${i}
    if [ ! -f ${pre}.d${i}.bed ]; then
      awk -F "\t" -v ix=$i '{print $1"\t"$2-ix"\t"$3+ix"\t"$4"\t"$5"\t"$6}' ${bed} |awk -F "\t"  '{if ($2<0) print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6; else print $0}' >${pre}.d${i}.raw.bed
      bedtools merge -i ${pre}.d${i}.raw.bed -c 4 -o collapse >${pre}.d${i}.merge.bed
      bedtools intersect -a ${pre}.d${i}.merge.bed -b hg19.mutations.border.bed >${pre}.d${i}.bed
      rm ${pre}.d${i}.raw.bed ${pre}.d${i}.merge.bed
    fi
    bedtools intersect -wb -a PairedGlioma.WGS.mutations.HMMYC.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 1>HM.MYCgain.mutations.chr.sample150000.seed${sdi}.${pre}.dist${i}.txt; 2>/dev/null
    bedtools intersect -wb -a PairedGlioma.WGS.mutations.HMnoMYC.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 1>HM.noMYCgain.mutations.chr.sample150000.seed${sdi}.${pre}.dist${i}.txt; 2>/dev/null
    bedtools intersect -wb -a PairedGlioma.WGS.mutations.NHM.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 1>nonHM.mutations.chr.sample150000.seed${sdi}.${pre}.dist${i}.txt; 2>/dev/null
    #bedtools intersect -wb -a PairedGlioma.WGS.mutations.HM.ix.bed -b ${pre}.d${i}.bed 1>HM.mutations.chr.${pre}.dist${i}.txt; 2>/dev/null
  done
done
