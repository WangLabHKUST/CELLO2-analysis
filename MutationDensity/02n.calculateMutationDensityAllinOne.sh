#! /bin/bash
bed=$1 #e.g. U87.MYC.SRX129069.05.bed



pre=${bed%.*}
#pre=${pre/U87./}
echo -e $pre"\t"$bed

echo -e "Type\tseed\tbin\tMuts\tWidth" >MutationDensity.per150000.${pre}.txt
for sdi in {123,42,5,666,7890,1,2,22,3,333}; do
  echo $sdi

  for ii in {1..9} $(seq 10 10 90) $(seq 100 100 900) $(seq 1000 1000 10000) $(seq 20000 20000 100000) $(seq 500000 500000 2000000) {5000000,10000000}; do
    i=$(printf "%.0f" $ii)
    echo $i
    if [ ! -f ${pre}.d${i}.bed ]; then
      awk -F "\t" -v ix=$i '{print $1"\t"$2-ix"\t"$3+ix"\t"$4"\t"$5"\t"$6}' ${bed} |awk -F "\t"  '{if ($2<0) print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6; else print $0}' >${pre}.d${i}.raw.bed
      bedtools sort -i ${pre}.d${i}.raw.bed |awk '{print $1"\t"$2"\t"$3"\trecord"NR}'>${pre}.d${i}.sorted.bed
      bedtools merge -i ${pre}.d${i}.sorted.bed -c 4 -o collapse >${pre}.d${i}.merge.bed
      bedtools intersect -a ${pre}.d${i}.merge.bed -b hg19.mutations.border.bed >${pre}.d${i}.bed
      rm ${pre}.d${i}.raw.bed ${pre}.d${i}.sorted.bed ${pre}.d${i}.merge.bed
    fi

    x=$(bedtools intersect -a PairedGlioma.WGS.mutations.HM.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 2>/dev/null |wc -l);
    y=$(awk '{sum+=$3-$2} END {print sum}' ${pre}.d${i}.bed);
    echo -e "HM\t"$sdi"\t"$i"\t"$x"\t"$y >>MutationDensity.per150000.${pre}.txt

    x=$(bedtools intersect -a PairedGlioma.WGS.mutations.HMMYC.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 2>/dev/null |wc -l);
    y=$(awk '{sum+=$3-$2} END {print sum}' ${pre}.d${i}.bed);
    echo -e "HMMYC\t"$sdi"\t"$i"\t"$x"\t"$y >>MutationDensity.per150000.${pre}.txt

    x=$(bedtools intersect -a PairedGlioma.WGS.mutations.HMnoMYC.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 2>/dev/null |wc -l);
    y=$(awk '{sum+=$3-$2} END {print sum}' ${pre}.d${i}.bed);
    echo -e "HMnoMYC\t"$sdi"\t"$i"\t"$x"\t"$y >>MutationDensity.per150000.${pre}.txt

    x=$(bedtools intersect -a PairedGlioma.WGS.mutations.NHM.ix.sample150000.seed${sdi}.bed -b ${pre}.d${i}.bed 2>/dev/null |wc -l);
    y=$(awk '{sum+=$3-$2} END {print sum}' ${pre}.d${i}.bed);
    echo -e "NHM\t"$sdi"\t"$i"\t"$x"\t"$y >>MutationDensity.per150000.${pre}.txt
  done
done
echo "Done. Results in MutationDensity.per150000.${pre}.txt"
echo "Run rm ${pre}.d[1-9]*.bed to remove unnecessary files."
