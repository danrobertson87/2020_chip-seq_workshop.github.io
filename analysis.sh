##DEEPTOOLS

##Correlation
multiBamSummary bins -b bwa_out/*/*uniq.bam -bs 10000 -o cov.matrix -e 200 --smartLabels
plotCorrelation -in cov.matrix -p heatmap -c spearman -o spearman_cor.png
##Fingerprint
#Leave commented for tutorial as takes a long time to run
#plotFingerprint -b bwa_out/*/*uniq.bam -e 200 --smartLabels -o fingerprint.png 
##Normalised ChIP vs Input scores
bamCompare -b1 bwa_out/Reb1_R1/Reb1_R1.uniq.bam -b2 bwa_out/Input_R1/Input_R1.uniq.bam -o Reb1_R1.norm.bw -e 200 --operation log2 --normalizeUsing BPM --scaleFactorsMethod None
bamCompare -b1 bwa_out/Reb1_R2/Reb1_R2.uniq.bam -b2 bwa_out/Input_R2/Input_R2.uniq.bam -o Reb1_R2.norm.bw -e 200 --operation log2 --normalizeUsing BPM --scaleFactorsMethod None
ln -s $PWD/*norm.bw visualisation
##Get gene annotations
cp /homes/genomes/s.cerevisiae/sacCer3/annotation/UCSC_sgdGene.bed .
##TSS plot
computeMatrix reference-point -R UCSC_sgdGene.bed -S Reb1_R1.norm.bw Reb1_R2.norm.bw -o TSS.matrix --referencePoint TSS --upstream 2000 --downstream 2000 -bs 100 --smartLabels
plotHeatmap -m TSS.matrix -o TSS.heatmap.png

##PEAK CALLING

##macs
macs2 callpeak -t bwa_out/Reb1_R1/Reb1_R1.uniq.bam bwa_out/Reb1_R2/Reb1_R2.uniq.bam -c bwa_out/Input_R1/Input_R1.uniq.bam bwa_out/Input_R2/Input_R2.uniq.bam -g 12000000 --nomodel -n Reb1
##Get sequences under peaks
fastaFromBed -fi /homes/genomes/s.cerevisiae/sacCer3/sacCer3.fa -bed Reb1_peaks.narrowPeak -fo Reb1_peaks.fasta
##MEME
meme-chip Reb1_peaks.fasta
