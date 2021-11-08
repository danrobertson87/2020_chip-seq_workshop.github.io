##PREPARATION

##Create a web directory to view output
ln -s $PWD ~/public_html/
##Create a folder for fastq files
mkdir fastq
##Copy the raw data to the fastq folder
cp /homes/library/training/ChIP-seq_workshop/data/*fq.gz fastq/.

##QC

##FastQ Screen
cat samples.txt | parallel -j 4 "fastq_screen  --conf /homes/genomes/tool_configs/fastq_screen/fastq_screen.conf fastq/{}.fq.gz --outdir fastq"
##FastQC
cat samples.txt | parallel -j 4 fastqc fastq/{}.fq.gz
##Multi QC
multiqc -o fastq fastq

##PREPROCESSING

##Cutadapt
cat samples.txt | parallel -j 4 "cutadapt -a AGATCGGAAGAG -q 20 --minimum-length 36 -o fastq/{}.trim.fq.gz fastq/{}.fq.gz > fastq/{}.trim.cutadapt_report.txt"


##FastQC
cat samples.txt | parallel -j 4 fastqc fastq/{}.trim.fq.gz
multiqc -f -o fastq fastq

##ALIGNMENT

##BWA
mkdir bwa_out
cat samples.txt | parallel -j 4 mkdir bwa_out/{}
cat samples.txt | parallel -j 4 "bwa mem -t 5 -a -R '@RG\tID:{}\tPL:ILLUMINA' -M /homes/genomes/s.cerevisiae/sacCer3/bwa_indexes/sacCer3 fastq/{}.trim.fq.gz > bwa_out/{}/{}.sam"
##Sort and index
cat samples.txt | parallel -j 4 "samtools sort bwa_out/{}/{}.sam -o bwa_out/{}/{}.bam -T bwa_out/{}/{} -O BAM; samtools index bwa_out/{}/{}.bam"
##Filter for uniquely mapped reads
cat samples.txt | parallel -j 4 "samtools view -b -q 20 bwa_out/{}/{}.bam -o bwa_out/{}/{}.uniq.bam; samtools index bwa_out/{}/{}.uniq.bam"
##FastQC on BAM files
cat samples.txt | parallel -j 4 fastqc bwa_out/{}/{}.uniq.bam
multiqc -o bwa_out bwa_out

##VISUALISATION

##Convert BAM to normalised bigWig
cat samples.txt | parallel -j 4 bamCoverage -bs 1 --normalizeUsing BPM -e 200 -b bwa_out/{}/{}.uniq.bam --outFileName bwa_out/{}/{}.uniq.bw

##Make visualisation directory
mkdir visualisation
ln -s $PWD/bwa_out/*/*bw visualisation
ln -s $PWD/bwa_out/*/*uniq.bam* visualisation

##Tidy up
rm fastq/*trim*fq.gz #Remove trimmed fastq temp files
rm bwa_out/*/*.sam #Remove sam files

