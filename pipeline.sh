##PREPARATION

##Create a web directory to view output
ln -s $PWD ~/public_html/
##Create a folder for fastq files
mkdir fastq
##Copy the raw data to the fastq folder
cp /homes/library/training/ChIP-seq_workshop/data/*fq.gz fastq/.


##QC

##FastQ Screen
fastq_screen  --conf /homes/genomes/tool_configs/fastq_screen/fastq_screen.conf fastq/$1.fq.gz --outdir fastq
##FastQC
fastqc fastq/$1.fq.gz


##PREPROCESSING

##Cutadapt
cutadapt -a AGATCGGAAGAG -q 20 --minimum-length 36 -o fastq/$1.trim.fq.gz fastq/$1.fq.gz > fastq/$1.trim.cutadapt_report.txt

##FastQC
fastqc fastq/$1.trim.fq.gz


##ALIGNMENT

##BWA
mkdir bwa_out
mkdir bwa_out/$1
bwa mem -t 5 -a -R '@RG\tID:{}\tPL:ILLUMINA' -M /homes/genomes/s.cerevisiae/sacCer3/bwa_indexes/sacCer3 fastq/$1.trim.fq.gz > bwa_out/$1/$1.sam
##Sort and index
samtools sort bwa_out/$1/$1.sam -o bwa_out/$1/$1.bam -T bwa_out/$1/$1 -O BAM; samtools index bwa_out/$1/$1.bam
##Filter for uniquely mapped reads
samtools view -b -q 20 bwa_out/$1/$1.bam -o bwa_out/$1/$1.uniq.bam; samtools index bwa_out/$1/$1.uniq.bam
##FastQC on BAM files
fastqc bwa_out/$1/$1.uniq.bam


##VISUALISATION

##Convert BAM to normalised bigWig
bamCoverage -bs 1 --normalizeUsing BPM -e 200 -b bwa_out/$1/$1.uniq.bam --outFileName bwa_out/$1/$1.uniq.bw
