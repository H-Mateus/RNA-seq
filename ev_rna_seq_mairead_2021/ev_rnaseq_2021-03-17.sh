#!/bin/bash
pwd
echo $PATH
# set number of cores to one less than the total
CORES=$(($(nproc) - 1))
echo ${CORES}

# First we set up an environment for the project
conda create -n rna_seq
# Then we can install the software we need
conda install python=3.8
conda install -c bioconda fastqc
conda install -c bioconda multiqc
conda install -c bioconda cutadapt
conda install -c bioconda trim-galore
conda install salmon=1.4.0

# Then we activate the environment
conda activate rna_seq
# Then set up our channels - the last channel added is the highest priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# check conda info
conda info
# list packages
conda list

# cd Raw
ls -l -h
  for dir in ./Raw/*/
  do
      echo $dir
      mv ${dir}*.gz ~/RNA_seq_mairead_2021-03-16/Raw
  done

ls -l -h

# generate fastqc reports
fastqc -o Raw/fastqc_reports Raw/*fastq.gz

# aggregate reports
multiqc -o Raw/fastqc_reports Raw/fastqc_reports

# trim reads - set overlap to 3 and length to 15 to match Liverpool settings
trim_galore -a " AGATCGGAAGAGCACACGTCT -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG -a CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA -a TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG" --stringency 3 --cores 8 --length 15 -o trimmed_manual Raw/*.fastq.gz

# make directory for trimming reports
mkdir trimmed_manual/trimming_reports
# move reports
mv trimmed_manual/*.txt trimmed_manual/trimming_reports

# generate fastqc reports on trimmed data
fastqc -o trimmed_manual/fastqc_reports trimmed_manual/*.fq.gz
# use multiqc
multiqc -o trimmed_manual/fastqc_reports trimmed_manual/fastqc_reports

# set up salmon index - using gencode v37 transcripts
cd reference_genes
salmon index -t gencode_human.v37.transcripts.fa -i gencode_v37_index --gencode
cd ..

# make .txt with the sample strings
cd trimmed_manual
command ls *.fq.gz | while read file;do echo $file | cut -dR -f1 >> samples.txt; done

# get a list of unique sample file strings
SAMPLE=$(cat samples.txt | sort | uniq)
echo $SAMPLE
# loop over sample files with salmon
for i in ${SAMPLE}
do
    salmon quant -i ../reference_genes/gencode_v37_index -l A -o ${i} -r ${i}*
done

# check that we have 17 directories as we expect
ls -d *L001_ | wc -l

cd 1-16081_Nor_Cells_AGTCAAAT_L001_
  command ls -lh
  # check counts and quant
  cat lib_format_counts.json
  head quant.sf
# check mapping rate
grep -i 'mapping rate' logs/salmon_quant.log
  cd ..

# get ensamble gene and transcript ids
grep -P -o 'ENST\d{11}' ../reference_genes/gencode_human.v37.transcripts.fa > enst.txt
grep -P -o 'ENSG\d{11}' ../reference_genes/gencode_human.v37.transcripts.fa > ensg.txt

# check
head enst.txt
paste -d ',' enst.txt ensg.txt | head
# past gene and transcipt ids
paste -d ',' enst.txt ensg.txt > ../../git_work/RNA-seq/ev_rna_seq_mairead_2021/gene_map.csv
