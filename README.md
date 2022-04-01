# Integrated-RNA-seq-and-ChIP-seq-analysis

Here, I provide the code used to perform an integrated RNA-seq and ChIP-seq analysis of data accessed through GEO, producing a list of transcription factor target genes that5 were used in the published paper, [Genes regulated by BCL11B during T-cell development are enriched for de novo mutations found in schizophrenia patients.](https://onlinelibrary.wiley.com/doi/abs/10.1002/ajmg.b.32811)

I am interested in the target genes of tracription factor bcl11b. For the analysis below, I used RNA-seq and ChIP-seq data from [Ha VL, Luong A, Li F, Casero D et al., 2017](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84678), carried out on human T-cells in DN3/DN4 stage at differentiation.

### RNA-seq Data

SRR5194388 - Control CD1a+ replicate 1

SRR5194390 - Knockdown CD1a+ replicate 1

SRR5194392 - Control CD1a+ replicate 2

SRR5194394 - Knockdown CD1a+ replicate 2

SRR5194396 - Control CD1a+ replicate 3

SRR5194398 - Knockdown CD1a+ replicate 3

CD1a+ = marker of DN3 and DN4 T-cell differentiation.

### ChIP-seq data.

SRR3938841.sra - CD34- cells-replicate 1-IP

SRR3938842.sra -  CD34- cells-replicate 1-control

SRR3938845.sra - CD34- cells-replicate 2-IP

SRR3938846.sra - CD34- cells-replicate 2-control

CD34- = marker of DN4 T-cell differentiation.


## RNA-seq Analysis

#### 1. Fetch SRA Files


```

for i in {1..6}
do
/home/nextgen2015/.aspera/connect/bin/ascp -i /home/nextgen2015/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -Tr -l200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP079/SRP079176/SRR393884$i/SRR393884$i.sra ./
done

```
#### 2. Convert SRA files to fastq format


```

fastq-dump --skip-technical --split-files --readids --read-filter pass --dumpbase --clip *.sra

```

#### 3. Carry out FASTQC and then trim based on fastqc results

```

fastqc *.fastq

```

```

for i in {88,90,92,94,96,98};
do
java -jar /home/liam/bin/trimmomatic-0.33.jar PE -phred33 SRR51943${i}_pass_1.fastq SRR51943${i}_pass_2.fastq SRR51943${i}_output_forward_paired.fq SRR51943${i}_forward_unpaired.fq SRR51943${i}_reverse_paired.fq SRR51943${i}_reverse_unpaired.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:15 MINLEN:70;
done

```

#### 4. Alignment using hisat2


##### Index

```

wget "ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

wget "ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"

hisat2-build Homo_sapiens.GRCh37.75.dna.primary_assembly.fa index

# make splicesites file:

python extract_splice_sites.py Homo_sapiens.GRCh37.75.gtf > splicesites.txt

```
##### Align back to genome

```

for i in {88,90,92,94,96,98}
do
/data4/laura/hisat_stranded/hisat2-2.0.5/hisat2 -x index --known-splicesite-infile splicesites.txt -p 12 --dta -1 SRR51943${i}_output_forward_paired.fq -2 SRR51943${i}_reverse_paired.fq -S SRR51943${i}.sam
done

```

#### 5. Post alignment processing using samtools

```

for i in *.sam
do
id=$(echo ${i} | sed 's/.sam//')
# Convert SAM files to BAM format
samtools view -Sb ${id}.sam > ${id}.bam
# Remove possible PCR duplicates
samtools rmdup ${id}.bam ${id}.rmdup.bam
# Sort BAM file
samtools sort ${id}.rmdup.bam ${id}.rmdup.sorted
# index the bam file to creat .bai files
samtools index ${id}.rmdup.sorted.bam
# Generate file containing mapping stats
samtools flagstat ${id}.rmdup.sorted.bam > ${id}_mappingstats.txt
done

```

#### 6. Generate count files using Stringtie

```

for i in {88,90,92,94,96,98}
do
/data4/laura/hisat_stranded/stringtie-1.3.3b.Linux_x86_64/stringtie SRR51943$i.rmdup.sorted.bam -p 12 -G Homo_sapiens.GRCh38.88.gtf -fr -e -o SRR51943$i.gtf
done

# Create file called gtf_list.txt in the format:
# SRR5194388 /data4/lfahey/rna_seq_ha/SRR5194388.gtf
# SRR5194390 /data4/lfahey/rna_seq_ha/SRR5194390.gtf
# SRR5194392 /data4/lfahey/rna_seq_ha/SRR5194392.gtf
# SRR5194394 /data4/lfahey/rna_seq_ha/SRR5194394.gtf
# SRR5194396 /data4/lfahey/rna_seq_ha/SRR5194396.gtf
# SRR5194398 /data4/lfahey/rna_seq_ha/SRR5194398.gtf

# Use this as input to prepDe.py

python prepDE.py -i gtf_list.txt

# Generates gene_count_matrix.csv and  transcript_count_matrix.csv

```
#### 7. Statistical analysis using EdgeR 

Run edgeR.R to output a list of differentially expressed genes as gene symbols using cutoff FDR < 0.01.

## ChIP-seq Analysis

#### 1. Fetch sra files

```

for i in {1,2,5,6};
do
/home/nextgen2015/.aspera/connect/bin/ascp -i /home/nextgen2015/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -Tr -l200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP079/SRP079176/SRR393884$i/SRR393884$i.sra ./;
done

```

#### 2. Covert to fastq files


```

fastq-dump --skip-technical --split-files --readids --read-filter pass --dumpbase --clip *.sra


```
#### 3. Carry out FASTQC and then trim based on fastqc results

```

fastqc *.fastq


```

```

for i in {41,42,45,46}
do
java -jar /home/liam/bin/trimmomatic-0.33.jar PE -phred33 SRR39388${i}_pass_1.fastq SRR39388${i}_pass_2.fastq SRR39388${i}_output_forward_paired.fq SRR39388${i}_forward_unpaired.fq SRR39388${i}_reverse_paired.fq SRR39388${i}_reverse_unpaired.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:15 MINLEN:70;
done

```

#### 4. Aignment using bowtie2

```

for i in {41,42,45,46}; do
bowtie2 -x HS_38 -1 SRR39388${i}_output_forward_paired.fq  -2 SRR39388${i}_reverse_paired.fq -S SRR39388$i.sam;
done

```

#### 5. Post-processing using samtools

```

for i in *.sam
do
id=$(echo ${i} | sed 's/.sam//')
# Convert SAM files to BAM format
samtools view -Sb ${id}.sam > ${id}.bam
# Remove possible PCR duplicates
samtools rmdup ${id}.bam ${id}.rmdup.bam
# Sort BAM file
samtools sort ${id}.rmdup.bam ${id}.rmdup.sorted
# index the bam file to creat .bai files
samtools index ${id}.rmdup.sorted.bam
# Generate file containing mapping stats
samtools flagstat ${id}.rmdup.sorted.bam > ${id}_mappingstats.txt
done

```

#### 6. Call peaks using macs:

```

macs2 callpeak -t SRR3330046.rmdup.sorted.bam SRR3330049.rmdup.sorted.bam -c SRR3330048.rmdup.sorted.bam SRR3330051.rmdup.sorted.bam -f BAMPE -g hs -n sep_reps -q 0.01 -B --tempdir /home/lfahey/

```
Edit macs .narrowpeak file output for beta:

```
sed -e '/^GL\|MT/d' macs_merged_reps_cd34-_peaks.narrowPeak | awk '{print "chr" $1, $2, $3, $4, $9}' > macs_merged_reps_cd34-_peaks.bed

```

### Integrative RNA-seq 

Integration of RNA-seq and ChIP-seq results using [BETA](http://cistrome.org/BETA/):

```

BETA_1.0.7/bin/BETA basic -p macs_xls_output.bed -e deGenes_sym.txt -k BSF -g hg19 --cutoff 1 -n basic --gname2 --output beta_ha_cd1a+


```
Cutoff was set to 1 because the regulatory potential of up-regulated or downregulated genes do not differ significantly from statically expressed genes.





