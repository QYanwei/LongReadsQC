# LongReadsQC
Quality control tool for long-read sequencing data




Input sequencing files type:
---
1. **For FASTQ files**
  Usage: fastq_read_parse.py [-h] --input INPUT [--stat STAT] [--filtering FILTERING] [--qualified_qscore_threshold QUALIFIED_QSCORE_THRESHOLD]
                           [--qualified_length_threshold QUALIFIED_LENGTH_THRESHOLD] [--cut_read_ends CUT_READ_ENDS] [--cut_head_length CUT_HEAD_LENGTH] [--cut_tail_length CUT_TAIL_LENGTH]
                           [--trim_ends_homopolybase TRIM_ENDS_HOMOPOLYBASE] [--trim_head_homopolybase TRIM_HEAD_HOMOPOLYBASE] [--trim_tail_homopolybase TRIM_TAIL_HOMOPOLYBASE]
                           [--downsampling DOWNSAMPLING] [--pickup_number_of_reads PICKUP_NUMBER_OF_READS] [--pickup_over_length_of_reads PICKUP_OVER_LENGTH_OF_READS]
                           [--pickup_less_length_of_reads PICKUP_LESS_LENGTH_OF_READS] [--downsample_genomesize DOWNSAMPLE_GENOMESIZE] [--downsample_depthx DOWNSAMPLE_DEPTHX] [--kmer KMER]
                           [--kmer_size KMER_SIZE] [--homopolymer HOMOPOLYMER] [--homopolymer_length_min HOMOPOLYMER_LENGTH_MIN] [--homopolymer_length_max HOMOPOLYMER_LENGTH_MAX] --output OUTPUT
                           [--txt TXT] [--json JSON]
  Help Doc for each fastq_read_parse parameters:
  ```
    optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         raw data: fastq/fq.gz
  --stat STAT           statistics the quality and length of reads in fastq file.
  --filtering FILTERING
                        filter short or low-quality read or not: 1 or 0.
  --qualified_qscore_threshold QUALIFIED_QSCORE_THRESHOLD
                        filtered quality of read (default: 12).
  --qualified_length_threshold QUALIFIED_LENGTH_THRESHOLD
                        filtered quality of read (default: 1000).
  --cut_read_ends CUT_READ_ENDS
                        cut the head or tail segment of read or not: 1 or 0.
  --cut_head_length CUT_HEAD_LENGTH
                        cut the head segment of read (default: 100) under --cut_double_ends_length with 1, only cut head segment setting --cut_tail_length to 0.
  --cut_tail_length CUT_TAIL_LENGTH
                        cut the head segment of read (default: 100) under --cut_double_ends_length with 1, only cut tail segment setting --cut_head_length to 0.
  --trim_ends_homopolybase TRIM_ENDS_HOMOPOLYBASE
                        trim the long homopoly bases segment of read or not: 1 or 0.
  --trim_head_homopolybase TRIM_HEAD_HOMOPOLYBASE
                        trim the head homopolybase of read (default: 100) under --trim_ends_homopolybase with 1, only cut head segment setting --trim_tail_homopolybase to 0.
  --trim_tail_homopolybase TRIM_TAIL_HOMOPOLYBASE
                        trim the head homopolybase of read (default: 100) under --trim_ends_homopolybase with 1, only cut tail segment setting --trim_head_homopolybase to 0.
  --downsampling DOWNSAMPLING
                        downsampling or not: 1 or 0.
  --pickup_number_of_reads PICKUP_NUMBER_OF_READS
                        pick up the number of read, parameter under --downsampling with 1.
  --pickup_over_length_of_reads PICKUP_OVER_LENGTH_OF_READS
                        pick up the minimum length of read, parameter under --downsampling with 1.
  --pickup_less_length_of_reads PICKUP_LESS_LENGTH_OF_READS
                        pick up the maximum length of read, parameter under --downsampling with 1.
  --downsample_genomesize DOWNSAMPLE_GENOMESIZE
                        the genome size of downsample data followed by sample species, binding with --downsample_depthx under --downsample with 1.
  --downsample_depthx DOWNSAMPLE_DEPTHX
                        the depth (x) fold of downsample data custom-made by user, binding with --downsample_genomesize under --downsample with 1.
  --kmer KMER           homopolymer analysis or not: 1 or 0.
  --kmer_size KMER_SIZE
                        observed the minimum length of kmer(default: 6).
  --homopolymer HOMOPOLYMER
                        homopolymer analysis or not: 1 or 0.
  --homopolymer_length_min HOMOPOLYMER_LENGTH_MIN
                        observed the minimum length of homopolymer(default: 5).
  --homopolymer_length_max HOMOPOLYMER_LENGTH_MAX
                        observed the maximum length of homopolymer(default: 50).
  --output OUTPUT       the sign of output file name, following by sample name is better.
  --txt TXT             output table file name, filetype=TXT.
  --json JSON           output table file name, filetype=JSON.

  ```
2. **For FASTQ files**
usage: fasta_read_parse.py [-h] --input INPUT [--stat STAT] [--filtering FILTERING] [--qualified_length_threshold QUALIFIED_LENGTH_THRESHOLD] [--cut_read_ends CUT_READ_ENDS]
                           [--cut_head_length CUT_HEAD_LENGTH] [--cut_tail_length CUT_TAIL_LENGTH] [--trim_ends_homopolybase TRIM_ENDS_HOMOPOLYBASE]
                           [--trim_head_homopolybase TRIM_HEAD_HOMOPOLYBASE] [--trim_tail_homopolybase TRIM_TAIL_HOMOPOLYBASE] [--downsampling DOWNSAMPLING]
                           [--pickup_number_of_reads PICKUP_NUMBER_OF_READS] [--pickup_over_length_of_reads PICKUP_OVER_LENGTH_OF_READS] [--pickup_less_length_of_reads PICKUP_LESS_LENGTH_OF_READS]
                           [--downsample_genomesize DOWNSAMPLE_GENOMESIZE] [--downsample_depthx DOWNSAMPLE_DEPTHX] [--kmer KMER] [--kmer_size KMER_SIZE] [--homopolymer HOMOPOLYMER]
                           [--homopolymer_length_min HOMOPOLYMER_LENGTH_MIN] [--homopolymer_length_max HOMOPOLYMER_LENGTH_MAX] --output OUTPUT [--xlsx XLSX] [--json JSON]

Help Doc for each fastq_read_parse parameters:
  ```
    optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         raw data: fasta/fa.gz
  --stat STAT           statistic the length of reads in fasta file.
  --filtering FILTERING
                        filter short or low-quality read or not: 1 or 0.
  --qualified_length_threshold QUALIFIED_LENGTH_THRESHOLD
                        filtered quality of read (default: 1000).
  --cut_read_ends CUT_READ_ENDS
                        cut the head or tail segment of read or not: 1 or 0.
  --cut_head_length CUT_HEAD_LENGTH
                        cut the head segment of read (default: 100) under --cut_double_ends_length with 1, only cut head segment setting --cut_tail_length to 0.
  --cut_tail_length CUT_TAIL_LENGTH
                        cut the head segment of read (default: 100) under --cut_double_ends_length with 1, only cut tail segment setting --cut_head_length to 0.
  --trim_ends_homopolybase TRIM_ENDS_HOMOPOLYBASE
                        trim the long homopoly bases segment of read or not: 1 or 0.
  --trim_head_homopolybase TRIM_HEAD_HOMOPOLYBASE
                        trim the head homopolybase of read (default: 100) under --trim_ends_homopolybase with 1, only cut head segment setting --trim_tail_homopolybase to 0.
  --trim_tail_homopolybase TRIM_TAIL_HOMOPOLYBASE
                        trim the head homopolybase of read (default: 100) under --trim_ends_homopolybase with 1, only cut tail segment setting --trim_head_homopolybase to 0.
  --downsampling DOWNSAMPLING
                        downsampling or not: 1 or 0.
  --pickup_number_of_reads PICKUP_NUMBER_OF_READS
                        pick up the number of read, parameter under --downsampling with 1.
  --pickup_over_length_of_reads PICKUP_OVER_LENGTH_OF_READS
                        pick up the minimum length of read, parameter under --downsampling with 1.
  --pickup_less_length_of_reads PICKUP_LESS_LENGTH_OF_READS
                        pick up the maximum length of read, parameter under --downsampling with 1.
  --downsample_genomesize DOWNSAMPLE_GENOMESIZE
                        the genome size of downsample data followed by sample species, binding with --downsample_depthx under --downsample with 1.
  --downsample_depthx DOWNSAMPLE_DEPTHX
                        the depth (x) fold of downsample data custom-made by user, binding with --downsample_genomesize under --downsample with 1.
  --kmer KMER           homopolymer analysis or not: 1 or 0.
  --kmer_size KMER_SIZE
                        observed the minimum length of kmer(default: 6).
  --homopolymer HOMOPOLYMER
                        homopolymer analysis or not: 1 or 0.
  --homopolymer_length_min HOMOPOLYMER_LENGTH_MIN
                        observed the minimum length of homopolymer(default: 5).
  --homopolymer_length_max HOMOPOLYMER_LENGTH_MAX
                        observed the maximum length of homopolymer(default: 50).
  --output OUTPUT       the sign of output file name, following by sample name is better.
  --xlsx XLSX           output table file name, filetype=XLSX
  --json JSON           output table file name, filetype=JSON

  ```
