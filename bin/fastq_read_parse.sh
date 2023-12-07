python fastq_read_parse.py --input ./data/wt_rta.fq.gz --trim_ends_homopolybase 1 --trim_head_homopolybase 4 --trim_tail_homopolybase 4  --output ./result/fastq_output 
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --cut_read_ends 1 --cut_head_length 100 --cut_tail_length 100 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --filtering 1 --qualified_qscore_threshold 10 --qualified_length_threshold 200 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --downsampling 1 --pickup_number_of_reads 100 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --downsampling 1 --pickup_over_length_of_reads 1000 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --downsampling 1 --pickup_less_length_of_reads 1000 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --downsampling 1 --downsample_genomesize 100000 --downsample_depthx 10 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --kmer 1 --kmer_size 16  --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --homopolymer 1 --homopolymer_length_min 7 --homopolymer_length_max 100 --output ./result/fastq_output
python fastq_read_parse.py --input ./data/wt_rta.fq.gz --trim_ends_homopolybase 1 --trim_head_homopolybase 4 --trim_tail_homopolybase 4 --cut_read_ends 1 --cut_head_length 100 --cut_tail_length 100 --filtering 1 --qualified_qscore_threshold 10 --qualified_length_threshold 200  --downsampling 1 --pickup_number_of_reads 100 --pickup_over_length_of_reads 1000 --pickup_less_length_of_reads 1000  --downsample_genomesize 100000 --downsample_depthx 10 --kmer 1 --kmer_size 16 --homopolymer 1 --homopolymer_length_min 7 --homopolymer_length_max 100 --output ./result/fastq_output > result/fastq_output_LR_windows_freq.txt
