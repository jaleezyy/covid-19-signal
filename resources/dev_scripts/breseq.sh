# Reference MN908947.3

mkdir Iran1-LA_MN
breseq --reference /workspace/sars-cov-2/ref_seq/MN908947_3.gbk \
--num-processors 80 \
--polymorphism-prediction \
--brief-html-output \
--output Iran1-LA_MN \
Iran1-LA-R1_paired.fastq.gz \
Iran1-LA-R2_paired.fastq.gz 
