kraken2 \
--db /home/raphenar/sars-cov-2/round1/human/Iran1-LA/kraken/db \
--threads 20 --quick \
--unclassified-out unclassified-sequences# \
--classified-out classified-sequences# \
--output output \
--paired \
--gzip-compressed \
*paired.fastq.gz \
--report report
