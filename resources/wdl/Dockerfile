# TAG dnastack/signal:1
FROM finlaymaguire/signal

MAINTAINER Heather Ward <heather@dnastack.com>

RUN git pull

# This is not optimal, but it works & stays up to date
RUN cp example_config.yaml config.yaml && \
	cp example_sample_table.csv sample_table.csv && \
	mkdir data && \
	touch data/MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz \
		data/MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz \
		data/MN908947_3.fasta && \
	conda run \
		-n snakemake \
		snakemake \
		--verbose \
		--use-conda \
		--create-envs-only \
		--cores 12 \
		--conda-prefix $HOME/.snakemake \
		-s Snakefile all && \
	rm -rf data config.yaml sample_table.csv
