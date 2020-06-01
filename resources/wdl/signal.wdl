version 1.0

workflow signal {
	input {
		String samplename
		Array [File] fastq_R1s
		Array [File] fastq_R2s

		File viral_reference_genome
		File human_reference

		## Trimgalore, Ivar
		Int min_length = 20
		Int min_qual = 20

		## Ivar
		File scheme_bed
		Float ivar_freq_threshold = 0.75
		Int ivar_min_coverage_depth = 10
		Float ivar_min_freq_threshold = 0.25
		Int ivar_min_variant_quality = 20

		## Samtools mpileup
		Int mpileup_depth = 100000

		## breseq
		File breseq_reference

		## kraken
		File kraken_db

		## quast
		File viral_reference_feature_coords
	}

	call run_signal {
		input:
			samplename = samplename,
			fastq_R1s = fastq_R1s,
			fastq_R2s = fastq_R2s,
			min_qual = min_qual,
			min_len = min_length,
			scheme_bed = scheme_bed,
			human_reference = human_reference,
			viral_reference_genome = viral_reference_genome,
			viral_reference_feature_coords = viral_reference_feature_coords,
			breseq_reference = breseq_reference,
			kraken2_db = kraken_db,
			mpileup_depth = mpileup_depth,
			ivar_freq_threshold = ivar_freq_threshold,
			ivar_min_coverage_depth = ivar_min_coverage_depth,
			ivar_min_freq_threshold = ivar_min_freq_threshold,
			ivar_min_variant_quality = ivar_min_variant_quality		
	}

	output {
		File sample_output = run_signal.sample_output
	}
}

task run_signal {
	input {
		String samplename
		Array [File] fastq_R1s
		Array [File] fastq_R2s
		Int min_qual
		Int min_len
		File scheme_bed
		File human_reference
		File viral_reference_genome
		File viral_reference_feature_coords
		File breseq_reference
		File kraken2_db
		Int mpileup_depth
		Float ivar_freq_threshold
		Int ivar_min_coverage_depth
		Float ivar_min_freq_threshold
		Int ivar_min_variant_quality
	}

	String kraken2_db_base = basename(kraken2_db, ".tar.gz") + "/db"
	String human_reference_base = basename(human_reference, ".tar.gz")
	Int num_paired_fastqs = length(fastq_R1s)
	Int threads = 8

	command {
		# generate config
		echo \
		"min_qual: ~{min_qual}
		min_len: ~{min_len}
		scheme_bed: '~{scheme_bed}'
		human_reference: '$(pwd)/~{human_reference_base}'
		viral_reference_genome: '~{viral_reference_genome}'
		viral_reference_feature_coords: '~{viral_reference_feature_coords}'
		breseq_reference: '~{breseq_reference}'
		kraken2_db: '$(pwd)/~{kraken2_db_base}'
		mpileup_depth: ~{mpileup_depth}
		ivar_freq_threshold: ~{ivar_freq_threshold}
		ivar_min_coverage_depth: ~{ivar_min_coverage_depth}
		ivar_min_freq_threshold: ~{ivar_min_freq_threshold}
		ivar_min_variant_quality: ~{ivar_min_variant_quality}
		samples: '$(pwd)/sample_table.csv'" | tr -d '\t' > config.yaml

		tar -zxvf ~{kraken2_db}
		tar -zxvf ~{human_reference}

		# generate sample_table.csv
		yes ~{samplename} | head -~{num_paired_fastqs} > samplename.tmp
		echo -e "~{sep='\n' fastq_R1s}" >> fastq_R1s.tmp
		echo -e "~{sep='\n' fastq_R2s}" >> fastq_R2s.tmp
		echo sample,r1_path,r2_path > sample_table.csv
		paste -d , samplename.tmp fastq_R1s.tmp fastq_R2s.tmp >> sample_table.csv && rm samplename.tmp fastq_R1s.tmp fastq_R2s.tmp

		conda run \
			-n snakemake \
			snakemake \
			--verbose \
			--use-conda \
			--conda-prefix $HOME/.snakemake \
			--cores ~{threads} \
			-s /covid-19-signal/Snakefile all

		tar -zcvf ~{samplename}.tar.gz ~{samplename}
	}

	output {
		File sample_output = "~{samplename}.tar.gz"
	}

	runtime {
		docker: "dnastack/signal:1"
		cpu: 8
		memory: "32 GB"
		disks: "local-disk 500 HDD"
		bootDiskSizeGb: 20
		preemptible: 2
	}
}
