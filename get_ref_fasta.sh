#!/bin/bash
ref=${1}
if [ $ref = "dmel" ]; then echo $PI_SCRATCH/REFERENCE/d_mel/fasta_all/all_dmel.fasta; else
	if [ $ref = "dmelmicro" ]; then echo $PI_SCRATCH/REFERENCE/d_mel/fasta_all/dmel_and_bacteria.fasta; else
		if [ $ref = "human" ]; then echo $PI_SCRATCH/REFERENCE/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna; else
			if [ $ref = "dsim" ]; then echo $PI_SCRATCH/REFERENCE/d_sim/Drosophila_simulans.ASM75419v3.dna.chromosome.3R.fa; else
				echo "ref not known"; fi
		fi
	fi
fi