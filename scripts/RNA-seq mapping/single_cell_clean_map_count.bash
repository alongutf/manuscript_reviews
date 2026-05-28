#!/bin/bash




export PATH=/mnt/Spinning1/Adir/cellranger/cellranger-7.0.0:$PATH

# create the genome refference:
cellranger mkref --genome=Ecoli_with_extras_Adi_Rotem --fasta=Ecoli_with_extras_Adi_Rotem.FASTA --genes=Ecoli_with_extras_Adi_Rotem_corrected2.gtf	

samples=(
NB_Adi_5-15a_CTAGGTGA_S5
NB_Adi_6-15b_CGCTATGT_S6
NB_Adi_4-13b_TATGATTC_S4
NB_Adi_3-13a_CAGTACTG_S3
NB_Adi_1-2a_GGTTTACT_S1
NB_Adi_2-2b_TTTCATGA_S2
)

# This should be the same as samples, just remove the _s# from the end
samples_ids_for_cell_ranger=(
NB_Adi_5-15a_CTAGGTGA
NB_Adi_6-15b_CGCTATGT
NB_Adi_4-13b_TATGATTC
NB_Adi_3-13a_CAGTACTG
NB_Adi_1-2a_GGTTTACT
NB_Adi_2-2b_TTTCATGA
)

i=0
for filename in "${samples[@]}";
	do 
	#filename="PL_ROSEN_9070_ACAGAGGT_S7";
	echo $filename;
	
	
	echo "${samples_ids_for_cell_ranger[$i]}"
: <<'END_COMMENT'	
	# replace 10x barcode with the PCR UMI from read2

	
	for f in "$filename"*_R2_*.fastq.gz;
		do
		r1str="${f/_R2_/_R1_}";
		r2str=$f;

		# last 12 bp are 10x molecule barcode that we need to replace. keep the first 16
		# -j 0: automaticlly detect availble cores. 
		~/.local/bin/cutadapt -l 16 -j 0 -o "R1cellUMI_${r1str%.fastq.gz}.fastq" $r1str;

		# make a fastq file with the probe UMIs
		# (in most reads...) first 100 chareters are illumina read1 P5 and the read data. After that there is a 12bp long probe UMI we need
		# ~/.local/bin/cutadapt -u 100 -l 12 -j 0 -o "R2probeUMI_${r2str%.fastq.gz}.fastq" $r2str;
		
		# The 12bp UMI (inserted by Adam) is right after the "Extender" sequence. We use cut adapt -g option to get a fasta 
		# with the 12bp following the extender:
		~/.local/bin/cutadapt -l 12 -j 0 -g CATAGTTTCTTCGAGCAAGCTT -o "R2probeUMI_${r2str%.fastq.gz}.fastq" $r2str;
		
		
		# concatenate the probe UMIs and the cell UMIs to create read1
		# (keeps the header of the first file):
		paste -d '\n' "R1cellUMI_${r1str%.fastq.gz}.fastq" "R2probeUMI_${r2str%.fastq.gz}.fastq"  | sed -n 'p;n;n;N;s/\n//p' > "probeUMI_${r1str%.fastq.gz}.fastq";

		# rename read2
		cp $r2str "probeUMI_${r2str%.fastq.gz}.fastq";
		r1str="probeUMI_${r1str%.fastq.gz}.fastq";
		r2str="probeUMI_${r2str%.fastq.gz}.fastq";
		echo $r1str;
		echo $r2str;
		
		# remove polyA tail (do not need any tail trimming if you trim after the extender)
		#~/.local/bin/cutadapt  -m 6 -j 0 -A "A{100}" -o "trimmed1_${r1str}" -p "trimmed1_${r2str}" $r1str $r2str;
		# remove polyG tail
		#~/.local/bin/cutadapt  -m 6 -j 0 -A "G{100}" -o "trimmed2_${r1str}" -p "trimmed2_${r2str}" "trimmed1_${r1str}" "trimmed1_${r2str}";
		# remove tail with low quality
		#~/.local/bin/cutadapt -q 20 -m 6 -j 0 -o "trimmed3_${r1str}" -p "trimmed3_${r2str}" "trimmed2_${r1str}" "trimmed2_${r2str}";
		
		
		# remove PCR 3' handle (every read2 has this)
		~/.local/bin/cutadapt -G GATGACCCGGTCCATACA -m 6 -j 0 -o "trimmed_${r1str}" -p "trimmed_${r2str}" "${r1str}" "${r2str}"; 
		
		
		# remove extender 5' and anything biond it (every read2 has this)
		~/.local/bin/cutadapt -a CATAGTTTCTTCGAGCAAGCTT -m 6 -j 0 -o "trimmed2_${r1str}" -p "trimmed2_${r2str}" "trimmed_${r1str}" "trimmed_${r2str}"; 
		
		cp "trimmed2_${r1str}" "trimmed_${r1str}";
		cp "trimmed2_${r2str}" "trimmed_${r2str}";
		
		# up to ~54-115 bp in the R2 there is a region of low quality, PCR adaptors,  UMI2 , polyA
		~/.local/bin/cutadapt -l 54 -j 0 -o "trimmed2_${r2str}" "trimmed_${r2str}"
		
		#cp "trimmed2_${r1str}" "trimmed_${r1str}";
		cp "trimmed2_${r2str}" "trimmed_${r2str}";
		
		# another cutadapt run just make sure the paired-end files R1 and R2 are compatible.
		# Only keep reads with 28bp. If read1 has less then that, its not compatible with our chemestry (V3)
		~/.local/bin/cutadapt -m 28 -j 0 -o "trimmed2_${r1str}" -p "trimmed2_${r2str}" "trimmed_${r1str}" "trimmed_${r2str}";
		
		cp "trimmed2_${r1str}" "trimmed_${r1str}";
		cp "trimmed2_${r2str}" "trimmed_${r2str}";
		
	done;


	# Remove temp files
	#rm trimmed5_*
	#rm trimmed4_*
	#rm trimmed3_*
	rm trimmed2_*
	#rm trimmed1_*


END_COMMENT

	# Run cellranger mapping:
	
	#cellranger count --id="trimmed6_probeUMI_${filename}" --r1-length 26 --transcriptome=/mnt/Spinning1/Adir/scRNAseq_CASP_SHX_bioreps_220814/Ecoli_with_extras2 --fastqs=/mnt/Spinning1/Adir/scRNAseq_CASP_SHX_bioreps_220814/ --sample="trimmed6_probeUMI_${filename:0:22}"
	#cellranger count --id="trimmed6_genome2_force10k_probeUMI_${filename}" --r1-length 26 --transcriptome=/mnt/Spinning1/Adir/scRNAseq_CASP_SHX_bioreps_220814/Ecoli_with_extras2 --force-cells 10000 --fastqs=/mnt/Spinning1/Adir/scRNAseq_CASP_SHX_bioreps_220814/ --sample="trimmed6_probeUMI_${filename:0:22}"

	cellranger count --id="trimmed_probeUMI_${filename}" --chemistry=SC3Pv3 --transcriptome=/mnt/Spinning1/Adir/scRNAseq_3samples_230305/Ecoli_with_extras_Adi_Rotem --fastqs=/mnt/Spinning1/Adir/scRNAseq_3samples_230305 --sample="trimmed_probeUMI_${samples_ids_for_cell_ranger[$i]}"
	

	# to get a UMI count of each probe: 
	cd "trimmed_probeUMI_${filename}/outs/filtered_feature_bc_matrix/"
	# Print line number along with contents of barcodes.tsv.gz and genes.tsv.gz 

	zcat barcodes.tsv.gz | awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1}' | sort -t, -k 1b,1 > numbered_barcodes.csv
	zcat features.tsv.gz | awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1,$2,$3}' | sort -t, -k 1b,1 > numbered_features.csv

	# Skip the header lines and sort matrix.mtx.gz
	zcat matrix.mtx.gz | tail -n +4 | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 1b,1 > feature_sorted_matrix.csv
	zcat matrix.mtx.gz | tail -n +4 | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 2b,2 > barcode_sorted_matrix.csv

	# Use join to replace line number with barcodes and genes
	join -t, -1 1 -2 1 numbered_features.csv feature_sorted_matrix.csv | cut -d, -f 2,3,4,5,6 | sort -t, -k 4b,4 | join -t, -1 1 -2 4 numbered_barcodes.csv - | cut -d, -f 2,3,4,5,6 > "counts_trimmed_probeUMI_${filename}.csv"

	# add column names to the csv (for more info on the final matrix see https://kb.10xgenomics.com/hc/en-us/articles/360023793031-How-can-I-convert-the-feature-barcode-matrix-from-Cell-Ranger-3-to-a-CSV-file-)
	sed -i 1i"10x Genomics cellular barcode,Feature ID,Feature name,Feature type,UMI count" "counts_trimmed_probeUMI_${filename}.csv"

	# Remove temp files
	rm -f barcode_sorted_matrix.csv feature_sorted_matrix.csv numbered_barcodes.csv numbered_features.csv


	echo "the counts are in (name sample)/outs/filtered_feature_bc_matrix/counts(name sample).csv".

	cp "counts_trimmed_probeUMI_${filename}.csv" ../../../
	
	cd ../
	
	cp web_summary.html "../../web_summary_${samples[$i]}.html"
	
	cd ../../
	i=$(($i + 1));
	
done;
