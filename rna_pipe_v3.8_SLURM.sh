#!/bin/bash
#!/usr/bin/awk
# programmer: Holden Liang， modified by Xiaoyu Zhuo

###################################################
# This is a pipline to perform initial analysis on RNA-seq data including de novo assemble and gene expression quantification, which will be ready to use for TEProf2 and DESeq2 directly, respectively.
###################################################
# how to run this script:
# bash rna_pipe_v3.8_SLURM.sh --sample rna_JHU029_NEG_BRep2 --genome hg38 --input $PWD
#####################################################################
# Update log:
# 1. pangenome library prep kit is rf. add an option for pangenome data (rf data)
# v3.8: allow single-end rna-seq data processing
#####################################################################
conda activate rna_pipe
# default settings
INPUT_DIR="." 
READ1EXTENSION="R1.fastq.gz"
STRAND="no" ## old kit (universal) we use is fr. But in most cases like pangenome data, encode data and data generated using NEB kit, it's rf
GENOME="hg38"
MAXJOBSALIGN=12 # Cores for each alignment
SAMPLE="mRNA_tumor1_AllCell"

# Arguments to bash script to decide on various parameters
while [[ $# > 1 ]]
do
key="$1"
case $key in
	-i|--input)
	INPUT_DIR="$2"
	shift # past argument
	;;
	-st|--strand)
	STRAND="$2"
	shift # past argument
	;;
	-s|--sample)
	SAMPLE="$2"
	shift # past argument
	;;
	-r|--read1extension)
	READ1EXTENSION="$2"
	shift # past argument
	;;
	-g|--genome)
	GENOME="$2"
	shift # past argument
	;;
	-m|--max)
	MAXJOBSALIGN="$2"
	shift # past argument
	;;  
	--default)
	DEFAULT=YES
	;;
	*)
	        # unknown option
	;;
esac
shift # past argument or value
done

echo "sample name: ${SAMPLE}"
echo "read 1 extension: ${READ1EXTENSION}"
echo "parallel job limit STAR: ${MAXJOBSALIGN}"
echo "fastqfile location: ${INPUT_DIR}/fastq"
echo "reference genome: ${GENOME}"
echo "library strandness: ${STRAND}"

if [[ $GENOME == "hg38" ]];then
	GENOME_DIR="/genome/STAR_v2.7.11a_index_hg38_gencodeV36/"
	GENOME_FA="/genome/GRCh38.primary_assembly.genome.fa"
elif [[ $GENOME == "mm39" ]]; then
	GENOME_DIR="/genome/STAR_index_mm39_gencodeVM34/"
	GENOME_FA="/genome/mm39.fa"
fi

if [ -d "./trimmed" ] 
then
    rm -rf trimmed aligned rseqc assembled quantification
fi

fastqdir=${INPUT_DIR}/fastq
debugdir=${INPUT_DIR}/debug
trimdir=${INPUT_DIR}/trimmed
aligndir=${INPUT_DIR}/aligned
rseqcdir=${INPUT_DIR}/rseqc
assembleddir=${INPUT_DIR}/assembled
quantificationdir=${INPUT_DIR}/quantification

mkdir ${debugdir} ${trimdir} ${aligndir} ${rseqcdir} ${assembleddir} ${quantificationdir}

	jid=`sbatch <<- VERSION | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/0_version-%j.out
	#SBATCH -e $debugdir/0_version-%j.err 
	#SBATCH -c 1
	#SBATCH --mem=1G
	#SBATCH --time=24:00:00
	#SBATCH -J "0_version_${SAMPLE}"
	#SBATCH --partition=general
	
	date

	conda activate rna_pipe

	## save softwares_version information
	echo "cutadapt version:"
	cutadapt --version
	echo "STAR version:"
	STAR --version
	echo "FastQC version:"
	fastqc --version
	echo "stringtie version:"
	stringtie --version
	echo "picard_EstimateLibraryComplexity version:"
	picard EstimateLibraryComplexity --version
	echo "samtools version:"
	samtools --version

	date
VERSION`

	jid=`sbatch <<- FASTQC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/1_fastqc-%j.out
	#SBATCH -e $debugdir/1_fastqc-%j.err
	#SBATCH -c 2
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "1_fastqc_${SAMPLE}"
	#SBATCH --partition=general

	date
	find ${fastqdir} -maxdepth 1 -name "*.gz"  | while read file ; do xbase=\\\$(basename \\\$file) ; mkdir ${fastqdir}/\\\${xbase}_fastqc; echo "fastqc -o "${fastqdir}/\\\${xbase}"_fastqc "\\\$file"" >> ${fastqdir}/1_fastqc_commands.txt; done ;
	parallel -j 2 < ${fastqdir}/1_fastqc_commands.txt
	date
FASTQC`

	jid=`sbatch <<- TRIMMING | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l 
	#SBATCH -o $debugdir/2_trimming-%j.out
	#SBATCH -e $debugdir/2_trimming-%j.err 
	#SBATCH -c 2
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "2_trimming_${SAMPLE}"
	#SBATCH --partition=general

	date

	conda activate rna_pipe

	find ${fastqdir} -maxdepth 1 -name "*${READ1EXTENSION}" | while read file ; do xbase=\\\$(basename \\\$file) ;
	echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50  -o ${trimdir}/Trimmed_"\\\$xbase" -p ${trimdir}/Trimmed_"\\\${xbase/R1/R2}" "\\\$file" "\\\${file/R1/R2}" > ${trimdir}/"\\\$xbase"_cutadapt.log" >> ${trimdir}/2_cutadaptcommands.txt ; done
	parallel -j 1 < ${trimdir}/2_cutadaptcommands.txt

	date
TRIMMING`

dependtrimming="afterok:$jid"

	jid=`sbatch <<- FASTQC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/3_fastqc-%j.out
	#SBATCH -e $debugdir/3_fastqc-%j.err 
	#SBATCH -c 2
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "3_fastqc_${SAMPLE}"
	#SBATCH --partition=general
	#SBATCH -d $dependtrimming

	date
	find ${trimdir} -maxdepth 1 -name "*.gz"  | while read file ; do xbase=\\\$(basename \\\$file) ; mkdir ${trimdir}/\\\${xbase}_fastqc; echo "fastqc -o "${trimdir}/\\\${xbase}"_fastqc "\\\$file"" >> ${trimdir}/3_fastqc_commands.txt; done ;
	parallel -j 2 < ${trimdir}/3_fastqc_commands.txt
	date
FASTQC`

	

	jid=`sbatch <<- ALIGN | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/4_alignment-%j.out
	#SBATCH -e $debugdir/4_alignment-%j.err
	#SBATCH -c $MAXJOBSALIGN
	#SBATCH --mem=128G
	#SBATCH --time=24:00:00
	#SBATCH -J "4_alignment_${SAMPLE}"
	#SBATCH -d $dependtrimming
	#SBATCH --partition=general

	date


	## this session follow the gdc guideline as i handle CPTAC RNA-seq data deposited on GDC
	## https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
	## they used gencode v36

	find ${trimdir} -name "*${READ1EXTENSION}" | while read file; do xbase=\\\$(basename \\\$file); mkdir ${aligndir}/\\\${xbase/R1.fastq.gz/1stpass_output}; echo "STAR --genomeDir ${GENOME_DIR} --readFilesIn "\\\$file" "\\\${file/R1/R2}" --runThreadN 4 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None --outFileNamePrefix ${aligndir}/\\\${xbase/R1.fastq.gz/1stpass_output}/ 2> ${aligndir}/1st_pass.err &> ${aligndir}/1st_pass.log" >> ${aligndir}/4a_alignCommands_1stpass.txt ; done ;
	echo "Aligning Reas with STAR (first pass)"
	parallel -j 1 < ${aligndir}/4a_alignCommands_1stpass.txt

	find ${trimdir} -name "*${READ1EXTENSION}" | while read file; do xbase=\\\$(basename \\\$file); echo "STAR --runMode genomeGenerate --genomeDir ${aligndir}/\\\${xbase/R1.fastq.gz/1stpass_output} --genomeFastaFiles ${GENOME_FA} --sjdbOverhang 100 --runThreadN 4 --sjdbFileChrStartEnd ${aligndir}/\\\${xbase/R1.fastq.gz/1stpass_output}/SJ.out.tab --outTmpDir ${aligndir}/\\\${xbase/R1.fastq.gz/genome_index_temp} 2> ${aligndir}/genome_indexing.err &> ${aligndir}/genome_indexing.log; mv ${INPUT_DIR}/Log.out ${aligndir}/genome_indexing.Log.out"  >> ${aligndir}/4b_alignCommands_generate_index.txt ; done ;
	echo "Aligning Reas with STAR (generating index)"
	parallel -j 1 < ${aligndir}/4b_alignCommands_generate_index.txt 

	find ${trimdir} -name "*${READ1EXTENSION}" | while read file; do xbase=\\\$(basename \\\$file); echo "STAR --runMode alignReads --genomeDir ${aligndir}/\\\${xbase/R1.fastq.gz/1stpass_output} --readFilesIn "\\\$file" "\\\${file/R1/R2}" --runThreadN 4 --outFileNamePrefix ${aligndir}/"\\\${xbase%.*}" --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 2> ${aligndir}/2nd_pass.err &> ${aligndir}/2nd_pass.log" >> ${aligndir}/4c_alignCommands_2ndpass.txt;done;
	echo "Aligning Reas with STAR (second pass)"
	parallel -j 1 < ${aligndir}/4c_alignCommands_2ndpass.txt

	date
ALIGN`

dependalignment="afterok:$jid"


	jid=`sbatch <<- INDEX | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/5_index-%j.out
	#SBATCH -e $debugdir/5_index-%j.err 
	#SBATCH -c 1
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "5_index_${SAMPLE}"
	#SBATCH --partition=general
	#SBATCH -d $dependalignment    

	date
	find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "samtools index \\\$file ; samtools idxstats \\\$file > \\\${file}_idxstats" >> ${aligndir}/5_indexingCommands.txt; done;
	parallel -j 1 < ${aligndir}/5_indexingCommands.txt
	date
INDEX`

dependindex="afterok:$jid"

	jid=`sbatch <<- BIGWIG | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/6_bigwig-%j.out
	#SBATCH -e $debugdir/6_bigwig-%j.err 
	#SBATCH -c 4
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "6_bigwig_${SAMPLE}"
	#SBATCH --partition=general
	#SBATCH -d $dependindex

	date

	conda activate rna_pipe

	if [[ $STRAND == "no" ]]; then
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/forward.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand forward 2> $debugdir/6_bam2bigwigCommand_forward.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/reverse.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand reverse --scaleFactor -1 2> $debugdir/6_bam2bigwigCommand_reverse.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/.bigwig} -of bigwig -p 3 --normalizeUsing RPKM 2> $debugdir/bam2bigwigCommand.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
	elif [[ $STRAND == "fr" ]]; then
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/reverse.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand forward --scaleFactor -1 2> $debugdir/6_bam2bigwigCommand_forward.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/forward.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand reverse 2> $debugdir/6_bam2bigwigCommand_reverse.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
	elif [[ $STRAND == "rf" ]] ; then
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/forward.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand forward 2> $debugdir/6_bam2bigwigCommand_forward.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file; do xbase=\\\$(basename \\\$file); echo "bamCoverage -b \\\$file -o \\\${file/bam/reverse.bigwig} -of bigwig -p 3 --normalizeUsing RPKM --filterRNAstrand reverse --scaleFactor -1 2> $debugdir6_/bam2bigwigCommand_reverse.log" >> ${aligndir}/6_bam2bigwigCommand.txt; done;
	fi

 	parallel -j 3 < ${aligndir}/6_bam2bigwigCommand.txt

	date
BIGWIG`


	jid=`sbatch <<- RSEQC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/7_rseqc-%j.out
	#SBATCH -e $debugdir/7_rseqc-%j.err 
	#SBATCH -c 3
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "7_rseqc_${SAMPLE}"
	#SBATCH --partition=general
	#SBATCH -d $dependindex

	date

	conda activate rna_pipe

	find ${aligndir} -name "*sortedByCoord.out.bam" | while read file ; do xbase=\\\$(basename \\\$file) ; \
	echo "geneBody_coverage.py -i "\\\$file" -o ${rseqcdir}/"\\\${xbase%.*}" -r /genome/rseqc/hg38.HouseKeepingGenes.bed" >> ${rseqcdir}/7_rseqQCcommands.txt ; \
	echo "read_distribution.py -i "\\\$file" -r /genome/rseqc/hg38_Gencode_V28.bed > ${rseqcdir}/"\\\${xbase%.*}".readdistribution.txt" >> ${rseqcdir}/7_rseqQCcommands.txt ; \
	echo "junction_saturation.py -i "\\\$file" -o ${rseqcdir}/"\\\${xbase%.*}" -r /genome/rseqc/hg38_Gencode_V28.bed" >> ${rseqcdir}/7_rseqQCcommands.txt ; \
	echo "picard EstimateLibraryComplexity I="\\\$file" O=${rseqcdir}/"\\\${xbase%.*}"_duplication_stats.txt; mv ${INPUT_DIR}/log.txt ${rseqcdir}/log.txt" >> ${rseqcdir}/7_rseqQCcommands.txt ; done ;
	parallel -j 4 < ${rseqcdir}/7_rseqQCcommands.txt

	date
RSEQC`


	jid=`sbatch <<- ASSEMBLE | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/8_assemble-%j.out
	#SBATCH -e $debugdir/8_assemble-%j.err 
	#SBATCH -c 1
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "8_assemble_${SAMPLE}"
	#SBATCH --partition=general
	#SBATCH -d $dependindex

	date

	if [[ $STRAND == "no" ]]; then
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "samtools view -q 255 -h "\\\$file" | stringtie - -o ${assembleddir}/"\\\${xbase%.*}".gtf -p 4 -m 100 -c 1" >> ${assembleddir}/8_assembleCommands.txt ; done ;
	elif [[ $STRAND == "fr" ]]; then
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "samtools view -q 255 -h "\\\$file" | stringtie - -o ${assembleddir}/"\\\${xbase%.*}".gtf -p 4 -m 100 -c 1 --fr" >> ${assembleddir}/8_assembleCommands.txt ; done ;
	elif [[ $STRAND == "rf" ]] ; then
		find ${aligndir} -name "*sortedByCoord.out.bam" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "samtools view -q 255 -h "\\\$file" | stringtie - -o ${assembleddir}/"\\\${xbase%.*}".gtf -p 4 -m 100 -c 1 --rf" >> ${assembleddir}/8_assembleCommands.txt ; done ;
	fi

	parallel -j 1 < ${assembleddir}/8_assembleCommands.txt

	date
ASSEMBLE`


	jid=`sbatch <<- QUANTIFICATION | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -o $debugdir/9_quantification-%j.out
	#SBATCH -e $debugdir/9_quantification-%j.err 
	#SBATCH -c 1
	#SBATCH --mem=20G
	#SBATCH --time=24:00:00
	#SBATCH -J "9_quantification_${SAMPLE}"
	#SBATCH --partition=general
	#SBATCH -d $dependindex

	date

	find ${aligndir} -name "*.sortedByCoord.out.bam" | while read file ; do xbase=\\\$(basename \\\$file) ; echo "samtools view -q 255 -h "\\\$file" | stringtie - -o ${quantificationdir}/"\\\${xbase%.*}".gtf -e -b ${quantificationdir}/"\\\${xbase%.*}".stats -p 4 -m 100 -c 1 -G /genome/gencode.v36.primary_assembly.annotation.gtf" >> ${quantificationdir}/9_quantification_commands.txt ; done ;
	parallel -j 1 < ${quantificationdir}/9_quantification_commands.txt

	date
QUANTIFICATION`





