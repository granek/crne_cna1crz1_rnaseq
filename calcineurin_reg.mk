# make all  -l 16 -j # for parallel make
# make --load-average 12 --jobs 12

# need bash to allow for process substitution for fastq-mcf: <(gunzip ...)
# might be able to use bash -c "COMMANDLINE" for that line specifically
SHELL := /bin/bash

.PHONY: clean test results
# .PRECIOUS: %.bam
# .SECONDARY: %.bam
.PRECIOUS:
.SECONDARY:
#--------------------------------------------------
GENOME_DIR=genome
INDEXBASE=$(GENOME_DIR)/h99_2
SEQ_URL="http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EASupercontigsFasta&sp=SCNA2&sp=S.gz"
GTF_URL="http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EATranscriptsGtf&sp=SCNA2&sp=S.gz"
SEQ= $(INDEXBASE).fa
GTF= $(INDEXBASE).gtf
# INDEXBASE=$(basename $(GENOME))
BT2_INDEX_FILES=$(addprefix $(INDEXBASE),.1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2)
#--------------------------------------------------
RAW_FASTQ_DIR=/nfs/gems_sata/heitman/calcineurin_reg/rnaseq/raw_fastqs
FINAL_FASTQ_DIR=final_fastqs
#--------------------------------------------------
RAW_FASTQS_FULLPATH := $(wildcard $(RAW_FASTQ_DIR)/*.fastq.gz)
NUMTHREADS=6
# RAW_FASTQS_FULLPATH := $(wildcard $(RAW_FASTQ_DIR)/SC-ECRNA14_*.fastq.gz)
# NUMTHREADS=10
#--------------------------------------------------
FASTQS := $(notdir $(RAW_FASTQS_FULLPATH))
FINAL_FASTQS := $(addprefix $(FINAL_FASTQ_DIR)/,$(FASTQS))
#--------------------------------------------------
TOPHAT_BASE_DIR=thout
FASTQ_SUFFIX=_R1_001.fastq.gz
# FINAL_BAMS := $(patsubst %$(FASTQ_SUFFIX),$(TOPHAT_BASE_DIR)/%/accepted_hits.bam,$(FASTQS))
FINAL_BAMS := $(patsubst %$(FASTQ_SUFFIX),$(TOPHAT_BASE_DIR)/%.bam,$(FASTQS))
FINAL_BAIS := $(addsuffix .bai,$(FINAL_BAMS))
#--------------------------------------------------
#--------------------------------------------------
COUNT_DIR=counts
READ_COUNTS := $(patsubst %$(FASTQ_SUFFIX),$(COUNT_DIR)/%_counts.tab,$(FASTQS))

INFO_DIR=info
CNEO_ANNOT=$(INFO_DIR)/Cneo_H99.AHRD.20131001.tab

TMP_DIR=tmp
FIG_DIR=figures

PILEUP_PLOTS = $(FIG_DIR)/CNAG_00156_pileup.pdf $(FIG_DIR)/CNAG_04796_pileup.pdf
PILEUP_PLOT_TABLE = $(INFO_DIR)/plot_pileup_table.csv

ILLUMINA_ADAPTERS=$(INFO_DIR)/illumina_adapters.fasta
ADAPTERS_USED=$(FINAL_FASTQ_DIR)/adapter.fa
INDICES=$(TMP_DIR)/adapter_indices.txt
VERSION_FILE=version.txt

DESEQ_RESULTS_DIR:=results

#--------------------------------------------------
TOPHAT=tophat
SAMTOOLS=samtools
#--------------------------------------------------
dir_guard=@mkdir -p $(@D)

all:  results todo version

results : $(READ_COUNTS) $(FINAL_BAIS) $(FINAL_BAMS)

under_development : $(READ_COUNTS) $(FINAL_BAIS) $(FINAL_BAMS) # $(BT2_INDEX_FILES) $(FINAL_FASTQS) 

todo:
	@echo "GET ANNOTATIONS FROM GUILHEM(sp?)"

test: $(TOPHAT_BASE_DIR)/SC-ECRNA10_TAGCTT_L003/test.bam
	@echo $(FINAL_BAMS)
	@echo $(FINAL_BAIS)
	@echo $(READ_COUNTS)

figs : $(PILEUP_PLOTS)

bams :  $(FINAL_BAMS)

version : $(VERSION_FILE)

$(VERSION_FILE) :
	$(TOPHAT) -v > $@
	$(SAMTOOLS) -v >> $@
	$(RSCRIPT) --version >> $@

#===============================================================================
#===============================================================================
# Download and merge reference genomes 
refgenome : $(SEQ) $(GTF)

$(SEQ) :
	$(dir_guard)
	# The following changes the chromosome names in the genome sequence file so they are compatible with the gtf file
	curl $(SEQ_URL) | zcat | sed s/Supercontig_2/Chromosome_2/ > $@.tmp
	mv $@.tmp $@

$(GTF) :
	$(dir_guard)
	curl $(GTF_URL) | zcat > $@.tmp
	mv $@.tmp $@
#--------------------------------------------------------------------------------
# Index reference genome
# $(BT2_INDEX_FILES) : $(SEQ)
# 	bowtie2-build $< $(INDEXBASE)
%.1.bt2 %.2.bt2 %.3.bt2 %.4.bt2 %.rev.1.bt2 %.rev.2.bt2 : %.fa
	bowtie2-build $< $*

%.bwt %.pac %.ann %.amb %.sa : %.fa
	bwa index $< -p $*
#--------------------------------------------------------------------------------
# Crude Filtering of Cassava failed reads
# from http://cancan.cshl.edu/labmembers/gordon/fastq_illumina_filter/

# $(FINAL_FASTQ_DIR)/%.fastq.gz : $(RAW_FASTQ_DIR)/%.fastq.gz
# 	$(dir_guard)
# 	zcat  $^ | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$$" | gzip -c > $@.tmp
# 	mv $@.tmp $@
#--------------------------------------------------------------------------------
# Run Tophat
#-----------
# -T/--transcriptome-only 	Only align the reads to the transcriptome and report only those mappings as genomic mappings.
# -x/--transcriptome-max-hits 	Maximum number of mappings allowed for a read, when aligned to the transcriptome (any reads found with more then this number of mappings will be discarded).
# -M/--prefilter-multihits 	When mapping reads on the transcriptome, some repetitive or low complexity reads that would be discarded in the context of the genome may appear to align to the transcript sequences and thus may end up reported as mapped to those genes only. This option directs TopHat to first align the reads to the whole genome in order to determine and exclude such multi-mapped reads (according to the value of the -g/--max-multihits option). 
# -G known_genes.gtf \
# --transcriptome-index=transcriptome_data/known \
#--------------------------------------------------
$(TOPHAT_BASE_DIR)/%/accepted_hits.bam : $(GTF) $(FINAL_FASTQ_DIR)/%$(FASTQ_SUFFIX) $(BT2_INDEX_FILES)
	$(eval OUTDIR :=  $(@D)_tmpthoutdir)
	mkdir -p $(OUTDIR)
	$(TOPHAT) --output-dir $(OUTDIR) --GTF $(word 1,$^) \
		--transcriptome-max-hits 1 --max-multihits 1 \
		--max-intron-length 4000 --num-threads $(NUMTHREADS) \
		--library-type fr-unstranded --no-coverage-search \
		$(INDEXBASE) $(word 2,$^)
	mv $(OUTDIR) $(@D)

$(TOPHAT_BASE_DIR)/%.bam : $(TOPHAT_BASE_DIR)/%/accepted_hits.bam
	ln -s $< $@

$(TOPHAT_BASE_DIR)/%/test.bam : $(RAW_FASTQ_DIR)/%$(FASTQ_SUFFIX) $(BT2_INDEX_FILES)
	$(eval OUTDIR :=  $(@D)_tmpthoutdir)
	mkdir -p $(OUTDIR)
	echo "blah" > $(OUTDIR)/$(@F)
	mv $(OUTDIR) $(@D)
#--------------------------------------------------------------------------------
# Index BAMs
%.bam.bai : %.bam
	$(SAMTOOLS) index $^
#--------------------------------------------------------------------------------
# Count Reads
#-----------
$(COUNT_DIR)/%_counts.tab : $(TOPHAT_BASE_DIR)/%/accepted_hits.bam $(GTF)
	$(dir_guard)
	$(eval TEMP_OUT := $@.tmp)
	$(eval ERR := $(basename $@).err)
	$(eval TEMP_ERR := $(ERR).tmp)
	$(SAMTOOLS) view $< | htseq-count --stranded=no --type=exon --idattr=gene_id - $(word 2,$^) 1> $(TEMP_OUT) 2> $(TEMP_ERR)
	mv $(TEMP_OUT) $@
	mv $(TEMP_ERR) $(ERR)

#===============================================================================
# Run Analysis in R
#=============================
## wget http://fungalgenomes.org/public/cryptococcus/CryptoDB/product_names/Cneo_H99.AHRD.20131001.tab

$(CNEO_ANNOT) :
	$(dir_guard)
	wget --no-directories --directory-prefix $(@D) http://fungalgenomes.org/public/cryptococcus/CryptoDB/product_names/Cneo_H99.AHRD.20131001.tab

deseq_test : $(CNEO_ANNOT)
	Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R
	printf '\a';printf '\a';printf '\a'
	printf '\a';printf '\a';printf '\a'

deseq_analysis : $(CNEO_ANNOT)
	# Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --label drops7_ --table info/calcineurin_sample_table_drop_crz1s7.csv --fc 2 --fdr 0.2
	# Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --label drops7_ --table info/calcineurin_sample_table_drop_crz1s7.csv --fc 0 --fdr 0.2
	# Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --label drops7_ --table info/calcineurin_sample_table_drop_crz1s7.csv --fc 0 --fdr 0.05
	Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --table info/calcineurin_sample_table.csv --fc 0 --fdr 0.2
	Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --table info/calcineurin_sample_table.csv --fc 0 --fdr 0.05
	printf '\a';printf '\a';printf '\a'
	printf '\a';printf '\a';printf '\a'

deseq_label_test : $(CNEO_ANNOT)
	## Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --usecwd
	Rscript --no-restore $(CNA)/calcineurin_reg_analysis.R --label test
	printf '\a';printf '\a';printf '\a'
	printf '\a';printf '\a';printf '\a'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#--------------------------------------------------------------------------------
# /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
# WORKS ABOVE HERE
# /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
#--------------------------------------------------------------------------------

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# BEGIN IN PROGRESS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#===============================================================================

$(FINAL_FASTQ_DIR)/adapter.fa : $(ILLUMINA_ADAPTERS) $(RAW_FASTQ_DIR)
	$(dir_guard)
	mkdir -p $(TMP_DIR)
	ls -1 $(word 2,$^) | cut -f2 -d_ | sort | uniq > $(INDICES)
	grep -B1 -f $(INDICES) $< > $@.tmp
	grep -B1 -A1 Universal $< >> $@.tmp
	mv $@.tmp $@

$(FINAL_FASTQ_DIR)/%.fastq.gz : $(ADAPTERS_USED) $(RAW_FASTQ_DIR)/%.fastq.gz
	$(dir_guard)
	fastq-mcf -o $(basename $@).TMP.gz $< <(gunzip -c $(word 2,$^);)
	mv $(basename $@).TMP.gz $@

#--------------------------------------------------
$(FIG_DIR)/%_pileup.pdf : $(PILEUP_PLOT_TABLE) $(FINAL_BAMS) $(FINAL_BAIS)
	$(dir_guard)
	# python2.7 $SCRIPTS/plot_pileup.py --gff  genome/h99_2.gtf --linewidth 1 --legendsize 12 --noxlabel --gene CNAG_00156 --subplot --table info/plot_pileup_table.csv -o $CNA/CNAG_00156_pileup.pdf --figsize 8 30 --independenty
	# python2.7 $SCRIPTS/plot_pileup.py --gff  genome/h99_2.gtf --linewidth 1 --legendsize 12 --noxlabel --gene CNAG_00156 --subplot --table info/plot_pileup_table.csv -o $CNA/CNAG_00156_pileup_y10.pdf --figsize 8 30 --ymax 10
	# python2.7 $SCRIPTS/plot_pileup.py --gff  genome/h99_2.gtf --linewidth 1 --legendsize 12 --noxlabel --gene CNAG_00156 --subplot --table info/plot_pileup_table.csv -o $CNA/CNAG_00156_pileup_y5.pdf --figsize 8 30 --ymax 5
	python2.7 $(SCRIPTS)/plot_pileup.py --gff  genome/h99_2.gtf --linewidth 0.5 --legendsize 12 --noxlabel --subplot --figsize 8 30 --table $< -o $@ --gene $*
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# END IN PROGRESS 
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#--------------------------------------------------------------------------------
# \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/ 
# REFERENCE CODE BELOW HERE
# \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/ 
#--------------------------------------------------------------------------------




#===============================================================================
# [josh@gemscompute01 raw_fastqs]$ pwd
# /nfs/gems_sata/heitman/calcineurin_reg/rnaseq/raw_fastqs
# [josh@gemscompute01 raw_fastqs]$ \ls -1
# SC-ECRNA10_TAGCTT_L003_R1_001.fastq.gz
# SC-ECRNA11_GGCTAC_L003_R1_001.fastq.gz
# SC-ECRNA12_CTTGTA_L003_R1_001.fastq.gz
# SC-ECRNA13_AGTCAA_L004_R1_001.fastq.gz
# SC-ECRNA14_AGTTCC_L004_R1_001.fastq.gz
# SC-ECRNA15_ATGTCA_L004_R1_001.fastq.gz
# SC-ECRNA16_CCGTCC_L004_R1_001.fastq.gz
# SC-ECRNA17_GTCCGC_L004_R1_001.fastq.gz
# SC-ECRNA18_GTGAAA_L004_R1_001.fastq.gz
# SC-ECRNA19_GTGGCC_L005_R1_001.fastq.gz
# SC-ECRNA1_ATCACG_L002_R1_001.fastq.gz
# SC-ECRNA20_GTTTCG_L005_R1_001.fastq.gz
# SC-ECRNA21_CGTACG_L005_R1_001.fastq.gz
# SC-ECRNA22_GAGTGG_L005_R1_001.fastq.gz
# SC-ECRNA23_ACTGAT_L005_R1_001.fastq.gz
# SC-ECRNA24_ATTCCT_L005_R1_001.fastq.gz
# SC-ECRNA25_ATCACG_L006_R1_001.fastq.gz
# SC-ECRNA26_CGATGT_L006_R1_001.fastq.gz
# SC-ECRNA27_TTAGGC_L006_R1_001.fastq.gz
# SC-ECRNA28_TGACCA_L006_R1_001.fastq.gz
# SC-ECRNA29_ACAGTG_L006_R1_001.fastq.gz
# SC-ECRNA2_CGATGT_L002_R1_001.fastq.gz
# SC-ECRNA30_GCCAAT_L006_R1_001.fastq.gz
# SC-ECRNA3_TTAGGC_L002_R1_001.fastq.gz
# SC-ECRNA4_TGACCA_L002_R1_001.fastq.gz
# SC-ECRNA5_ACAGTG_L002_R1_001.fastq.gz
# SC-ECRNA6_GCCAAT_L002_R1_001.fastq.gz
# SC-ECRNA7_CAGATC_L003_R1_001.fastq.gz
# SC-ECRNA8_ACTTGA_L003_R1_001.fastq.gz
# SC-ECRNA9_GATCAG_L003_R1_001.fastq.gz
#===============================================================================

#--------------------------------------------------------------------------------
# Index BAMs


$(TOPHAT_BASE_DIR)/%.bam : $(TOPHAT_BASE_DIR)/%/accepted_hits.bam
	$(dir_guard)
	ln -s $(subst $(TOPHAT_BASE_DIR)/,,$<) $@
#===============================================================================
#===============================================================================
#===============================================================================

# Cleanup
#-----------
clean :
	rm -rf $(DESEQ_RESULTS_DIR)

veryclean :
	rm -rf $(FINAL_FASTQ_DIR) $(COUNT_DIR) $(TOPHAT_BASE_DIR)

pristine:
	find . -type d -delete
tidy:
	rm -rf $(STAGE_DIR)/*

rm_th_results:
	rm -rf $(TOPHAT_BASE_DIR)/*
