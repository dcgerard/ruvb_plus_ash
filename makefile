## The data output from gtex_extract_tissues_v6p.R
tissue_dat = ./output/gtex_tissue_gene_reads_v6p/adiposetissue.csv \
	     ./output/gtex_tissue_gene_reads_v6p/bladder.csv \
             ./output/gtex_tissue_gene_reads_v6p/bloodvessel.csv \
             ./output/gtex_tissue_gene_reads_v6p/breast.csv \
             ./output/gtex_tissue_gene_reads_v6p/colon.csv \
             ./output/gtex_tissue_gene_reads_v6p/fallopiantube.csv \
             ./output/gtex_tissue_gene_reads_v6p/kidney.csv \
             ./output/gtex_tissue_gene_reads_v6p/lung.csv \
             ./output/gtex_tissue_gene_reads_v6p/nerve.csv \
             ./output/gtex_tissue_gene_reads_v6p/pancreas.csv \
             ./output/gtex_tissue_gene_reads_v6p/prostate.csv \
             ./output/gtex_tissue_gene_reads_v6p/skin.csv \
             ./output/gtex_tissue_gene_reads_v6p/spleen.csv \
             ./output/gtex_tissue_gene_reads_v6p/testis.csv \
             ./output/gtex_tissue_gene_reads_v6p/uterus.csv \
             ./output/gtex_tissue_gene_reads_v6p/adrenalgland.csv \
             ./output/gtex_tissue_gene_reads_v6p/blood.csv \
             ./output/gtex_tissue_gene_reads_v6p/brain.csv \
             ./output/gtex_tissue_gene_reads_v6p/cervixuteri.csv \
             ./output/gtex_tissue_gene_reads_v6p/esophagus.csv \
             ./output/gtex_tissue_gene_reads_v6p/heart.csv \
             ./output/gtex_tissue_gene_reads_v6p/liver.csv \
             ./output/gtex_tissue_gene_reads_v6p/muscle.csv \
             ./output/gtex_tissue_gene_reads_v6p/ovary.csv \
             ./output/gtex_tissue_gene_reads_v6p/pituitary.csv \
             ./output/gtex_tissue_gene_reads_v6p/salivarygland.csv \
             ./output/gtex_tissue_gene_reads_v6p/smallintestine.csv \
             ./output/gtex_tissue_gene_reads_v6p/stomach.csv \
             ./output/gtex_tissue_gene_reads_v6p/thyroid.csv \
             ./output/gtex_tissue_gene_reads_v6p/vagina.csv

all: sims

## run simulations
.PHONY : sims
sims : $(tissue_dat) ./code/sample_sims.R
	mkdir -p output/sims_out
	Rscript ./code/sample_sims.R
