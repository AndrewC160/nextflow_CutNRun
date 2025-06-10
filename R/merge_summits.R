#!/bin/Rscript
# Flatten summits, prioritizing those with the highest fold enrichment in cases where multiple 
# summits overlap.

options(warn=-1,show.error.messages = FALSE)

library(dplyr)
library(magrittr)
library(data.table)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)

# args <- c("/N/p/asclab/ASC-cutNrun/25APR26-antibody_survey/work/96/02e04d2080889eaa4788e62b0a0864/OS186_IgG_NONE_2_all_summits.bed","/N/p/asclab/ASC-cutNrun/25APR26-antibody_survey/work/96/02e04d2080889eaa4788e62b0a0864/OS186_IgG_NONE_2_summit_regions.bed","1001")

fl_nm <- args[1]
fl_out<- args[2]
reg_width <- args[3] %>% as.integer

gr_sums <- fread(fl_nm,col.names = c("seqnames","start","end","name","score")) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
  resize(width=reg_width,fix='center')

gr_flat <- reduce(gr_sums)
mcols(gr_flat)$idx <- c(1:length(gr_flat))

olaps <- findOverlaps(gr_sums,gr_flat)
gr_out<- gr_sums[queryHits(olaps)]
mcols(gr_out)$idx <- mcols(gr_flat)$idx[subjectHits(olaps)]

gr_out %>%
  as_tibble %>%
  group_by(idx) %>%
  arrange(desc(score)) %>%
  filter(row_number() == 1) %>%
  ungroup %>%
  select(seqnames,start,end,width,name,score) %>%
  arrange(seqnames,start,end) %>%
  fwrite(sep="\t",file=fl_out,quote=FALSE,row.names=FALSE,col.names = FALSE) %>%
  invisible

