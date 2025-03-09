#!/bin/Rscript

library(tidyverse)
library(magrittr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

#args <- Sys.glob("/N/p/asclab/ASC-cutNrun/25FEB27-mromanmoreno/src/nextflow_cutNrun/test_nf/replicates/*/align/*_spike.tsv")

lapply(args,function(fl_nm) {
  fread(fl_nm,col.names=c("type","count")) %>%
    as_tibble %>%
    mutate(type = case_when(row_number() == 1 ~ "hg38",
                            row_number() == 2 ~ "sac3",
                            TRUE ~ type) %>%
             factor(levels=c("hg38","sac3","ambiguous","unmapped")))
}) %>% do.call(rbind,.) %>%
  group_by(type) %>%
  summarize(count = sum(count),.groups='drop') %>%
  apply(1,function(row) cat(paste0(row[1],"\t",row[2],"\n"))) %>% invisible
