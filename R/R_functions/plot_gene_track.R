#' @title
#' Plot gene track.
#'
#' @description
#' Plot a single-panel gene track.
#'
#' @param grange_win GRange object describing the genomic region to plot. Must be one contiguous region.
#' @param gene_cache_file Cache file location for gene data...DOES NOT cache plot-specific subset of genes, but all genes found by genes_and_exons, so for multiple plots use the same cache file.
#' @param overwrite_cache Should gene cache be overwritten? Defaults to FALSE.
#' @param plot_biotypes List of gene biotypes to plot, defaults to all biotypes. Biotypes not featured in this list are not plotted (labels or otherwise).
#' @param point_size Size of geom_points in gene track. Defaults to 6.
#' @param text_lab_size Size of text labels on plot. Defaults to 4.
#' @param gene_arrows Plot arrows to show gene strand? Defaults to TRUE.
#' @param label_genes List of genes to label with gene names (if present). If "all" is included in this list, all genes are labeled. If "oncogenes" is included, all genes that arelisted in the Bushman lab's oncogene list are included. Can also be "protein_coding".
#' @param min_gene_gap Minimum gap between genes before they are incremented to another track; defaults to NULL (i.e. 1bp apart).
#'
#' @import ggplot2
#' @import tidyr
#' @import GenomicRanges
#' @import tibble
#' @import clugPac
#' @import data.table
#' @import magrittr
#' @import dplyr
#' @import scales
#' @import ggrepel
#' @import ggsci
#'
#' @export

plot_gene_track  <- function(grange_win,gene_cache_file=NULL,overwrite_cache=FALSE,
                             plot_biotypes=NULL,point_size=2,text_lab_size=4,gene_arrows=TRUE,
                             gene_height_frac=0.7,label_genes="all",min_gene_gap=NULL){
  if(length(grange_win) > 1) stop("<grange_win> should be a GRange of length 1.")
  x_rng   <- c(start=start(grange_win),end=end(grange_win))
  tb_gns  <- genes_and_exons(grange_win = grange_win,clip_features=TRUE,
                             cache_gene_tsv = gene_cache_file,
                             overwrite_cache = overwrite_cache)
  tb_exn  <- tb_gns$exons
  tb_gns  <- tb_gns$genes

  #Get palette for biotypes.
  gene_pal<- biotype_palette()
  gene_pal<- gene_pal[names(gene_pal) %in% unique(tb_gns$gene_biotype)]
  g_theme <-theme(
    panel.border = element_rect(linewidth=0.25,color="black",fill=NA),
    panel.background = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linewidth=0.25,color="black",linetype = "dotted"),
    plot.margin = unit(c(0,0,-0.2,0),"lines"),
    plot.background = element_rect(color=NA,fill="white"),
    legend.direction = "horizontal",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle=0,hjust=1,vjust=0.5))

  if(nrow(tb_gns) == 0){
    p <- ggplot() +
      scale_x_continuous(name = grange_desc(grange_win),
                         limits=x_rng,expand=c(0,0),
                         labels = comma) +
      scale_y_continuous(expand = c(0,0)) +
      scale_color_manual(name="Biotype",values=gene_pal) +
      scale_fill_manual(name="Biotype",values=gene_pal) +
      g_theme
  }else{
    if("all" %in% label_genes){
      lab_genes <- unique(tb_gns$gene_name)
    }else{
      lab_genes   <- label_genes
      if("oncogenes" %in% label_genes){
        lab_genes <- c(lab_genes,
                       filter(tb_gns,bushman_onco) %>% select(gene_name) %>% unlist %>% unique)
      }
      if("protein_coding" %in% label_genes){
        lab_genes <- c(lab_genes,
                       filter(tb_gns,gene_biotype == "protein_coding") %>% select(gene_name) %>% unlist %>% unique)
      }
    }
    if(is.null(plot_biotypes)){
      plot_biotypes   <- unique(as.character(tb_gns$gene_biotype))
    }

    tb_gns  <- tb_gns %>%
      filter(gene_biotype %in% plot_biotypes) %>%
      filter(end > start) %>%
      mutate(ymax = -pile_coords(start_vals = start,end_vals = end,min_gap = min_gene_gap),
             ymin = ymax + 1,
             gn_start = ifelse(strand == "+",start,end),
             gn_end = ifelse(strand == "+",end,start),
             start = gn_start,
             end = gn_end,
             label = ifelse(gene_name %in% lab_genes | gene_id %in% lab_genes,gene_name,"")) %>%
      select(-gn_start,-gn_end) %>%
      mutate(ymin = ymin - (1-gene_height_frac/2), ymax=ymax + (1-gene_height_frac/2))

    if(!inherits(tb_exn,"data.frame")){
      geom_exn<- NULL
    }else{
      tb_exn  <- tb_exn %>%
        mutate(arr_start = ifelse(strand == "+",start + 0.1 * width,end - 0.1 * width),
               arr_end = ifelse(strand == "+",end - 0.1 * width,start + 0.1 * width)) %>%
        select(seqnames,start,end,arr_start,arr_end,width,strand,gene_id) %>%
        inner_join(by="gene_id",select(tb_gns,gene_name,gene_id,ymin,ymax,gene_biotype))
      geom_exn <- geom_rect(data=tb_exn,alpha=0.4)
    }

    p <- ggplot(tb_gns,
          aes(xmin=start,xmax=end,
              x=start,xend=end,
              ymin=ymin,ymax=ymax,
              y=(ymin + ymax)/2,yend=(ymin+ymax)/2,
              fill=gene_biotype,color=gene_biotype)) +
      scale_x_continuous(name = grange_desc(grange_win),
                         limits=x_rng,expand=c(0,0),
                         labels = comma) +
      scale_y_continuous(name = "Genes",expand = c(0,0)) +
      scale_color_manual(name="Biotype",values=gene_pal) +
      scale_fill_manual(name="Biotype",values=gene_pal) +
      geom_rect(color=NA,alpha=0.2,show.legend=FALSE) +
      geom_segment(show.legend=TRUE) +
      geom_exn +
      #geom_rect(data=tb_exn,alpha=0.4) +
      g_theme
    tb_labs <- filter(tb_gns,label!="")
    if(nrow(tb_labs) > 0){
      p <- p + geom_text_repel(tb_labs,mapping=aes(x=(start+end)/2,y=(ymin + ymax)/2,label=label),
                               size=text_lab_size,color="black",nudge_y = -0.4,min.segment.length = 0.01,show.legend=FALSE)
    }

    if(gene_arrows){
      tb_g  <- filter(tb_gns,strand == "+")
      if(nrow(tb_g) > 0){
        p <- p +
          geom_point(data=tb_g,mapping=aes(x=(start + end)*0.5,y=(ymin+ymax)/2,color=gene_biotype),
                     show.legend=FALSE,shape = "\u25BA",size=point_size)
      }
      tb_g  <- filter(tb_gns,strand == "-")
      if(nrow(tb_g) > 0){
        p <- p +
          geom_point(data=tb_g,mapping=aes(x=(start + end)*0.5,y=(ymin+ymax)/2,color=gene_biotype),
                     show.legend=FALSE,shape = "\u25C4",size=point_size)
      }
    }
  }
  return(p)
}
