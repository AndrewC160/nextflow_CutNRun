#' @title
#' Plot summary region.
#'
#' @description
#' Plot treat/control BedGraph signal along with peak and summit annotations for
#' files annotated within a summary.tsv file. A pre-selected subset of a 
#' bedGraph table (from read_summary_bedgraphs()) can be included instead, but 
#' if neither is provided the function will exit with an error. Peak annotations
#' are also read from the <summary_file>, but if a GRanges object formatted as 
#' is the output of read_summary_peaks() is provided this will be used instead.
#' If no peaks are to be illustratred, provide an empty GRanges object.
#' 
#' Gene annotations can be provided as a GTF file, which will be read using 
#' gtf_to_genes(). A cache file (<gene_tsv>) can be provided here to speed up 
#' the retrieval process. If neither is provided, a gene facet will not be 
#' included.
#'
#' "Chr1: 1,000,000 - 2,000,000 (1MB)"
#'
#' @param grange_in GRange to describe, should be of length one (if not, only the first will be described.)
#' @param force_unit Enforce a desired unit (bp, kb, MB) for width readout.
#' @param digits Number of digits to display, defaults to 1.
#' @param sep Separation between width value and units, defaults to "".
#' @param append_ending Insert a string between the width printout and the closing parenthesis; intended for bin size if relevant (for instance, "(1MB, 2kb bin size)").
#'
#' @import stringr
#' @import GenomicRanges
#'
#' @export

plot_summary_region <- function(gr_window,summary_file=NULL,bdg_table=NULL,peak_gr=NULL,gene_gr=NULL,gene_gtf=NULL,gene_tsv=NULL){
  if(is.null(bdg_table)){
    if(is.null(summary_file)) stop("At minimum, either a summary.tsv file or a bedGraph table from read_summary_bedgraphs() is required.")
    bdg_table <- read_summary_bedgraphs(summary_file)
  }
  if(is.null(peak_gr)){
    gr_pks <- read_summary_peaks(summary_file)
  }else{
    gr_pks <- peak_gr
  }
  if(is.null(gene_gr)){
    if(!is.null(gene_tsv) | !is.null(gene_gtf)){
      gene_gr <- gtf_to_genes(gene_gtf,gns_tsv,overwrite_cache = FALSE)
    }
  }
  x_rng   <- c(start(gr_window),end(gr_window))
  tb_sig  <- bdg_table %>%
    apply(1,function(row_in){
      scan_bdg(row_in[['file']],gr_regions = gr_window) %>%
        mutate(name = row_in[['name']],
               set = row_in[['set']],
               cond = row_in[['cond']],
               spike = as.double(row_in[['spike']]))
    }) %>% do.call(rbind,.) %>%
    mutate(score = ifelse(cond == "control",-score,score),
           score_adj = ifelse(cond == "treat",score * spike,score)) %>%
    group_by(name) %>%
    mutate(x = (start + end)/2,
           x = case_when(start == min(start) ~ x_rng[1],
                         end == max(end) ~ x_rng[2],
                         TRUE ~ x))
  
  tb_scales <- tb_sig %>%
    group_by(name) %>%
    summarize(xmin=x_rng[1],
              xmax=x_rng[2],
              ymin=1.05 * min(score),
              ymin = ifelse(ymin >0,0,ymin),
              ymax=1.1*max(score),.groups='drop') %>%
    mutate(ystep = 0.1 * (ymax - ymin))
  
  tb_pks <- gr_pks %>%
    subset_and_trim_gr(gr_window) %>%
    as_tibble %>%
    select(-strand) %>%
    left_join(select(tb_scales,name,ymin,ymax,ystep),by="name") %>%
    mutate(y1 = case_when(set == "summits" ~ ymax,
                          set == "narrowPeaks" ~ ymax + ystep,
                          set == "broadPeaks" ~ ymax + 2 * ystep),
           y2 = y1 + ystep)
  tb_sums <- filter(tb_pks,set == "summits")
  tb_pks  <- filter(tb_pks,set != "summits")
  
  geom_pks <- NULL
  geom_sums <- NULL
  if(nrow(tb_sums) > 0){
    geom_sums <- geom_point(data=tb_sums,mapping=aes(x=(start+end)/2,y=y1),color="black",alpha=1,pch=8,inherit.aes=FALSE)
  }
  if(nrow(tb_pks) > 0){
    geom_pks  <- geom_rect(data=tb_pks,mapping=aes(xmin=start,xmax=end,ymin=y1,ymax=y2,fill=set),color='black',alpha=1,linewidth=0.1,inherit.aes=FALSE)
  }
  
  p <- ggplot(tb_sig,aes(x=x,ymin=0,ymax=score,alpha=cond)) +
    facet_wrap(.~name,scales="free_y",strip.position = "right",ncol=1) +
    scale_x_continuous(name= grange_desc(gr_window),expand=c(0,0),limits=x_rng,labels=prettyBP) +
    scale_y_continuous(name= "Pileup",labels=function(x) ifelse(x>0,prettyNumbers(x),""),n.breaks = 3) +
    scale_alpha_manual(values=c(treat=1,control=0.3)) +
    geom_ribbon(fill="gray25",color='black') +
    geom_pks +
    geom_sums +
    new_scale_fill() +
    geom_rect(data=tb_scales,mapping=aes(ymax=ymax),xmin=-Inf,xmax=Inf,ymin=0,fill=NA,color=NA,inherit.aes=FALSE) +
    theme(plot.background = element_rect(fill="white",color=NA),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black"),
          panel.grid = element_blank(),
          panel.spacing.y = unit(0,"lines"),
          strip.text.y.right = element_text(angle=0,hjust=0),
          strip.background = element_blank(),
          legend.title = element_blank()) %>% suppressWarnings()
  
  if(!is.null(gene_gr)){
    p_gns <- plot_genes(genomic_region = gr_window,gr_genes = gene_gr) +
      theme(legend.position = 'none',
            panel.grid = element_blank())
    
    p <- cowplot::plot_grid(p_gns,p,ncol=1,align="v",axis='lr',rel_heights=c(1,5))
  }
  return(p)
}
