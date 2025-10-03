#' @title Plot Gene Track from GTF
#'
#' @description
#' Given a genomic range and a GTF file of feature annotations, plots a simple
#' gene track using gene annotations from gtf_extract_features(). If no gtf is
#' provided (gtf_file=NULL), uses default GTF file from gtf_extract_features().
#' Returns a GGplot2 object unless <return_tables> is TRUE, in which case a
#' list of tables for features, scales, and gene labels will be returned
#' instead.
#'
#' Gene annotations include exons, and each gene represents the highest
#' confidence transcript for that gene/ensembl_id combination. In cases where
#' an exon's width is less than 1% of the x-range, that exon is drawn with a
#' vertical line to ensure all exons are still displayed.
#'
#' @param genomic_region GRange denoting the region to be illustrated.
#' @param gtf_file GTF file to illustrate. If NULL (default), uses default GTF for gtf_extract_features().
#' @param max_support_level Max transcript support level allowed (1 == TSL1 == highest confidence). Defaults to 7, in which case all genes are shown.
#' @param focus_genes Gene names / EnsemblIDs to highlight by giving all other features a background alpha. Defaults to NULL, in which case all genes have an alpha of 1.
#' @param background_alpha Alpha value to use for background (non-<focus-genes>). Defaults to 0.25.
#' @param intron_fill Color to fill introns with. Defaults to lightblue.
#' @param exon_fill Color to fill exons with. Defaults to darkblue.
#' @param label_color Color to draw segments/errows/text for gene labels with. Defaults to gray20.
#' @param segment_linewidth Linewidths to draw segments/arrows with. Defaults to 0.25.
#' @param text_size Text size for gene names. Defaults to 1.5.
#' @param draw_arrows Boolean; should arrows be included? Defaults to TRUE.
#' @param arrow_sweep Distance to extend arrows upstream from TSS when drawing looped direction arrows. Units of X-range, defaults to 0.005 (0.5% of x-range).
#' @param arrow_size Size in cm for arrows to be drawn (arrow(length=unit(<arrow_size>,"cm"))). Defaults to 0.25.
#' @param min_y_ax Minimum y-axis height. Ensures that in cases where all annotations are on one level the plot area doesn't look too stretched. Defaults to 3, in which case no gene track will cover more than 30% of the panel's width.
#' @param min_gap Minimum gap between genes before they are moved to separate Y-tracks. Defaults to 0.1 (10% of the X-axis).
#' @param feature_colors Colors to use for various feature types.
#' @param include_genes Names/EnsIDs of specific genes to include in the specified region. Useful for subsetting genes in the GTF file.
#' @param exclude_genes Names/EnsIDs of specific genes to exclude in the specified region. Useful for subsetting genes in the GTF file.
#' @param return_tables Boolean; should tables of data and label positions be returned instead of a plot? Defaults to FALSE.
#'
#' @import ggplot2
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import ggrepel
#' @import scales
#'
#' @export

plot_genes_from_gtf <- function(genomic_region,gtf_file=NULL,max_support_level=7,focus_genes=NULL,background_alpha=0.25,
                                intron_fill="lightblue",exon_fill="darkblue",label_color="gray20",segment_linewidth=0.5,
                                text_size=1.5,draw_arrows=TRUE,arrow_sweep=0.005,arrow_size=0.25,min_y_ax = 3,min_gap=0.1,
                                include_genes=NULL,exclude_genes=NULL,return_tables=FALSE){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  arrange <- dplyr::arrange
  first   <- dplyr::first

  x_rng   <- c(start=GenomicRanges::start(genomic_region),end=GenomicRanges::end(genomic_region))

  scale_x <- scale_x_continuous(name = grange_desc(genomic_region),
                                oob = scales::oob_keep,
                                limits=x_rng,expand=c(0,0),labels = comma)

  if(length(genomic_region) > 1) stop("plot_transcripts() is intended for contiguous regions only.")

  gr_feats<- gtf_extract_features(genomic_region,gtf_file = gtf_file,retain_attributes=c("gene_name","gene_id","transcript_id","exon_id","exon_number","transcript_support_level","gene_version"))
  if(nrow(as_tibble(gr_feats)) == 0){
    p_out <- ggplot() + theme_void()
  }else{
    seqlevelsStyle(gr_feats) <- "UCSC"
    gr_feats<- clugPac::clip_granges(grange_in = gr_feats,grange_window = genomic_region)
  
    focus_gns <- focus_genes %||% c(unique(gr_feats$gene_name),unique(gr_feats$gene_id))
    inclu_gns <- include_genes %||% unique(gr_feats$gene_id)
    exclu_gns <- exclude_genes %||% ""
  
    tb_feats  <- gr_feats %>%
      as_tibble %>%
      filter(!is.na(transcript_id)) %>%
      filter( gene_name %in% inclu_gns | gene_id %in% inclu_gns) %>%
      filter(!(gene_name %in% exclu_gns) & !(gene_id %in% exclu_gns)) %>%
      rename(support=transcript_support_level) %>%
      mutate(support=ifelse(support %in% as.character(1:5),support,"6"),
             support=as.integer(support)) %>%
      group_by(seqnames,gene_name,gene_id) %>%
      arrange(support) %>%
      filter(transcript_id == dplyr::first(transcript_id)) %>%
      ungroup %>%
      select(-source,-score,-frame,-gene_version,-start_adj,-end_adj)
  
    tb_sum <- tb_feats %>%
      group_by(seqnames,gene_name,gene_id,transcript_id) %>%
      summarize(gene_min = min(start),
                gene_max = max(end),
                support = min(support),
                .groups='keep') %>%
      mutate(width = gene_max - gene_min) %>%
      ungroup %>%
      filter(support <= max_support_level) %>%
      arrange(gene_name,support,desc(width),gene_min,gene_max) %>%
      mutate(y = pile_coords(gene_min,gene_max,min_gap = min_gap * diff(x_rng)),
             focus = gene_name %in% focus_gns | gene_id %in% focus_gns) %>%
      select(gene_name,gene_id,transcript_id,gene_min,gene_max,support,focus,y)
  
    # feat_cols <- feature_colors %||% c(transcript="gray70",
    #                                    exon = "black",
    #                                    CDS = "blue",
    #                                    three_prime_utr="red",
    #                                    five_prime_utr ="green",
    #                                    stop_codon = "darkred",
    #                                    start_codon="darkgreen") %>% enframe(name="feature",value="color")
  
    tb_p <- tb_feats %>%
      as_tibble %>%
      select(-support) %>%
      right_join(tb_sum,by=c("gene_name","gene_id","transcript_id")) %>%
      # inner_join(feat_cols,by="feature") %>%
      mutate(# feature = factor(feature,levels = feat_cols$feature),
             x_frac = width / diff(x_rng)) %>%
      filter(!is.na(feature))
  
    tb_labs <- tb_p %>%
      filter(feature == "transcript") %>%
      select(seqnames,start,end,y,width,strand,gene_name,gene_id,focus) %>%
      mutate(x1 = ifelse(strand == "+",start,end),
             x2 = ifelse(strand == "+",start-arrow_sweep*diff(x_rng),end+arrow_sweep*diff(x_rng)),
             x3 = ifelse(strand == "+",start+arrow_sweep*diff(x_rng),end-arrow_sweep*diff(x_rng)),
             y1 = y,
             y2 = ifelse(strand=="+",y + 0.2,y-0.2),
             y3 = ifelse(strand=="+",y + 0.25,y-0.25)) %>%
      mutate(y3 = case_when(x2 < x_rng[1] ~ y2,
                            x2 > x_rng[2] ~ y2,
                            TRUE ~ y3),
             x2 = case_when(x2 < x_rng[1] ~ x1,
                            x2 > x_rng[2] ~ x1,
                            TRUE ~ x2))
  
    tb_scale <- tibble(xmin=x_rng[1],xmax=x_rng[2],
                       ymin=1,ymax=max(c(max(tb_p$y+1),min_y_ax)))
    if(return_tables){
      p_out <- list(features=tb_p,labels=tb_labs,scales=tb_scale)
    }else{
      geoms_arrows <- NULL
      if(draw_arrows){
        geoms_arrows <- list(
          geom_segment(data=tb_labs,mapping=aes(x=x1,xend=x2,y=y1,yend=y2,alpha=focus),
                       inherit.aes=FALSE,
                       color=label_color,linewidth=segment_linewidth),
          geom_segment(data=tb_labs,mapping=aes(x=x2,xend=x2,y=y2,yend=y3,alpha=focus),
                       inherit.aes=FALSE,
                       color=label_color,linewidth=segment_linewidth),
          geom_segment(data=filter(tb_labs,strand == "+"),
                       mapping=aes(x=x2,xend=x3,y=y3,yend=y3,alpha=focus),
                       inherit.aes=FALSE,
                       color=label_color,
                       linewidth=segment_linewidth,
                       arrow = arrow(type="closed",length=unit(arrow_size,"cm"))),
          geom_segment(data=filter(tb_labs,strand == "-"),
                       mapping=aes(x=x2,xend=x3,y=y3,yend=y3,alpha=focus),
                       inherit.aes=FALSE,
                       color=label_color,
                       linewidth=segment_linewidth,
                       arrow = arrow(type="closed",length=unit(arrow_size,"cm")))
        )
      }
      p_out <- ggplot(tb_p,aes(xmin=start,xmax=end,ymin=y-0.2,ymax=y+0.2,
                      x=(start+end)/2,xend=(start+end)/2,
                      y=y-0.2,yend=y+0.2,
                      alpha=focus,label=gene_name,fill=color,color=color)) +
        scale_x +
        scale_y_continuous(name="Genes",expand=expansion(mult=c(0.3,0.3))) +
        scale_fill_identity() +
        scale_color_identity() +
        scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=background_alpha)) +
        geom_rect(data=tb_scale,fill=NA,color=NA,
                  mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes=FALSE) +
        geoms_arrows +
        geom_rect(data=filter(tb_p,feature == "transcript"),
                  mapping=aes(ymin=y-0.1,ymax=y+0.1),
                  fill=intron_fill,color=NA) +
        geom_rect(data=filter(tb_p,feature == "exon"),
                  fill=exon_fill,color=NA) +
        # Use segment to make sure narrow exons (< 1% of X-axis) are still represented.
        geom_segment(data=filter(tb_p,feature == "exon" & x_frac < 0.001),
                     mapping=aes(x=(start+end)/2,xend=(start+end)/2,
                                 y=y-0.2,yend=y+0.2,alpha=focus),
                     color=exon_fill,inherit.aes=FALSE) +
        geom_text_repel(data=filter(tb_labs,strand=="+"),mapping=aes(x=x2,y=y3,label=gene_name,alpha=focus),
                  hjust=0,vjust=1,direction = "y",inherit.aes=FALSE,color=label_color,
                  size=text_size) +
        geom_text_repel(data=filter(tb_labs,strand=="-"),mapping=aes(x=x2,y=y3,label=gene_name,alpha=focus),
                  hjust=1,vjust=0,direction = "y",inherit.aes=FALSE,color=label_color,
                  size=text_size) +
        guides(alpha="none",
               color="none") +
        theme(plot.background = element_rect(fill="white",color=NA),
              plot.margin = unit(c(0,0,0,0),"lines"),
              panel.background = element_blank(),
              panel.border = element_rect(linewidth=0.5,color="black",fill=NA),
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(linewidth=0.5,color="gray",linetype = "dotted"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())
    }
  }
  return(p_out)
}
