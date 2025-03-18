#' @title
#' Plot gene track
#'
#' @description
#' Given a genomic range, plots a simple gene track using gene annotations from
#' gtf_to_genes(). Returns a GGplot2 object.
#'
#' @param genomic_region GRange denoting the region to be illustrated.
#' @param gr_genes GRanges object containing gene annotations. If not provided, gtf_to_genes() is used, but this can be time consuming: better to use gtf_to_genes() earlier and provide the output GR to this function directly.
#' @param label_genes Gene names to be labeled. Defaults to NULL, in which case no genes are labeled. Supersedes text_genes.
#' @param label_size Size of text in labels' font. Defaults to 3.
#' @param label_seed Seed number to use for ggrepel labels (used to maintain consistent label positions). Defaults to NULL.
#' @param focus_genes Gene names to highlight (i.e. give an alpha of 1). If given, genes not in this list have their alpha reduced to <background_alpha>. Defaults to NULL, in which case all genes are alpha 1.
#' @param background_alpha Alpha value to use for background (non-<focus-genes>). Defaults to 0.25.
#' @param text_genes Gene names to be labeled with text. Defaults to NULL, in which case no genes are labeled. Genes must have non-empty gene name columns.
#' @param text_biotypes Vector of gene biotypes to be labeled with text. Defaults to "protein_coding" and "lncRNA." If "all" is listed, all biotypes will be labeled.
#' @param text_size Size of text font. Defaults to 2.
#' @param text_bins X-axis is split into <text_bins> "bins" and labels are stacked within each. Fewer bins result in higher stacks of labels, while lower bin numbers may result in more overlaps with labels/geometry. Defaults to 10.
#' @param text_nudge_y Distance above each text bin's highest gene annotation above which to start adding text labels in units of gene layers (one layer = one gene annotation). Defaults to 0.25.
#' @param text_line_height Height of text labels, useful to adjust in cases where larger text sizes are used causing them to overlap. Units of layers, defaults to 0.2.
#' @param text_nudge_x Distance to nudge labels in the negative x-direction when staggering stacked labels. Units of fractions of the X-axis, defaults to 0.005 (text is nudged 0.05 of the x-axis per stacked label).
#' @param x_trace_genes Gene names/IDs which should have traces drawn to the x-axis. Can be individual names/IDs, as well as "all", "focus" (foreground-only), and "label" (labeled genes only).
#' @param x_trace_alpha Opacity of x-axis traces. Defaults to 0.3, can be between zero and 1.
#' @param arrow_genes Should arrows be drawn for gene directions? Defaults to TRUE.
#' @param gr_tads GRanges of TADs to be drawn as diamonds in the background; defaults to NULL (none), and will be plotted in series separated by "name" column.
#' @param tad_height Y-position at which to center TAD diamonds. Defaults to 1.5 (centered on minimum gene level).
#' @param tad_alpha Opacity of TAD diamonds; defaults to 0.1.
#' @param tad_fill Fill of TAD diamonds; defaults to "gray50."
#' @param tad_color Color of TAD diamonds; defaults to NA (no outline).
#' @param tad_linetype Linetype of TADs; defaults to "solid."
#' @param tad_linewidth Linewidth of TADs; defaults to 0.5.
#' @param tad_y_frac Fraction of the y-axis to scale TADs to; defaults to depend on x-range.
#'
#' @import ggplot2
#' @import magrittr
#' @import GenomicRanges
#' @import dplyr
#' @import ggrepel
#' @import scales
#' @import ggsci
#'
#' @export

# gr_genes <- gtf_to_genes() %>% keepStandardChromosomes(pruning.mode="coarse")
# genomic_region <- gr_genes[gr_genes$gene_name == "DLX5"] %>% resize(width=5E6,fix="center")
# text_biotypes="all"
# label_genes <- "DLX5"

plot_genes<- function(genomic_region,gr_genes=NULL,label_genes=NULL,label_size=3,label_seed=NULL,
                      focus_genes=NULL,background_alpha=0.25,
                      text_genes=NULL,text_biotypes=c("protein_coding","lncRNA"),
                      text_size=2,text_bins=10,text_nudge_y=0.25,text_line_height=0.2,text_nudge_x = 0.005,
                      x_trace_genes = "label",x_trace_alpha=0.3,
                      arrow_genes=TRUE,
                      gr_tads=NULL,tad_height=1.5,tad_alpha=0.1,tad_fill="gray50",tad_color=NA,tad_linetype="solid",tad_linewidth=0.5,tad_y_frac=NULL,
                      gr_eigs){
  rename  <- dplyr::rename
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  x_rng   <- c(start=GenomicRanges::start(genomic_region),end=GenomicRanges::end(genomic_region))
  if(is.null(gr_genes)){
    message("Gathering genes from GFF file, provide a GR object to 'gr_genes' in order to speed this up in the future.")
    gr_genes <- gtf_to_genes() %>% keepStandardChromosomes(pruning.mode="coarse")
  }
  gr <- gr_genes[queryHits(findOverlaps(gr_genes,genomic_region))]
  scale_x <- scale_x_continuous(name = grange_desc(genomic_region),
                                oob = scales::oob_keep,
                                limits=x_rng,expand=c(0,0),labels = comma)
  base_plot   <- ggplot(data=NULL) +
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

  if(length(gr) > 0){
    if(is.null(label_genes))label_genes <- ""
    if(is.null(text_genes)) text_genes  <- ""
    if(is.null(focus_genes)) focus_genes<- unique(gr$gene_name)
    if("all" %in% tolower(text_biotypes)) text_biotypes  <- unique(gr$gene_biotype)
    bt_levs <- c("protein_coding","miRNA","lncRNA")
    bt_levs <- c(bt_levs,setdiff(levels(gr_genes$gene_biotype),bt_levs))
    tb_p  <- gr %>%
      as_tibble %>%
      mutate(ymin = pile_coords(start_vals = start,end_vals = end,seqname_vals = seqnames,min_gap=0.001 * diff(x_rng)),
             ymax = ymin + 1,
             start_adj = ifelse(start < x_rng[1],x_rng[1],start),
             end_adj = ifelse(end > x_rng[2],x_rng[2],end),
             start = ifelse(strand == "-",end_adj,start_adj),
             end = ifelse(strand == "-",start_adj,end_adj),
             gene_biotype = factor(gene_biotype,levels=bt_levs)) %>%
      select(-start_adj,-end_adj) %>%
      mutate(start2 = start + 0.2 * (end - start),
             end2 = end - 0.2 * (end - start),
             mid = (start + end) / 2,
             focus_type = ifelse(gene_name %in% focus_genes,"foreground","background"),
             lab_type = case_when(gene_name == "" ~ "none",
                                  gene_name %in% label_genes ~ "label",
                                  gene_name %in% text_genes ~ "text",
                                  gene_biotype %in% text_biotypes ~ "text",
                                  TRUE ~ "none"))

    tb_labs <- filter(tb_p,lab_type != "none")
    if(nrow(tb_labs) > 0){
      tb_labs<- tb_labs %>%
        arrange(desc(start),desc(end)) %>%
        mutate(text_bin = cut(start,breaks=text_bins,dig.lab = 10)) %>%
        mutate(text_bin_start = str_match(as.character(text_bin),"^[\\[\\(]([[:digit:]\\.]+)")[,2] %>% as.double %>% ceiling,
               text_bin_end = str_match(as.character(text_bin),"([[:digit:]\\.]+)[\\]\\)]$")[,2] %>% as.double %>% ceiling) %>%
        group_by(text_bin) %>%
        mutate(nudge_y_lower = row_number() * (text_nudge_y),
               text_y1 = max(ymax) + text_nudge_y + 0.2 * max(tb_p$ymax),
               text_y2 = text_y1 + text_line_height * row_number() * max(ymax),
               text_bin_start = text_bin_start - (row_number() * text_nudge_x * diff(x_rng)),
               text_bin_start = ifelse(text_bin_start < x_rng[1],x_rng[1],text_bin_start)) %>%
        ungroup %>%
        mutate(text_y1 = text_y1 - nudge_y_lower)
    }else{
      tb_labs <- mutate(tb_labs,text_y1 = 0,text_y2 = 0)
    }

    # Scale y up for very large windows to represent "zoomed out" scale.
    min_y_scale <- 6 + ceiling(width(genomic_region)/1E6)/2
    y_rng <- c(0,max(c(min_y_scale,tb_labs$text_y2,tb_p$ymax)))
    scale_y <- scale_y_continuous(name = "Genes",limits=c(-1,y_rng[2]+1),expand=c(0,0),oob = scales::oob_keep)

    #Gene text.
    tb_l <- filter(tb_labs,lab_type == "text")
    gene_text <- NULL
    if(nrow(tb_l) > 0){
      gene_text <- list(
        geom_segment(data=tb_l,
                     mapping=aes(x=text_bin_start,
                                 xend=text_bin_start,
                                 y=text_y2,yend=text_y1,
                                 alpha=focus_type),
                     color="gray75"),
        geom_segment(data=tb_l,
                     mapping=aes(x=text_bin_start,
                                 xend=start,
                                 y=text_y1,yend=ymax,
                                 alpha=focus_type),
                     color="gray75"),
        geom_segment(data=tb_l,
                     mapping=aes(x=start,
                                 xend=start,
                                 y=ymax,yend=(ymin+ymax)/2,
                                 alpha=focus_type),
                     color="gray75"),
        geom_text(data=tb_l,
                  mapping=aes(x=text_bin_start,
                              y=text_y2,
                              label=gene_name,
                              alpha=focus_type),
                  hjust=0,vjust=0,
                  size=text_size)
      )
    }

    #TADs.
    geoms_tads  <- NULL
    if(!is.null(gr_tads)){
      if(is.null(tad_y_frac)){
        gr_wid <- width(genomic_region)
        tad_y_frac<-
          case_when(gr_wid < 1E5 ~ 1,
                    gr_wid < 1E6 ~ 0.8,
                    gr_wid < 1E7 ~ 0.5,
                    gr_wid < 1E8 ~ 0.2,
                    TRUE ~ 0.05)
      }
      tb_tads <- gr_tads %>%
        subsetByOverlaps(genomic_region) %>%
        as_tibble %>%
        select(seqnames,start,end,name) %>%
        mutate(tad_idx=row_number()) %>%
        mutate(x1=start,
               x2=(start+end)/2,
               x3=end,
               x4=x2,
               y1=0,
               y2=(end - start)/2,
               y3=0,
               y4=-y2) %>%
        pivot_longer(cols=c(x1,x2,x3,x4,y1,y2,y3,y4),
                     names_to = c("dim","num"),
                     names_pattern="(.)(.)",
                     values_to = "coord") %>%
        pivot_wider(id_cols = c(seqnames,start,end,name,tad_idx,num),
                    names_from = dim,values_from = coord) %>%
        rowwise %>%
        mutate(off_start = x < x_rng[1],
               off_end = x > x_rng[2],
               y = case_when(off_start ~ list(c(-1,1) * (x_rng[1] - x)/2),
                             off_end ~ list(c(-1,1) * (x - x_rng[2])/2),
                             TRUE ~ list(y)),
               x = case_when(off_start ~ list(rep(x_rng[1],length.out=2)),
                             off_end ~ list(rep(x_rng[2],length.out=2)),
                             TRUE ~ list(x))) %>%
        ungroup %>%
        group_by(tad_idx) %>%
        mutate(off_start = cumsum(off_start),
               off_end = cumsum(off_end)) %>%
        mutate(keep = off_start <= max(off_start) & off_end <= 1) %>%
        ungroup %>%
        filter(keep) %>%
        select(-off_start,-off_end,-keep) %>%
        unnest(c(x,y)) %>%
        group_by(name,tad_idx) %>%
        mutate(num = row_number()) %>%
        ungroup %>%
        mutate(y = scale_to(y,-y_rng[2]*tad_y_frac,y_rng[2]*tad_y_frac),
               y = y + tad_height)

      # base_plot +
      #   scale_y +
      #   scale_x +
      #   geom_polygon(data=tb_tads,mapping=aes(x=x,y=y,color=name,fill=name,group=tad_idx),alpha=0.1,color=NA) +
      #   theme(legend.position='right')
      if(nrow(tb_tads) > 0){
        geoms_tads <-
          geom_polygon(data=tb_tads,
                       mapping=aes(x=x,y=y,color=name,fill=name,group=tad_idx),
                       alpha=tad_alpha,color=tad_color,fill=tad_fill,linetype=tad_linetype,linewidth=tad_linewidth,
                       show.legend=FALSE)
      }
    }

    #X-axis traces.
    poly_gns <- x_trace_genes
    if("all" %in% tolower(x_trace_genes)){
      poly_gns<- c(poly_gns,unique(tb_p$gene_id))
    }
    if("label" %in% tolower(x_trace_genes)){
      poly_gns<- c(poly_gns,filter(tb_p,lab_type != "none")$gene_id)
    }
    if("focus" %in% tolower(x_trace_genes)){
      poly_gns<- c(poly_gns,filter(tb_p,focus_type == "foreground")$gene_id)
    }
    geom_segs <- NULL
    geom_polys<- NULL
    tb_poly <- tb_p %>%
      filter(gene_name %in% unique(poly_gns) | gene_id %in% unique(poly_gns)) %>%
      mutate(gene_biotype = factor(as.character(gene_biotype),levels=levels(tb_p$gene_biotype)))

    if(nrow(tb_poly) > 0){
      tb_poly <- tb_poly %>%
        mutate(x_frac = width/width(genomic_region),
               x_polygon = x_frac > 0.01,
               idx = row_number()) %>%
        select(gene_name,gene_id,gene_biotype,start,end,mid,x_polygon,ymin,focus_type,lab_type,idx)

      tb_segs <- filter(tb_poly,!x_polygon)
      if(nrow(tb_segs) > 0){
        geom_segs <- geom_segment(
          data=tb_segs,
          mapping=aes(x=mid,xend=mid,y=ymin+0.1,yend=-1),
          color = "gray75",alpha=x_trace_alpha)
      }
      tb_poly <- filter(tb_poly,x_polygon)
      if(nrow(tb_poly) > 0){
        tb_poly <- tb_poly %>%
          rowwise %>%
          mutate(x_pos = list(c(start,end,mid)),
                 y_pos = list(c(ymin+0.1,ymin+0.1,-1))) %>%
          ungroup %>%
          unnest(c(x_pos,y_pos))
        geom_polys <- geom_polygon(data=tb_poly,mapping=aes(x=x_pos,y=y_pos,fill=gene_biotype,group=idx),alpha=x_trace_alpha)
      }
    }

    base_plot  <- base_plot +
      geoms_tads +
      gene_text +
      scale_alpha_manual(values=c(foreground=1,background=background_alpha)) +
      scale_color_manual(values=c(foreground="black",background=NA)) +
      scale_x +
      scale_y +
      geom_polys +
      geom_segs +
      geom_rect(data=tb_p,
                mapping=aes(xmin = start,xmax=end,
                            ymin = ymin+0.1, ymax=ymax-0.1,
                            fill = gene_biotype,
                            alpha=focus_type,
                            color=focus_type),
                linewidth=0.25)

    if(arrow_genes){
      base_plot <- base_plot +
      geom_segment(data=filter(tb_p,strand != "*" & focus_type != "background" & lab_type != "none"),
                   mapping=aes(x=start2,xend = end2,
                               y=(ymin + ymax)/2,yend=(ymin + ymax) / 2,
                               alpha=focus_type,
                               color=focus_type),
                   linewidth=1,
                   arrow = arrow(length = unit(0.5,"lines"),type = "closed"))
    }

    #Gene labels.
    tb_l <- filter(tb_labs,lab_type == "label")
    if(nrow(tb_l) > 0){
      base_plot <- base_plot +
        geom_label_repel(
          data=tb_l,
          seed = label_seed,
          mapping=aes(x=(start + end) / 2,
                      y=ymin+0.1,
                      label = gene_name,
                      alpha=focus_type),
          min.segment.length = 0.01,
          size=label_size,ylim = c(0,-0.5))
    }
  }
  return(base_plot)
}
