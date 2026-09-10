#!/usr/bin/env Rscript
# Purpose: original full-cohort pair-distance boxplots and bubble plot.
# Inputs: stage-03 pairs and stage-04 summaries.
# Outputs: PDFs plus bubble_plot_data.csv auditing rounded plotting p-values.
args_all <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=",args_all)][1]))
source(file.path(dirname(script),"R","common.R"))

suppressPackageStartupMessages({library(ggplot2);library(ggpubr)})
cfg <- cli_config(); record_stage(cfg,"05")
x <- readRDS(file.path(cfg$paths$output,"pairwise_distances.rds"))
s <- readRDS(file.path(cfg$paths$output,"similarity_summary.rds"))
s <- s[s$run_id=="FULL",]; s <- s[match(cfg$cell_type_order,s$cell_type),]
# Recover historical loop order from network index, distinct from display order.
idx <- readRDS(file.path(cfg$paths$output,"network_index.rds"))
loop_order <- idx$cell_type[idx$run_id=="FULL"]
p_annotation <- character(length(loop_order))
for (k in seq_along(loop_order)) {
  ct <- loop_order[k]; d <- x[x$run_id=="FULL" & x$cell_type==ct,]
  d$activity <- factor(c(RR="L_L",R_NR="H_L",NRNR="H_H")[d$group],levels=c("L_L","H_L","H_H"))
  p <- ggboxplot(d,x="activity",y="distance",fill="activity",outlier.colour=NA) +
    ylab("euclidean distance") + ggtitle(paste0("euclidean distance of ",ct," networks")) +
    theme_classic() + scale_fill_manual(values=c("#D00000","#FFBA08","#3F88C5")) +
    stat_compare_means(comparisons=list(c("H_H","L_L")),method="wilcox.test",label="p.format",
      hide.ns=TRUE,bracket.size=0.6,tip.length=0.01,label.y=max(d$distance)*1.05)
  p_annotation[k] <- as.character(ggplot_build(p)$data[[2]]$annotation[1])
  ggsave(file.path(cfg$paths$output,paste0("response_revised_p_value_New_euclidean_",sanitize(ct),"_disease_activity_personal_network.pdf")),p,width=7,height=7)
}
s_loop <- s[match(loop_order,s$cell_type),]
b <- data.frame(cell_type=rep(loop_order,each=3),group=rep(c("RR","NRR","NRNR"),length(loop_order)),
 distance=as.vector(t(as.matrix(s_loop[,c("median_RR","median_R_NR","median_NRNR")]))))
if (cfg$plot$bubble_size=="legacy_rounded_p") {
  # Original bracket layer repeats the rounded annotation on its three segments.
  b$pval <- rep(p_annotation,each=3)
  b$pval_source_cell_type <- rep(loop_order,each=3)
  b$size_val <- -log10(suppressWarnings(as.numeric(b$pval)))
} else {
  b$pval <- rep(s_loop$wilcox_p_NRNR_vs_RR,each=3)
  b$pval_source_cell_type <- b$cell_type
  b$size_val <- -log10(as.numeric(b$pval))
}
write_csv(b,file.path(cfg$paths$output,"bubble_plot_data.csv"))
b$cell_type <- factor(b$cell_type,levels=cfg$cell_type_order)
b$group <- factor(b$group,levels=c("RR","NRR","NRNR"))
g <- ggplot(b,aes(cell_type,group)) + geom_point(aes(fill=distance,size=size_val),shape=21,color="grey20",alpha=0.9) +
 scale_fill_viridis_c(name="Euclidean distance",option="magma") + scale_size(range=c(2,10),name=expression(-log[10](p))) +
 scale_x_discrete(expand=expansion(add=0.6)) + scale_y_discrete(expand=expansion(add=0.6)) +
 coord_fixed() + labs(x=NULL,y=NULL,title="Network similarity across cell types",
 subtitle=if(cfg$plot$bubble_size=="legacy_rounded_p")"Color: median distance; size: rounded NRNR vs RR p-value" else "Color: median distance; size: NRNR vs RR p-value per cell type") +
 theme_minimal(base_size=11) + theme(plot.margin=margin(12,12,12,80),panel.grid=element_blank(),panel.border=element_rect(fill=NA,color="black",size=0.6),axis.ticks=element_blank(),
 axis.text.x=element_text(angle=45,hjust=1,vjust=1,face="bold"),axis.text.y=element_text(face="bold"),plot.title=element_text(face="bold")) +
 guides(size=guide_legend(order=1,override.aes=list(shape=21,fill="grey80")),fill=guide_colorbar(order=2))
ggsave(file.path(cfg$paths$output,"response_revised_bubble_network_similarity.pdf"),g,width=12,height=6)
record_stage(cfg,"05")
