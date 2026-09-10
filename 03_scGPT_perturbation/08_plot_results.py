"""Generate the four Python manuscript plot types from computed summary tables.

Input: SRC top-20/filtered rank XLSX and kNN k=100 filtered metrics XLSX.
Output: figures/ ranked SPC, ranked kNN, top-20 barplot, gene-distance-change histogram PDF/PNG and rank summaries.
Plot functions, rank ordering, ties and empirical p-values are retained from source.
Run all manifest perturbations before interpreting background ranks.
"""
from pipeline.common import cli, load_config, output_dir, input_path


def main():
    args=cli(__doc__).parse_args()
    cfg=load_config(args.config)
    from pipeline import figure_plots as plots, top20_plot
    out=output_dir(cfg,'figures')
    src=output_dir(cfg,'SRC')
    knn=output_dir(cfg,'kNN')
    k=cfg['knn']['k']
    specs=[
        plots.PlotSpec(src/f"spc_mean_matrix_true_count_ge{cfg['src']['min_original_nonzero']}.xlsx",
            'top20_mean_spc_euclidean_mean','selected_genes_spc_top20_mean_euclidean_mean',
            'Selected genes in SPC mean matrix','Top20 mean SPC euclidean mean'),
        plots.PlotSpec(knn/f'knn_label_enrichment_filtered_metrics_no_gate_nr_genezero_paired_k{k}.xlsx',
            'knn_nr_fraction_shift_mean_euclidean','selected_genes_knn_shift_mean_euclidean',
            'Selected genes in kNN enrichment','kNN NR fraction shift mean (euclidean)')]
    rows=[]
    for spec in specs:
        work, targets=plots.load_metric_table(spec.input_path,spec.metric,cfg['figure5']['targets'])
        summary=plots.compute_summaries(work,targets,spec.metric)
        plots.save_rank_plot(work,targets,spec.metric,spec,out,summary)
        for row in summary.to_dict('records'):
            rows.append({'source_file':spec.input_path.name,'metric':spec.metric,**row})
    plots.save_summary_table(rows,out)
    top20_plot.INPUT_XLSX=src/'spc_euclidean_top20_selected_genes.xlsx'
    top20_plot.OUTPUT_PDF=out/'spc_euclidean_top20_selected_genes_barplot_custom.pdf'
    top20_plot.OUTPUT_PNG=out/'spc_euclidean_top20_selected_genes_barplot_custom.png'
    top20_plot.main()
    from pipeline.distance_change_plot import plot_distance_change
    plot_distance_change(input_path(cfg, "distance_change_scores"), out)


if __name__=='__main__':
    main()
