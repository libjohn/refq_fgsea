# README

<!-- badges: start -->

<!-- badges: end -->

Using tidyverse principles, demonstrate iteration of *fast gene enrichment* as seen in the {[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html)} package.

This annotated reproducible code uses the example data from the {fgsea} ([vignette](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)) package after considering a biosquid [fgsea-tutorial](https://biostatsquid.com/fgsea-tutorial-gsea/). The main goal of this code is to share a reproducible example with accessible data. Code for an applied example is also shown. The data for the applied example is **not shown** because I am unaware of the data's provenance.

## Basic steps

1.  Import gmt file as a data frame
    1.  Pivot to tall so each pathway has a vector of gene symbols
    2.  `drop_na`() on gene_symbols
    3.  `group_by()` pathways ; summarize the gene symbols as a list
    4.  Add a list column as a named-vector of rankings from a separate computation
    5.  iterate rowwise (via `map2()` ) to create plots via `fgsea::plotEnrichment()`

```{my_gmt_file <- read_table("data/my_symbols.gmt",}
                                       col_names = FALSE) |> 
  janitor::clean_names() |> 
  select(-x2) 

my_df <- read_csv("data/my_model_ranking.csv", show_col_types = FALSE) |> 
  janitor::clean_names() |> 
  rename(gene_symbol = x1, rank_avg_log = avg_log2fc) 

my_rank_vector <- my_df |> 
  drop_na(rank_avg_log) |>
  arrange(rank_avg_log) |> 
  select(gene_symbol, rank_avg_log) |> 
  summarise(my_ranked_vector = list(setNames(rank_avg_log, gene_symbol))) |> 
  pull(my_ranked_vector) |> 
  unlist()

my_rank_vector |> head()

my_enriched_plots_df <- my_gmt_file |> 
  rename(pathways = x1) |> 
  pivot_longer(x3:last_col(), values_to = "genes_in_pathway") |> 
  drop_na(genes_in_pathway) |> 
  group_by(pathways) |> 
  summarise(gene_symbols = list(genes_in_pathway)) |> 
  mutate(my_ranks = list(my_rank_vector)) |> 
  mutate(my_plot = map2(gene_symbols, my_ranks, \(x, y) 
                        plotEnrichment(x, y)))

my_enriched_plots_df

my_enriched_plots_df |> 
  pull(my_plot)
```
