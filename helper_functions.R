# DESeq2 count normalisation across replicates
normalize_counts = function(raw.counts, cond.table, diff.levels, dds.filename) {
  raw.counts = round(raw.counts)
  
  row.names(cond.table) = names(raw.counts)
  
  cond.table$condition = factor(cond.table$condition, levels = diff.levels)
  
  dds = DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cond.table,
                               design = ~ condition)
  
  dds = DESeq(dds)
  
  saveRDS(dds, dds.filename)
  
  dds = estimateSizeFactors(dds)
  
  return(counts(dds, normalized = T))
}

# Wrapper function for DESeq2 expression normalisation
generate_norm_gene_counts = function(master.table, domain.name) {
  p.gene.sample.names = unlist(stringr::str_match_all(names(master.table), paste0("WT_D[179]+_", domain.name, ".*")))
  
  p.raw.counts = master.table %>%
    dplyr::select(all_of(p.gene.sample.names))
  
  p.cond = data.frame(condition = stringr::str_extract(names(p.raw.counts), "D(7|9|11)"))
  
  p.diff.levels = c("D7", "D9", "D11")
  
  p.norm.gene.counts = as.data.frame(normalize_counts(p.raw.counts, 
                                                      p.cond, 
                                                      p.diff.levels,
                                                      paste0("../r_results/predict_correlated_expressed_gene/",
                                                             domain.name, "_gene_dds.rds"))) %>%
    mutate(gene_names = gene.names) %>%
    dplyr::select(gene_names, all_of(names(.)[names(.) != "gene_names"]))
  
  saveRDS(p.norm.gene.counts,
          file = paste0("../r_results/predict_correlated_expressed_gene/",
                        domain.name, "_norm_counts_expression.rds"))
  
  return(p.norm.gene.counts)
}

# Calculate z-scores per feature across samples
generate_scaled_counts = function(count.df, annot.col.name) {
  rownames(count.df) = NULL
  
  count.df.clean = count.df %>%
    column_to_rownames(var = annot.col.name)
  
  count.df.clean.scaled = as.data.frame(t(apply(count.df.clean, 1, scale)))
  
  names(count.df.clean.scaled) = names(count.df.clean)
  
  count.df.clean.final = count.df.clean.scaled %>%
    rownames_to_column(var = annot.col.name)
  
  return(count.df.clean.final)
}

# Wrapper function for DESeq2 accessibility normalisation
generate_norm_chrom_counts = function(master.table, region.annot, domain.name) {
  p.region.sample.names = unlist(stringr::str_match_all(names(master.table), paste0("WT_D[179]+_", domain.name, ".*")))
  
  p.raw.counts = master.table %>%
    dplyr::select(all_of(p.region.sample.names))
  
  p.cond = data.frame(condition = stringr::str_extract(names(p.raw.counts), "D(7|9|11)"))
  
  p.diff.levels = c("D7", "D9", "D11")
  
  p.norm.region.counts = as.data.frame(normalize_counts(p.raw.counts, 
                                                        p.cond, 
                                                        p.diff.levels,
                                                        paste0("../r_results/predict_correlated_expressed_gene/", 
                                                               domain.name, "_region_dds.rds"))) %>%
    mutate(region_names = region.annot$Geneid) %>%
    dplyr::select(region_names, all_of(names(.)[names(.) != "region_names"]))
  
  saveRDS(p.norm.region.counts,
          file = paste0("../r_results/predict_correlated_expressed_gene/",
                        domain.name, "_norm_counts_accessibility.rds"))
  
  return(p.norm.region.counts)
}

# Select expressed genes
select_expressed_genes = function(p.norm.gene.counts, expression.median) {
  sample.names = names(p.norm.gene.counts)[names(p.norm.gene.counts) != "gene_names"]
  
  return(p.norm.gene.counts %>% 
           rowwise() %>%
           mutate(median_expression = median(c_across(all_of(sample.names)))) %>% 
           ungroup() %>%
           filter(median_expression > expression.median) %>%
           pull(gene_names))
}

# Select genes by expression amplitude
select_genes_by_amplitude = function(p.norm.gene.counts, expression.amplitude.remove) {
  sample.names = names(p.norm.gene.counts)[names(p.norm.gene.counts) != "gene_names"]
  
  p.norm.gene.counts %<>% 
    rowwise() %>%
    mutate(min_z = min(c_across(all_of(sample.names)))) %>% 
    mutate(max_z = max(c_across(all_of(sample.names)))) %>% 
    mutate(ampl = max_z - min_z) %>% 
    ungroup() %>%
    arrange(-ampl)

  cutoff.rownumber = round((1 - expression.amplitude.remove) * nrow(p.norm.gene.counts))

  return (p.norm.gene.counts %>%
            head(cutoff.rownumber))
}

# Generate TSS coordinates from transcript coordinates
generate_tss = function(tx.ranges, tss.bed.filename) {
  tss.ranges = GRanges(seqnames = seqnames(tx.ranges),
                       ranges = IRanges(start = ifelse(strand(tx.ranges) == "+",
                                                       start(ranges(tx.ranges)),
                                                       end(ranges(tx.ranges)) - 1),
                                        end = ifelse(strand(tx.ranges) == "+",
                                                     start(ranges(tx.ranges)) + 1,
                                                     end(ranges(tx.ranges))),
                                        names = names(ranges(tx.ranges))),
                       strand = strand(tx.ranges),
                       gene_name = tx.ranges$gene_name)

  export(tss.ranges, tss.bed.filename)

  return(tss.ranges)
}

# Generate a TSS annotation
generate_tss_annot = function(p.genes, mm10.annot.genes, p.tss.filename) {
  p.genes.annot = mm10.annot.genes[mm10.annot.genes$gene_name %in% p.genes, ]

  p.genes.ranges = GRanges(seqnames = p.genes.annot$chr,
                           ranges = IRanges(start = p.genes.annot$start,
                                            end   = p.genes.annot$stop,
                                            names = p.genes.annot$tx_name),
                           strand = strand(p.genes.annot$strand),
                           gene_name = p.genes.annot$gene_name)

  p.tss.ranges = generate_tss(p.genes.ranges, p.tss.filename)

  return(p.tss.ranges)
}

# Calculate and plot the distribution of distance between NFIA-dependent region and the closest TSS (strand ignored)
calc_dist_distributions = function(dep.regions, expr.tss.ranges, domain.name) {
  dist.hits.obj = distanceToNearest(dep.regions, expr.tss.ranges, ignore.strand = T)
  
  dist.vector = mcols(dist.hits.obj)$distance
  
  cat("The number of NFIA-dependent regions in promoters (<= 1 kbp from a TSS; strand ignored):", 
      length(dist.vector[dist.vector <= 1000]),
      "(", round(length(dist.vector[dist.vector <= 1000]) / length(dist.vector) * 100, 2), "% )\n")
}

# Generate the region-gene assignment area around a region as a vicinity of a fixed radius (in kbp)
generate_vicinity_radius = function(dep.regions, vicinity.radius, vicinity.ranges.bed.filename) {
  vicinity.radius = vicinity.radius * 1000

  dep.regions.df = as.data.frame(dep.regions) %>%
    mutate(seqnames = as.character(seqnames)) %>%
    mutate(strand = as.character(strand)) %>%
    rowwise() %>%
    mutate(., max_position = mm10.chr.sizes[seqnames] - 1) %>%
    ungroup()

  vicinity.ranges = GRanges(seqnames = dep.regions.df$seqnames,
                            ranges = IRanges(start = ifelse(dep.regions.df$start - vicinity.radius > 0,
                                                            dep.regions.df$start - vicinity.radius,
                                                            1),
                                             end = ifelse(dep.regions.df$end + vicinity.radius < dep.regions.df$max_position,
                                                          dep.regions.df$end + vicinity.radius,
                                                          dep.regions.df$max_position),
                                             names = dep.regions.df$name),
                            strand = Rle(strand("*")))

  export(vicinity.ranges, vicinity.ranges.bed.filename)

  return(vicinity.ranges)
}

# Add a "_gene" suffix
add_gene_suffix = function(column.names) {
  return(paste0(column.names, "_gene"))
}

# Add a "_region" suffix
add_region_suffix = function(column.names) {
  return(paste0(column.names, "_region"))
}

# Find Pearson correlation coefficients between region accessibility and gene expression
calc_correlations = function(vicinity.gr, dep.tss,
                             norm.region.counts, norm.gene.counts,
                             region.sample.names, gene.sample.names,
                             domain.name, corr.rds.filename) {
  region.tss.shared.sample.names = intersect(region.sample.names, gene.sample.names)
  
  tss.vicinity.hits = GenomicRanges::findOverlaps(query = dep.tss,
                                                  subject = vicinity.gr,
                                                  type = "within",
                                                  select = "all",
                                                  ignore.strand = T)

  tss.hits = data.frame(tss_num = queryHits(tss.vicinity.hits)) %>%
    left_join(as.data.frame(dep.tss) %>%
                rownames_to_column(var = "tss_id") %>%
                rownames_to_column(var = "tss_num") %>%
                mutate(tss_num = as.integer(tss_num)) %>%
                dplyr::select(tss_num,
                              gene_name),
              by = c("tss_num" = "tss_num")) %>%
    left_join(norm.gene.counts,
              by = c("gene_name" = "gene_names")) %>%
    dplyr::select(tss_num,
                  gene_name,
                  all_of(sort(names(.)[names(.) %in% region.tss.shared.sample.names]))) %>%
    dplyr::rename_with(add_gene_suffix, all_of(region.tss.shared.sample.names))

  vicinity.hits = data.frame(vicinity_num = subjectHits(tss.vicinity.hits)) %>%
    left_join(as.data.frame(vicinity.gr) %>%
                rownames_to_column(var = "region_id") %>%
                mutate(region_names = unlist(stringr::str_replace(region_id, "\\.[1-9]+$", ""))) %>%
                rownames_to_column(var = "vicinity_num") %>%
                mutate(vicinity_num = as.integer(vicinity_num)) %>%
                dplyr::select(vicinity_num,
                              region_names,
                              region_id),
              by = c("vicinity_num" = "vicinity_num")) %>%
    left_join(norm.region.counts,
              by = c("region_names" = "region_names")) %>%
    dplyr::select(vicinity_num,
                  region_id,
                  all_of(sort(names(.)[names(.) %in% region.tss.shared.sample.names]))) %>%
    dplyr::rename_with(add_region_suffix, all_of(region.tss.shared.sample.names))

  vicinity.tss.corr = vicinity.hits %>%
    bind_cols(tss.hits) %>%
    dplyr::select(-vicinity_num,
                  -tss_num) %>%
    distinct() %>%
    rowwise() %>%
    mutate(pcc = cor(x = c_across(ends_with("_region")),
                     y = c_across(ends_with("_gene")),
                     method = "pearson")) %>%
    ungroup()

  vicinity.tss.corr.final = vicinity.tss.corr %>%
    dplyr::select(region_id,
                  gene_name,
                  pcc) %>%
    mutate(abs_pcc = abs(pcc))

  cat(domain.name, ":\n")

  cat("Number of region-gene associations:", nrow(vicinity.tss.corr.final), "\n")

  cat("Number of unique regions          :", length(unique(vicinity.tss.corr.final$region_id)), "\n")

  cat("Number of unique genes            :", length(unique(vicinity.tss.corr.final$gene_name)), "\n")

  cat("---\n")

  saveRDS(vicinity.tss.corr.final, corr.rds.filename)

  return(vicinity.tss.corr.final)
}

# Calc all pairwise correlations between regions and genes
calc_correlations_pairs_all = function(p.dep.regions, 
                                       p.norm.region.counts, p.norm.gene.counts,
                                       p.region.sample.names, p.gene.sample.names,
                                       domain.name, corr.rds.filename) {
  p.names = intersect(p.region.sample.names, p.gene.sample.names)
  
  p.norm.region.counts.dep = p.norm.region.counts %>%
    filter(region_names %in% p.dep.regions$name) %>%
    dplyr::select(all_of(c("region_names", p.names))) %>%
    column_to_rownames(var = "region_names")
  
  p.norm.gene.counts.shared = p.norm.gene.counts %>%
    dplyr::select(all_of(c("gene_names", p.names))) %>%
    column_to_rownames(var = "gene_names")
  
  p.bkgd.radius.df = as.data.frame(cor(x = t(p.norm.region.counts.dep), y = t(p.norm.gene.counts.shared))) %>%
    rownames_to_column(var = "region_names") %>%
    gather("gene_names", "pcc", -region_names) %>%
    dplyr::rename("region_id" = "region_names",
                  "gene_name" = "gene_names") %>%
    mutate(abs_pcc = abs(pcc)) %>%
    tibble()
  
  cat(domain.name, ":\n")
  
  cat("Number of region-gene associations:", nrow(p.bkgd.radius.df), "\n")
  
  cat("Number of unique regions          :", length(unique(p.bkgd.radius.df$region_id)), "\n")
  
  cat("Number of unique genes            :", length(unique(p.bkgd.radius.df$gene_name)), "\n")
  
  cat("---\n")
  
  saveRDS(p.bkgd.radius.df, corr.rds.filename)
  
  return(p.bkgd.radius.df)
}

# Shuffle region names across vicinities
# shuffle_regions_across_vicinities = function(p.vicinity.radius, output.filename) {
#   names(p.vicinity.radius) = sample(names(p.vicinity.radius))
#   
#   export(p.vicinity.radius, output.filename)
#   
#   return(p.vicinity.radius)
# }

# Calculate empirical p-values for target PCCs
calc_empirical_pvalue = function(abs.pcc, abs_pcc_sorted, bkgd.df_nrow) {
  return(sum(abs_pcc_sorted >= abs.pcc) / bkgd.df_nrow)
}

# Calculate and plot pairwise PCCs between NFIA-dependent regions
calc_and_plot_dep_region_pccs = function(dep.regions,
                                         norm.region.counts,
                                         domain.name) {
  dep.region.norm.counts.transp = norm.region.counts %>% 
    filter(region_names %in% dep.regions$name) %>%
    column_to_rownames(var = "region_names") %>% 
    t()
  
  norm.region.counts.transp.cor = cor(dep.region.norm.counts.transp, method = "pearson")
  
  dep.region.cor.plot = data.frame(region_pcc = norm.region.counts.transp.cor[upper.tri(norm.region.counts.transp.cor)]) %>%
    ggplot(aes(x = region_pcc)) +
    geom_density(fill = "grey") +
    scale_x_continuous(limits = c(-1, 1)) + #,
                       #breaks = breaks_pretty(n = 20)) +
    theme_classic()
  
  ggsave(filename = paste0("../r_results/predict_correlated_expressed_gene/plots/", 
                           domain.name, "_dep_regions_corr",
                           "_regions_fdr", fdr, "_min-l2fc", min.l2fc,
                           "_min-baseMean", min.baseMean, "_density.pdf"),
         plot = dep.region.cor.plot)
  
  cat("The number of pairwise PCCs between NFIA-dependent regions in", domain.name, ":", 
      length(norm.region.counts.transp.cor[upper.tri(norm.region.counts.transp.cor)]), "\n")
}

# Count expressed genes whose TSS(s) are inside a radius-defined vicinity of a least one NFIA-dependent region
count_expr_genes_inside_vicinity = function(dep.regions.vicinity.radius, expr.tss.ranges) {
  overlaps.hits.obj = findOverlaps(expr.tss.ranges, 
                                   dep.regions.vicinity.radius,
                                   type = "within",
                                   select = "all",
                                   ignore.strand = T)
  
  expr.genes.inside = unique(expr.tss.ranges$gene_name[unique(queryHits(overlaps.hits.obj))])
  
  return(length(expr.genes.inside))
}

find_expr_genes_inside_vicinity = function(dep.regions.vicinity.radius, expr.tss.ranges) {
  overlaps.hits.obj = findOverlaps(expr.tss.ranges, 
                                   dep.regions.vicinity.radius,
                                   type = "within",
                                   select = "all",
                                   ignore.strand = T)
  
  expr.genes.inside = unique(expr.tss.ranges$gene_name[unique(queryHits(overlaps.hits.obj))])
  
  return(expr.genes.inside)
}

# Make sampling.n copies of each vicinity in a GRanges object p1.vicinity.radius
# multiply_vicinities = function(vicinity.gr, sampling.n, output.filename) {
#   vicinity.gr.mult = unlist(as(lapply(1:sampling.n,
#                                       function(i) {
#                                         names(vicinity.gr) = paste0(names(vicinity.gr), ".", i)
#                                         return(vicinity.gr)
#                                       }), 
#                                "GRangesList"))
#   
#   export(vicinity.gr.mult, output.filename)
#   
#   return(vicinity.gr.mult)
# }

# Re-format a vicinity GRanges read from a BED file
trim_vicinity_gr = function(vicinity.gr) {
  names(vicinity.gr) = vicinity.gr$name
  
  vicinity.gr$name = NULL
  
  vicinity.gr$score = NULL
  
  return(vicinity.gr)
}
