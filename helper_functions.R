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

# Select expressed genes (normalised count > 5 in at least 2 samples)
select_expressed_genes = function(p.norm.gene.counts) {
  sample.names = names(p.norm.gene.counts)[names(p.norm.gene.counts) != "gene_names"]
  
  return(p.norm.gene.counts %>% 
           rowwise() %>%
           mutate(pass_count = sum(c_across(all_of(sample.names)) > expression.cutoff)) %>% 
           ungroup() %>%
           filter(pass_count > pass.count.cutoff) %>%
           pull(gene_names))
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
                rownames_to_column(var = "vicinity_num") %>%
                mutate(vicinity_num = as.integer(vicinity_num)) %>%
                dplyr::select(vicinity_num,
                              region_id),
              by = c("vicinity_num" = "vicinity_num")) %>%
    left_join(norm.region.counts,
              by = c("region_id" = "region_names")) %>%
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

# Sample each region sampling.n times with random coordinates
randomize_region_coords = function(regions.gr, mm10.chr.sizes, sampling.n) {
  chr.names.all = names(mm10.chr.sizes)
  
  chr.names.detected = stringr::str_detect(chr.names.all, "chr[0-9X]+$")
  
  chr.names = chr.names.all[chr.names.detected]
  
  chr.count = length(chr.names)
  
  random.regions.df = bind_rows(lapply(1:(length(regions.gr)),
                                       function(i) {
                                         region.length = end(regions.gr[i]) - start(regions.gr[i]) + 1
                                         region.name = regions.gr[i]$name
                                         return(bind_rows(lapply(1:sampling.n,
                                                                 function(j) {
                                                                   region.chr = sample(chr.names, size = 1)
                                                                   region.start = sample(0:(mm10.chr.sizes[region.chr] - region.length), size = 1)
                                                                   region.end = region.start + region.length - 1
                                                                   return(data.frame(region.chr = region.chr,
                                                                                     region.start = region.start,
                                                                                     region.end = region.end,
                                                                                     region.name = paste0(region.name, "_", j)))
                                                                 })))
                                       }))
  
  return(GRanges(seqnames = random.regions.df$region.chr,
                 ranges = IRanges(start = random.regions.df$region.start,
                                  end = random.regions.df$region.end),
                 strand = rep(strand(regions.gr), sampling.n),
                 name = random.regions.df$region.name,
                 score = rep(regions.gr$score, sampling.n)))
}

# Define a function to calculate empirical p-values for PCCs:
#   
#   ```{r, include=T}
# calc_empirical_pvalue = function(abs.pcc, bkgd.df) {
#   return(nrow(bkgd.df %>% filter(abs_pcc > abs.pcc)) / nrow(bkgd.df))
# }
# ```
# 
# Define a function to find significant NFIA-dependent region-gene assignments:
#   
#   ```{r, include=T}
# assign_regions_to_genes = function(domain.name, area.type) {
#   cat(domain.name, ",", area.type, ": fdr =", fdr, ";", "min.l2fc =", min.l2fc, ";", "min.baseMean =", min.baseMean, "\n")
#   
#   dep.regions = import(paste0("../r_results/select_diff_regions/", domain.name, "_dep_ranges", 
#                               "_fdr", fdr, "_min-l2fc", min.l2fc, "_min-baseMean", min.baseMean,
#                               ".bed"))
#   
#   if (area.type == "radius") {
#     vicinity.file.name = paste0("../r_results/predict_correlated_deg/", domain.name, "_vicinities_radius_", 
#                                 vicinity.radius, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc, 
#                                 "_min-baseMean", min.baseMean, ".bed")
#     
#     vicinity = generate_vicinity_radius(dep.regions, 
#                                         vicinity.radius,
#                                         vicinity.file.name)
#     
#   } else if (area.type == "tad") {
#     vicinity.file.name = paste0("../r_results/predict_correlated_deg/", domain.name, "_vicinities_largest_tad_", 
#                                 min.tad.size, "-", max.tad.size, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                                 "_min-baseMean", min.baseMean, ".bed")
#     
#     vicinity = generate_vicinity_tad(min.tad.size,
#                                      max.tad.size,
#                                      dep.regions,
#                                      domain.name,
#                                      vicinity.file.name)
#   }
#   
#   if (domain.name == "p1") {
#     dep.tss.ranges = p1.dep.tss.ranges
#     
#     norm.region.counts = p1.norm.region.counts
#     
#     norm.gene.counts = p1.norm.gene.counts
#     
#     gene.sample.names = p1.gene.sample.names
#     
#     region.sample.names = p1.region.sample.names
#     
#   } else if (domain.name == "p2") {
#     dep.tss.ranges = p2.dep.tss.ranges
#     
#     norm.region.counts = p2.norm.region.counts
#     
#     norm.gene.counts = p2.norm.gene.counts
#     
#     gene.sample.names = p2.gene.sample.names
#     
#     region.sample.names = p2.region.sample.names
#     
#   } else if (domain.name == "pM") {
#     dep.tss.ranges = pM.dep.tss.ranges
#     
#     norm.region.counts = pM.norm.region.counts
#     
#     norm.gene.counts = pM.norm.gene.counts
#     
#     gene.sample.names = pM.gene.sample.names
#     
#     region.sample.names = pM.region.sample.names
#   }
#   
#   if (area.type == "radius") {
#     corr.filename = paste0("../r_results/predict_correlated_deg/", domain.name, "_corr_in_radius_",
#                            vicinity.radius, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                            "_min-baseMean", min.baseMean, ".rds")
#     
#     corr.plot.filename = paste0("../r_results/predict_correlated_deg/plots/", domain.name, "_corr_in_radius_",
#                                 vicinity.radius, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                                 "_min-baseMean", min.baseMean, "_hist.pdf")
#     
#     bkgd.filename = paste0("../r_results/predict_correlated_deg/", domain.name, "_bkgd_in_radius_",
#                            vicinity.radius, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                            "_min-baseMean", min.baseMean, ".rds")
#     
#     bkgd.plot.filename = paste0("../r_results/predict_correlated_deg/plots/", domain.name, "_bkgd_in_radius_",
#                                 vicinity.radius, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                                 "_min-baseMean", min.baseMean, "_hist.pdf")
#     
#   } else if (area.type == "tad") {
#     corr.filename = paste0("../r_results/predict_correlated_deg/", domain.name, "_corr_in_largest_tad_",
#                            min.tad.size, "-", max.tad.size, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                            "_min-baseMean", min.baseMean, ".rds")
#     
#     corr.plot.filename = paste0("../r_results/predict_correlated_deg/plots/", domain.name, "_corr_in_largest_tad_",
#                                 min.tad.size, "-", max.tad.size, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                                 "_min-baseMean", min.baseMean, "_hist.pdf")
#     
#     bkgd.filename = paste0("../r_results/predict_correlated_deg/", domain.name, "_bkgd_in_largest_tad_",
#                            min.tad.size, "-", max.tad.size, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                            "_min-baseMean", min.baseMean, ".rds")
#     
#     bkgd.plot.filename = paste0("../r_results/predict_correlated_deg/plots/", domain.name, "_bkgd_in_largest_tad_",
#                                 min.tad.size, "-", max.tad.size, "kbp_dep_regions_fdr", fdr, "_min-l2fc", min.l2fc,
#                                 "_min-baseMean", min.baseMean, "_hist.pdf")
#   }
#   
#   corr = calc_correlations(vicinity, dep.tss.ranges,
#                            norm.region.counts, norm.gene.counts,
#                            region.sample.names, gene.sample.names,
#                            domain.name,
#                            "target",
#                            corr.filename)
#   
#   p = ggplot(corr) +
#     geom_histogram(aes(x = abs_pcc), 
#                    binwidth = 0.05) +
#     theme_classic()
#   
#   ggsave(filename = corr.plot.filename,
#          plot = p)
#   
#   bkgd = calc_correlations(vicinity, dep.tss.ranges,
#                            norm.region.counts, norm.gene.counts,
#                            region.sample.names, gene.sample.names,
#                            domain.name,
#                            "background",
#                            bkgd.filename)
#   
#   p = ggplot(bkgd) +
#     geom_histogram(aes(x = abs_pcc), 
#                    binwidth = 0.05) +
#     theme_classic()
#   
#   ggsave(filename = bkgd.plot.filename,
#          plot = p)
#   
#   corr.sign = corr %>%
#     rowwise() %>%
#     mutate(pvalue = calc_empirical_pvalue(abs_pcc, bkgd)) %>%
#     ungroup() %>%
#     mutate(padj = p.adjust(pvalue, method = "BH")) %>%
#     mutate(fdr = pcc.fdr) %>%
#     mutate(sign.padj = (padj < fdr))
#   
#   cat(domain.name, ":", area.type, "\n")
#   
#   cat("Number of significant NFIA-dependent region-gene pairs ( FDR =", pcc.fdr, ") :", nrow(corr.sign %>% filter(sign.padj)), "\n")
#   
#   cat("---\n")
#   
#   if (nrow(corr.sign %>% filter(sign.padj)) > 0) {
#     if (domain.name == "p1") {
#       strictest.dep.regions = p1.dep.regions$name
#       
#     } else if (domain.name == "p2") {
#       strictest.dep.regions = p2.dep.regions$name
#       
#     } else if (domain.name == "pM") {
#       strictest.dep.regions = pM.dep.regions$name
#     }
#     
#     sign.assigned.regions = corr.sign %>% 
#       filter(sign.padj) %>%
#       pull(region_id) %>%
#       unique()
#     
#     cat("Number of the strictest regions among all regions in significant assignments:", 
#         length(intersect(sign.assigned.regions, strictest.dep.regions)),
#         "(", length(intersect(sign.assigned.regions, strictest.dep.regions)) / length(sign.assigned.regions), ")\n")
#   }
#   
#   return(list("vicinity" = vicinity,
#               "corr" = corr,
#               "bkgd" = bkgd,
#               "corr.sign" = corr.sign))
# }
# ```