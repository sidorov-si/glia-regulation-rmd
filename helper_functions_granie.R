generate_sample_meta = function(p_gene_expression,
                                p_region_accessibility) {
  return(tibble(sample_id = c(unlist(stringr::str_match_all(c(names(p_gene_expression), 
                                                              names(p_region_accessibility)),
                                                            "^WT_.+")))) %>%
     mutate(condition   = unlist(stringr::str_match_all(sample_id, "^..|^...")),
            day         = unlist(stringr::str_match_all(sample_id, "D[0-9]+")),
            nfia_status = unlist(stringr::str_match_all(sample_id, "NFIAp|NFIAn")),
            rep_id      = unlist(stringr::str_match_all(sample_id, "R[0-9]+$")),
            domain      = unlist(stringr::str_match_all(sample_id, "p[12M]"))))
}

