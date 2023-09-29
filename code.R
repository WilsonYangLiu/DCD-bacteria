#
setwd(dir = '/storage1/wilson/Amplicon/DCD-mapping-16Slong_Wilson')

#
require(ggplot2)
# .libPaths('/Library/Frameworks/R.framework/Versions/library-3.5')
require(microbiome)
require(ComplexHeatmap)
require(circlize)
require(dendextend)

mapping <- c('k', 'p', 'c', 'o', 'f', 'g')
names(mapping) = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')

get_lv <- function(x, usedLv, mapping = c('k', 'p', 'c', 'o', 'f', 'g', 's')) {
  if (usedLv == 'k') {
    return(x[usedLv])
  }
  if (x[usedLv] == '' || endsWith(x = x[usedLv], suffix = '__')) {
    usedLv <- mapping[which(mapping %in% usedLv) - 1]
    return(get_lv(x, usedLv))
  } 
  return(x[usedLv])
}

taxa_chr2df_sp.silva <- function(x) {
  mapping <- rep('', 7)
  names(mapping) <- c('D_0', 'D_1', 'D_2', 'D_3', 'D_4', 'D_5', 'D_6')
  taxa_raw <- sapply(X = x, FUN = function(t) {
    t_new <- rep('', 7)
    t <- stringi::stri_split(str = t, regex = '\\s*[;|]\\s*')[[1]]
    level <- gsub(pattern = '^(D_.)__.*$', replacement = '\\1', x = t)
    idx <- match(level, names(mapping))
    if (all(is.na(idx))) {
      return(t_new)
    } else {
      t_new[idx[!is.na(idx)]] <- t[!is.na(idx)]
      t_new = gsub('^D_[0-9]__', '', t_new)
      return(t_new)
    }
  })
  taxa_raw <- t(taxa_raw)
  colnames(taxa_raw) <- names(mapping)
  
  return(taxa_raw)
}

taxa_chr2df_sp <- function(x) {
  mapping <- rep('', 7)
  names(mapping) <- c('k', 'p', 'c', 'o', 'f', 'g', 's')
  taxa_raw <- sapply(X = x, FUN = function(t) {
    t_new <- rep('', 7)
    t <- stringi::stri_split(str = t, regex = '\\s*[;|]\\s*')[[1]]
    level <- gsub(pattern = '^(.)__.*$', replacement = '\\1', x = t)
    idx <- match(level, names(mapping))
    if (all(is.na(idx))) {
      return(t_new)
    } else {
      t_new[idx[!is.na(idx)]] <- t[!is.na(idx)]
      return(t_new)
    }
  })
  taxa_raw <- t(taxa_raw)
  colnames(taxa_raw) <- names(mapping)
  
  return(taxa_raw)
}


# -0 Load the data
if (FALSE) {
  load(file = 'Result/Saved-tab.RData') # tax, tax_wSpecies, seqtab.nochim
  
  df_tax <- as.data.frame(tax, stringsAsFactors = FALSE);rm(tax)
  table(is.na(df_tax$Species))
  # FALSE  TRUE 
  # 1648 78045
  
  df_tax_wSpecies <- as.data.frame(tax_wSpecies, stringsAsFactors = FALSE);rm(tax_wSpecies) # <<<------ use this
  table(is.na(df_tax_wSpecies$Species))
  # FALSE  TRUE 
  # 28813 50880
  
  if (FALSE) {
    #
    df_tax_blast <- read.csv(file = gzfile('Result/df_tax.blast_silva_132_99_16S.tsv.gz'), 
                             header = FALSE, sep = '\t', check.names = FALSE, comment.char = "#",
                             strip.white = TRUE, stringsAsFactors = FALSE)
    colnames(df_tax_blast) = c('query', 'subject', '% identity', 'alignment length', 'mismatches', 'gap opens', 
                               'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score')
    
    #
    query_seq = seqinr::read.fasta(file = 'Result/rep_seq.fa')
    df_tax_blast$`q. len` <- sapply(df_tax_blast$query, function(x) length(query_seq[[x]]))
    
    tmp = sapply(query_seq, FUN = function(x) toupper(seqinr::getSequence(x, as.string = TRUE)[[1]]))
    df_tax_blast$`q. seq` <- tmp[df_tax_blast$query]
    
    df_tax_blast$`q. taxonomy` = apply(df_tax_wSpecies[df_tax_blast$`q. seq`, ], 1, function(x) {
      x[is.na(x)] = ''
      return(paste(x, collapse = ';'))
    })
    
    ## -
    subject_tax = read.csv(file = 'Result/SILVA_132_QIIME_release/taxonomy_7_levels.txt', 
                           header = FALSE, sep = '\t', check.names = FALSE, comment.char = "#", 
                           strip.white = TRUE, stringsAsFactors = FALSE)
    colnames(subject_tax) = c('ID', 'taxonomy')
    rownames(subject_tax) = subject_tax$ID
    
    tmp = as.data.frame(taxa_chr2df_sp.silva(subject_tax$taxonomy), stringsAsFactors = FALSE)
    rownames(tmp) = subject_tax$ID
    colnames(tmp) = colnames(df_tax_wSpecies)
    subject_tax = tmp; rm(tmp)
    
    df_tax_blast$`s. taxonomy` = apply(subject_tax[df_tax_blast$subject, ], 1, paste, collapse = ';')
    # df_tax_blast = cbind(df_tax_blast, tmp[df_tax_blast$subject, ])
    
    #
    # ratio of aligned part of query sequence
    df_tax_blast$Ratio <- as.numeric(df_tax_blast$`alignment length`) / as.numeric(df_tax_blast$`q. len`)
    table(df_tax_blast$Ratio > 0.99)# FALSE:TRUE, 516511:267468 
    table(df_tax_blast$Ratio > 0.95)# FALSE:TRUE,  238950:545029
    table(df_tax_blast$Ratio > 0.9)# FALSE:TRUE,  88187:695792 
    
    write.csv(x = df_tax_blast, file = gzfile('Result/df_tax.blast_silva_132_99_16S.curated.csv.gz'))
  }
  # save(df_tax_blast, file = 'Result/df_tax.blast_silva_132_99_16S.curated.RData')
  load(file = 'Result/df_tax.blast_silva_132_99_16S.curated.RData') # df_tax_blast
  
  
  ##
  query_seq = seqinr::read.fasta(file = 'Result/rep_seq.fa')
  tmp = sapply(query_seq, FUN = function(x) toupper(seqinr::getSequence(x, as.string = TRUE)[[1]]))
  query_seq = names(tmp); names(query_seq) = tmp
  
  rownames(df_tax) = unname(query_seq[rownames(df_tax)])
  rownames(df_tax_wSpecies) = unname(query_seq[rownames(df_tax_wSpecies)])
  colnames(seqtab.nochim) = unname(query_seq[colnames(seqtab.nochim)])
  
  
  # -
  metadat = read.csv(file = 'metadata_dcd_PacBio.tsv', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, 
                     strip.white = TRUE, stringsAsFactors = FALSE)
  add_info = read.csv(file = 'Result/summary.csv', header = TRUE, row.names = 1, check.names = FALSE, 
                      strip.white = TRUE, stringsAsFactors = FALSE)
  rownames(add_info) = gsub('-filtered.fq.gz', '', rownames(add_info))
  
  metadat = cbind(sample = rownames(metadat), add_info[rownames(metadat), ], metadat); rm(add_info)
  
  ## 
  microbe_count <- as.data.frame(t(seqtab.nochim)); rm(seqtab.nochim)
  colnames(microbe_count) <- gsub('-filtered.fq.gz', '', colnames(microbe_count))
  microbe_count <- microbe_count[, rownames(metadat)]
  
  microbe_feature = df_tax_wSpecies # 
  all(rownames(microbe_feature) == rownames(microbe_count))
  
  # -
  save(df_tax, df_tax_blast, df_tax_wSpecies, microbe_count, microbe_feature, metadat, 
       file = 'Result/00.Saved.input.RData')
}

# -0 collect sample type: sample_type.csv
if (FALSE) { 
  # df_tax, df_tax_blast, df_tax_wSpecies, microbe_count, microbe_feature, metadat
  load(file = 'Result/00.Saved.input.RData') 
  
  map_site_orgam = c('NEG', 
                     'SKINLP', 'SKINUP', 
                     'ORALmix', 
                     'ESOM', 'ZA1', 'Z', 'ZB1', 
                     'SF', 'SB', 'GJ', 'SA', 'PY', 
                     'DUS', 'DUPA', 'TRE', 'JEJ100', 'JEJ200', 'IIE300', 'IIE400', 'IICP', 'IIC', 
                     'AP', 
                     'TOC', 'ASCM', 'TRC', 'DECM', 'SICM', 'RECU', 'ANAL')
  
  metadat_used <- metadat
  metadat_used$sampletype[metadat_used$sampletype == 'U'] = 'F'
  sample_type <- matrix(data = 0, nrow = length(unique(metadat_used$site)), ncol = 3, 
                        dimnames = list(map_site_orgam, c('M', 'S', 'F')))
  tmp <- apply(as.data.frame(table(paste(metadat_used$site, metadat_used$sampletype, sep = ':'))), 1, function(item) {
    tmp <- stringi::stri_split(item[1], regex = ':')[[1]]
    sample_type[tmp[1], tmp[2]] <<- as.numeric(item[2])
    
    return(NA)
  })
  
  sample_loc = unique(metadat_used[, c("site", "MajorBS")])
  rownames(sample_loc) = sample_loc$site
  
  sample_type <- cbind(sample_type, sample_loc[rownames(sample_type), ])
  write.csv(sample_type, file = 'Result/sample_type.csv')
}

# -0 negative control
if (FALSE) {
  # df_tax, df_tax_blast, df_tax_wSpecies, microbe_count, microbe_feature, metadat
  load(file = 'Result/00.Saved.input.RData') 
  
  # load and set the data: metadat
  if (TRUE) {
    metadat$Antibiotics[metadat$subject == 'S006'] <- 'Y'
    
    map_site_orgam <- c(
      'NEG' = 'Negtive control'
      ,
      "SKINLP" = 'Body skin',
      "SKINUP" = 'Body skin'
      ,
      'ORALmix' = 'Oral cavity'
      ,
      'ESOM' = 'Esophageal',
      'ZA1' = 'Esophageal',
      'Z' = 'Esophageal',
      'ZB1' = 'Esophageal'
      ,
      'SF' = 'Stomach',
      'SB' = 'Stomach',
      "GJ" = 'Stomach',
      "SA" = 'Stomach',
      "PY" = 'Stomach'
      ,
      "DUS" = 'SI',
      'DUPA' = 'SI',
      'TRE' = 'SI',
      "JEJ100" = 'SI',
      "JEJ200" = 'SI',
      "IIE300" = 'SI',
      "IIE400" = 'SI',
      "IICP" = 'SI',
      "IIC" = 'SI'
      ,
      'AP' = 'Appendix'
      ,
      'Cecum' = 'LI',
      "ASCM" = 'LI',
      "TRC" = 'LI',
      "DECM" = 'LI',
      "SICM" = 'LI',
      "RECU" = 'LI',
      'ANAL' = 'LI'
    )
    
    metadat$site[metadat$site %in% c('IICD', 'TOC')] <- 'Cecum'
    metadat$site <- map_site_orgam[metadat$site]
    
    metadat$sampletype[metadat$sampletype == 'U'] = 'F'
  }
  
  # load and filter the data: rel in all samples > 0.1%
  if (FALSE) {
    # filter
    microbe_rel <- sweep(microbe_count, 2, apply(microbe_count, 2, sum), '/')
    isRM <- apply(X = microbe_rel, MARGIN = 1, FUN = function(ar) {
      all(ar < 0.001) # 0.1%
    }); table(isRM) # FALSE  TRUE: 42189 37504
    microbe_rel <- microbe_rel[!isRM, ]
    microbe_count <- microbe_count[rownames(microbe_rel), ]
    
    # rarefy
    if (TRUE) {
      if (!file.exists('Result/rarefy/before_rm_contm/microbe_count_rare.tsv')) {
        write.table(x = cbind(ID = rownames(microbe_count), microbe_count), 
                    file = 'Result/rarefy/before_rm_contm/microbe_count.tsv', quote = FALSE, sep = '\t', 
                    row.names = FALSE)
      }
      
      ## after process
      if (file.exists('Result/rarefy/before_rm_contm/microbe_count_rare.tsv')) {
        tmp <- read.csv(file = 'Result/rarefy/before_rm_contm/microbe_count_rare.tsv', header = TRUE, 
                        skip = 1, sep = '\t', row.names = 1, check.names = FALSE, 
                        strip.white = TRUE, stringsAsFactors = FALSE)
        
        all(colnames(microbe_count) == colnames(tmp))
        all(rownames(tmp) %in% rownames(microbe_count))
        microbe_rare <- tmp
      }
    }
  }
  # save(microbe_count, microbe_rel, microbe_rare, microbe_feature, file = 'Result/00.Saved.input.flt_rarefy.RData')
  load(file = 'Result/00.Saved.input.flt_rarefy.RData')
  
  ## 
  Saved.contamdf = list()
  sapply(list(list(Site = 'LI', negType = c('M')), 
              list(Site = 'Appendix', negType = c('S')),
              list(Site = 'SI', negType = c('M')), # M, S
              list(Site = 'Stomach', negType = c('M', 'S')), # M, S
              list(Site = 'Esophageal', negType = c('M')),
              list(Site = 'Oral cavity', negType = c('S'))
              ), 
         function (x) {
           # Site = 'LI'; negType <- c('M') 
           Site = x$Site; negType = x$negType; print(sprintf('[Info] Process %s ...', Site))
           pheno.used = metadat[metadat$site == 'Negtive control' |
                                  (metadat$site == Site & metadat$sampletype %in% negType), ]
           genus.used = microbe_rare[, rownames(pheno.used)]
           
           #
           phyl <- phyloseq(otu_table(genus.used, taxa_are_rows = TRUE),
                            sample_data(pheno.used)
                            , tax_table(as.matrix(microbe_feature[rownames(genus.used), ]))
                            , read_tree(treefile = 'Result/tree-rooted.nwk')
           )
           
           sample_data(phyl)$is.neg <- sample_data(phyl)$site == "Negtive control"
           
           ## 02, use this if only sample info (neg and real, for prev) is available
           contamdf.prev <- isContaminant(phyl, method="prevalence", neg="is.neg", threshold = 0.1, detailed = TRUE, normalize = TRUE)
           contamdf.prev <- cbind(contamdf.prev, microbe_feature[rownames(contamdf.prev), ])
           table(contamdf.prev$contaminant)
           
           # write.csv(contamdf.prev[order(contamdf.prev$contaminant, decreasing = TRUE), ], 
           #           file = sprintf('Result/contamdf/contamdf.prev-%s.csv', Site))
           
           ### the abundance of these taxa in real samples
           if (TRUE) {
             taxa_conts <- rownames(contamdf.prev)[contamdf.prev$contaminant] 
             genus.rel_contamdf.pos <- microbe_rel[taxa_conts, 
                                                   rownames(pheno.used)[pheno.used$site != 'Negtive control']]
             genus.rel_contamdf.neg <- microbe_rel[taxa_conts, 
                                                   rownames(pheno.used)[pheno.used$site == 'Negtive control']]
             
             genus.rel_contamdf.pos <- asin(sqrt(genus.rel_contamdf.pos) )
             genus.rel_contamdf.neg <- asin(sqrt(genus.rel_contamdf.neg) )
             
             df.pa <- data.frame(mean.pos = rowMeans(genus.rel_contamdf.pos), 
                                 mean.neg = rowMeans(genus.rel_contamdf.neg), 
                                 sd.pos = apply(genus.rel_contamdf.pos, 1, sd),
                                 sd.neg = apply(genus.rel_contamdf.neg, 1, sd),
                                 median.pos = apply(genus.rel_contamdf.pos, 1, median), 
                                 median.neg = apply(genus.rel_contamdf.neg, 1, median),
                                 IQR.pos = apply(genus.rel_contamdf.pos, 1, IQR),
                                 IQR.neg = apply(genus.rel_contamdf.neg, 1, IQR))
             
             # write.csv(df.pa[order(df.pa$mean.pos, decreasing = TRUE), ], 
             #           file = sprintf('Result/contamdf/df.pa-%s.csv', Site))
             
             g <- ggplot(data=df.pa, aes(x=mean.neg, y=mean.pos)) + 
               geom_errorbarh(aes(xmin=mean.neg-sd.neg/sqrt(dim(genus.rel_contamdf.neg)[2]), 
                                  xmax=mean.neg+sd.neg/sqrt(dim(genus.rel_contamdf.neg)[2])), color = 'grey', linetype = 'solid', size = 0.5, height = 0.0001) +
               geom_errorbar(aes(ymin=mean.pos-sd.pos/sqrt(dim(genus.rel_contamdf.pos)[2]), 
                                 ymax=mean.pos+sd.pos/sqrt(dim(genus.rel_contamdf.pos)[2])), color = 'grey', linetype = 'solid', size = 0.5, width = 0.001) +
               geom_hline(yintercept = asin(sqrt(0.001) ), color = '#ffa700', linetype = 'dashed', size = 0.5) +
               geom_vline(xintercept = asin(sqrt(0.001) ), color = '#ffa700', linetype = 'dashed', size = 0.5) +
               geom_point(color = 'black', fill = 'red', alpha = 0.5, shape = 21) + 
               labs(title = sprintf('%s', Site)) + 
               xlab("Normalized Abundance\n(Negative Controls)") + 
               ylab("Normalized Abundance\n(True Samples)")
             g <- g + theme_bw() +
               theme(legend.position = "none", 
                     panel.grid = element_blank(), 
                     axis.text.x = element_text(size = 10,colour = 'black', angle = 0, hjust = 0.5, vjust = 0.5), # element_text(size = 15), 
                     axis.text.y = element_text(size = 10, colour = 'black'), 
                     axis.title = element_text(size = 11, colour = 'black'), 
                     panel.border = element_rect(colour = 'black')
               ) + 
               scale_x_continuous(expand = expand_scale(mult = c(0.05, 0.1))) +
               scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.1)))
             # g
             # ggsave(filename = sprintf('Result/contamdf/Abun_real_negCont.%s.pdf', Site), plot = g, width = 3.53, height = 3.08)
           }
           
           Saved.contamdf[[Site]] <<- df.pa
           
           # -. get the clean data
           print(colSums(microbe_count[, rownames(metadat)[metadat$site == 'Negtive control'] ]))
           microbe_count[rownames(df.pa), rownames(pheno.used)[pheno.used$site != 'Negtive control']] <<- 0
           print(colSums(microbe_count[, rownames(metadat)[metadat$site == 'Negtive control'] ]))
           
           return(NA)
        })
  save(Saved.contamdf, file = 'Result/contamdf/00.Saved.contamdf.RData')
  
  
  # -. get the clean data
  ## - phyla
  if (TRUE) {
    tmp <- microbe_feature[rownames(microbe_count), ]
    tmp$ID <- apply(X = tmp[, c('Kingdom', 'Phylum')], MARGIN = 1, FUN = function(ar) paste(ar, collapse = ';'))
    tmp.p <- unique(as.character(tmp$ID))
    phyla <- do.call('rbind', lapply(X = tmp.p, FUN = function(i) {
      i <- rownames(tmp)[tmp$ID == i]
      if (length(i) > 1) {
        r = colSums(microbe_count[i, ])
      } else {
        r = as.numeric(microbe_count[i, ])
      }
      return(r)
    }))
    rownames(phyla) <- tmp.p
    phyla <- as.data.frame(phyla)
  }
  
  ## - genus
  if (TRUE) {
    tmp <- microbe_feature[rownames(microbe_count), ]
    tmp$ID <- apply(X = tmp[, names(mapping)], MARGIN = 1, FUN = function(ar) paste(ar, collapse = ';'))
    tmp.g <- unique(as.character(tmp$ID))
    genus <- do.call('rbind', lapply(X = tmp.g, FUN = function(i) {
      i <- rownames(tmp)[tmp$ID == i]
      if (length(i) > 1) {
        r = colSums(microbe_count[i, ])
      } else {
        r = as.numeric(microbe_count[i, ])
      }
      return(r)
    }))
    rownames(genus) <- tmp.g
    genus <- as.data.frame(genus)
  }
  
  ## - species
  if (TRUE) {
    tmp <- microbe_feature[rownames(microbe_count), ]
    tmp$ID <- apply(X = tmp[, c(names(mapping), 'Species')], MARGIN = 1, FUN = function(ar) paste(ar, collapse = ';'))
    tmp.s <- unique(as.character(tmp$ID))
    species <- do.call('rbind', lapply(X = tmp.s, FUN = function(i) {
      i <- rownames(tmp)[tmp$ID == i]
      if (length(i) > 1) {
        r = colSums(microbe_count[i, ])
      } else {
        r = as.numeric(microbe_count[i, ])
      }
      return(r)
    }))
    rownames(species) <- tmp.s
    species <- as.data.frame(species)
  }
  
  save(df_tax, df_tax_blast, df_tax_wSpecies, 
       microbe_count, microbe_feature, phyla, genus, species, 
       file = 'Result/00.Data(asv, p, g, s), after rm contm.RData')
}


####### Start here ######## <-----------------------------------------
# - contamination
load(file = 'Result/contamdf/00.Saved.contamdf.RData')


# - load the data
# df_tax, df_tax_blast, df_tax_wSpecies, 
# microbe_count(raw), microbe_feature, metadat
load(file = 'Result/00.Saved.input.RData'); tmp = microbe_count

metadat$Antibiotics[metadat$subject == 'S006'] <- 'Y'
metadat$site[metadat$site %in% c('IICD', 'TOC')] <- 'Cecum'
metadat$sampletype[metadat$sampletype == 'U'] = 'F'

# df_tax, df_tax_blast, df_tax_wSpecies, 
# microbe_count(after_rm_contm), microbe_feature, phyla, genus, species
load(file = 'Result/00.Data(asv, p, g, s), after rm contm.RData')
if (FALSE) {
  write.csv(x = data.frame(smaple = colnames(tmp), ori = colSums(tmp), rm.contm = colSums(microbe_count)), 
            file = 'Result/00.Data, depth.csv', row.names = FALSE)
}
if (TRUE) {
  if (!file.exists('Result/rarefy/after_rm_contm/microbe_count_rare.tsv')) {
    write.table(x = cbind(ID = rownames(microbe_count), microbe_count), 
                file = 'Result/rarefy/after_rm_contm/microbe_count.tsv', quote = FALSE, sep = '\t', 
                row.names = FALSE)
  }
  
  ## after process
  if (file.exists('Result/rarefy/after_rm_contm/microbe_count_rare.tsv')) {
    tmp <- read.csv(file = 'Result/rarefy/after_rm_contm/microbe_count_rare.tsv', header = TRUE, 
                    skip = 1, sep = '\t', row.names = 1, check.names = FALSE, 
                    strip.white = TRUE, stringsAsFactors = FALSE)
    
    all(colnames(microbe_count) == colnames(tmp))
    all(rownames(tmp) %in% rownames(microbe_count))
    microbe_rare <- tmp
  }
}
if (TRUE) {
  if (!file.exists('Result/rarefy/after_rm_contm/species_rare.tsv')) {
    write.table(x = cbind(ID = rownames(species), species), 
                file = 'Result/rarefy/after_rm_contm/species.tsv', quote = FALSE, sep = '\t', 
                row.names = FALSE)
  }
  
  ## after process
  if (file.exists('Result/rarefy/after_rm_contm/species_rare.tsv')) {
    tmp <- read.csv(file = 'Result/rarefy/after_rm_contm/species_rare.tsv', header = TRUE, 
                    skip = 1, sep = '\t', row.names = 1, check.names = FALSE, 
                    strip.white = TRUE, stringsAsFactors = FALSE)
    
    all(colnames(species) == colnames(tmp))
    all(rownames(tmp) %in% rownames(species))
    species_rare <- tmp
  }
}


# - seq id and seq
if (FALSE) {
  query_seq = seqinr::read.fasta(file = 'Result/rep_seq.fa')
  tmp = sapply(query_seq, FUN = function(x) toupper(seqinr::getSequence(x, as.string = TRUE)[[1]]))
  query_seq = tmp; names(query_seq) = names(tmp)
  
  save(query_seq, file = 'Result/rep_seq.fa.RData')
}
load(file = 'Result/rep_seq.fa.RData') # query_seq
#######    END    ######## <-----------------------------------------


if (TRUE) {
  ### Lumen-Mucosa, Only SI
  if (FALSE) {
    lv <- c(
      'Peri-duodenum' = 'DUPA',
      'Jejunum 1m' = "JEJ100",
      'Ileum 3m' = "IIE300",
      'Ileocecal' = "IIC")
    
    metadat_used <- metadat[metadat$MajorBS == 'SmallIntestine', ]
    metadat_used <- metadat_used[metadat_used$site %in% lv, ]
    
    metadat_used$sampletype[metadat_used$sampletype %in% c('F', 'S')] <- 'L'
    metadat_used$sampletype <- factor(metadat_used$sampletype, c('L', 'M'))
    
    metadat_used$Site_belonging <- metadat_used$site
    metadat_used$site <- factor(metadat_used$site, lv)
  }
  
  ### Lumen-Mucosa, Only LI
  if (FALSE) {
    lv <- c(#
      'Cecum' = 'Cecum',
      'Ascending colon' = "ASCM",
      'Transverse colon' = "TRC",
      'Descending colon' = "DECM",
      'Signoid colon' = "SICM",
      'Rectum' = "RECU"
    )
    
    metadat_used <- metadat[metadat$MajorBS == 'LargeIntestine', ]
    metadat_used <- metadat_used[metadat_used$site %in% lv, ]
    
    metadat_used$sampletype[metadat_used$sampletype %in% c('F', 'S')] <- 'L'
    metadat_used$sampletype <- factor(metadat_used$sampletype, c('L', 'M'))
    
    metadat_used$Site_belonging <- metadat_used$site
    metadat_used$site <- factor(metadat_used$site, lv)
  }
  
  
  ### -- logic & beta-bionominal regression
  library(lme4)
  library(glmmML)
  library(aod)
  
  used <- lv
  pheno.used <- metadat_used
  pheno.used$group <- factor(pheno.used$Site_belonging, used)
  
  ##### - Only SI : detection = 1/1000, prevalence = 10/100 and fdr < 0.05 for L-M
  if (FALSE) {
    ## get rarefied data
    genus_rarefy <- species_rare
    table(colSums(genus_rarefy))
    
    ##
    L_M.SI <- list()
    for (grp in used) {
      # grp <- used[1]
      pheno.used_LM <- pheno.used[pheno.used$group == grp, ]
      pheno.used_LM.patient <- unique(pheno.used_LM[, c("subject", "SEX", "AGE", "Antibiotics", "BMI")])
      rownames(pheno.used_LM.patient) <- pheno.used_LM.patient$subject
      
      genus.used <- genus_rarefy[, rownames(pheno.used_LM)]
      genus.used <- genus.used[rowSums(genus.used) > 0, ]
      
      ## find core microbes
      pseq <- phyloseq(otu_table(genus.used, taxa_are_rows = TRUE))
      pseq.comp <- microbiome::transform(pseq, "compositional")
      core_taxa <- core_members(pseq.comp, detection = 1/1000, prevalence = 10/100)
      
      genus.used <- genus.used[core_taxa, ]
      Stats <- do.call('rbind', lapply(X = core_taxa, FUN = function(otu) {
        # otu <- core_taxa[1]
        print(otu)
        tab <- data.frame(otu = as.integer(genus.used[otu, ]), 
                          others = 6000 - as.integer(genus.used[otu, ]),
                          abun = as.integer(genus.used[otu, ]) / 6000,
                          Type = pheno.used_LM[colnames(genus.used), 'sampletype'],
                          sex = pheno.used_LM[colnames(genus.used), 'SEX'],
                          age = as.numeric(pheno.used_LM[colnames(genus.used), 'AGE']),
                          antibiotic = pheno.used_LM[colnames(genus.used), 'Antibiotics'],
                          patient = pheno.used_LM[colnames(genus.used), 'subject'])
      
        ##
        corr <- as.data.frame(matrix(data = NA, nrow = length(unique(tab$patient)), ncol = 2, 
                                     dimnames = list(unique(tab$patient), c('M', 'L'))))
        tmp <- apply(tab[, c('patient', 'Type', 'abun')], 1, function(ar) {
          corr[ar[1], ar[2]] <<- as.numeric(ar[3])
          return(NA)
        })
        corr <- corr[!(is.na(corr$M) | is.na(corr$L)), ]
        if (TRUE) {
          tmp <- cor.test(corr$M, corr$L, method = 'spearman')
          rd <- c(sp.pval = tmp$p.value, sp = tmp$estimate)
          
          tmp <- cor.test(corr$M, corr$L, method = 'pearson')
          rd <- c(rd, pr.pval = tmp$p.value, pr = tmp$estimate)
          
          tmp = cbind(corr, pheno.used_LM.patient[rownames(corr), ])
          tmp$r1 = tmp$M; tmp$r2 = tmp$L
          tmp$AGE = as.numeric(tmp$AGE); tmp$BMI = as.numeric(tmp$BMI)
          tag <- tryCatch(
            # Specifying expression
            expr = {                     
              fit <- PResiduals::partial_Spearman(formula = r1|r2~AGE+BMI, 
                                                  data = tmp, 
                                                  fit.x = 'orm', fit.y = 'orm')
              c(psp.pval = fit$TS$TB$pval, psp.rho = fit$TS$TB$ts)
            },
            # Specifying error message
            error = function(e){         
              c(psp.pval = NA, psp.rho = NA)
            },
            warning = function(w){      
              c(psp.pval = NA, psp.rho = NA)
            }
          )
          rd <- c(rd, tag)
        }
        
        
        ## 
        tab$Type <- as.factor(tab$Type)
        fit <- betabin(cbind(otu, others) ~ Type + age, random = ~ Type, data = tab)
        fit_coef <- summary(fit)@Coef
        
        c(beta = fit_coef['TypeM', 'Estimate'], pval = fit_coef['TypeM', 'Pr(> |z|)'], bonf = NA, fdr = NA, 
          rd)
      } ) )
      Stats <- as.data.frame(Stats); rownames(Stats) <- core_taxa
      Stats$bonf <- p.adjust(Stats$pval, 'bonferroni')
      Stats$fdr <- p.adjust(Stats$pval, 'fdr')
      
      taxa = as.data.frame(do.call('rbind', lapply(core_taxa, function(x) {
        x = stringi::stri_split(x, regex = ';')[[1]]
        paste(c(mapping,"s"), x, sep = '__')
      })), stringsAsFactors = FALSE)
      colnames(taxa) = c(mapping, 's')
      rownames(taxa) = core_taxa
      
      Stats <- cbind(Stats, taxa)
      L_M.SI[[grp]] <- Stats
    }
    save(L_M.SI, file = 'Result/Beta-bin/beta-bin.SI.Rds')
    
    
    ## filtering
    # load(file = 'Result/Beta-bin/beta-bin.SI.Rds')
    core_taxa <- lapply(X = used, FUN = function(grp) {
      Stats <- L_M.SI[[grp]]
      Stats <- Stats[!is.na(Stats$pval), ]
      rownames(Stats)[Stats$fdr < 0.05]
    })
    core_taxa <- Reduce(union, core_taxa)
    
    taxa = as.data.frame(do.call('rbind', lapply(core_taxa, function(x) {
      x = stringi::stri_split(x, regex = ';')[[1]]
      paste(c(mapping,"s"), x, sep = '__')
    })), stringsAsFactors = FALSE)
    colnames(taxa) = c(mapping, 's')
    rownames(taxa) = core_taxa
    taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
      f <- get_lv(x = x, usedLv = 'f')
      g <- get_lv(x = x, usedLv = 'g')
      s <- get_lv(x = x, usedLv = 's')
      if (startsWith(s, 's__')) {
        s <- paste(g, s, sep = ';')
      }
      if (s == 'g__NA;s__NA') return(paste(f, s, sep = ';')) else return(s)
    })
    
    tab <- do.call('rbind', lapply(X = core_taxa, function(otu) {
      sapply(X = used, FUN = function(grp) {
        Stats <- L_M.SI[[grp]]
        Stats <- Stats[!is.na(Stats$pval), ]
        if (otu %in% rownames(Stats) && Stats[otu, "fdr"] < 0.05) {
          return(Stats[otu, "beta"])
        } else {
          return(NA)
        }
      })
    }))
    tab <- as.data.frame(tab)
    rownames(tab) <- core_taxa

    tab$stats <- apply(X = tab, MARGIN = 1, FUN = function(r) {
      r <- r[!is.na(r)]
      if (all(sign(r) == 1)) {
        return('pos')
      } else if (all(sign(r) == -1)) {
        return('neg')
      } else {
        return('conflict')
      }
    })
    # write.csv(cbind(tab, taxa[rownames(tab), ]), file = 'Result/Beta-bin/beta-bin.SI.csv')
    
    
    ## draw
    col_phyla <- c("p__Proteobacteria" = '#1b9e77',
                   "p__Actinobacteriota" = '#d95f02',
                   "p__Firmicutes" = '#7570b3',       # 1
                   "p__Bacteroidota" = '#386cb0',
                   "p__Fusobacteria" = '#66a61e',
                   "p__Verrucomicrobia" = '#e6ab02',
                   'p__Deinococcota' = 'black')
    lab <- taxa$Name; names(lab) <- rownames(taxa)
    
    tab.flt <- tab[tab$stats != 'conflict', colnames(tab) != 'stats']
    tab.flt_stat <- data.frame(Var1 = rownames(tab.flt), 
                               mean = apply(X = tab.flt, MARGIN = 1, FUN = function(r) mean(r[!is.na(r)])), 
                               median = apply(X = tab.flt, MARGIN = 1, FUN = function(r) median(r[!is.na(r)])), 
                               sd = apply(X = tab.flt, MARGIN = 1, FUN = function(r) {
                                 val <- sd(r[!is.na(r)])
                                 if (!is.na(val)) return(val) else return(0)
                               }))
    
    tab.flt <- reshape2::melt(as.matrix(tab.flt))
    tab.flt <- tab.flt[!is.na(tab.flt$value), ]
    tab.flt$NAME <- taxa[as.character(tab.flt$Var1), 'Name']
    # tab.flt <- tab.flt[order(tab.flt_stat[as.character(tab.flt$Var1), 'mean'], decreasing = FALSE), ]
    tab.flt$Var1 <- factor(tab.flt$Var1, rownames(tab.flt_stat)[order(tab.flt_stat$mean, decreasing = FALSE)] )
    
    p <- ggplot(tab.flt, aes(x = Var1, y = value, col = Var2)) +
      geom_hline(yintercept = 0, colour = '#99aab5', lwd = 0.5, lty = 'dashed') +
      geom_point(aes(fill = Var2), shape = 21, col = 'black', size = 1.5) +
      geom_point(aes(y = mean), data = tab.flt_stat, shape = '|', colour = 'black') +
      labs(x = '', y = "Beta", title = '')
    p <- p + theme_bw() +
      theme(legend.position = "none", 
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 7, color = 'black'), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid = element_blank(), 
            axis.title = element_text(size = 12)
      ) + 
      scale_x_discrete(labels = lab) +
      scale_fill_manual(values = c('Peri-duodenum' = '#377eb8', 
                                   'Jejunum 1m' = '#984ea3',
                                   'Ileum 3m' = '#ffff33',
                                   'Ileocecal 1cm' = '#f781bf',
                                   'Ileocecal' = '#999999')) +
      scale_shape_manual(values = c('Peri-duodenum' = 19, 
                                    'Jejunum 1m' = 15,
                                    'Ileum 3m' = 23,
                                    'Ileocecal 1cm' = 24,
                                    'Ileocecal' = 25)) +
      coord_flip()
    p
    ggsave(filename = 'Result/Beta-bin/beta-bin.SI.pdf', 
           plot = p, width = 5.5, height = 5.5)
    stop('')
  }
  ##### - Only LI : detection = 1/1000, prevalence = 10/100 and fdr < 0.05 for L-M
  if (FALSE) {
    ## get rarefied data
    genus_rarefy <- species_rare
    table(colSums(genus_rarefy))
    
    ##
    L_M.LI <- list()
    for (grp in used) {
      # grp <- used[1]
      pheno.used_LM <- pheno.used[pheno.used$group == grp, ]
      pheno.used_LM.patient <- unique(pheno.used_LM[, c("subject", "SEX", "AGE", "Antibiotics", "BMI")])
      rownames(pheno.used_LM.patient) <- pheno.used_LM.patient$subject
      
      genus.used <- genus_rarefy[, rownames(pheno.used_LM)]
      genus.used <- genus.used[rowSums(genus.used) > 0, ]
      
      ## find core microbes
      pseq <- phyloseq(otu_table(genus.used, taxa_are_rows = TRUE))
      pseq.comp <- microbiome::transform(pseq, "compositional")
      core_taxa <- core_members(pseq.comp, detection = 1/1000, prevalence = 10/100)
      
      genus.used <- genus.used[core_taxa, ]
      Stats <- do.call('rbind', lapply(X = core_taxa, FUN = function(otu) {
        # otu <- core_taxa[1]
        print(otu)
        tab <- data.frame(otu = as.integer(genus.used[otu, ]), 
                          others = 6000 - as.integer(genus.used[otu, ]),
                          abun = as.integer(genus.used[otu, ]) / 6000,
                          Type = pheno.used_LM[colnames(genus.used), 'sampletype'],
                          sex = pheno.used_LM[colnames(genus.used), 'SEX'],
                          age = as.numeric(pheno.used_LM[colnames(genus.used), 'AGE']),
                          antibiotic = pheno.used_LM[colnames(genus.used), 'Antibiotics'],
                          patient = pheno.used_LM[colnames(genus.used), 'subject'])
        
        ##
        corr <- as.data.frame(matrix(data = NA, nrow = length(unique(tab$patient)), ncol = 2, 
                                     dimnames = list(unique(tab$patient), c('M', 'L'))))
        tmp <- apply(tab[, c('patient', 'Type', 'abun')], 1, function(ar) {
          corr[ar[1], ar[2]] <<- as.numeric(ar[3])
          return(NA)
        })
        corr <- corr[!(is.na(corr$M) | is.na(corr$L)), ]
        if (TRUE) {
          tmp <- cor.test(corr$M, corr$L, method = 'spearman')
          rd <- c(sp.pval = tmp$p.value, sp = tmp$estimate)
          
          tmp <- cor.test(corr$M, corr$L, method = 'pearson')
          rd <- c(rd, pr.pval = tmp$p.value, pr = tmp$estimate)
          
          tmp = cbind(corr, pheno.used_LM.patient[rownames(corr), ])
          tmp$r1 = tmp$M; tmp$r2 = tmp$L
          tmp$AGE = as.numeric(tmp$AGE); tmp$BMI = as.numeric(tmp$BMI)
          tag <- tryCatch(
            # Specifying expression
            expr = {                     
              fit <- PResiduals::partial_Spearman(formula = r1|r2~AGE+BMI, 
                                                  data = tmp, 
                                                  fit.x = 'orm', fit.y = 'orm')
              c(psp.pval = fit$TS$TB$pval, psp.rho = fit$TS$TB$ts)
            },
            # Specifying error message
            error = function(e){         
              c(psp.pval = NA, psp.rho = NA)
            },
            warning = function(w){      
              c(psp.pval = NA, psp.rho = NA)
            }
          )
          rd <- c(rd, tag)
        }
        
        
        ## 
        tab$Type <- as.factor(tab$Type)
        fit <- betabin(cbind(otu, others) ~ Type + age, random = ~ Type, data = tab)
        fit_coef <- summary(fit)@Coef
        
        c(beta = fit_coef['TypeM', 'Estimate'], pval = fit_coef['TypeM', 'Pr(> |z|)'], bonf = NA, fdr = NA, 
          rd)
      } ) )
      Stats <- as.data.frame(Stats); rownames(Stats) <- core_taxa
      Stats$bonf <- p.adjust(Stats$pval, 'bonferroni')
      Stats$fdr <- p.adjust(Stats$pval, 'fdr')
      
      taxa = as.data.frame(do.call('rbind', lapply(core_taxa, function(x) {
        x = stringi::stri_split(x, regex = ';')[[1]]
        paste(c(mapping,"s"), x, sep = '__')
      })), stringsAsFactors = FALSE)
      colnames(taxa) = c(mapping, 's')
      rownames(taxa) = core_taxa
      
      Stats <- cbind(Stats, taxa)
      L_M.LI[[grp]] <- Stats
    }
    save(L_M.LI, file = 'Result/Beta-bin/beta-bin.LI.Rds')
    
    
    ## filtering
    # load(file = 'Result/Beta-bin/beta-bin.LI.Rds')
    core_taxa <- lapply(X = used, FUN = function(grp) {
      Stats <- L_M.LI[[grp]]
      Stats <- Stats[!is.na(Stats$pval), ]
      rownames(Stats)[Stats$fdr < 0.05]
    })
    core_taxa <- Reduce(union, core_taxa)
    
    taxa = as.data.frame(do.call('rbind', lapply(core_taxa, function(x) {
      x = stringi::stri_split(x, regex = ';')[[1]]
      paste(c(mapping,"s"), x, sep = '__')
    })), stringsAsFactors = FALSE)
    colnames(taxa) = c(mapping, 's')
    rownames(taxa) = core_taxa
    taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
      f <- get_lv(x = x, usedLv = 'f')
      g <- get_lv(x = x, usedLv = 'g')
      s <- get_lv(x = x, usedLv = 's')
      if (startsWith(s, 's__')) {
        s <- paste(g, s, sep = ';')
      }
      if (s == 'g__NA;s__NA') return(paste(f, s, sep = ';')) else return(s)
    })
    
    tab <- do.call('rbind', lapply(X = core_taxa, function(otu) {
      sapply(X = used, FUN = function(grp) {
        print(grp)
        Stats <- L_M.LI[[grp]]
        Stats <- Stats[!is.na(Stats$pval), ]
        if (otu %in% rownames(Stats) && Stats[otu, "fdr"] < 0.05) {
          return(Stats[otu, "beta"])
        } else {
          return(NA)
        }
      })
    }))
    tab <- as.data.frame(tab)
    rownames(tab) <- core_taxa
    
    tab$stats <- apply(X = tab, MARGIN = 1, FUN = function(r) {
      r <- r[!is.na(r)]
      if (all(sign(r) == 1)) {
        return('pos')
      } else if (all(sign(r) == -1)) {
        return('neg')
      } else {
        return('conflict')
      }
    })
    write.csv(cbind(tab, taxa[rownames(tab), ]), file = 'Result/Beta-bin/beta-bin.LI.csv')
    
    
    ## draw
    col_phyla <- c("p__Proteobacteria" = '#1b9e77',
                   "p__Actinobacteriota" = '#d95f02',
                   "p__Firmicutes" = '#7570b3',       # 1
                   "p__Bacteroidota" = '#386cb0',
                   "p__Fusobacteria" = '#66a61e',
                   "p__Verrucomicrobia" = '#e6ab02',
                   'p__Deinococcota' = 'black')
    lab <- taxa$Name; names(lab) <- rownames(taxa)
    
    tab.flt <- tab[tab$stats != 'conflict', colnames(tab) != 'stats']
    tab.flt_stat <- data.frame(Var1 = rownames(tab.flt), 
                               mean = apply(X = tab.flt, MARGIN = 1, FUN = function(r) mean(r[!is.na(r)])), 
                               median = apply(X = tab.flt, MARGIN = 1, FUN = function(r) median(r[!is.na(r)])), 
                               sd = apply(X = tab.flt, MARGIN = 1, FUN = function(r) {
                                 val <- sd(r[!is.na(r)])
                                 if (!is.na(val)) return(val) else return(0)
                               }))
    
    tab.flt <- reshape2::melt(as.matrix(tab.flt))
    tab.flt <- tab.flt[!is.na(tab.flt$value), ]
    tab.flt$NAME <- taxa[as.character(tab.flt$Var1), 'Name']
    # tab.flt <- tab.flt[order(tab.flt_stat[as.character(tab.flt$Var1), 'mean'], decreasing = FALSE), ]
    tab.flt$Var1 <- factor(tab.flt$Var1, rownames(tab.flt_stat)[order(tab.flt_stat$mean, decreasing = FALSE)] )
    
    p <- ggplot(tab.flt, aes(x = Var1, y = value, col = Var2)) +
      geom_hline(yintercept = 0, colour = '#99aab5', lwd = 0.5, lty = 'dashed') +
      geom_point(aes(fill = Var2), shape = 21, col = 'black', size = 1.5) +
      geom_point(aes(y = mean), data = tab.flt_stat, shape = '|', colour = 'black') +
      labs(x = '', y = "Beta", title = '')
    p <- p + theme_bw() +
      theme(legend.position = "none", 
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 7, color = 'black'), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid = element_blank(), 
            axis.title = element_text(size = 12)
      ) + 
      scale_x_discrete(labels = lab) +
      scale_fill_manual(values = c('Cecum' = '#377eb8', 
                                   'Ascending colon' = '#4daf4a',
                                   'Transverse colon' = '#984ea3',
                                   'Descending colon' = '#ff7f00',
                                   'Signoid colon' = '#ffff33',
                                   'Rectum' = '#a65628')) +
      scale_shape_manual(values = c('Cecum' = 19, 
                                    'Ascending colon' = 15,
                                    'Transverse colon' = 23,
                                    'Descending colon' = 24,
                                    'Signoid colon' = 25,
                                    'Rectum' = 10)) +
      coord_flip()
    p
    ggsave(filename = 'Result/Beta-bin/beta-bin.LI.pdf', 
           plot = p, width = 5.5, height = 5.5)
    stop('')
  }
  
  #####-
  L_M.SI_stat <- read.csv(file = 'Result/Beta-bin/beta-bin.SI.csv', header = TRUE,
                          row.names = 1, check.names = FALSE, stringsAsFactors = FALSE,
                          strip.white = TRUE)
  L_M.LI_stat <- read.csv(file = 'Result/Beta-bin/beta-bin.LI.csv', header = TRUE,
                          row.names = 1, check.names = FALSE, stringsAsFactors = FALSE,
                          strip.white = TRUE)
  
  ID = intersect(rownames(L_M.SI_stat)[L_M.SI_stat$stats == 'pos'], rownames(L_M.LI_stat)[L_M.LI_stat$stats == 'pos'])
  cat(L_M.SI_stat[ID, 'Name'], sep = '\n')
  
  ID = intersect(rownames(L_M.SI_stat)[L_M.SI_stat$stats == 'neg'], rownames(L_M.LI_stat)[L_M.LI_stat$stats == 'neg'])
  L_M.SI_stat[ID, 'Name']
}


if (TRUE) {
  ## obtain the H. pylriy lv
  if (FALSE) {
    map_stoamch = c('SF' = 'Stomach',
                    'SB' = 'Stomach',
                    "GJ" = 'Stomach',
                    "SA" = 'Stomach',
                    "PY" = 'Stomach')
    
    metadat_stomach = metadat[metadat$site %in% names(map_stoamch), ]
    
    tmp = rownames(microbe_feature)[!is.na(microbe_feature$Species) & 
                                      !is.na(microbe_feature$Genus) & 
                                      microbe_feature$Genus == 'Helicobacter' & 
                                      microbe_feature$Species == 'pylori']
    Hp = microbe_rare[tmp[tmp %in% rownames(microbe_rare)], rownames(metadat_stomach)]
    Hp = colSums(Hp)
    
    metadat_stomach$Hp = Hp[rownames(metadat_stomach)] / 6000
    
    Hp = sapply(unique(metadat_stomach$subject), function(x) {
      mean(metadat_stomach$Hp[metadat_stomach$subject == x])
    })
    summary(Hp)
    
    save(Hp, file = 'Result/otu_cores/H_pylori.RData')
  }
  load(file = 'Result/otu_cores/H_pylori.RData') # Hp
  
  taxonomy = microbe_feature
  genus_rarefy <- microbe_rare
  #all(rownames(genus_rarefy) %in% rownames(taxonomy))
  
  
  ## - get oral gut axis at ind level
  if (FALSE) {
    map_site_orgam <- c(
      'ORALmix' = 'Oral cavity'
      ,
      'AP' = 'Appendix'
      ,
      'Cecum' = 'Cecum',
      "ASCM" = 'ASCM',
      "TRC" = 'TRC',
      "DECM" = 'DECM',
      "SICM" = 'SICM',
      "RECU" = 'RECU',
      'ANAL' = 'ANAL'
    )
    
    
    metadat_core <- metadat
    metadat_core <- metadat_core[(metadat_core$site %in% names(map_site_orgam)), ]
    metadat_core$site <- map_site_orgam[metadat_core$site]
    table(metadat_core$site)
    
    ## - 
    genus_rarefy <- genus_rarefy[, rownames(metadat_core)]
    
    otu_cores_ind <- list()
    tmp <- lapply(X = unique(metadat_core$subject), FUN = function(ind) {
      print(ind)
      metadat_core.ind <- metadat_core[metadat_core$subject == ind, ]
      otu_core <- genus_rarefy[, rownames(metadat_core.ind)]
      otu_core <- otu_core[rowSums(otu_core) > 0, ]
      
      NUM_site = dim(metadat_core.ind)[1]
      cutoff = 1/1000 # abun > 0.1%, prev > 0
      depth = 6000
      Ind_site <- lapply(X = unique(metadat_core.ind$site), FUN = function(site) {
        idx <- which(metadat_core.ind$site == site)
        if (length(idx) >= 1) {
          if (length(idx) == 1) {
            otu_core.site <- as.numeric(otu_core[, rownames(metadat_core.ind)[idx] ]); names(otu_core.site) <- rownames(otu_core)
            otu_core.site = otu_core.site / depth
            otu_core.site <- otu_core.site[otu_core.site > cutoff]; 
            
            otu_core_summ <- data.frame(Mean = otu_core.site, 
                                        Prevalence = as.integer(otu_core.site > cutoff), 
                                        row.names = names(otu_core.site), 
                                        stringsAsFactors = FALSE)
            otu_core_summ$Index = (1 / NUM_site) * otu_core_summ$Prevalence
            otu_core_summ <- otu_core_summ[order(otu_core_summ$Prevalence, decreasing = TRUE), ]
          } else {
            otu_core.site <- otu_core[, rownames(metadat_core.ind)[idx] ]
            otu_core.site <- otu_core.site[rowSums(otu_core.site) > 0, ]
            
            pseq <- phyloseq(otu_table(otu_core.site, taxa_are_rows = TRUE))
            pseq.comp <- microbiome::transform(pseq, "compositional")
            taxa <- core_members(pseq.comp, detection = cutoff, prevalence = 0); 
            otu_core.site <- prune_taxa(taxa, pseq)@.Data
            
            otu_core_summ <- apply(X = otu_core.site, MARGIN = 1, FUN = function(otu) {
              otu <- otu / depth
              c('Mean' = mean(otu), 'Prevalence' = sum(otu > cutoff) / length(otu))
            })
            otu_core_summ <- as.data.frame(t(otu_core_summ))
            otu_core_summ$Mean <- as.numeric(as.character(otu_core_summ$Mean))
            otu_core_summ$Prevalence <- as.numeric(as.character(otu_core_summ$Prevalence))
            otu_core_summ$Index = (length(idx) / NUM_site) * otu_core_summ$Prevalence
            otu_core_summ <- otu_core_summ[order(otu_core_summ$Prevalence, decreasing = TRUE), ]
          }
        } else {
          return(NA)
        }
      })
      names(Ind_site) <- unique(metadat_core.ind$site)
      
      taxa <- Reduce(union, lapply(Ind_site, function(site) {
        rownames(site)
      }))
      otu_core_summ <- matrix(data = NA, nrow = length(taxa), ncol = length(Ind_site), 
                              dimnames = list(taxa, names(Ind_site)))
      tmp <- sapply(names(Ind_site), function(n) {
        # print(n)
        site <- Ind_site[[n]]
        otu_core_summ[rownames(site), n] <<- as.numeric(site$Mean) # Index, Prevalence, Mean
        return(NA)
      })
      
      taxa <- taxonomy[rownames(otu_core_summ), ]
      otu_core_summ <- cbind(taxa, otu_core_summ)
      write.csv(otu_core_summ, file = sprintf('Result/otu_cores/oral_gut/Ind-%s.csv', ind))
      
      otu_cores_ind[[ind]] <<- otu_core_summ
    })
    # save(otu_cores_ind, file = 'Result/otu_cores/oral_gut/otu_cores_ind.Rds')
    
    
    ##-
    load(file = 'Result/otu_cores/oral_gut/otu_cores_ind.Rds') # otu_cores_ind
    
    metadat_core.ind <- read.csv(file = 'metadata/patient_info-new.tsv', header = TRUE, sep = '\t', row.names = 1, 
                                 check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
    metadat_core.ind$Hospitalization[grepl('^>', metadat_core.ind$Hospitalization)] <- 10
    
    
    ## -
    if (TRUE) {
      cutoff = 0.001 # prev: 0; abun: 0.001
      Ind_summ.All <- list()
      Ind_summ.Oral_gut <- list()
      ratio.oral_gut <- list()
      tmp <- lapply(X = names(otu_cores_ind), FUN = function(n) {
        ind = otu_cores_ind[[n]]
        taxa <- ind[, colnames(taxonomy)]
        taxa <- taxa[taxa$Kingdom != 'k__Archaea', ]
        ind <- ind[rownames(taxa), ]
        
        ID_sites <- list(oral = rownames(ind)[!is.na(ind$`Oral cavity`) & ind$`Oral cavity` >= cutoff],
                         Appendix = rownames(ind)[!is.na(ind$Appendix) & ind$Appendix >= cutoff],
                         Cecum = rownames(ind)[!is.na(ind$Cecum) & ind$Cecum >= cutoff],
                         ASCM = rownames(ind)[!is.na(ind$ASCM) & ind$ASCM >= cutoff],
                         TRC = rownames(ind)[!is.na(ind$TRC) & ind$TRC >= cutoff],
                         DECM = rownames(ind)[!is.na(ind$DECM) & ind$DECM >= cutoff],
                         SICM = rownames(ind)[!is.na(ind$SICM) & ind$SICM >= cutoff],
                         RECU = rownames(ind)[!is.na(ind$RECU) & ind$RECU >= cutoff],
                         ANAL = rownames(ind)[!is.na(ind$ANAL) & ind$ANAL >= cutoff] )
        
        ID = Reduce(union, ID_sites)
        
        summ <- matrix(data = 0, nrow = length(ID), ncol = length(ID_sites), 
                       dimnames = list(ID, names(ID_sites)))
        tmp <- sapply(names(ID_sites), function(n) {
          # print(n)
          summ[ID_sites[[n]], n] <<- 1
          return(NA)
        })
        summ <- as.data.frame(summ)
        
        Ind_summ.All[[n]] <<- summ
        
        ##
        stat <- rowSums(summ)
        Oral_gut <- ifelse(summ$oral == 1, rowSums(summ), 0)
        
        ratio.oral_gut[[n]] <<- sum(Oral_gut >= 2) / sum(Oral_gut >= 1)
        
        summ$stat <- stat
        summ$Oral_gut <- Oral_gut
        
        Ind_summ.Oral_gut[[n]] <<- summ[summ$Oral_gut >= 2, ]
        
        return(NA)
      })
      # save(Ind_summ.All, Ind_summ.Oral_gut,
      #      file = 'Result/rare/core_bugs/Oral_gut/Ind_summ.RData')
      
      ## -
      summary(unlist(ratio.oral_gut)*100)
      sd(unlist(ratio.oral_gut)*100, na.rm = TRUE)
      #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's      sd
      # 0.00000 0.02241 0.05572 0.05558 0.07740 0.16000       1  0.0395
      
      
      ## -
      if (TRUE) {
        Ind_summ = Ind_summ.Oral_gut
        type = 'Oral_gut'
        ##
        ID <- Reduce(union, lapply(Ind_summ, function(x) rownames(x)))
        summ <- matrix(data = 0, nrow = length(ID), ncol = length(Ind_summ), 
                       dimnames = list(ID, names(Ind_summ)))
        tmp <- sapply(names(Ind_summ), function(n) {
          # print(n)
          ind = Ind_summ[[n]]
          summ[rownames(ind), n] <<- 1
          return(NA)
        })
        summ <- as.data.frame(summ)
        
        ##
        taxa = taxonomy[rownames(summ), c(names(mapping), 'Species')]
        taxa[is.na(taxa)] <- ''
        taxa$Name = apply(taxa, 1, paste, collapse = ';')
        tmp.s = unique(taxa$Name)
        
        sp.summ <- matrix(data = 0, nrow = length(tmp.s), ncol = dim(summ)[2], 
                          dimnames = list(tmp.s, colnames(summ)))
        sp.summ_val <- matrix(data = '', nrow = length(tmp.s), ncol = dim(summ)[2], 
                              dimnames = list(tmp.s, colnames(summ)))
        tmp <- sapply(tmp.s, function(s) {
          idx = rownames(taxa)[taxa$Name == s]
          if (length(idx) == 1) {
            sp.summ[s, ] <<- as.numeric(summ[idx, ])
            sp.summ_val[s, ] <<- sapply(as.numeric(summ[idx, ]), function(x) ifelse(x == 1, idx, ''))
          } else if (length(idx) > 1) {
            sp.summ[s, ] <<- colSums(summ[idx, ])
            sp.summ_val[s, ] <<- apply(summ[idx, ], 2, function(x) paste(idx[x == 1], collapse = '//'))
          }
          return(NA)
        })
        sp.summ <- as.data.frame(sp.summ)
        sp.summ_val <- as.data.frame(sp.summ_val)
        sp.summ = sp.summ[rownames(sp.summ) != "Bacteria;;;;;;", ]
        sp.summ_val = sp.summ_val[rownames(sp.summ), ]
        
        # save(summ, sp.summ, sp.summ_val, file = sprintf('Result/rare/core_bugs/Oral_gut/sp_summ_%s.RData', type))
        
        
        ## -
        sp.summ = sp.summ[, sort(colnames(sp.summ))]
        
        taxa = as.data.frame(do.call('rbind', lapply(rownames(sp.summ), function(x) {
          x = stringi::stri_split(x, regex = ';')[[1]]
          paste(c(mapping,"s"), x, sep = '__')
        })), stringsAsFactors = FALSE)
        colnames(taxa) = c(mapping, 's')
        rownames(taxa) = rownames(sp.summ)
        taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
          # f <- get_lv(x = x, usedLv = 'f')
          g <- get_lv(x = x, usedLv = 'g')
          s <- get_lv(x = x, usedLv = 's')
          if (startsWith(s, 's__')) {
            s <- paste(g, s, sep = ';')
          }
          # return(paste(f, s, sep = ';'))
          return(s)
        })
        
        metadat_core.ind <- metadat_core.ind[colnames(sp.summ), ]
        load(file = 'Result/otu_cores/H_pylori.RData') # Hp
        metadat_core.ind$Hp = Hp[rownames(metadat_core.ind)]
        
        ## -
        ## - draw heat map
        if (TRUE) {
          if (TRUE) {
            wanted = c('g__Klebsiella;s__pneumoniae', 
                       'g__Dialister;s__pneumosintes', 
                       'g__Veillonella;s__parvula', 
                       'g__Streptococcus;s__salivarius', 
                       'g__Streptococcus;s__anginosus', 
                       'g__Parvimonas;s__micra', 
                       'g__Fusobacterium;s__nucleatum')
            
            taxa = taxa[taxa$Name %in% wanted, ]
            
            sp.summ = sp.summ[rownames(taxa), ]
          }
          sp.summ = sp.summ[order(rowSums(sp.summ != 0), decreasing = TRUE), ]
          
          Fun_row_anno <- rowAnnotation(
            Text = anno_text(paste(' ', taxa[rownames(sp.summ), 'Name'], " ", sep = ''), just = 'right', location = unit(1, 'npc'),
                             #rownames(taxonomy[rownames(sp.summ), ]),
                             gp = gpar(cex = 0.5, col = 'black', fontface = 'bold')),
            # status = taxa$status,
            col = list(status = c('up' = '#ff0000', 'down' = '#028900')),
            border = TRUE,
            gap = unit(2, "points"),
            simple_anno_size = unit(0.25, "cm"),
            gp = gpar(col = NA, width = 0.1),
            show_annotation_name = TRUE,
            annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                           nrow = 9,
                                           border = 'black',
                                           grid_width = unit(5, "mm"))
          )
          
          Fun_row_right <- rowAnnotation(
            Ratio = anno_barplot(rowSums(sp.summ != 0) / dim(sp.summ)[2], border = FALSE),
            border = FALSE,
            gap = unit(2, "points"),
            simple_anno_size = unit(0.25, "cm"),
            gp = gpar(col = NA, width = 0.1),
            show_annotation_name = TRUE,
            annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                           nrow = 9,
                                           border = 'black',
                                           grid_width = unit(5, "mm"))
          )
          
          Fun_col_anno <- HeatmapAnnotation(
            which = 'column',
            # ICU = as.character(metadat_core.ind$`ICU (Yes: 1; No: 0)`),
            # time = as.numeric(metadat_core.ind$`Sampling time after death  (Min)`),
            # Antibiotic = as.character(metadat_core.ind$`Antibiotic uses (Yes: 1; No: 0)`),
            # Hospitalization = as.numeric(metadat_core.ind$Hospitalization),
            BMI = as.numeric(metadat_core.ind$BMI),
            Age = as.numeric(metadat_core.ind$Age),
            Hp = as.numeric(metadat_core.ind$Hp),
            # Gender = metadat_core.ind$Gender,
            col = list(ICU = c('0' = 'white', '1' = 'black'),
                       time = colorRamp2(c(70, 100), c("white", "black")),
                       Antibiotic = c('0' = 'white', '1' = 'black'),
                       Hospitalization = colorRamp2(c(1, 10), c("white", "black")),
                       BMI = colorRamp2(c(0, 30), c("white", "black")),
                       Age = colorRamp2(c(0, 100), c("white", "black")),
                       Gender = c(Male = '#3399ff', Female = '#ff93ac'),
                       Hp = colorRamp2(c(0, 0.001, 1), c("white", "white", "black"))),
            border = TRUE,
            show_annotation_name = TRUE,
            annotation_name_gp = gpar(cex = 0.6),
            gap = unit(2, "points"),
            simple_anno_size = unit(0.2, "cm"), 
            # gp = gpar(col = "grey", width = 0.1),
            annotation_legend_param = list(legend_direction = "vertical",
                                           nrow = 9,
                                           border = 'black',
                                           grid_width = unit(5, "mm"))
          )
          
          
          hp <- Heatmap(matrix = as.matrix(sp.summ), 
                        col = colorRamp2(breaks = c(0, 2), colors = c("white", "red"), space = 'sRGB'),
                        # col = c('0' = 'white', '1' = 'red'),
                        
                        border = 'black', 
                        # rect_gp = gpar(color = '#f0f8ff', lwd = 0.01),
                        
                        cluster_rows = FALSE,
                        # row_order = ordering_row,
                        cluster_columns = FALSE,
                        # column_order = ordering_col,
                        show_row_dend = FALSE,
                        show_column_dend = FALSE,
                        
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        column_names_gp = gpar(cex = 0.6),
                        
                        split = factor(c(rep('A', 1), rep('B', 4), rep('C', 2)), levels = c('A', 'B', 'C')),
                        row_title = NULL,
                        row_gap = unit(0, "mm"),
                        
                        right_annotation = Fun_row_right,
                        left_annotation = Fun_row_anno,
                        top_annotation = Fun_col_anno,
                        # bottom_annotation = Fun_col_anno,
                        
                        show_heatmap_legend = TRUE, 
                        heatmap_legend_param = list(legend_direction = "ver", title = ''),
                        column_title = paste('', '', sep = ''))
          
          draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
               heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
          
          
          pdf(sprintf('Result/rare/core_bugs/Oral_gut/ind(with ID)_%s.pdf', type), width = 5.04, height = 5.5) # Oral_gut
          pdf(sprintf('Result/rare/core_bugs/Oral_gut/ind(with ID)_%s-.pdf', type), width = 5.04, height = 1.6) # Oral_gut
          
          draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
               heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
          dev.off()
        }
        
        ## -
        ## -
        if (TRUE) {
          taxa = taxa[rownames(sp.summ), ]
          sp.summ_val = sp.summ_val[rownames(sp.summ), colnames(sp.summ)]
          
          sp.tab = lapply(rownames(sp.summ_val), function(sp) {
            # sp = rownames(sp.summ_val)[1]
            tab = do.call(rbind, lapply(colnames(sp.summ_val), function(ind) {
              # ind = 'S001'
              idx = stringi::stri_split(sp.summ_val[sp, ind], regex = '//')[[1]]
              if (length(idx) == 1 && idx == "") {
                y = rep(NA, length(unique(metadat_core$site)))
              } else {
                y = c()
                for (region in unique(metadat_core$site)) {
                  # region = 'Appendix'
                  idx_r = which(metadat_core$subject == ind & metadat_core$site == region)
                  if (length(idx_r) > 1) {
                    y = c(y, sum(rowMeans(genus_rarefy[idx, rownames(metadat_core)[idx_r]])))
                  } else if (length(idx_r) == 1) {
                    y = c(y, sum(genus_rarefy[idx, rownames(metadat_core)[idx_r]]))
                  } else {
                    y = c(y, 0)
                  }
                }
                names(y) = unique(metadat_core$site)
              }
              y
            }))
            rownames(tab) = colnames(sp.summ_val)
            tab = as.data.frame(tab)
            
            tab
          })
          names(sp.tab) = rownames(sp.summ_val)
          
          save(sp.tab, file = 'Result/rare/core_bugs/Oral_gut/abun.sp_collapse_from_asv.RData')
          
          ## -
          load(file = 'Result/rare/core_bugs/Oral_gut/abun.sp_collapse_from_asv.RData')
          
          abun = sp.tab$`Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;pneumoniae`
          abun.m = reshape2::melt(as.matrix(abun))
          abun.m$Var2 = factor(abun.m$Var2, map_site_orgam)
          abun.m = abun.m[abun.m$Var2 != 'Appendix', ]
          
          abun.m = abun.m[!is.na(abun.m$value), ]
          abun.m$value = abun.m$value / 6000
          
        }
        
        ## -
        ## -
        if (TRUE) {
          taxa = taxa[rownames(sp.summ), ]
          species.rel = sweep(species, 2, apply(species, 2, sum), '/')
          
          sp.tab = lapply(rownames(taxa), function(sp) {
            # sp = rownames(taxa)[1]
            tab = do.call(rbind, lapply(colnames(sp.summ), function(ind) {
              # ind = 'S001'
              y = c()
              for (region in unique(metadat_core$site)) {
                # region = 'ANAL'
                idx_r = rownames(metadat_core)[which(metadat_core$subject == ind & metadat_core$site == region)]
                if (length(idx_r) > 1) {
                  y = c(y, mean(as.numeric(species.rel[sp, idx_r])))
                } else if (length(idx_r) == 1) {
                  y = c(y, as.numeric(species.rel[sp, idx_r]))
                } else {
                  y = c(y, 0)
                }
              }
              names(y) = unique(metadat_core$site)
              
              y
            }))
            rownames(tab) = colnames(sp.summ_val)
            tab = as.data.frame(tab)
            
            tab
          })
          names(sp.tab) = rownames(taxa)
          
          save(sp.tab, file = 'Result/rare/core_bugs/Oral_gut/abun.sp.RData')
          
          ## -
          load(file = 'Result/rare/core_bugs/Oral_gut/abun.sp.RData')
          
          abun = sp.tab$`Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;pneumoniae`
          abun.m = reshape2::melt(as.matrix(abun))
          abun.m$Var2 = factor(abun.m$Var2, map_site_orgam)
          abun.m = abun.m[abun.m$Var2 != 'Appendix', ]
          
          ## -
          if (FALSE) {
            c(wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("Oral cavity", "Cecum"), ], paired = TRUE)$p.value,
              wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("Cecum", "ASCM"), ], paired = TRUE)$p.value,
              wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("ASCM", "TRC"), ], paired = TRUE)$p.value,
              wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("TRC", "DECM"), ], paired = TRUE)$p.value,
              wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("DECM", "SICM"), ], paired = TRUE)$p.value,
              wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("SICM", "RECU"), ], paired = TRUE)$p.value,
              wilcox.test(value~Var2, data = abun.m[abun.m$Var2 %in% c("RECU", "ANAL"), ], paired = TRUE)$p.value)
            
            cor.test(abun$`Oral cavity`, abun$Cecum, method = 'spearman')$p.value
            cor.test(abun$`Oral cavity`, abun$ASCM, method = 'spearman')$p.value
            cor.test(abun$`Oral cavity`, abun$TRC, method = 'spearman')$p.value
            cor.test(abun$`Oral cavity`, abun$DECM, method = 'spearman')$p.value
            cor.test(abun$`Oral cavity`, abun$SICM, method = 'spearman')$p.value
            cor.test(abun$`Oral cavity`, abun$RECU, method = 'spearman')$p.value
            
            c(cor.test(abun$`Oral cavity`, abun$Cecum, method = 'spearman')$p.value,
              cor.test(abun$Cecum, abun$ASCM, method = 'spearman')$p.value,
              cor.test(abun$ASCM, abun$TRC, method = 'spearman')$p.value,
              cor.test(abun$TRC, abun$DECM, method = 'spearman')$p.value,
              cor.test(abun$DECM, abun$SICM, method = 'spearman')$p.value,
              cor.test(abun$SICM, abun$RECU, method = 'spearman')$p.value,
              cor.test(abun$RECU, abun$ANAL, method = 'spearman')$p.value)
            
            
            
            g <- ggplot(abun.m, aes(x = Var2, y = sqrt(asin(value)))) + 
              geom_boxplot(outlier.shape = NA, #outlier.size = 0.5, outlier.colour = '#d62d20',
                           position = position_dodge(0.5),
                           width = 0.3, lwd = 0.3, notch = FALSE) +
              stat_summary(fun.y = mean, geom = "point", shape = 3, size = 0.7, col = 'white',
                           position = position_dodge(0.5)) +
              geom_point(color = 'black', pch = 21, size = 1) +
              geom_line(aes(group = Var1, color = Var1), lty = 'dashed', lwd = 0.2) +
              labs(x = '', y = "Normalized Abundance")
            g <- g + theme_minimal() + 
              theme(panel.grid = element_line(color = 'white', linetype = 1), 
                    legend.position = 'none', #, c(1.0, 1.0), 
                    legend.spacing.x = unit(x = 6, units = 'pt'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8),
                    legend.key.size = unit(x = 8, units = 'pt'),
                    legend.background = element_rect(fill="white", size=.5, linetype="solid", colour = 'black'), 
                    axis.title = element_text(size = 8), 
                    axis.text.x = element_text(size = 8, color = 'black', angle = 60, vjust = 1, hjust = 1),
                    axis.text.y = element_text(size = 8, color = 'black'),
                    axis.line.x.bottom = element_line(colour = 'black'),
                    axis.ticks.x.bottom = element_line(colour = 'black'),
                    axis.line.y.left = element_line(colour = 'black'),
                    axis.ticks.y.left = element_line(colour = 'black'),
                    axis.ticks.length = unit(.05, "cm")) +
              scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
            g
            ggsave(filename = sprintf('Result/rare/core_bugs/Oral_gut/%s.pdf', 'K.pneumoniae'), plot = g, width = 3.79, height = 2.20)
          }
          
        }
        
        
        stop('')
      }
      
      
    }
  }
  
  
  ## - get core ASV at ind level
  if (FALSE) {
    map_site_orgam <- c(
      'NEG' = 'Negtive control'
      ,
      "SKINLP" = 'Body skin',
      "SKINUP" = 'Body skin'
      ,
      'ORALmix' = 'Oral cavity'
      ,
      'ESOM' = 'Esophageal',
      'ZA1' = 'Esophageal',
      'Z' = 'Esophageal',
      'ZB1' = 'Esophageal'
      ,
      'SF' = 'Stomach',
      'SB' = 'Stomach',
      "GJ" = 'Stomach',
      "SA" = 'Stomach',
      "PY" = 'Stomach'
      ,
      "DUS" = 'duodenum',
      'DUPA' = 'duodenum',
      'TRE' = 'duodenum',
      "JEJ100" = 'SI',
      "JEJ200" = 'SI',
      "IIE300" = 'SI',
      "IIE400" = 'SI',
      "IICP" = 'SI',
      "IIC" = 'SI'
      ,
      'AP' = 'Appendix'
      ,
      'Cecum' = 'LI',
      "ASCM" = 'LI',
      "TRC" = 'LI',
      "DECM" = 'LI',
      "SICM" = 'LI',
      "RECU" = 'LI',
      'ANAL' = 'LI'
    )
    
    metadat_core <- metadat
    metadat_core$site <- map_site_orgam[metadat_core$site]
    
    metadat_core <- metadat_core[!(metadat_core$site %in% c('Negtive control', 'Body skin')), ]
    table(metadat_core$site)
    
    ## - 
    genus_rarefy <- genus_rarefy[, rownames(metadat_core)]
    
    otu_cores_ind <- list()
    tmp <- lapply(X = unique(metadat_core$subject), FUN = function(ind) {
      print(ind)
      metadat_core.ind <- metadat_core[metadat_core$subject == ind, ]
      otu_core <- genus_rarefy[, rownames(metadat_core.ind)]
      otu_core <- otu_core[rowSums(otu_core) > 0, ]
      
      NUM_site = dim(metadat_core.ind)[1]
      cutoff = 1/1000 # abun > 0.1%, prev > 0
      depth = 6000
      Ind_site <- lapply(X = unique(metadat_core.ind$site), FUN = function(site) {
        idx <- which(metadat_core.ind$site == site)
        if (length(idx) >= 1) {
          if (length(idx) == 1) {
            otu_core.site <- as.numeric(otu_core[, rownames(metadat_core.ind)[idx] ]); names(otu_core.site) <- rownames(otu_core)
            otu_core.site = otu_core.site / depth
            otu_core.site <- otu_core.site[otu_core.site > cutoff]; 
            
            otu_core_summ <- data.frame(Mean = otu_core.site, 
                                        Prevalence = as.integer(otu_core.site > cutoff), 
                                        row.names = names(otu_core.site), 
                                        stringsAsFactors = FALSE)
            otu_core_summ$Index = (1 / NUM_site) * otu_core_summ$Prevalence
            otu_core_summ <- otu_core_summ[order(otu_core_summ$Prevalence, decreasing = TRUE), ]
          } else {
            otu_core.site <- otu_core[, rownames(metadat_core.ind)[idx] ]
            otu_core.site <- otu_core.site[rowSums(otu_core.site) > 0, ]
            
            pseq <- phyloseq(otu_table(otu_core.site, taxa_are_rows = TRUE))
            pseq.comp <- microbiome::transform(pseq, "compositional")
            taxa <- core_members(pseq.comp, detection = cutoff, prevalence = 0); 
            otu_core.site <- prune_taxa(taxa, pseq)@.Data
            
            otu_core_summ <- apply(X = otu_core.site, MARGIN = 1, FUN = function(otu) {
              otu <- otu / depth
              c('Mean' = mean(otu), 'Prevalence' = sum(otu > cutoff) / length(otu))
            })
            otu_core_summ <- as.data.frame(t(otu_core_summ))
            otu_core_summ$Mean <- as.numeric(as.character(otu_core_summ$Mean))
            otu_core_summ$Prevalence <- as.numeric(as.character(otu_core_summ$Prevalence))
            otu_core_summ$Index = (length(idx) / NUM_site) * otu_core_summ$Prevalence
            otu_core_summ <- otu_core_summ[order(otu_core_summ$Prevalence, decreasing = TRUE), ]
          }
        } else {
          return(NA)
        }
      })
      names(Ind_site) <- unique(metadat_core.ind$site)
      
      taxa <- Reduce(union, lapply(Ind_site, function(site) {
        rownames(site)
      }))
      otu_core_summ <- matrix(data = NA, nrow = length(taxa), ncol = length(Ind_site), 
                              dimnames = list(taxa, names(Ind_site)))
      tmp <- sapply(names(Ind_site), function(n) {
        # print(n)
        site <- Ind_site[[n]]
        otu_core_summ[rownames(site), n] <<- as.numeric(site$Mean) # Index, Prevalence, Mean
        return(NA)
      })
      
      taxa <- taxonomy[rownames(otu_core_summ), ]
      otu_core_summ <- cbind(taxa, otu_core_summ)
      write.csv(otu_core_summ, file = sprintf('Result/otu_cores/Ind/Ind-%s.csv', ind))
      
      otu_cores_ind[[ind]] <<- otu_core_summ
    })
    # save(otu_cores_ind, file = 'Result/otu_cores/Ind/otu_cores_ind.Rds')
    
    
    ##-
    load(file = 'Result/otu_cores/Ind/otu_cores_ind.Rds') # otu_cores_ind
    
    metadat_core.ind <- read.csv(file = 'metadata/patient_info-new.tsv', header = TRUE, sep = '\t', row.names = 1, 
                                 check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
    metadat_core.ind$Hospitalization[grepl('^>', metadat_core.ind$Hospitalization)] <- 10
    
    
    ##-
    ## - 
    if (TRUE) {
      cutoff = 0.001 # prev: 0; abun: 0.001
      Ind_summ.UD <- list()
      Ind_summ.LD <- list()
      Ind_summ.Oral_gut <- list()
      Ind_summ.Uniq <- list()
      Ind_summ.All <- list()
      tmp <- lapply(X = names(otu_cores_ind), FUN = function(n) {
        ind = otu_cores_ind[[n]]
        taxa <- ind[, colnames(taxonomy)]
        taxa <- taxa[taxa$Kingdom != 'k__Archaea', ]
        ind <- ind[rownames(taxa), ]
        
        ID_sites <- list(skin = rownames(ind)[!is.na(ind$`Body skin`) & ind$`Body skin` >= cutoff],
                         oral = rownames(ind)[!is.na(ind$`Oral cavity`) & ind$`Oral cavity` >= cutoff],
                         esophagus = rownames(ind)[!is.na(ind$Esophageal) & ind$Esophageal >= cutoff],
                         stomach = rownames(ind)[!is.na(ind$Stomach) & ind$Stomach >= cutoff],
                         duodenum = rownames(ind)[!is.na(ind$duodenum) & ind$duodenum >= cutoff],
                         SI = rownames(ind)[!is.na(ind$SI) & ind$SI >= cutoff],
                         appendix = rownames(ind)[!is.na(ind$Appendix) & ind$Appendix >= cutoff],
                         LI = rownames(ind)[!is.na(ind$LI) & ind$LI >= cutoff] )
        
        ID = Reduce(union, ID_sites)
        
        summ <- matrix(data = 0, nrow = length(ID), ncol = length(ID_sites), 
                       dimnames = list(ID, names(ID_sites)))
        tmp <- sapply(names(ID_sites), function(n) {
          # print(n)
          summ[ID_sites[[n]], n] <<- 1
          return(NA)
        })
        summ <- as.data.frame(summ)
        
        Ind_summ.All[[n]] <<- summ
        
        ##
        stat <- rowSums(summ)
        # UD <- rowSums(summ[, c("oral", "esophagus", "stomach", "duodenum")])
        # LD <- rowSums(summ[, c("SI", "appendix", "LI")])
        UD <- rowSums(summ[, c("esophagus", "stomach")])
        LD <- rowSums(summ[, c("duodenum", "SI", "appendix", "LI")])
        Oral_gut <- rowSums(summ[, c("oral", "LI")])
        
        summ$stat <- stat
        summ$UD <- UD
        summ$LD <- LD
        summ$Oral_gut <- Oral_gut
        
        # summ = summ[summ$UD == 4 | summ$LD == 3 | summ$Oral_gut == 2, ]
        Ind_summ.Uniq[[n]] <<- summ[summ$stat == 1, ]
        
        Ind_summ.UD[[n]] <<- summ[summ$UD == 2, ]
        
        if (sum(summ$appendix) == 0) 
          Ind_summ.LD[[n]] <<- summ[summ$LD >= 2, ] 
        else 
          Ind_summ.LD[[n]] <<- summ[summ$LD >= 3, ]
        
        Ind_summ.Oral_gut[[n]] <<- summ[summ$Oral_gut == 2, ]
        
        ## at least present in at least 3 organs
        # return(summ[summ$stat >= 3, ])
        return(NA)
      })
      # save(Ind_summ.LD, Ind_summ.UD, Ind_summ.Oral_gut, 
      #      file = 'Result/rare/core_bugs/Ind/Ind_summ.RData')
      # save(Ind_summ.Uniq, 
      #      file = 'Result/rare/core_bugs/Organ.Uniq/Ind_summ.Uniq.RData')
      # save(Ind_summ.All,
      #      file = 'Result/rare/core_bugs/All/Ind_summ.All.RData')
      
      
      ## -. presence_of_asv_among_inds: Ind_summ.All
      if (TRUE) {
        load(file = 'Result/rare/core_bugs/All/Ind_summ.All.RData')
        
        Ind_summ = Ind_summ.All
        type = 'All'
        ##
        ID <- Reduce(union, lapply(Ind_summ, function(x) rownames(x)))
        summ <- matrix(data = '', nrow = length(ID), ncol = length(Ind_summ), 
                       dimnames = list(ID, names(Ind_summ)))
        tmp <- sapply(names(Ind_summ), function(n) {
          # print(n)
          ind = Ind_summ[[n]]
          summ[rownames(ind), n] <<- apply(ind, 1, function(x) {
            organs = colnames(ind)[which(x != 0)]
            
            is_UD = all(c("esophagus", "stomach") %in% organs)
            if (sum(ind$appendix) == 0) {
              is_LD = sum(c("duodenum", "SI", "LI") %in% organs) >= 2
            } else {
              is_LD = sum(c("duodenum", "SI", 'appendix', "LI") %in% organs) >= 3
            }
            is_OG = all(c("oral", "LI") %in% organs)
            
            if (is_UD) organs = c(organs, 'UD')
            if (is_LD) organs = c(organs, 'LD')
            if (is_OG) organs = c(organs, 'OG')
            
            paste(organs, collapse = ';')
            
          })
          return(NA)
        })
        summ <- as.data.frame(summ)
        
        ## taxa
        if (TRUE) {
          taxa = taxonomy[rownames(summ), c(names(mapping), 'Species')]
          colnames(taxa) = c(mapping, 's')
          taxa[is.na(taxa)] <- ''
          taxa = do.call(rbind, lapply(1:dim(taxa)[1], function(i) {
            paste(colnames(taxa), as.character(taxa[i, ]), sep = '__')
          }))
          taxa = as.data.frame(taxa, stringsAsFactors = FALSE)
          colnames(taxa) = c(mapping, 's')
          rownames(taxa) = rownames(summ)
          
          Name = apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
            g <- get_lv(x = x, usedLv = 'g')
            s <- get_lv(x = x, usedLv = 's')
            if (startsWith(s, 's__')) {
              s <- paste(g, s, sep = ';')
            }
            return(s)
          })
          taxa$Name = paste(sapply(names(Name), function(n) {
            stringi::stri_split(n, regex = ';')[[1]][1]
          }), Name, sep = '-')
        }
        # save(summ, taxa, file = 'Result/rare/core_bugs/ASV/Organ.All(presence_of_asv_among_inds).RData')
        
        
        ## -- 
        asv.sp = 'g__Bacteroides.s__fragilis'
        asv = summ[rownames(taxa)[taxa$g == 'g__Bacteroides' & taxa$s == 's__fragilis'], ]
        # View(asv)
        
        organs = c("oral", "esophagus", "stomach", "duodenum", "SI", "appendix", "LI", 'UD', 'LD', 'OG')
        asv.organs <- list()
        asv.stat = do.call('rbind', lapply(organs, function(organ) {
          # organ = 'LD' 
          tmp = asv[apply(asv, 1, function(x) { # the asvs belonging to the `sp` in `organ`
            x = as.character(x)
            any(grepl(pattern = organ, x = x))
          }), ]
          asv.organs[[organ]] <<- rownames(tmp)
          tmp.prev = apply(tmp, 2, function(x) { # the prevalence of that set of asvs among population
            x = as.character(x)
            sum(grepl(pattern = organ, x = x))
          })
          # c("# asv" = dim(tmp)[1], "prev" = sum(tmp.prev != 0))
          # return(c(tmp.prev, "# asv" = dim(tmp)[1]))
          return(tmp.prev)
        }) )
        rownames(asv.stat) = organs
        
        ### draw heatmap
        if (TRUE) {
          asv.stat = asv.stat[, sort(colnames(asv.stat))]
          metadat_core.ind <- metadat_core.ind[colnames(asv.stat), ]
          
          Fun_row_anno <- rowAnnotation(
            Text = anno_text(paste(' ', rownames(asv.stat), " ", sep = ''), just = 'right', location = unit(1, 'npc'),
                             gp = gpar(cex = 0.5, col = 'black', fontface = 'bold')),
            # status = taxa$status,
            col = list(oral = colorRamp2(c(0, 33), c("white", "black")),
                       status = c('up' = '#ff0000', 'down' = '#028900')),
            border = TRUE,
            gap = unit(2, "points"),
            simple_anno_size = unit(0.25, "cm"),
            gp = gpar(col = NA, width = 0.1),
            show_annotation_name = TRUE,
            annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                           nrow = 9,
                                           border = 'black',
                                           grid_width = unit(5, "mm"))
          )
          
          Fun_row_right <- rowAnnotation(
            `# asv` = anno_barplot(sapply(asv.organs, length), border = FALSE, gp = gpar(fill = 'red')),
            border = FALSE,
            gap = unit(2, "points"),
            simple_anno_size = unit(0.25, "cm"),
            gp = gpar(col = NA, width = 0.1),
            show_annotation_name = TRUE,
            annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                           nrow = 9,
                                           border = 'black',
                                           grid_width = unit(5, "mm"))
          )
          
          Fun_col_anno <- HeatmapAnnotation(
            which = 'column',
            # ICU = as.character(metadat_core.ind$`ICU (Yes: 1; No: 0)`),
            # time = as.numeric(metadat_core.ind$`Sampling time after death  (Min)`),
            # Antibiotic = as.character(metadat_core.ind$`Antibiotic uses (Yes: 1; No: 0)`),
            # Hospitalization = as.numeric(metadat_core.ind$Hospitalization),
            BMI = as.numeric(metadat_core.ind$BMI),
            Age = as.numeric(metadat_core.ind$Age),
            # Gender = metadat_core.ind$Gender,
            col = list(ICU = c('0' = 'white', '1' = 'black'),
                       time = colorRamp2(c(70, 100), c("white", "black")),
                       Antibiotic = c('0' = 'white', '1' = 'black'),
                       Hospitalization = colorRamp2(c(1, 10), c("white", "black")),
                       BMI = colorRamp2(c(0, 30), c("white", "black")),
                       Age = colorRamp2(c(0, 100), c("white", "black")),
                       Gender = c(Male = '#3399ff', Female = '#ff93ac')),
            border = TRUE,
            show_annotation_name = TRUE,
            annotation_name_gp = gpar(cex = 0.6),
            gap = unit(2, "points"),
            simple_anno_size = unit(0.2, "cm"), 
            # gp = gpar(col = "grey", width = 0.1),
            annotation_legend_param = list(legend_direction = "vertical",
                                           nrow = 9,
                                           border = 'black',
                                           grid_width = unit(5, "mm"))
          )
          
          
          hp <- Heatmap(matrix = as.matrix(asv.stat), 
                        col = colorRamp2(breaks = c(0, 2), colors = c("white", "red"), space = 'sRGB'),
                        # col = c('0' = 'white', '1' = 'red'),
                        
                        border = 'black', 
                        # rect_gp = gpar(color = '#f0f8ff', lwd = 0.01),
                        
                        cluster_rows = FALSE,
                        # row_order = ordering_row,
                        cluster_columns = FALSE,
                        # column_order = ordering_col,
                        show_row_dend = FALSE,
                        show_column_dend = FALSE,
                        
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        column_names_gp = gpar(cex = 0.6),
                        
                        split = c(rep(1, 7), 2, 3, 4),
                        row_title = NULL,
                        row_gap = unit(0, "mm"),
                        
                        right_annotation = Fun_row_right,
                        left_annotation = Fun_row_anno,
                        # top_annotation = Fun_col_anno,
                        # bottom_annotation = Fun_col_anno,
                        
                        show_heatmap_legend = TRUE, 
                        heatmap_legend_param = list(legend_direction = "ver", title = ''),
                        column_title = sprintf('Total ASV of %s: %s', asv.sp, dim(asv)[1]))
          
          draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
               heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
          
          
          pdf(sprintf('Result/rare/core_bugs/ASV/%s.pdf', asv.sp), width = 3.7, height = 1.50)
          
          draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
               heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
          dev.off()
        }
        
        ### got the sequences
        if (TRUE) {
          asv.seqs = query_seq[rownames(asv)]
          seqinr::write.fasta(sequences = as.list(asv.seqs), names = as.list(names(asv.seqs)), 
                              file.out = sprintf('Result/rare/core_bugs/ASV/%s.fasta', asv.sp))
          
          organ_col = c('oral' = '#008000',
                        'esophagus' = '#830303',
                        'stomach' = '#B967FF',
                        'duodenum' = '#7f7fff',
                        'SI' = '#0000FF',
                        'LI' = '#C51B7D',
                        'appendix' = '#00ced1',
                        'LI' = 'black',
                        'UD' = '#008000',
                        'LD' = '#830303',
                        'OG' = '#B967FF')
          asv.seqs_anno = matrix(data = '#FFFFFF', nrow = length(asv.seqs), ncol = length(organs)+1, 
                                 dimnames = list(names(asv.seqs), c(organs, 'prev')))
          tmp = sapply(names(asv.organs), function(n) {
            organ = asv.organs[[n]]
            if (length(organ) > 0) {
              asv.seqs_anno[organ, n] <<- organ_col[n]
            }
            return(NA)
          })
          asv.seqs_anno = as.data.frame(asv.seqs_anno)
          asv.seqs_anno$prev = apply(asv[rownames(asv.seqs_anno), ], 1, function(x) sum(x != ''))
          write.csv(asv.seqs_anno, file = sprintf('Result/rare/core_bugs/ASV/%s.anno.csv', asv.sp))
        }
        
        
        stop('')
      }
      
      
      ## - 01. organ specific and All(data)
      if (TRUE) {
        load(file = 'Result/rare/core_bugs/All/Ind_summ.All.RData')
        
        Ind_summ = Ind_summ.All # Ind_summ.Uniq; Ind_summ.All
        type = 'All' # # '.'; 'All'
        ##
        Organ.Uniq <- list()
        for (organ in c("oral", "esophagus", "stomach", "duodenum", "SI", "appendix", "LI")) {
          # organ = 'oral' # "oral", "esophagus", "stomach", "duodenum", "SI", "appendix", "LI"
          print(organ)
          ID <- Reduce(union, lapply(Ind_summ, function(x) rownames(x)))
          summ <- matrix(data = 0, nrow = length(ID), ncol = length(Ind_summ), 
                         dimnames = list(ID, names(Ind_summ)))
          tmp <- sapply(names(Ind_summ), function(n) {
            # print(n)
            ind = Ind_summ[[n]]
            summ[rownames(ind)[ind[, organ] == 1], n] <<- 1
            return(NA)
          })
          summ <- as.data.frame(summ)
          summ = summ[rowSums(summ) > 0, ]
          
          ## ASV
          if (FALSE) {
            summ = summ[order(rowSums(summ != 0), decreasing = TRUE), ]
            
            taxa = taxonomy[rownames(summ), c(names(mapping), 'Species')]
            colnames(taxa) = c(mapping, 's')
            taxa[is.na(taxa)] <- ''
            taxa = do.call(rbind, lapply(1:dim(taxa)[1], function(i) {
              paste(colnames(taxa), as.character(taxa[i, ]), sep = '__')
            }))
            taxa = as.data.frame(taxa, stringsAsFactors = FALSE)
            colnames(taxa) = c(mapping, 's')
            rownames(taxa) = rownames(summ)
            
            Name = apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
              g <- get_lv(x = x, usedLv = 'g')
              s <- get_lv(x = x, usedLv = 's')
              if (startsWith(s, 's__')) {
                s <- paste(g, s, sep = ';')
              }
              return(s)
            })
            taxa$Name = paste(sapply(names(Name), function(n) {
              stringi::stri_split(n, regex = ';')[[1]][1]
            }), Name, sep = '-')
            taxa$stat = rowSums(summ != 0)
          }
          
          ## ASV collapse to sp
          if (TRUE) {
            taxa = taxonomy[rownames(summ), c(names(mapping), 'Species')]
            colnames(taxa) = c(mapping, 's')
            taxa[is.na(taxa)] <- ''
            taxa = do.call(rbind, lapply(1:dim(taxa)[1], function(i) {
              paste(colnames(taxa), as.character(taxa[i, ]), sep = '__')
            }))
            taxa = as.data.frame(taxa, stringsAsFactors = FALSE)
            colnames(taxa) = c(mapping, 's')
            rownames(taxa) = rownames(summ)
            
            taxa$Name = apply(taxa, 1, paste, collapse = ';')
            tmp.s = unique(taxa$Name)
            
            sp.summ <- matrix(data = 0, nrow = length(tmp.s), ncol = dim(summ)[2], 
                              dimnames = list(tmp.s, colnames(summ)))
            tmp <- sapply(tmp.s, function(s) {
              idx = rownames(taxa)[taxa$Name == s]
              if (length(idx) == 1) {
                sp.summ[s, ] <<- as.numeric(summ[idx, ])
              } else if (length(idx) > 1) {
                sp.summ[s, ] <<- colSums(summ[idx, ])
              }
              return(NA)
            })
            sp.summ <- as.data.frame(sp.summ)
            sp.summ = sp.summ[order(rowSums(sp.summ != 0), decreasing = TRUE), ]
            
            
            ##
            taxa = as.data.frame(taxa_chr2df_sp(rownames(sp.summ)), stringsAsFactors = FALSE)
            taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
              g <- get_lv(x = x, usedLv = 'g')
              s <- get_lv(x = x, usedLv = 's')
              if (startsWith(s, 's__')) {
                s <- paste(g, s, sep = ';')
              }
              return(s)
            })
            taxa$stat = rowSums(sp.summ != 0)
          }
          
          Organ.Uniq[[organ]] <- taxa
        }
        # save(Organ.Uniq, file = 'Result/rare/core_bugs/Organ.Uniq/Organ.Uniq(asv).RData')
        # save(Organ.Uniq, file = 'Result/rare/core_bugs/Organ.Uniq/Organ.Uniq(asv_collapse_to_sp).RData')
        # save(Organ.Uniq, file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).RData')
        
        ##
        load(file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).RData')
        
        ID = do.call(c, lapply(Organ.Uniq, function(x) {
          x$Name[!(x$s %in% c('', 's__')) & x$stat >= 16]
        }))
        ID = unname(ID)
        ID = unique(ID)
        summ = matrix(data = 0, nrow = length(ID), ncol = length(Organ.Uniq), 
                      dimnames = list(ID, names(Organ.Uniq) ))
        
        tmp = sapply(names(Organ.Uniq), function(organ) {
          tab = Organ.Uniq[[organ]]
          idx = match(ID, tab$Name)
          if (sum(!is.na(idx)) > 0) summ[which(!is.na(idx)), organ] <<- tab[idx[!is.na(idx)], 'stat']
          
          return(NA)
        })
        
        # write.csv(x = cbind(summ >= 16, stat = rowSums(summ >= 16)), 
        #           file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).csv')
        # write.csv(x = summ, 
        #           file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp)-ori.csv')
        
        
        ## From here
        # summ = summ[, colnames(summ) != 'duodenum']
        if (FALSE) {
          unwanted = c('g__Pseudomonas;s__putida', 
                       'g__Acinetobacter;s__johnsonii', 
                       'g__Bacillus;s__anthracis', 
                       'g__Photobacterium;s__damselae', 
                       'g__Lysinibacillus;s__xylanilyticus', 
                       'g__Psychrobacter;s__celer', 
                       'g__Pantoea;s__vagans', 
                       'g__Kluyvera;s__intermedia')
          summ = summ[!(rownames(summ) %in% unwanted), ]
          
          ID = rownames(summ)
          # write.table(x = ID, file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).ID.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
          # save(ID, file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).ID.RData')
          
          ##
          if (FALSE) {
            comb_cor = as.data.frame(t(combn(colnames(summ), 2)))
            comb_cor$status = 0
            comb_cor$status = apply(comb_cor, 1, function(ar) {
              x = summ[, ar[1:2]]
              sum(apply(x, 1, function(i) all(i>=16)))
            })
            
            comb_cor.r = matrix(data = NA, nrow = 7, ncol = 7, 
                                dimnames = list(colnames(summ), colnames(summ)))
            comb_cor.r = as.data.frame(comb_cor.r)
            tmp = apply(comb_cor, 1, function(x) {
              comb_cor.r[x[1], x[2]] <<- as.numeric(x[3]) / dim(summ)[1]
              comb_cor.r[x[2], x[1]] <<- as.numeric(x[3]) / dim(summ)[1]
              
              return(NA)
            })
            
            pdf('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).ratio_in_ID.pdf', width = 4.5, height = 4.5)
            M_order <- corrplot(corr = as.matrix(comb_cor.r), method="square", type = 'lower', 
                                addCoef.col = 'black', number.cex = 1,
                                col = colorRampPalette(c('#0057e7', "#f5f5f5", "#d62d20"))(20), is.corr=TRUE,
                                diag = FALSE, bg="#fbfdfb", addgrid.col = 'black',
                                tl.pos = 'l', tl.col="black", tl.srt = 30, tl.cex = 0.8, tl.offset = 0.5,
                                cl.pos = 'n', 
                                # p.mat = as.matrix(heatdata.p), sig.level = 0.06, insig = 'label_sig', pch.cex = 1.2, pch.col = '#ffa700',
                                mar = c(1, 1, 1, 1))
            dev.off()
            
          }
        }
        
        Fun_row_anno <- rowAnnotation(
          Text = anno_text(paste(' ', rownames(summ), " ", sep = ''), just = 'right', location = unit(1, 'npc'),
                           gp = gpar(cex = 0.5, col = 'black', fontface = 'bold')),
          border = TRUE,
          gap = unit(2, "points"),
          simple_anno_size = unit(0.25, "cm"),
          gp = gpar(col = NA, width = 0.1),
          show_annotation_name = TRUE,
          annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                         nrow = 9,
                                         border = 'black',
                                         grid_width = unit(5, "mm"))
        )
        
        
        hp <- Heatmap(matrix = as.matrix(summ), 
                      col = colorRamp2(breaks = c(0, 33), colors = c("white", "blue"), space = 'sRGB'),
                      # col = c('0' = 'white', '1' = 'red'),
                      
                      border = 'black', 
                      # rect_gp = gpar(color = '#f0f8ff', lwd = 0.01),
                      
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        if (summ[i, j] >= 16) {
                          grid.points(x, y, pch = 18, size = unit(2, "mm"), gp = gpar(fill = "red", col = "red"))
                        }
                      },
                      
                      cluster_rows = FALSE,
                      # row_order = ordering_row,
                      cluster_columns = FALSE,
                      # column_order = ordering_col,
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                      
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      column_names_gp = gpar(cex = 0.6),
                      
                      left_annotation = Fun_row_anno,
                      
                      show_heatmap_legend = TRUE, 
                      heatmap_legend_param = list(legend_direction = "ver", title = ''),
                      column_title = paste('', '', sep = ''))
        
        draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
             heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        
        
        # pdf(sprintf('Result/rare/core_bugs/Organ.Uniq/Organ.Uniq(asv).pdf'), width = 3.50, height = 4.15)  # dir: Organ.Uniq
        pdf(sprintf('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pdf'), width = 3.50, height = 8) # dir: All
        
        draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
             heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
      }
      
      
      ## - 01. organ All(data) for corr. analysis
      if (TRUE) {
        load(file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).ID.RData') # ID
        
        ## -> species
        taxa = as.data.frame(do.call('rbind', lapply(rownames(species), function(x) {
          x = stringi::stri_split(x, regex = ';')[[1]]
          paste(c(mapping,"s"), x, sep = '__')
        })), stringsAsFactors = FALSE)
        colnames(taxa) = c(mapping, 's')
        rownames(taxa) = rownames(species)
        taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
          # f <- get_lv(x = x, usedLv = 'f')
          g <- get_lv(x = x, usedLv = 'g')
          s <- get_lv(x = x, usedLv = 's')
          if (startsWith(s, 's__')) {
            s <- paste(g, s, sep = ';')
          }
          # return(paste(f, s, sep = ';'))
          return(s)
        })
        taxa = taxa[taxa$Name %in% ID, ]
        
        ## -
        species.rel = sweep(species, 2, apply(species, 2, sum), '/')
        sp.cor = list()
        sp.tab = lapply(rownames(taxa), function(sp) {
          # sp = rownames(taxa)[1]
          inds = sort(unique(metadat_core$subject))
          tab = do.call(rbind, lapply(inds, function(ind) {
            # ind = 'S001'
            y = do.call(c, lapply(unique(metadat_core$site), function(region) {
              idx_r = rownames(metadat_core)[which(metadat_core$subject == ind & metadat_core$site == region)]
              if (length(idx_r) > 1) {
                return(mean(as.numeric(species.rel[sp, idx_r]))) # mean, median
              } else if (length(idx_r) == 1) {
                return(as.numeric(species.rel[sp, idx_r]))
              } else {
                return(0)
              }
            }))
            names(y) = unique(metadat_core$site)
            
            y
          }))
          rownames(tab) = inds
          tab = as.data.frame(tab)
          
          
          ##
          comps <- as.data.frame(t(combn(colnames(tab), 2)))
          tmp = do.call('rbind', lapply(1:dim(comps)[1], function(idx) {
            region = as.character(comps[idx, ])
            
            tmp <- cor.test(tab[, region[1]], tab[, region[2]], method = 'spearman')
            rd <- c(sp.pval = tmp$p.value, sp = tmp$estimate)
            
            tmp <- cor.test(tab[, region[1]], tab[, region[2]], method = 'pearson')
            rd <- c(rd, pr.pval = tmp$p.value, pr = tmp$estimate)
            
            tmp = cbind(tab, metadat_core.ind[rownames(tab), ])
            tmp$r1 = tmp[, region[1]]
            tmp$r2 = tmp[, region[2]]
            tag <- tryCatch(
              # Specifying expression
              expr = {                     
                fit <- PResiduals::partial_Spearman(formula = r1|r2~Age+BMI, 
                                                    data = tmp, 
                                                    fit.x = 'orm', fit.y = 'orm')
                c(psp.pval = fit$TS$TB$pval, psp = fit$TS$TB$ts)
              },
              # Specifying error message
              error = function(e){         
                c(psp.pval = NA, psp = NA)
              },
              warning = function(w){      
                c(psp.pval = NA, psp = NA)
              }
            )
            rd <- c(rd, tag)
            rd['sp.pval'] <- sign(rd['sp.rho'])*rd['sp.pval']
            rd['pr.pval'] <- sign(rd['pr.cor'])*rd['pr.pval']
            rd['psp.pval'] <- sign(rd['psp'])*rd['psp.pval']
            
            rd
          }))
          comps = cbind(comps, tmp)
          comps$status = apply(comps[,grep('pval$', colnames(comps), value = TRUE)], 1, function(x) sum(!is.na(x) & x > 0 & x < 0.05))
          
          sp.cor[[sp]] <<- comps
          
          ## return the tab
          tab
        })
        names(sp.tab) = rownames(taxa)
        # save(sp.cor, sp.tab, file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr.RData')
        
        
        load(file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr.RData')
        ## - corrplot for overall
        organs = c("Oral cavity", "Esophageal", "Stomach", "duodenum", "SI", "Appendix", 'LI')
        if (TRUE) {
          load(file = 'Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr.RData') # sp.cor, sp.tab
          
          ## -
          comb_cor = as.data.frame(t(combn(unique(metadat_core$site), 2)))
          comb_cor$status = 0
          tmp = lapply(sp.cor, function(x) {
            comb_cor$status <<- comb_cor$status + as.integer(x$status > 0)
            
            return(NA)
          })
          
          ## -
          comb_cor.r = matrix(data = NA, nrow = 7, ncol = 7, 
                              dimnames = list(organs, organs))
          comb_cor.r = as.data.frame(comb_cor.r)
          tmp = apply(comb_cor, 1, function(x) {
            comb_cor.r[x[1], x[2]] <<- as.numeric(x[3]) / length(sp.cor)
            comb_cor.r[x[2], x[1]] <<- as.numeric(x[3]) / length(sp.cor)
            
            return(NA)
          })
          
          pdf('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr_ratio_in_ID.pdf', width = 4.5, height = 4.5)
          # jpeg('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr_ratio_in_ID.jpg', width = 400, height = 400)
          M_order <- corrplot(corr = as.matrix(comb_cor.r), method="square", type = 'lower', 
                              addCoef.col = 'black', number.cex = 1,
                              col = colorRampPalette(c('#0057e7', "#f5f5f5", "#d62d20"))(20), is.corr=TRUE,
                              diag = FALSE, bg="#fbfdfb", addgrid.col = 'black',
                              tl.pos = 'l', tl.col="black", tl.srt = 30, tl.cex = 0.8, tl.offset = 0.5,
                              cl.pos = 'n', 
                              # p.mat = as.matrix(heatdata.p), sig.level = 0.06, insig = 'label_sig', pch.cex = 1.2, pch.col = '#ffa700',
                              mar = c(1, 1, 1, 1))
          dev.off()
        }
        
        
        ## - corrplot and abundance change of each sp
        if (TRUE) {
          sp = 'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;pneumoniae'
          
          ## - corrplot for each sp
          for (sp in names(sp.cor)) {
            if (sum(sp.cor[[sp]]$status) > 0) {
              comb_cor.r = matrix(data = NA, nrow = 7, ncol = 7, 
                                  dimnames = list(organs, organs))
              comb_cor.r = as.data.frame(comb_cor.r)
              tmp = apply(sp.cor[[sp]], 1, function(x) {
                if (as.numeric(x['status']) > 0) {
                  # print(as.numeric(x['status']))
                  val = as.numeric(x[c('sp.rho', 'pr.cor', 'psp')])[which.min(abs(as.numeric(x[c('sp.pval', 'pr.pval', 'psp.pval')])) )]
                } else {
                  val = 0
                }
                
                comb_cor.r[x[1], x[2]] <<- val
                comb_cor.r[x[2], x[1]] <<- val
                
                return(NA)
              })
              
              pdf(sprintf('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr_of_%s.pdf', taxa[sp, 'Name']), width = 4.5, height = 4.5)
              # jpeg('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).pos_corr_ratio_in_ID.jpg', width = 400, height = 400)
              M_order <- corrplot(corr = as.matrix(comb_cor.r), method="ellipse", type = 'lower', 
                                  # addCoef.col = 'black', number.cex = 1,
                                  col = colorRampPalette(c('#0057e7', "#f5f5f5", "#d62d20"))(20), is.corr=TRUE,
                                  diag = FALSE, bg="#fbfdfb", addgrid.col = 'black',
                                  tl.pos = 'l', tl.col="black", tl.srt = 30, tl.cex = 0.8, tl.offset = 0.5,
                                  cl.pos = 'n', 
                                  p.mat = -as.matrix(comb_cor.r), sig.level = 0, insig = 'label_sig', pch.cex = 1.5, pch.col = '#0057e7',
                                  mar = c(1, 1, 1, 1))
              dev.off()
            }
          }
          
          ## - abundance change of each sp
          if (TRUE) {
            abun = sp.tab[[sp]]
            abun.m = reshape2::melt(as.matrix(abun))
            abun.m$Var2 = factor(abun.m$Var2, organs)
            abun.m = abun.m[abun.m$Var2 != 'Appendix', ]
            abun.m = abun.m[abun.m$Var2 != 'duodenum', ]
            # abun.m = abun.m[abun.m$Var1 != 'S007', ]
            
            g <- ggplot(abun.m, aes(x = Var2, y = sqrt(asin(value)))) + 
              geom_hline(yintercept = sqrt(asin(0.001)), lty = 'dashed', color = 'grey') +
              geom_boxplot(outlier.shape = NA, #outlier.size = 0.5, outlier.colour = '#d62d20',
                           position = position_dodge(0.5),
                           width = 0.3, lwd = 0.3, notch = FALSE) +
              stat_summary(fun.y = mean, geom = "point", shape = 3, size = 0.7, col = 'white',
                           position = position_dodge(0.5)) +
              geom_point(color = 'black', pch = 21, size = 1) +
              # geom_line(aes(group = Var1, color = Var1), data = abun.m[abun.m$Var2 %in% c("Oral cavity", "Esophageal", 'Stomach', 'SI', 'LI'), ],
              #           lty = 'dashed', lwd = 0.2) +
              geom_line(aes(group = Var1, color = Var1), data = abun.m[abun.m$Var2 %in% c("Esophageal", 'Stomach', 'SI', 'LI'), ],
                        lty = 'dashed', lwd = 0.2) +
              labs(x = '', y = "Normalized Abundance")
            g <- g + theme_minimal() + 
              theme(panel.grid = element_line(color = 'white', linetype = 1), 
                    legend.position = 'none', #, c(1.0, 1.0), 
                    legend.spacing.x = unit(x = 6, units = 'pt'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8),
                    legend.key.size = unit(x = 8, units = 'pt'),
                    legend.background = element_rect(fill="white", size=.5, linetype="solid", colour = 'black'), 
                    axis.title = element_text(size = 8), 
                    axis.text.x = element_text(size = 8, color = 'black', angle = 60, vjust = 1, hjust = 1),
                    axis.text.y = element_text(size = 8, color = 'black'),
                    axis.line.x.bottom = element_line(colour = 'black'),
                    axis.ticks.x.bottom = element_line(colour = 'black'),
                    axis.line.y.left = element_line(colour = 'black'),
                    axis.ticks.y.left = element_line(colour = 'black'),
                    axis.ticks.length = unit(.05, "cm")) +
              scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
              # scale_y_continuous(limits = c(0,0.4), expand = expansion(mult = c(0, 0.1)))
            g
            ggsave(filename = sprintf('Result/rare/core_bugs/All/Organ.All.(asv_collapse_to_sp).abun_%s.pdf', taxa[sp, 'Name']), 
                   plot = g, width = 2.5, height = 2.20)
            
            View(sp.cor[[sp]])
          }
        }
        
        
        stop('')
      }
      
      
      ## - 02. microbiome core
      if (TRUE) {
        Ind_summ = Ind_summ.UD # Ind_summ.LD, Ind_summ.UD, Ind_summ.Oral_gut
        type = 'UD'
        ##
        ID <- Reduce(union, lapply(Ind_summ, function(x) rownames(x)))
        summ <- matrix(data = 0, nrow = length(ID), ncol = length(Ind_summ), 
                       dimnames = list(ID, names(Ind_summ)))
        summ.oral <- summ
        tmp <- sapply(names(Ind_summ), function(n) {
          # print(n)
          ind = Ind_summ[[n]]
          summ[rownames(ind), n] <<- 1
          summ.oral[rownames(ind)[ind$oral == 1], n] <<- 1
          return(NA)
        })
        summ <- as.data.frame(summ)
        summ.oral <- as.data.frame(summ.oral)
        
        # table(rowSums(summ))
        # summ <- summ[rowSums(summ) > 1, ]
        
        ##
        taxa = taxonomy[rownames(summ), c(names(mapping), 'Species')]
        taxa[is.na(taxa)] <- ''
        taxa$Name = apply(taxa, 1, paste, collapse = ';')
        tmp.s = unique(taxa$Name)
        
        sp.summ <- matrix(data = 0, nrow = length(tmp.s), ncol = dim(summ)[2], 
                          dimnames = list(tmp.s, colnames(summ)))
        sp.summ.oral <- sp.summ
        tmp <- sapply(tmp.s, function(s) {
          idx = rownames(taxa)[taxa$Name == s]
          if (length(idx) == 1) {
            sp.summ[s, ] <<- as.numeric(summ[idx, ])
            sp.summ.oral[s, ] <<- as.numeric(summ.oral[idx, ])
          } else if (length(idx) > 1) {
            sp.summ[s, ] <<- colSums(summ[idx, ])
            sp.summ.oral[s, ] <<- colSums(summ.oral[idx, ])
          }
          return(NA)
        })
        sp.summ <- as.data.frame(sp.summ)
        sp.summ = sp.summ[rownames(sp.summ) != "Bacteria;;;;;;", ]
        sp.summ.oral <- as.data.frame(sp.summ.oral)
        sp.summ.oral = sp.summ.oral[rownames(sp.summ), ]
        
        save(summ, sp.summ, sp.summ.oral,
             file = sprintf('Result/rare/core_bugs/Ind/sp_summ_%s.RData', type))
        
        
        ## From Here
        if (FALSE) {
          ## LD
          # summ, sp.summ, sp.summ.oral
          load(file = sprintf('Result/rare/core_bugs/Ind/Ind_summ_LD.RData'))
          type = 'LD'
          if (FALSE) {
            taxa = taxonomy[rownames(summ), c(names(mapping), 'Species')]
            taxa[is.na(taxa)] <- ''
            taxa$Name = apply(taxa, 1, paste, collapse = ';')
            taxa = taxa[taxa$Name != 'Bacteria;;;;;;', ]
            
            # K.pneumoniae
            idx = rownames(taxa)[grepl('Klebsiella;pneumoniae', taxa$Name)]
            tmp = summ[idx, ]
            
            tmp = query_seq[idx]
            seqinr::write.fasta(sequences = as.list(tmp), names = as.list(names(tmp)), 
                                file.out = sprintf('Result/rare/core_bugs/Ind/K.pneumoniae_%s.fasta', type))
          }
          if (TRUE) {
            taxa = as.data.frame(do.call('rbind', lapply(rownames(sp.summ), function(x) {
              x = stringi::stri_split(x, regex = ';')[[1]]
              paste(c(mapping,"s"), x, sep = '__')
            })), stringsAsFactors = FALSE)
            colnames(taxa) = c(mapping, 's')
            rownames(taxa) = rownames(sp.summ)
            taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
              f <- get_lv(x = x, usedLv = 'f')
              g <- get_lv(x = x, usedLv = 'g')
              s <- get_lv(x = x, usedLv = 's')
              if (startsWith(s, 's__')) {
                s <- paste(g, s, sep = ';')
              }
              return(paste(f, s, sep = ';'))
            })
            
            map.v3v4_full <- list()
            map.v3v4_full[['f__Enterobacteriaceae']] = rownames(taxa)[grepl(pattern = 'f__Enterobacteriaceae', x = taxa$Name) & 
                                                                        !grepl(pattern = 'g__Klebsiella', x = taxa$Name)]
            map.v3v4_full[['g__Klebsiella']] = rownames(taxa)[grep(pattern = 'g__Klebsiella', x = taxa$Name)]
            map.v3v4_full[['g__[Ruminococcus];s__gnavus']] = rownames(taxa)[grep(pattern = 'g__\\[Ruminococcus\\].*gnavus*', x = taxa$Name)]
            map.v3v4_full[['g__Enterococcus']] = rownames(taxa)[grep(pattern = 'g__Enterococcus', x = taxa$Name)]
            # map.v3v4_full[['f__Lachnospiraceae']] = rownames(taxa)[grepl(pattern = 'f__Lachnospiraceae', x = taxa$Name) & 
            #                                                             !grepl(pattern = '(g__Blautia)|(g__\\[Ruminococcus\\].*gnavus*)|(g__\\[Ruminococcus\\].*torques*)', x = taxa$Name)]
            map.v3v4_full[['g__Bacteroides']] = rownames(taxa)[grepl(pattern = 'g__Bacteroides', x = taxa$Name) & 
                                                                 !grepl(pattern = 's__fragilis|s__caccae', x = taxa$Name)]
            # map.v3v4_full[['g__Blautia']] = rownames(taxa)[grep(pattern = 'g__Blautia', x = taxa$Name)]
            map.v3v4_full[['g__Bacteroides;s__fragilis']] = rownames(taxa)[grep(pattern = 'g__Bacteroides;s__fragilis', x = taxa$Name)]
            map.v3v4_full[['g__Bacteroides;s__caccae']] = rownames(taxa)[grep(pattern = 'g__Bacteroides;s__caccae', x = taxa$Name)]
            # map.v3v4_full[['g__Parabacteroides;s__distasonis']] = rownames(taxa)[grep(pattern = 'g__Parabacteroides;s__distasonis', x = taxa$Name)]
            # map.v3v4_full[['g__[Ruminococcus];s__torques']] = rownames(taxa)[grep(pattern = 'g__\\[Ruminococcus\\].*torques*', x = taxa$Name)]
            map.v3v4_full[['g__Faecalibacterium;s__prausnitzii']] = rownames(taxa)[grep(pattern = 'g__Faecalibacterium;s__prausnitzii', x = taxa$Name)]
          }
          
          ## UD
          # summ, sp.summ, sp.summ.oral
          load(file = sprintf('Result/rare/core_bugs/Ind/Ind_summ_UD.RData'))
          type = 'UD'
          if (TRUE) {
            taxa = as.data.frame(do.call('rbind', lapply(rownames(sp.summ), function(x) {
              x = stringi::stri_split(x, regex = ';')[[1]]
              paste(c(mapping,"s"), x, sep = '__')
            })), stringsAsFactors = FALSE)
            colnames(taxa) = c(mapping, 's')
            rownames(taxa) = rownames(sp.summ)
            taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
              f <- get_lv(x = x, usedLv = 'f')
              g <- get_lv(x = x, usedLv = 'g')
              s <- get_lv(x = x, usedLv = 's')
              if (startsWith(s, 's__')) {
                s <- paste(g, s, sep = ';')
              }
              return(paste(f, s, sep = ';'))
            })
            
            map.v3v4_full <- list()
            map.v3v4_full[['g__Klebsiella']] = rownames(taxa)[grep(pattern = 'g__Klebsiella', x = taxa$Name)]
            map.v3v4_full[['g__Streptococcus']] = rownames(taxa)[grepl(pattern = 'g__Streptococcus', x = taxa$Name) & 
                                                                   !grepl(pattern = 's__anginosus', x = taxa$Name)]
            map.v3v4_full[['g__Enterococcus']] = rownames(taxa)[grep(pattern = 'g__Enterococcus', x = taxa$Name)]
            map.v3v4_full[['g__Helicobacter;s__pylori']] = rownames(taxa)[grep(pattern = 'g__Helicobacter;s__pylori', x = taxa$Name)]
            map.v3v4_full[['f__Enterobacteriaceae']] = rownames(taxa)[grepl(pattern = 'f__Enterobacteriaceae', x = taxa$Name) & 
                                                                        !grepl(pattern = 'g__Klebsiella', x = taxa$Name)]
            map.v3v4_full[['g__Peptostreptococcus']] = rownames(taxa)[grep(pattern = 'g__Peptostreptococcus', x = taxa$Name)]
            # map.v3v4_full[['g__Lactobacillus']] = rownames(taxa)[grep(pattern = 'g__Lactobacillus', x = taxa$Name)]
            map.v3v4_full[['g__Streptococcus;s__anginosus']] = rownames(taxa)[grep(pattern = 'g__Streptococcus;s__anginosus', x = taxa$Name)]
            map.v3v4_full[['g__Haemophilus;s__parainfluenzae']] = rownames(taxa)[grep(pattern = 'g__Haemophilus;s__parainfluenzae', x = taxa$Name)]
            # map.v3v4_full[['g__Clostridium;s__perfringens']] = rownames(taxa)[grep(pattern = 'g__Clostridium.*perfringens.*', x = taxa$Name)]
            # map.v3v4_full[['g__Stenotrophomonas']] = rownames(taxa)[grep(pattern = 'g__Stenotrophomonas', x = taxa$Name)]
            map.v3v4_full[['g__Prevotella;s__melaninogenica']] = rownames(taxa)[grep(pattern = 'g__Prevotella.*;s__melaninogenica', x = taxa$Name)]
            # map.v3v4_full[['g__Corynebacterium']] = rownames(taxa)[grep(pattern = 'g__Corynebacterium', x = taxa$Name)]
          }
          
          ## Oral_gut
          # summ, sp.summ, sp.summ.oral
          load(file = sprintf('Result/rare/core_bugs/Ind/Ind_summ_Oral_gut.RData'))
          type = 'Oral_gut'
          if (TRUE) {
            taxa = as.data.frame(do.call('rbind', lapply(rownames(sp.summ), function(x) {
              x = stringi::stri_split(x, regex = ';')[[1]]
              paste(c(mapping,"s"), x, sep = '__')
            })), stringsAsFactors = FALSE)
            colnames(taxa) = c(mapping, 's')
            rownames(taxa) = rownames(sp.summ)
            taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
              f <- get_lv(x = x, usedLv = 'f')
              g <- get_lv(x = x, usedLv = 'g')
              s <- get_lv(x = x, usedLv = 's')
              if (startsWith(s, 's__')) {
                s <- paste(g, s, sep = ';')
              }
              return(paste(f, s, sep = ';'))
            })
            
            map.v3v4_full <- list()
            map.v3v4_full[['g__Klebsiella']] = rownames(taxa)[grep(pattern = 'g__Klebsiella', x = taxa$Name)]
            map.v3v4_full[['g__Enterococcus']] = rownames(taxa)[grep(pattern = 'g__Enterococcus', x = taxa$Name)]
            map.v3v4_full[['g__Streptococcus']] = rownames(taxa)[grep(pattern = 'g__Streptococcus', x = taxa$Name)]
            map.v3v4_full[['g__Dialister']] = rownames(taxa)[grep(pattern = 'g__Dialister', x = taxa$Name)]
          }
          
          ##
          sp.summ = sp.summ[order(rowSums(sp.summ != 0), decreasing = TRUE), ]
          sp.summ = sp.summ[Reduce(union, map.v3v4_full), ]
        }
        # # sp.summ = sp.summ[rowSums(sp.summ != 0) > 1, ]
        # sp.summ = sp.summ[order(rowSums(sp.summ != 0), decreasing = TRUE), ]
        sp.summ = sp.summ[, sort(colnames(sp.summ))]
        
        sp.summ.oral <- sp.summ.oral[rownames(sp.summ), colnames(sp.summ)]
        
        taxa = as.data.frame(do.call('rbind', lapply(rownames(sp.summ), function(x) {
          x = stringi::stri_split(x, regex = ';')[[1]]
          paste(c(mapping,"s"), x, sep = '__')
        })), stringsAsFactors = FALSE)
        colnames(taxa) = c(mapping, 's')
        rownames(taxa) = rownames(sp.summ)
        taxa$Name <- apply(X = taxa[, c(mapping, 's')], MARGIN = 1, FUN = function(x) {
          # f <- get_lv(x = x, usedLv = 'f')
          g <- get_lv(x = x, usedLv = 'g')
          s <- get_lv(x = x, usedLv = 's')
          if (startsWith(s, 's__')) {
            s <- paste(g, s, sep = ';')
          }
          # return(paste(f, s, sep = ';'))
          return(s)
        })
        
        metadat_core.ind <- metadat_core.ind[colnames(sp.summ), ]
        
        Fun_row_anno <- rowAnnotation(
          Text = anno_text(paste(' ', taxa[rownames(sp.summ), 'Name'], " ", sep = ''), just = 'right', location = unit(1, 'npc'),
                           #rownames(taxonomy[rownames(sp.summ), ]),
                           gp = gpar(cex = 0.5, col = 'black', fontface = 'bold')),
          oral = rowSums(sp.summ.oral != 0),
          # status = taxa$status,
          col = list(oral = colorRamp2(c(0, 33), c("white", "black")),
                     status = c('up' = '#ff0000', 'down' = '#028900')),
          border = TRUE,
          gap = unit(2, "points"),
          simple_anno_size = unit(0.25, "cm"),
          gp = gpar(col = NA, width = 0.1),
          show_annotation_name = TRUE,
          annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                         nrow = 9,
                                         border = 'black',
                                         grid_width = unit(5, "mm"))
        )
        
        Fun_row_right <- rowAnnotation(
          Ratio = anno_barplot(rowSums(sp.summ != 0) / dim(sp.summ)[2], border = FALSE),
          border = FALSE,
          gap = unit(2, "points"),
          simple_anno_size = unit(0.25, "cm"),
          gp = gpar(col = NA, width = 0.1),
          show_annotation_name = TRUE,
          annotation_legend_param = list(legend_direction = "hor", #legend_direction = "vertical",
                                         nrow = 9,
                                         border = 'black',
                                         grid_width = unit(5, "mm"))
        )
        
        Fun_col_anno <- HeatmapAnnotation(
          which = 'column',
          # ICU = as.character(metadat_core.ind$`ICU (Yes: 1; No: 0)`),
          # time = as.numeric(metadat_core.ind$`Sampling time after death  (Min)`),
          # Antibiotic = as.character(metadat_core.ind$`Antibiotic uses (Yes: 1; No: 0)`),
          # Hospitalization = as.numeric(metadat_core.ind$Hospitalization),
          BMI = as.numeric(metadat_core.ind$BMI),
          Age = as.numeric(metadat_core.ind$Age),
          # Gender = metadat_core.ind$Gender,
          col = list(ICU = c('0' = 'white', '1' = 'black'),
                     time = colorRamp2(c(70, 100), c("white", "black")),
                     Antibiotic = c('0' = 'white', '1' = 'black'),
                     Hospitalization = colorRamp2(c(1, 10), c("white", "black")),
                     BMI = colorRamp2(c(0, 30), c("white", "black")),
                     Age = colorRamp2(c(0, 100), c("white", "black")),
                     Gender = c(Male = '#3399ff', Female = '#ff93ac')),
          border = TRUE,
          show_annotation_name = TRUE,
          annotation_name_gp = gpar(cex = 0.6),
          gap = unit(2, "points"),
          simple_anno_size = unit(0.2, "cm"), 
          # gp = gpar(col = "grey", width = 0.1),
          annotation_legend_param = list(legend_direction = "vertical",
                                         nrow = 9,
                                         border = 'black',
                                         grid_width = unit(5, "mm"))
        )
        
        
        hp <- Heatmap(matrix = as.matrix(sp.summ), 
                      col = colorRamp2(breaks = c(0, 2), colors = c("white", "red"), space = 'sRGB'),
                      # col = c('0' = 'white', '1' = 'red'),
                      
                      border = 'black', 
                      # rect_gp = gpar(color = '#f0f8ff', lwd = 0.01),
                      
                      cluster_rows = FALSE,
                      # row_order = ordering_row,
                      cluster_columns = FALSE,
                      # column_order = ordering_col,
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                      
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      column_names_gp = gpar(cex = 0.6),
                      
                      split = factor(do.call(c, lapply(names(map.v3v4_full), function(x) {
                        rep(x, length(map.v3v4_full[[x]]))
                      })), names(map.v3v4_full)),
                      row_title = NULL,
                      row_gap = unit(0, "mm"),
                      
                      right_annotation = Fun_row_right,
                      left_annotation = Fun_row_anno,
                      # top_annotation = Fun_col_anno,
                      # bottom_annotation = Fun_col_anno,
                      
                      show_heatmap_legend = TRUE, 
                      heatmap_legend_param = list(legend_direction = "ver", title = ''),
                      column_title = paste('', '', sep = ''))
        
        draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
             heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        
        
        pdf(sprintf('Result/rare/core_bugs/Ind/ind(with ID)_%s-.pdf', type), width = 5.04, height = 3.52) # LD
        pdf(sprintf('Result/rare/core_bugs/Ind/ind(with ID)_%s-.pdf', type), width = 5.04, height = 4.15) # UD
        pdf(sprintf('Result/rare/core_bugs/Ind/ind(with ID)_%s-.pdf', type), width = 5.04, height = 3.14) # Oral_gut
        
        draw(hp, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
             heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
      }
      
    }
    
    stop('')
  }
  
}


