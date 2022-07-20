#####Functions to be called to dada2_workflow.Rmd

#==========AllOrients==========
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

#==========PrimerHits==========
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  require(ShortRead)
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


###### Functions to be called to phyloseq ACTUAL.Rmd

#==========Fill NAs in taxa==========
seqtab_noNAs <- function(seqtab, taxa_levels = c("silva-genus", "silva-species", "pr2")) {
  if (taxa_levels == "silva-genus") {
    tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  } else if (taxa_levels == "silva-species") {
    tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  } else if (taxa_levels == "pr2" ) {
    tax_levels <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
  } else {
    tax_levels <- taxa_levels
  }
  
  for (i in 1:(length(tax_levels)-1)){
    ps_na <- which(is.na(seqtab[,(i+1)]))
    seqtab[ps_na,(i+1)] <- str_c(seqtab[ps_na,i],"_X")
  }
  return(seqtab)
}





#==========Normalise===============
ps_normalize_median <- function(ps, title) {
  ps_median = median(sample_sums(ps))
  cat(sprintf("\nThe median number of reads used for normalization of %s is  %.0f", 
              title, ps_median))
  normalize_median = function(x, t = ps_median) (if (sum(x) > 0) {
    t * (x/sum(x))
  } else {
    x
  })
  ps = transform_sample_counts(ps, normalize_median)
  cat(str_c("\nPhyloseq ", title, "\n========== \n"))
  print(ps)
}

#==============find abundant taxa============
ps_abundant <- function(ps, contrib_min = 0.1, title) {
  total_per_sample <- max(sample_sums(ps))
  ps <- filter_taxa(ps, function(x) sum(x > total_per_sample * contrib_min) > 0, TRUE)
  ps <- ps_normalize_median(ps, title)
}

#=========convert phyloseq object to long format==========
ps_to_long <- function(ps) {
  otu_df <- data.frame(otu_table(ps),stringsAsFactors =FALSE) %>%
    rownames_to_column(var = "otu_id")
  taxo_df <- data.frame(tax_table(ps),stringsAsFactors =FALSE) %>%
    rownames_to_column(var = "otu_id")
  otu_df <- left_join(taxo_df, otu_df)
  otu_df <- gather(otu_df, "sample", "n_seq", contains(c("WDL", "SBW", "STL", "SNB")))
  metadata_df <- data.frame(sample_data(ps), stringsAsFactors =FALSE) %>% 
    rownames_to_column(var = "sample")
  otu_df <- left_join(otu_df, metadata_df)
}

#==========make treemap==========
treemap_gg_dv <- function(df, group1, group2, title) {
  
  df <- df %>% 
    group_by({{group1}}, {{group2}}) %>% 
    summarise(n_seq = sum(n_seq)) %>% 
    na.omit()
  
  g_treemap <- ggplot(df, aes(area = n_seq, fill = {{group2}}, label = {{group2}}, subgroup = {{group1}})) + 
    ggtitle(title) + 
    treemapify::geom_treemap() + 
    treemapify::geom_treemap_subgroup_border() + 
    treemapify::geom_treemap_text(colour = "black", place = "topleft", reflow = T,
                                  padding.x = grid::unit(3, "mm"),
                                  padding.y = grid::unit(3, "mm")) + 
    treemapify::geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5,
                                           colour = "white", fontface = "italic", min.size = 0) + 
    scale_fill_viridis(discrete = TRUE) + 
    theme(legend.position = "none", plot.title = element_text(size = 16, face = "bold"))
  
  return(g_treemap)
}

#==========Theme for DESeq2 plot==========
theme_dviz_grid <- function(font_size = 14, font_family = "") {
  color = "grey90"
  line_size = 0.5
  
  # Starts with theme_cowplot and then modify some parts
  cowplot::theme_cowplot(font_size = font_size, font_family = font_family) %+replace%
    theme(
      # make horizontal grid lines
      panel.grid.major   = element_line(colour = color,
                                        size = line_size),
      
      # adjust axis tickmarks
      axis.ticks        = element_line(colour = color, size = line_size),
      
      # no x or y axis lines
      axis.line.x       = element_blank(),
      axis.line.y       = element_blank()
    )
}

#==========Phyloseq Object to DESeq2 Plot==========
ps_do_deseq <- function(deseq_list) {
  require(DESeq2)
  plot_array <- list()
  
  for (i in 1:length(deseq_list$deseq)) {
    deseq_data <- deseq_list$deseq[[i]]
    
    diagdds <- phyloseq_to_deseq2(deseq_data, ~ Monsoon)
    
    #============Estimate geometric mean with zeros============
    gm_mean <- function(x, na.rm = TRUE) {
      exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
    }
    
    geoMeans <- apply(counts(diagdds), 1, gm_mean)
    diagdds <- DESeq2::estimateSizeFactors(diagdds, geoMeans = geoMeans)
    #===========================end============================
    diagdds = DESeq(diagdds, fitType = "local")
    
    res = results(diagdds, cooksCutoff = FALSE)
    alpha = 0.005
    sigtab = res[which(res$padj < alpha), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(deseq_list$deseq[[i]])[rownames(sigtab), ], "matrix"))
    sigtablen <- nrow(sigtab)
    sigtabgen <- length(unique(sigtab$Genus))
    cat(sigtablen, "ASVs in", sigtabgen, "unique genus that are significantly different. \n")
    
    sigtab_plot <- sigtab %>% 
      rownames_to_column(var = "otu_id") %>% 
      arrange(desc(abs(log2FoldChange))) %>% 
      top_n(-100, wt = otu_id) %>% 
      mutate(label = case_when(Kingdom == "Eukaryota" ~ str_c(otu_id, Class, Genus, sep = "_"),
                               TRUE ~ str_c(otu_id, Class, Genus, sep = "_"))) %>%
      rowid_to_column()
    
    g <- ggplot(sigtab_plot,
                aes(x = reorder(label, log2FoldChange),
                    y = log2FoldChange,
                    color = Kingdom, alpha = log(padj, base = 10))) + 
      theme_dviz_grid() + 
      geom_point(size = 3) + 
      theme(axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.title.x = element_text(size = 7),
            axis.text.y = element_text(size = 7, angle = 0, hjust = 0, vjust = 0.5),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 15, hjust = 0.5)) + 
      ggtitle("Significantly Different ASVs") +
      guides(alpha = FALSE) +
      ylab(str_c(deseq_list$title[i])) + 
      xlab("") + 
      scale_color_viridis_d(option = "C", name = "Kingdom") + 
      coord_flip()
    
    plot_array[[i]] <- g
  }
  return(plot_array)
}

#==========Phyloseq Object to DESeq2 Plot==========
ps_do_deseq_genus <- function(deseq_list) {
  require(DESeq2)
  plot_array <- list()
  
  for (i in 1:length(deseq_list$deseq)) {
    deseq_data <- deseq_list$deseq[[i]]
    
    diagdds <- phyloseq_to_deseq2(deseq_data, ~ type)
    
    #============Estimate geometric mean with zeros============
    gm_mean <- function(x, na.rm = TRUE) {
      exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
    }
    
    geoMeans <- apply(counts(diagdds), 1, gm_mean)
    diagdds <- DESeq2::estimateSizeFactors(diagdds, geoMeans = geoMeans)
    #===========================end============================
    diagdds = DESeq(diagdds, fitType = "local")
    
    res = results(diagdds, cooksCutoff = FALSE)
    alpha = 0.01
    sigtab = res[which(res$padj < alpha), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(deseq_list$deseq[[i]])[rownames(sigtab), ], "matrix"))
    
    sigtab_plot <- sigtab %>% 
      rownames_to_column(var = "otu_id") %>%
      mutate(label = case_when(Kingdom == "Eukaryota" ~ str_c(Class, Genus, sep = "_"),
                               TRUE ~ str_c(Class, Genus, sep = "_"))) %>%
      rowid_to_column()
    
    sigtab_plot_genus <- sigtab_plot %>% 
      group_by(Genus) %>% 
      summarise(across(c(4, 8), median))
    
    sigtab_plot_collapsed <- sigtab_plot[match(unique(sigtab_plot$Genus), sigtab_plot$Genus) , c(9:17)] %>% 
      left_join(sigtab_plot_genus, by = "Genus") %>% 
      rowid_to_column("GenusNumber") %>% 
      mutate(GenusNumber = sprintf("genus%03d", GenusNumber))
    sigtabgen <- nrow(sigtab_plot_collapsed)
    cat(sigtabgen, "unique genus that are significantly different. \n")
    
    g <- ggplot(sigtab_plot_collapsed,
                aes(x = reorder(label, log2FoldChange),
                    y = log2FoldChange,
                    color = Kingdom, alpha = log(padj, base = 10))) + 
      theme_dviz_grid() + 
      geom_point(size = 3) + 
      theme(axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.title.x = element_text(size = 7),
            axis.text.y = element_text(size = 7, angle = 0, hjust = 0, vjust = 0.5),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 15, hjust = 0.5)) + 
      ggtitle("Significantly Different Genus") +
      guides(alpha = FALSE) +
      ylab(str_c(deseq_list$title[i])) + 
      xlab("") + 
      scale_color_viridis_d(option = "C", name = "Kingdom") + 
      coord_flip()
    
    plot_array[[i]] <- g
  }
  return(plot_array)
}

################################################################################
# Define an internal function for accessing and orienting the OTU table
# in a fashion suitable for vegan functions
# @keywords internal
# https://rdrr.io/bioc/phyloseq/src/R/ordination-methods.R
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#==========Phyloseq to Dataframe==========

phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  
  # require(phyloseq)
  
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }
  
  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }
  
  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }
    
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }
  
  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
      
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
      
    } # end of addmaxrank
  }   # end of addtax
  
  ## Reorder OTUs
  if(!is.null(sorting)){
    
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }
    
    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }
  
  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }
  
  rownames(res) <- NULL
  return(res)
}

#============convrt phyloseq object into sth vegan can digest================
veganotu <- function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}


#============Make gephi object================
ps_to_gephi <- function(ps, max_dist = 0.7, file_name = "./ps_to_gephi.gexf"){
  require("rgexf")
  require("phyloseq")
  require("dplyr")
  require("igraph")
  ############################################################################################
  # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
  taxa_fam <- as(tax_table(ps), "matrix")
  network_fam <- simplify(make_network(ps, type = "taxa", distance = "bray", max.dist = max_dist))
  
  # Print number of nodes and edges
  cat("\n Number of vertices:")
  print(vcount(network_fam))
  cat("\n Number of edges:")
  print(ecount(network_fam))
  
  #Get names
  names_fam <- rownames(as.matrix(V(network_fam)))
  #names_fam
  
  #Edit names
  taxa_fam <- data.frame(taxa_fam, stringsAsFactors = F) %>% 
    rownames_to_column()
  names_fam_taxa <- data.frame(rowname = names_fam) %>% 
    left_join(taxa_fam)
  
  network_fam <- set.vertex.attribute(network_fam, "name", value = names_fam_taxa$Family)
  V(network_fam)
  
  ##############################################################################################
  # Create a dataframe nodes: 1st column - node ID, 2nd column -node name
  nodes_df <- data.frame(ID = c(1:vcount(network_fam)), NAME = V(network_fam)$name)
  # Create a dataframe edges: 1st column - source node ID, 2nd column -target node ID
  edges_df <- as.data.frame(get.edges(network_fam, c(1:ecount(network_fam))))
  
  write.gexf(nodes = nodes_df, edges = edges_df, defaultedgetype = "undirected", output = file_name)
}


#============Ps to sparCC input==============

ps_to_SparCCinput <- function(ps, normalise = "median", min_sample_reads = 1e4, min_mean_reads = 2, min_occurrence = 0.1, level = c("Family", "ASV")){
  # data preparation for SparCC input
  # Removing samples with low read counts
  # Removing taxa with low reads on average
  # Removing taxa with minimum 10% occurrence across samples
  # Normalise samples to median (default) or centered log-ratio (in-progress)
  require(tidyverse)
  require(phyloseq)
  require(SpiecEasi)
  
  asvtable <- otu_table(ps) # get table of abundances from phyloseq object
  asvtable <- data.frame(asvtable@.Data)
  
  temp_tokeep <- colSums(asvtable) > min_sample_reads
  asvtable <- asvtable[ , temp_tokeep] # remove samples with low read counts
  
  temp_tokeep <- apply(asvtable, 1, mean) > min_mean_reads
  asvtable <- asvtable[temp_tokeep, ] # remove families with <2 reads on average
  
  temp_tokeep <- rowSums(asvtable != 0) > (min_occurrence * nrow(asvtable))
  asvtable <- asvtable[temp_tokeep, ] # remove families with occurences in all samples <10%
  
  if(normalise == "median"){
    # normalise to median read count!
    asvtable_norm <- asvtable
    m <- median(colSums(asvtable_norm))
    for(c in 1:ncol(asvtable_norm)) {asvtable_norm[ , c] <- asvtable[ , c]/colSums(asvtable)[c]*m}
  } # else if(normalise == "clr"){
  #   # normalise to centered log-ratio
  #   asvtable_norm <- asvtable %>% 
  #     clr()
  # }
  
  
  # get corresponding family name & asvno
  TAXtable <- tax_table(ps) 
  TAXtable <- data.frame(TAXtable@.Data) %>% 
    rownames_to_column(var = "ASVno") %>% 
    select(ASVno, Family)

  # make table for SparCC input  
  if(level == "Family") {
    asvtable_input <- asvtable_norm %>% 
      rownames_to_column(var = "ASVno") %>% 
      left_join(TAXtable, by = "ASVno") %>% 
      arrange(Family) %>% 
      column_to_rownames("Family") %>% 
      select(where(is.numeric)) %>% 
      t()
  }else if(level == "ASV"){
      asvtable_input <- asvtable_norm %>% 
        rownames_to_column(var = "ASVno") %>% 
        left_join(TAXtable, by = "ASVno") %>% 
        arrange(Family) %>% 
        column_to_rownames("ASVno") %>% 
        select(where(is.numeric)) %>% 
        t()
  }

  return(asvtable_input)
}

#============Running SparCC & tidying data==============

SparCC_inout <- function(FAMtable_input, seed = 42, boot_perm = 1000, threads = 2){
  # ======================================================
  #  SparCC correlation coefficients and p-values calculation
  #  Adapted from https://biovcnet.github.io/_pages/NetworkScience_SparCC.nb.html
  #  and https://rachaellappan.github.io/16S-analysis/correlation-between-otus-with-sparcc.html
  # ======================================================
  require(tidyverse)
  require(SpiecEasi)
  
  set.seed(seed)
  out_sparcc <- sparcc(FAMtable_input) # calculate sparcc correlation values
  
  rownames(out_sparcc$Cor) <- colnames(FAMtable_input) # re-add column/row names with the family
  colnames(out_sparcc$Cor) <- colnames(FAMtable_input)
  rownames(out_sparcc$Cov) <- colnames(FAMtable_input)
  colnames(out_sparcc$Cov) <- colnames(FAMtable_input)
  
  out_sparboot <- sparccboot(FAMtable_input, R = boot_perm, ncpus = threads) # calculate bootstrapped/pseudo p-values. 1000 permutations are recommended. can do (R = 1000, ncpus = 8) in the server later
  saveRDS(out_sparboot, "21_out_sparboot.RDS")
  
  out_sparempval <- pval.sparccboot(out_sparboot)  # calculate empirical p-values from bootstrapped p-values.
  # Outputs the lower part of a diagonal matrix.
  
  # Rescue using this:
  cors <- out_sparempval$cors
  pvals <- out_sparempval$pvals
  #correlation
  out_cors <- diag(0.5, nrow = dim(out_sparcc$Cor)[1], ncol = dim(out_sparcc$Cor)[1])
  out_cors[upper.tri(out_cors, diag=FALSE)] <- cors # fill upper triangle of out_cors with cors
  out_cors <- out_cors + t(out_cors) # fill lower triangle with the same thing
  #p-value
  out_pval <- diag(0.5, nrow = dim(out_sparcc$Cor)[1], ncol = dim(out_sparcc$Cor)[1])
  out_pval[upper.tri(out_pval, diag=FALSE)] <- pvals #fill upper triangle of out_pval with pvals
  out_pval <- out_pval + t(out_pval) #fill lower triangle with the same thing
  
  rownames(out_cors) <- colnames(FAMtable_input)
  colnames(out_cors) <- colnames(FAMtable_input)
  rownames(out_pval) <- colnames(FAMtable_input)
  colnames(out_pval) <- colnames(FAMtable_input)
  
  return(list(cors = out_cors, pval = out_pval))
}

#============SparCC output to igraph object==============

SparCC_to_igraph <- function(list_cor_pval, min_corr = 0.3, max_pval = 0.05){
  # Remove the links that correspond to absolute correlation <0.3  and for which p-value is > 0.05
  # Daniel removed every corr <0.18 including negative values, and uses pseudo p-value
  require(tidyverse)
  require(igraph)
  
  graph_cor <- list_cor_pval$cors
  graph_cor [abs(graph_cor) < min_corr] <- 0
  graph_cor [list_cor_pval$pval > max_pval] <- 0 
  
  # create the graph
  my_graph_sparcc <- graph.adjacency(graph_cor, weighted=TRUE, mode="upper")
  
  # Remove the links that correspond to self nodes
  my_graph <- simplify(my_graph_sparcc)
  
  # --- add edge attributes ---
  my_graph$layout <- layout_with_fr #The Fruchterman-Reingold layout algorithm
  E(my_graph)$width <- 50^abs(E(my_graph)$weight) # Make the width of the edges proportional to the absolute correlation. )weight = corr coef!!
  E(my_graph)$color <- ifelse(E(my_graph)$weight > 0,'blue', 'red') #Make the colour of the edges red if -ve correlation, blue if +ve correlation
  
  return(my_graph)
}

#============Running SpiecEasi MB & tidying data==============

SpiecEasi_inout <- function(tab_graph_input, method = "mb", repetitions = 100, threads = 8){
  set.seed(42)
  spieceasi.net <- spiec.easi(tab_graph_input, method = method, lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = repetitions, ncores = threads)) # reps have to increase for real data
  
  spieceasi.matrix.dsc<- symBeta(getOptBeta(spieceasi.net), mode = "maxabs") # get optimal coefficient matrix, and make it symmetrical by taking the max value between the two values
  spieceasi.matrix <- as.matrix(spieceasi.matrix.dsc)
  
  spieceasi.stability.dsc<- getOptMerge(spieceasi.net) # get a symmetric matrix of edge-wise stability. Works like [1-p.value], bigger number = greater stability.
  spieceasi.stability <- as.matrix(spieceasi.stability.dsc)
  
  colnames(spieceasi.matrix) <- rownames(spieceasi.matrix) <- colnames(spieceasi.stability) <- rownames(spieceasi.stability) <- colnames(tab_graph_input) # rename columns & row with ASV number
  
  return(list(cors = spieceasi.matrix, stab = spieceasi.stability))
}

#============SpiecEasi output to igraph object==============

SpiecEasi_to_igraph <- function(SpiecEasi_out_list, min_coef = 0, min_stab = 0){
  spieceasi.matrix <- SpiecEasi_out_list$cors
  spieceasi.stabil <- SpiecEasi_out_list$stab
  
  spieceasi.matrix[abs(spieceasi.matrix) < min_coef] <- 0 # remove anything with coef lower than min_coef
  spieceasi.matrix[spieceasi.stabil < min_stab] <- 0      # remove anything with stab lower than min_stab
  
  # create the graph
  otu_graph_mb <- graph.adjacency(spieceasi.matrix, mode = "upper", weighted = TRUE, diag = FALSE)
  
  # Remove the links that correspond to self nodes
  otu_graph_mb <- simplify(otu_graph_mb)
  
  # --- add edge attributes ---
  otu_graph_mb$layout <- layout_with_fr #The Fruchterman-Reingold layout algorithm
  E(otu_graph_mb)$width <- 50^abs(E(otu_graph_mb)$weight) # Make the width of the edges proportional to the absolute correlation. )weight = corr coef!!
  E(otu_graph_mb)$color <- ifelse(E(otu_graph_mb)$weight > 0,'blue', 'red') #Make the colour of the edges red if -ve correlation, blue if +ve correlation
  E(otu_graph_mb)$lty <- ifelse(E(otu_graph_mb)$weight > 0,'solid', 'dashed') #Make the linetype of the edges dashed if -ve correlation, solid if +ve correlation
  
  return(otu_graph_mb)
}


#============Saving igraph object for gephi input==============

igraph_to_gephi_jts <- function(my_graph, output_file = "my_graph.gexf"){
  require(igraph)
  require(rgexf)
  require(tidyverse)
  # Create a dataframe for nodes: 1st column node ID, 2nd column node name
  nodes_df <- data.frame(ID = c(1:vcount(my_graph)),
                         NAME = V(my_graph)$name)
  # Create a dataframe for edges: 1st column source node ID, 2nd column target node ID
  edges_df <- as.data.frame(get.edges(my_graph, c(1:ecount(my_graph))))
  # Create a dataframe for edge attribute: 1st column width of line, 2nd column color of line
  edge_attr_df <- data.frame(WEIGHT = E(my_graph)$weight,
                             ABSWEIGHT = abs(E(my_graph)$weight),
                             WIDTH = E(my_graph)$width,
                             COLOR = E(my_graph)$color)
  
  # remove nodes with no links
  temp_remove <- setdiff(nodes_df$ID, as.vector(c(edges_df$V1, edges_df$V2)))
  if(length(temp_remove) != 0){nodes_df <- nodes_df[ -temp_remove, ]}
  
  #write gexf for gephi input
  rgexf::write.gexf(nodes = nodes_df, edges = edges_df, edgesAtt=edge_attr_df, defaultedgetype = "undirected", output = output_file)
  cat("Import to gephi now :)")
}