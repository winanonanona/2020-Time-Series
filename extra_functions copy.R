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
