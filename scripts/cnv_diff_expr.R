suppressMessages(library('DESeq2'))
suppressMessages(library('edgeR'))
suppressMessages(library('foreach'))
suppressMessages(library('doMC'))



# ---------------- FUNCTION definitions --------------- #
perform_limma  <- function(data_matrix, samples, controls, gene){
    
    print(paste('Performing LIMMA-VOOM for gene', gene))
    
    row.names(data_matrix) <- data_matrix$Hugo_Symbol
    data_matrix$Hugo_Symbol <- NULL
  
     # get lengths of samples and counts for design matrix
    nsample   <- length(samples)
    ncontrols <- length(controls)
    totals    <- nsample + ncontrols
  
    # Set samples and controls as characters for indexing
    samples  <- as.character(samples)
    controls <- as.character(controls)

    design <- model.matrix(~0 + factor(c(rep(0, ncontrols), rep(1, nsample))))
    rownames(design) <- c(controls, samples)
    colnames(design) <- c('neutral', 'amplified')
    
    # Perform LIMMA-VOOM
    dge  <- DGEList(counts = data_matrix[,c(controls,samples)])
    keep <- filterByExpr(dge, design = design)
    dge  <- dge[keep,,keep.lib.sizes=FALSE]
    dge  <- calcNormFactors(dge)
    v    <- voom(dge, design = design)
    fit  <- lmFit(v, design = design)
    fit  <- eBayes(fit)
    results <- topTable(fit, coef=ncol(design), number = Inf)
   
  # check if our gene of interest has passed the test and it's logFC is positive
    p_value <- results[gene, 'adj.P.Val']
    log2FC  <- results[gene, 'logFC']
    
    return (c(p_value, log2FC))
}

perform_deseq2 <- function(data_matrix, samples, controls, gene){
  
    print(paste('Performing DESEq2 for gene', gene))
    
    # Set probes, samples and controls 
    tested_gene <- gene 

    # Set probes as row names
    row.names(data_matrix) <- data_matrix$Hugo_Symbol
    data_matrix$Hugo_Symbol <- NULL
    samples_controls <- c(as.character(samples), as.character(controls))
 
    # Select samples and controls of interest
    data_matrix <- data_matrix[,samples_controls]
    data_matrix <- as.matrix(data_matrix)

    # Create design matrix
    sample_names  <- rep('cnv_gain', times=length(samples))
    control_names <- rep('neutral', times=length(controls))
    designmatrix  <- data.frame('condition' = c(sample_names,control_names),
                               row.names=samples_controls)

    # RUN DESEQ2
    dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                                  colData = designmatrix,
                                  design = ~ condition)
   
    dds$condition <- relevel(dds$condition, ref= 'neutral')
    dds2 <- DESeq(dds, quiet = TRUE)
    res <- results(dds2)
    
    # check if our gene of interest has passed the test and it's logFC is positive
    p_value <- res[tested_gene, 'padj']
    log2FC  <- res[tested_gene, 'log2FoldChange']

    return (c(p_value, log2FC))
}

save_dummy_dataframe <-function(target)
{
    # Generating empty dataframe for now, so it won't crash downstream
    empty_data <- data.frame('Hugo_Symbol'= character(), 
                             'adj_pval'   = character(),
                             'log2fc'     = character(),
                             'method'     = character(),
                             'cases'      = character(),
                             'bias'       = character())
    
    write.csv(x=empty_data, file=where_to_save, row.names = FALSE)
}


# --------- SNAKEMAKE I/O -------------- # 
cases_table     = read.csv(snakemake@input[["cases_table"]], sep=',')
expr_matrix     = read.csv(snakemake@input[["expression_table"]], sep=',', check.names = FALSE)
where_to_save   = snakemake@output[['tested_cnv']]

# Setup available cores for parallel
registerDoMC(cores=snakemake@threads)
print(paste('Registered', snakemake@threads, 'threads'))


# Select cases with cnv_gains, GoF, high gscore and RNA-seq samples available
has_expr_data <- cases_table$case_id %in% colnames(expr_matrix)
cases_of_interest <- cases_table[has_expr_data,]

cnv_gain      <- cases_of_interest$copy_number == 'cnv_gain'
gof           <- cases_of_interest$Consequence == 'GoF'     
not_other     <- cases_of_interest$Variant_Classification == 'None' 
high_gscore   <- cases_of_interest$gscore >= 0.4 

genes_of_interest <- cases_of_interest[high_gscore &
                                       cnv_gain &
                                       gof & 
                                       not_other,
                                       c('Hugo_Symbol', 'case_id', 'Consequence')]

# Check if the dataframe is empty after filtering genes of interest
if(dim(genes_of_interest)[1] > 0)
{
    # get a table of genes with unique samples and sample count
    genes_to_test <- aggregate(case_id ~ Hugo_Symbol,
                               genes_of_interest,
                               function(x) length(unique(x)))
    
    colnames(genes_to_test) <- c('gene','count') 
    genes_to_test <- genes_to_test[order(-genes_to_test$count),] 
    
    # Select genes with sample count > 5 
    genes_to_test <- genes_to_test[genes_to_test$count >= 5, 'gene']
    
    # Check whether we have enough samples for DE
    if (length(genes_to_test) > 0){  
        
        print(paste(length(genes_to_test), 'genes will be tested for cnv DE'))
        
        # Loop over genes of interest to perform deseqs
        gtests <- foreach(gene=genes_to_test,
                          .packages = c('DESeq2', 'edgeR'),
                          .inorder = FALSE,
                          .combine = 'rbind') %dopar% {
                              
                              # cases are those patients sharing a cnv_gain event
                              cases <- unique(genes_of_interest[genes_of_interest$Hugo_Symbol == gene,'case_id'])
                              
                              # get cases not carrying our gene of interest as a GoF or LoF
                              not_cases        <- cases_of_interest[!cases_of_interest$case_id %in% cases,]
                              valid_gene_alt   <- (not_cases$Hugo_Symbol == gene & not_cases$Consequence == 'Unknown')
                              gene_not_alt     <- not_cases$Hugo_Symbol != gene
                              controls         <- unique(not_cases[valid_gene_alt | gene_not_alt, 'case_id'])
                              bias             <- 0 # Whether cases were coerced to controls or not
                              
                              if(length(cases) > length(controls)){
                                  cases <- sample(cases, size=length(controls))
                                  bias  <- 1
                              }else{
                                  controls = sample(controls, size = length(cases))
                              }
                              
                              total_samples <- length(cases) + length(controls)
                              
                              analysis_class <- 'unknown'
                              test_results   <- c(0,0)
                              
                              if(total_samples <= 10){
                                  analysis_class <- 'deseq2'  
                                  test_results <- perform_deseq2(data_matrix = expr_matrix,
                                                                 samples = cases,
                                                                 controls = controls,
                                                                 gene = gene)
                              }else{
                                  analysis_class <- 'limma'
                                  test_results <- perform_limma(data_matrix = expr_matrix,
                                                                samples = cases,
                                                                controls = controls,
                                                                gene = gene)
                              }
                              # Return a row to later merge in results dataframe. 
                              gene <- as.character(gene)
                              gene_de_sum <- c(gene, test_results, analysis_class, length(cases), bias)
                              return(gene_de_sum)
                          }
        
        colnames(gtests) <- c('Hugo_Symbol', 'adj_pval', 'log2fc', 'method', 'cases', 'bias')
        write.csv(x=as.data.frame(gtests), file=where_to_save, row.names = FALSE)
        
            #TODO: Refactor this and add reasons. Ex: "Could not perform DEG because of xxx"
            }else{ save_dummy_dataframe(where_to_save) }
}else{ save_dummy_dataframe(where_to_save) }






