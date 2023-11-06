#!/usr/bin/env nextflow

/*
Coriell Institute for Medical Research
RNAseq Pipeline. Started November 2023.

Contributors:
Anthony Pompetti <apompetti@coriell.org>

Methodology adapted from:
Gennaro Calendo
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Define local params 
*/
params.outdir = "./results"
params.pubdir = "gen_se"
gtf = params.genomes ? params.genomes[ params.genome ].gtf ?:false : false
grl = params.genomes ? params.genomes[ params.genome ].grl ?:false : false
rd = params.genomes ? params.genomes[ params.genome ].rd ?:false : false

/*
Run R script to generate Summarized Experiement using salmon quants
*/
process GEN_SE {
    maxForks 4
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/gen_se", mode: 'copy'

    input:
    path('quants/*')

    output:
    tuple path("*.rds"), emit: se

    script:
    """
    #!/usr/bin/env Rscript
    
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("data.table library is necessary")
    }
    suppressPackageStartupMessages(library(data.table))

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR library is necessary")
    }
    if (!requireNamespace("tximport", quietly = TRUE)) {
        stop("tximport library is necessary")
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("rtracklayer library is necessary")
    }
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("GenomicRanges library is necessary")
    }
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
        stop("GenomicFeatures library is necessary")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment library is necessary")
    }

    quants <- "./quants/"
    gtf_file <- ${gtf}
    grl_file <- ${grl}
    rwd_file <- ${rd}
    outdir <- "."

    message("Loading transcript GTF file: ", gtf_file)
    gtf <- rtracklayer::import(gtf_file)
    message("Loading GRangesList file: ", grl_file)
    te_grl <- readRDS(grl_file)
    message("Loading TE rowData file: ", rwd_file)
    te_rowdata <- fread(rwd_file)
    message("Creating txdb from GTF...")
    txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gtf))
    message("Creating a GRangesList object for all transcripts by genes...")
    gene_grl <- GenomicFeatures::transcriptsBy(txdb, by = "gene")

    message("Importing gene-level counts with `tximport::tximport()`...")
    gtf_dt <- as.data.table(data.frame(gtf))
    tx2gene <- unique(gtf_dt[, .(transcript_id, gene_id)])[!is.na(transcript_id)]

    importTxQuants <- function(x) {
        dt <- fread(x)
        dt[Name %like% "^ENS"]
    }

    quant_files <- list.files(
        path = quants,
        pattern = "quant.sf",
        recursive = TRUE,
        full.names = TRUE
    )

    # Make colnames of tximport data same as edgeR::catchSalmon import
    names(quant_files) <- gsub("/quant.sf", "", quant_files, fixed = TRUE)

    txi <- tximport::tximport(
    files = quant_files,
    type = "salmon",
    countsFromAbundance = "scaledTPM",
    tx2gene = tx2gene,
    importer = importTxQuants,
    dropInfReps = TRUE
    )
    genes <- txi$counts
    genes <- genes[rowSums(genes) > 0, ]

    # Filter the gene grl for only those genes in the count matrix
    gene_grl <- gene_grl[rownames(genes)]

    # Check alignment
    stopifnot("Failed creating GRangesList for transcripts by genes" = 
            all(rownames(genes) == names(gene_grl)))

    # TE counts ---------------------------------------------------------------

    message("Importing TE-loci counts with `edgeR::catchSalmon()`...")
    dirs <- list.dirs(
        path = quants, 
        full.names = TRUE,
        recursive = FALSE
    )

    catch <- edgeR::catchSalmon(dirs, verbose = FALSE)
    scaled_counts <- with(catch, counts / annotation$Overdispersion)
    scaled_counts <- scaled_counts[rowSums(scaled_counts) > 0, ]

    message("Annotating TE counts...")
    te <- scaled_counts[!grepl("^ENS", rownames(scaled_counts)), ]

    # Subset the GRangesList and the rowData for only these hashes
    te_df <- te_rowdata[Hash %chin% rownames(te)]
    te_grl <- te_grl[rownames(te)]

    stopifnot("Failed subsetting GRangesList for TEs" = 
            all(rownames(te) == names(te_grl)))

    # Construct SummarizedExperiments -----------------------------------------

    # Check that the genes and te matrices have the same colnames
    stopifnot("colnames of gene matrix does not match colnames of TEs matrix" = 
                all(colnames(te) == colnames(genes)))

    message("Combining gene and hash-level counts...")
    counts <- rbind(genes, te)

    message("Combining GRanges for genes and hashes...")
    rowranges <- c(gene_grl, te_grl)
    stopifnot("Failed combining GRangesList objects" = 
            all(names(rowranges) == rownames(counts)))

    message("Combining rowData feature-level annotations...")
    gene_df <- data.table(gene_id = rownames(genes))
    message("Joining gene names onto Ensembl IDs...")
    gene_df <- merge(gene_df, unique(gtf_dt[, .(gene_id, gene_name)]), by = "gene_id", all.x = TRUE) 
    coding <- gtf_dt[type == "gene" & gene_type == "protein_coding", unique(gene_id)]
    gene_df[, isCoding := fifelse(gene_id %chin% coding, TRUE, FALSE)]
    message("Binding rowData annotations...")
    rowdata <- rbindlist(list(gene = gene_df, TE = te_df), idcol = "FeatureType", fill = TRUE)
    rowdata[, FeatureID := fcoalesce(gene_id, Hash)]
    setcolorder(rowdata, "FeatureID")
    setDF(rowdata, rownames = rowdata$FeatureID)
    message("Subsetting final row annotations for features present in count matrix...")
    rowdata <- rowdata[rownames(counts), ]
    stopifnot("Failed combining rowData information for genes and hashes" = 
            all(rownames(counts) == rownames(rowdata)))

    message("Creating SummarizedExperiment objects...")
    # RangedSummarizedExperiment
    se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowRanges = rowranges,
    metadata = list(catch = catch)
    )

    # SummarizedExperiment
    se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowdata,
    metadata = list(catch = catch)
    )

    message("Saving RangedSummarizedExperiment to: ", file.path(outdir, "RangedSummarizedExperiment.rds"))
    saveRDS(se, file.path(outdir, "RangedSummarizedExperiment.rds"))
    message("Saving SummarizedExperiment to: ", file.path(outdir, "SummarizedExperiment.rds"))
    saveRDS(se2, file.path(outdir, "SummarizedExperiment.rds"))
    message("Done.")
    """
}