# %%
library(optparse)
library(arrow)
library(data.table)
library(dplyr)
library(yaml)

source("watershed.R")

# %%

######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/fset@watershed/config.yaml', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/fset@watershed/data.parquet', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/fset@watershed/data.parquet.done', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/subset/dna_only.parquet', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/subset/dna_only.parquet.done', "featureset_config" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/fset@watershed/config.yaml', "data_pq" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/fset@watershed/data.parquet', "data_pq_done" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/fset@watershed/data.parquet.done', "subset_pq" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/subset/dna_only.parquet', "subset_pq_done" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/subset/dna_only.parquet.done'),
    output = list('/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/train.parquet', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/test.parquet', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/model.joblib', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/features.yaml', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/done', "train_data_dir" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/train.parquet', "test_data_dir" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/test.parquet', "model_joblib" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/model.joblib', "features_yaml" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/features.yaml', "touch_file" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm/done'),
    params = list('/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/feature_sets', '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models', 'train_simplecv.py', "output_basedir" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models/dna_only/watershed@train_simplecv.py#lightgbm', "feature_dir" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/feature_sets', "model_dir" = '/s/project/rep/processed/training_results_v14/gtex_v8_old_dna/models', "nb_script" = 'train_simplecv.py'),
    wildcards = list('gtex_v8_old_dna', 'dna_only', 'watershed', 'lightgbm', "ds_dir" = 'gtex_v8_old_dna', "subset" = 'dna_only', "feature_set" = 'watershed', "model_type" = 'lightgbm'),
    threads = 1,
    log = list(),
    resources = list('tmpdir', 'ntasks', 'mem_mb', 'mem_mib', "tmpdir" = '/scratch/tmp/hoelzlwi', "ntasks" = 1, "mem_mb" = 120000, "mem_mib" = 114441),
    config = list("projectTitle" = 'REP', "n_splits" = 6, "dirs" = list("CACHE_DIR" = '/s/project/rep/cache', "RAW_DATA_DIR" = '/s/project/rep/raw', "PROCESSED_DATA_DIR" = '/s/project/rep/processed', "MODEL_DIR" = '/s/project/rep/processed/training_results_v14'), "canonical_transcript_pq" = '/s/project/rep/processed/training_results_v14/general/gtex_canonical_transcript.parquet', "expressed_genes_pq" = '/s/project/rep/processed/training_results_v14/gtex_v8/gene_is_expressed_per_subtissue.parquet', "system" = list("MANE_select_url" = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz', "gencode_basepath" = '/s/genomes/Gencode/Gencode_human/release_{gencode_version}', "vep" = list("version" = 108, "vep_bin" = 'vep', "perl_bin" = 'perl', "vep_cache_dir" = '/s/raw/vep/{vep_version}', "cadd_dir" = '/s/raw/cadd/v1.6/{human_genome_assembly}', "loftee_data_dir" = '/s/raw/loftee/{human_genome_assembly}', "loftee_src_path" = '/s/raw/loftee/{human_genome_assembly}_src'), "absplice" = list("cache" = c('/s/project/rep/processed/absplice_cache/{human_genome_version}/veff.snp.parquet')))),
    rule = 'models__train_simplecv',
    bench_iteration = as.numeric(NA),
    scriptdir = '/data/nasif12/home_if12/hoelzlwi/Projects/REP/rep_scripts/smk/ds_dir/models',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########


# %%
snakemake@input

# %%
snakemake@output

# %%
featureset_config = read_yaml(snakemake@input$featureset_config)
featureset_config

# %%
data_dt = (
    open_dataset(snakemake@input$subset_pq)
    %>% select(`individual`, `gene`)
    %>% unique()
    %>% collect()
    %>% left_join(
        open_dataset(snakemake@input$data_pq) %>% collect(),
        by=c("individual", "gene")
    )
    # %>% select(`subtissue`, `tissue`, `hilo_padj`)
    # %>% rename(tissue_type=tissue, tissue=subtissue)
    %>% as.data.table()
)
data_dt

# %%
# remove infinite values
data_dt$`feature.vep.features.sift_score.pval_max_significant`[is.infinite(data_dt$`feature.vep.features.sift_score.pval_max_significant`)] = 200

# %%
data_dt[, unique(`feature.vep.features.TF_binding_site_variant.max`)]

# %%
data_dt[, unique(`feature.vep.features.regulatory_region_variant.max`)]

# %%
data_dt[, unique(`feature.vep.features.intergenic_variant.max`)]

# %%
data_dt = data_dt[, !c(
    "feature.vep.features.TF_binding_site_variant.max",
    "feature.vep.features.regulatory_region_variant.max",
    "feature.vep.features.intergenic_variant.max"
)]

# %%
# Load in and parse Watershed input file
load_watershed_data <- function(data_dt, use_adjusted_pvalues=FALSE, pvalue_threshold=0.05, num_dims=NULL, num_features=NULL) {
    target_dt = data_dt %>% select(all_of(colnames(data_dt)[startsWith(colnames(data_dt), "y_truecat")]))
    pval_dt = data_dt %>% select(all_of(colnames(data_dt)[startsWith(colnames(data_dt), "pval")]))
    padj_dt = data_dt %>% select(all_of(colnames(data_dt)[startsWith(colnames(data_dt), "padj")]))
    feat_dt = data_dt %>% select(all_of(colnames(data_dt)[startsWith(colnames(data_dt), "feature")]))
    
    if (is.null(num_dims)) {
        number_of_dimensions = length(colnames(pval_dt))
    } else{
        number_of_dimensions = num_dims
        target_dt = target_dt[, 1:num_dims]
        pval_dt = pval_dt[, 1:num_dims]
        padj_dt = padj_dt[, 1:num_dims]
    }
    if (! is.null(num_features)) {
        feat_dt = feat_dt[, 1:num_features]
    }

    rownames = data_dt[,paste(`individual`, ":", `gene`, sep="")]

    feat = as.matrix(feat_dt)
    feat[is.na(feat)] = 0
    feat[is.nan(feat)] = 0
    rownames(feat) <- rownames

    if (use_adjusted_pvalues) {
        print("Using adjusted p-values...")
        # use adjusted p-values
        outlier_pvalues = as.matrix(padj_dt)
    } else {
        # use nominal p-values
        print("Using nominal p-values...")
        outlier_pvalues = as.matrix(pval_dt)
    }
    rownames(outlier_pvalues) <- rownames

    outliers_binary <- ifelse(abs(outlier_pvalues) <= pvalue_threshold, 1, 0)
    # Convert outlier status into discretized random variables
    outliers_discrete <- get_discretized_outliers(outlier_pvalues)
    # Put all data into compact data structure
    data_input <- list(
        feat=as.matrix(feat),
        outliers_binary=as.matrix(outliers_binary),
        outliers_discrete=outliers_discrete,
        number_of_dimensions=number_of_dimensions
    )
    return(data_input)
}


# %%
predict_watershed = function(watershed_object, prediction_data_input) {
    # Parse data_input for relevent data
    predictions_feat <- prediction_data_input$feat
    predictions_discretized_outliers <- prediction_data_input$outliers_discrete
    # Scale prediction features (according to mean and standard deviation from training data)
    predictions_feat <- scale(predictions_feat, center=watershed_object$mean_feat, scale=watershed_object$sd_feat)
    number_dimensions <- prediction_data_input$number_of_dimensions


    ########################
    ## Inference to compute Watershed posterior probabilities
    ########################
    # watershed_info <- update_marginal_posterior_probabilities(predictions_feat, predictions_discretized_outliers, watershed_object$model_params)
    # watershed_posteriors <- watershed_info$probability  # Marginal posteriors

    # We do not need Expression for inference, therefore replace `predictions_discretized_outliers` with an empty `matrix()`.
    watershed_z_given_g = update_conditional_z_given_g_probabilities(predictions_feat, matrix(), watershed_object$model_params)
    watershed_z_given_g_posteriors = watershed_z_given_g$probability
    colnames(watershed_z_given_g_posteriors) = watershed_object$pred_names

    # ########################
    # # Add row names and column names to posterior predictions matrix
    # ########################
    # watershed_z_given_g_posteriors_mat <- cbind(rownames(predictions_feat), watershed_z_given_g_posteriors)
    # colnames(watershed_z_given_g_posteriors_mat) = c("sample_names", paste0("watershed_z_given_g_posterior_outlier_signal_", 1:number_dimensions))
    
    return(watershed_z_given_g_posteriors)
}

# %%
row_any <- function(dt, cols) {
  return(dt[, Reduce(`|`, .SD), .SDcols=cols])
}

# %%
row_all_na <- function(dt, cols) {
  return(dt[, Reduce(`&`, lapply(.SD, is.na)), .SDcols=cols])
}

# %%
folds = data_dt[["fold"]] %>% unique()
folds

# %%
fold_dirs = paste(snakemake@params$output_basedir, "fold", folds, sep="/")
fold_dirs

# %%
predict_fold = function(target_fold) {
    pvalue_threshold = 0.2
    fold_dir = paste(snakemake@params$output_basedir, "fold", target_fold, sep="/")
    
    print(paste0("Reading fold ", target_fold, "..."))

    testing_data = data_dt %>% filter(`fold` == .env$target_fold)
    testing_watershed_data = load_watershed_data(
        testing_data,
        use_adjusted_pvalues=FALSE,
        pvalue_threshold=pvalue_threshold
    )

    #######################################
    ## Train Watershed model on training data
    #######################################
    print("Loading watershed object...")
    watershed_object <- readRDS(paste(fold_dir, "model.rds", sep="/"))

    print("Predicting watershed scores...")
    test_predictions = cbind(
        testing_data %>% select(`individual`, `gene`),
        as.data.table(predict_watershed(watershed_object, testing_watershed_data))
    ) %>% mutate(fold = target_fold)
    # test_predictions
    
    return(test_predictions)
}

# %%
predictions_dt = rbindlist(lapply(folds, predict_fold))
predictions_dt

# %%
write_dataset(predictions_dt, paste0(snakemake@params$output_basedir, "/validation_predictions.parquet"), partitioning="fold")

# %% [raw]
# predictions_dt = (
#     open_dataset(paste0(snakemake@params$output_basedir, "/validation_predictions.parquet"))
#     %>% collect()
#     %>% as.data.table()
# )
# predictions_dt

# %%
colnames(predictions_dt)

# %%
column_mapping <- c(
    "pval.Adipose.20..20Subcutaneous" = "Adipose - Subcutaneous",
    "pval.Adipose.20..20Visceral.20.28Omentum.29" = "Adipose - Visceral (Omentum)",
    "pval.Adrenal.20Gland" = "Adrenal Gland",
    "pval.Artery.20..20Aorta" = "Artery - Aorta",
    "pval.Artery.20..20Coronary" = "Artery - Coronary",
    "pval.Artery.20..20Tibial" = "Artery - Tibial",
    "pval.Brain.20..20Amygdala" = "Brain - Amygdala",
    "pval.Brain.20..20Anterior.20cingulate.20cortex.20.28BA24.29" = "Brain - Anterior cingulate cortex (BA24)",
    "pval.Brain.20..20Caudate.20.28basal.20ganglia.29" = "Brain - Caudate (basal ganglia)",
    "pval.Brain.20..20Cerebellar.20Hemisphere" = "Brain - Cerebellar Hemisphere",
    "pval.Brain.20..20Cerebellum" = "Brain - Cerebellum",
    "pval.Brain.20..20Cortex" = "Brain - Cortex",
    "pval.Brain.20..20Frontal.20Cortex.20.28BA9.29" = "Brain - Frontal Cortex (BA9)",
    "pval.Brain.20..20Hippocampus" = "Brain - Hippocampus",
    "pval.Brain.20..20Hypothalamus" = "Brain - Hypothalamus",
    "pval.Brain.20..20Nucleus.20accumbens.20.28basal.20ganglia.29" = "Brain - Nucleus accumbens (basal ganglia)",
    "pval.Brain.20..20Putamen.20.28basal.20ganglia.29" = "Brain - Putamen (basal ganglia)",
    "pval.Brain.20..20Spinal.20cord.20.28cervical.20c.1.29" = "Brain - Spinal cord (cervical c-1)",
    "pval.Brain.20..20Substantia.20nigra" = "Brain - Substantia nigra",
    "pval.Breast.20..20Mammary.20Tissue" = "Breast - Mammary Tissue",
    "pval.Cells.20..20Cultured.20fibroblasts" = "Cells - Cultured fibroblasts",
    "pval.Cells.20..20EBV.transformed.20lymphocytes" = "Cells - EBV-transformed lymphocytes",
    "pval.Colon.20..20Sigmoid" = "Colon - Sigmoid",
    "pval.Colon.20..20Transverse" = "Colon - Transverse",
    "pval.Esophagus.20..20Gastroesophageal.20Junction" = "Esophagus - Gastroesophageal Junction",
    "pval.Esophagus.20..20Mucosa" = "Esophagus - Mucosa",
    "pval.Esophagus.20..20Muscularis" = "Esophagus - Muscularis",
    "pval.Heart.20..20Atrial.20Appendage" = "Heart - Atrial Appendage",
    "pval.Heart.20..20Left.20Ventricle" = "Heart - Left Ventricle",
    "pval.Liver" = "Liver",
    "pval.Lung" = "Lung",
    "pval.Minor.20Salivary.20Gland" = "Minor Salivary Gland",
    "pval.Muscle.20..20Skeletal" = "Muscle - Skeletal",
    "pval.Nerve.20..20Tibial" = "Nerve - Tibial",
    "pval.Ovary" = "Ovary",
    "pval.Pancreas" = "Pancreas",
    "pval.Pituitary" = "Pituitary",
    "pval.Prostate" = "Prostate",
    "pval.Skin.20..20Not.20Sun.20Exposed.20.28Suprapubic.29" = "Skin - Not Sun Exposed (Suprapubic)",
    "pval.Skin.20..20Sun.20Exposed.20.28Lower.20leg.29" = "Skin - Sun Exposed (Lower leg)",
    "pval.Small.20Intestine.20..20Terminal.20Ileum" = "Small Intestine - Terminal Ileum",
    "pval.Spleen" = "Spleen",
    "pval.Stomach" = "Stomach",
    "pval.Testis" = "Testis",
    "pval.Thyroid" = "Thyroid",
    "pval.Uterus" = "Uterus",
    "pval.Vagina" = "Vagina",
    "pval.Whole.20Blood" = "Whole Blood"
)


# %%
melted_predictions_dt <- melt(
    predictions_dt,
    id.vars = c("individual", "gene", "fold"),
    variable.name = "subtissue",
    value.name = "y_pred"
)
# Rename the 'subtissue' column using the mapping defined above
melted_predictions_dt$subtissue <- column_mapping[melted_predictions_dt$subtissue]
# convert p-values to -log10(p)
melted_predictions_dt$y_pred_proba <- -log10(melted_predictions_dt$y_pred)
melted_predictions_dt

# %%
column_mapping_y_true <- c(
    "y_truecat.Adipose.20..20Subcutaneous" = "Adipose - Subcutaneous",
    "y_truecat.Adipose.20..20Visceral.20.28Omentum.29" = "Adipose - Visceral (Omentum)",
    "y_truecat.Adrenal.20Gland" = "Adrenal Gland",
    "y_truecat.Artery.20..20Aorta" = "Artery - Aorta",
    "y_truecat.Artery.20..20Coronary" = "Artery - Coronary",
    "y_truecat.Artery.20..20Tibial" = "Artery - Tibial",
    "y_truecat.Brain.20..20Amygdala" = "Brain - Amygdala",
    "y_truecat.Brain.20..20Anterior.20cingulate.20cortex.20.28BA24.29" = "Brain - Anterior cingulate cortex (BA24)",
    "y_truecat.Brain.20..20Caudate.20.28basal.20ganglia.29" = "Brain - Caudate (basal ganglia)",
    "y_truecat.Brain.20..20Cerebellar.20Hemisphere" = "Brain - Cerebellar Hemisphere",
    "y_truecat.Brain.20..20Cerebellum" = "Brain - Cerebellum",
    "y_truecat.Brain.20..20Cortex" = "Brain - Cortex",
    "y_truecat.Brain.20..20Frontal.20Cortex.20.28BA9.29" = "Brain - Frontal Cortex (BA9)",
    "y_truecat.Brain.20..20Hippocampus" = "Brain - Hippocampus",
    "y_truecat.Brain.20..20Hypothalamus" = "Brain - Hypothalamus",
    "y_truecat.Brain.20..20Nucleus.20accumbens.20.28basal.20ganglia.29" = "Brain - Nucleus accumbens (basal ganglia)",
    "y_truecat.Brain.20..20Putamen.20.28basal.20ganglia.29" = "Brain - Putamen (basal ganglia)",
    "y_truecat.Brain.20..20Spinal.20cord.20.28cervical.20c.1.29" = "Brain - Spinal cord (cervical c-1)",
    "y_truecat.Brain.20..20Substantia.20nigra" = "Brain - Substantia nigra",
    "y_truecat.Breast.20..20Mammary.20Tissue" = "Breast - Mammary Tissue",
    "y_truecat.Cells.20..20Cultured.20fibroblasts" = "Cells - Cultured fibroblasts",
    "y_truecat.Cells.20..20EBV.transformed.20lymphocytes" = "Cells - EBV-transformed lymphocytes",
    "y_truecat.Colon.20..20Sigmoid" = "Colon - Sigmoid",
    "y_truecat.Colon.20..20Transverse" = "Colon - Transverse",
    "y_truecat.Esophagus.20..20Gastroesophageal.20Junction" = "Esophagus - Gastroesophageal Junction",
    "y_truecat.Esophagus.20..20Mucosa" = "Esophagus - Mucosa",
    "y_truecat.Esophagus.20..20Muscularis" = "Esophagus - Muscularis",
    "y_truecat.Heart.20..20Atrial.20Appendage" = "Heart - Atrial Appendage",
    "y_truecat.Heart.20..20Left.20Ventricle" = "Heart - Left Ventricle",
    "y_truecat.Liver" = "Liver",
    "y_truecat.Lung" = "Lung",
    "y_truecat.Minor.20Salivary.20Gland" = "Minor Salivary Gland",
    "y_truecat.Muscle.20..20Skeletal" = "Muscle - Skeletal",
    "y_truecat.Nerve.20..20Tibial" = "Nerve - Tibial",
    "y_truecat.Ovary" = "Ovary",
    "y_truecat.Pancreas" = "Pancreas",
    "y_truecat.Pituitary" = "Pituitary",
    "y_truecat.Prostate" = "Prostate",
    "y_truecat.Skin.20..20Not.20Sun.20Exposed.20.28Suprapubic.29" = "Skin - Not Sun Exposed (Suprapubic)",
    "y_truecat.Skin.20..20Sun.20Exposed.20.28Lower.20leg.29" = "Skin - Sun Exposed (Lower leg)",
    "y_truecat.Small.20Intestine.20..20Terminal.20Ileum" = "Small Intestine - Terminal Ileum",
    "y_truecat.Spleen" = "Spleen",
    "y_truecat.Stomach" = "Stomach",
    "y_truecat.Testis" = "Testis",
    "y_truecat.Thyroid" = "Thyroid",
    "y_truecat.Uterus" = "Uterus",
    "y_truecat.Vagina" = "Vagina",
    "y_truecat.Whole.20Blood" = "Whole Blood"
)


# %%
y_truecat_dt = (
    data_dt %>% select(all_of(c("individual", "gene", names(column_mapping_y_true))))
    %>% melt(
        id.vars = c("individual", "gene"),
        variable.name = "subtissue",
        value.name = "y_truecat"
    )
    %>% filter(! is.na(`y_truecat`))
)
# Rename the 'subtissue' column using the mapping defined above
y_truecat_dt$subtissue <- column_mapping_y_true[y_truecat_dt$subtissue]
y_truecat_dt

# %%
joint_predictions_dt = (
    y_truecat_dt
    %>% left_join(melted_predictions_dt, by=c("individual", "gene", "subtissue"))
)
joint_predictions_dt

# %%
print(paste0(snakemake@params$output_basedir, "/data.parquet"))

# %%
write_dataset(melted_predictions_dt, paste0(snakemake@params$output_basedir, "/data.parquet"), partitioning="fold")

# %%
