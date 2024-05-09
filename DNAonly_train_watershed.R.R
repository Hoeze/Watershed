# %%
library(optparse)
library(arrow)
library(data.table)
library(dplyr)
library(yaml)

source("watershed.R")

# %%
# trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
# trace(what = "cat", where = getNamespace("base"), exit = flush.console, print = FALSE)

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
#########################################
# Command line arguments
#########################################
# arguments <- parse_args(OptionParser(usage = "%prog [options]", description="Watershed command line args",
# 	option_list=list(
# 		make_option(c("-t","--training_input"), default = NULL, help="The Watershed input file containing instances used to train the model [default %default]"),
# # 		make_option(c("-i","--prediction_input"), default = NULL, help="The Watershed input file containing instances to predict on"),
# 		make_option(c("-d","--number_dimensions"), default=1, help="The number of outlier variables."),
# 		make_option(c("-m","--model_name"), default="Watershed_exact", help="Name of model. Options are Watershed_exact, Watershed_approximate, and RIVER"),
# 		make_option(c("-p","--dirichlet_prior_parameter"), default=10, help="Parameter defining Dirichlet distribution the acts as a prior a Phi (the model parameters defining E|Z"),
# 		make_option(c("-l","--l2_prior_parameter"), default=.01, help="Parameter defining L2 (gaussian) distribution the acts as a prior on the parameters of the conditional random Field (the model parameters defining Z|G"),
# 		make_option(c("-o","--output_file"), default="watershed_prediction_object.rds", help="file name to save results to"),
# 		make_option(c("-b","--binary_pvalue_threshold"), default=.1, help="Absolute p-value threshold used to create binary outliers used for Genomic Annotation Model")
# 	)
# ))
# process command line args
# training_input_file <- arguments$training_input
# prediction_input_file <- arguments$prediction_input
# number_of_dimensions <- arguments$number_dimensions
model_name <- "Watershed_approximate" # arguments$model_name
pseudoc <- 10 # arguments$dirichlet_prior_parameter
lambda_init <- .01 # arguments$l2_prior_parameter
if (lambda_init == "NA") {
	lambda_init <- NA
}
# # output_stem <- arguments$output_prefix
# binary_pvalue_threshold <- arguments$binary_pvalue_threshold

# %%
#########################################
# Fixed Parameters
#########################################
# If arguments$l2_prior_parameter == NA, perform grid search over the following values of lambda to determine optimal lambda
lambda_costs <- c(.1, .01, 1e-3)
# If arguments$l2_prior_parameter == NA, Number of folds to be used in K-fold cross validation for Genomic annotation model ()
nfolds <- 5
# Parameters used for Variational Optimization (only applies if arguments$model_name=="Watershed_approximate")
vi_step_size <- .8
vi_threshold <- 1e-8

# %%
# Set seed for reproducability (feel free to change!)
set.seed(1)

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
learn_watershed_model_parameters_from_training_data <- function(data_input, model_name, pseudoc, lambda_init, binary_pvalue_threshold, lambda_costs, nfolds, vi_step_size, vi_threshold) {
	# Parse data_input for relevent data
	feat_all <- data_input$feat
	discrete_outliers_all <- data_input$outliers_discrete
	binary_outliers_all <- data_input$outliers_binary
    number_of_dimensions <- data_input$number_of_dimensions

	#######################################
	## Standardize Genomic Annotations (features)
	#######################################
	mean_feat <- apply(feat_all, 2, mean)
	sd_feat <- apply(feat_all, 2, sd)
 	feat_all <- scale(feat_all, center=mean_feat, scale=sd_feat)

  	#######################################
	## Fit Genomic Annotation Model (GAM)
	#######################################
    # print(feat_all)
	gam_data <- logistic_regression_genomic_annotation_model_cv(feat_all, binary_outliers_all, nfolds, lambda_costs, lambda_init)
	# Report optimal lambda learned from cross-validation data (if applicable)
	if (is.na(lambda_init)) {
		cat(paste0(nfolds,"-fold cross validation on GAM yielded optimal lambda of ", gam_data$lambda, "\n"))
	}

	#######################################
	### Initialize phi using GAM
	#######################################
	# Compute GAM Predictions on data via function in  CPP file ("independent_crf_exact_updates.cpp")
	gam_posterior_obj <- update_independent_marginal_probabilities_exact_inference_cpp(
        feat_all,
        binary_outliers_all,
        gam_data$gam_parameters$theta_singleton,
        gam_data$gam_parameters$theta_pair,
        gam_data$gam_parameters$theta,
        matrix(0,2,2),
        matrix(0,2,2),
        number_of_dimensions,
        choose(number_of_dimensions, 2),
        FALSE
    )
	gam_posteriors <- gam_posterior_obj$probability
	# Initialize Phi using GAM posteriors
	# ie. Compute MAP estimates of the coefficients defined by P(outlier_status| FR)
	phi_init <- map_phi_initialization(discrete_outliers_all, gam_posteriors, number_of_dimensions, pseudoc)

	#######################################
	### Fit Watershed Model
	#######################################
	watershed_model <- train_watershed_model(
        feat_all,
        discrete_outliers_all,
        phi_init,
        gam_data$gam_parameters$theta_pair,
        gam_data$gam_parameters$theta_singleton,
        gam_data$gam_parameters$theta,
        pseudoc,
        gam_data$lambda,
        number_of_dimensions,
        model_name,
        vi_step_size,
        vi_threshold
    )

	return(list(
        mean_feat=mean_feat,
        sd_feat=sd_feat,
        model_params=watershed_model,
        gam_model_params=gam_data
    ))
}

# %%
predict_watershed = function(watershed_object, prediction_data_input) {
    # Parse data_input for relevent data
    predictions_feat <- prediction_data_input$feat
    # predictions_discretized_outliers <- prediction_data_input$outliers_discrete
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
# watershed_data = load_watershed_data(
#     data_dt,
#     use_adjusted_pvalues=TRUE
#     pvalue_threshold=0.2
# )

# %%
dir.create(paste(snakemake@params$output_basedir, "fold", sep="/"), showWarnings=FALSE)

# %%
data_dt[, unique(`fold`)]

# %%
row_any <- function(dt, cols) {
  return(dt[, Reduce(`|`, .SD), .SDcols=cols])
}

# %%
row_all_na <- function(dt, cols) {
  return(dt[, Reduce(`&`, lapply(.SD, is.na)), .SDcols=cols])
}

# %%
train_fold = function(target_fold, subsample_training_negative_class=NULL) {
    pvalue_threshold=0.2
    use_adjusted_pvalues = FALSE
    # # limit dimensions and features for testing
    # num_dims=3
    # num_features=3
    num_dims=NULL
    num_features=NULL
    
    print(paste0("Reading fold ", target_fold, "..."))

    training_data = data_dt %>% filter(`fold` != .env$target_fold)
    if (! is.null(subsample_training_negative_class)) {
        print(paste0("Subsampling non-outliers of training data by '", subsample_training_negative_class, "x' the number of outliers..."))
        training_data$has_outlier=row_any(training_data, cols=colnames(training_data)[startsWith(colnames(training_data), "y_truecat")])
        num_has_outlier = length(which(training_data$has_outlier))
        num_sample = min(subsample_training_negative_class * num_has_outlier, length(training_data$has_outlier) - num_has_outlier)
        
        training_data = rbind(
            training_data[(`has_outlier` != TRUE) | is.na(`has_outlier`)][sample(.N, num_sample)],
            training_data[`has_outlier` == TRUE]
        )[order(`individual`, `gene`)]
    }
    training_watershed_data = load_watershed_data(
        training_data,
        use_adjusted_pvalues=use_adjusted_pvalues,
        pvalue_threshold=pvalue_threshold,
        num_features=num_features,
        num_dims=num_dims
    )

    testing_data = data_dt %>% filter(`fold` == .env$target_fold)
    testing_watershed_data = load_watershed_data(
        testing_data,
        use_adjusted_pvalues=use_adjusted_pvalues,
        pvalue_threshold=pvalue_threshold,
        num_features=num_features,
        num_dims=num_dims
    )

    #######################################
    ## Train Watershed model on training data
    #######################################
    print("Training watershed...")
    watershed_object <- learn_watershed_model_parameters_from_training_data(training_watershed_data, model_name, pseudoc, lambda_init, binary_pvalue_threshold, lambda_costs, nfolds, vi_step_size, vi_threshold)
    watershed_object$pred_names = colnames(training_watershed_data$outliers_binary)

    #######################################
    ## Calculate predictions
    #######################################
    train_predictions = cbind(
        training_data %>% select(`individual`, `gene`),
        as.data.table(predict_watershed(watershed_object, training_watershed_data))
    )
    # train_predictions

    test_predictions = cbind(
        testing_data %>% select(`individual`, `gene`),
        as.data.table(predict_watershed(watershed_object, testing_watershed_data))
    )
    # test_predictions

    #######################################
    ## Save predictions
    #######################################
    print(paste0("Writing to fold ", target_fold, "..."))
    
    fold_dir = paste(snakemake@params$output_basedir, "fold", target_fold, sep="/")
    dir.create(fold_dir, showWarnings=FALSE)
    
    saveRDS(watershed_object, paste(fold_dir, "model.rds", sep="/"))
    write_parquet(train_predictions, paste(fold_dir, "train_predictions.parquet", sep="/"))
    write_parquet(test_predictions, paste(fold_dir, "test_predictions.parquet", sep="/"))
    
    return(list(
        training_data=training_data,
        testing_data=testing_data,
        watershed_object=watershed_object,
        training_predictions=train_predictions,
        testing_predictions=test_predictions
    ))
}

# %%
subsample_training_negative_class = 10

# %%
folds = data_dt[["fold"]] %>% unique()
folds

# %%
train_fold('0', subsample_training_negative_class=subsample_training_negative_class)

# %%
train_fold('1', subsample_training_negative_class=subsample_training_negative_class)

# %%
train_fold('2', subsample_training_negative_class=subsample_training_negative_class)

# %%
train_fold('3', subsample_training_negative_class=subsample_training_negative_class)

# %%
train_fold('4', subsample_training_negative_class=subsample_training_negative_class)

# %%
train_fold('test', subsample_training_negative_class=subsample_training_negative_class)

# %%
