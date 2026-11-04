# =============================================================================
# MaAsLin 3 CLINICAL PIPELINE - TEMPLATE UNIVERSALE E GENERALIZZABILE
# =============================================================================
# Questo script implementa un'analisi longitudinale robusta e pubblicabile
# per dati di conteggio metagenomici (EC, KO, generi, specie, pathways) normalizzati
# o da normalizzare (es. su conteggi ORF o profondità di sequenziamento).
#
# Modello: MaAsLin 3 (Hurdle model: Abundance log2 lineare + Prevalence logit).
#
# QUESITI CLINICI SUPPORTATI (attivabili/disattivabili a monte):
#   - REQ 1 : Verifica sbilanciamento al basale (T0)
#   - REQ 5 : Interazione Gruppo x Tempo (Trial principale: Treated vs Control nel tempo)
#   - REQ 6 : Contrasto continuo tra due dosi/interventi (es. A vs B all'interno dei trattati)
#   - REQ 7 : Deconvoluzione multifattoriale (più interventi simultanei aggiustati reciprocamente)
#   - REQ 9 : Ciclo di sensibilità su covariate/farmaci opzionali (loop su pattern di colonne)
# =============================================================================

suppressPackageStartupMessages({
  library(maaslin3)
  library(phyloseq)
  library(tidyverse)
})

# =============================================================================
# ---- 0. PARAMETRI DI CONFIGURAZIONE (USER INPUT - MODIFICA SOLO QUI) --------
# =============================================================================

# --- 0.1 FILE DI INPUT E FORMATO DEI DATI ---
# File phyloseq (.Rdata con oggetto 'phyloseq' o data.frame sample_data)
# Se NULL, cerca automaticamente l'unico file .Rdata nella directory corrente
input_rdata_file <- NULL
phyloseq_obj_name <- "phy.ok"  # Nome dell'oggetto all'interno del file .Rdata (o NULL per primo oggetto)

# Matrice delle feature (campioni sulle righe, feature sulle colonne; oppure tsv/csv tab-delimited)
features_file <- "final_summary_identity_70.tsv"

# File opzionale di normalizzazione (es. ORF counts per campione, depth totale, ecc.)
# Impostare su NULL se i dati in features_file sono già normalizzati (es. abbondanze relative)
norm_counts_file <- "ORF_counts.txt"
norm_count_col   <- "ORF_count" # Nome della colonna dei conteggi di profondità
min_norm_depth   <- 1000        # Soglia minima di profondità per considerare valido il campione

# Cartella di output generale
out_dir <- "Maaslin3_General_Analysis_Output"

# --- 0.2 VARIABILI DI DISEGNO LONGITUDINALE ---
# Nome della colonna identificativa del soggetto per intercetta casuale: (1 | Subject)
col_subject <- "Subject"

# Nome della colonna del timepoint e livelli temporali ordinati (il primo è il baseline T0)
col_time <- "Time"
time_levels <- c("T0", "T1", "T2", "T3", "T4")

# Timepoint pre-specificato di fine intervento (per endpoint primario e washout)
primary_timepoint <- "T2"

# Variabile del braccio principale (per Req 1 e Req 5)
col_treatment <- "Treatment"
treatment_ref <- "Control"     # Livello di riferimento (controllo/placebo)
treatment_act <- "Treatment"     # Livello attivo/trattamento

# --- 0.3 COVARIATE CLINICHE DI AGGIUSTAMENTO ---
# Colonne cliniche da includere in tutti i modelli (es. "Age", "Sex", "BMI", "Smoking")
base_covariates <- c("Age", "Sex", "Bristol", "City", "PPI")

# Aggiungere il log10 della profondità come covariata di aggiustamento (raccomandato da MaAsLin 3)?
include_depth_covariate <- TRUE

# Covariata opzionale aggiuntiva per analisi di sensibilità
sensitivity_covariate <- "burden_now"

# Imputazione LOCF (Last Observation Carried Forward) per colonne con pattern specifico (es. farmaci)
# Impostare su NULL se non si desidera effettuare LOCF
locf_col_pattern <- "_L4$"

# --- 0.4 CONFIGURAZIONE RICHIESTE CLINICHE (SWITCH ON/OFF) ---
run_main_models <- TRUE  # Esegue Req 1 (basale) e Req 5 (Trial principale); il Dose-model ha il proprio switch (run_dose_model)
run_req6        <- TRUE  # Req 6: Contrasto continuo o bilancio tra 2 interventi nei trattati
run_req7        <- TRUE  # Req 7: Deconvoluzione multifattoriale (più interventi simultanei)
run_req9_loop   <- FALSE # Req 9: Ciclo di robustezza su covariate extra (es. classi di farmaci)
run_dose_model     <- TRUE  # Dose-response model (Time * composite adherence dose)
run_sensitivity    <- TRUE  # Sensitivity analyses of the primary contrast (prevalence, covariate, no median, T0-T2 only)
run_evidence_tiers <- TRUE  # Evidence tiers A/B/C per EC and model

# --- 0.5 PARAMETRI SPECIFICI PER REQ 6 (CONSTRASTO CONTINUO DOSE/INTERVENTO) ---
# Compilare solo se run_req6 = TRUE
# Due variabili numeriche continue di assunzione (es. grammi o aderenza di due componenti)
req6_component_A <- "Carciofini"
req6_component_B <- "Pane"
req6_subset_condition <- list(column = col_treatment, value = treatment_act)

# --- 0.6 PARAMETRI SPECIFICI PER REQ 7 (DECONVOLUZIONE MULTIFATTORIALE) ---
# Compilare solo se run_req7 = TRUE
# Elenco delle variabili dicotomiche o continue dei singoli interventi e relativi riferimenti
req7_interventions <- list(
  Pane_dic               = list(type = "factor", ref = "NoDiet"),
  Carciofini_dic         = list(type = "factor", ref = "NoDiet"),
  Training_cognitivo_dic = list(type = "factor", ref = "NoTraining"),
  Training_fisico_dic    = list(type = "factor", ref = "NoTraining")
)

# --- 0.7 PARAMETRI SPECIFICI PER REQ 9 (LOOP COVARIATE EXTRA) ---
req9_covar_pattern  <- "_L4$"
req9_min_pos_cases  <- 10

# --- 0.8 SOGLIE STATISTICHE E DI FILTRO ---
min_prevalence_primary <- 0.20 # Prevalenza minima per il modello primario
min_prevalence_sens    <- 0.40 # Prevalenza minima per sensitivity
alpha_primary          <- 0.05 # Soglia FDR primaria (BH)
alpha_support          <- 0.10 # Soglia FDR di supporto per sensitivity
top_n_plot             <- 15   # Numero di feature da includere nei barplot

# Filtri di stabilità per il ramo logistico di MaAsLin 3
se_max_logit   <- 2.0
coef_max_logit <- 5.0
max_prev_logit <- 0.95

# Parametri operativi / debug
max_features_test <- 50      # Impostare ad un intero (es. 50) per un test run rapido, altrimenti NULL
random_seed       <- 2026

set.seed(random_seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# ---- 1. CARICAMENTO E PRE-PROCESSING DINAMICO DEI METADATI E DATI -----------
# =============================================================================
cat(">>> [1/7] Inizializzazione e caricamento dati...\n")

# Caricamento Rdata / phyloseq
if (is.null(input_rdata_file)) {
  found_rdata <- list.files(".", pattern = "\\.rdata$", ignore.case = TRUE)
  if (length(found_rdata) != 1) {
    stop("Impossibile rilevare automaticamente un unico file .Rdata. Specificare 'input_rdata_file' nella Sezione 0.")
  }
  input_rdata_file <- found_rdata[1]
}

temp_env <- new.env()
load(input_rdata_file, envir = temp_env)
if (!is.null(phyloseq_obj_name) && phyloseq_obj_name %in% ls(temp_env)) {
  phy_obj <- temp_env[[phyloseq_obj_name]]
} else {
  phy_obj <- temp_env[[ls(temp_env)[1]]]
}
rm(temp_env)

# Estrazione data.frame dei metadati
if (inherits(phy_obj, "phyloseq")) {
  meta <- as(sample_data(phy_obj), "data.frame")
} else if (is.data.frame(phy_obj)) {
  meta <- phy_obj
} else {
  stop("L'oggetto caricato non è né un phyloseq né un data.frame valido.")
}
meta$Sample <- rownames(meta)

if ("Bristol" %in% names(meta)) {
  meta$Bristol <- as.numeric(as.character(meta$Bristol))
}
if ("Age" %in% names(meta)) {
  meta$Age <- as.numeric(as.character(meta$Age))
}
# numeric columns can be stored as factor/character inside sample_data: coerce explicitly
for (v in intersect(c("Pane", "Carciofini", "burden_now"), names(meta))) meta[[v]] <- as.numeric(as.character(meta[[v]]))
for (v in grep("_L4$", names(meta), value = TRUE)) meta[[v]] <- as.numeric(as.character(meta[[v]]))

# Caricamento tabelle conteggi
ec_raw <- read.delim(features_file, row.names = 1, check.names = FALSE)

# Normalizzazione su fattore di scala/profondità (se fornito)
if (!is.null(norm_counts_file) && file.exists(norm_counts_file)) {
  norm_df <- read.delim(norm_counts_file, row.names = 1, check.names = FALSE)
  valid_s <- rownames(norm_df)[norm_df[[norm_count_col]] >= min_norm_depth]
  common_s <- Reduce(intersect, list(rownames(ec_raw), rownames(meta), rownames(norm_df), valid_s))
  
  meta   <- meta[common_s, , drop = FALSE]
  meta   <- meta[meta[[col_time]] %in% time_levels, , drop = FALSE]
  common_s <- rownames(meta)
  
  scale_vector <- norm_df[common_s, norm_count_col]
  data_in <- as.data.frame(as.matrix(ec_raw[common_s, , drop = FALSE]) / scale_vector)
  meta$log10_Depth <- log10(scale_vector)
} else {
  common_s <- intersect(rownames(ec_raw), rownames(meta))
  meta <- meta[common_s, , drop = FALSE]
  meta <- meta[meta[[col_time]] %in% time_levels, , drop = FALSE]
  common_s <- rownames(meta)
  data_in <- as.data.frame(as.matrix(ec_raw[common_s, , drop = FALSE]))
  meta$log10_Depth <- 0
}

cat("Campioni validi allineati:", length(common_s), "\n")

# Sanificazione nomi delle feature
safe_names <- make.names(colnames(data_in), unique = TRUE)
name_dict  <- setNames(colnames(data_in), safe_names)
colnames(data_in) <- safe_names

if (!is.null(max_features_test)) {
  data_in <- data_in[, seq_len(min(max_features_test, ncol(data_in))), drop = FALSE]
}
cat("Feature microbiologiche/funzionali in analisi:", ncol(data_in), "\n")

# Formattazione fattori temporali e bracci
meta[[col_time]] <- factor(meta[[col_time]], levels = intersect(time_levels, unique(meta[[col_time]])))
later_timepoints <- levels(meta[[col_time]])[-1]

if (col_treatment %in% names(meta)) {
  meta[[col_treatment]] <- factor(
    ifelse(meta[[col_treatment]] == treatment_ref, treatment_ref, treatment_act),
    levels = c(treatment_ref, treatment_act)
  )
}

# composite dose (0-1): mean adherence to bread and artichokes per subject, controls = 0
# (used by the dose-response model)
if (all(c("Pane", "Carciofini") %in% names(meta)) && col_treatment %in% names(meta)) {
  is_ctrl <- meta[[col_treatment]] == treatment_ref
  pn <- ifelse(is_ctrl, 0, meta$Pane); ca <- ifelse(is_ctrl, 0, meta$Carciofini)
  meta$dose <- (ave(pn, meta[[col_subject]], FUN = function(x) mean(x, na.rm = TRUE)) +
                ave(ca, meta[[col_subject]], FUN = function(x) mean(x, na.rm = TRUE))) / 200
}

# LOCF opzionale
if (!is.null(locf_col_pattern) && nchar(locf_col_pattern) > 0) {
  matched_locf <- grep(locf_col_pattern, colnames(meta), value = TRUE)
  if (length(matched_locf) > 0) {
    ord <- order(meta[[col_subject]], meta[[col_time]])
    for (ac in matched_locf) {
      x <- meta[[ac]]
      for (i in ord[-1]) {
        j <- ord[match(i, ord) - 1]
        if (is.na(x[i]) && meta[[col_subject]][i] == meta[[col_subject]][j]) x[i] <- x[j]
      }
      meta[[ac]] <- x
    }
  }
}

# Configurazione finale delle covariate
selected_covs <- base_covariates[base_covariates %in% names(meta)]
if (include_depth_covariate && "log10_Depth" %in% names(meta) && var(meta$log10_Depth, na.rm = TRUE) > 0) {
  selected_covs <- c(selected_covs, "log10_Depth")
}
covs_txt <- paste(selected_covs, collapse = " + ")

formula_arm <- paste("~", col_time, "*", col_treatment, "+", covs_txt, "+ (1 |", col_subject, ")")
cat("Formula principale (Trial):", formula_arm, "\n")

# =============================================================================
# ---- 2. HELPERS STATISTICI E MOTORE CONTRASTI --------------------------------
# =============================================================================
bh_na <- function(p) { 
  q <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  q[ok] <- p.adjust(p[ok], "BH")
  q 
}

# Testing FAMILIES (defined a priori; BH is applied within a family, pooled over the two models):
#   primary                  : ix at the primary time point (the only confirmatory family)
#   secondary_timecourse     : ix at the other time points + washout contrasts
#   secondary_crosssectional : exposure effect at each time point (time-stratified)
#   baseline_check           : exposure effect at T0
#   one primary + one secondary family per component in multi-component models (Req 7)
fam_of <- function(test) {
  tk_of <- function(x) sub("^[a-z]+_(T[0-9]+).*$", "\\1", x)
  ifelse(test == "base_diff", "baseline_check",
  ifelse(grepl("^ix_T[0-9]+_", test),
         paste0(sub("^ix_T[0-9]+_", "", test), ifelse(tk_of(test) == primary_timepoint, "|primary", "|secondary")),
  ifelse(grepl("^ix_", test), ifelse(tk_of(test) == primary_timepoint, "primary", "secondary_timecourse"),
  ifelse(grepl("^diff_", test), "secondary_crosssectional",
  ifelse(grepl("^wash_", test), "secondary_timecourse", "other")))))
}

add_flags <- function(ct, analysis = "analysis") {
  ct$unstable <- ct$model == "logistic" &
    (!is.na(ct$error) | is.na(ct$stderr) | ct$stderr > se_max_logit | abs(ct$coef) > coef_max_logit |
       (ct$N_not_zero / ct$N) > max_prev_logit)
  ct$unstable[is.na(ct$unstable)] <- TRUE
  p <- ifelse(ct$unstable, NA_real_, ct$pval_individual)
  ct$family   <- paste(analysis, fam_of(ct$test), sep = "|")
  ct$q_within <- ave(p, paste(ct$test, ct$model), FUN = bh_na)   # within contrast and model (lenient)
  ct$q_clean  <- ave(p, ct$test, FUN = bh_na)                    # within contrast, models pooled
  ct$q_family <- ave(p, ct$family, FUN = bh_na)                  # within family, models pooled (use this)
  # Inferential status (a priori): only the primary contrast of the main trial model is confirmatory;
  # sensitivity analyses are supportive; everything else is exploratory (no control across analyses)
  ct$inference <- ifelse(analysis == "main" & ct$test == paste0("ix_", primary_timepoint), "confirmatory",
                  ifelse(grepl("^sens_", analysis), "sensitivity", "exploratory"))
  ct
}

run_fit <- function(formula_str, run_name, min_prev, md = meta, data = data_in) {
  maaslin3(input_data = data[rownames(md), , drop = FALSE], input_metadata = md,
           output = file.path(out_dir, run_name), formula = formula_str,
           normalization = "NONE", transform = "LOG",
           augment = TRUE, standardize = TRUE,
           max_significance = alpha_support,
           median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
           min_prevalence = min_prev, save_models = TRUE,
           plot_summary_plot = FALSE, plot_associations = FALSE,
           max_pngs = 20, verbosity = "WARN")
}

coef_names <- function(fit) {
  unique(c(fit$fit_data_abundance$results$name, fit$fit_data_prevalence$results$name))
}

spec_to_matrix <- function(S) {
  cn <- unique(unlist(lapply(S, names)))
  M <- matrix(0, nrow = length(S), ncol = length(cn), dimnames = list(names(S), cn))
  for (r in names(S)) M[r, names(S[[r]])] <- S[[r]]
  M
}

run_contrasts <- function(fit, S, median_abund = TRUE, analysis = "analysis") {
  ct <- maaslin_contrast_test(maaslin3_fit = fit, contrast_mat = spec_to_matrix(S),
                              max_significance = alpha_support, correction = "BH",
                              median_comparison_abundance = median_abund,
                              median_comparison_prevalence = FALSE)
  if (!is.data.frame(ct)) ct <- rbind(ct$fit_data_abundance$results, ct$fit_data_prevalence$results)
  ct <- as.data.frame(ct)
  ct <- add_flags(ct, analysis)
  ct$feature_name <- name_dict[ct$feature]
  ct
}

# Generatore automatico della matrice contrasti per Termine * Time
build_generic_spec <- function(cn, exposure_term, later_tps, time_col_name = col_time) {
  need <- c(exposure_term)
  for (tk in later_tps) {
    term_time <- paste0(time_col_name, tk)
    term_int1 <- paste0(term_time, ":", exposure_term)
    term_int2 <- paste0(exposure_term, ":", term_time)
    if (!term_int1 %in% cn && !term_int2 %in% cn) {
      stop("Interazione non trovata nei coefficienti per: ", term_time, " e ", exposure_term)
    }
  }
  
  # Req 1: Effetto al baseline (T0)
  S <- list(base_diff = setNames(1, exposure_term))
  
  # Req 5: Interazioni longitudinali (variazione da T0 a Tk)
  for (tk in later_tps) {
    term_time <- paste0(time_col_name, tk)
    term_int <- if (paste0(term_time, ":", exposure_term) %in% cn) {
      paste0(term_time, ":", exposure_term)
    } else {
      paste0(exposure_term, ":", term_time)
    }
    S[[paste0("ix_", tk)]] <- setNames(1, term_int)
    S[[paste0("diff_", tk)]] <- setNames(c(1, 1), c(exposure_term, term_int))   # exposure effect at tk
  }
  
  # Washout relativo all'endpoint primario
  if (primary_timepoint %in% later_tps) {
    later_after_prim <- intersect(c("T3", "T4", "T5"), later_tps)
    for (tk in later_after_prim) {
      t_prim <- paste0(time_col_name, primary_timepoint)
      t_wash <- paste0(time_col_name, tk)
      int_p <- if (paste0(t_prim, ":", exposure_term) %in% cn) paste0(t_prim, ":", exposure_term) else paste0(exposure_term, ":", t_prim)
      int_w <- if (paste0(t_wash, ":", exposure_term) %in% cn) paste0(t_wash, ":", exposure_term) else paste0(exposure_term, ":", t_wash)
      S[[paste0("wash_", tk, "_vs_", primary_timepoint)]] <- setNames(c(1, -1), c(int_w, int_p))
    }
  }
  S
}

# =============================================================================
# ---- 3. HELPER EXPORT TSV E GRAFICI GGPLOT2 (STILE GASLINI) -----------------
# =============================================================================
export_tsv_and_plots <- function(ct_df, raw_mat, sub_meta, out_subdir, title_prefix = "Analysis") {
  dir_tsv  <- file.path(out_subdir, "tsv_results")
  dir_plot <- file.path(out_subdir, "plots")
  dir.create(dir_tsv,  showWarnings = FALSE, recursive = TRUE)
  
  # Calcolo statistiche conteggi reali
  raw_sub  <- raw_mat[rownames(sub_meta), , drop = FALSE]
  raw_n_pos <- colSums(raw_sub > 0, na.rm = TRUE)
  raw_tot   <- colSums(raw_sub, na.rm = TRUE)
  raw_mean  <- round(colMeans(raw_sub, na.rm = TRUE), 2)
  
  ct_df$N_Samples_Scenario <- nrow(sub_meta)
  ct_df$N_Samples_Pos      <- raw_n_pos[ct_df$feature_name]
  ct_df$Total_Raw_Counts   <- raw_tot[ct_df$feature_name]
  ct_df$Mean_Raw_Counts    <- raw_mean[ct_df$feature_name]
  ct_df$qval_contrast      <- ct_df$q_within
  ct_df$qval_pooled        <- ct_df$q_clean
  ct_df$qval_family        <- ct_df$q_family   # q-value to report; read together with the 'inference' column
  
  write.table(ct_df, file.path(out_subdir, paste0("MASTER_", title_prefix, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Export singoli TSV e grafici per ciascun test e modello
  for (tt in unique(ct_df$test)) {
    for (mm in unique(ct_df$model)) {
      sub_df <- ct_df[ct_df$test == tt & ct_df$model == mm & !ct_df$unstable, ]
      if (nrow(sub_df) == 0) next
      sub_df <- sub_df[order(sub_df$qval_family), ]
      
      write.table(sub_df, file.path(dir_tsv, paste0("results_", tt, "_", mm, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      res_sig <- sub_df[!is.na(sub_df$qval_family) & sub_df$qval_family <= alpha_primary, ]
      if (nrow(res_sig) == 0) next      # plots only when at least one EC is significant
      
      pdata <- head(res_sig, top_n_plot)
      mod_label <- if (mm == "linear") "Abundance (log2 FC)" else "Prevalence (log-odds)"
      
      dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)   # created only when a plot is produced
      pdf_path <- file.path(dir_plot, paste0(tt, "_", mm, "_plots.pdf"))
      dyn_h <- max(4, 1.2 + (nrow(pdata) * 0.35))
      if (dyn_h > 15) dyn_h <- 15
      
      pdf(pdf_path, width = 9, height = dyn_h)
      
      # --- VOLCANO PLOT ---
      sub_df$color_status <- "NS"
      sub_df$color_status[sub_df$qval_family <= alpha_primary & sub_df$coef > 0] <- "Up"
      sub_df$color_status[sub_df$qval_family <= alpha_primary & sub_df$coef < 0] <- "Down"
      
      p_volc <- ggplot(sub_df, aes(x = coef, y = -log10(qval_family))) +
        geom_point(aes(color = color_status), alpha = 0.65, size = 2) +
        scale_color_manual(values = c("Up" = "#084594", "Down" = "#41B6C4", "NS" = "grey70"),
                           labels = c("Up" = "Positivamente associato", "Down" = "Negativamente associato", "NS" = "Non Sig."),
                           name = "Status") +
        geom_hline(yintercept = -log10(alpha_primary), linetype = "dashed", color = "red", alpha = 0.7) +
        labs(title = paste0("MaAsLin 3 Volcano: ", title_prefix, " [", tt, "]"),
             subtitle = paste0("Modello: ", mod_label, " | ", toupper(unique(sub_df$inference)[1])),
             x = paste0("Coefficiente (", mod_label, ")"), y = "-log10(q-value, BH within family)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
              plot.subtitle = element_text(hjust = 0.5, size = 8.5, face = "italic"),
              legend.position = "bottom")
      print(p_volc)
      
      if (nrow(pdata) > 0) {
        # --- BARPLOT ---
      pdata$feature_name <- factor(pdata$feature_name, levels = unique(pdata$feature_name[order(pdata$coef)]))
      p_bar <- ggplot(pdata, aes(x = feature_name, y = coef)) +
        geom_bar(stat = "identity", aes(fill = coef > 0), width = 0.7, color = "black", linewidth = 0.2) +
        geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2, color = "black", alpha = 0.7) +
        scale_fill_manual(values = c("TRUE" = "#084594", "FALSE" = "#41B6C4"),
                          labels = c("TRUE" = "Enriched in Target/Active", "FALSE" = "Enriched in Reference"),
                          name = "Arricchimento") +
        coord_flip() +
        labs(title = paste0("Top Features: ", title_prefix, " [", tt, "]"),
             subtitle = paste0("MaAsLin 3 (", mod_label, ")"),
             x = NULL, y = paste0("Coefficiente +/- SE (", mod_label, ")")) +
        theme_classic() +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_text(size = 8.5, face = "italic", color = "black"),
              plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 8.5, face = "italic"),
              legend.position = "bottom")
      print(p_bar)
      
      dev.off()
    }
  }
}
}   # chiusura di export_tsv_and_plots (mancante: errore di parsing)

# =============================================================================
# ---- 4. MODELLO PRINCIPALE (TRIAL LONGITUDINALE: REQ 1 & REQ 5) -------------
# =============================================================================
if (run_main_models) {
  cat("\n>>> [4/7] Esecuzione Modello Principale (Req 1 Basale & Req 5 Interazione)...\n")
  fit_arm <- run_fit(formula_arm, "main_arm", min_prevalence_primary)
  
  exp_term_name <- paste0(col_treatment, treatment_act)
  S_arm <- build_generic_spec(coef_names(fit_arm), exp_term_name, later_timepoints)
  ct_arm <- run_contrasts(fit_arm, S_arm, analysis = "main")
  
  export_tsv_and_plots(ct_arm, ec_raw, meta, file.path(out_dir, "req1_req5_trial_output"), "Trial_Arm")
  saveRDS(list(arm = ct_arm), file.path(out_dir, "main_trial_results.rds"))   # contrasts only: the full fit object can be several GB
  cat("Req 1 & Req 5 completate con successo.\n")
}

# =============================================================================
# ---- 4b. DOSE-RESPONSE MODEL (TIME * DOSE COMPOSITA DI ADERENZA) -------------
# =============================================================================
# dose (0-1) = aderenza media a pane e carciofini per soggetto, controlli = 0.
# Contrasti: interazione Time x dose (coefficiente per +1 SD, standardize = TRUE).
if (run_dose_model) {
  cat("\n>>> [4b] Dose-response model (Time * dose)...\n")
  if ("dose" %in% names(meta)) {
    tryCatch({
      meta_dose <- meta[!is.na(meta$dose), , drop = FALSE]
      f_dose <- paste("~", col_time, "* dose +", covs_txt, "+ (1 |", col_subject, ")")
      cat("Formula Dose:", f_dose, "\n")
      fit_dose <- run_fit(f_dose, "dose_response", min_prevalence_primary, md = meta_dose)
      S_dose   <- build_generic_spec(coef_names(fit_dose), "dose", later_timepoints)
      ct_dose  <- run_contrasts(fit_dose, S_dose, analysis = "dose")
      export_tsv_and_plots(ct_dose, ec_raw, meta_dose, file.path(out_dir, "dose_response_output"), "Dose_Response")
      saveRDS(list(dose = ct_dose), file.path(out_dir, "dose_response_results.rds"))   # contrasts only
      cat("Dose-response completato.\n")
    }, error = function(e) message("[Dose] failed: ", conditionMessage(e)))
  } else {
    cat("[SKIP Dose] Colonna 'dose' non disponibile (servono Pane e Carciofini).\n")
  }
}


# =============================================================================
# ---- 4c. ANALISI DI SENSIBILITA' DEL CONTRASTO PRIMARIO ----------------------
# =============================================================================
# Ricalcola solo il contrasto di interazione al timepoint primario in 4 varianti:
#   prev40   : prevalenza minima piu' stringente (min_prevalence_sens)
#   cov      : aggiustamento per sensitivity_covariate (casi completi)
#   nomedian : senza confronto con la mediana dei coefficienti (abbondanza)
#   T0_prim  : solo baseline e timepoint primario
# Richiede il modello principale (usa ct_arm, fit_arm e S_arm).
sens_list <- list()
if (run_sensitivity) {
  if (!exists("ct_arm")) {
    cat("[SKIP Sensitivity] Richiede run_main_models = TRUE.\n")
  } else {
    cat("\n>>> [4c] Analisi di sensibilita' del contrasto primario...\n")
    prim_id <- paste0("ix_", primary_timepoint)

    run_sens <- function(label, formula_str, min_prev, md, later_tps) {
      tryCatch({
        fit_s <- run_fit(formula_str, paste0("sens_", label), min_prev, md = md)
        S_s   <- build_generic_spec(coef_names(fit_s), exp_term_name, later_tps)
        ct_s  <- run_contrasts(fit_s, S_s[prim_id], analysis = paste0("sens_", label))
        ct_s$sensitivity <- label
        ct_s
      }, error = function(e) { message("[Sensitivity ", label, "] failed: ", conditionMessage(e)); NULL })
    }

    # (a) prevalenza minima piu' stringente
    sens_list[["prev40"]] <- run_sens("prev40", formula_arm, min_prevalence_sens, meta, later_timepoints)

    # (b) covariata aggiuntiva (se presente, non gia' nel modello e con variabilita')
    if (!is.null(sensitivity_covariate) && sensitivity_covariate %in% names(meta) &&
        !sensitivity_covariate %in% selected_covs) {
      md_cov <- meta[!is.na(meta[[sensitivity_covariate]]), , drop = FALSE]
      cat("Sensitivity cov:", sensitivity_covariate, "- campioni con dato:", nrow(md_cov), "di", nrow(meta), "\n")
      if (length(unique(md_cov[[sensitivity_covariate]])) > 1) {
        f_cov <- paste("~", col_time, "*", col_treatment, "+", sensitivity_covariate, "+", covs_txt, "+ (1 |", col_subject, ")")
        sens_list[["cov"]] <- run_sens("cov", f_cov, min_prevalence_primary, md_cov, later_timepoints)
      }
    }

    # (c) senza confronto con la mediana (riuso del fit principale)
    sens_list[["nomedian"]] <- tryCatch({
      ct_nm <- run_contrasts(fit_arm, S_arm[prim_id], median_abund = FALSE, analysis = "sens_nomedian")
      ct_nm$sensitivity <- "nomedian"
      ct_nm
    }, error = function(e) { message("[Sensitivity nomedian] failed: ", conditionMessage(e)); NULL })

    # (d) solo baseline + timepoint primario
    md_t2 <- meta[meta[[col_time]] %in% c(levels(meta[[col_time]])[1], primary_timepoint), , drop = FALSE]
    md_t2[[col_time]] <- droplevels(md_t2[[col_time]])
    sens_list[["T0_prim"]] <- run_sens("T0_prim", formula_arm, min_prevalence_primary, md_t2, primary_timepoint)

    sens_list <- Filter(Negate(is.null), sens_list)
    if (length(sens_list) > 0) {
      sens_all <- do.call(rbind, sens_list)
      rownames(sens_all) <- NULL
      write.table(sens_all, file.path(out_dir, "sensitivity_primary_contrast.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Sensibilita' completata. Varianti riuscite:", paste(names(sens_list), collapse = ", "), "\n")
    }
  }
}


# =============================================================================
# ---- 4d. LIVELLI DI EVIDENZA A/B/C PER EC E MODELLO --------------------------
# =============================================================================
# Contrasto primario del modello principale, per feature e modello (linear/logistic):
#   A : q_family <= alpha_primary e, in tutte le varianti di sensibilita' stimabili
#       (almeno una), stesso segno e q_family <= alpha_support
#   B : q_family <= alpha_primary ma non confermato in tutte le varianti
#   C : alpha_primary < q_family <= alpha_support (esplorativo)
#   NS: altrimenti (incluse stime instabili)
assign_evidence_tiers <- function(main_ct, sens_cts, test_id) {
  key <- function(d) paste(d$feature, d$model, sep = "|")
  m <- main_ct[main_ct$test == test_id, c("feature", "feature_name", "model", "coef", "stderr", "q_family"), drop = FALSE]
  names(m)[names(m) == "q_family"] <- "q_main"
  n_test <- integer(nrow(m)); n_conf <- integer(nrow(m))
  for (s in sens_cts) {
    s <- s[s$test == test_id, , drop = FALSE]
    i <- match(key(m), key(s))
    q <- s$q_family[i]; cf <- s$coef[i]
    testable <- !is.na(q)
    conf <- testable & q <= alpha_support & sign(cf) == sign(m$coef)
    conf[is.na(conf)] <- FALSE
    n_test <- n_test + testable
    n_conf <- n_conf + conf
  }
  m$n_sens_testable  <- n_test
  m$n_sens_confirmed <- n_conf
  sig <- !is.na(m$q_main) & m$q_main <= alpha_primary
  m$tier <- ifelse(sig & n_test >= 1 & n_conf == n_test, "A",
            ifelse(sig, "B",
            ifelse(!is.na(m$q_main) & m$q_main <= alpha_support, "C", "NS")))
  m[order(m$tier, m$q_main), ]
}

if (run_evidence_tiers) {
  if (exists("ct_arm") && length(sens_list) > 0) {
    cat("\n>>> [4d] Livelli di evidenza A/B/C...\n")
    tiers <- assign_evidence_tiers(ct_arm, sens_list, paste0("ix_", primary_timepoint))
    write.table(tiers, file.path(out_dir, "evidence_tiers_primary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    print(table(tiers$model, tiers$tier))
  } else {
    cat("[SKIP Evidence tiers] Richiedono modello principale e almeno una sensibilita' riuscita.\n")
  }
}


# =============================================================================
# ---- 5. REQ 6: CONTRASTO CONTINUO BILANCIO TRA DUE COMPONENTI (DOSE-BALANCE) -
# =============================================================================
# Valuta se la composizione della dieta/trattamento (Componente A vs B) ha
# un andamento temporale divergente a parità di dose totale assunta.
if (run_req6) {
  cat("\n>>> [5/7] Esecuzione Req 6 (Contrasto Continuo tra Componenti nei Trattati)...\n")
  
  # Isolamento sottopopolazione (es. solo i trattati)
  sub_meta6 <- if (!is.null(req6_subset_condition)) {
    meta[meta[[req6_subset_condition$column]] == req6_subset_condition$value, , drop = FALSE]
  } else {
    meta
  }
  
  if (req6_component_A %in% names(sub_meta6) && req6_component_B %in% names(sub_meta6)) {
    sub_meta6$Component_Total   <- as.numeric(sub_meta6[[req6_component_A]]) + as.numeric(sub_meta6[[req6_component_B]])
    sub_meta6$Component_Balance <- as.numeric(sub_meta6[[req6_component_A]]) - as.numeric(sub_meta6[[req6_component_B]])
    
    f6 <- paste("~", col_time, "* Component_Balance + Component_Total +", covs_txt, "+ (1 |", col_subject, ")")
    cat("Formula Req 6:", f6, "\n")
    
    fit6 <- run_fit(f6, "req6_dose_balance", min_prevalence_primary, md = sub_meta6)
    S6 <- build_generic_spec(coef_names(fit6), "Component_Balance", later_timepoints)
    ct6 <- run_contrasts(fit6, S6, analysis = "req6")
    
    export_tsv_and_plots(ct6, ec_raw, sub_meta6, file.path(out_dir, "req6_contrast_output"), "Req6_Balance")
    cat("Req 6 completata con successo.\n")
  } else {
    cat("[SKIP Req 6] Variabili dei componenti A/B non trovate nei metadati.\n")
  }
}

# =============================================================================
# ---- 6. REQ 7: DECONVOLUZIONE MULTIFATTORIALE (INTERVENTI MULTIPLI AGGIUSTATI) -
# =============================================================================
# DESIGN NOTE. In this study every treated subject receives ALL four components (bread,
# artichokes, cognitive and physical training) at a Low or High level, and controls receive
# none. The 8 dummies Pane_dic x {Low, High}, ... therefore have rank 5, not 8: "Low vs No" of
# each component cannot be separated from the overall treated-vs-control effect. The identifiable
# parametrisation is:
#     Treatment (all components, at Low level)  +  "High vs Low" extra effect of each component.
# Req 7 below uses it (interactions with time). Collinearity among the High indicators is
# reported as VIF (values above ~5 would mean the components cannot be separated). With only
# 39 treated subjects, read the component-specific terms as exploratory.
if (run_req7) {
  cat("\n>>> [6/7] Req 7 (multifactorial: Treatment + High-vs-Low gradient per component)...\n")
  meta7 <- meta
  hi_vars <- character(0)
  for (iv in intersect(names(req7_interventions), colnames(meta7))) {
    ref <- req7_interventions[[iv]]$ref
    lv  <- setdiff(unique(as.character(meta7[[iv]])), c(ref, NA, "None", "NA", ""))
    high <- lv[grepl("High", lv, ignore.case = TRUE)]
    if (length(high) != 1) { cat("[Req 7] No unique 'High' level for", iv, "- skipped\n"); next }
    nm <- paste0(iv, "_isHigh")
    meta7[[nm]] <- factor(ifelse(as.character(meta7[[iv]]) == high, "yes", "no"), levels = c("no", "yes"))
    hi_vars <- c(hi_vars, nm)
  }

  if (length(hi_vars) >= 2) {
    # collinearity among the High indicators (treated subjects, one row per subject)
    sub7 <- meta7[meta7[[col_treatment]] == treatment_act & !duplicated(meta7[[col_subject]]), hi_vars, drop = FALSE]
    X7 <- sapply(sub7, function(v) as.numeric(v == "yes"))
    vif <- sapply(seq_along(hi_vars), function(j) {
      r2 <- summary(lm(X7[, j] ~ X7[, -j, drop = FALSE]))$r.squared; 1 / (1 - r2) })
    names(vif) <- hi_vars
    cat("VIF of the High indicators among treated subjects (>5 = strong collinearity):\n"); print(round(vif, 2))
    write.csv2(data.frame(variable = hi_vars, VIF = vif), file.path(out_dir, "req7_VIF.csv"), row.names = FALSE)

    terms7 <- c(paste0(col_time, " * ", col_treatment), paste0(col_time, " * ", hi_vars))
    f7 <- paste("~", paste(terms7, collapse = " + "), "+", covs_txt, "+ (1 |", col_subject, ")")
    cat("Formula Req 7:", f7, "\n")
    fit7 <- run_fit(f7, "req7_multifactorial", min_prevalence_primary, md = meta7)
    cn7  <- coef_names(fit7)

    pick <- function(a, b) if (a %in% cn7) a else if (b %in% cn7) b else NULL
    S7 <- list()
    for (tk in later_timepoints) {
      tr <- paste0(col_treatment, treatment_act)
      hit <- pick(paste0(col_time, tk, ":", tr), paste0(tr, ":", col_time, tk))
      if (!is.null(hit)) S7[[paste0("ix_", tk, "_Treatment_overall")]] <- setNames(1, hit)
      for (hv in hi_vars) {
        cb <- paste0(hv, "yes")
        hit <- pick(paste0(col_time, tk, ":", cb), paste0(cb, ":", col_time, tk))
        if (!is.null(hit)) S7[[paste0("ix_", tk, "_", hv, "_HighVsLow")]] <- setNames(1, hit)
      }
    }
    dropped <- setdiff(c(paste0(col_treatment, treatment_act), paste0(hi_vars, "yes")), cn7)
    if (length(dropped) > 0) cat("[WARNING Req 7] coefficients dropped by the model (aliased):", paste(dropped, collapse = ", "), "\n")

    if (length(S7) > 0) {
      ct7 <- run_contrasts(fit7, S7, analysis = "req7")
      export_tsv_and_plots(ct7, ec_raw, meta7, file.path(out_dir, "req7_multifactorial_output"), "Req7_Multifactorial")
      cat("Req 7 completed.\n")
    } else {
      cat("[WARNING Req 7] No estimable contrast extracted.\n")
    }
  } else {
    cat("[SKIP Req 7] Fewer than two intervention components available.\n")
  }
}

# =============================================================================
# ---- 7. REQ 9: LOOP DI SENSIBILITÀ SU COVARIATE EXTRA (ATC / FARMACI) --------
# =============================================================================
if (run_req9_loop) {
  cat("\n>>> [7/7] Esecuzione Req 9 (Ciclo di robustezza su covariate extra)...\n")
  candidate_cols <- grep(req9_covar_pattern, colnames(meta), value = TRUE)
  tested_cols <- candidate_cols[vapply(candidate_cols, function(cn) {
    v <- meta[[cn]]; v <- v[!is.na(v)]
    length(unique(v)) > 1 && sum(v != 0) >= req9_min_pos_cases
  }, logical(1))]
  
  cat("Covariate idonee identificate per il ciclo:", length(tested_cols), "\n")
  loop_results <- list()
  for (cv in tested_cols) {
    res_cv <- tryCatch({
      f_cv <- paste("~", col_time, "*", col_treatment, "+", cv, "+", covs_txt, "+ (1 |", col_subject, ")")
      fit_cv <- run_fit(f_cv, paste0("loop_", make.names(cv)), min_prevalence_primary)

      exp_term_name <- paste0(col_treatment, treatment_act)
      S_cv <- build_generic_spec(coef_names(fit_cv), exp_term_name, later_timepoints)

      # primary interaction contrast only
      prim_contrast_id <- paste0("ix_", primary_timepoint)
      if (prim_contrast_id %in% names(S_cv)) {
        ct_cv <- run_contrasts(fit_cv, S_cv[prim_contrast_id], analysis = "req9")
        ct_cv$extra_covariate <- cv
        ct_cv
      } else NULL
    }, error = function(e) { message("[Req 9] covariate ", cv, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(res_cv)) loop_results[[cv]] <- res_cv
    saveRDS(loop_results, file.path(out_dir, "req9_progress.rds"))   # progress is kept if the run is interrupted
  }

  if (length(loop_results) > 0) {
    all_loop_df <- do.call(rbind, loop_results)
    write.table(all_loop_df, file.path(out_dir, "req9_extra_covariates_primary_summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\n=============================================================================\n")
cat("[COMPLETATO] Analisi completata con successo! Risultati in:", out_dir, "\n")
cat("=============================================================================\n")
