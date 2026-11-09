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
base_covariates <- c("Age", "Sex", "Bristol_T0", "City", "PPI")

# Covariate rilevate a ogni visita che possono essere influenzate dall'intervento (es. Bristol) NON vanno
# usate col valore corrente nel modello principale: si usa il valore BASALE (primo timepoint disponibile
# per soggetto), salvato in una nuova colonna <nome>_T0. Il valore corrente entra solo nella sensibilita'
# "bristol_now". Le covariate note solo al basale (es. PPI) sono sostituite dal loro valore basale.
baseline_vars    <- c("Bristol")   # crea Bristol_T0
baseline_replace <- c("PPI")       # sostituisce PPI col valore basale
bristol_current  <- "Bristol"      # colonna per la sensibilita' "Bristol corrente"
batch_col        <- "DNABatch"     # solo diagnostica: verifica se il batch e' confuso con tempo/braccio

# Aggiungere il log10 della profondità come covariata di aggiustamento (raccomandato da MaAsLin 3)?
include_depth_covariate <- TRUE

# Covariata opzionale aggiuntiva per analisi di sensibilità
sensitivity_covariate <- "burden_now"

# Imputazione LOCF (Last Observation Carried Forward) per colonne con pattern specifico (es. farmaci)
# Impostare su NULL se non si desidera effettuare LOCF
locf_col_pattern <- "_L4$"

# --- 0.4 CONFIGURAZIONE RICHIESTE CLINICHE (SWITCH ON/OFF) ---
run_main_models <- TRUE   # Esegue Req 1 (basale) e Req 5 (Trial principale); il Dose-model ha il proprio switch (run_dose_model)
run_req6        <- TRUE   # Req 6: Contrasto continuo o bilancio tra 2 interventi nei trattati
run_req7        <- TRUE   # Req 7: Deconvoluzione multifattoriale (più interventi simultanei)
run_req9_loop   <- FALSE # Req 9: Ciclo di robustezza su covariate extra (es. classi di farmaci)
run_dose_model     <- TRUE   # Dose-response model (Time * composite adherence dose)
run_sensitivity    <- TRUE   # Sensitivity analyses of the primary contrast (prevalence, covariate, no median, T0-T2 only)
run_evidence_tiers <- TRUE   # Evidence tiers A/B/C per EC and model

# Optional command-line override: Rscript script.R main|dose|req6|req7|req9
.run_sel <- commandArgs(trailingOnly = TRUE)
if (length(.run_sel) > 0) {
  run_main_models <- "main" %in% .run_sel
  run_sensitivity <- run_evidence_tiers <- "main" %in% .run_sel
  run_dose_model  <- "dose" %in% .run_sel
  run_req6        <- "req6" %in% .run_sel
  run_req7        <- "req7" %in% .run_sel
  run_req9_loop   <- "req9" %in% .run_sel
}

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

# Escludere (flag 'unstable') anche le stime del ramo LINEARE con warning/errore di fitting lme4
# (es. mancata convergenza, Hessiana degenere)? I 'boundary (singular) fit' NON vengono flaggati.
flag_linear_convergence <- TRUE

# MaAsLin 3 segnala "<4 average observations per random effect group often inflates coefficients..."
# quando ci sono in media meno di 4 osservazioni per soggetto (sempre vero nella variante T0_prim,
# che ha solo 2 tempi). E' un avviso di cautela, non un fallimento del fit.
# TRUE  = non conta come stima instabile (resta registrato nella colonna warn_few_obs);
#         le stime comunque inaffidabili sono ancora escluse dalle soglie su SE/coef e da "Fit is NA",
#         "large coefficient", mancata convergenza.
# FALSE = qualsiasi messaggio in 'error' esclude la stima (molto conservativo).
ignore_few_obs_warning <- TRUE

# Escludere gli EC parziali (con trattino, es. "1.1.1.-"): non sono funzioni enzimatiche specifiche.
# Decisione a priori, indipendente dai risultati; l'elenco degli esclusi e' salvato in out_dir.
# Si tengono i codici completi a 4 livelli, compresi quelli provvisori (es. 3.4.21.n1).
exclude_partial_ec <- TRUE

# Parametri operativi / debug
max_features_test <- NULL    # Impostare ad un intero (es. 50) per un test run rapido, altrimenti NULL
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

if (exclude_partial_ec) {
  is_full <- grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[A-Za-z]?[0-9]+$", colnames(ec_raw))
  if (sum(is_full) == 0) stop("Nessun EC completo trovato: controllare il formato dei nomi delle feature.")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  excl <- colnames(ec_raw)[!is_full]
  cat("EC totali:", ncol(ec_raw), "| completi tenuti:", sum(is_full), "| esclusi (parziali/altro):", length(excl), "\n")
  if (length(excl) > 0) {
    cat("Esempi di feature escluse:", paste(head(excl, 6), collapse = ", "), "\n")
    write.table(data.frame(feature = excl), file.path(out_dir, "excluded_features_partial_EC.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  ec_raw <- ec_raw[, is_full, drop = FALSE]
}

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

# ---- Covariate al basale ------------------------------------------------------------------
baseline_value <- function(v) {
  ord <- order(meta[[col_subject]], match(as.character(meta[[col_time]]), time_levels))
  first <- tapply(meta[[v]][ord], meta[[col_subject]][ord], function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else x[NA_integer_][1] })
  unname(first[as.character(meta[[col_subject]])])
}
for (v in baseline_vars) {
  if (v %in% names(meta)) {
    meta[[paste0(v, "_T0")]] <- baseline_value(v)
    cat("Covariata basale creata:", paste0(v, "_T0"), "| soggetti senza dato:",
        sum(is.na(tapply(meta[[paste0(v, "_T0")]], meta[[col_subject]], function(x) x[1]))), "\n")
  }
}
for (v in baseline_replace) {
  if (v %in% names(meta)) {
    new_v <- baseline_value(v)
    n_chg <- length(unique(meta[[col_subject]][!is.na(meta[[v]]) & !is.na(new_v) & meta[[v]] != new_v]))
    cat("Covariata", v, "sostituita col valore basale | soggetti con valori diversi nel tempo:", n_chg, "\n")
    meta[[v]] <- new_v
  }
}

# Diagnostica batch: se coincide con il tempo non e' separabile da esso (e non puo' essere aggiustato)
if (!is.null(batch_col) && batch_col %in% names(meta)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  bt <- as.data.frame(table(Treatment = meta[[col_treatment]], Time = meta[[col_time]], Batch = meta[[batch_col]]),
                      responseName = "N_samples")
  bt <- bt[bt$N_samples > 0, ]
  cat("\nBatch (", batch_col, ") per braccio x tempo:\n", sep = ""); print(bt, row.names = FALSE)
  write.table(bt, file.path(out_dir, "batch_by_arm_time.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  tb <- tapply(as.character(meta[[batch_col]]), as.character(meta[[col_time]]), function(x) length(unique(x)))
  if (all(tb == 1) && length(unique(as.character(meta[[batch_col]]))) > 1) {
    cat("[ATTENZIONE] Il batch", batch_col, "e' interamente determinato dal tempo: non separabile dall'effetto tempo.\n",
        "I contrasti fra bracci restano interpretabili (stesso batch nello stesso tempo); le variazioni entro braccio no.\n")
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

# Livelli di riferimento delle variabili categoriche del modello (il PRIMO livello e' il riferimento:
# gli altri livelli sono confrontati con quello). Tempo e braccio hanno ordine esplicito (factor);
# le altre covariate categoriche seguono l'ordine alfabetico, salvo che siano gia' factor.
# Gli altri confronti (diff_Tk, wash_*, High vs Low) NON richiedono di cambiare riferimento ne' di
# rifare il modello: sono combinazioni lineari dei coefficienti (maaslin_contrast_test).
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ref_tab <- do.call(rbind, lapply(unique(c(col_time, col_treatment, selected_covs)), function(v) {
  x <- meta[[v]]
  if (is.null(x) || is.numeric(x)) return(NULL)
  lv <- if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
  data.frame(variable = v, reference_level = lv[1], other_levels = paste(lv[-1], collapse = " | "),
             n_levels = length(lv),
             order_source = if (is.factor(x)) "factor (ordine esplicito)" else "carattere (ordine alfabetico)")
}))
if (!is.null(ref_tab)) {
  cat("\nLivelli di riferimento delle variabili categoriche:\n"); print(ref_tab, row.names = FALSE)
  write.table(ref_tab, file.path(out_dir, "model_reference_levels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}


# ---- Sample size check per arm x time (diagnostic only) ----
min_n_cell <- 10
n_cell <- as.data.frame(table(Treatment = meta[[col_treatment]], Time = meta[[col_time]]),
                        responseName = "N_samples")
n_subj <- aggregate(meta[[col_subject]],
                    by = list(Treatment = meta[[col_treatment]], Time = meta[[col_time]]),
                    FUN = function(x) length(unique(x)))
names(n_subj)[3] <- "N_subjects"
n_cell <- merge(n_cell, n_subj, by = c("Treatment", "Time"))
cat("\nSample size per arm x time:\n"); print(n_cell, row.names = FALSE)
write.table(n_cell, file.path(out_dir, "sample_size_arm_time.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
low <- n_cell[n_cell$N_samples < min_n_cell, ]
if (nrow(low) > 0) {
  cat("[WARNING] Cells below", min_n_cell, "samples:\n"); print(low, row.names = FALSE)
}


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
  # NB: nei risultati di MaAsLin 3 la colonna 'model' vale "abundance" oppure "prevalence"
  # (NON "linear"/"logistic"): i confronti sotto usano questi valori.
  if (!all(ct$model %in% c("abundance", "prevalence"))) {
    warning("[add_flags] valori inattesi in ct$model: ", paste(unique(ct$model), collapse = ", "))
  }
  
  # conv_warning: la feature ha un messaggio di errore/warning di fitting (colonna 'error' di MaAsLin 3)
  if (!"error" %in% names(ct)) {
    message("[add_flags] colonna 'error' assente nei risultati: il flag di convergenza non e' applicabile.")
    err_txt <- rep(NA_character_, nrow(ct))
  } else {
    # i messaggi di lme4 possono contenere a-capo e tab: li sostituisco con spazi,
    # altrimenti le righe dei TSV esportati si spezzano e le colonne si sfasano
    err_txt <- gsub("[\r\n\t]+", " ", as.character(ct$error))
    ct$error <- err_txt
  }
  has_msg   <- !is.na(err_txt) & nzchar(trimws(err_txt))
  few_obs   <- has_msg & grepl("average observations per random effect group", err_txt, ignore.case = TRUE)
  hard_fail <- has_msg & grepl("Fit is NA|large coefficient|converge|Hessian|gradient", err_txt, ignore.case = TRUE)
  ct$warn_few_obs <- few_obs
  ct$conv_warning <- if (ignore_few_obs_warning) has_msg & !(few_obs & !hard_fail) else has_msg
  if (any(has_msg)) {
    cat(sprintf("[add_flags:%s] messaggi 'error' piu' frequenti (primi 60 caratteri):\n", analysis))
    print(head(sort(table(substr(err_txt[has_msg], 1, 60)), decreasing = TRUE), 5))
  }
  
  prev_hi <- if (all(c("N", "N_not_zero") %in% names(ct))) {
    (ct$N_not_zero / ct$N) > max_prev_logit
  } else rep(FALSE, nrow(ct))
  
  # ramo prevalence (logistico): errore di fitting, SE mancante/alto, coefficiente estremo, prevalenza ~100%
  ct$unstable <- ct$model == "prevalence" &
    (ct$conv_warning | is.na(ct$stderr) | ct$stderr > se_max_logit | abs(ct$coef) > coef_max_logit | prev_hi)
  # ramo abundance (lineare): escluso se c'e' un warning/errore di convergenza
  if (flag_linear_convergence) {
    ct$unstable <- ct$unstable | (ct$model == "abundance" & ct$conv_warning)
  }
  ct$unstable[is.na(ct$unstable)] <- TRUE
  cat(sprintf("[add_flags:%s] stime escluse come instabili: %d su %d (abundance: %d, prevalence: %d)\n",
              analysis, sum(ct$unstable), nrow(ct),
              sum(ct$unstable & ct$model == "abundance"), sum(ct$unstable & ct$model == "prevalence")))
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
  ct <- add_effect_sizes(ct)
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
      mod_label <- if (mm == "abundance") "Abundance (log2 FC)" else "Prevalence (log-odds)"
      
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
# ---- 3b. EFFECT SIZE, CONTEGGI PER GRUPPO E TABELLE FINALI ------------------
# =============================================================================
# Effect size aggiustato per le covariate del modello (il coefficiente del contrasto e' gia'
# aggiustato per Age, Sex, ecc.), con intervallo di confidenza al 95% (approssimazione normale).
#   abundance : fold change = 2^coef            (coef su scala log2)
#   prevalence: odds ratio  = exp(coef)         (coef = log-odds)
# Per i contrasti di interazione (ix_*, wash_*) il rapporto e' un rapporto di rapporti
# (variazione nel braccio attivo / variazione nel braccio di controllo).
# Per analisi con variabili continue standardizzate (dose, bilancio componenti) il coefficiente
# e il rapporto sono per +1 SD della variabile.
add_effect_sizes <- function(ct) {
  z <- qnorm(0.975)
  is_ab <- ct$model == "abundance"
  conv  <- function(x) ifelse(is_ab, 2^x, exp(x))
  ct$coef_CI95_low          <- ct$coef - z * ct$stderr
  ct$coef_CI95_high         <- ct$coef + z * ct$stderr
  ct$Effect_ratio_type      <- ifelse(is_ab, "fold_change (2^coef)", "odds_ratio (exp(coef))")
  ct$Effect_ratio           <- conv(ct$coef)
  ct$Effect_ratio_CI95_low  <- conv(ct$coef_CI95_low)
  ct$Effect_ratio_CI95_high <- conv(ct$coef_CI95_high)
  ct
}

# Conteggi grezzi per gruppo "test" e "ref" di ciascun contrasto del modello principale (Treatment x Time).
# Convenzione: il coefficiente e' sempre TEST meno REF.
#   base_diff   : test = attivo a T0           | ref = controllo a T0
#   diff_Tk     : test = attivo a Tk           | ref = controllo a Tk
#   ix_Tk       : test = variazione attivo T0->Tk | ref = variazione controllo T0->Tk
#   wash_Tk_vs_P: test = variazione attivo P->Tk  | ref = variazione controllo P->Tk  (P = timepoint primario)
# Per i contrasti di variazione (ix_*, wash_*) sono riportati anche i valori ai due tempi (colonne *_a e *_b).
# Colonne: N_* campioni, Subj_* soggetti, Pos_* campioni con conteggio > 0, Prev_*_pct prevalenza (%),
# Counts_* somma dei conteggi grezzi, Mean_* media dei conteggi grezzi per campione.
add_group_counts <- function(ct, raw_mat, sub_meta) {
  tm <- as.character(sub_meta[[col_time]])
  ar <- as.character(sub_meta[[col_treatment]])
  t0 <- levels(sub_meta[[col_time]])[1]
  ids_of  <- function(arm, times) rownames(sub_meta)[ar %in% arm & tm %in% times]
  n_subj  <- function(ids) length(unique(sub_meta[ids, col_subject]))
  st_of   <- function(ids, feats) {
    m <- raw_mat[ids, feats, drop = FALSE]
    list(pos = colSums(m > 0, na.rm = TRUE), tot = colSums(m, na.rm = TRUE),
         mean = if (length(ids) > 0) colMeans(m, na.rm = TRUE) else setNames(rep(NA_real_, length(feats)), feats))
  }
  parse_test <- function(tt) {
    if (tt == "base_diff") return(list(kind = "cross", time = t0))
    if (grepl("^diff_T[0-9]+$", tt)) return(list(kind = "cross", time = sub("^diff_", "", tt)))
    if (grepl("^ix_T[0-9]+$", tt)) return(list(kind = "change", a = t0, b = sub("^ix_", "", tt)))
    m <- regmatches(tt, regexec("^wash_(T[0-9]+)_vs_(T[0-9]+)$", tt))[[1]]
    if (length(m) == 3) return(list(kind = "change", a = m[3], b = m[2]))
    NULL
  }
  chr_cols <- c("Test_group", "Ref_group", "Time_a", "Time_b")
  num_cols <- c("N_test", "N_ref", "Subj_test", "Subj_ref", "Pos_test", "Pos_ref",
                "Prev_test_pct", "Prev_ref_pct", "Counts_test", "Counts_ref", "Mean_test", "Mean_ref",
                "Pos_test_a", "Pos_test_b", "Pos_ref_a", "Pos_ref_b",
                "Mean_test_a", "Mean_test_b", "Mean_ref_a", "Mean_ref_b")
  for (cn in chr_cols) ct[[cn]] <- NA_character_
  for (cn in num_cols) ct[[cn]] <- NA_real_
  
  for (tt in unique(ct$test)) {
    sp <- parse_test(tt)
    if (is.null(sp)) next
    rows <- which(ct$test == tt & !is.na(ct$feature_name) & ct$feature_name %in% colnames(raw_mat))
    if (length(rows) == 0) next
    f <- ct$feature_name[rows]
    feats <- unique(f)
    if (sp$kind == "cross") {
      ids_t <- ids_of(treatment_act, sp$time); ids_r <- ids_of(treatment_ref, sp$time)
      ct$Test_group[rows] <- paste0(treatment_act, "@", sp$time)
      ct$Ref_group[rows]  <- paste0(treatment_ref, "@", sp$time)
    } else {
      ids_t <- ids_of(treatment_act, c(sp$a, sp$b)); ids_r <- ids_of(treatment_ref, c(sp$a, sp$b))
      ct$Test_group[rows] <- paste0(treatment_act, ": ", sp$b, " vs ", sp$a)
      ct$Ref_group[rows]  <- paste0(treatment_ref, ": ", sp$b, " vs ", sp$a)
      ct$Time_a[rows] <- sp$a; ct$Time_b[rows] <- sp$b
      ta <- st_of(ids_of(treatment_act, sp$a), feats); tb <- st_of(ids_of(treatment_act, sp$b), feats)
      ra <- st_of(ids_of(treatment_ref, sp$a), feats); rb <- st_of(ids_of(treatment_ref, sp$b), feats)
      ct$Pos_test_a[rows]  <- ta$pos[f];  ct$Pos_test_b[rows]  <- tb$pos[f]
      ct$Pos_ref_a[rows]   <- ra$pos[f];  ct$Pos_ref_b[rows]   <- rb$pos[f]
      ct$Mean_test_a[rows] <- round(ta$mean[f], 3); ct$Mean_test_b[rows] <- round(tb$mean[f], 3)
      ct$Mean_ref_a[rows]  <- round(ra$mean[f], 3); ct$Mean_ref_b[rows]  <- round(rb$mean[f], 3)
    }
    st <- st_of(ids_t, feats); sr <- st_of(ids_r, feats)
    ct$N_test[rows] <- length(ids_t); ct$N_ref[rows] <- length(ids_r)
    ct$Subj_test[rows] <- n_subj(ids_t); ct$Subj_ref[rows] <- n_subj(ids_r)
    ct$Pos_test[rows] <- st$pos[f]; ct$Pos_ref[rows] <- sr$pos[f]
    ct$Prev_test_pct[rows] <- if (length(ids_t) > 0) round(100 * st$pos[f] / length(ids_t), 1) else NA
    ct$Prev_ref_pct[rows]  <- if (length(ids_r) > 0) round(100 * sr$pos[f] / length(ids_r), 1) else NA
    ct$Counts_test[rows] <- st$tot[f]; ct$Counts_ref[rows] <- sr$tot[f]
    ct$Mean_test[rows] <- round(st$mean[f], 3); ct$Mean_ref[rows] <- round(sr$mean[f], 3)
  }
  ct
}

# Tabella finale per un contrasto: solo stime stabili, ordinate per q-value, con effect size,
# conteggi per gruppo, conteggi totali dell'EC e (se disponibili) livelli di evidenza A/B/C.
build_publication_table <- function(ct, raw_mat, sub_meta, test_id, tiers = NULL) {
  d <- ct[ct$test == test_id & !ct$unstable, , drop = FALSE]
  tot <- colSums(raw_mat[rownames(sub_meta), , drop = FALSE], na.rm = TRUE)
  d$Total_Raw_Counts <- tot[d$feature_name]
  if (!is.null(tiers)) {
    k <- match(paste(d$feature, d$model), paste(tiers$feature, tiers$model))
    d$Evidence_tier     <- tiers$tier[k]
    d$n_sens_testable   <- tiers$n_sens_testable[k]
    d$n_sens_confirmed  <- tiers$n_sens_confirmed[k]
  }
  d <- d[order(d$q_family, d$pval_individual), , drop = FALSE]
  cols <- c("feature_name", "model", "Test_group", "Ref_group",
            "coef", "stderr", "coef_CI95_low", "coef_CI95_high",
            "Effect_ratio_type", "Effect_ratio", "Effect_ratio_CI95_low", "Effect_ratio_CI95_high",
            "pval_individual", "q_family", "Evidence_tier", "n_sens_testable", "n_sens_confirmed",
            "Total_Raw_Counts", "N_test", "N_ref", "Subj_test", "Subj_ref",
            "Pos_test", "Pos_ref", "Prev_test_pct", "Prev_ref_pct",
            "Counts_test", "Counts_ref", "Mean_test", "Mean_ref",
            "Time_a", "Time_b", "Pos_test_a", "Pos_test_b", "Pos_ref_a", "Pos_ref_b",
            "Mean_test_a", "Mean_test_b", "Mean_ref_a", "Mean_ref_b",
            "warn_few_obs", "inference")
  d <- d[, intersect(cols, names(d)), drop = FALSE]
  names(d)[names(d) == "feature_name"] <- "EC"
  rownames(d) <- NULL
  d
}

# =============================================================================
# ---- 4. MODELLO PRINCIPALE (TRIAL LONGITUDINALE: REQ 1 & REQ 5) -------------
# =============================================================================
if (run_main_models) {
  cat("\n>>> [4/7] Esecuzione Modello Principale (Req 1 Basale & Req 5 Interazione)...\n")
  fit_arm <- run_fit(formula_arm, "main_arm", min_prevalence_primary)
  
  exp_term_name <- paste0(col_treatment, treatment_act)
  S_arm <- build_generic_spec(coef_names(fit_arm), exp_term_name, later_timepoints)
  cat("Coefficienti del modello principale:", paste(coef_names(fit_arm), collapse = ", "), "\n")
  ct_arm <- run_contrasts(fit_arm, S_arm, analysis = "main")
  ct_arm <- add_group_counts(ct_arm, ec_raw, meta)   # conteggi per gruppo test/ref di ogni contrasto
  
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
      
      rm(fit_dose, S_dose); invisible(gc())
      
      cat("Dose-response completato.\n")
    }, error = function(e) message("[Dose] failed: ", conditionMessage(e)))
  } else {
    cat("[SKIP Dose] Colonna 'dose' non disponibile (servono Pane e Carciofini).\n")
  }
}


# =============================================================================
# ---- 4c. ANALISI DI SENSIBILITA' DEL CONTRASTO PRIMARIO ----------------------
# =============================================================================
# Ricalcola solo il contrasto di interazione al timepoint primario in 5 varianti:
#   prev40   : prevalenza minima piu' stringente (min_prevalence_sens)
#   cov      : aggiustamento per sensitivity_covariate (casi completi)
#   nomedian : senza confronto con la mediana dei coefficienti (abbondanza)
#   T0_prim  : solo baseline e timepoint primario
#   bristol_now : Bristol corrente (ogni visita) al posto del valore basale
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
    
    # (c) senza confronto con la mediana (riuso del fit principale)
    sens_list[["nomedian"]] <- tryCatch({
      ct_nm <- run_contrasts(fit_arm, S_arm[prim_id], median_abund = FALSE, analysis = "sens_nomedian")
      ct_nm$sensitivity <- "nomedian"
      ct_nm
    }, error = function(e) { message("[Sensitivity nomedian] failed: ", conditionMessage(e)); NULL })
    
    rm(fit_arm, S_arm); invisible(gc())
    
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
    
    # (e) Bristol CORRENTE (valore di ogni visita) al posto del valore basale
    bristol_base <- paste0(bristol_current, "_T0")
    if (bristol_current %in% names(meta) && bristol_base %in% selected_covs) {
      covs_now <- paste(replace(selected_covs, selected_covs == bristol_base, bristol_current), collapse = " + ")
      f_now <- paste("~", col_time, "*", col_treatment, "+", covs_now, "+ (1 |", col_subject, ")")
      cat("Sensitivity bristol_now - formula:", f_now, "\n")
      sens_list[["bristol_now"]] <- run_sens("bristol_now", f_now, min_prevalence_primary, meta, later_timepoints)
    }
    
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
# Contrasto primario del modello principale, per feature e modello (abundance/prevalence):
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

if (exists("fit_arm")) rm(fit_arm, S_arm); invisible(gc())

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
    
    rm(fit6, S6); invisible(gc())
    
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
  
  if (length(hi_vars) > 0) {
    # ---- Req 7 group size check (diagnostic only) ----
    min_n_level <- 10
    d7 <- meta7[meta7[[col_treatment]] == treatment_act, , drop = FALSE]
    n_hi <- do.call(rbind, lapply(hi_vars, function(hv) {
      g <- as.data.frame(table(Time = d7[[col_time]], Level = d7[[hv]]), responseName = "N_samples")
      cbind(variable = hv, g)
    }))
    cat("\n[Req 7] Treated samples per component level x time (no = Low, yes = High):\n")
    print(n_hi, row.names = FALSE)
    write.table(n_hi, file.path(out_dir, "req7_group_sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    low7 <- n_hi[n_hi$N_samples < min_n_level, ]
    if (nrow(low7) > 0) { cat("[WARNING Req 7] Cells below", min_n_level, "samples:\n"); print(low7, row.names = FALSE) }
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
      
      rm(fit7); invisible(gc())
      
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
    ok <- !is.na(meta[[cn]])
    length(unique(meta[[cn]][ok])) > 1 &&
      length(unique(meta[[col_subject]][ok & meta[[cn]] != 0])) >= req9_min_pos_cases
  }, logical(1))]
  
  cat("Covariate idonee identificate per il ciclo:", length(tested_cols), "\n")
  cat("Tested covariates:", paste(tested_cols, collapse = ", "), "\n")
  
  loop_results <- list()
  exp_term_name <- paste0(col_treatment, treatment_act)
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
    
    if (exists("fit_cv")) rm(fit_cv); invisible(gc())
    
  }
  
  if (length(loop_results) > 0) {
    all_loop_df <- do.call(rbind, loop_results)
    write.table(all_loop_df, file.path(out_dir, "req9_extra_covariates_primary_summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# =============================================================================
# ---- 8. TABELLE FINALI PER PUBBLICAZIONE -------------------------------------
# =============================================================================
if (exists("ct_arm")) {
  cat("\n>>> [8] Tabelle finali per pubblicazione...\n")
  pub_id    <- paste0("ix_", primary_timepoint)
  tiers_obj <- if (run_evidence_tiers && exists("tiers")) tiers else NULL
  pub <- build_publication_table(ct_arm, ec_raw, meta, pub_id, tiers_obj)
  write.table(pub, file.path(out_dir, paste0("PUBLICATION_primary_", pub_id, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pub_sig <- pub[!is.na(pub$q_family) & pub$q_family <= alpha_primary, , drop = FALSE]
  cat("Contrasto primario", pub_id, ": stime stabili =", nrow(pub),
      "| significative (q_family <=", alpha_primary, ") =", nrow(pub_sig), "\n")
  if (nrow(pub_sig) > 0) {
    write.table(pub_sig, file.path(out_dir, paste0("PUBLICATION_primary_", pub_id, "_significant.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # tutti i contrasti del trial in un unico file (stime stabili), con conteggi per gruppo
  all_main <- ct_arm[!ct_arm$unstable, , drop = FALSE]
  all_main <- all_main[order(all_main$test, all_main$q_family), , drop = FALSE]
  write.table(all_main, file.path(out_dir, "PUBLICATION_all_contrasts_stable.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# ---- Cleanup: remove large model files (contrast results are already saved as TSV/RDS) ----
drop_model_files <- function(min_mb = 50, keep = "main_arm") {
  dirs <- list.dirs(out_dir, recursive = FALSE, full.names = FALSE)
  run_dirs <- dirs[dirs %in% c("main_arm", "dose_response", "req6_dose_balance", "req7_multifactorial") |
                     grepl("^(sens_|loop_)", dirs)]
  for (rd in setdiff(run_dirs, keep)) {
    f <- list.files(file.path(out_dir, rd), pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
    big <- f[file.size(f) > min_mb * 1024^2]
    if (length(big) > 0) {
      cat("[cleanup]", rd, ": removing", length(big), "model file(s),",
          round(sum(file.size(big)) / 1024^3, 2), "GB\n")
      unlink(big)
    }
  }
}
drop_model_files()

cat("\n=============================================================================\n")
cat("[COMPLETATO] Analisi completata con successo! Risultati in:", out_dir, "\n")
cat("=============================================================================\n")