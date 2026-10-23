#!/usr/bin/env Rscript
# =====================================================================================================
# BayesPrism benchmark with single-cell ground truth (pseudo-bulk), colorectal cancer.
#
# AGENT NOTE: self-contained, written for Linux. Run `Rscript bayesprism_benchmark.R --help` for options.
# Outputs are listed below; report.md summarizes the run.
#
# Linux setup: R >= 4.2 plus build libraries (Debian/Ubuntu:
#   sudo apt-get install -y r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev)
#   Rscript bayesprism_benchmark.R --install=TRUE --help   # installs missing packages into R_LIBS_USER, then exits
# or manually in R:
#   install.packages(c("BiocManager","remotes","data.table","R.utils","Matrix","snowfall","NMF","gplots"))
#   BiocManager::install(c("scran","BiocParallel","Biobase")); remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")
# Full run, detached from the terminal (all cores minus one):
#   nohup Rscript bayesprism_benchmark.R --out=bp_benchmark > bp_benchmark.log 2>&1 &
# Quick check:
#   Rscript bayesprism_benchmark.R --ref_cap=60 --test_cells=200 --n_mix=6 --mix_cells=400
#
# What it does
#   1. Downloads Lee et al. 2020 colorectal 10x scRNA-seq (GEO GSE132465: raw UMI counts + cell annotation)
#      unless already present in --data_dir.
#   2. Labels cells. cell type: Tumor (epithelial cells from tumor samples), Normal_epithelium (epithelial
#      cells from normal samples), T cells, B cells, Myeloids, Stromal cells, Mast cells. cell state: the
#      authors' Cell_subtype within each type; states with < --min_cells cells are merged into <type>|other.
#   3. Splits PATIENTS (not cells) into reference and held-out test sets (--test_frac), so no cell or
#      patient contributes to both the reference and the pseudo-bulk being deconvolved.
#   4. Builds pseudo-bulk samples from held-out cells, with known truth:
#        natural : one per held-out sample (tumor or normal), cells subsampled (--test_cells), so the real
#                  composition of that sample is kept
#        random  : --n_mix mixtures; for each, one held-out patient and Dirichlet(1) cell-type proportions,
#                  --mix_cells cells drawn from that patient (wider range of fractions)
#      Truth per pseudo-bulk: UMI fraction and cell fraction per cell type and state, and the summed counts
#      of each cell type (true cell-type-specific expression).
#   5. Runs BayesPrism (reference from reference patients, --ref_cap cells per state, key = "Tumor",
#      gene filtering as recommended by the authors; --genes pc = all protein-coding, markers = marker genes).
#   6. Scores
#        fractions  : per cell type and state, Pearson, Spearman, RMSE and bias of BayesPrism theta versus
#                     the true UMI fraction (theta estimates the fraction of reads) and versus cell fraction
#        expression : per cell type, BayesPrism Z versus true counts: per-gene Pearson across pseudo-bulks
#                     (log1p), per-sample Pearson across genes (log1p), normalized absolute error
#                     sum|Z - truth| / sum(truth), and the share of marker-gene reads assigned to Tumor
#                     (estimated vs true, a leakage check)
#
# Outputs (--out)
#   truth_fractions_type.tsv, truth_fractions_state.tsv, est_fractions_type.tsv, est_fractions_state.tsv,
#   metrics_fractions.tsv, metrics_expression.tsv, leakage.tsv, fractions_scatter.png, report.md, run_info.txt
#
# Caveat: reference and pseudo-bulk come from the same 10x platform, so this measures the best case. Real
# bulk RNA-seq (e.g. TCGA) adds platform differences that BayesPrism can only partly correct.
# =====================================================================================================

defaults <- list(data_dir = "bp_bench_data", out = "bp_benchmark", n_cores = "auto", ref_cap = 200, test_cells = 600,
                 n_mix = 30, mix_cells = 1000, genes = "pc", seed = 1, test_frac = 1 / 3, min_cells = 20, install = "FALSE")
args <- commandArgs(trailingOnly = TRUE)
show_help <- any(args %in% c("-h", "--help"))
args <- setdiff(args, c("-h", "--help"))
opt <- defaults
for (a in args) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2 || !kv[1] %in% names(opt)) stop("Unknown or malformed argument: ", a)
  opt[[kv[1]]] <- if (is.numeric(defaults[[kv[1]]])) as.numeric(kv[2]) else kv[2]
}
if (show_help) {
  cat("Usage: Rscript bayesprism_benchmark.R [--key=value ...]\nOptions (defaults):\n")
  for (k in names(defaults)) cat(sprintf("  --%-11s %s\n", k, defaults[[k]]))
  cat("  n_cores=auto uses all cores minus one; install=TRUE installs missing R packages into R_LIBS_USER\n")
  cat("Quick check: Rscript bayesprism_benchmark.R --ref_cap=60 --test_cells=200 --n_mix=6 --mix_cells=400\n")
}

# ---------------------------------------------------------------- packages (optional install)
lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }
if (isTRUE(as.logical(opt$install))) {
  repos <- "https://cloud.r-project.org"
  for (p in c("BiocManager", "remotes", "data.table", "R.utils", "Matrix", "snowfall", "NMF", "gplots"))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib, repos = repos)
  BiocManager::install(c("scran", "BiocParallel", "Biobase"), lib = lib, ask = FALSE, update = FALSE)
  if (!requireNamespace("BayesPrism", quietly = TRUE))
    remotes::install_github("Danko-Lab/BayesPrism/BayesPrism", lib = lib, upgrade = "never")
}
if (show_help) quit(status = 0)
req <- c("BayesPrism", "Matrix", "data.table", "R.utils")
missing_pkgs <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), ". Re-run with --install=TRUE or see the header.")
suppressPackageStartupMessages({ library(BayesPrism); library(Matrix); library(data.table) })
opt$n_cores <- if (identical(opt$n_cores, "auto")) max(1L, parallel::detectCores() - 1L) else as.integer(opt$n_cores)
options(timeout = 3600)
t0 <- Sys.time()
dir.create(opt$data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---------------------------------------------------------------- 1. data
geo <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/"
f_ann <- file.path(opt$data_dir, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
f_mat <- file.path(opt$data_dir, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
for (f in c(f_ann, f_mat)) if (!file.exists(f)) { msg("downloading ", basename(f)); download.file(paste0(geo, basename(f)), f, mode = "wb") }

# ---------------------------------------------------------------- 2. labels
ann <- fread(f_ann)
ann[, cell_type := Cell_type]
ann[Cell_type == "Epithelial cells" & Class == "Tumor", cell_type := "Tumor"]
ann[Cell_type == "Epithelial cells" & Class == "Normal", cell_type := "Normal_epithelium"]
ann[, cell_state := paste0(gsub(" ", "_", cell_type), "|", gsub(" ", "_", Cell_subtype))]
small <- ann[, .N, by = cell_state][N < opt$min_cells, cell_state]
ann[cell_state %in% small, cell_state := paste0(gsub(" ", "_", cell_type), "|other")]

# ---------------------------------------------------------------- 3. patient split + cell selection
set.seed(opt$seed)
pts <- sort(unique(ann$Patient))
test_pts <- sort(sample(pts, max(1, round(length(pts) * opt$test_frac))))
ref <- ann[!Patient %in% test_pts][, .SD[sample(.N, min(.N, opt$ref_cap))], by = cell_state]
tst <- ann[Patient %in% test_pts][, .SD[sample(.N, min(.N, opt$test_cells))], by = Sample]
msg(length(pts), " patients: ", length(pts) - length(test_pts), " reference, ", length(test_pts), " held out (",
    paste(test_pts, collapse = ","), "); reference cells ", nrow(ref), ", test cells ", nrow(tst))

# ---------------------------------------------------------------- read only the selected cells
hdr <- names(fread(f_mat, nrows = 0))
need <- unique(c(ref$Index, tst$Index))
pos <- match(need, hdr)
if (anyNA(pos)) stop(sum(is.na(pos)), " selected cells not found in the count matrix header")
msg("reading ", length(need), " cells from the UMI matrix (this parses the whole file; a few minutes)")
dt <- fread(f_mat, select = c(1L, pos))
genes <- dt[[1]]
blocks <- lapply(split(seq_along(need) + 1L, ceiling(seq_along(need) / 1000)),
                 function(ix) Matrix(as.matrix(dt[, ix, with = FALSE]), sparse = TRUE))
M <- do.call(cbind, blocks); rm(dt, blocks); invisible(gc())
dimnames(M) <- list(genes, need)                        # genes x cells
M <- M[Matrix::rowSums(M[, ref$Index] > 0) >= 5, ]      # genes detected in >= 5 reference cells
msg("count matrix: ", nrow(M), " genes x ", ncol(M), " cells")

# ---------------------------------------------------------------- 4. pseudo-bulk with truth
types <- sort(unique(ann$cell_type)); states <- sort(unique(ann$cell_state))
make_bulk <- function(cells, id) {
  lab <- tst[match(cells, tst$Index)]
  X <- M[, cells, drop = FALSE]
  umi <- Matrix::colSums(X)
  list(id = id, counts = Matrix::rowSums(X),
       type_umi = tapply(umi, factor(lab$cell_type, types), sum, default = 0),
       type_cell = tapply(rep(1, length(cells)), factor(lab$cell_type, types), sum, default = 0),
       state_umi = tapply(umi, factor(lab$cell_state, states), sum, default = 0),
       type_expr = sapply(types, function(t) { ix <- lab$cell_type == t; if (any(ix)) Matrix::rowSums(X[, ix, drop = FALSE]) else rep(0, nrow(X)) }))
}
bulks <- list()
for (s in sort(unique(tst$Sample))) bulks[[length(bulks) + 1]] <- make_bulk(tst[Sample == s, Index], paste0("natural_", s))
for (k in seq_len(opt$n_mix)) {
  p <- sample(test_pts, 1); pool <- tst[Patient == p]
  avail <- intersect(types, unique(pool$cell_type))
  w <- rgamma(length(avail), 1); w <- w / sum(w)
  n <- as.vector(rmultinom(1, opt$mix_cells, w))
  cells <- unlist(mapply(function(t, m) { ids <- pool[cell_type == t, Index]; if (m == 0) character(0) else ids[sample.int(length(ids), m, replace = TRUE)] },
                         avail, n, SIMPLIFY = FALSE), use.names = FALSE)
  bulks[[length(bulks) + 1]] <- make_bulk(cells, sprintf("random_%02d_%s", k, p))
}
ids <- sapply(bulks, `[[`, "id")
bk <- t(sapply(bulks, `[[`, "counts")); dimnames(bk) <- list(ids, rownames(M))
frac <- function(field) { x <- t(sapply(bulks, `[[`, field)); rownames(x) <- ids; sweep(x, 1, rowSums(x), "/") }
truth_type_umi <- frac("type_umi"); truth_type_cell <- frac("type_cell"); truth_state_umi <- frac("state_umi")
msg(nrow(bk), " pseudo-bulks (", sum(grepl("^natural", ids)), " natural, ", sum(grepl("^random", ids)), " random)")

# ---------------------------------------------------------------- 5. BayesPrism
sc <- as.matrix(t(M[, ref$Index]))
sc <- cleanup.genes(input = sc, input.type = "count.matrix", species = "hs",
                    gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)
sc <- select.gene.type(sc, gene.type = "protein_coding")
if (opt$genes == "markers") {
  st <- get.exp.stat(sc.dat = sc[, colSums(sc > 0) > 3], cell.type.labels = ref$cell_type, cell.state.labels = ref$cell_state,
                     pseudo.count = 0.1, cell.count.cutoff = 50, n.cores = opt$n_cores)
  sc <- select.marker(sc.dat = sc, stat = st, pval.max = 0.01, lfc.min = 0.1)
}
msg("reference: ", nrow(sc), " cells x ", ncol(sc), " genes; running BayesPrism on ", opt$n_cores, " cores")
prism <- new.prism(reference = sc, mixture = bk, input.type = "count.matrix", cell.type.labels = ref$cell_type,
                   cell.state.labels = ref$cell_state, key = "Tumor", outlier.cut = 0.01, outlier.fraction = 0.1)
res <- run.prism(prism = prism, n.cores = opt$n_cores)
est_type <- get.fraction(bp = res, which.theta = "final", state.or.type = "type")
est_state <- get.fraction(bp = res, which.theta = "first", state.or.type = "state")

# ---------------------------------------------------------------- 6. scoring
wt <- function(x, f) write.table(data.frame(sample = rownames(x), x, check.names = FALSE), file.path(opt$out, f), sep = "\t", quote = FALSE, row.names = FALSE)
wt(truth_type_umi, "truth_fractions_type.tsv"); wt(truth_state_umi, "truth_fractions_state.tsv")
wt(est_type, "est_fractions_type.tsv"); wt(est_state, "est_fractions_state.tsv")
score <- function(est, tru, level, truth_kind) {
  out <- list()
  for (c in intersect(colnames(est), colnames(tru))) {
    e <- est[ids, c]; t <- tru[ids, c]
    if (sd(t) == 0 || sd(e) == 0) next
    out[[c]] <- data.frame(level = level, truth = truth_kind, component = c, n = length(e), mean_truth = mean(t), mean_est = mean(e),
                           pearson = cor(e, t), spearman = cor(e, t, method = "spearman"), rmse = sqrt(mean((e - t)^2)), bias = mean(e - t))
  }
  do.call(rbind, out)
}
mf <- rbind(score(est_type, truth_type_umi, "type", "UMI fraction"), score(est_type, truth_type_cell, "type", "cell fraction"),
            score(est_state, truth_state_umi, "state", "UMI fraction"))
write.table(mf, file.path(opt$out, "metrics_fractions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

me <- list(); leak <- list()
z_all <- list()
for (ct in colnames(est_type)) {
  Z <- get.exp(bp = res, state.or.type = "type", cell.name = ct)
  z_all[[ct]] <- Z
  g <- colnames(Z); k <- match(ct, types)
  T <- t(sapply(bulks, function(b) b$type_expr[g, k])); rownames(T) <- ids
  Zs <- Z[ids, g]
  keep <- colMeans(T) >= 10
  gcor <- sapply(which(keep), function(j) suppressWarnings(cor(log1p(Zs[, j]), log1p(T[, j]))))
  pres <- rowSums(T) > 0
  scor <- sapply(which(pres), function(i) cor(log1p(Zs[i, keep]), log1p(T[i, keep])))
  me[[ct]] <- data.frame(cell_type = ct, genes = length(g), genes_scored = sum(keep), median_gene_pearson = median(gcor, na.rm = TRUE),
                         frac_genes_r_gt_0.5 = mean(gcor > 0.5, na.rm = TRUE), median_sample_pearson = median(scor, na.rm = TRUE),
                         normalized_abs_error = sum(abs(Zs - T)) / sum(T))
}
markers <- c("EPCAM", "CEACAM5", "KRT20", "PTPRC", "CD3E", "CD8A", "MS4A1", "LYZ", "COL1A1", "PECAM1")
tum_true <- t(sapply(bulks, function(b) b$type_expr[, "Tumor"])); colnames(tum_true) <- rownames(M); rownames(tum_true) <- ids
for (m in markers) {
  if (!m %in% colnames(z_all[["Tumor"]])) next
  tot_est <- sum(sapply(z_all, function(Z) sum(Z[ids, m])))
  leak[[m]] <- data.frame(gene = m, share_to_tumor_estimated = sum(z_all[["Tumor"]][ids, m]) / tot_est,
                          share_to_tumor_true = sum(tum_true[, m]) / sum(bk[, m]))
}
me <- do.call(rbind, me); leak <- do.call(rbind, leak)
write.table(me, file.path(opt$out, "metrics_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(leak, file.path(opt$out, "leakage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(opt$out, "fractions_scatter.png"), width = 1400, height = 800, res = 120)
cts <- colnames(est_type); par(mfrow = c(2, ceiling(length(cts) / 2)), mar = c(4, 4, 2, 1))
for (c in cts) {
  plot(truth_type_umi[ids, c], est_type[ids, c], pch = ifelse(grepl("^natural", ids), 17, 16),
       col = ifelse(grepl("^natural", ids), "firebrick", "grey30"), xlab = "true UMI fraction", ylab = "BayesPrism", main = c,
       xlim = c(0, 1), ylim = c(0, 1)); abline(0, 1, lty = 2)
}
dev.off()

el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
ft <- subset(mf, level == "type" & truth == "UMI fraction")
fmt <- function(d, cols) paste(apply(d[, cols], 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), collapse = "\n")
num <- function(d) { for (c in names(d)) if (is.numeric(d[[c]])) d[[c]] <- signif(d[[c]], 3); d }
rep <- c("# BayesPrism benchmark (single-cell pseudo-bulk ground truth)", "",
         sprintf("Data: GSE132465 (Lee et al. 2020, colorectal 10x). Held-out patients: %s. Pseudo-bulks: %d (%d natural, %d random). Genes: %s (%d). Runtime %.1f min on %d cores.",
                 paste(test_pts, collapse = ", "), length(ids), sum(grepl("^natural", ids)), sum(grepl("^random", ids)), opt$genes, ncol(sc), el, opt$n_cores), "",
         "## Cell-type fractions (BayesPrism theta vs true UMI fraction)", "",
         "| cell type | mean truth | mean est | Pearson | Spearman | RMSE | bias |", "|---|---|---|---|---|---|---|",
         fmt(num(ft), c("component", "mean_truth", "mean_est", "pearson", "spearman", "rmse", "bias")), "",
         "## Cell-type expression (BayesPrism Z vs true counts)", "",
         "| cell type | genes scored | median per-gene r | genes r > 0.5 | median per-sample r | normalized abs error |", "|---|---|---|---|---|---|",
         fmt(num(me), c("cell_type", "genes_scored", "median_gene_pearson", "frac_genes_r_gt_0.5", "median_sample_pearson", "normalized_abs_error")), "",
         "## Leakage check (share of each marker gene's reads assigned to Tumor)", "",
         "| gene | estimated | true |", "|---|---|---|", fmt(num(leak), c("gene", "share_to_tumor_estimated", "share_to_tumor_true")), "",
         "Best case: reference and pseudo-bulk share the 10x platform. Expect lower accuracy on real bulk RNA-seq.")
writeLines(rep, file.path(opt$out, "report.md"))
writeLines(c(paste0(names(opt), ": ", unlist(opt)), sprintf("minutes: %.1f", el)), file.path(opt$out, "run_info.txt"))
cat(paste(rep, collapse = "\n"), "\n")
msg("done in ", sprintf("%.1f", el), " min -> ", opt$out)
