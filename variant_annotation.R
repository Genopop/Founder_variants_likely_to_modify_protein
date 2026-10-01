
#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# Automatic annotation of variants (format chr:pos:allele1:allele2, GRCh38)
# ------------------------------------------------------------------------
# For each variant, this script gets :
#   - gene (via Ensembl VEP)
#   - reference allele
#   - gene function (via NCBI Gene, "Summary" - same source as UCSC)
#   - associated diseases according to OMIM identification (via NCBI OMIM)
#   - ClinVar status 
#
# INSTALLATION :
#   install.packages(c("httr", "jsonlite"))
#   # on a cluster: module load r  (or equivalent), and run with the following command
#
# UTILISATION :
#   Rscript variant_annotation.R list_of_variants.txt
#
# list_of_variant : one variant per lign, format chr:pos:A1:A2
#   ex: chr10:100243705:A:C
#
# Output : variants_annotes.tsv
# ------------------------------------------------------------------------

suppressMessages({
  library(httr)
  library(jsonlite)
})

# ---------------------------- CONFIGURATION -----------------------------
NCBI_API_KEY   <- ""                 # optionnal 
NCBI_EMAIL     <- "you@example.com"  
VEP_BATCH_SIZE <- 200                
SEQ_BATCH_SIZE <- 50                
SLEEP_NCBI     <- if (nchar(NCBI_API_KEY) > 0) 0.11 else 0.34

ENSEMBL_SERVER <- "https://rest.ensembl.org"
NCBI_EUTILS    <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
G2P_FILE <- "/path/to/G2P.tsv" #Add phenotypes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript variant_annotation.R list_of_variants.txt")
}
input_path <- args[1]

# ------------------------------------------------------------------------
# 1. READING AND PARSING
# ------------------------------------------------------------------------
lines <- readLines(input_path, warn = FALSE)
lines <- trimws(lines)
lines <- lines[nchar(lines) > 0]

if (!grepl(":", lines[1])) {
  lines <- lines[-1]
}

parse_variant <- function(v) {
  parts <- strsplit(v, ":")[[1]]
  chrom <- sub("^chr", "", parts[1])
  pos   <- as.integer(parts[2])
  a1    <- parts[3]
  a2    <- parts[4]
  list(chrom = chrom, pos = pos, a1 = a1, a2 = a2)
}

variants_raw <- lines
parsed <- lapply(variants_raw, parse_variant)

cat(sprintf("Fichier lu: %d variants\n", length(parsed)))

# ------------------------------------------------------------------------
# 2. REFERENCE ALLELE
# ------------------------------------------------------------------------
ensembl_post <- function(path, body, retries = 3) {
  for (i in seq_len(retries)) {
    resp <- POST(
      url = paste0(ENSEMBL_SERVER, path),
      body = body, encode = "json",
      content_type("application/json"),
      add_headers(Accept = "application/json"),
      timeout(300)
    )
    if (status_code(resp) == 200) {
      return(fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                       simplifyVector = FALSE))
    }
    Sys.sleep(1)
  }
  warning(sprintf("Error requete Ensembl %s after %d tries (code %s)",
                   path, retries, status_code(resp)))
  return(NULL)
}

get_ref_alleles_batch <- function(parsed_batch) {
  regions <- sapply(parsed_batch, function(v) sprintf("%s:%d-%d", v$chrom, v$pos, v$pos))
  res <- ensembl_post("/sequence/region/human", list(regions = as.list(regions)))
  if (is.null(res)) return(rep(NA_character_, length(parsed_batch)))
  sapply(res, function(x) if (!is.null(x$seq)) toupper(x$seq) else NA_character_)
}

cat("1/3 - Determination of reference allele...\n")
ref_bases <- character(length(parsed))
for (i in seq(1, length(parsed), by = SEQ_BATCH_SIZE)) {
  idx <- i:min(i + SEQ_BATCH_SIZE - 1, length(parsed))
  cat(sprintf("  sequence: variants %d a %d / %d\n", min(idx), max(idx), length(parsed)))
  ref_bases[idx] <- get_ref_alleles_batch(parsed[idx])
}

for (i in seq_along(parsed)) {
  ref_base <- ref_bases[i]
  a1 <- parsed[[i]]$a1; a2 <- parsed[[i]]$a2
  if (!is.na(ref_base) && ref_base == a1) {
    parsed[[i]]$ref <- a1; parsed[[i]]$alt <- a2
  } else if (!is.na(ref_base) && ref_base == a2) {
    parsed[[i]]$ref <- a2; parsed[[i]]$alt <- a1
  } else {
    # aucun des deux alleles ne correspond au genome de reference
    # (multi-nucleotide, indel, ou erreur) -> on garde l'ordre original
    parsed[[i]]$ref <- a1; parsed[[i]]$alt <- a2
    parsed[[i]]$orientation_incertaine <- TRUE
  }
}

# ------------------------------------------------------------------------
# 3. ANNOTATION VEP (gene + ClinVar status if known)
# ------------------------------------------------------------------------
to_vep_region <- function(v) {
  end <- v$pos + nchar(v$ref) - 1
  sprintf("%s %d %d %s/%s 1", v$chrom, v$pos, end, v$ref, v$alt)
}

extract_gene_and_clinvar <- function(entry) {
  gene <- NA_character_
  tcs <- entry$transcript_consequences
  if (!is.null(tcs)) {
    for (tc in tcs) {
      if (!is.null(tc$gene_symbol)) { gene <- tc$gene_symbol; break }
    }
  }
  rsid <- NA_character_
  clin_sig <- NA_character_
  covs <- entry$colocated_variants
  if (!is.null(covs)) {
    for (cv in covs) {
      if (!is.null(cv$id) && startsWith(cv$id, "rs")) rsid <- cv$id
      if (!is.null(cv$clin_sig)) {
        clin_sig <- paste(unlist(cv$clin_sig), collapse = ", ")
      }
    }
  }
  list(gene = gene, rsid = rsid, clinvar_status = clin_sig)
}

cat("2/3 - Annotation VEP...\n")
vep_results <- vector("list", length(parsed))
for (i in seq(1, length(parsed), by = VEP_BATCH_SIZE)) {
  idx <- i:min(i + VEP_BATCH_SIZE - 1, length(parsed))
  cat(sprintf("  VEP: variants %d a %d / %d\n", min(idx), max(idx), length(parsed)))
  regions <- lapply(parsed[idx], to_vep_region)
  res <- ensembl_post("/vep/human/region", list(variants = regions, pick = 1))
  if (is.null(res)) {
    for (j in idx) vep_results[[j]] <- list(gene = NA, rsid = NA, clinvar_status = NA)
  } else {
    for (k in seq_along(idx)) {
      vep_results[[idx[k]]] <- extract_gene_and_clinvar(res[[k]])
    }
  }
  Sys.sleep(0.2)
}

# ------------------------------------------------------------------------
# 4. GENE FUNCTION + OMIM (NCBI Gene / OMIM)
# ------------------------------------------------------------------------
ncbi_get <- function(path, query) {
  query$retmode <- "json"
  query$email <- NCBI_EMAIL
  if (nchar(NCBI_API_KEY) > 0) query$api_key <- NCBI_API_KEY
  resp <- GET(paste0(NCBI_EUTILS, path), query = query, timeout(30))
  stop_for_status(resp)
  fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE)
}

get_gene_id <- function(symbol) {
  res <- tryCatch(
    ncbi_get("/esearch.fcgi", list(db = "gene", term = sprintf("%s[sym] AND human[orgn]", symbol))),
    error = function(e) NULL
  )
  ids <- res$esearchresult$idlist
  if (length(ids) == 0) return(NA_character_)
  ids[1]
}

get_gene_summary_and_omim <- function(gene_id) {
  res <- tryCatch(
    ncbi_get("/esummary.fcgi", list(db = "gene", id = gene_id)),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(summary = "", mim_ids = character(0)))
  doc <- res$result[[gene_id]]
  summary <- if (!is.null(doc$summary)) doc$summary else ""
  mim_ids <- if (!is.null(doc$mim)) unlist(doc$mim) else character(0)
  list(summary = summary, mim_ids = mim_ids)
}

get_omim_titles <- function(mim_ids) {
  if (length(mim_ids) == 0) return(character(0))
  res <- tryCatch(
    ncbi_get("/esummary.fcgi", list(db = "omim", id = paste(mim_ids, collapse = ","))),
    error = function(e) NULL
  )
  if (is.null(res)) return(character(0))
  titles <- c()
  for (mid in mim_ids) {
    doc <- res$result[[mid]]
    if (!is.null(doc$title) && nchar(doc$title) > 0) {
      titles <- c(titles, sprintf("%s (OMIM:%s)", doc$title, mid))
    }
  }
  titles
}

cat("3/3 - Gene functions et OMIM data...\n")
genes <- sapply(vep_results, function(x) x$gene)
unique_genes <- sort(unique(na.omit(genes)))

gene_cache <- list()
for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  cat(sprintf("  Gene %d/%d: %s\n", i, length(unique_genes), gene))
  gid <- get_gene_id(gene)
  Sys.sleep(SLEEP_NCBI)
  if (is.na(gid)) {
    gene_cache[[gene]] <- list(summary = "", diseases = "")
    next
  }
  info <- get_gene_summary_and_omim(gid)
  Sys.sleep(SLEEP_NCBI)
  titles <- get_omim_titles(info$mim_ids)
  Sys.sleep(SLEEP_NCBI)
  gene_cache[[gene]] <- list(summary = info$summary, diseases = paste(titles, collapse = "; "))
}

################################Add phenotypes#################################
# ------------------------------------------------------------------------
# 4b. GENE2PHENOTYPE (G2P)
# ------------------------------------------------------------------------

cat("Loading Gene2Phenotype data...\n")

g2p <- read.csv(
  G2P_FILE,
  sep = ",",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(sprintf("G2P charge: %d associations, %d colonnes\n",
            nrow(g2p), ncol(g2p)))
print(colnames(g2p))

g2p <- g2p[g2p$confidence %in%
            c("definitive", "strong", "moderate"), ]


get_g2p <- function(gene){

  res <- g2p[g2p$`gene symbol` == gene, ]

  if(nrow(res) == 0){
    return("")
  }

  annotations <- paste0(
    res$`disease name`,
    " [",
    res$`allelic requirement`,
    "; ",
    res$confidence,
    "; ",
    res$`molecular mechanism`,
    "]"
  )

  paste(unique(annotations), collapse="; ")
}
###########################################Add phenotypes####################################################
# ------------------------------------------------------------------------
# 5. FINAL TABLE
# ------------------------------------------------------------------------
cat("Final table...\n")
out_rows <- lapply(seq_along(parsed), function(i) {
  v <- parsed[[i]]
  vres <- vep_results[[i]]
  gene <- vres$gene
  ginfo <- if (!is.na(gene) && gene %in% names(gene_cache)) gene_cache[[gene]] else list(summary = "", diseases = "")
#  data.frame(
#    Variant = variants_raw[i],
#    Ref_reel = v$ref,
#    Alt_reel = v$alt,
#    Orientation_incertaine = isTRUE(v$orientation_incertaine),
#    Gene = ifelse(is.na(gene), "", gene),
#    rsID = ifelse(is.na(vres$rsid), "", vres$rsid),
#    Fonction_du_gene = ginfo$summary,
#    Maladies_associees_OMIM = ginfo$diseases,
#    Statut_ClinVar = ifelse(is.na(vres$clinvar_status), "Non trouve dans ClinVar", vres$clinvar_status),
#    stringsAsFactors = FALSE
#  )
data.frame(
    Variant = variants_raw[i],
    Ref_reel = v$ref,
    Alt_reel = v$alt,
    Orientation_incertaine = isTRUE(v$orientation_incertaine),
    Gene = ifelse(is.na(gene), "", gene),
    rsID = ifelse(is.na(vres$rsid), "", vres$rsid),
    Fonction_du_gene = ginfo$summary,
    Maladies_associees_OMIM = ginfo$diseases,
    Maladies_G2P = get_g2p(gene),
    Statut_ClinVar = ifelse(is.na(vres$clinvar_status),
                            "Non trouve dans ClinVar",
                            vres$clinvar_status),
    stringsAsFactors = FALSE
)
})

out_df <- do.call(rbind, out_rows)
out_path <- "variants_annotes_ad.tsv"
write.table(out_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nTermine ! Resultat ecrit dans : %s\n", out_path))
