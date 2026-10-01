library(data.table)

annot_file      <- "/path/to/variants_annotes.tsv"
constraint_file <- "/path/to/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

### =========================
### LOAD YOUR FINAL ANNOTATED TABLE (1 gene/variant, already picked)
### =========================
annot <- fread(annot_file, sep = "\t", header = TRUE)

### =========================
### LOAD pLI (gnomAD constraint file)
### =========================
constraint <- fread(constraint_file)
constraint_small <- constraint[, .(gene, pLI)]

constraint_small <- unique(constraint_small, by = "gene")

### =========================
### MERGE BY GENE SYMBOL
### =========================
result <- merge(
  annot,
  constraint_small,
  by.x = "Gene",
  by.y = "gene",
  all.x = TRUE
)

### =========================
### CHECK FOR UNMATCHED GENES 
### =========================
missing_genes <- unique(result[is.na(pLI), Gene])
if (length(missing_genes) > 0) {
  cat("Gene without pLI found (", length(missing_genes), "):\n")
  print(missing_genes)
  cat("-> check for HGNC synonyms or updated symbols.\n")
}

### =========================
### SAVE
### =========================
fwrite(result, "variants_annotes_with_pLI.tsv", sep = "\t", quote = FALSE)

cat("DONE\n")
cat("Lignes totales :", nrow(result), "\n")
cat("Lignes avec pLI trouvé :", sum(!is.na(result$pLI)), "\n")

