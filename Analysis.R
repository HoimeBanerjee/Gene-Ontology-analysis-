# ==============================================================================
# Step 1: Load Required Packages (Prevent Source Compilation)
# ==============================================================================
options(install.packages.compile.from.source = "never")
options(pkgType = "binary")

required_cran <- c("readxl", "openxlsx", "ggplot2", "stringr", "BiocManager")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, type = "binary")
  }
}

required_bioc <- c("clusterProfiler", "org.Dr.eg.db", "enrichplot")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(readxl)
library(openxlsx)
library(clusterProfiler)
library(org.Dr.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)

# ==============================================================================
# Step 2: Read Excel File & Clean Accession Numbers
# ==============================================================================
input_path          <- "E:/APS 3.5/miRNA computational data/Downregulated match.xlsx"
output_path         <- "E:/APS 3.5/miRNA computational data/Downreg.xlsx"
excel_analysis_path <- "E:/APS 3.5/miRNA computational data/Down reg analysis.xlsx"

# Read original Excel file
df <- read_excel(input_path)

# Clean accession numbers in Column B ("Representative 3' UTR") e.g. .1 suffix
df$`Representative 3' UTR` <- sub("\\.\\d+$", "", df$`Representative 3' UTR`)
write.xlsx(df, file = output_path)

ensembl_ids <- na.omit(df$`Representative 3' UTR`)
cat("\n[SUCCESS] Cleaned Ensembl IDs saved to:\n", output_path, "\n\n")

# ==============================================================================
# Step 3: Run GO Analysis & Remove Redundant Parent-Child Terms
# ==============================================================================
# 1. Full GO Biological Process Run
go_res_full <- enrichGO(
  gene          = ensembl_ids,
  OrgDb         = org.Dr.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",            # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 1.0,             # Keep all terms to compute cutoff statistics
  qvalueCutoff  = 1.0
)

# 2. Translate Ensembl IDs to readable Gene Symbols
go_res_full <- setReadable(go_res_full, OrgDb = org.Dr.eg.db, keyType = "ENSEMBL")

# 3. Simplify Redundant GO Terms (Merge terms sharing identical gene sets)
go_res_simplified <- simplify(
  go_res_full,
  cutoff     = 0.7,               # Semantic similarity cutoff
  by         = "pvalue",
  select_fun = min
)

all_simplified_df <- as.data.frame(go_res_simplified)

# Separate Included (Nominal p < 0.05) vs Excluded (p >= 0.05)
p_cutoff    <- 0.05
included_df <- all_simplified_df[all_simplified_df$pvalue < p_cutoff, ]
excluded_df <- all_simplified_df[all_simplified_df$pvalue >= p_cutoff, ]

# Select Top 15 non-redundant terms sorted by nominal p-value
included_df <- head(included_df[order(included_df$pvalue), ], 15)

# Format label text:
# Line 1: Count + Gene names on a SINGLE straight line
# Line 2: Nominal p-value on its own line underneath
included_df$Gene_Single_Line <- gsub("/", ", ", included_df$geneID)
included_df$Label_Text <- paste0(
  included_df$Count, " [", included_df$Gene_Single_Line, "]\n",
  "p = ", sprintf("%.4f", included_df$pvalue)
)

# ==============================================================================
# Step 4: Streamlined Significance-Ranked Plot (X-Axis = 0 to 5)
# ==============================================================================
if (nrow(included_df) > 0) {
  
  # reorder(Description, -pvalue) ranks terms by significance (most significant at top)
  p_simplified <- ggplot(included_df, aes(x = Count, y = reorder(Description, -pvalue), fill = pvalue)) +
    geom_col(width = 0.68, color = "black", linewidth = 0.2) +
    # Draw gene labels + p-values cleanly next to bars
    geom_text(
      aes(label = Label_Text), 
      hjust = 0, 
      nudge_x = 0.08, 
      size = 2.6, 
      fontface = "bold", 
      color = "#1A252C",
      lineheight = 0.95
    ) +
    scale_fill_gradient(low = "#e41a1c", high = "#377eb8", name = "Nominal p-val") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    # Set X-axis scale strictly from 0 to 5
    scale_x_continuous(
      limits = c(0, 5), 
      breaks = 0:5, 
      expand = c(0, 0)
    ) + 
    labs(
      title = "Non-Redundant GO Biological Process Clusters",
      subtitle = paste0("Ranked by Significance (Similarity Cutoff = 0.7 | Nominal p < ", p_cutoff, ")"),
      x = "Number of Target Genes (Count)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, italic = TRUE, color = "grey30"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 9.5, face = "bold"),
      legend.position = "right",
      # Left padding shows full cluster titles; right padding keeps legend tight
      plot.margin = margin(t = 12, r = 20, b = 12, l = 170, unit = "pt")
    )
  
  dev.new()
  print(p_simplified)
  
  # Save PNG with 13.5-inch width for clean rendering
  plot_path <- "E:/APS 3.5/miRNA computational data/GO_NonRedundant_Barplot.png"
  ggsave(plot_path, plot = p_simplified, width = 13.5, height = 7.5, dpi = 300)
  cat("[SUCCESS] Saved significance-ranked plot to:\n", plot_path, "\n\n")
}

# ==============================================================================
# Step 5: Export 'Down reg analysis.xlsx' Workbook
# ==============================================================================
wb <- createWorkbook()

# Sheet 1: Overview
addWorksheet(wb, "Overview")
overview_data <- data.frame(
  Metric = c("Organism", "Input Gene File", "Ontology Evaluated", 
             "Significance Cutoff Applied", "Redundancy Filter (simplify)", 
             "Total Input Genes", "Included Non-Redundant Clusters (p < 0.05)", 
             "Excluded Clusters (p >= 0.05)"),
  Value = c("Danio rerio (Zebrafish)", "Downreg.xlsx", "Biological Process (BP)",
            "Nominal p-value < 0.05", "Cutoff = 0.7 (Semantic Similarity)",
            length(ensembl_ids), nrow(included_df), nrow(excluded_df))
)
writeData(wb, "Overview", overview_data)

# Sheet 2: Enriched Clusters
addWorksheet(wb, "Enriched Clusters")
writeData(wb, "Enriched Clusters", 
          included_df[, c("ID", "Description", "Count", "GeneRatio", "pvalue", "p.adjust", "geneID")])

# Sheet 3: Excluded Clusters
addWorksheet(wb, "Excluded Clusters")
if (nrow(excluded_df) > 0) {
  excluded_df$Exclusion_Reason <- paste0("Nominal p-value (", round(excluded_df$pvalue, 4), ") >= 0.05 cutoff")
  writeData(wb, "Excluded Clusters", 
            excluded_df[, c("ID", "Description", "Count", "pvalue", "p.adjust", "Exclusion_Reason", "geneID")])
}

saveWorkbook(wb, excel_analysis_path, overwrite = TRUE)
cat("[SUCCESS] Created analysis workbook:\n", excel_analysis_path, "\n")