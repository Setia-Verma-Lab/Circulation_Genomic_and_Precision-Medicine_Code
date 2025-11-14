library(dplyr)
library(PheWAS)
library(tidyr)
library(logistf)
library(meta)

pheno <- read.delim(file = "Deidentified_Phecode_File.txt",
                    sep = '\t', header = TRUE, check.names = FALSE) %>%
  rename(PMBB_ID = id)

no_geno_variability <- c()

filter <- function(geno, geno_summary) {
  for (c in 7:ncol(geno)){
    x <- as.matrix(geno[c])
    geno_summary[c-6,1] <- colnames(x)
    x_table <- as.data.frame(table(x))
    geno_summary[c-6,2:(nrow(x_table)+1)] <- x_table$Freq
    # Remove variants that have more than 5% NA for genotype
    if (sum(is.na(x)) > 2087){
      geno[c] <- 0
      geno_summary[c-6,2] <- 'REMOVED_NA'
    }
    # Remove non-rare variants (>1% MAF)
    if (sum(x, na.rm = TRUE) > 835){
      geno[c] <- 0
      geno_summary[c-6,2] <- 'REMOVED_Common'
    }
  }
  geno_summary[is.na(geno_summary)] <- 0
  return(geno_summary)
}

for (file in c('genotypes_8:105802051:A:G.raw')){
  # gene_name <- substr(unlist(strsplit(file,"genotypes_"))[1], 1, -4)
  gene_name <- substr(file, 11, nchar(file)-4)
  print(paste0("Processing ", gene_name, "..."))

  geno <- read.delim(paste0(file),
                     sep = " ", header = TRUE, check.names = FALSE) %>%
    rename(PMBB_ID = FID)
  geno_summary <- data.frame(variant = numeric(), n_0 = numeric(), n_1 = numeric(), n_2 = numeric())

  geno_summary <- filter(geno,geno_summary)

  # Calculate gene burden
  geno <- geno %>%
    mutate(lof=rowSums(.[7:length(colnames(.))], na.rm=T)) %>%
    dplyr::select(PMBB_ID,lof)
  ###PHEWAS###
  print("Performing PheWAS...")

  results_all <- phewas(phenotypes = pheno, genotypes = geno)
  results_all$study<-'all'

  print("Printing plots...")
  if (sum(is.na(results_all$beta)) == nrow(results_all)){
    no_geno_variability <- append(no_geno_variability, paste0(gene_name,"_ALL"))
  } else {
    pdf(file=paste0("PheWasOut/",gene_name,"_phewas_plot_rollup_all.pdf"),
        paper='USr', width=11, height=8.5)
    myplot <- phewasManhattan(results_all, OR.direction=T, annotate.angle=0, annotate.size=4,annotate.level=0.05) +
      ggtitle(paste0(gene_name," PheWAS")) +
      theme(plot.title = element_text(hjust=0.5, size=20),
            axis.text.x = element_text(angle=45, hjust=1),
            plot.margin = margin(0,0,0,20))
    print(myplot)
    dev.off()
  }
  results_all2 <- addPhecodeInfo(results_all)
  write.csv(results_all2, file=paste0("PheWasOut/",gene_name,"_phewas_rollup_all.csv"))
}
