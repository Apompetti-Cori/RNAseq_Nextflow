md <- read.csv("~/data/03Aug2023_RRBS_20230703/129svev_dss_sample_metadata.csv")
sample_table_template <- read_csv("~/data/03Aug2023_RRBS_20230703/sample_table_template.csv")
sample_table_template[1:nrow(md),] <- NA
sample_table_template$sample <- md$sample_id
sample_table_template <- as.data.frame(sample_table_template)
rownames(sample_table_template) <- sample_table_template$sample

fastq <- list.files(path = "~/data/03Aug2023_RRBS_20230703/fastq", full.names = T)

for(sample in sample_table_template$sample){
  hit <- print(grep(pattern = sample, x = fastq, value = TRUE))
  sample_table_template[sample,]$r1_L1 <- hit[1]
  sample_table_template[sample,]$r2_L1 <- hit[2]
}

sample_table_template[is.na(sample_table_template)] <- ""

write.csv(sample_table_template, file = "~/data/03Aug2023_RRBS_20230703/sample_table.csv", quote = F, row.names = F)
