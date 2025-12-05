## code to prepare curatedMotifs dataset

library(universalmotif)
library(methods)

# 读取 universalMotif 格式数据
curated_motifs <- readRDS("data-raw/curatedMotif.rds")
raw_motifs <- readRDS("data-raw/rawMotif.rds")

# 读取元数据（移除 row.names = 1，因为CSV通常没有行名）
curated_motifs_meta <- read.csv("data-raw/curatedMotif_meta.csv")
raw_motifs_meta <- read.csv("data-raw/rawMotif_meta.csv")

# 读取peak
peaks <- read.table("data-raw/H3K4me3_Peaks.bed", header = F)

# 确保motif名称一致性
if ("motif_id" %in% names(curated_motifs_meta)) {
  names(curated_motifs) <- curated_motifs_meta$motif_id
} else {
  # 如果没有motif_id列，使用元数据的行名或其他唯一标识
  names(curated_motifs) <- curated_motifs_meta[,1]  # 使用第一列
}

if ("motif_id" %in% names(raw_motifs_meta)) {
  names(raw_motifs) <- raw_motifs_meta$motif_id
} else {
  names(raw_motifs) <- raw_motifs_meta[,1]
}

# 验证数据完整性
cat("Curated motifs - Number of motifs:", length(curated_motifs), "\n")
cat("Curated motifs - Number of metadata entries:", nrow(curated_motifs_meta), "\n")
cat("Raw motifs - Number of motifs:", length(raw_motifs), "\n")
cat("Raw motifs - Number of metadata entries:", nrow(raw_motifs_meta), "\n")
cat("peaks - Number of peaks:", nrow(peaks), "\n")

# 检查是否是universalMotif对象
cat("Curated motifs class:", class(curated_motifs[[1]]), "\n")
cat("Raw motifs class:", class(raw_motifs[[1]]), "\n")

# 保存为包数据
usethis::use_data(curated_motifs, overwrite = TRUE)
usethis::use_data(curated_motifs_meta, overwrite = TRUE)
usethis::use_data(raw_motifs, overwrite = TRUE)
usethis::use_data(raw_motifs_meta, overwrite = TRUE)
usethis::use_data(peaks, overwrite = TRUE)


