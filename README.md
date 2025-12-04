[CuratedMotifs.md](https://github.com/user-attachments/files/23922706/CuratedMotifs.md)
---
curatedMotif
---





- 3255 High-quality Non-Redundant DNA binding motifs for ~1400 mouse/human genes.
- supported by a transparent curation workflow. 
- Seamless integration with sequencing data, such as ChIP-seq, ATAC-seq, and snATAC-seq.



Data Sources and Workflow



1. Original motifs were obtained from cisBP, JASPAR, Encode, and homer

   - `cisBP_v3`:  
     1. Only human and mouse motifs were retained. 
     2. HocoMoco and previous JASPAR versions motifs were removed.

   - `JASPAR2024`: vertebrates core
   - `Encode`: https://compbio.mit.edu/encode-motifs/motifs.txt
   - `homer`:  vertebrates

   

2. Gene symbols were updated to the most current standard using: Google, GeneCards, and  GENCODE annotations 

   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.metadata.HGNC.gz

   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.metadata.MGI.gz

   

3. Motif Processing Workflow (using the R package universalmotif)

   ​       filter ➡️ trim ➡️ merge ➡️ curation

   1. **Quality Filtering:** We filtered out motifs with low total information content (IC) scores to retain only high-quality, informative sequences.
   2. **Flank Trimming:** We trimmed low-information-content bases from both ends of each motif.
   3. **Similarity-based Merging:** We clustered and merged redundant motifs associated with the same gene to generate a non-redundant set.
   4. **Expert Curation:** We manually curated the results of each merge to validate and refine the final motif models.



## motif Database

**Requirements for Motif Processing:**

1. **Motif Naming**

   - Each motif must have a **unique name**.

   - Motif Naming Guidelines:
     1. Use special characters sparingly to avoid errors during single-cell analysis.
     2. Permitted characters: `.` and `-` are acceptable.
     3. Avoid `_`: They are automatically converted to `-` by chromVAR, which may cause inconsistent motif naming.

2. **Metadata Table**

   Format: a data frame or table with exactly two columns:

   - `motif_name`: Must **exactly match** the names in the motif list.
   - `gene`: Each motif must be annotated with its corresponding **official gene symbol**.  

```r
library(tidyverse)
library(universalmotif)
library(ggpubr)

outDir <- 'E:/motif/curatedMotif'
figDir <- paste0(outDir, '/png')
dir.create(figDir, recursive = T)
setwd(outDir)

theme1 <- theme(legend.position = 'none', axis.text.x = element_text(size = 6)  )
```

### cisBP

**meta中gene为OTP的不是真正的gene Symbol**, 部分可以注释; 但大部分为fly gene, 直接删去了这部分难以注释的motif

```r
dirPWM <- 'E:/motif/cisbpV3/pwms'
aa <- list.files(dirPWM, full.names = T)
nameX <- gsub('_.*$', '', basename(aa))
names(aa) <- nameX

cisBP <- list()
for (i in nameX) {
  bb <-  read.table(aa[i])
  if (nrow(bb) > 2 ) { 
      cisBP[[i]] <-  read_cisbp(aa[i]) 
      cisBP[[i]]@name <- i  } } # cisBP有的motif是空的

metaBP <- data.table::fread("E:/motif/cisbpV3/TF_Information_all_motifs.txt")
metaBP <- metaBP[metaBP$Motif_ID != '.', ]
metaBP <- metaBP[metaBP$TF_Species %in% c( 'Homo_sapiens', 'Mus_musculus'), ]
metaBP$Motif_ID <- gsub('_3.00', '', metaBP$Motif_ID)
metaBP <- metaBP[!duplicated(metaBP$Motif_ID), ]
metaBP <- metaBP[order(metaBP$gene), ]
metaBP <- as.data.frame(metaBP)
metaBP <- metaBP[, c("Motif_ID", "gene", "MSource_Identifier")]
metaBP <- metaBP[!metaBP$MSource_Identifier %in% c('HocoMoco', 'JASPAR'), ] # 去掉'HocoMoco', 'JASPAR'
rownames(metaBP) <- metaBP$Motif_ID

cisBP <- cisBP[metaBP$Motif_ID]
```

### JASPAR

```r
jas2024 <- read_jaspar( 'E:/motif/JASPAR2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt')
jas2024 <- convert_type(motifs = jas2024, type = 'PPM')

motifID <- vector()
geneSymbol <- vector()

for (i in 1:length(jas2024) ) {
  motifID[i] <- jas2024[[i]]@name
  geneSymbol[i] <- jas2024[[i]]@altname
}

metaJAS <- data.frame('Motif_ID'= motifID, 
                      'gene'= geneSymbol, 
                      "MSource_Identifier" = "JASPAR" )
rownames(metaJAS) <- metaJAS$Motif_ID
names(jas2024) <- metaJAS$Motif_ID
```

### encode motif

```r
# encode motif 标准命名, 但部分错误  https://compbio.mit.edu/encode-motifs/
motif <- read_matrix('E:/motif/encode/encode2013.txt', type = 'PPM', positions = 'row', alphabet = 'DNA', sep = NULL, headers = '>')
# 不要使用 chromVAR包的encode_pwms, 矩阵的数值偏小
orgName <- sapply(motif, function(x) x@name)
meta <- as.data.frame(orgName)
meta <- separate(meta, 'orgName', remove = F, into = c('orgMotif', 'source'), sep = '\t')
rownames(meta) <- meta$orgMotif
meta$order <- 1:nrow(meta)
meta$chromVAR <- names(encode_pwms@listData)
meta$geneVAR <- sapply(encode_pwms@listData, function(x) x@name )

meta$gene <- meta$geneVAR
aa <- meta[ !grepl('::', meta$gene ) & !toupper(meta$gene  ) %in% toupper(c(human$V2, mouse$V2)) ,] 
aa$gene <- gsub('_.*$', '', aa$sourceVAR)

meta[ !grepl('::', meta$geneVAR) & !toupper(meta$geneVAR ) %in% toupper(c(human$V2, mouse$V2)) ,]$gene <- aa$gene

aa <- meta[ !grepl('::', meta$gene ) & !toupper(meta$gene  ) %in% toupper(c(human$V2, mouse$V2)) ,] 
meta[meta$gene=='NR1H',]$gene <- 'NR1H2' # 许多来自transfac的motif只能指导protein family, 无法指导具体gene

meta$motif <- make.unique(meta$gene, sep = '.e')
meta[-grep('\\.e', meta$motif), ]$motif <- paste0( meta[-grep('\\.e', meta$motif), ]$motif, '.e0')
rownames(meta) <- meta$motif
names(motif) <- meta$motif

for (i in names(motif) ) {
  motif[[i]]@extrainfo <- gsub('\t', ' ', motif[[i]]@name)
  motif[[i]]@name <- i
}

encode$meta <- meta
encode$motif <- motif
```

### homer motif

```
motif <- read_homer('E:/motif/homer/allVerbrate.txt' )
geneHomer <- sapply(motif, function(x) x@name)
meta <- as.data.frame(geneHomer)
meta$gene <- meta$geneHomer
```

## combine all collections

```r
motif <- c( cisBP3$motif, JASPAR2024$motif, homer$motif, encode$motif)

meta <- rbind(cisBP3$meta[, c('motif', 'gene', 'source')],  
                         JASPAR2024$meta[, c('motif', 'gene', 'source')], 
                         homer$meta[, c('motif', 'gene', 'source')], 
                         encode$meta[, c('motif', 'gene', 'source')] )

for (i in names(preCurated$motif) ) { preCurated$motif[[i]]@name <- i }

# 排序
meta <- meta[order(meta$gene),] 
motif <- motif[rownames(meta)]

preCurated <- list()
preCurated$meta <- meta
preCurated$motif <- motif
```

## official gene symbol

**手动注释 名为`OTP`的motif** 

```r
human <- read.delim("E:/motif/gencode.v49.metadata.HGNC", header=FALSE)
mouse <- read.delim("E:/motif/gencode.vM38.metadata.MGI", header=FALSE)

#找出需要校正的gene
meta$order <- 1:nrow(meta)
aa <- meta[ !grepl('--', meta$gene) & !toupper(meta$gene ) %in% toupper(c(human$V2, mouse$V2)) ,] 
# 逐个校正gene name
meta$gene <- ifelse( grepl('ARNTL', meta$gene )  , 'BMAL1', meta$gene )
```

## motif quality

```r
# 必须先过滤，如果直接trim， 有一些会全部修剪掉的motif会导致出错
for (i in names(motif) ) {motif[[i]]@name <- i }

motif <- filter_motifs( motifs = motif, icscore = 6 ) 
motif <- trim_motifs(motifs = motif, min.ic = 0.4, trim.from ='both' ) #
motif <- filter_motifs( motifs = motif, icscore = 6 ) 

meta <- meta[names(motif),] 


keep <- meta[ grepl('--', meta$gene) | toupper(meta$gene ) %in% toupper(c(human$V2, mouse$V2)) ,]$motif
motif <- motif[keep]
meta <- meta[keep,]
```

## First-round merging

Fixed-parameter-based merging of similar motifs within the same genes.

<img src="https://typora-1329573677.cos.ap-shanghai.myqcloud.com/img/Nr1h2_motif_rP.png" alt="Nr1h2_motif_rP" style="zoom:25%;" />

```r
figDir1 <- paste0(outDir, '/png1')
figDir2 <- paste0(outDir, '/png2')
dir.create(figDir1, recursive = T)
dir.create(figDir2, recursive = T)

motifMerge <- list( )
for (i in unique(meta$gene) ) {
  gene1 = i
  motifx <- motif[ meta[meta$gene==i,]$motif]
  mergex <- list()
  if (length(motifx) == 1) {  # 该gene只有一个motif

    mergex <- motifx
    names(mergex) <- paste0(gene1, '.m0'  )
    plotM <- view_motifs(motifs = mergex ) +
      theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )# 颜色顺序ACGT
    ggsave( paste0( gene1, "_motif.png"), plot = plotM, path = figDir1, scale = 1.5, width = 10 , 
            height = 2, units = "cm", dpi = 150, limitsize = F )
    
  } else {# 该gene有多个motif
    names(motifx) <- paste0("N", 1:length(motifx), "_", names(motifx))
    for (i in names(motifx)) { motifx[[i]]@name <- i }
    mergex1 <- merge_similar(motifs = motifx, threshold = 0.6, method = 'PCC', normalise.scores = T, min.overlap = 0.6, tryRC = T) # 
    mergex1 <- trim_motifs(motifs = mergex1, min.ic = 0.4, trim.from = 'both' )
    if (length(mergex1)==1) { # 合并后只有一个motif
      mergex[[ paste0(gene1, '.m0') ]] <- mergex1 
    } else { # 合并后有都多个motif
      mergex <- mergex1
      names(mergex) <-  paste0(gene1, '.m', 0:(length(mergex)-1) )
    }
    plot1 <- view_motifs(motifs = motifx) +
      theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )
    plot2 <- view_motifs(motifs = mergex) +
      theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") ) 
    plotM <- ggarrange(plot1, ggarrange(plot2, 'blank', ncol = 1, heights = c(length(mergex), (length(motifx) - length(mergex))*0.5 )) ,
                       nrow = 1, widths = c(1,1.5) )
    ggsave( paste0( gene1, "_motif.png"), plot = plotM, path = figDir2, scale = 1.5, width = 20, 
          height = length(motifx)*1.2 , units = "cm", dpi = 150, limitsize = F )
  }
  
  motifMerge <- c(motifMerge, mergex)
}
```

## Manual curation

Based on the image inspection, parameters will be adjusted and re-merged for unreasonable clustering results.

<img src="https://typora-1329573677.cos.ap-shanghai.myqcloud.com/img/Olig2_motif_1.png" alt="Olig2_motif_1" style="zoom:25%;" />

```r
gene1 <- 'DMRT1' # 逐个校对
motifx <- motif[meta[ grepl(pattern = paste0('^', gene1, '$'), meta$gene, ignore.case = T ), ]$motif ]
names(motifx) <- paste0("N", 1:length(motifx), "_", names(motifx))
for (i in names(motifx)) { motifx[[i]]@name <- i }
plot1 <- view_motifs(motifs = motifx, method = 'PCC', tryRC = T, min.overlap = 0.5, normalise.scores =T ) +
  theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )
# 去掉低质量的motif,调整threshold等参数
mergex1 <- merge_similar(motifs = motifx[ -c(4:6) ], 
                         threshold = 0.7, method = 'PCC', normalise.scores = T, min.overlap = 0.5, tryRC = T) # 
mergex1 <- trim_motifs(motifs = mergex1, min.ic = 0.6, trim.from = 'both' ) # 

plot2 <- view_motifs(motifs = mergex1, method = 'PCC', tryRC = T, min.overlap = 0.5, normalise.scores =T) +
  theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") ) # 
plotM <- ggarrange(plot1, 
                   ggarrange(plot2, 'blank', ncol = 1, heights = c(length(mergex1), 0.5*(length(motifx) - length(mergex1)) ) ), 
                   nrow = 1, widths = c(1,1.5) )
 
ggsave( paste0( gene1, "_motif.png"), plot = plotM, path = figDir2, scale = 1.5, width = 12, 
        height = length(motifx), units = "cm", dpi = 150, limitsize = F )

motifMerge[grep(pattern = paste0('^', gene1, ".m[0-9]" ), names(motifMerge), ignore.case = T)]  <- NULL     #删除不合理的合并

mergex <- list()
if (length(mergex1)==1) {
  mergex[[ paste0(gene1, '.m0')  ]] <- mergex1 
} else{
  mergex <- mergex1
  names(mergex) <-  paste0(gene1, '.m', 0:(length(mergex)-1) )
}

motifMerge <- c(motifMerge, mergex)
grep(pattern = paste0('^', gene1, ".m[0-9]" ), names(motifMerge), ignore.case = T)
```

```r
mergex1[[6]] <- merge_motifs( motifx[c(5,6)],method = 'PCC', normalise.scores = T, min.overlap = 0.5, tryRC = T )  
mergex1[[4]] <- merge_motifs(motifx[sapply(mergex1[c(4,5)], function(x) x@name) %>% strsplit('/') %>% unlist()] )
mergex1[c(5 )] <- NULL

motifx[sapply(mergex1[c(1 ) ], function(x) x@name) %>% strsplit('/') %>% unlist()]  %>%  
  view_motifs( method = 'PCC', tryRC = T, min.overlap = 0.5, normalise.scores =T ) +
  theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )

merge_motifs( motifx[c( 1,4,6,8 )], method = 'PCC', normalise.scores = T, min.overlap = 0.8, tryRC = T )   %>%  
  view_motifs( method = 'PCC', tryRC = T, min.overlap = 0.5, normalise.scores =T ) +
  theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") )
```

## rename motif

```r
motifMerge <- motifMerge[!duplicated(names(motifMerge))]

curatedMotif <- motifMerge[order(names(motifMerge))]

metaMerge <- data.frame('motif'=names(curatedMotif))
rownames(metaMerge) <- metaMerge$motif
for (i in rownames(metaMerge)) { metaMerge[i, 'sourceMotif'] <- curatedMotif[[i]]@name}

for (i in names(curatedMotif)) {
  curatedMotif[[i]]@extrainfo <- curatedMotif[[i]]@name
  curatedMotif[[i]]@name <- i
}
```

## final QC

```r
curatedMotif <- trim_motifs(curatedMotif, min.ic = 0.6, trim.from = 'both')
curatedMotif <- filter_motifs(curatedMotif, width = 4, icscore = 6)
```

## Versus other databases

<img src="https://typora-1329573677.cos.ap-shanghai.myqcloud.com/img/Sox2_motif2.png" alt="Sox2_motif2" style="zoom:25%;" />

```r
library(chromVARmotifs)
data(mouse_pwms_v2) #

chromVAR2 <- convert_motifs(motifs = mouse_pwms_v2@listData )
chromVAR2 <- convert_type(motifs = chromVAR2, type = 'PPM')

metaV2 <- data.frame()
for (i in names(chromVAR2) ) {
  metaV2[i, 'gene'] <- chromVAR2[[i]]@name
  metaV2[i, 'motif'] <- paste0( gsub( paste0( '^.*', chromVAR2[[i]]@name), '', i ) )
  metaV2[i, 'ensembl'] <- chromVAR2[[i]]@extrainfo
}

aa <- metaV2[grep('LINE', metaV2$gene),]$motif # 错误的gene 
aa <- gsub('NP_', 'NP', aa ) # 名称为NP的
aa <- sub('_', '', aa)
geneV2 <- sub('_.*$', '', aa ) 
motifV2 <- sub('^[^_]*', '', aa ) 
metaV2[grep('LINE', metaV2$gene),]$motif <- motifV2
metaV2[grep('LINE', metaV2$gene),]$gene <- geneV2
metaV2$motif <- paste0(metaV2$gene, '.v2', metaV2$motif)
 
names(chromVAR2) 
for (i in names(chromVAR2))  { chromVAR2[[i]]@name <- i }

# 
chromVAR1 <- convert_motifs(motifs = mouse_pwms_v1@listData )
chromVAR1 <- convert_type(motifs = chromVAR1, type = 'PPM')


HocoMoco <- read_jaspar('E:/motif/HocoMocoV13/H13CORE_jaspar_format.txt')
motifName <- vector()
for (i in 1:length(HocoMoco)) { motifName[i] <- HocoMoco[[i]]@name }
names(HocoMoco) <-  motifName
HocoMoco <- convert_type(motifs = HocoMoco, type = 'PPM')

JASPAR2024X <- RSQLite::dbConnect(RSQLite::SQLite(),JASPAR2024::db(JASPAR2024::JASPAR2024()))
JASPAR2024 <- TFBSTools::getMatrixSet( x = JASPAR2024X, opts = list(species = c('Homo sapiens','Mus musculus' ) , all_versions = F) )
JASPAR2024 <- JASPAR2024 %>% convert_motifs() %>% convert_type(type = 'PPM')

for (i in 1:length(JASPAR2024)) { names(JASPAR2024 )[i] <- paste0(JASPAR2024[[i]]@name, '.', names(JASPAR2024 )[i] ) }
for (i in names(JASPAR2024)) {JASPAR2024[[i]]@name <- i}


motif <- c(curatedMotif$motif, HocoMoco, chromVAR1, chromVAR2, JASPAR2024)

gene1 <-  'Sox2'
motifx <- motif[grep( paste0("(^|_)", gene1, '[._]') , names(motif), ignore.case = T) ]
plotM <- view_motifs(motifs = motifx, method ='PCC' , tryRC = T, min.overlap = 0.5, normalise.scores =T ) +
  theme_void() + theme1 + scale_fill_manual(values = c( "#50C878", "#FF7185", "#ffbf00", "#57ABFF") ) 


ggsave( paste0( gene1, "_motif.pdf"), plot = plotM, path = figDir, scale = 1.3, width = 8, 
        height = length(motifx), units = "cm", dpi = 150, limitsize = F )
```

## Implementation on your data

### using with Signac (snATAC-seq)

```r
library(universalmotif) # BiocManager::install("universalmotif")
# library("motifStack") # BiocManager::install('motifStack')
library('chromVAR') # BiocManager::install('chromVAR')
library(motifmatchr) # BiocManager::install('motifmatchr')
library('BSgenome.Mmusculus.UCSC.mm39') 
library(TFBSTools)
library(Signac) # 本身也是利用motifmatchr的功能进行扫描

motifTFBS <- convert_motifs( motifs = curatedMotif, class = 'TFBSTools-PFMatrix')
motifTFBS <- do.call(TFBSTools::PFMatrixList, motifTFBS)  
brain <- AddMotifs(object = brain, genome = BSgenome.Mmusculus.UCSC.mm39, pfm = motifTFBS )

brain <- RunChromVAR( object = brain, genome = BSgenome.Mmusculus.UCSC.mm39)
# 此过程中motif名字中的特殊字符 | _ 会被转换成- 会导致ATAC assay和chromvar assay中的motif名字不一致
DefaultAssay(brain) <- 'chromvar' 
```

### scan any genomic region

```r
library(universalmotif) # BiocManager::install("universalmotif")
# library("motifStack") # BiocManager::install('motifStack')
library('chromVAR') # BiocManager::install('chromVAR')
library(motifmatchr) # BiocManager::install('motifmatchr')
library('BSgenome.Mmusculus.UCSC.mm39') 
library(TFBSTools)

motifTFBS <- universalmotif::convert_motifs( motifs = curatedMotif, class = 'TFBSTools-PFMatrix')
motifTFBS <- do.call(TFBSTools::PFMatrixList, motifTFBS) 
motifMtx <- Signac::CreateMotifMatrix(features = Signac::StringToGRanges(regions = peak, sep = c("-", "-")), 
                                       pwm = motifTFBS, score = T, use.counts = F,
                                       genome = 'BSgenome.Mmusculus.UCSC.mm39', sep = c("-", "-") ) 
motifMtx <- as.matrix(motifMtx)
```

### with peak counts

|                      | Lhx2cKO_1 | Lhx2cKO_2 | WT_1 | WT_2 |
| -------------------- | --------- | --------- | :--- | ---- |
| chr1:3164708-3165648 | 59        | 64        | 117  | 108  |
| chr1:3189930-3190957 | 80        | 94        | 60   | 55   |
| chr1:3469987-3470570 | 19        | 30        | 55   | 63   |
| chr1:3740787-3741419 | 57        | 60        | 48   | 39   |

```r
library(universalmotif)
library(motifStack)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm39)
library(tidyverse)
library(motifmatchr)

curatedMotif <- readRDS(paste0(OneDrive, '/code/R-code/curatedMotif.rds') )
# counts <- counts[grep(pattern = '^chr', rownames(counts)),] 
# stander chromatin, error evolced from mismatch chr scoffold between data and BSgenome
counts <- countX[, c(mutantCol, controlCol)]
countRC <- 1e+07*t(t(counts)/colSums(counts)) %>% as.data.frame() # normalization sequencing depth

scanMotifs <- function(motifX = curatedMotif, peakCounts = countRC, 
                       genome = 'BSgenome.Mmusculus.UCSC.mm39', method = 'matches'  ) {
  motifTFBS <- convert_motifs( motifs = motifX, class = 'TFBSTools-PFMatrix')
  motifTFBS <- do.call(TFBSTools::PFMatrixList, motifTFBS)
  motifMtx <- matchMotifs(pwms = motifTFBS, genome = genome, 
                          subject = StringToGRanges(regions = rownames(peakCounts), sep = c(":", "-")), 
                          out = "scores" ) # bg = background
  
  if (method == 'matches') {
    bind <- motifmatchr::motifMatches(motifMtx) # 是否出现  
    bind2 <- as.integer(bind)
    dim(bind2) <- dim(bind)
    rownames(bind2) <- paste(motifMtx@rowRanges@seqnames, motifMtx@rowRanges@ranges, sep = ':' ) 
    colnames(bind2) <- motifMtx@colData$name
  }
  
  if (method == 'counts') {
    bind2 <- motifmatchr::motifCounts(motifMtx)  
    rownames(bind2) <- paste(motifMtx@rowRanges@seqnames, motifMtx@rowRanges@ranges, sep = ':' )
  }
  
  motifCounts <- data.frame(row.names = motifMtx@colData$name)
  for (i in colnames(peakCounts)) {
    motifCounts[[i]] <- colSums( bind2*peakCounts[[i]])
  }
  
  return(motifCounts)
}

Lhx2cKOmotif <- scanMotifs(motifX = curatedMotif, peakCounts = countRC, 
                           genome = 'BSgenome.Mmusculus.UCSC.mm39', method = 'matches' )
```

```r
projectX <- 'Lhx2cKO_ATAC'
outDir <- 'E:/Lhx2cKO/Lhx2cKO_ATAC'
figDir <- 'E:/Lhx2cKO/Lhx2cKO_ATAC/png'
dir.create(figDir,recursive = T)

library(edgeR)
library(tidyverse)
 
y <- DGEList(count= motifMtx, remove.zeros= T, group=groupX$group ) # 构建EdgeR对象
# y$samples$group <- rep(c('Lhx2cKO', 'control'), each =2)
y <- calcNormFactors(y)  # 计算标准化因子 默认使用TMM方法
y <- estimateDisp(y) # 同时计算common disperse 和 tagwise disperse
DiffExp <- exactTest(y, pair = c('control', 'Lhx2cKO')) # 对照在前
DiffExp$table$FDR <- p.adjust(DiffExp$table$PValue, method="fdr", n=length(DiffExp$table$PValue)) ## 校正p值
dataE <- DiffExp$table    ##导出全部结果
dataE$motif <- rownames(dataE)
dataE$logCPM <- NULL
dataE$'nlogQ' <- -log((dataE$FDR + 1.028776e-323), base = 10) 
dataE$change <- ifelse(dataE$FDR > 0.05, 'NS', ifelse(dataE$logFC <0, 'down', 'up'))
# dataE$log2FC <- log2(rowMeans( motifMtx[,1:2])  / rowMeans(motifMtx[,3:4]) )
dataE$diff <- rowMeans(motifMtx[,1:2]) - rowMeans(motifMtx[,3:4])
                        
lab1 <- labs(title ='E15Ctx ATAC, Lhx2cKO vs WT', 
             subtitle = paste0("total = ", nrow(dataE), " motifs", ", up = ", 
                               nrow(dataE[dataE$change == "up", ]), ", down = ", nrow(dataE[dataE$change == "down", ]))) 

theme1 <- theme(plot.title=element_text(face="bold", size= 20, hjust=0.5,vjust=0.5, angle=360,lineheight=113),
                plot.subtitle = element_text(face="bold", color="steelblue",size= 16, hjust=1), 
                legend.key= element_blank(), legend.key.size=unit(0.8,'cm'), 
                axis.title=element_text(face="bold", color="black",size=14),
                axis.title.y = element_text(angle = 90),
                axis.line=element_line(colour='black', linewidth = 0.3),
                panel.background=element_rect(fill = "transparent"), 
                plot.background = element_rect(fill = "transparent", color = 'transparent')  )

ggplot(data = dataE,  mapping = aes( x = diff, y = nlogQ, color = change, fill = change, shape=, size= nlogQ)) + 
  geom_point(shape=21, alpha = 0.5, stroke = 0, color = 'transparent' )  + lab1 + theme_void() + theme1  +
  scale_fill_manual(values = c( '#4169E1', 'gray50',  '#FF0000'), na.value = 'gray15' ) +
  geom_hline(yintercept = 0, linetype = 4 ) + geom_vline(xintercept = 0, linetype = 4 ) +
  geom_text(data = dataE[dataE$change != "NS",], aes(label = motif ), 
            show.legend = F, check_overlap = T, color = 'black', alpha = 1 )  +
  guides(color = guide_legend(ncol = 1, bycol=T, override.aes = list(size = 6, stroke = 2, alpha = 0.8)), 
         fill = guide_legend(ncol = 1, bycol=T, override.aes = list(size = 6, stroke = 0)),
         size = guide_legend(ncol = 1, bycol=T, override.aes = list(fill = 'black', alpha = 0.8)))

ggsave(paste0( projectX, '_DEmotif.png'), plot = last_plot(), path = figDir, scale = 1.2, width = 16, height =  12, units = "cm", dpi = 144, limitsize = T)
```

```
topA <- DEmotif[DEmotif$change=='up',] %>% top_n(n=50, nlogQ) %>% rownames()
motifA <- curatedMotif[topA] %>% convert_motifs(class = "motifStack-pfm")
for (i in names(motifA)) { motifA[[i]]@color <- c("#50C878", "#FF7185", "#ffbf00", "#57ABFF") }
color2 <- scales::seq_gradient_pal(RColorBrewer::brewer.pal(8, "Dark2"))(seq(0, 1, length.out = length(motifA)))
png(file = paste0( figDir, "/up", '_motif.png'), units="in", width=10, height=10, res=300) # 创建画布
motifStack(pfms = motifA, layout="radialPhylog", 
           circle = 0.95, cleaves = 0.8, clabel.leaves = 0.6, col.bg=color2, col.bg.alpha=0.2,
           col.leaves=color2, col.inner.label.circle=color2, inner.label.circle.width=0.0,  
           # col.outer.label.circle=color2, outer.label.circle.width=0.02, motifScale="logarithmic",
           circle.motif= 1.25, angle=180, font = "Arial" ) 
dev.off()
```



## other merging methods 

1. `motifStack::mergeMotifs` merge multiple motifs by calculate mean of each position 

   -  对齐    

     motif是矩阵，需要先进行矩阵对齐, 再相似性打分(distance)并聚类.

     motifstack::matalign(Matrix Aligner)矩阵对齐  由 Matalign-v4a修改而来， 每次对比2个位置特异的矩阵

     motifstack::DNAmotifAlignment() # align DNA motifs for plotting motifs stack

   -  计算相似性

     计算distance, 并进行hclust等聚类  motifstack::motifHclust

     **层次聚类分析 Hierarchical Cluster Analysis (HCA)**: 最开始由一个数据点作为集群，然后针对每个集群，基于同样的标准合并集群。一直进行着个机构，最后只剩下一个集群，就找到了集群的层次结构。



2. `universalmotif::merge_motifs` 

   从读入的第一个motif开始，每次merge一个motif 

   -  对齐     Aligns the motifs using **compare_motifs()**, then averages the motif PPMs. 
   -  计算相似性  **compare_motifs()**支持多种算法 method= PCC, EUCL, SW, KL, ALLR, BHAT, HELL, SEUCL, MAN, ALLR_LL, WEUCL, WPCC; 
   - 其中ALLR and ALLR_LL 不支持merge_motif,仅可作图使用

   

3. `universalmotif::merge_similar`

   对于给定的motif列表，merge_similar() 函数会先通过 compare_motifs() 识别相似的motif，再使用 merge_motifs() 进行合并。

   

4.  `TOM聚类`

   1. 相似性打分, 生成***拓扑异构矩阵（TOM）*** 

      Topological overlap（TO，拓扑重叠）# 

      邻接/相关性矩阵（Adjacency Matrix）每两个motif之间的相似性分数 homer -matrix

   2. 用TOM矩阵进行聚类，从而对基因分块，得到不同的模块（module）。

      用1减去TOM，从而得到对应的基因不相似性的矩阵（dissTOM）。

      dissTOM转换成距离矩阵，越相似的motif距离越小，越容易聚集成类

      ```r
      write_homer(univMotif,'JAScis_homer.txt')
      compareMotifs.pl JAScis_homer.txt -reduceThresh 0.5 -matrix mtx.txt -cpu 8
      
      mtx <- read.delim("F:/ATACsn/motif/mtx.txt", row.names=1, header=T)
      mtxSmi <- mtx[names(ic4), names(ic4)]
      gfaps <- CreateSeuratObject(counts = mtxSmi)
      gfaps@assays$RNA@scale.data <- as.matrix(mtxSmi)
      gfaps <- RunPCA(gfaps, features = rownames(gfaps))
      gfaps <- RunTSNE(gfaps, dims = 1:50, reduction = 'pca', check_duplicates = F ) 
      # gfaps <- RunUMAP(gfaps, dims = 1:50)
      gfaps <- FindNeighbors(gfaps)
      gfaps <- FindClusters(gfaps, res=0.5)
      
      # gfaps@meta.data <- metaMotif[rownames(gfaps@meta.data),]
      
      theme1 <- theme(axis.title = element_blank(),  axis.ticks = element_blank(), axis.text = element_blank())
      p2 <- DimPlot(object = gfaps, label = T, label.size = 5, pt.size = 1, cols = color$col27) + guides(color=guide_legend(ncol = 1, bycol=T, override.aes = list(size = 5))) + ggtitle('louvain cluster') & theme1
      gene1 <- 'Lhx2'
      DimPlot(gfaps, pt.size = 1, label = T, cells.highlight = metaMotif[metaMotif$Symbol==gene1,]$motif ) +ggtitle(gene1) & theme1 & NoLegend()
      ggsave( paste0( 'motif_', gene1, '.png'), plot = last_plot(), path = getwd(), scale = 1, width = 10, height = 10, units = "cm", dpi = 300, limitsize = T )
      
      gfaps$'tSEN_1' <-  gfaps@reductions$tsne@cell.embeddings[1]
      gfaps$'tSEN_2' <-  gfaps@reductions$tsne@cell.embeddings[2]
      ```

      ```r
      # 每个基因按类合并
      ## meta: Symbol geneName; Cluster, 不能只是数字; motif, 与motifList中的motif的名字保持一致;
      mergeMotifByGene <- function( geneList = NULL, motifMeta = meta, motifList = stacMotif) {
        tfMotif <- list()
        for (i in geneList) {
          meta <- motifMeta[motifMeta$Symbol == i,]
          N <- table(meta$cluster) %>% as.data.frame()
          rownames(N) <- N$Var1
          for (ii in N$Var1) {
            motifX <- meta[meta$cluster == ii,]$motif
            if(length(motifX) >1){
              tfMotif[[paste0(i,".", ii, 'N', N[ii,]$Freq ) ]] <- mergeMotifs(stacMotif[motifX])
            }
            else{tfMotif[[paste0(i,".", ii, 'N', N[ii,]$Freq ) ]] <- stacMotif[[motifX]] }
          }
        }
        return(tfMotif)
      }
      
      tfMotif <- mergeMotifByGene(geneList = c('Lhx2', 'Olig2'), motifMeta = meta, motifList = stacMotif )
      tfMotif <- filter_motifs(motifs = tfMotif, icscore = 5)
      ```

      







#### **stacMotif合并**

```R
TFname = 'Lhx2'
pfms = stacMotif[metaMotif[metaMotif$gene==TFname,]$motif ]

# 聚类
hc <- clusterMotifs(pfms)
## convert the hclust to phylog object

library(ade4)
phylog <- ade4::hclust2phylog(hc)
## reorder the pfms by the order of hclust

pfms <- pfms[names(phylog$leaves)]
## extract the motif signatures
motifSig <- motifSignature(pfms, phylog, cutoffPval=0.0001, min.freq=1)
sig <- signatures(motifSig)

pfmsAligned <- DNAmotifAlignment(pfms)
motifPiles(phylog=phylog, pfms=pfmsAligned, cleaves = 0.2, clabel.leaves = 3, r.tree = 0.1, motifScale="logarithmic", plotIndex=T )

plotMotifStackWithRadialPhylog(phylog=phylog, pfms=sig, 
  circle=1, cleaves = 0.5, clabel.leaves = 0.6, col.bg=color2, col.bg.alpha=0.3, col.leaves=color2,
  # col.inner.label.circle=gpCol,  inner.label.circle.width=0, 
  angle=358, circle.motif=1.8, motifScale="logarithmic")
```



# 
