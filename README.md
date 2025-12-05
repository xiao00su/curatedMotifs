



# curatedMotifs



- 3230 High-quality Non-Redundant motifs for 1410 mouse/human genes.
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

      

## Installation

You can install the development version of curatedMotifs like so:

``` r
remotes::install_github('xiao00su/curatedMotifs')
```



## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(curatedMotifs)
data("curated_motifs")
data("curated_motifs_meta")
data("raw_motifs")
data("raw_motifs_meta")

plot_rawVScurated( gene = 'Lhx2')
```



### signac

```r
library(universalmotif) # BiocManager::install("universalmotif")
library(curatedMotifs) 
library('BSgenome.Mmusculus.UCSC.mm39') # your BSgenome
library(TFBSTools)
library(Signac)  

data("curated_motifs")
motifTFBS <- convert_motifs( motifs = curated_motifs, class = 'TFBSTools-PFMatrix')
motifTFBS <- do.call(TFBSTools::PFMatrixList, motifTFBS)  

brain <- AddMotifs(object = brain, genome = BSgenome.Mmusculus.UCSC.mm39, pfm = motifTFBS )
brain <- RunChromVAR( object = brain, genome = BSgenome.Mmusculus.UCSC.mm39)
```



### scan any genomic region

```r
library(curatedMotifs)
library(universalmotif) # BiocManager::install("universalmotif")
library(Signac)  
library('BSgenome.Mmusculus.UCSC.mm39') 
library(TFBSTools)
library(GenomicRanges)

data("curated_motifs")
motifTFBS <- universalmotif::convert_motifs( motifs = curated_motifs, class = 'TFBSTools-PFMatrix')
motifTFBS <- do.call(TFBSTools::PFMatrixList, motifTFBS) 

data("peaks")
colnames(peaks)[1:3] <- c('chr', 'start', 'end')
GR <- GenomicRanges::makeGRangesFromDataFrame(df = peaks )
```

```r
# motif  scan using signac
motifMtx <- Signac::CreateMotifMatrix(features = GR, pwm = motifTFBS, score = T, use.counts = F,
                                       genome = 'BSgenome.Mmusculus.UCSC.mm39', sep = c("-", "-") ) 
motifMtx <- as.matrix(motifMtx)
```

```R
# motif scan using motifmatchr
motifMtx <- motifmatchr::matchMotifs(subject = GR, pwms = motifTFBS, genome = genome, 
                                     sep = c(":", "-")),  out = "scores" ) 
bind <- motifmatchr::motifMatches(motifMtx)
rownames(bind) <- paste(motifMtx@rowRanges@seqnames, motifMtx@rowRanges@ranges, sep = ':' ) 
```



### MEME motif

To obtain the curated DNA binding motif dataset, simply download the file `curatedMotif_meme.txt`—no package installation is required.

```bash
curated_motifs=PATH/curatedMotif_meme.txt
centrimo XXXX_summits.fa  ${curated_motifs} --oc outDir 
```



## Versus other databases

Compared to conventional motif databases (chromVar, HocoMoco13, JASPAR2024, and Cisbp3), curatedMotifs provides a comprehensive, non‑redundant, and manually curated motif collection that combines wide transcription‑factor coverage with high‑quality motifs, offering researchers a clean, analysis‑ready resource that minimizes noise and maximizes specificity in transcription‑factor binding‑site analyses.

<img src="https://typora-1329573677.cos.ap-shanghai.myqcloud.com/img/Sox2_motif2.png" alt="Sox2_motif2" style="zoom:25%;" />

