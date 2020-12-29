library(tidyverse)
library(ShortRead)

r <- c('201012_M03249_0104_000000000-JBCL8', '201014_M03249_0105_000000000-J78D3', '201028_M03249_0110_000000000-G675M')
d <- lapply(r, function(x){
       o <- readFastq(paste0('/data/sequencing/Illumina-archive/', x, '/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz'))
       o <- trimTailw(o, 2, '?', 5)
       o[width(o) >= 60]
      })

d <- Reduce('append', d)
ids <- as.character(d@id)
d <- d@sread
names(d) <- ids

sort(table(as.character(subseq(d, 1, 8))), decreasing = TRUE)[1:5]

# LTR reads are expected to GGGTCTTT.
s <- as.character(subseq(d, 1, 8))
d <- d[s == 'GGGTCTTT']
d <- unique(d)
sort(table(as.character(subseq(d, 1, 8))), decreasing = TRUE)

cluster <- parallel::makeCluster(20)
if(! file.exists('alignments.rds')){
  alignments <- bind_rows(parallel::parLapply(cluster, split(d, ntile(1:length(d), 20)), function(x){
         library(ShortRead)
         source('lib.R')
         alignReads.BLAT(x, '/home/everett/data/sequenceDatabases/BLAT/hg38/hg38.2bit')
       }))
  parallel::stopCluster(cluster)
  saveRDS(alignments, 'alignments.rds')
} else {
  alignments <- readRDS('alignments.rds')
}

alignments <- 
  filter(alignments, 
         alignmentPercentID >= 95,
         blockCount  == 1,
         tNumInsert  <= 1, 
         qNumInsert  <= 1,
         tBaseInsert <= 2,
         qBaseInsert <= 2) 


# Only consider alignments with 50 more more letters and few unaligned letters at the end of the query.
alignments <-  filter(alignments, (qSize - qEnd) <= 3, matches >= 50)

alignments <- group_by(alignments, qName) %>%
              mutate(n = n()) %>%
              ungroup() %>%
              filter(n == 1)



View(subset(alignments, qStart == 21 & tStart == 49678803))


