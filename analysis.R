library(gt23)
library(tidyverse)
library(RMySQL)
library(RColorBrewer)
library(gtools) 
library(wordcloud) 
library(gridExtra)

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- unname(unlist(dbGetQuery(dbConn, 'select specimenAccNum from gtsp where Trial="MND-ADA-LTFU"')))

if(! file.exists('sites.rds')){
  intSites <- getDBgenomicFragments(samples, 'specimen_management', 'intsites_miseq') %>%
              stdIntSiteFragments() %>%
              collapseReplicatesCalcAbunds() %>%
              annotateIntSites()
  saveRDS(intSites, 'sites.rds')
} else {
  intSites <- readRDS('sites.rds')
}

d <- data.frame(intSites)
d$patient <- sub('^p', '', d$patient)
d$cellType <- gsub('^\\s+|\\s+$', '', toupper(d$cellType))

d <- dplyr::group_by(d, timePoint, cellType) %>%
     dplyr::mutate(readsPerSample = sum(reads)) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(timePoint, cellType, posid) %>%
     dplyr::mutate(readsRelAbund = (sum(reads) / readsPerSample[1])*100) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(cellType, timePoint, GTSP) %>%  
     dplyr::mutate(totalSampleFrags = sum(estAbund)) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(cellType, timePoint) %>% 
     dplyr::mutate(include = ifelse(totalSampleFrags == max(totalSampleFrags), 'yes', 'no')) %>%
     dplyr::ungroup()


# Create data frame needed to generate relative abundance plots.
numClones <- 10

d$posidLabel <- d$nearestFeature
#d$posidLabel <- paste0(d$nearestFeature, '\n', paste0(d$seqnames, d$strand, d$start))

abundantClones <- bind_rows(lapply(split(d, d$patient), function(x){
  
  # Adjust the number of clones to return based on the number of sites per cell type.
  if(nrow(x) < numClones) numClones <- nrow(x)
  
  # Sort nearest genes by abundance.
  x <- x[order(x$estAbund, decreasing = TRUE),]
  
  # Select clones to report.
  topClones <-  unique(x$posidLabel)[1:numClones]
  
  # For each time point, create a data frame for relative abundance plots
  bind_rows(lapply(split(x, x$timePoint), function(x2){
    
    lowAbundData <- dplyr::mutate(x2,
                                  patient = x$patient[1],
                                  posidLabel = 'LowAbund',
                                  totalCells = sum(estAbund),
                                  relAbund   = 100) %>%
      dplyr::slice(1) %>% 
      dplyr::select(patient, cellType, timePoint, posidLabel, totalCells, timePointDays, relAbund)
    
    x3 <- subset(x2, posidLabel %in% topClones)
    if(nrow(x3) == 0) return(lowAbundData)
    x3$totalCells <- sum(x2$estAbund)
    
    lowAbundData$relAbund <- 100 - sum(x3$relAbund)
    bind_rows(lowAbundData,  dplyr::select(x3, patient, cellType, timePoint, posidLabel, totalCells, timePointDays, relAbund))
  }))
}))



# Create named color vector for unique clones.
cloneColorsVector <- setNames(c('#eeeeee', colorRampPalette(brewer.pal(12, "Paired"))(n_distinct(abundantClones$posidLabel))),  c('LowAbund', unique(abundantClones$posidLabel)))

abundantClonesPlots <- lapply(split(abundantClones, abundantClones$patient), function(x){
  o <- subset(x, posidLabel != 'LowAbund')
  o <- o[order(o$relAbund, decreasing = TRUE),]
  
  x$posidLabel <- factor(x$posidLabel, levels = c('LowAbund', unique(o$posidLabel)))
  x <- x[order(x$timePointDays),]

  x$timePoint  <- factor(x$timePoint, levels = (unique(x$timePoint))) 
 
  totalCellLabel <- unname(unlist(lapply(split(x, x$timePoint), function(x) ppNum(x$totalCells[1]))))
  
  ggplot(x) +
    theme_bw() +
    scale_x_discrete(drop=FALSE) + 
    geom_bar(aes(timePoint, relAbund/100, fill=posidLabel), stat='identity', color = 'black', size = 0.20) + 
    scale_fill_manual(name = 'Clones', values = cloneColorsVector) +
    scale_shape_manual(values = c(16, 17, 15), drop = FALSE) +
    labs(x = '', y = '') +
    ggtitle(x$patient[1]) +
    guides(fill=guide_legend(title.position = "top", ncol=1, keyheight=0.35, default.unit="inch")) +
    scale_y_continuous(labels = scales::percent) + 
    annotate('text', x=1:length(totalCellLabel), y=1.04, label=totalCellLabel, size=2.7, angle=45, hjust=0.5) +
    theme(axis.text.x = element_text(angle = 315, hjust = 0),
          text = element_text(size=20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none")
})




pdf(height = 18, width = 4, file = 'relativeAbundancePlots_2x5.pdf')
do.call("grid.arrange", c(abundantClonesPlots, ncol=2))
dev.off()

pdf(height = 10, width = 15, file = 'relativeAbundancePlots_5x2.pdf')
do.call("grid.arrange", c(abundantClonesPlots, ncol=5))
dev.off()

pdf(height = 10, width = 10, file = 'relativeAbundancePlots_5x5.pdf')
do.call("grid.arrange", c(abundantClonesPlots, ncol=5))
dev.off()


maxWordCloudWords <- 100
x <- d
x <- x[order(x$estAbund, decreasing = TRUE),]


if(nrow(x) < maxWordCloudWords) maxWordCloudWords <- nrow(x)
x <- x[1:maxWordCloudWords,]

w <- setNames(x$estAbund, x$nearestFeature)
w <- w[! is.na(w)]

png(file = 'wordCloud_allSubjects_estAbund.png', res = 150)
wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWordCloudWords), 
          rot.per=0, max.words = maxWordCloudWords, scale = c(2.5, 0.15))
dev.off()

pdf(file = 'wordCloud_allSubjects_estAbund.pdf')
wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWordCloudWords), 
          rot.per=0, max.words = maxWordCloudWords, scale = c(2.5, 0.15))
dev.off()



x <- group_by(d, nearestFeature) %>%
     summarise(n = n_distinct(posid)) %>%
     ungroup() %>%
     arrange(desc(n)) %>%
     dplyr::slice(1:100)

w <- setNames(x$n, x$nearestFeature)

names(w) <- sub('\\-AS1', '', names(w))
     
png(file = 'wordCloud_allSubjects_nSites.png', res = 150)
wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWordCloudWords), 
          rot.per=0, max.words = maxWordCloudWords, scale = c(2.5, 0.25))
dev.off()

pdf(file = 'wordCloud_allSubjects_nSites.pdf')
wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWordCloudWords), 
          rot.per=0, max.words = maxWordCloudWords, scale = c(2.5, 0.25))
dev.off()




