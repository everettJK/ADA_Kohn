library(gt23)
library(tidyverse)
library(RMySQL)
library(RColorBrewer)
library(gtools) 
library(wordcloud) 
library(gridExtra)

# Create intSite data objects if it is not present on disk.

if(! file.exists('data/sites.rds') | ! file.exists('data/sites_reps.rds')){
  dbConn  <- dbConnect(MySQL(), group='specimen_management')
  samples <- unname(unlist(dbGetQuery(dbConn, 'select specimenAccNum from gtsp where Trial="MND-ADA-LTFU"')))
  
  intSites <- getDBgenomicFragments(samples, 'specimen_management', 'intsites_miseq') %>%
              stdIntSiteFragments() %>%
              collapseReplicatesCalcAbunds() %>%
              annotateIntSites()
  
  intSites_reps <- getDBgenomicFragments(samples, 'specimen_management', 'intsites_miseq') %>%
                   stdIntSiteFragments() %>%
                   calcSampleAbunds() %>%
                   annotateIntSites()
  
  
  saveRDS(intSites, 'data/sites.rds')
  saveRDS(intSites_reps, 'data/sites_reps.rds')
} else {
  intSites <- readRDS('data/sites.rds')
  intSites_reps <- readRDS('data/sites_reps.rds')
}

# Standardize patient and cellType identifiers.
d <- data.frame(intSites)
d$patient  <- sub('^p', '', d$patient)
d$cellType <- gsub('^\\s+|\\s+$', '', toupper(d$cellType))

# Expand intSite data object by adding reads per sample and fragments (cells) per sample.
d <- dplyr::group_by(d, timePoint, cellType) %>%
     dplyr::mutate(readsPerSample = sum(reads)) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(timePoint, cellType, posid) %>%
     dplyr::mutate(readsRelAbund = (sum(reads) / readsPerSample[1])*100) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(cellType, timePoint, GTSP) %>%  
     dplyr::mutate(totalSampleFrags = sum(estAbund)) %>%
     dplyr::ungroup() 


# MECOM analysis 
#--------------------------------------------------------------------------------------------------
bind_rows(lapply(split(intSites, paste(intSites$patient, intSites$timePoint)), function(x){
  tibble(patient = x$patient[1], timePoint = x$timePoint[1], 
         MECOM_sites = n_distinct(subset(x, nearestFeature == 'MECOM')$posid),
         replicates = n_distinct(subset(intSites_reps, patient == x$patient[1] & timePoint == x$timePoint[1])$sampleName))
  
})) %>% openxlsx::write.xlsx('output/MECOM_sites.xlsx')


MECOM_replicates <- bind_rows(lapply(split(intSites_reps, intSites_reps$patient), function(x){
  
  x <- group_by(data.frame(x), sampleName) %>% mutate(sampleCells = sum(estAbund)) %>% 
       ungroup() %>%
       filter(sampleCells >= 50)
  
  if(n_distinct(x$timePoint) != 2) return(tibble())
  
  a <- subset(x, timePointDays == min(x$timePointDays) & nearestFeature == 'MECOM')
  b <- subset(x, timePointDays == max(x$timePointDays) & nearestFeature == 'MECOM')
  if(nrow(a) == 0 | nrow(b) == 0) return(tibble())
  
  sharedSites <- base::intersect(a$posid, b$posid)
  if(length(sharedSites) == 0) return(tibble())
  
  # Site in later time point with greatest relative abundance -- may be more than one.
  bind_rows(lapply(sharedSites, function(site){
    a <- select(subset(a, posid == site), patient, sampleName, timePoint, timePointDays, posid, nearestFeatureDist, sampleCells, relAbund) 
    b <- select(subset(b, posid == site), patient, sampleName, timePoint, timePointDays, posid, nearestFeatureDist, sampleCells, relAbund) 
  
    bind_rows(a, b) %>% mutate(Wilcox_pVal = wilcox.test(a$relAbund, b$relAbund)$p.value) %>% select(-timePointDays)
  })) 
})) 

openxlsx::write.xlsx(MECOM_replicates, 'output/MECOM_replicates.xlsx')



r <- bind_rows(lapply(split(intSites, intSites$patient), function(x){
  x <- group_by(data.frame(x), GTSP) %>%
       mutate(cellsPerSample = sum(estAbund)) %>%
       ungroup() %>%
       filter(cellsPerSample >= 50)
  
  a <- subset(x, timePointDays == min(x$timePointDays) & nearestFeature == 'MECOM')
  b <- subset(x, timePointDays == max(x$timePointDays) & nearestFeature == 'MECOM')
  if(nrow(a) == 0 | nrow(b) == 0) return(tibble())
  
  sharedSites <- base::intersect(a$posid, b$posid)
  if(length(sharedSites) == 0) return(tibble())
  
  bind_rows(lapply(sharedSites, function(site){
    data.frame(patient = x$patient[1], site = site, timePointDays = c(a$timePointDays[1], b$timePointDays[1]),
               timePoint = c(a$timePoint[1], b$timePoint[1]), 
               relAbund = c(subset(a, posid == site)$relAbund, subset(b, posid == site)$relAbund))
    
  }))
}))
                             
r <- arrange(r, timePointDays)
r$timePoint <- factor(as.character(r$timePoint), levels = unique(r$timePoint))
r$g <- paste(r$patient, r$site)                             

colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(r$g))

MECOM_sample_relAbund_plot <- 
  ggplot(r, aes(timePoint, relAbund/100, group = g, color = g)) + 
  theme_bw() +
  scale_color_manual(name = 'Clones', values = colors) +
  geom_point(size = 3) + 
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Time point', y = 'Relative Abundance') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(MECOM_sample_relAbund_plot, file = 'output/MECOM_sample_relAbund_plot.pdf', units = 'in', width = 8)


# Fishers exact tests
#--------------------------------------------------------------------------------------------------
o <- bind_rows(lapply(split(intSites, intSites$patient), function(x){
  x$late <- ifelse(x$timePointDays == min(x$timePointDays), FALSE, TRUE)   
  data.frame(x)
}))

m <- matrix(c(n_distinct(subset(o, late == FALSE & nearestFeature != 'MECOM')$posid),
              n_distinct(subset(o, late == FALSE & nearestFeature == 'MECOM')$posid),
              n_distinct(subset(o, late == TRUE  &  nearestFeature != 'MECOM')$posid),
              n_distinct(subset(o, late == TRUE  &  nearestFeature == 'MECOM')$posid)),
            byrow = FALSE, ncol = 2, dimnames = list(c('Not MECOM', 'MECOM'), c('Early', 'Late')))

pEarly <- (m[2,1] / sum(m[,1]))*100
pLate  <- (m[2,2] / sum(m[,2]))*100

fisher.test(m, alternative = 'two.sided')$p.value
fisher.test(m, alternative = 'greater')$p.value



# Relative abundance analyses 
#--------------------------------------------------------------------------------------------------

# Create data frame needed to generate relative abundance plots.
numClones <- 10


# Add nearest feature flags.
d <- d %>%
  mutate(labeledNearestFeature = paste0(nearestFeature, ' ')) %>% 
  mutate(labeledNearestFeature = ifelse(inFeature, paste0(labeledNearestFeature, '*'), labeledNearestFeature)) 

if('nearestOncoFeatureDist' %in% names(d))
  d <- mutate(d, labeledNearestFeature = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(labeledNearestFeature, '~'), labeledNearestFeature))

if('nearestlymphomaFeatureDist' %in% names(d))
  d <- mutate(d, labeledNearestFeature = ifelse(abs(nearestlymphomaFeatureDist) <= 50000, paste0(labeledNearestFeature, '!'), labeledNearestFeature)) 

d$posidLabel <- paste0(d$labeledNearestFeature, '\n', d$posid)

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

pdf(height = 20, width = 8, file = 'output/relativeAbundancePlots_2x5.pdf')
do.call("grid.arrange", c(abundantClonesPlots, ncol=2))
dev.off()

pdf(height = 10, width = 15, file = 'output/relativeAbundancePlots_5x2.pdf')
do.call("grid.arrange", c(abundantClonesPlots, ncol=5))
dev.off()

pdf(height = 10, width = 12, file = 'output/relativeAbundancePlots_5x5.pdf')
do.call("grid.arrange", c(abundantClonesPlots, ncol=5))
dev.off()


# Wordclouds
#--------------------------------------------------------------------------------------------------

nearestGeneSitesWordClounds <- function(x, maxWordCloudWords, outputFile){
  x <- group_by(x, labeledNearestFeature) %>%
    summarise(n = n_distinct(posid)) %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    dplyr::slice(1:maxWordCloudWords)
  
  w <- setNames(x$n, x$labeledNearestFeature)
  names(w) <- sub('\\-AS1', '', names(w))
  
  pdf(height = 4, width = 4, file = file.path('output', outputFile))
  wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWordCloudWords), 
            rot.per=0, max.words = maxWordCloudWords, scale = c(2.5, 0.25), min.freq = 1)
  dev.off()
}


nearestGeneSitesWordClounds(d, 100, 'wordCloud_allSubjects_nSites.pdf')

invisible(lapply(split(d, d$patient), function(x){
  nearestGeneSitesWordClounds(x, 100, paste0('wordCloud_subject_', x$patient[1], '_nSites.pdf'))
}))


