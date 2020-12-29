
tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }


alignReads.BLAT <- function(x, db, outputDir = '.', command.blat = ' ~/ext/blat'){
  f <- tmpFile()
  writeFasta(x, file = file.path(outputDir, 'tmp', paste0(f, '.fasta')))
  
  comm <- paste0(command.blat, ' ', db, ' ', 
                 file.path(outputDir, 'tmp', paste0(f, '.fasta')), ' ', 
                 file.path(outputDir, 'tmp', paste0(f, '.psl')),  
                 ' -tileSize=11 -stepSize=9 -minIdentity=90 -out=psl -noHead')
  system(comm)
  file.remove(file.path(outputDir, 'tmp', paste0(f, '.fasta')))
  
  if(file.exists(file.path(outputDir, 'tmp', paste0(f, '.psl')))){
    b <- parseBLAToutput(file.path(outputDir, 'tmp', paste0(f, '.psl')))
    file.remove(file.path(outputDir, 'tmp', paste0(f, '.psl')))
    b
  } else {
    return(tibble())
  }
}

# BLAT output parser.
parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = readr::cols())
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID       <- (b$matches/b$qSize)*100
  b$tAlignmentWidth      <- (b$tEnd - b$tStart) + 1
  b$queryWidth           <- (b$qEnd - b$qStart) + 1
  b$alignmentPercentID   <- (b$matches/b$tAlignmentWidth)*100
  b$percentQueryCoverage <- (b$queryWidth/b$qSize)*100
  b$qStarts              <- as.character(b$qStarts)
  b$tStarts              <- as.character(b$tStarts)
  b
}