##==========================================================================================================================================================
##===========
## FUNCTIONS
##===========

require(tools)

##
removeDuplicated <- function(referenceDataDF) {
  ## Verify if searchTerms are uniques
  if (TRUE %in% duplicated(referenceDataDF$searchTerm)) {
    duplicatedsIndexes <- which(duplicated(referenceDataDF$searchTerm))
    duplicatedsSearchTerms <- referenceDataDF$searchTerm[duplicatedsIndexes]
    for (duplicatedSearchTerm in duplicatedsSearchTerms) {
      duplicatedRelativeAbundanceSum <- referenceDataDF %>%
        filter(searchTerm == duplicatedSearchTerm) %$%
        sum(relativeAbundance)
      duplicatedSearchTermIndexes <- which(referenceDataDF$searchTerm == duplicatedSearchTerm)
      referenceDataDF$relativeAbundance[duplicatedSearchTermIndexes] <- duplicatedRelativeAbundanceSum
    }
    
    referenceDataDF <- referenceDataDF[-duplicatedsIndexes,]
  }
  
  return(referenceDataDF)
}

##
checkComposition <- function(referenceDataDF) {
  hasSearchTerm <- "searchTerm" %in% names(referenceDataDF)
  if (hasSearchTerm == FALSE) {
    stop("Ops! Your dataset should have a searchTerm field.")
  }
  
  ## Verify if searchTerms are uniques
  if (TRUE %in% duplicated(referenceDataDF$searchTerm)) {
    duplicatedsIndexes <- which(duplicated(referenceDataDF$searchTerm))
    duplicatedsSearchTerms <- referenceDataDF$searchTerm[duplicatedsIndexes]
    stop("Ops! Search term is duplicated." %>% paste("(", duplicatedsSearchTerms ,")"))
  }
  
  hasRelativeAbundance <- "relativeAbundance" %in% names(referenceDataDF)
  if (hasRelativeAbundance == FALSE) {
    ## add relative abundances
    
    randomRelativeAbundances <- NA
    maxRelativeAbundance <- 100
    N <- nrow(referenceDataDF)
    while(maxRelativeAbundance > 40) {
      randomAbundances <- poweRlaw::rpldis(N, 1, 1.5)
      randomAbundancesSum <- randomAbundances %>% sum
      randomRelativeAbundances <- randomAbundances / randomAbundancesSum * 100
      maxRelativeAbundance <- max(randomRelativeAbundances)
    }
    referenceDataDF$relativeAbundance <- randomRelativeAbundances %>% sort
  }
  
  if (hasRelativeAbundance == TRUE &&
      referenceDataDF$relativeAbundance %>% sum != 100) {
    stop("Ops! Your relative abundances sum should be 100. It was: " %>% paste0(referenceDataDF$relativeAbundance %>% sum))
  }
  
  referenceDataDF <- referenceDataDF %>% dplyr::select(searchTerm, relativeAbundance)
  referenceDataDF$downloadStatus <- F
  return(referenceDataDF)
}


## Function to retrieve assembly information from NCBI
retrieveAssemblySummary <- function(term = NA) {
  assembly.search.result <- rentrez::entrez_search("assembly", term)
  if (assembly.search.result$count == 0) {
    return(NA)
  }
  assembly.summary.result <- rentrez::entrez_summary("assembly", assembly.search.result$ids)
  return(assembly.summary.result)
}

retrieveAssemblySummaryThroughGenomeSummaryAssemblyName <- function(term = NA) {
  genome.search.result <- rentrez::entrez_search("genome", term)
  if (genome.search.result$count == 0) {
    return(NA)
  }
  genome.summary.result <- rentrez::entrez_summary("genome", genome.search.result$ids)
  if (genome.summary.result$assembly_name == "") {
    return(NA)
  }
  assembly.summary.result <- retrieveAssemblySummary(genome.summary.result$assembly_name)
  return(assembly.summary.result)
}

retrieveAssemblySummaryThroughGenomeSummaryAssemblyId <- function(term = NA) {
  genome.search.result <- rentrez::entrez_search("genome", term)
  if (genome.search.result$count == 0) {
    return(NA)
  }
  genome.summary.result <- rentrez::entrez_summary("genome", genome.search.result$ids)
  if (genome.summary.result$assemblyid == "") {
    return(NA)
  }
  assembly.summary.result <- rentrez::entrez_summary("assembly", genome.summary.result$assemblyid)
  return(assembly.summary.result)
}

retrieveAssemblySummaryThroughGenomeSummaryOrganismName <- function(term = NA) {
  genome.search.result <- rentrez::entrez_search("genome", term)
  if (genome.search.result$count == 0) {
    return(NA)
  }
  genome.summary.result <- rentrez::entrez_summary("genome", genome.search.result$ids)
  if (genome.summary.result$organism_name == "") {
    return(NA)
  }
  assembly.summary.result <- retrieveAssemblySummary(genome.summary.result$organism_name)
  return(assembly.summary.result)
}


##
getGcaGcf <- function(gcaGcfPath = NA) {
  gca_gcf.strsplit <- strsplit(gcaGcfPath %>% as.character, "/")
  gca_gcf.length <- gca_gcf.strsplit[[1]] %>% length
  gca_gcf <- gca_gcf.strsplit[[1]][gca_gcf.length]
  return(gca_gcf)
}
## getGcaGcf("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/171/635/GCF_000171635.1_ASM17163v1")

##
getGenomeUrlInAssemblyOrRefseq <- function(assemblyOrRefseqFtpPath) {
  gcaGcf <- getGcaGcf(assemblyOrRefseqFtpPath)
  fastaUrl <- assemblyOrRefseqFtpPath %>%
    paste0("/") %>%
    paste0(gcaGcf) %>%
    paste0("_genomic.fna.gz")
  return(fastaUrl)
}
## getGenomeUrlInAssemblyOrRefseq("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/171/635/GCF_000171635.1_ASM17163v1")

##
getCdssUrlInAssemblyOrRefseq <- function(assemblyOrRefseqFtpPath) {
  gcaGcf <- getGcaGcf(assemblyOrRefseqFtpPath)
  fastaUrl <- assemblyOrRefseqFtpPath %>%
    paste0("/") %>%
    paste0(gcaGcf) %>%
    paste0("_cds_from_genomic.fna.gz")
  return(fastaUrl)
}
## getCdssUrlInAssemblyOrRefseq("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/171/635/GCF_000171635.1_ASM17163v1")


##
downloadGenome <- function(destinationFolder = NA, genomeIdentifier = NA, assemblyOrRefseqFtpPath = NA) {
  ## Create destination folder if not exists
  dir.create(destinationFolder, showWarnings = F, recursive = T)
  
  ## Get file url in server
  fastaUrl <- getGenomeUrlInAssemblyOrRefseq(assemblyOrRefseqFtpPath)
  
  ## Get file destination in local
  destinationFile <- destinationFolder %>% paste0("/", genomeIdentifier %>% paste0(".fna.gz"))
  
  ## Download
  download.file(fastaUrl, destinationFile, method = "wget", quiet = T)
  
  ## Extract
  system("gunzip -k -f " %>% paste0(destinationFile))
  
  ## Apply a new unique identifier for all chromosomes and plasmid sequences
  uncompressedDestinationFile <- destinationFolder %>% paste0("/", genomeIdentifier %>% paste0(".fna"))
  system("sed -i 's/>/>REF /' FNA_FILENAME" %>%
           str_replace("REF", genomeIdentifier) %>%
           str_replace("FNA_FILENAME", uncompressedDestinationFile))
}

##
checkGenome <- function(destinationFolder = NA, genomeIdentifier = NA, assemblyOrRefseqFtpPath = NA) {
  ## Get file destination in local
  destinationFile <- destinationFolder %>% paste0("/", genomeIdentifier %>% paste0(".fna.gz"))
  
  ## Compute md5
  destinationFileMd5 <- md5sum(destinationFile) %>% as.vector
  
  ## Download files md5
  filesMd5DownloadUrl <- assemblyOrRefseqFtpPath %>% paste0("/md5checksums.txt")
  filesMd5Destination <- destinationFolder %>% paste0("/", genomeIdentifier ,".md5checksums.txt")
  download.file(filesMd5DownloadUrl, filesMd5Destination, method = "wget", quiet = T)
  
  ## Load md5s
  genomeMd5ChecksumsDF <- read.delim(filesMd5Destination, header = F, sep = " ")
  
  ## Check md5
  isMd5Ok <- destinationFileMd5 %in% genomeMd5ChecksumsDF$V1
  
  ## And return check result
  return(isMd5Ok)
}


##
downloadCdss <- function(destinationFolder = NA, genomeIdentifier = NA, assemblyOrRefseqFtpPath = NA) {
  ## Create destination folder if not exists
  dir.create(destinationFolder, showWarnings = F, recursive = T)
  
  ## Get file url in server
  fastaUrl <- getCdssUrlInAssemblyOrRefseq(assemblyOrRefseqFtpPath)
  
  ## Get file destination in local
  destinationFile <- destinationFolder %>% paste0("/", genomeIdentifier %>% paste0(".fna.gz"))
  
  ## Download
  download.file(fastaUrl, destinationFile, method = "wget", quiet = T)
  
  ## Extract
  system("gunzip -k -f " %>% paste0(destinationFile))
  
  ## Apply a new unique identifier for all chromosomes and plasmid sequences
  uncompressedDestinationFile <- destinationFolder %>% paste0("/", genomeIdentifier %>% paste0(".fna"))
  system("sed -i 's/>/>REF /' FNA_FILENAME" %>%
           str_replace("REF", genomeIdentifier) %>%
           str_replace("FNA_FILENAME", uncompressedDestinationFile))
}

##
checkCdss <- function(destinationFolder = NA, genomeIdentifier = NA, assemblyOrRefseqFtpPath = NA) {
  ## Get file destination in local
  destinationFile <- destinationFolder %>% paste0("/", genomeIdentifier %>% paste0(".fna.gz"))
  
  ## Compute md5
  destinationFileMd5 <- md5sum(destinationFile) %>% as.vector
  
  ## Download files md5
  filesMd5DownloadUrl <- assemblyOrRefseqFtpPath %>% paste0("/md5checksums.txt")
  filesMd5Destination <- destinationFolder %>% paste0("/", genomeIdentifier ,".md5checksums.txt")
  download.file(filesMd5DownloadUrl, filesMd5Destination, method = "wget", quiet = T)
  
  ## Load md5s
  genomeMd5ChecksumsDF <- read.delim(filesMd5Destination, header = F, sep = " ")
  
  ## Check md5
  isMd5Ok <- destinationFileMd5 %in% genomeMd5ChecksumsDF$V1
  
  ## And return check result
  return(isMd5Ok)
}

##
createGenomesReferenceFileForGrinder <- function(genomesFolder, outputFolder) {
  genomes.cat <- "cat GENOMES_FOLDER/*.fna > OUTPUT_FOLDER/genomes.fna" %>%
    str_replace("GENOMES_FOLDER", genomesFolder) %>%
    str_replace("OUTPUT_FOLDER", outputFolder)
  system(genomes.cat)
}

##
createCdssReferenceFileForGrinder <- function(cdssFolder, outputFolder) {
  cdss.cat <- "cat CDSS_FOLDER/*.fna > OUTPUT_FOLDER/cdss.fna" %>%
    str_replace("CDSS_FOLDER", cdssFolder) %>%
    str_replace("OUTPUT_FOLDER", outputFolder)
  system(cdss.cat)
}

##
#V1
# getSearchSummary <- function(searchTerm = NA) {
#   searchSummary <- NA
# 
#   searchSummary <- retrieveAssemblySummary(searchTerm)
# 
#   if (searchSummary %>% is.na) {
#     searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryAssemblyId(searchTerm)
#   }
# 
#   if (searchSummary %>% is.na) {
#     searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryOrganismName(searchTerm)
#   }
# 
#   if (searchSummary %>% is.na) {
#     stop("Ops! We did not find the NCBI's Assembly information for search term: " %>% paste0(searchTerm))
#   }
# 
#   return(searchSummary)
# }

# V2
getSearchSummary <- function(searchTerm = NA) {
  searchSummary <- NA
  
  searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryAssemblyName(searchTerm)
  
  if (searchSummary %>% is.na) {
    searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryAssemblyId(searchTerm)
  }
  
  if (searchSummary %>% is.na) {
    searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryOrganismName(searchTerm)
  }
  
  if (searchSummary %>% is.na) {
    searchSummary <- retrieveAssemblySummary(searchTerm)
  }
  
  if (searchSummary %>% is.na) {
    stop("Ops! We did not find the NCBI's Assembly information for search term: " %>% paste0(searchTerm))
  }
  
  return(searchSummary)
}

## V3
getSearchSummary <- function(searchTerm = NA) {
  searchSummary <- NA
  
  searchSummary <- retrieveAssemblySummary(searchTerm)
  
  if (searchSummary %>% is.na) {
    searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryAssemblyName(searchTerm)
  }
  
  if (searchSummary %>% is.na) {
    searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryAssemblyId(searchTerm)
  }
  
  if (searchSummary %>% is.na) {
    searchSummary <- retrieveAssemblySummaryThroughGenomeSummaryOrganismName(searchTerm)
  }
  
  if (searchSummary %>% is.na) {
    stop("Ops! We did not find the NCBI's Assembly information for search term: " %>% paste0(searchTerm))
  }
  
  return(searchSummary)
}
## txid1904.searchSummary <- getSearchSummary("txid1515")
## noSenseSummary <- getSearchSummary("txidnoSenseString")
## rm(txid1904.searchSummary)
## rm(noSenseSummary)

##
isSingleSummary <- function(searchSummary) {
  hasUidField <- "uid" %in% (searchSummary %>% names)
  return(hasUidField)
}

##
exportReferenceMetada <- function(destinationFolder, referenceDataDF) {
  referenceDataDF <- referenceDataDF %>% dplyr::select(searchTerm, relativeAbundance,
                                                taxid, organism, speciestaxid, speciesname, ftppath_genbank, ftppath_refseq, 
                                                usedFtpPathForGenome, usedFtpPathForCdss,
                                                genomeCheck, cdsCheck, downloadStatus, prodigalCds,
                                                primerCompatibleGenome, relativeAbundanceForAmplicons)
  write.table(referenceDataDF, destinationFolder %>% paste0("/referenceMetadata.tsv"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}

##
getListOfCompatibleGenomesForPrimers <- function(projectFolder, primersFile, ampliconSearchScriptFile) {
  ampliconSearch.command <- "perl AMPLICON_SEARCH_PL --fasta PROJECT_FOLDER/genomes.fna --primers PRIMERS_FILE > PROJECT_FOLDER/amplicons.txt" %>%
    str_replace_all("AMPLICON_SEARCH_PL", ampliconSearchScriptFile) %>%
    str_replace_all("PROJECT_FOLDER", projectFolder) %>%
    str_replace_all("PRIMERS_FILE", primersFile)
  system(ampliconSearch.command)
  ampliconsFile <- "PROJECT_FOLDER/amplicons.txt" %>% str_replace_all("PROJECT_FOLDER", projectFolder)
  compatibleGenomesList <- read.table(ampliconsFile, header = F)[,1]
  return(compatibleGenomesList)
}

##
addAdjustedRelativeAbundanceForAmplicons <- function(referenceDataDF) {
  ## Add a field for relative abundance for amplicons
  referenceDataDF$relativeAbundanceForAmplicons <- referenceDataDF$relativeAbundance
  
  ## Turn into zero the relative abundances of the genomes which had no compatibility with the primers
  toTurnIntoZero <- which(referenceDataDF$primerCompatibleGenome == F)
  referenceDataDF$relativeAbundanceForAmplicons[toTurnIntoZero] <- 0
  
  ## Calculate a normalize factor
  normalizationFactor <- 100/sum(referenceDataDF$relativeAbundanceForAmplicons)
  referenceDataDF$relativeAbundanceForAmplicons <- referenceDataDF$relativeAbundanceForAmplicons * normalizationFactor
  
  return(referenceDataDF)
}

##
addAmpliconsInformation <- function(projectFolder, primersFile, ampliconsSearchScriptFile, referenceDataDF) {
  genomesList <- getListOfCompatibleGenomesForPrimers(projectFolder, primersFile, ampliconsSearchScriptFile)
  referenceDataDF$primerCompatibleGenome <- referenceDataDF$searchTerm %in% genomesList
  referenceDataDF <- addAdjustedRelativeAbundanceForAmplicons(referenceDataDF)
  file.copy(primersFile, projectFolder %>% paste0("/primers.fna"))
  return(referenceDataDF)
}


##
checkSequencesFiles <- function(
  sequencesFolder,
  referenceDataDF) {
  
  filesSizes <- sequencesFolder %>% paste0("/") %>% paste0(composition.df$searchTerm) %>% paste0(".fna") %>% file.size
  
  filesSizesChecks <- !(filesSizes == 0 | filesSizes %>% is.na)
  
  return(filesSizesChecks)
}

##
exportRanks <- function(projectFolder, referenceDataDF) {
  write.table(referenceDataDF %>% dplyr::filter(genomeCheck) %>% dplyr::select(searchTerm, relativeAbundance),
              projectFolder %>% paste0("/sh-mt-ranks.txt"),
              col.names = F, row.names = F, quote = F, sep = "\t")
  
  write.table(referenceDataDF %>% dplyr::filter(primerCompatibleGenome) %>% dplyr::select(searchTerm, relativeAbundanceForAmplicons),
              projectFolder %>% paste0("/amp-ranks.txt"),
              col.names = F, row.names = F, quote = F, sep = "\t")
}

##
# getGenomesAndCdss <- function(PROJECT_FOLDER, GENOMES_FOLDER, CDSS_FOLDER, composition.df) {
#   composition.df$taxid <- NA
#   composition.df$organism <- NA
#   composition.df$speciestaxid <- NA
#   composition.df$speciesname <- NA
#   composition.df$ftppath_refseq <- NA
#   composition.df$ftppath_genbank <- NA
#   composition.df$usedFtpPathForGenome <- NA
#   composition.df$usedFtpPathForCdss <- NA
#   for (i in 1:nrow(composition.df)) {
#     ## Get serach term
#     searchTerm <- composition.df$searchTerm[i]
#     
#     ## Get summaries by search term
#     summaryResults <- getSearchSummary(searchTerm)
#     if (isSingleSummary(summaryResults)) {
#       summaryResults <- list(summaryResults)
#     }
#     
#     ## Get Genomes and Cdss
#     usedFtpPathForGenome <- NA
#     usedFtpPathForCdss <- NA
#     currentSummaryResultIndex <- 0
#     summaryResultsNumber <- summaryResults %>% length
#     areGenomesAndCdsOk <- FALSE
#     while (areGenomesAndCdsOk == FALSE &&
#            currentSummaryResultIndex < summaryResultsNumber) {
#       ## Get a summary result
#       currentSummaryResultIndex <- currentSummaryResultIndex + 1
#       summaryResult <- summaryResults[[currentSummaryResultIndex]]
#       
#       ## Try Genbank First (GENOME)
#       isGenomeOk <- FALSE
#       usedFtpPathForGenome <- summaryResult$ftppath_genbank
#       if ((usedFtpPathForGenome %>% is.na == FALSE) &&
#           (usedFtpPathForGenome == "") == FALSE) {
#         
#         try(downloadGenome(GENOMES_FOLDER, searchTerm, usedFtpPathForGenome))
#         isGenomeOk <- checkGenome(GENOMES_FOLDER, searchTerm, usedFtpPathForGenome)
#       }
#       
#       ## If not ok, try Refseq (GENOME)
#       if (isGenomeOk == FALSE) {
#         usedFtpPathForGenome <- summaryResult$ftppath_refseq
#         if ((usedFtpPathForGenome %>% is.na == FALSE) &&
#             (usedFtpPathForGenome == "") == FALSE) {
#           try(downloadGenome(GENOMES_FOLDER, searchTerm, usedFtpPathForGenome))
#           isGenomeOk <- checkGenome(GENOMES_FOLDER, searchTerm, usedFtpPathForGenome)
#         }
#       }
#       
#       ## Try Genbank First (CDSs)
#       isCdssOk <- FALSE
#       usedFtpPathForCdss <- summaryResult$ftppath_genbank
#       if ((usedFtpPathForCdss %>% is.na == FALSE) &&
#           (usedFtpPathForCdss == "") == FALSE) {
#         try(downloadCdss(CDSS_FOLDER, searchTerm, usedFtpPathForCdss))
#         isCdssOk <- checkCdss(CDSS_FOLDER, searchTerm, usedFtpPathForCdss)
#       }
#       
#       ## If not ok, try Refseq (CDSs)
#       if (isCdssOk == FALSE) {
#         usedFtpPathForCdss <- summaryResult$ftppath_refseq
#         if ((usedFtpPathForCdss %>% is.na == FALSE) &&
#             (usedFtpPathForCdss == "") == FALSE) {
#           try(downloadCdss(CDSS_FOLDER, searchTerm, usedFtpPathForCdss))
#           isCdssOk <- checkCdss(CDSS_FOLDER, searchTerm, usedFtpPathForCdss)
#         }
#       }
#       
#       ## Verification flag, TRUE if both Genome and CDSs is OK
#       areGenomesAndCdsOk <- isGenomeOk && isCdssOk
#     }
#     
#     if (areGenomesAndCdsOk == FALSE) {
#       stop("Ops! We could not found genomes and CDSs for search term: " %>% paste0(searchTerm))
#     }
#     
#     ## Keep information
#     composition.df$usedFtpPathForGenome[i] <- usedFtpPathForGenome
#     composition.df$usedFtpPathForCdss[i] <- usedFtpPathForCdss
#     
#     summaryResult <- summaryResults[[currentSummaryResultIndex]]
#     composition.df$taxid[i] <- summaryResult$taxid
#     composition.df$organism[i] <- summaryResult$organism
#     composition.df$speciestaxid[i] <- summaryResult$speciestaxid
#     composition.df$speciesname[i] <- summaryResult$speciesname
#     composition.df$ftppath_refseq[i] <- summaryResult$ftppath_refseq
#     composition.df$ftppath_genbank[i] <- summaryResult$ftppath_genbank
#   }
#   
#   return(composition.df)
# }

getFtpPath <- function(summaryResult, ncbiDatabase = "genbank") {
  ftpPath <- ""
  
  if (ncbiDatabase == "genbank") {
    ftpPath <- summaryResult$ftppath_genbank
  } else {
    ftpPath <- summaryResult$ftppath_refseq
  }
  
  if (ftpPath == "") {
    gcaGcf <- ""
    if (ncbiDatabase == "genbank") {
      gcaGcf <- summaryResult$synonym$genbank
    } else {
      gcaGcf <- summaryResult$synonym$refseq
    }
    
    filteredGcaGcf <- gcaGcf %>% str_replace_all("\\_", "") %>% str_replace_all("\\.1", "")
    gcaCgf_1 <- filteredGcaGcf %>% str_sub(1,3)
    gcaCgf_2 <- filteredGcaGcf %>% str_sub(4,6)
    gcaCgf_3 <- filteredGcaGcf %>% str_sub(7,9)
    gcaCgf_4 <- filteredGcaGcf %>% str_sub(10,12)
    
    ftpPath <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all" %>%
      paste0("/", gcaCgf_1) %>%
      paste0("/", gcaCgf_2) %>%
      paste0("/", gcaCgf_3) %>%
      paste0("/", gcaCgf_4) %>%
      paste0("/", gcaGcf) %>%
      paste0("_", summaryResult$assemblyname)
  }
  
  return(ftpPath)
}

getGenomesAndCdss <- function(genomesFolder, cdssFolder, referenceDataDF, checkDownloadStatus = T) {
  if (checkDownloadStatus == FALSE) {
    referenceDataDF$assemblyname <- NA
    referenceDataDF$taxid <- NA
    referenceDataDF$organism <- NA
    referenceDataDF$speciestaxid <- NA
    referenceDataDF$speciesname <- NA
    referenceDataDF$ftppath_refseq <- NA
    referenceDataDF$ftppath_genbank <- NA
    referenceDataDF$usedFtpPathForGenome <- NA
    referenceDataDF$usedFtpPathForCdss <- NA
    referenceDataDF$downloadStatus <- FALSE
  }
  
  for (i in 1:nrow(referenceDataDF)) {
    if(referenceDataDF$downloadStatus[i] == F) {
      
      ## Get serach term
      searchTerm <- referenceDataDF$searchTerm[i]
      
      ## Get summaries by search term
      summaryResults <- getSearchSummary(searchTerm)
      if (isSingleSummary(summaryResults)) {
        summaryResults <- list(summaryResults)
      }
      
      if (summaryResults %>% length > 20) {
        summaryResults <- getSearchSummary(searchTerm %>% paste0(' "complete genome"'))
      }
      
      ## Get Genomes and Cdss
      usedFtpPathForGenome <- NA
      usedFtpPathForCdss <- NA
      currentSummaryResultIndex <- 0
      summaryResultsNumber <- summaryResults %>% length
      areGenomesAndCdsOk <- FALSE
      while (areGenomesAndCdsOk == FALSE &&
             currentSummaryResultIndex < summaryResultsNumber) {
        ## Get a summary result
        currentSummaryResultIndex <- currentSummaryResultIndex + 1
        summaryResult <- summaryResults[[currentSummaryResultIndex]]
        
        ## Try Genbank First (GENOME)
        isGenomeOk <- FALSE
        usedFtpPathForGenome <- getFtpPath(summaryResult, "genbank")
        if ((usedFtpPathForGenome %>% is.na == FALSE) &&
            (usedFtpPathForGenome == "") == FALSE) {
          
          try(downloadGenome(genomesFolder, searchTerm, usedFtpPathForGenome))
          try(isGenomeOk <- checkGenome(genomesFolder, searchTerm, usedFtpPathForGenome))
        }
        
        ## If not ok, try Refseq (GENOME)
        if (isGenomeOk == FALSE) {
          usedFtpPathForGenome <- getFtpPath(summaryResult, "refseq")
          if ((usedFtpPathForGenome %>% is.na == FALSE) &&
              (usedFtpPathForGenome == "") == FALSE) {
            try(downloadGenome(genomesFolder, searchTerm, usedFtpPathForGenome))
            try(isGenomeOk <- checkGenome(genomesFolder, searchTerm, usedFtpPathForGenome))
          }
        }
        
        ## Try Genbank First (CDSs)
        isCdssOk <- FALSE
        usedFtpPathForCdss <- getFtpPath(summaryResult, "genbank")
        if ((usedFtpPathForCdss %>% is.na == FALSE) &&
            (usedFtpPathForCdss == "") == FALSE) {
          try(downloadCdss(cdssFolder, searchTerm, usedFtpPathForCdss))
          try(isCdssOk <- checkCdss(cdssFolder, searchTerm, usedFtpPathForCdss))
        }
        
        ## If not ok, try Refseq (CDSs)
        if (isCdssOk == FALSE) {
          usedFtpPathForCdss <- getFtpPath(summaryResult, "refseq")
          if ((usedFtpPathForCdss %>% is.na == FALSE) &&
              (usedFtpPathForCdss == "") == FALSE) {
            try(downloadCdss(cdssFolder, searchTerm, usedFtpPathForCdss))
            try(isCdssOk <- checkCdss(cdssFolder, searchTerm, usedFtpPathForCdss))
          }
        }
        
        ##
        print(paste(searchTerm, usedFtpPathForGenome, usedFtpPathForCdss, isGenomeOk, isCdssOk))
        
        ## Verification flag, TRUE if both Genome and CDSs is OK
        areGenomesAndCdsOk <- isGenomeOk && isCdssOk
      }
      
      # if (areGenomesAndCdsOk == FALSE) {
      #   stop("Ops! We could not found genomes and CDSs for search term: " %>% paste0(searchTerm))
      # }
      ## Remove downloaded data
      # if (areGenomesAndCdsOk == FALSE) {
      #   removeGenomeFilesCommand <- "rm GENOMES_FOLDER/" %>% paste0(searchTerm, ".*") %>% str_replace("GENOMES_FOLDER", genomesFolder)
      #   removeCdssFilesCommand <- "rm CDSS_FOLDER/" %>% paste0(searchTerm, ".*") %>% str_replace("CDSS_FOLDER", cdssFolder)
      #   system(removeGenomeFilesCommand)
      #   system(removeCdssFilesCommand)
      # }
      referenceDataDF$downloadStatus[i] <- areGenomesAndCdsOk
      
      ## Keep information
      referenceDataDF$usedFtpPathForGenome[i] <- usedFtpPathForGenome
      referenceDataDF$usedFtpPathForCdss[i] <- usedFtpPathForCdss
      
      summaryResult <- summaryResults[[currentSummaryResultIndex]]
      referenceDataDF$assemblyname[i] <- summaryResult$assemblyname
      referenceDataDF$taxid[i] <- summaryResult$taxid
      referenceDataDF$organism[i] <- summaryResult$organism
      referenceDataDF$speciestaxid[i] <- summaryResult$speciestaxid
      referenceDataDF$speciesname[i] <- summaryResult$speciesname
      referenceDataDF$ftppath_refseq[i] <- summaryResult$ftppath_refseq
      referenceDataDF$ftppath_genbank[i] <- summaryResult$ftppath_genbank
    }
  }
  
  return(referenceDataDF)
}

predictCdss <- function(searchTerm, genomesFolder, cdssFolder, prodigalLocation) {
  inputGenome <- genomesFolder %>% paste0("/") %>% paste0(searchTerm) %>% paste0(".fna")
  outputCdss <- cdssFolder %>% paste0("/") %>% paste0(searchTerm) %>% paste0(".fna")
  outputGff <- cdssFolder %>% paste0("/") %>% paste0(searchTerm) %>% paste0(".gff")
  
  prodigalCommand <- "PRODIGAL_LOCATION -i INPUT_GENOME -d OUTPUT_CDSS -f gff -o GFF_OUTPUT" %>%
    str_replace_all("PRODIGAL_LOCATION", prodigalLocation) %>%
    str_replace_all("INPUT_GENOME", inputGenome) %>%
    str_replace_all("OUTPUT_CDSS", outputCdss) %>%
    str_replace_all("GFF_OUTPUT", outputGff)
  system(prodigalCommand)
  
  sedCommand <- "sed -i 's/>/>REF /' FNA_FILENAME" %>%
    str_replace("REF", searchTerm) %>%
    str_replace("FNA_FILENAME", outputCdss)
  system(sedCommand)
  
  return(outputCdss)
}

predictMissingCdss <- function(genomesFolder, cdssFolder, referenceDataDF, prodigalLocation) {
  referenceDataDF$prodigalCds <- FALSE
  for (i in 1:nrow(referenceDataDF)) {
    if (referenceDataDF$downloadStatus[i] == FALSE &&
        referenceDataDF$genomeCheck[i] == TRUE) {
      predictCdss(referenceDataDF$searchTerm[i], genomesFolder, cdssFolder, prodigalLocation)
      referenceDataDF$prodigalCds[i] <- TRUE
    }
  }
  return(referenceDataDF)
}

createGrinderScripts <- function(projectFolder,
                                 shReadsNumber, mtReadsNumber, ampReadsNumber,
                                 readSize, readSizeStandardDeviation,
                                 ampChimeraPercentual) {
  BASE_GRINDER_COMMAND_FOR_SH <-
    "grinder -reference_file PROJECT_FOLDER/genomes.fna" %>%
    paste("-total_reads READS_NUMBER") %>%
    paste("-read_dist READ_SIZE") %>%
    paste("normal STANDARD_DEVIATION") %>%
    paste("-length_bias 1 -random_seed 1 -qual_levels 30 10 -fastq_output 1") %>%
    paste("-output_dir PROJECT_FOLDER/sh") %>%
    paste("-md poly4 3e-3 3.3e-8 -mr 80 20") %>%
    paste("-af PROJECT_FOLDER/sh-mt-ranks.txt")
  
  BASE_GRINDER_COMMAND_FOR_MT <-
    "grinder -reference_file PROJECT_FOLDER/cdss.fna" %>%
    paste("-total_reads READS_NUMBER") %>%
    paste("-read_dist READ_SIZE") %>%
    paste("normal STANDARD_DEVIATION") %>%
    paste("-length_bias 1 -random_seed 1 -qual_levels 30 10 -fastq_output 1") %>%
    paste("-output_dir PROJECT_FOLDER/mt") %>%
    paste("-md poly4 3e-3 3.3e-8 -mr 80 20") %>%
    paste("-af PROJECT_FOLDER/sh-mt-ranks.txt")
  
  BASE_GRINDER_COMMAND_FOR_AMP <-
    "grinder -reference_file PROJECT_FOLDER/genomes.fna" %>%
    paste("-total_reads READS_NUMBER") %>%
    paste("-read_dist READ_SIZE") %>%
    paste("normal STANDARD_DEVIATION") %>%
    paste("-copy_bias 1 -random_seed 1 -qual_levels 30 10 -fastq_output 1") %>%
    paste("-output_dir PROJECT_FOLDER/amp -md poly4 3e-3 3.3e-8 -mr 80 20") %>%
    paste("-af PROJECT_FOLDER/amp-ranks.txt -length_bias 0 -unidirectional 1") %>%
    paste("-chimera_perc CHIMERA_PERCENTUAL") %>%
    paste("-fr PROJECT_FOLDER/primers.fna")
  
  shCommand <- BASE_GRINDER_COMMAND_FOR_SH %>% str_replace_all("PROJECT_FOLDER", projectFolder)
  shCommand <- shCommand %>% str_replace_all("READS_NUMBER", shReadsNumber)
  shCommand <- shCommand %>% str_replace_all("READ_SIZE", readSize)
  shCommand <- shCommand %>%  str_replace_all("STANDARD_DEVIATION", readSizeStandardDeviation)
  
  mtCommand <- BASE_GRINDER_COMMAND_FOR_MT %>% str_replace_all("PROJECT_FOLDER", projectFolder)
  mtCommand <- mtCommand %>% str_replace_all("READS_NUMBER", mtReadsNumber)
  mtCommand <- mtCommand %>% str_replace_all("READ_SIZE", readSize)
  mtCommand <- mtCommand %>% str_replace_all("STANDARD_DEVIATION", readSizeStandardDeviation)
  
  ampCommand <- BASE_GRINDER_COMMAND_FOR_AMP %>% str_replace_all("PROJECT_FOLDER", projectFolder)
  ampCommand <- ampCommand %>% str_replace_all("READS_NUMBER", ampReadsNumber)
  ampCommand <- ampCommand %>% str_replace_all("READ_SIZE", readSize)
  ampCommand <- ampCommand %>% str_replace_all("STANDARD_DEVIATION", readSizeStandardDeviation)
  ampCommand <- ampCommand %>% str_replace_all("CHIMERA_PERCENTUAL", ampChimeraPercentual)
  
  sink(projectFolder %>% paste0("/grinderScriptForSh.sh"))
  cat(shCommand)
  sink(NULL)
  
  sink(projectFolder %>% paste0("/grinderScriptForMt.sh"))
  cat(mtCommand)
  sink(NULL)
  
  sink(projectFolder %>% paste0("/grinderScriptForAmp.sh"))
  cat(ampCommand)
  sink(NULL)
}

createGrinderJobs <- function(projectFolder, jobsBaseName = "Grinder", email = "", group = "", condaEnviroment = "") {
  dir.create(projectFolder %>% paste0("/jobs"), showWarnings = F, recursive = T)
  
  ## Sh
  grinderShJobLocation <- projectFolder %>% paste0("/grinderSh.job")
  sink(grinderShJobLocation)
  cat("#!/bin/bash\n")
  cat("#PBS -N JOB_NAME\n" %>% str_replace("JOB_NAME", paste0(jobsBaseName, "Sh")))
  cat("#PBS -S /bin/bash\n")
  cat("#PBS -l nodes=1:ppn=1,walltime=60:00:00,vmem=16gb\n")
  cat("#PBS -M EMAIL\n" %>% str_replace("EMAIL", email))
  cat("#PBS -m bea\n")
  cat("#PBS -d PROJECT_FOLDER\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -o PROJECT_FOLDER/jobs/sh.stdout\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -e PROJECT_FOLDER/jobs/sh.stderr\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -W group_list=GROUP\n" %>% str_replace("GROUP", group))
  cat("#PBS -V\n")
  cat("\n")
  cat("## Configure\n")
  if (condaEnviroment != "") {
    cat("source activate CONDA_ENVIROMENT\n" %>% str_replace("CONDA_ENVIROMENT", condaEnviroment))
  }
  cat("sh grinderScriptForSh.sh\n")
  if (condaEnviroment != "") {
    cat("source deactivate\n")
  }
  sink(NULL)
  
  ## Mt
  grinderMtJobLocation <- projectFolder %>% paste0("/grinderMt.job")
  sink(grinderMtJobLocation)
  cat("#!/bin/bash\n")
  cat("#PBS -N JOB_NAME\n" %>% str_replace("JOB_NAME", paste0(jobsBaseName, "Mt")))
  cat("#PBS -S /bin/bash\n")
  cat("#PBS -l nodes=1:ppn=1,walltime=30:00:00,vmem=16gb\n")
  cat("#PBS -M EMAIL\n" %>% str_replace("EMAIL", email))
  cat("#PBS -m bea\n")
  cat("#PBS -d PROJECT_FOLDER\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -o PROJECT_FOLDER/jobs/mt.stdout\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -e PROJECT_FOLDER/jobs/mt.stderr\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -W group_list=GROUP\n" %>% str_replace("GROUP", group))
  cat("#PBS -V\n")
  cat("\n")
  cat("## Configure\n")
  if (condaEnviroment != "") {
    cat("source activate CONDA_ENVIROMENT\n" %>% str_replace("CONDA_ENVIROMENT", condaEnviroment))
  }
  cat("sh grinderScriptForMt.sh\n")
  if (condaEnviroment != "") {
    cat("source deactivate\n")
  }
  sink(NULL)
  
  ## Amp
  grinderAmpJobLocation <- projectFolder %>% paste0("/grinderAmp.job")
  sink(grinderAmpJobLocation)
  cat("#!/bin/bash\n")
  cat("#PBS -N JOB_NAME\n" %>% str_replace("JOB_NAME", paste0(jobsBaseName, "Amp")))
  cat("#PBS -S /bin/bash\n")
  cat("#PBS -l nodes=1:ppn=1,walltime=30:00:00,vmem=16gb\n")
  cat("#PBS -M EMAIL\n" %>% str_replace("EMAIL", email))
  cat("#PBS -m bea\n")
  cat("#PBS -d PROJECT_FOLDER\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -o PROJECT_FOLDER/jobs/amp.stdout\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -e PROJECT_FOLDER/jobs/amp.stderr\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -W group_list=GROUP\n" %>% str_replace("GROUP", group))
  cat("#PBS -V\n")
  cat("\n")
  cat("## Configure\n")
  if (condaEnviroment != "") {
    cat("source activate CONDA_ENVIROMENT\n" %>% str_replace("CONDA_ENVIROMENT", condaEnviroment))
  }
  cat("sh grinderScriptForAmp.sh\n")
  if (condaEnviroment != "") {
    cat("source deactivate\n")
  }
  sink(NULL)
  
}

submitGrinderJobs <- function(projectFolder) {
  grinderShJobLocation <- projectFolder %>% paste0("/grinderSh.job")
  grinderMtJobLocation <- projectFolder %>% paste0("/grinderMt.job")
  grinderAmpJobLocation <- projectFolder %>% paste0("/grinderAmp.job")
  
  sh.command <- "qsub " %>% paste0(grinderShJobLocation)
  mt.command <- "qsub " %>% paste0(grinderMtJobLocation)
  amp.command <- "qsub " %>% paste0(grinderAmpJobLocation)
  
  system(sh.command)
  system(mt.command)
  system(amp.command)
}

createKrakenScripts <- function(projectFolder, sickleProgramFolder, krakenProgramFolder, krakenDatabaseFolder) {
  ## SH
  krakenShFileLocation <- projectFolder %>% paste0("/krakenScriptForSh.sh")
  sink(krakenShFileLocation)
  ## Run sickle
  cat("SICKLE_FOLDER/sickle se -f  PROJECT_FOLDER/sh/grinder-reads.fastq -t sanger -o PROJECT_FOLDER/sh/qc_filtered.fastq\n" %>%
        str_replace_all("SICKLE_FOLDER", sickleProgramFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Run kraken
  cat("KRAKEN_FOLDER/kraken --db KRAKEN_DATABASE --threads 12 PROJECT_FOLDER/sh/qc_filtered.fastq --fastq-input --output PROJECT_FOLDER/sh/kraken_classifications.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Add scores
  cat("KRAKEN_FOLDER/kraken-filter --db KRAKEN_DATABASE PROJECT_FOLDER/sh/kraken_classifications.tsv > PROJECT_FOLDER/sh/kraken_scores.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Generate report
  cat("KRAKEN_FOLDER/kraken-report --db KRAKEN_DATABASE PROJECT_FOLDER/sh/kraken_classifications.tsv > PROJECT_FOLDER/sh/kraken_report.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Generate labels
  cat("KRAKEN_FOLDER/kraken-translate --db KRAKEN_DATABASE PROJECT_FOLDER/sh/kraken_classifications.tsv --mpa-format > PROJECT_FOLDER/sh/kraken_labels.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  
  sink(NULL)
  
  ## MT
  krakenMtFileLocation <- projectFolder %>% paste0("/krakenScriptForMt.sh")
  sink(krakenMtFileLocation)
  ## Run sickle
  cat("SICKLE_FOLDER/sickle se -f  PROJECT_FOLDER/mt/grinder-reads.fastq -t sanger -o PROJECT_FOLDER/mt/qc_filtered.fastq\n" %>%
        str_replace_all("SICKLE_FOLDER", sickleProgramFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Run kraken
  cat("KRAKEN_FOLDER/kraken --db KRAKEN_DATABASE --threads 12 PROJECT_FOLDER/mt/qc_filtered.fastq --fastq-input --output PROJECT_FOLDER/mt/kraken_classifications.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Add scores
  cat("KRAKEN_FOLDER/kraken-filter --db KRAKEN_DATABASE PROJECT_FOLDER/mt/kraken_classifications.tsv > PROJECT_FOLDER/mt/kraken_scores.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Generate report
  cat("KRAKEN_FOLDER/kraken-report --db KRAKEN_DATABASE PROJECT_FOLDER/mt/kraken_classifications.tsv > PROJECT_FOLDER/mt/kraken_report.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  ## Generate labels
  cat("KRAKEN_FOLDER/kraken-translate --db KRAKEN_DATABASE PROJECT_FOLDER/mt/kraken_classifications.tsv --mpa-format > PROJECT_FOLDER/mt/kraken_labels.tsv\n" %>%
        str_replace_all("KRAKEN_FOLDER", krakenProgramFolder) %>%
        str_replace_all("KRAKEN_DATABASE", krakenDatabaseFolder) %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  
  sink(NULL)
}

createKrakenJobs <- function(projectFolder, jobsBaseName = "Kraken", email = "", group = "") {
  ## Sh
  krakenShJobLocation <- projectFolder %>% paste0("/krakenSh.job")
  sink(krakenShJobLocation)
  cat("#!/bin/bash\n")
  cat("#PBS -N JOB_NAME\n" %>% str_replace("JOB_NAME", paste0(jobsBaseName, "Sh")))
  cat("#PBS -S /bin/bash\n")
  cat("#PBS -l nodes=1:ppn=16,walltime=30:00:00,vmem=150gb\n")
  cat("#PBS -M EMAIL\n" %>% str_replace("EMAIL", email))
  cat("#PBS -m bea\n")
  cat("#PBS -d PROJECT_FOLDER\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -o PROJECT_FOLDER/jobs/krakenSh.stdout\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -e PROJECT_FOLDER/jobs/krakenSh.stderr\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -W group_list=GROUP\n" %>% str_replace("GROUP", group))
  cat("#PBS -V\n")
  cat("\n")
  cat("sh krakenScriptForSh.sh\n")
  sink(NULL)
  
  ## Mt
  krakenMtJobLocation <- projectFolder %>% paste0("/krakenMt.job")
  sink(krakenMtJobLocation)
  cat("#!/bin/bash\n")
  cat("#PBS -N JOB_NAME\n" %>% str_replace("JOB_NAME", paste0(jobsBaseName, "Mt")))
  cat("#PBS -S /bin/bash\n")
  cat("#PBS -l nodes=1:ppn=16,walltime=30:00:00,vmem=150gb\n")
  cat("#PBS -M EMAIL\n" %>% str_replace("EMAIL", email))
  cat("#PBS -m bea\n")
  cat("#PBS -d PROJECT_FOLDER\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -o PROJECT_FOLDER/jobs/krakenMt.stdout\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -e PROJECT_FOLDER/jobs/krakenMt.stderr\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -W group_list=GROUP\n" %>% str_replace("GROUP", group))
  cat("#PBS -V\n")
  cat("\n")
  cat("sh krakenScriptForMt.sh\n")
  sink(NULL)
}

createQiimeScript <- function(projectFolder, sampleName) {
  sink(projectFolder %>% paste0("/qiimeParameters.txt"))
  cat("pick_otus:otu_picking_method\tuclust\nassign_taxonomy:assignment_method\trdp")
  sink(NULL)
  
  sink(projectFolder %>% paste0("/qiimeScriptForAmp.sh"))
  cat("split_libraries_fastq.py -i PROJECT_FOLDER/amp/grinder-reads.fastq -o PROJECT_FOLDER/qiime/slout -q 19 --barcode_type 'not-barcoded' --sample_ids SAMPLE_NAME --phred_offset 33" %>%
        str_replace_all("PROJECT_FOLDER", projectFolder) %>%
        str_replace_all("SAMPLE_NAME", sampleName))
  cat("\n")
  cat("\n")
  cat("pick_open_reference_otus.py -i PROJECT_FOLDER/qiime/slout/seqs.fna -o PROJECT_FOLDER/qiime/open_reference_otus" %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  cat("\n")
  cat("\n")
  # cat("assign_taxonomy.py -i PROJECT_FOLDER/qiime/open_reference_otus/rep_set.fna -o PROJECT_FOLDER/qiime/open_reference_otus/rdp_assigned_taxonomy -m rdp -c 0.1" %>%
  #       str_replace_all("PROJECT_FOLDER", projectFolder))
  # cat("\n")
  # cat("\n")
  # cat("make_otu_table.py -i PROJECT_FOLDER/qiime/open_reference_otus/final_otu_map_mc2.txt -t PROJECT_FOLDER/qiime/open_reference_otus/rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o PROJECT_FOLDER/qiime/open_reference_otus/otu_table_mc2_w_tax.biom" %>%
  #       str_replace_all("PROJECT_FOLDER", projectFolder))
  # cat("\n")
  # cat("\n")
  cat("summarize_taxa.py -i PROJECT_FOLDER/qiime/open_reference_otus/otu_table_mc2_w_tax.biom -o PROJECT_FOLDER/qiime/open_reference_otus/taxa_summary -L 1,2,3,4,5,6,7" %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  cat("\n")
  cat("\n")
  cat("biom convert -i PROJECT_FOLDER/qiime/open_reference_otus/otu_table_mc2_w_tax.biom -o PROJECT_FOLDER/qiime/open_reference_otus/otu_table.txt --to-tsv --header-key taxonomy" %>%
        str_replace_all("PROJECT_FOLDER", projectFolder))
  cat("\n")
  sink(NULL)
}

createQiimeJob <- function(projectFolder, jobBaseName, email, group) {
  qiimeAmpJobLocation <- projectFolder %>% paste0("/qiimeAmp.job")
  sink(qiimeAmpJobLocation)
  cat("#!/bin/bash\n")
  cat("#PBS -N JOB_NAME\n" %>% str_replace("JOB_NAME", paste0(jobBaseName, "Amp")))
  cat("#PBS -S /bin/bash\n")
  cat("#PBS -l nodes=1:ppn=16,walltime=6:00:00,vmem=16gb\n")
  cat("#PBS -M EMAIL\n" %>% str_replace("EMAIL", email))
  cat("#PBS -m bea\n")
  cat("#PBS -d PROJECT_FOLDER\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -o PROJECT_FOLDER/jobs/qiimeAmp.stdout\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -e PROJECT_FOLDER/jobs/qiimeAmp.stderr\n" %>% str_replace("PROJECT_FOLDER", projectFolder))
  cat("#PBS -W group_list=GROUP\n" %>% str_replace("GROUP", group))
  cat("#PBS -V\n")
  cat("\n")
  cat("sh qiimeScriptForAmp.sh\n")
  sink(NULL)
}
##==========================================================================================================================================================