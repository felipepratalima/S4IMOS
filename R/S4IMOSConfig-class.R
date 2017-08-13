## Required packages
require(magrittr)


## Helper Constants
DEFAULT_APPLICATION_FOLDER <- "~/s4imos_data"
DEFAULT_GENOMES_FOLDER <- "~/s4imos_data/genomes"


#'@export
S4IMOSConfig <- setClass(
  "S4IMOSConfig",
  
  representation = representation(
    applicationFolder = "character",
    genomesFolder = "character",
    useDocker = "logical"
  ),
  
  prototype = prototype(
    applicationFolder = DEFAULT_APPLICATION_FOLDER,
    genomesFolder = DEFAULT_GENOMES_FOLDER,
    useDocker = TRUE
  ),
  
  validity = function(object) {
    return(TRUE)
  }
)


## Setters and getters

# --

setGeneric("getApplicationFolder", function(self) {
  standardGeneric("getApplicationFolder")
})

setMethod("getApplicationFolder", "S4IMOSConfig", function(self) {
  return(self@applicationFolder)
})

setGeneric("setApplicationFolder", function(self, applicationFolder) {
  standardGeneric("setApplicationFolder")
})

setMethod("setApplicationFolder", "S4IMOSConfig", function(self, applicationFolder = DEFAULT_APPLICATION_FOLDER) {
  if (applicationFolder %>% is.null) {
    cat("applicationFolder should not be NULL\n\n")
    return(self)
  }
  
  if (applicationFolder %>% is.na) {
    cat("applicationFolder should not be NA\n\n")
    return(self)
  }
  
  if (applicationFolder %>% trimws == "") {
    cat("applicationFolder should not be empty\n\n")
    return(self)
  }
  
  applicationFolder <- applicationFolder %>% path.expand
  self@applicationFolder <- applicationFolder
  return(self)
})

# --

setGeneric("getGenomesFolder", function(self) {
  standardGeneric("getGenomesFolder")
})
setMethod("getGenomesFolder", "S4IMOSConfig", function(self) {
  return(self@genomesFolder)
})

setGeneric("setGenomesFolder", function(self, genomesFolder) {
  standardGeneric("setGenomesFolder")
})
setMethod("setGenomesFolder", "S4IMOSConfig", function(self, genomesFolder = DEFAULT_GENOMES_FOLDER) {
  if (genomesFolder %>% is.null) {
    cat("genomesFolder should not be NULL\n\n")
    return(self)
  }
  
  if (genomesFolder %>% is.na) {
    cat("genomesFolder should not be NA\n\n")
    return(self)
  }
  
  if (genomesFolder %>% trimws == "") {
    cat("genomesFolder should not be empty\n\n")
    return(self)
  }
  
  genomesFolder <- genomesFolder %>% path.expand
  self@genomesFolder <- genomesFolder
  return(self)
})

# --

setGeneric("isUseDocker", function(self) {
  standardGeneric("isUseDocker")
})
setMethod("isUseDocker", "S4IMOSConfig", function(self) {
  return(self@useDocker)
})

setGeneric("setUseDocker", function(self, useDocker) {
  standardGeneric("setUseDocker")
})
setMethod("setUseDocker", "S4IMOSConfig", function(self, useDocker = TRUE) {
  if (useDocker %>% is.null) {
    cat("genomesFolder should not be NULL\n\n")
    return(self)
  }
  
  if (useDocker %>% is.na) {
    cat("genomesFolder should not be NA\n\n")
    return(self)
  }
  
  self@useDocker <- useDocker
  return(self)
})
