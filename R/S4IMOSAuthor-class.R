## Class definition
#'@exportClass S4IMOSAuthor
S4IMOSAuthor <- setClass(
  "S4IMOSAuthor",
  
  representation = representation(
    name = "character",
    email = "character"
  ),
  
  prototype = list(
    name = "",
    email = ""
  ),
  
  validity = function(object) {
    return(TRUE)
  }
)

## Setters and getters

##--

setGeneric("getName", function(self) {
  standardGeneric("getName")
})

setMethod("getName", "S4IMOSAuthor", function(self) {
  return(self@name)
})

setGeneric("setName", function(self, name) {
  standardGeneric("setName")
})

setMethod("setName", "S4IMOSAuthor", function(self, name){
  if (name %>% is.null) {
    cat("name should not be NULL\n\n")
    return(self)
  }
  
  if (name %>% is.na) {
    cat("name should not be NA\n\n")
    return(self)
  }
  
  if (name %>% trimws == "") {
    cat("name should not be empty\n\n")
    return(self)
  }
  
  self@name <- name
  return(self)
})

##--

setGeneric("getEmail", function(self) {
  standardGeneric("getEmail")
})

setMethod("getEmail", "S4IMOSAuthor", function(self) {
  return(self@email)
})

setGeneric("setEmail", function(self, email) {
  standardGeneric("setEmail")
})

setMethod("setEmail", "S4IMOSAuthor", function(self, email){
  if (email %>% is.null) {
    cat("email should not be NULL\n\n")
    return(self)
  }
  
  if (email %>% is.na) {
    cat("email should not be NA\n\n")
    return(self)
  }
  
  if (email %>% trimws == "") {
    cat("email should not be empty\n\n")
    return(self)
  }
  
  self@email <- email
  return(self)
})
