
# GO dataset to list ------------------------------------------------------

go <- readLines("Data/go.obo.txt")
go <- go[30:(grep("Typedef", go)[1]-1)]

go <- gsub(pattern = "\\[Term\\]", replacement = "", x = go)
go <- go[!go == ""]

n <- grep("id: GO", go)
go1 <- list()
for(i in 1:length(n)){
  if(i == length(n)){
    go1[[as.character(go[n[i]])]] <- as.vector(go[(n[i]+1):length(go)])
  } else {
    go1[[as.character(go[n[i]])]] <- as.vector(go[(n[i]+1):(n[i+1]-1)])
  }
}

go1[grep("name: biological_process", go1)]


# biocLite("gageData")
library(gageData)
data(go.sets.hs)
go.names <- names(go.sets.hs)[unlist(lapply(go.sel$term[go.sel$geneSet == 4], grep, names(go.sets.hs)))]

nam <- unlist(lapply(go.names, function(x) substr(x, start = 1, stop = 10)))


## example ->  looking for parents term!! 
for(name in nam){
  term <- go1[paste0("id: ", name)]
  if (names(term) == "id: GO:0008150"){
    stop
  }
  is_a <- term[[1]][grep("is_a", term[[1]])]
  unlist(lapply(is_a, function(x) substr(x, start = 7, stop = 16)))
}

