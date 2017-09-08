## useful isa functions
# Modules informations -----------------------------------------------------
isaModules <- function(data.isa, type = c("isa", "ppa")) {
  if(type == "isa"){
    colGroups <- matrix(NA, ncol(data.isa$rows))
    rowGroups <- matrix(NA, ncol(data.isa$rows))
    for (i in 1:ncol(data.isa$rows)) {
      colGroups[i] <- sum(data.isa$columns[,i] != 0)
      rowGroups[i] <- sum(data.isa$rows[,i] != 0)
    }
    modules <- data.frame(colGroups, rowGroups, rob = data.isa$seeddata$rob,  
                          thr.row = data.isa$seeddata$thr.row, thr.col = data.isa$seeddata$thr.col)
    colnames(modules) <- c("colGroups", "rowGroups", "rob", "thr.row", "thr.col")
    rownames(modules) <- paste0("Module", 1:nrow(modules))
  }
  if(type == "ppa"){
    colGroups <- matrix(NA, ncol(data.isa$rows1))
    row1Groups <- matrix(NA, ncol(data.isa$rows1))
    row2Groups <- matrix(NA, ncol(data.isa$rows2))
    for (i in 1:ncol(data.isa$rows1)) {
      colGroups[i] <- sum(data.isa$columns[,i] != 0)
      row1Groups[i] <- sum(data.isa$rows1[,i] != 0)
      row2Groups[i] <- sum(data.isa$rows2[,i] != 0)
    }
    modules <- data.frame(colGroups, row1Groups, row2Groups, rob = data.isa$seeddata$rob,  
                          thr.row1 = data.isa$seeddata$thr.row1, 
                          thr.row2 = data.isa$seeddata$thr.row2, thr.col = data.isa$seeddata$thr.col)
    colnames(modules) <- c("colGroups", "row1Groups","row2Groups", "rob", "thr.row1", "thr.row2", "thr.col")
    rownames(modules) <- paste0("Module", 1:nrow(modules))
  }
  print(modules)
}
# Extract module informations ----------------------------------------------------
isaRowNames <- function(data, data2 = NULL, data.isa, type = "isa", n){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- as.matrix(data[isaRow, isaCol, drop=FALSE])
    return(rownames(module.n))
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- t(data)[isaRow1, isaCol, drop=FALSE]
    module.n2 <- t(data2)[isaRow2, isaCol, drop=FALSE]
    return(list(data1 = colnames(module.n1), data2 = colnames(module.n2)))
  }
}

isaColNames <- function(data, data2 = NULL, data.isa, type = "isa", n){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- as.matrix(data[isaRow, isaCol, drop=FALSE])
    return(colnames(module.n))
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- t(data)[isaRow1, isaCol, drop=FALSE]
    module.n2 <- t(data2)[isaRow2, isaCol, drop=FALSE]
    return(list(data1 = rownames(module.n1), data2 = rownames(module.n2)))
  }
}

isaScore <- function(data, data.isa, Row = FALSE, Col = FALSE, n){
    isaRow = data.isa$rows[, n] != 0
    isaNames <- rownames(data)[isaRow]
    scoreRow <- data.frame(RowScore = data.isa$row[isaRow, n], row.names = isaNames)

    isaCol = data.isa$columns[, n] != 0
    isaNames <- colnames(data)[isaCol]
    scoreCol <- data.frame(ColScore = data.isa$column[isaCol, n], row.names = isaNames)

  return(list(scoreRow, scoreCol))
}

# Module visualization ----------------------------------------------------

isa2image <- function(data, data2 = NULL, data.isa, type = "isa", n, name1 = NULL, name2 = NULL, cex = 0.6, color1 = "red", color2 = "yellow", all = FALSE){
  if(type == "isa"){
    colors <- colorRampPalette(c(color1, color2))(n = 10000)
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
    ColorUsed <- colors[round(1+(min(module.n)-min(data))*10000/(max(data)-min(data))) : round( (max(module.n)-min(data))*10000/(max(data)-min(data)) )]
    image(module.n, axes = F, main = paste("Module",  n), col=ColorUsed)
    mtext(text=colnames(module.n), side=2, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = cex)
    mtext(text=rownames(module.n), side=1, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = cex)
    # if(all == TRUE){
    #   allCol <- c(colnames(module.n), rownames(data)[!rownames(data) %in% colnames(module.n)])
    #   allRow <- c(rownames(module.n), colnames(data)[!colnames(data) %in% rownames(module.n)])
    #   image(as.matrix(t(data[allCol, allRow])), axes = F)
    # }
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- t(data)[isaRow1, isaCol, drop=FALSE]
    module.n2 <- t(data2)[isaRow2, isaCol, drop=FALSE]
    par(mfrow = c(1,2))
    image(module.n1, axes = F, main = paste("Module",  n, name1))
    mtext(text=colnames(module.n1), side=2, line=0.3, at=seq(0,1,l=ncol(module.n1)), las=2, cex = cex)
    mtext(text=rownames(module.n1), side=1, line=0.3, at=seq(0,1,l=nrow(module.n1)), las=2, cex = cex)
    image(module.n2, axes = F, main = paste("Module",  n, name2))
    mtext(text=colnames(module.n2), side=2, line=0.3, at=seq(0,1,l=ncol(module.n2)), las=2, cex = cex)
    mtext(text=rownames(module.n2), side=1, line=0.3, at=seq(0,1,l=nrow(module.n2)), las=2, cex = cex)
    par(mfrow = c(1,1))
  }
  if(all == TRUE){
    isaRow <- c(rownames(data)[data.isa$rows[, n] != 0], rownames(data)[data.isa$rows[, n] == 0])
    isaCol <- c(colnames(data)[data.isa$columns[, n] != 0], colnames(data)[data.isa$columns[, n] == 0])
    module.n   <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
    image(module.n, axes = F, main = paste("Module",  n))
    mtext(text=colnames(module.n), side=2, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = cex)
    mtext(text=rownames(module.n), side=1, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = cex)
  }
}



# Identity function -------------------------------------------------------

## Identity function
identity <- function(data1, data2, data.isa1, data.isa2, modules1, modules2, sel = 0, Col = FALSE){
  id.matr <- matrix(NA, nrow = nrow(modules1), ncol = nrow(modules2))
  if(Col == TRUE){
    for(i in 1:nrow(modules1)){
      id1 <- isaColNames(data = data1, type = "isa", data.isa = data.isa1, n = i)
      for(j in 1:nrow(modules2)){
        id2 <- isaColNames(data = data2, type = "isa", data.isa = data.isa2, n = j)
        perc <- round(ifelse(length(id1) > length(id2), sum(id1 %in% id2)/length(id1), sum(id2 %in% id1)/length(id2)), 2)
        id.matr[i,j] <- perc
      }
    }
  } else {
    for(i in 1:nrow(modules1)){
      id1 <- isaRowNames(data = data1, type = "isa", data.isa = data.isa1, n = i)
      for(j in 1:nrow(modules2)){
        id2 <- isaRowNames(data = data2, type = "isa", data.isa = data.isa2, n = j)
        perc <- round(ifelse(length(id1) > length(id2), sum(id1 %in% id2)/length(id1), sum(id2 %in% id1)/length(id2)), 2)
        id.matr[i,j] <- perc
      }
    }
  }
  rownames(id.matr) <- paste0("mod", seq(1,nrow(modules1)))
  colnames(id.matr) <- paste0("mod", seq(1,nrow(modules2)))
  
  identity.sel <- id.matr[apply(id.matr, MARGIN = 1, function(x) any(x %in% seq(sel, 0.99, by = 0.001))), ]
  identity.sel <- identity.sel[,apply(identity.sel, MARGIN = 2, function(x) any(x %in% seq(sel, 0.99, by = 0.001)))]
  
  return(identity.sel)
  
}

