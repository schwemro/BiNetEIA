#' @title Degree centrality in a weighted network
#'
#' @description Modified function of package tnet. Like for the two-mode, every node is now considered in the returned data frame.
#' @param net see \code{\link{degree_w}}
#' @param measure see \code{\link{degree_w}}
#' @param type see \code{\link{degree_w}}
#' @param alpha see \code{\link{degree_w}}
#' @param nn Number of actions or environmental components, respectively.
#' @seealso \code{\link{degree_w}}

deg_w <-
  function(net, measure = c("degree", "output"), type = "out", alpha=1, nn){
    # Ensure that the network conforms to the tnet standard
    if (is.null(attributes(net)$tnet))
      net <- as.tnet(net, type = "weighted one-mode tnet")
    if (attributes(net)$tnet != "weighted one-mode tnet")
      stop("Network not loaded properly")

    # Reverse data if calculating in-degrees
    if(type == "in") {
      net <- data.frame(i=net[,2], j=net[,1], w=net[,3])
      net <- net[order(net[,"i"],net[,"j"]),]
    }
    # Create an index for each node
    index <- cumsum(!duplicated(net[,1]))
    # Create output object
    k.list <- cbind(unique(net[,1]), NaN, NaN, NaN)
    # Assign names
    dimnames(k.list)[[2]] <- c("node", "degree", "output","alpha")
    # Calculating degree
    if(sum(c("degree","alpha") %in% measure)>0)
      k.list[,"degree"] <- tapply(net[, "w"], index, length)
    # Calculating strength
    if(sum(c("output","alpha") %in% measure)>0)
      k.list[,"output"] <- tapply(net[, "w"], index, sum)
    # Calculating alpha parameter
    if("alpha" %in% measure)
      k.list[,"alpha"] <- k.list[,"degree"]*((k.list[,"output"]/k.list[,"degree"])^alpha)
    # Add rows to the output object if isolates exists
    if(nn != nrow(k.list)) {
      k.list <- rbind(k.list, cbind(node=1:nn, 0, 0, 0))
      k.list <- k.list[order(k.list[, "node"]), ]
      k.list <- k.list[!duplicated(k.list[, "node"]), ]
    }
    # Extract just relevant columns, and return
    return(k.list[,c("node", measure)])
  }
