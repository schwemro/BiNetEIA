#' @title Visualization of the network and the basic network statistic
#'
#' @description The bipartite network from the interaction matrix is build, visualized and the main network statistics are computed
#' @param mat.int Interaction matrix as matrix
#' @param negative.values Defines if interaction matrix contains negative interactions. Default is set to TRUE.
#' @param vcat A vector which contains the labels of the vertice categories.
#' @param vcat.env.comps A vector which contains the vertice categories of the environmental components.
#' @param vcat.actions A vector which contains the vertice categories of the actions.
#' @export
#' @details ...
#' @return Two data frames. \code{bas.stats} contains the results of the basic network statistic and \code{deg} includes the weighted vertice degrees and mean degrees.
#' @author Robin Schwemmle \email{robin.schwemmle@venus.uni-freiburg.de}
#' @seealso \code{\link{assortativity_nominal}}, \code{\link{degree_tm}}, \code{\link{degree_w}}, \code{\link{distance_tm}}, \code{\link{distance_w}}, \code{\link{clustering_tm}}
#' @examples
#' data('tallahala')
#' vcat <- c("Water", "Land", "Biology", "Socioeconomy", "Infrastructure")
#' vcat.env.comps <- c(1,1,1,1,2,2,1,1,2,2,3,3,3,3,3,3,3,4,4,4,4)
#' vcat.actions <- c(1,1,1,5,5,5,1,3,3,5)
#' visstats(tallahala, negative.values = TRUE, vcat = vcat, vcat.env.comps = vcat.env.comps,
#' vcat.actions = vcat.actions)
#'
#' data('eastkarun')
#' vcat <- c("Water", "Land", "Biology", "Socioeconomy")
#' vcat.env.comps <- c(1,2,2,2,2,1,1,3,3,3,3,3,3,4,4,4,4,4,4,4)
#' vcat.actions <- c(1,1,1,3,1)
#' visstats(eastkarun, negative.values = TRUE, vcat = vcat, vcat.env.comps = vcat.env.comps,
#' vcat.actions = vcat.actions)
#'
#' data('kladovo')
#' vcat <- c("Water", "Land", "Biology", "Socioeconomy", "Infrastructure")
#' vcat.env.comps <- c(1,2,2,2,2,2,3,3,3,3,3,2,2,4,4,4)
#' vcat.actions <- c(5,5,5,5,5,5,5,5,5)
#' visstats(kladovo, negative.values = FALSE, vcat = vcat, vcat.env.comps = vcat.env.comps,
#' vcat.actions = vcat.actions)
#'
#' @import igraph
#' @import tnet
#' @import stringr
#' @import moments

visstats <- function(mat.int, negative.values = TRUE, vcat, vcat.env.comps, vcat.actions){

  if (negative.values == TRUE){

    # short numeric labels for plotting
    rownames(mat.int) <- seq(1,nrow(mat.int),1)
    colnames(mat.int) <- seq(1,ncol(mat.int),1)

    # separating positive and negative interactions into two matrices
    mat.int.pos <- mat.int
    mat.int.pos[mat.int.pos < 0] <- 0
    mat.int.neg <- mat.int
    mat.int.neg[mat.int.neg > 0] <- 0
    mat.int.neg <- mat.int.neg*(-1)

    # create igraph and tnet objects from interaction matrix
    ig.int <- graph_from_incidence_matrix(mat.int, weighted = T)
    net.int <- as.tnet(mat.int, type = "weighted two-mode tnet")
    # create igraph and tnet objects from positve interaction matrix
    ig.int.pos <- graph_from_incidence_matrix(mat.int.pos, weighted = T)
    # check if network is bipartite
    if (is_bipartite(ig.int.pos) == FALSE){
      stop("Network is not bipartite!")
      }
    net.int.pos <- as.tnet(mat.int.pos, type = "weighted two-mode tnet")
    # create igraph and tnet objects from negative interaction matrix
    ig.int.neg <- graph_from_incidence_matrix(mat.int.neg, weighted = T)
    # check if network is bipartite
    if (is_bipartite(ig.int.neg) == FALSE){
      stop("Network is not bipartite!")
      }
    net.int.neg <- as.tnet(mat.int.neg, type = "weighted two-mode tnet")

    ## using tnet for one-mode-projection and to create igraph object
    ## multi-edges are collapsed to single edge where the sum of the weights of the multi-edges is assigned
    ## create igraph and tnet objects for one-mode-projection from positive interaction matrix
    net.int.row.pos <- projecting_tm(net.int.pos, method = "sum")
    net.int.col.pos <- projecting_tm(net.int.pos[,c(2,1,3)], method = "sum")
    ig.int.row.pos <- graph_from_edgelist(as.matrix(net.int.row.pos[,c(1,2)]), directed = F)
    E(ig.int.row.pos)$weight <- net.int.row.pos[,3]
    ig.int.row.pos  <- simplify(ig.int.row.pos, edge.attr.comb="sum")
    ig.int.row.pos <- check_vertices(ig.int.row.pos, no_v=nrow(mat.int))
    ig.int.col.pos <- graph_from_edgelist(as.matrix(net.int.col.pos[,c(1,2)]), directed = F)
    E(ig.int.col.pos)$weight <- net.int.col.pos[,3]
    ig.int.col.pos <- simplify(ig.int.col.pos, edge.attr.comb="sum")
    ig.int.col.pos <- check_vertices(ig.int.col.pos, no_v=ncol(mat.int))

    ## using tnet for one-mode-projection and to create igraph object
    ## multi-edges are collapsed to single edge where the sum of the weights of the multi-edges is assigned
    ## create igraph and tnet objects for one-mode-projection from negative interaction matrix
    net.int.row.neg <- projecting_tm(net.int.neg, method = "sum")
    net.int.col.neg <- projecting_tm(net.int.neg[,c(2,1,3)], method = "sum")
    ig.int.row.neg <- graph_from_edgelist(as.matrix(net.int.row.neg[,c(1,2)]), directed = F)
    E(ig.int.row.neg)$weight <- net.int.row.neg[,3]
    ig.int.row.neg <- simplify(ig.int.row.neg, edge.attr.comb="sum")
    ig.int.row.neg <- check_vertices(ig.int.row.neg, no_v=nrow(mat.int))
    ig.int.col.neg <- graph_from_edgelist(as.matrix(net.int.col.neg[,c(1,2)]), directed = F)
    E(ig.int.col.neg)$weight <- net.int.col.neg[,3]
    ig.int.col.neg <- simplify(ig.int.col.neg, edge.attr.comb="sum")
    ig.int.col.neg <- check_vertices(ig.int.col.neg, no_v=ncol(mat.int))

    ## data frame with basic statistics for two-mode and one-mode
    bas.stats <- as.data.frame(matrix(data = NA, nrow = 9, ncol = 8))
    colnames(bas.stats) <- c(' ', 'Two-mode', ' ', ' ',
                             ' ', 'One-mode',' ', ' ')
    rownames(bas.stats) <- c(' ', '  ', 'Number_of_vertices', 'Number_of_edges','Network_asymmetry',
                             'Connectance', 'Mean_distance', 'Cluster_coefficient',  'Assortativity_coefficient')
    bas.stats[1,] <- c('positive interactions', ' ', 'negative interactions', ' ',
                       'positive interactions', ' ', 'negative interactions', ' ')
    bas.stats[2,] <- c('Actions', 'Environmental components', 'Actions', 'Environmental components',
                       'Actions', 'Environmental components','Actions', 'Environmental components')

    ## number of vertices
    bas.stats[3,1]  <- ncol(mat.int)
    bas.stats[3,2]  <- nrow(mat.int)
    bas.stats[3,3]  <- ncol(mat.int)
    bas.stats[3,4]  <- nrow(mat.int)
    bas.stats[3,5]  <- ncol(mat.int)
    bas.stats[3,6]  <- nrow(mat.int)
    bas.stats[3,7]  <- ncol(mat.int)
    bas.stats[3,8]  <- nrow(mat.int)

    ## number of edges
    bas.stats[4,1]  <- gsize(ig.int.pos)
    bas.stats[4,2]  <- gsize(ig.int.pos)
    bas.stats[4,3]  <- gsize(ig.int.neg)
    bas.stats[4,4]  <- gsize(ig.int.neg)
    bas.stats[4,5]  <- gsize(ig.int.col.pos)
    bas.stats[4,6]  <- gsize(ig.int.row.pos)
    bas.stats[4,7]  <- gsize(ig.int.col.neg)
    bas.stats[4,8]  <- gsize(ig.int.row.neg)

    ## network asymmetry
    bas.stats[5,1] <- (nrow(mat.int)-ncol(mat.int))/(nrow(mat.int)+ncol(mat.int))
    bas.stats[5,5] <- 0

    ## connectance
    bas.stats[6,1] <- gsize(ig.int.pos)/(nrow(mat.int)*ncol(mat.int))
    bas.stats[6,3] <- gsize(ig.int.neg)/(nrow(mat.int)*ncol(mat.int))
    bas.stats[6,5] <- gsize(ig.int.col.pos)/(ncol(mat.int)*ncol(mat.int))
    bas.stats[6,6] <- gsize(ig.int.row.pos)/(nrow(mat.int)*nrow(mat.int))
    bas.stats[6,7] <- gsize(ig.int.col.neg)/(ncol(mat.int)*ncol(mat.int))
    bas.stats[6,8] <- gsize(ig.int.row.neg)/(nrow(mat.int)*nrow(mat.int))

    ## shortest paths and mean distance
    # two-mode
    dist.int.row.pos.2 <- distance_tm(net.int.pos[,c(1,2,3)], projection.method="sum", seed=42)
    dist.int.col.pos.2 <- distance_tm(net.int.pos[,c(2,1,3)], projection.method="sum", seed=42)
    dist.int.row.neg.2 <- distance_tm(net.int.neg[,c(1,2,3)], projection.method="sum", seed=42)
    dist.int.col.neg.2 <- distance_tm(net.int.neg[,c(2,1,3)], projection.method="sum", seed=42)

    bas.stats[7,1] <- mean(dist.int.col.pos.2, na.rm = T)
    bas.stats[7,2] <- mean(dist.int.row.pos.2, na.rm = T)
    bas.stats[7,3] <- mean(dist.int.col.neg.2, na.rm = T)
    bas.stats[7,4] <- mean(dist.int.row.neg.2, na.rm = T)

    # one-mode
    dist.int.row.pos.1 <- distance_w(net.int.row.pos, directed=F)
    dist.int.col.pos.1 <- distance_w(net.int.col.pos, directed=F)
    dist.int.row.neg.1 <- distance_w(net.int.row.neg, directed=F)
    dist.int.col.neg.1 <- distance_w(net.int.col.neg, directed=F)

    bas.stats[7,5] <- mean(dist.int.col.pos.1, na.rm = T)
    bas.stats[7,6] <- mean(dist.int.row.pos.1, na.rm = T)
    bas.stats[7,7] <- mean(dist.int.col.neg.1, na.rm = T)
    bas.stats[7,8] <- mean(dist.int.row.neg.1, na.rm = T)

    ## clustering
    # global cluster coefficient
    # two mode
    # in which way are weights considered?
    bas.stats[8,1] <- clustering_tm(net.int.pos[,c(2,1,3)])[3] #geometric mean
    bas.stats[8,2] <- clustering_tm(net.int.pos[,c(1,2,3)])[3]
    bas.stats[8,3] <- clustering_tm(net.int.neg[,c(2,1,3)])[3]
    bas.stats[8,4] <- clustering_tm(net.int.neg[,c(1,2,3)])[3]

    # one mode
    bas.stats[8,5] <- clustering_w(net.int.col.pos, measure = "gm") #geometric mean
    bas.stats[8,6] <- clustering_w(net.int.row.pos, measure = "gm")
    bas.stats[8,7] <- try(clustering_w(net.int.col.neg, measure = "gm"), TRUE)
    bas.stats[8,8] <- try(clustering_w(net.int.row.neg, measure = "gm"), TRUE)

    ## nominal assortativity by categories
    # check if category vectors are empty
    if (is.null(vcat) == TRUE){
      stop("Vector with the labels for the vertice categories is empty!")
    }
    if (is.null(vcat.env.comps) == TRUE){
      stop("Vector with categories for environmental components is empty!")
    }
    if (is.null(vcat.actions) == TRUE){
      stop("Vector with categories for actions is empty!")
    }

    # check if string categories and numeric categories have the same dimension
    vcat.num <- unique(c(vcat.env.comps, vcat.actions))
    if (length(vcat) != length(vcat.num)){
      stop("Vectors with category labels and numerical categories do not have the same dimension!")
    }

    # assign categories as vertice attributes
    V(ig.int.pos)$cat <- c(vcat.env.comps, vcat.actions)
    V(ig.int.neg)$cat <- c(vcat.env.comps, vcat.actions)
    V(ig.int.row.pos)$cat <- vcat.env.comps
    V(ig.int.row.neg)$cat <- vcat.env.comps
    V(ig.int.col.pos)$cat <- vcat.actions
    V(ig.int.col.neg)$cat <- vcat.actions
    # two-mode
    bas.stats[9,1] <- assortativity_nominal(ig.int.pos, types=V(ig.int.pos)$cat, directed = F)
    bas.stats[9,3] <- assortativity_nominal(ig.int.neg, types=V(ig.int.neg)$cat, directed = F)

    # one-mode
    bas.stats[9,5] <- assortativity_nominal(ig.int.col.pos, types=V(ig.int.col.pos)$cat, directed = F)
    bas.stats[9,6] <- assortativity_nominal(ig.int.row.pos, types=V(ig.int.row.pos)$cat, directed = F)
    bas.stats[9,7] <- assortativity_nominal(ig.int.col.neg, types=V(ig.int.col.neg)$cat, directed = F)
    bas.stats[9,8] <- assortativity_nominal(ig.int.row.neg, types=V(ig.int.row.neg)$cat, directed = F)

    bas.stats[is.na(bas.stats)] <- " "
    rownames(bas.stats) <- str_replace_all(rownames(bas.stats), "_", " ")

    ## node degree & mean degree
    deg <- matrix(data = NA, nrow = nrow(mat.int)+2, ncol = 10)
    colnames(deg) <- c(' ', 'Two-mode', ' ', 'One-mode', ' ',
                       'One-mode', ' ','Two-mode', ' ',' ')
    deg[1,] <- c(' ', 'positive interactions', 'negative interactions', 'positive interactions', 'negative interactions',
                 'negative interactions', 'positive interactions', 'negative interactions', 'positive interactions', ' ')
    deg[2:(ncol(mat.int)+1),1] <- seq(1, ncol(mat.int), 1)
    deg[2:(nrow(mat.int)+1),10] <- seq(1, nrow(mat.int), 1)

    deg[2:(ncol(mat.int)+1),2] <- deg_tm(net.int.pos[,c(2,1,3)], measure="output", nn = ncol(mat.int))[,2]
    mean.deg2.act.pos <- round(sum(as.numeric(deg[2:(ncol(mat.int)+1),2]))/ncol(mat.int))
    deg[(ncol(mat.int)+2),2] <- mean.deg2.act.pos
    deg[2:(ncol(mat.int)+1),3] <- deg_tm(net.int.neg[,c(2,1,3)], measure="output", nn = ncol(mat.int))[,2]
    mean.deg2.act.neg <- round(sum(as.numeric(deg[2:(ncol(mat.int)+1),3]))/ncol(mat.int))
    deg[(ncol(mat.int)+2),3] <- mean.deg2.act.neg

    deg[2:(nrow(mat.int)+1),8] <- deg_tm(net.int.neg, measure="output", nn = nrow(mat.int))[,2]
    mean.deg2.env.comp.neg <- round(sum(as.numeric(deg[2:(nrow(mat.int)+1),8]))/nrow(mat.int))
    deg[(nrow(mat.int)+2),8] <- mean.deg2.env.comp.neg
    deg[2:(nrow(mat.int)+1),9] <- deg_tm(net.int.pos, measure="output", nn = nrow(mat.int))[,2]
    mean.deg2.env.comp.pos <- round(sum(as.numeric(deg[2:(nrow(mat.int)+1),9]))/nrow(mat.int))
    deg[(nrow(mat.int)+2),9] <- mean.deg2.env.comp.pos

    # one-mode
    deg[2:(ncol(mat.int)+1),4] <- deg_w(net.int.col.pos, measure="output", nn = ncol(mat.int))[,2]
    deg[(ncol(mat.int)+2),4] <- sum(as.numeric(deg[2:(ncol(mat.int)+1),4]))/ncol(mat.int)
    deg[2:(ncol(mat.int)+1),5] <- deg_w(net.int.col.neg, measure="output", nn = ncol(mat.int))[,2]
    deg[(ncol(mat.int)+2),5] <- sum(as.numeric(deg[2:(ncol(mat.int)+1),5]))/ncol(mat.int)

    deg[2:(nrow(mat.int)+1),6] <- deg_w(net.int.row.neg, measure="output", nn = nrow(mat.int))[,2]
    deg[(nrow(mat.int)+2),6] <- sum(as.numeric(deg[2:(nrow(mat.int)+1),6]))/nrow(mat.int)
    deg[2:(nrow(mat.int)+1),7] <- deg_w(net.int.row.pos, measure="output", nn = nrow(mat.int))[,2]
    deg[(nrow(mat.int)+2),7] <- sum(as.numeric(deg[2:(nrow(mat.int)+1),7]))/nrow(mat.int)

    deg[is.na(deg)] <- " "

    ## degree distribution with relative frequencies
    deg.dist.act.pos <- deg_tm(net.int.pos[,c(2,1,3)], measure="output", nn = ncol(mat.int))[,2]
    deg.dist.act.pos.1 <- as.data.frame(table(deg.dist.act.pos))
    colnames(deg.dist.act.pos.1) <- c('degree','count')
    deg.dist.act.pos.1[,1] <- as.numeric(as.character(deg.dist.act.pos.1[,1]))
    deg.dist.act.pos.1[,2] <- as.numeric(deg.dist.act.pos.1[,2])

    deg.dist.act.neg <- deg_tm(net.int.neg[,c(2,1,3)], measure="output", nn = ncol(mat.int))[,2]
    deg.dist.act.neg.1 <- as.data.frame(table(deg.dist.act.neg))
    colnames(deg.dist.act.neg.1) <- c('degree','count')
    deg.dist.act.neg.1[,1] <- as.numeric(as.character(deg.dist.act.neg.1[,1]))
    deg.dist.act.neg.1[,2] <- as.numeric(deg.dist.act.neg.1[,2])

    deg.dist.env.pos <- deg_tm(net.int.pos[,c(1,2,3)], measure="output", nn = nrow(mat.int))[,2]
    deg.dist.env.pos.1 <- as.data.frame(table(deg.dist.env.pos))
    colnames(deg.dist.env.pos.1) <- c('degree','count')
    deg.dist.env.pos.1[,1] <- as.numeric(as.character(deg.dist.env.pos.1[,1]))
    deg.dist.env.pos.1[,2] <- as.numeric(deg.dist.env.pos.1[,2])

    deg.dist.env.neg <- deg_tm(net.int.neg[,c(1,2,3)], measure="output", nn = nrow(mat.int))[,2]
    deg.dist.env.neg.1 <- as.data.frame(table(deg.dist.env.neg))
    colnames(deg.dist.env.neg.1) <- c('degree','count')
    deg.dist.env.neg.1[,1] <- as.numeric(as.character(deg.dist.env.neg.1[,1]))
    deg.dist.env.neg.1[,2] <- as.numeric(deg.dist.env.neg.1[,2])

    deg.hist.act.pos <- hist(deg.dist.act.pos, breaks = seq(0, max(c(max(deg.dist.act.pos.1[,1]),max(deg.dist.act.neg.1[,1])))+3, by = 3), plot = F)
    deg.hist.act.pos$counts <- deg.hist.act.pos$counts/sum(deg.hist.act.pos$counts)
    deg.hist.act.neg <- hist(deg.dist.act.neg, breaks = seq(0, max(c(max(deg.dist.act.pos.1[,1]),max(deg.dist.act.neg.1[,1])))+3, by = 3), plot = F)
    deg.hist.act.neg$counts <- deg.hist.act.neg$counts/sum(deg.hist.act.neg$counts)
    deg.hist.env.pos <- hist(deg.dist.env.pos, breaks = seq(0, max(c(max(deg.dist.env.pos.1[,1]),max(deg.dist.env.neg.1[,1])))+3, by = 3), plot = F)
    deg.hist.env.pos$counts <- deg.hist.env.pos$counts/sum(deg.hist.env.pos$counts)
    deg.hist.env.neg <- hist(deg.dist.env.neg, breaks = seq(0, max(c(max(deg.dist.env.pos.1[,1]),max(deg.dist.env.neg.1[,1])))+3, by = 3), plot = F)
    deg.hist.env.neg$counts <- deg.hist.env.neg$counts/sum(deg.hist.env.neg$counts)

    ## statistical moments of degree distribution
    # standard deviation
    sd.deg.act.pos <- round(sd(as.numeric(deg[2:(ncol(mat.int)+1),2])), 1)
    sd.deg.act.neg <- round(sd(as.numeric(deg[2:(ncol(mat.int)+1),3])), 1)
    sd.deg.env.comp.neg <- round(sd(as.numeric(deg[2:(nrow(mat.int)+1),8])), 1)
    sd.deg.env.comp.pos <- round(sd(as.numeric(deg[2:(nrow(mat.int)+1),9])), 1)

    # skewness
    skew.deg.act.pos <- round(skewness(as.numeric(deg[2:(ncol(mat.int)+1),2])), 1)
    skew.deg.act.neg <- round(skewness(as.numeric(deg[2:(ncol(mat.int)+1),3])), 1)
    skew.deg.env.comp.neg <- round(skewness(as.numeric(deg[2:(nrow(mat.int)+1),8])), 1)
    skew.deg.env.comp.pos <- round(skewness(as.numeric(deg[2:(nrow(mat.int)+1),9])), 1)

    # plotting the histograms
    op <- par(yaxs = "i", xaxs = "i", mfrow = c(2,2), oma = c(1,1,1,0) + 0.1, mar = c(4,4,2,1) + 0.1)
    plot(deg.hist.act.pos, col = "green", main = "", xlab = expression("k2"["j"]), xaxt = "n", yaxt = "n",
         ylim = c(0, 1),
         xlim = c(0, max(deg.hist.act.pos$breaks)), cex.lab=1.4, ylab = "Relative frequency")
    axis(1, at = deg.hist.act.pos$breaks, cex.axis = 1.2)
    axis(2, at = seq(0, 1, by = 0.2), cex.axis = 1.2, las = 1)
    box()
    mtext(bquote(mu == .(mean.deg2.act.pos)), 2, line = -22, padj=-11.6, las = 1, font = 1, cex = 1)
    mtext(bquote(sigma == .(sd.deg.act.pos)), 2, line = -28, padj=-15, las = 1, font = 1, cex = 1)
    mtext(bquote(upsilon == .(skew.deg.act.pos)), 2, line = -34, padj=-15, las = 1, font = 1, cex = 1)

    plot(deg.hist.act.neg, col = "red", main = "", xlab = expression("k2"["j"]), xaxt = "n", yaxt = "n",
         ylim = c(0, 1),
         xlim = c(0, max(deg.hist.act.neg$breaks)), cex.lab=1.4, ylab = "Relative frequency")
    axis(1, at = deg.hist.act.neg$breaks, cex.axis = 1.2)
    axis(2, at = seq(0, 1, by = 0.2), cex.axis = 1.2, las = 1)
    box()
    mtext(bquote(mu == .(mean.deg2.act.neg)), 2, line = -22, padj=-11.6, las = 1, font = 1, cex = 1)
    mtext(bquote(sigma == .(sd.deg.act.neg)), 2, line = -28, padj=-15, las = 1, font = 1, cex = 1)
    mtext(bquote(upsilon == .(skew.deg.act.neg)), 2, line = -34, padj=-15, las = 1, font = 1, cex = 1)

    plot(deg.hist.env.pos, col = "green", main = "", xlab = expression("k2"["i"]), xaxt = "n", yaxt = "n",
         ylim = c(0, 1),
         xlim = c(0, max(deg.hist.env.pos$breaks)), cex.lab=1.4, ylab = "Relative frequency")
    axis(1, at = deg.hist.env.pos$breaks, cex.axis = 1.2)
    axis(2, at = seq(0, 1, by = 0.2), cex.axis = 1.2, las = 1)
    box()
    mtext(bquote(mu == .(mean.deg2.env.comp.pos)), 2, line = -22, padj=-11.6, las = 1, font = 1, cex = 1)
    mtext(bquote(sigma == .(sd.deg.env.comp.pos)), 2, line = -28, padj=-15, las = 1, font = 1, cex = 1)
    mtext(bquote(upsilon == .(skew.deg.env.comp.pos)), 2, line = -34, padj=-15, las = 1, font = 1, cex = 1)

    plot(deg.hist.env.neg, col = "red", main = "", xlab = expression("k2"["i"]), xaxt = "n", yaxt = "n",
         ylim = c(0, 1),
         xlim = c(0, max(deg.hist.env.neg$breaks)), cex.lab=1.4, ylab = "Relative frequency")
    axis(1, at = deg.hist.env.neg$breaks, cex.axis = 1.2)
    axis(2, at = seq(0, 1, by = 0.2), cex.axis = 1.2, las = 1)
    box()
    mtext(bquote(mu == .(mean.deg2.env.comp.neg)), 2, line = -22, padj=-11.6, las = 1, font = 1, cex = 1)
    mtext(bquote(sigma == .(sd.deg.env.comp.neg)), 2, line = -28, padj=-15, las = 1, font = 1, cex = 1)
    mtext(bquote(upsilon == .(skew.deg.env.comp.neg)), 2, line = -34, padj=-15, las = 1, font = 1, cex = 1)
    par(op)

    ## plotting the bipartite network
    # category colors for vertices
    V(ig.int)$color <- c(vcat.env.comps, vcat.actions)
    V(ig.int.row.pos)$color <- vcat.env.comps
    V(ig.int.row.neg)$color <- vcat.env.comps
    V(ig.int.col.pos)$color <- vcat.actions
    V(ig.int.col.neg)$color <- vcat.actions
    cl <- c("steelblue", "tan3", "orchid", "yellow",
            "slategray1", "lightcyan", "tan", "thistle", "linen")
    cl <- cl[1:length(unique(c(vcat.env.comps, vcat.actions)))]
    for (i in 1:length(cl)){
      V(ig.int)$color[V(ig.int)$color==i] <- cl[i]
      V(ig.int.row.pos)$color[V(ig.int.row.pos)$color==i] <- cl[i]
      V(ig.int.row.neg)$color[V(ig.int.row.neg)$color==i] <- cl[i]
      V(ig.int.col.pos)$color[V(ig.int.col.pos)$color==i] <- cl[i]
      V(ig.int.col.neg)$color[V(ig.int.col.neg)$color==i] <- cl[i]
    }

    # layout of bipartite network
    E(ig.int)$color <- ifelse(E(ig.int)$weight > 0, "green", "red")
    E(ig.int)$weight <- abs(E(ig.int)$weight)
    plt.x <- c((nrow(mat.int)+1):2, seq(nrow(mat.int)+1,2, length.out=ncol(mat.int)))
    plt.y <- c(rep(2,nrow(mat.int)),rep(4,ncol(mat.int)))
    lay <- as.matrix(cbind(plt.x,plt.y))
    shapes <- c("circle","square")
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot(ig.int, vertex.shape=shapes[V(ig.int)$type+1], vertex.label.family="Helvetica", vertex.label.font=4,
         vertex.size=9, vertex.label.color="black", vertex.label.cex=0.7, edge.width=E(ig.int)$weight, layout=lay)
    mtext("Actions", 3, 0)
    mtext("Environmental Components", 1, 0)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(19, length(cl)), inset=c(-0.355,0.1), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    legend("topright", legend = rep(" / ", length(cl)),
           col = unique(cl), pch = rep(15, length(cl)), inset=c(-0.1,0.1), bty = "n", y.intersp=0.9,
           x.intersp=0.3, xjust=0, yjust=0)
    par(op)

    ## plotting the one-mode-projections
    # environmental components
    set.seed(42)
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot.igraph(ig.int.row.pos, vertex.size=9,
                vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label.font=4, vertex.label.color="black",
                edge.width=sqrt(E(ig.int.row.pos)$weight), edge.color = "green", curved=F)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(19, length(cl)), inset=c(-0.35,0.05), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    par(op)

    set.seed(42)
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot.igraph(ig.int.row.neg, vertex.size=9,
                vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label.font=4, vertex.label.color="black",
                edge.width=sqrt(E(ig.int.row.neg)$weight), edge.color = "red", curved=F)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(19, length(cl)), inset=c(-0.35,0.05), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    par(op)

    # actions
    set.seed(42)
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot.igraph(ig.int.col.pos, vertex.shape = "square", vertex.size=9,
                vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label.font=4, vertex.label.color="black",
                edge.width=sqrt(E(ig.int.col.pos)$weight), edge.color = "green", curved=T)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(15, length(cl)), inset=c(-0.35,0.05), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    par(op)

    set.seed(42)
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot.igraph(ig.int.col.neg, vertex.shape = "square", vertex.size=9,
                vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label.font=4, vertex.label.color="black",
                edge.width=sqrt(E(ig.int.col.neg)$weight), edge.color = "red", curved=F)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(15, length(cl)), inset=c(-0.35,0.05), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    par(op)

    assign("bas.stats", bas.stats, envir = .GlobalEnv)
    assign("deg", deg, envir = .GlobalEnv)

  } else if (negative.values == FALSE){

    # short numeric labels for plotting
    rownames(mat.int) <- seq(1,nrow(mat.int),1)
    colnames(mat.int) <- seq(1,ncol(mat.int),1)

    # create igraph and tnet objects from weighted interaction matrix
    ig.int <- graph_from_incidence_matrix(mat.int, weighted = T)
    if (is_bipartite(ig.int) == FALSE)
      stop("Network is not bipartite!")
    net.int <- as.tnet(mat.int, type = "weighted two-mode tnet")

    ## using tnet for one-mode-projection and to create igraph object
    ## multi-edges are collapsed to single edge where the sum of the weights of the multi-edges is assigned
    ## create igraph and tnet objects for one-mode-projection from weighted interaction matrix
    net.int.row <- projecting_tm(net.int, method = "sum")
    net.int.col <- projecting_tm(net.int[,c(2,1,3)], method = "sum")
    ig.int.row <- graph_from_edgelist(as.matrix(net.int.row[,c(1,2)]), directed = F)
    E(ig.int.row)$weight <- net.int.row[,3]
    ig.int.row  <- simplify(ig.int.row, edge.attr.comb="sum")
    ig.int.row <- check_vertices(ig.int.row, no_v=nrow(mat.int))
    ig.int.col <- graph_from_edgelist(as.matrix(net.int.col[,c(1,2)]), directed = F)
    E(ig.int.col)$weight <- net.int.col[,3]
    ig.int.col <- simplify(ig.int.col, edge.attr.comb="sum")
    ig.int.col <- check_vertices(ig.int.col, no_v=ncol(mat.int))

    ## data frame with basic statistics for two-mode and one-mode
    bas.stats <- as.data.frame(matrix(data = NA, nrow = 8, ncol = 4))
    colnames(bas.stats) <- c(' ', 'Two-mode',
                             ' ', 'One-mode')
    rownames(bas.stats) <- c(' ', 'Number_of_vertices', 'Number_of_edges','Network_asymmetry',
                             'Connectance', 'Mean_distance', 'Cluster_coefficient',  'Degree_correlation_coefficient')
    bas.stats[1,] <- c('Actions', 'Environmental components', 'Actions', 'Environmental components')

    ## number of vertices
    bas.stats[2,1]  <- ncol(mat.int)
    bas.stats[2,2]  <- nrow(mat.int)
    bas.stats[2,3]  <- ncol(mat.int)
    bas.stats[2,4]  <- nrow(mat.int)

    ## number of edges
    bas.stats[3,1]  <- gsize(ig.int)
    bas.stats[3,2]  <- gsize(ig.int)
    bas.stats[3,3]  <- gsize(ig.int.col)
    bas.stats[3,4]  <- gsize(ig.int.row)

    ## network asymmetry
    bas.stats[4,1] <- (nrow(mat.int)-ncol(mat.int))/(nrow(mat.int)+ncol(mat.int))
    bas.stats[4,3] <- 0

    ## connectance
    bas.stats[5,1] <- gsize(ig.int)/(nrow(mat.int)*ncol(mat.int))
    bas.stats[5,3] <- gsize(ig.int.col)/(ncol(mat.int)*ncol(mat.int))
    bas.stats[5,4] <- gsize(ig.int.row)/(nrow(mat.int)*nrow(mat.int))

    ## shortest paths and mean distance
    # two-mode
    dist.int.row.2 <- distance_tm(net.int[,c(1,2,3)], projection.method="sum")
    dist.int.col.2 <- distance_tm(net.int[,c(2,1,3)], projection.method="sum")

    bas.stats[6,1] <- mean(dist.int.col.2, na.rm = T)
    bas.stats[6,2] <- mean(dist.int.row.2, na.rm = T)

    # one-mode
    dist.int.row.1 <- distance_w(net.int.row, directed=F)
    dist.int.col.1 <- distance_w(net.int.col, directed=F)

    bas.stats[6,3] <- mean(dist.int.col.1, na.rm = T)
    bas.stats[6,4] <- mean(dist.int.row.1, na.rm = T)

    ## clustering
    # global cluster coefficient
    # two mode
    # in which way are weights considered?
    bas.stats[7,1] <- clustering_tm(net.int[,c(2,1,3)])[3] #geometric mean
    bas.stats[7,2] <- clustering_tm(net.int[,c(1,2,3)])[3]

    # one mode
    bas.stats[7,3] <- clustering_w(net.int.col, measure = "gm") #geometric mean
    bas.stats[7,4] <- clustering_w(net.int.row, measure = "gm")

    ## nominal assortativity by categories
    # check if category vectors are empty
    if (is.null(vcat) == TRUE){
      stop("Vector with the labels for the vertice categories is empty!")
    }
    if (is.null(vcat.env.comps) == TRUE){
      stop("Vector with categories for environmental components is empty!")
    }
    if (is.null(vcat.actions) == TRUE){
      stop("Vector with categories for actions is empty!")
    }

    # check if category labels and the number of unique numerical categories have the same dimension
    vcat.num <- unique(c(vcat.env.comps, vcat.actions))
    if (length(vcat) != length(vcat.num)){
      stop("Vectors with category labels and numerical categories do not have the same dimension!")
    }

    # assign categories as vertice attributes
    V(ig.int)$cat <- c(vcat.env.comps, vcat.actions)
    V(ig.int.row)$cat <- vcat.env.comps
    V(ig.int.col)$cat <- vcat.actions
    # two-mode
    bas.stats[8,1] <- assortativity_nominal(ig.int, types=V(ig.int)$cat, directed = F)

    # one-mode
    bas.stats[8,3] <- assortativity_nominal(ig.int.col, types=V(ig.int.col)$cat, directed = F)
    bas.stats[8,4] <- assortativity_nominal(ig.int.row, types=V(ig.int.row)$cat, directed = F)

    bas.stats[is.na(bas.stats)] <- " "
    rownames(bas.stats) <- str_replace_all(rownames(bas.stats), "_", " ")

    ## node degree & mean degree
    deg <- matrix(data = NA, nrow = nrow(mat.int)+1, ncol = 6)
    colnames(deg) <- c(' ', 'Two-mode', 'One-mode', 'One-mode', 'Two-mode', ' ')
    deg[1:(ncol(mat.int)), 1] <- seq(1, ncol(mat.int), 1)
    deg[1:(nrow(mat.int)), 6] <- seq(1, nrow(mat.int), 1)

    # two-mode
    deg[1:(ncol(mat.int)),2] <- deg_tm(net.int[,c(2,1,3)], measure="output", ncol(mat.int))[,2]
    mean.deg2.act <- round(sum(as.numeric(deg[1:(ncol(mat.int)),2]))/ncol(mat.int))
    deg[(ncol(mat.int)+1),2] <- mean.deg2.act

    deg[1:(nrow(mat.int)),5] <- deg_tm(net.int, measure="output", nrow(mat.int))[,2]
    mean.deg2.env.comp <- round(sum(as.numeric(deg[1:(nrow(mat.int)),5]))/nrow(mat.int))
    deg[(nrow(mat.int)+1),5] <- mean.deg2.env.comp

    # one-mode
    deg[1:(ncol(mat.int)),3] <- deg_w(net.int.col, measure="output", nn=ncol(mat.int))[,2]
    deg[(ncol(mat.int)+1),3] <- sum(as.numeric(deg[1:(ncol(mat.int)),3]))/ncol(mat.int)

    deg[1:(nrow(mat.int)),4] <- deg_w(net.int.row, measure="output", nn=nrow(mat.int))[,2]
    deg[(nrow(mat.int)+1),4] <- sum(as.numeric(deg[1:(nrow(mat.int)),4]))/nrow(mat.int)

    deg[is.na(deg)] <- " "

    ## degree distribution with relative frequencies
    deg.dist.act <- deg_tm(net.int[,c(2,1,3)], measure="output", ncol(mat.int))[,2]
    deg.dist.act.1 <- as.data.frame(table(deg.dist.act))
    colnames(deg.dist.act.1) <- c('degree','count')
    deg.dist.act.1[,1] <- as.numeric(as.character(deg.dist.act.1[,1]))
    deg.dist.act.1[,2] <- as.numeric(deg.dist.act.1[,2])

    deg.dist.env <- deg_tm(net.int[,c(1,2,3)], measure="output", nrow(mat.int))[,2]
    deg.dist.env.1 <- as.data.frame(table(deg.dist.env))
    colnames(deg.dist.env.1) <- c('degree','count')
    deg.dist.env.1[,1] <- as.numeric(as.character(deg.dist.env.1[,1]))
    deg.dist.env.1[,2] <- as.numeric(deg.dist.env.1[,2])

    deg.hist.act <- hist(deg.dist.act, breaks = seq(0, max(deg.dist.act.1[,1])+1, by = 3), plot = F)
    deg.hist.act$counts <- deg.hist.act$counts/sum(deg.hist.act$counts)
    deg.hist.env <- hist(deg.dist.env, breaks = seq(0, max(deg.dist.env.1[,1])+1, by = 3), plot = F)
    deg.hist.env$counts <- deg.hist.env$counts/sum(deg.hist.env$counts)

    ## statistical moments of histogram
    # standard deviation
    sd.deg.act <- round(sd(as.numeric(deg[2:(ncol(mat.int)+1),2])), 1)
    sd.deg.env.comp <- round(sd(as.numeric(deg[2:(nrow(mat.int)+1),5])), 1)

    # skewness
    skew.deg.act <- round(skewness(as.numeric(deg[2:(ncol(mat.int)+1),2])), 1)
    skew.deg.env.comp <- round(skewness(as.numeric(deg[2:(nrow(mat.int)+1),5])), 1)

    # plotting the histograms
    op <- par(yaxs = "i", xaxs = "i", mfrow = c(1,2), oma = c(1,1,1,0) + 0.1, mar = c(4,4,2,1) + 0.1)
    plot(deg.hist.act, col = "grey",
         main = "", xlab = expression("k2"["j"]), xaxt = "n", yaxt = "n", ylim = c(0, 1),
         xlim = c(0, max(deg.hist.act$breaks)), cex.lab=1.4, ylab = "Relative frequency")
    axis(1, at = deg.hist.act$breaks, cex.axis = 1.2)
    axis(2, at = seq(0, 1, 0.2), cex.axis = 1.2, las = 1)
    box()
    mtext(bquote(mu == .(mean.deg2.act)), 2, line = -18, padj=-23.4, las = 1, font = 1, cex = 1.2)
    mtext(bquote(sigma == .(sd.deg.act)), 2, line = -24, padj=-30, las = 1, font = 1, cex = 1.2)
    mtext(bquote(upsilon == .(skew.deg.act)), 2, line = -30, padj=-30, las = 1, font = 1, cex = 1.2)

    plot(deg.hist.env, col = "grey",
         main = "", xlab = expression("k2"["i"]), xaxt = "n", yaxt = "n", ylim = c(0, 1),
         xlim = c(0, max(deg.hist.env$breaks)), cex.lab=1.4, ylab = "Relative frequency")
    axis(1, at = deg.hist.env$breaks, cex.axis = 1.2)
    axis(2, at = seq(0, 1, 0.2), cex.axis = 1.2, las = 1)
    box()
    mtext(bquote(mu == .(mean.deg2.env.comp)), 2, line = -18, padj=-23.7, las = 1, font = 1, cex = 1.2)
    mtext(bquote(sigma == .(sd.deg.env.comp)), 2, line = -24, padj=-30, las = 1, font = 1, cex = 1.2)
    mtext(bquote(upsilon == .(skew.deg.env.comp)), 2, line = -30, padj=-30, las = 1, font = 1, cex = 1.2)
    par(op)

    ## plotting the bipartite network
    # category colors for vertices
    V(ig.int)$color <- c(vcat.env.comps, vcat.actions)
    V(ig.int.row)$color <- vcat.env.comps
    V(ig.int.col)$color <- vcat.actions
    cl <- c("steelblue", "tan3", "orchid", "yellow",
            "slategray1", "lightcyan", "tan", "thistle", "linen")
    cl <- cl[1:length(unique(c(vcat.env.comps, vcat.actions)))]
    for (i in 1:length(cl)){
      V(ig.int)$color[V(ig.int)$color==i] <- cl[i]
      V(ig.int.row)$color[V(ig.int.row)$color==i] <- cl[i]
      V(ig.int.col)$color[V(ig.int.col)$color==i] <- cl[i]
    }

    # layout of bipartite network
    plt.x <- c((nrow(mat.int)+1):2, seq(nrow(mat.int)+1,2, length.out=ncol(mat.int)))
    plt.y <- c(rep(2,nrow(mat.int)),rep(4,ncol(mat.int)))
    lay <- as.matrix(cbind(plt.x,plt.y))
    shapes <- c("circle","square")
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot(ig.int, vertex.shape=shapes[V(ig.int)$type+1], vertex.label.family="Helvetica", vertex.label.font=4,
         vertex.size=9, vertex.label.color="black", vertex.label.cex=0.7, edge.width=(E(ig.int)$weight)^2, edge.color="grey", layout=lay)
    mtext("Actions", 3, 0)
    mtext("Environmental Components", 1, 0)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(19, length(cl)), inset=c(-0.355,0.1), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    legend("topright", legend = rep(" / ", length(cl)),
           col = unique(cl), pch = rep(15, length(cl)), inset=c(-0.1,0.1), bty = "n", y.intersp=0.9,
           x.intersp=0.3, xjust=0, yjust=0)
    par(op)

    ## plotting the one-mode-projection
    # environmental components
    set.seed(42)
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot.igraph(ig.int.row, vertex.size=9,
                vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label.font=4, vertex.label.color="black",
                edge.width=sqrt(E(ig.int.row)$weight), edge.color = "grey", curved=F)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(19, length(cl)), inset=c(-0.35,0.05), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    par(op)

    # actions
    set.seed(42)
    op <- par(oma = c(1,1,1,1) + 0.1, mar = c(1,0,1,5) + 0.1)
    plot.igraph(ig.int.col, vertex.shape = "square", vertex.size=9,
                vertex.label.cex=0.7, vertex.label.family="Helvetica", vertex.label.font=4, vertex.label.color="black",
                edge.width=sqrt(E(ig.int.col)$weight), edge.color = "grey", curved=T)
    legend("topright", legend = vcat,
           col = unique(cl), pch = rep(15, length(cl)), inset=c(-0.35,0.05), bty = "n", y.intersp=0.9,
           x.intersp=0.8, xjust=0, yjust=0)
    par(op)

    assign("bas.stats", bas.stats, envir = .GlobalEnv)
    assign("deg", deg, envir = .GlobalEnv)
  }
}
