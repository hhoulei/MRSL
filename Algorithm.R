
Open_path_adj <- function(MM, start0, end0, adj) {

    locname <- which(MM != 0, arr.ind = T)
    dagstatement <- "dag {"
    for (i in 1:nrow(locname)) {
        dagstatement <- paste0(dagstatement, "D", locname[i, 1], " -> ", "D", locname[i, 2], " ; ")
    }
    dagstatement <- substr(dagstatement, 1, nchar(dagstatement) - 3)
    dagstatement <- paste0(dagstatement, "}")

    g <- dagitty(dagstatement)

    start0 = paste0("D", start0)
    end0 = paste0("D", end0)

    if (adj == 2) {

        node_adj <- NULL

        pp <- paths(g, start0, end0, directed = F)  ##list the paths
        if (length(pp) != 0) {
            openpath <- pp$paths[pp$open]
            nodes_all <- strsplit(openpath, " ")

            if (length(nodes_all) != 0) {
                for (oi in 1:length(nodes_all)) {
                  # cat('ppR2_',oi,'\n')
                  node_once <- nodes_all[[oi]]
                  # cat('ppR3','\n')
                  if (length(node_once) > 3) {
                    node_loc <- seq(3, length(node_once) - 2, 2)
                    node_adj <- c(node_adj, node_once[node_loc])
                  }
                }
            }

        }


    } else if (adj == 3) {

        loc_adj0 <- adjustmentSets(g, start0, end0)

        if (length(loc_adj0) == 0) {
            node_adj <- NULL
        } else {
            node_adj <- loc_adj0[[1]]
        }


        pp2 = paths(g, start0, end0, directed = TRUE)$paths

        if (length(pp2) != 0) {

            nodes_all <- strsplit(pp2, " ")

            for (oi in 1:length(nodes_all)) {
                node_once <- nodes_all[[oi]]
                if (length(node_once) > 3) {
                  node_loc <- seq(3, length(node_once) - 2, 2)
                  node_adj <- c(node_adj, node_once[node_loc])
                }
            }

        }


    }


    if (length(node_adj) != 0) {
        node_adj <- unique(node_adj)
        node_adj <- as.numeric(substr(node_adj, 2, nchar(node_adj)))
        node_adj <- node_adj - 1
    }

    node_adj <- as.numeric(node_adj)

    return(node_adj)
}
