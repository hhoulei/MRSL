
Clump_LD <- function(nodex1, nodeadjust, bet1, l19) {

    nodex11 = nodex1 + 1
    nodeadjust1 = nodeadjust + 1
    l191 = l19 + 1
    # cat('l191=',length(l191),'\n') cat('nodex11=',length(nodex11),'\n')
    # cat('nodeadjust1=',length(nodeadjust1),'\n')

    nonadj <- c(1:ncol(bet1))[-c(nodex11, nodeadjust1)]
    # cat('nonadj=',length(nonadj),'\n')
    if (length(nonadj) == 1) {
        l20 <- which(bet1[, nonadj] == 0)
    } else {
        l20 <- which(apply(bet1[, nonadj], 1, sum) == 0)
    }
    # cat('l20=',length(l20),'\n')
    l21 <- l191[l191 %in% l20]
    # cat('l21=',length(l21),'\n')

    load("beta.Rdata")
    snpsname <- rownames(beta)[l21]

    # cat('snpsname=',length(snpsname),'\n')

    load("panel.Rdata")
    dt1 <- panel[panel$SNP %in% snpsname, ]
    cat("SNPs1=", nrow(dt1), "\n")
    dt2 <- clump_data(dt1)
    cat("SNPs2=", nrow(dt2), "\n")

    l22 <- which(rownames(beta) %in% dt2$SNP)
    l22 <- l22 - 1

    return(l22)

}
