suppressPackageStartupMessages({library(data.table); library(dplyr); library(purrr); library(parallel)})
set.seed(20261023)
source("AcadianGY_12.3.6.r")
spp <- c(PIRU="RS",TSCA="EH",ACRU="RM",THOC="WC",ABBA="BF",PIST="WP",BEPA="PB",BEPO="GB",BEAL="YB",POGR="BT",POTR="QA",FAGR="AB",UNK="OH")
JENK <- list(spruce=c(-2.0773,2.3323),fir_hem=c(-2.5384,2.4814),pine=c(-2.5356,2.4349),cedar=c(-2.0336,2.2592),
             aspen=c(-2.2094,2.3867),smb=c(-1.9123,2.3651),mixhw=c(-2.4800,2.4835),hmo=c(-2.0127,2.4342))
grp <- c(RS="spruce",BS="spruce",WS="spruce",EH="fir_hem",BF="fir_hem",WP="pine",WC="cedar",RM="smb",PB="smb",GB="smb",
         YB="smb",BT="aspen",QA="aspen",AB="hmo",SM="hmo",OH="mixhw")
jenk_kg <- function(sp,dbh){g<-grp[sp];g[is.na(g)]<-"mixhw";b<-do.call(rbind,JENK[g]);exp(b[,1]+b[,2]*log(dbh))}
td <- fread("../trees.tab")
td[, PLOT := sub("^([A-Z]+_[0-9]+).*", "\\1", Plot)]
td <- td[grepl("^HL", PLOT) & (is.na(Snag) | Snag == "") & DBH >= 10]
td <- td[, .SD[Year == max(Year)], by = PLOT]
cat("plots", uniqueN(td$PLOT), "trees", nrow(td), "\n")
tree0 <- data.frame(STAND=1, YEAR=2025, PLOT=td$PLOT, TREE=seq_len(nrow(td)), SP=unname(spp[td$Species]), DBH=td$DBH,
                    HT=NA_real_, HCB=NA_real_, EXPF=10000/250, Form=NA, Risk=NA, Region="ME", dDBH.mult=1, dHt.mult=1,
                    mort.mult=1, max.dbh=999, max.height=999, stringsAsFactors=FALSE)
tree0$SP[is.na(tree0$SP)] <- "OH"
ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=2.54, CutPoint=0.95, usedHTCap=TRUE,
            rtnVars=c("STAND","YEAR","PLOT","TREE","SP","DBH","HT","HCB","EXPF","pHT","pHCB","Form","Risk",
                      "dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height","RD","SDImax"))
run_plot <- function(tr, csi, ncyc = 50) {
  stand <- list(CSI = csi, ELEV = 60); cur <- tr; out <- list()
  for (i in 1:ncyc) {
    prev <- cur
    nxt <- tryCatch(AcadianGYOneStand(cur, stand = stand, ops = ops), error = function(e) NULL)
    if (is.null(nxt)) break
    nxt$YEAR <- 2025 + i
    for (v in c("dDBH.mult","dHt.mult","mort.mult")) nxt[[v]][is.na(nxt[[v]])] <- 1
    for (v in c("max.dbh","max.height")) nxt[[v]][is.na(nxt[[v]])] <- 999
    m <- merge(prev[,c("TREE","SP","DBH","EXPF")], nxt[,c("TREE","DBH","EXPF")], by="TREE", suffixes=c("0","1"))
    a0 <- jenk_kg(m$SP, m$DBH0); a1 <- jenk_kg(m$SP, m$DBH1)
    growth <- sum((a1 - a0) * m$EXPF1) / 1000           # Mg/ha, survivors
    mort   <- sum(a1 * (m$EXPF0 - m$EXPF1)) / 1000       # Mg/ha lost to mortality
    agb    <- sum(jenk_kg(nxt$SP, nxt$DBH) * nxt$EXPF) / 1000
    ba     <- sum(pi * (nxt$DBH/200)^2 * nxt$EXPF)
    out[[i]] <- data.table(PLOT=tr$PLOT[1], CSI=csi, year=2025+i, agb_J=agb, ba=ba, tph=sum(nxt$EXPF),
                           gross_growth=growth, mortality=mort,
                           RD=if ("RD" %in% names(nxt)) mean(nxt$RD, na.rm=TRUE) else NA_real_,
                           SDImax=if ("SDImax" %in% names(nxt)) mean(nxt$SDImax, na.rm=TRUE) else NA_real_)
    cur <- nxt
  }
  rbindlist(out)
}
plots <- split(tree0, tree0$PLOT)
res <- rbindlist(lapply(c(10, 12, 14), function(csi)
  rbindlist(mclapply(plots, run_plot, csi = csi, mc.cores = 12))))
a0 <- tree0 %>% group_by(PLOT) %>% summarise(agb0 = sum(jenk_kg(SP, DBH) * EXPF) / 1000)
fwrite(res, "csi_trajectories.csv")
s <- as.data.table(res)[year %in% c(2033, 2075), .(agb = mean(agb_J), growth_cum = NA_real_, RD = mean(RD, na.rm=TRUE),
                                                 SDImax = mean(SDImax, na.rm=TRUE), n = .N), by = .(CSI, year)]
g <- as.data.table(res)[year <= 2033, .(gross8 = sum(gross_growth), mort8 = sum(mortality)), by = .(CSI, PLOT)][
  , .(gross_8yr = mean(gross8), mortality_8yr = mean(mort8)), by = CSI]
cat("initial mean AGB (Jenkins biomass Mg/ha):", mean(a0$agb0), "\n")
print(s); print(g)
fwrite(merge(s, g, by = "CSI"), "csi_summary.csv")
