library(openxlsx)

## load enrichment function
source("./EnrichrFunctions.R")

OE75.up.degs <- read.xlsx("./APM_OE75_DEGS_5donor.xlsx", sheet=1)
OE75.down.degs <- read.xlsx("./APM_OE75_DEGS_5donor.xlsx", sheet=2)

OE7_5.up.degs <- read.xlsx("./APM_OE7.5_DEGS_12donor.xlsx", sheet=1)
OE7_5.down.degs <- read.xlsx("./APM_OE7.5_DEGS_12donor.xlsx", sheet=2)

DEGs <- list()

DEGs[[1]] <- OE75.up.degs
DEGs[[1]]$comparison <- "OE75_up"
DEGs[[2]] <- OE75.down.degs
DEGs[[2]]$comparison <- "OE75_down"
DEGs[[3]] <- OE7_5.up.degs
DEGs[[3]]$comparison <- "OE7.5_up"
DEGs[[4]] <- OE7_5.down.degs
DEGs[[4]]$comparison <- "OE7.5_down"

DEGs.df <- do.call(rbind, DEGs)

# call enrichr
doEnrichOneAtATime(DEGs.df, dataset="OE")

