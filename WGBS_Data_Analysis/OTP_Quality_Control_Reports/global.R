## L I B R A R I E S ###################################################################################################

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
theme_set(theme_bw())

## Load all data to be visualized
#tbls <- readRDS("Z:/Clarissa/Clarissa-WGBS-quality/otp_metrics_prostate.RDS")
tbls <- readRDS("/home/yoann/C010-Projects/Computational/Projects/2018-03-Hematopoietic-DNMT1/data/otp_metrics_twgbs.RDS")

## F U N C T I O N S ###################################################################################################

gg_scales <- function(x, xstep, xticks = xstep, add.ormore = FALSE) {
	xrange <- range(x, na.rm = TRUE)
	xrange <- c(floor(xrange[1] / xstep), ceiling(xrange[2] / xstep)) * xstep
	if (xrange[1] == 1) {
		xbreaks <- seq(0, xrange[2], by = xticks)
		xbreaks[1] <- 1
		xrange <- xrange + 0.5 * c(-1, 1)
	} else {
		xbreaks <- seq(xrange[1], xrange[2], by = xticks)
	}
	scale.params <- list(breaks = xbreaks, limits = xrange, expand = c(0, 0))
	if (add.ormore) {
		lbs <- as.character(xbreaks)
		lbs[length(lbs)] <- paste(lbs[length(lbs)], "or more")
		scale.params$labels <- lbs
	}
	list('x' = do.call(scale_x_continuous, scale.params), 'y' = do.call(scale_y_continuous, scale.params))
}

########################################################################################################################

preprocess.methylation.vs <- function(dframe, multiplicative.factor) {
	dframe$mate <- factor(dframe$mate, levels = 1:2)
	levels(dframe$mate) <- c("first", "second")
	dframe$CG.F <- dframe$CG.mC / (dframe$CG.mC + dframe$CG.C) * multiplicative.factor
	dframe$CH.F <- dframe$CH.mC / (dframe$CH.mC + dframe$CH.C) * multiplicative.factor
	dframe$CG.mC <- NULL
	dframe$CG.C <- NULL
	dframe$CH.mC <- NULL
	dframe$CH.C <- NULL
	dframe
}

## G L O B A L S #######################################################################################################

## Preprocess the table of insert sizes
tbls$insertsizes$Frequency <- tbls$insertsizes$Frequency / 1000
attr(tbls$insertsizes, "in.scales") <- gg_scales(tbls$insertsizes$Size, 10, 50)
attr(tbls$insertsizes, "fr.scales") <- gg_scales(tbls$insertsizes$Frequency, 1000)

## Preprocess the table of coverage
tbls$coverage$CG <- tbls$coverage$CG / 1000
tbls$coverage$CH <- tbls$coverage$CH / 1000000
attr(tbls$coverage, "co.scales") <- gg_scales(tbls$coverage$cov, 10, 50, add.ormore = TRUE)
attr(tbls$coverage, "CG.scales") <- gg_scales(tbls$coverage$CG, 10, 20)
attr(tbls$coverage, "CH.scales") <- gg_scales(tbls$coverage$CH, 0.1, 0.2)

## Preprocess the tables 'methylation vs. position' and 'methylation vs. baseQ'
tbls$`methylation vs. position` <- preprocess.methylation.vs(tbls$`methylation vs. position`, 100)
tbls$`methylation vs. baseQ` <- preprocess.methylation.vs(tbls$`methylation vs. baseQ`, 100)
attr(tbls$`methylation vs. position`, "sc.scales") <- gg_scales(tbls$`methylation vs. position`$pos, 1, 10)
attr(tbls$`methylation vs. position`, "CG.scales") <- gg_scales(tbls$`methylation vs. position`$CG.F, 1, 10)
attr(tbls$`methylation vs. position`, "CH.scales") <- gg_scales(tbls$`methylation vs. position`$CH.F, 1, 10)
attr(tbls$`methylation vs. baseQ`, "sc.scales") <- gg_scales(tbls$`methylation vs. baseQ`$baseQ, 1)
attr(tbls$`methylation vs. baseQ`, "CG.scales") <- gg_scales(tbls$`methylation vs. baseQ`$CG.F, 1, 10)
attr(tbls$`methylation vs. baseQ`, "CH.scales") <- gg_scales(tbls$`methylation vs. baseQ`$CH.F, 1, 10)

## Remove unused functions
rm(gg_scales, preprocess.methylation.vs)

PLOT.TYPES <- c(
	"summary" = "summary",
	"insert sizes" = "isizes",
	"methylation vs. position" = "methpos",
	"methylation vs. baseQ" = "methbas",
	"coverage" = "coverage")
