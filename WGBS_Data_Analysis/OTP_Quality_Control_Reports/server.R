## F U N C T I O N S ###################################################################################################

gg_theme_blank <- function(...) {
	theme.elements <- as.character(list(...))
	names(theme.elements) <- theme.elements
	theme.elements <- lapply(theme.elements, function(x) { element_blank() })
	do.call(theme, theme.elements)
}

########################################################################################################################

gg_nodata <- function(txt = "No data available") {
	ggplot(data.frame(x = 1, y = 1, z = txt), aes_string(x = "x", y = "y", label = "z")) +
		geom_text(color = "grey50") + gg_theme_blank('axis.line', 'axis.title', 'axis.text', 'axis.ticks') +
		gg_theme_blank('panel.border', 'panel.grid', 'panel.background', 'plot.background')
}

########################################################################################################################

plot.samplestats <- function(dframe = tbls[["samples"]]) {
	cnames <- c("State", "Metrics")
	tbl <- tapply(1:nrow(dframe), dframe[, cnames], length)
	tbl[is.na(tbl)] <- 0L
	tbl <- data.frame(
		'x' = factor(rep(colnames(tbl), each = nrow(tbl)), levels = colnames(tbl)),
		'y' = factor(rep(rownames(tbl), ncol(tbl)), levels = rev(rownames(tbl))),
		'z' = as.vector(tbl))
	ggplot(tbl, aes_string(x = 'x', y = 'y', fill = 'z', label = 'z')) + geom_tile(color = "#FFFFFF") + geom_text() +
		scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
		scale_fill_continuous(low = "#FFFFFF", high = "#808080") +
		labs(x = cnames[2], y = cnames[1], z = 'Samples') +
		gg_theme_blank('axis.ticks', 'axis.line', 'panel.border', 'panel.grid', 'plot.background') +
		theme(legend.position = "none")
}

########################################################################################################################

plot.insertsizes <- function(sample.id, dframe = tbls[["insertsizes"]]) {
	tbl <- dframe[dframe$`Sample` == sample.id, c("Size", "Frequency")]
	tbl<- tbl[which(tbl$Frequency >= max(tbl$Frequency)/100),]
	if (nrow(tbl) == 0) {
		return(gg_nodata())
	}
	
	ggplot(tbl, aes_string(x = 'Size', y = 'Frequency')) +
		geom_col(width = 1, colour = NA, fill = "#000000") +
		attr(dframe, "in.scales")[['x']] + attr(dframe, "fr.scales")[['y']] +
		labs(x = "Insert size (bp)", y = "Frequency (thousand)") +
	  xlim(0,max(tbl$Size)) +
	  ylim(0,tbl$Frequency[which(tbl$Frequency == max(tbl$Frequency))])
}

########################################################################################################################

plot.meth <- function(sample.id, chrom.name, cytosine.context, dframe) {
	if ("pos" %in% colnames(dframe)) {
		cnames <- "pos"
		xlabel <- "Position in read"
	} else {
		cnames <- "baseQ"
		xlabel <- "Quality score"
	}
	cnames <- c("mate", cnames, paste0(cytosine.context, ".F"))
	tbl <- dframe[dframe$`Sample` == sample.id & dframe$`Chromosome` == chrom.name, cnames]
	tbl <- tbl[!is.na(tbl[, 3]), , drop = FALSE]
	ggplot(tbl, aes_string(x = colnames(tbl)[2], y = colnames(tbl)[3], color = 'mate')) + geom_line() + geom_point() +
		attr(dframe, "sc.scales")[['x']] + attr(dframe, paste0(cytosine.context, ".scales"))[['y']] +
		labs(x = xlabel, y = 'Methylated cytosines (%)')
}

########################################################################################################################

plot.coverage <- function(sample.id, chrom.name, cytosine.context, dframe = tbls[["coverage"]]) {
	tbl <- dframe[dframe$`Sample` == sample.id & dframe$`Chromosome` == chrom.name, c("cov", cytosine.context)]
	tbl<- tbl[which(tbl[cytosine.context] >= max(tbl[cytosine.context])/100),]
	yscaling <- ifelse(cytosine.context == "CG", "thousand", "million")
	pp <- ggplot(tbl, aes_string(x = 'cov', y = cytosine.context)) +
		geom_col(width = 1, colour = NA, fill = "#000000") +
		attr(dframe, "co.scales")[['x']] + attr(dframe, paste0(cytosine.context, ".scales"))[['y']] +
		labs(x = "Coverage", y = paste0("Frequency (", yscaling, ")")) +
		theme(plot.margin = unit(0.1 + c(0, 0.3, 0, 0), "in")) +
	  xlim(0,max(tbl$cov)) +
	  ylim(0,tbl$CG[which(tbl$CG == max(tbl$CG))])
	print(pp)
}

## S E R V E R #########################################################################################################

shinyServer(function(input, output) {

	observe({
		if (input$selPlot == PLOT.TYPES[1]) {
			output$plotMain <- renderPlot(plot.samplestats())
		} else if (input$selPlot == PLOT.TYPES[2]) {
			output$plotMain <- renderPlot(plot.insertsizes(input$selSample))
		} else if (input$selPlot == PLOT.TYPES[3]) {
			dframe <- tbls[["methylation vs. position"]]
			output$plotMain <- renderPlot(plot.meth(input$selSample, input$selChrom, input$selContext, dframe))
		} else if (input$selPlot == PLOT.TYPES[4]) {
			dframe <- tbls[["methylation vs. baseQ"]]
			output$plotMain <- renderPlot(plot.meth(input$selSample, input$selChrom, input$selContext, dframe))
		} else { # input$selPlot == PLOT.TYPES[5]
			output$plotMain <- renderPlot(plot.coverage(input$selSample, input$selChrom, input$selContext))
		}
	})
})
