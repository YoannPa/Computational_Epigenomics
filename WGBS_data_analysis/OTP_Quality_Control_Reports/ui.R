## U I #################################################################################################################

shinyUI(fluidPage(

	titlePanel("OTP Quality Control Reports"),

	sidebarLayout(
		sidebarPanel(
			selectInput("selPlot", "Plot type", PLOT.TYPES, selectize = FALSE),
			conditionalPanel(paste0("input.selPlot != '", PLOT.TYPES[1], "'"),
				br(),
				selectInput("selSample", "Sample", levels(tbls$samples$Sample), selectize = FALSE)),
			conditionalPanel(paste0("input.selPlot == '", PLOT.TYPES[3], "' || input.selPlot == '", PLOT.TYPES[4], "'"),
				br(),
				selectInput("selChrom", "Chromosome", levels(tbls$coverage$Chromosome), selectize = FALSE)),
			conditionalPanel(paste0("input.selPlot == '", PLOT.TYPES[4], "'"),
				br(),
				selectInput("selContext", "Cytosine context", c("CG", "CH"), selectize = FALSE))
		),

		mainPanel(plotOutput("plotMain"))
)))
