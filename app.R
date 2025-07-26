library(admisc)
library(date)
library(dplyr)
library(fresh)
library(geoR)
library(ggplot2)
library(kableExtra)
library(knitr)
library(MASS)
library(mgcv)
library(modeest)
library(msm)
library(nlme)
library(qgam)
library(refineR)
library(reflimR)
library(rhandsontable)
library(scales)
library(shiny)
library(shinyalert)
library(shinybusy)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(snpar) # installation: "install.packages('devtools')" and "devtools::install_github('debinqiu/snpar')"
library(stringr)
library(tidykosmic) # installation: "install.packages('devtools')" and "devtools::install_github('divinenephron/tidykosmic')"

# Colors
uoltheme <- create_theme(
  adminlte_color(light_blue = "#003B73"),
  adminlte_sidebar(width = "250px", dark_bg = "#40668d", dark_hover_bg = "#40668d", dark_color = "#FFF"),
  adminlte_global(content_bg = "#EEE", box_bg = "#FFF", info_box_bg = "#FFF" )
)

buttoncolors1 <- "color: #FFF; background-color: #003B73; border-color: #FFF; width: 90%;"
buttoncolors2 <- "color: #000; background-color: #BBBBBB; border-color: #000;"

# Variables
malelist <- c('male', 'männlich', 'Mann', 'M', 'm')
femalelist <- c('female', 'weiblich', 'Frau', 'F', 'f', 'W', 'w')
sexlist <- c(malelist, femalelist, 'D', 'X')
tablesize <- 200000

# Path for TMC and TML library
if (Sys.info()["sysname"] == "Windows") {
  source("Z:/R/kc.uol.de/referenceranger/TMC.R")
  source("Z:/R/kc.uol.de/referenceranger/TML.R")
} else {
  source("/srv/shiny-server/referenceranger/TMC.R")
  source("/srv/shiny-server/referenceranger/TML.R")
}

resultvalidator <- "function(value, callback) {
  setTimeout(function() {
   if (value === null || value === '') {
    callback(true);
   } else {
    var regex = /^<?\\d+(?:[\\.,]\\d+)?$/;
    callback(regex.test(value.trim()));
   }
  }, 100);}"


# Function for plotting the RI limits for comparison
plotcomparisonlimits <- function (estimatedlimits.low, estimatedlimits.high, referencelimits.low, referencelimits.high) {
  pU <- permissible_uncertainty(referencelimits.low, referencelimits.high)
  if (estimatedlimits.low > pU[1] && estimatedlimits.low < pU[2]) col1 = '#10D010'
  else col1 = 'red'
  if (estimatedlimits.high > pU[3] && estimatedlimits.high < pU[4]) col2 = '#10D010'
  else col2 = 'red'
  rect(pU[1], 0, pU[2], 99999, col = adjustcolor(col1, alpha.f = 0.05), border = FALSE)
  rect(pU[3], 0, pU[4], 99999, col = adjustcolor(col2, alpha.f = 0.05), border = FALSE)
  abline(v = referencelimits.low, col = adjustcolor(col1, alpha.f = 0.40), lwd = 2, lty = 3)
  abline(v = referencelimits.high, col = adjustcolor(col2, alpha.f = 0.40), lwd = 2, lty = 3)
  abline(v = estimatedlimits.low, col = "#10D010", lwd = 3)
  abline(v = estimatedlimits.high, col = "#10D010", lwd = 3)
}

# Dashboard and sidebar 
ui <-
  dashboardPage(
    dashboardHeader(title = tags$span("ReferenceRangeR", class = "header-title")),
    
    dashboardSidebar(
      tags$head(
        tags$style(HTML("
      .header-title {
          font-size: 16px;
          font-style: italic;
          font-weight: bold;
          text-align: left;
        }
      .sidebar {
          display: flex;
          flex-direction: column;
          height: 100vh;
          overflow: hidden;
        }
      .sidebar-scroll {
          flex: 1;
          overflow-y: auto;
          padding-bottom: 100px;
        }
      .sidebar-image {
          padding-bottom: 50px;
          text-align: center;
        }
      .sidebar-image img {
          padding: 0px;
          display: block;
          margin-left: auto;
          margin-right: auto;
          height: 50px;
        }
"))
      ),
      use_theme(uoltheme),
      tags$div(class = "sidebar-scroll",
               sidebarMenu(
                 br(),
                 div(img(src = 'rrr.webp', width = '110px'), style = 'text-align: center;'),
                 menuItem(tags$div(
                   style = "text-align:center; font-size: 0.8em; font-style: italic;",
                   tags$a(href = "https://kc.uol.de/referenceranger/help.pdf", "help", target = "_blank"),
                   "  |  ",
                   tags$a(href = "https://github.com/gubruol/ReferenceRangeR/", "github", target = "_blank"),
                   "  |  ",
                   tags$a(href = "https://kc.uol.de/disclaimer/", "disclaimer", target = "_blank")
                 )),
                 br(),
                 menuItem("select sex", startExpanded = FALSE,
                          br(),
                          actionButton("sexbox", HTML("visualize data"), style = buttoncolors1, width = '100%'),
                          radioButtons("sexradio", "", c("all" = "A", "male" = "M", "female" = "F", "non-binary" = "D")),
                          conditionalPanel(condition = "output.pregnancymode == 1", radioButtons("trimester", "", c("1. trimester" = 1, "2. trimester" = 2, "3. trimester" = 3))),
                          br(),
                          conditionalPanel(condition = "output.pregnancymode != 1", actionButton("pregnancy_button", HTML("add trimester column"), style = buttoncolors2, width = '90%')),
                          br()
                 ),
                 menuItem("select age", startExpanded = FALSE,
                          br(),
                          HTML("<div style='text-align:center;width:100%;font-size:80%'><i>from / to</i></div>"),
                          fluidRow(column(6, numericInput("agell", NULL, "0")), column(6, numericInput("ageul", NULL, "0"))),
                          actionButton("drift", HTML("drift check"), style = buttoncolors1, width = '100%'),
                          conditionalPanel(condition = "output.showstratslider != '0'", sliderInput("stratnum", "preferred no. of groups", min = 2, max = 12, value = 3, step = 1,ticks = FALSE)),
                          br()
                 ),
                 menuItem("select method", startExpanded = FALSE,
                          radioButtons("methodradio", "", choiceNames = list(
                            HTML("<b>refineR</b> <a href = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8346497/pdf/41598_2021_Article_95301.pdf'>(Lit.)</a>"),
                            HTML("<b>TMC</b> <a href = 'https://www.degruyter.com/document/doi/10.1515/cclm-2018-1341/html'>(Lit.)</a>"),
                            HTML("<b>TML</b> <a href = 'https://www.degruyter.com/document/doi/10.1515/CCLM.2007.249/html'>(Lit.)</a>"),
                            HTML("<b>kosmic</b> <a href = 'https://www.nature.com/articles/s41598-020-58749-2'>(Lit.)</a>"),
                            HTML("<b>reflimR</b> <a href = 'https://doi.org/10.1515/labmed-2023-0042'>(Lit.)</a>")
                          ),
                          choiceValues = list("refiner", "tmc", "tml", "kosmic", "reflimr")
                          ),
                          conditionalPanel(condition = "input.methodradio == 'tml'", checkboxInput(
                            "fasttml", HTML("Fast mode <i style = 'font-size:80%;'>(3 significant figures)</i>"),
                            value = TRUE
                          )),
                          br()
                 ),
                 menuItem("select limits for comparison", startExpanded = FALSE,
                          HTML("<br><div style='text-align:center;width:100%;font-size:80%'><i>lower limit / upper limit </i>"),
                          fluidRow(
                            column(6, numericInput("referencelimits.low", NULL, "0.0", step = 0.1)),
                            column(6, numericInput("referencelimits.high", NULL, "0.0", step = 0.1))),
                          HTML("<div style='text-align:center;width:100%;font-size:80%'>Permissible uncertainty <a href = 'https://www.degruyter.com/document/doi/10.1515/cclm-2014-0874/html'>(Lit.1 </a> and <a href = 'https://doi.org/10.1515/labmed-2023-0042'>Lit.2)</a></div>"),
                          br()
                 ),
                 br(),
                 actionButton("calc", HTML("<b>calculate</b>"), style = buttoncolors1, width = '100%'),
                 br(),
                 div(style = "display: flex; justify-content: space-between;",
                     actionButton("clear", "clear data", style = buttoncolors2, width = '100%'),
                     actionButton("demo", "demo data", style = buttoncolors2, style = 'margin-right: 10px;', width = '100%'))
                 
               )
      ),
      tags$div(class = "sidebar-image",
               tags$img(src = "umo.svg"))
    ),
    
    dashboardBody(
      add_busy_spinner(spin = "folding-cube", color = "#003B73", position = "full-page", height = "150px", width = "150px", onstart = FALSE, timeout = 2000),
      box(title = "plot", id = "boxplot", width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", plotOutput("plot"), htmlOutput("placeholder", style = "font-size:80%;")),
      box(title = "table", id = "boxtable", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary", rHandsontableOutput("table"))
    )
  )

# Server initialization
server <- function(input, output, session) {
  setwd(tempdir())
  dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist))
  dataframe$result = as.character(dataframe$result)
  dataframe$age = as.numeric(dataframe$age)
  output$table <- renderRHandsontable(rhandsontable(dataframe, width = '100%', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                        hot_col("result", validator = resultvalidator) %>%
                                        hot_col("sex", allowInvalid = TRUE)
  )
  
  output$showstratslider <- renderText({ 
    '0'
  })
  outputOptions(output, "showstratslider", suspendWhenHidden=FALSE)
  
  output$pregnancymode <- renderText({
    '0'
  })
  outputOptions(output, "pregnancymode", suspendWhenHidden=FALSE)
  
  observeEvent(input$table, {
    dataframe = hot_to_r(input$table)
  })
  
  # Clear data and table
  observeEvent(input$clear, {
    dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist))
    dataframe$age = as.numeric(dataframe$age)
    dataframe$result = as.character(dataframe$result)
    updateNumericInput(session, "referencelimits.low", value = 0)
    updateNumericInput(session, "referencelimits.high", value = 0)
    updateNumericInput(session, "agell", value = "")
    updateNumericInput(session, "ageul", value = "")
    output$showstratslider <- renderText({ '0' })
    output$pregnancymode <- renderText({ '0' })
    if (input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '100%', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                        hot_col("result", validator = resultvalidator) %>%
                                        hot_col("sex", allowInvalid = TRUE)
  )
  })
  
  # Generate demo data
  observeEvent(input$demo, {
    dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist))
    dataframe$result = as.numeric(dataframe$result)
    dataframe$age = as.numeric(dataframe$age)
    demosamplesize <- tablesize * 0.5
    dataframe$age[1:demosamplesize] <- round(runif(demosamplesize, min = 1, max = 99), digits = 2)
    dataframe$sex[1:demosamplesize] <- factor(sample(c("M", "F"), demosamplesize, replace = TRUE), levels = sexlist)
    dataframe$result[1:demosamplesize] <- SSlogis(dataframe$age[1:demosamplesize], 5, 60, 6) + rnorm(demosamplesize, mean = 37.5, sd = 3.75)
    pathsize <- round(demosamplesize / 100, digits = 0) # 10% pathological values
    dataframe$result[1:pathsize] <- SSlogis(dataframe$age[1:pathsize], 5, 60, 6) + rnorm(pathsize, 20, 5)
    dataframe$result[(1 + pathsize):(3 * pathsize)] <- SSlogis(dataframe$age[(1 + pathsize):(3 * pathsize)], 5, 60, 6) + rnorm(2 * pathsize, 28, 3)
    dataframe$result[(1 + 3 * pathsize):(8 * pathsize)] <- SSlogis(dataframe$age[(1 + 3 * pathsize):(8 * pathsize)], 5, 60, 6) + rnorm(5 * pathsize, 50, 4.5)
    dataframe$result[(1 + 8 * pathsize):(10 * pathsize)] <- SSlogis(dataframe$age[(1 + 8 * pathsize):(10 * pathsize)], 5, 60, 6) + rnorm(2 * pathsize, 60, 10)
    dataframe$result <- round(dataframe$result, 2)
    dataframe$result = as.character(dataframe$result)
    updateNumericInput(session, "referencelimits.low", value = 30)
    updateNumericInput(session, "referencelimits.high", value = 45)
    updateNumericInput(session, "agell", value = 18)
    updateNumericInput(session, "ageul", value = 45)
    output$showstratslider <- renderText({ '0' })
    output$pregnancymode <- renderText({ '0' })
    if (input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '100%', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                        hot_col("result", validator = resultvalidator) %>%
                                        hot_col("sex", allowInvalid = TRUE)
  )
  })
  
  # Add trimester column
  observeEvent(input$pregnancy_button, {
    dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist), trimester = factor(rep(NA, tablesize), levels = c(1, 2, 3)))
    dataframe$result = as.character(dataframe$result)
    dataframe$age = as.numeric(dataframe$age)
    updateNumericInput(session, "referencelimits.low", value = 0)
    updateNumericInput(session, "referencelimits.high", value = 0)
    updateNumericInput(session, "agell", value = "")
    updateNumericInput(session, "ageul", value = "")
    output$showstratslider <- renderText({ '0' })
    output$pregnancymode <- renderText({ '1' })
    if (input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '100%', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                        hot_col("result", validator = resultvalidator) %>%
                                        hot_col("sex", allowInvalid = TRUE)
  )
  })
  
  # Visualize data for sex differences
  observeEvent(input$sexbox, {
    dataframe = hot_to_r(input$table)
    dataframe$result <- as.numeric(dataframe$result)
    dataframe <- dataframe[dataframe$result > 0, ]
    dataframe <- dataframe[!is.na(dataframe$result), ]
    agell <- isolate(input$agell)
    ageul <- isolate(input$ageul)
    agelimitsvalid <- (ageul > 0 && (ageul > agell) && !is.na(agell) && !is.na(ageul))  
    if (agelimitsvalid) dataframe <- dataframe[(dataframe$age >= agell) & (dataframe$age <= ageul), ]
    dataframe <- dataframe %>% mutate(sex = ifelse(sex %in% femalelist, 'F', ifelse(sex %in% malelist, 'M', ifelse(sex == 'D', 'D', 'X'))))
    dataframe <- dataframe[!is.na(dataframe$sex), ]
    cases <- sum(dataframe$result > 0 & !is.na(dataframe$result) & dataframe$sex != "" & !is.na(dataframe$sex))
    dataframe <- dataframe %>% group_by(sex) %>% filter(n() >= 100)
    casesfiltered <- sum(dataframe$result > 0 & !is.na(dataframe$result) & dataframe$sex != "" & !is.na(dataframe$sex))
    if (casesfiltered < cases) dataremoved = TRUE else dataremoved = FALSE
    if (casesfiltered < 500) shinyalert("insufficient data...", paste("The selected dataset contains ", casesfiltered, " results.\nPlease increase the sample size (minimum n=500)."), type = "error")
    else {
      q10 <- quantile(dataframe$result, probs = 0.1)
      q90 <- quantile(dataframe$result, probs = 0.9)
      ylim.min <- q10 - (q90 - q10) / 1.3
      ylim.max <- q90 + (q90 - q10) / 1.3
      pu_percent <- 2.39*(-0.25 + 100*(-1+exp(((log(reflim(dataframe$result)$limits[2])-log(reflim(dataframe$result)$limits[1]))/3.92)^2))^0.5)^0.5
      n_groups <- length(unique(dataframe$sex))
      pvalue <- -1
      test_used <- ""
      if (n_groups == 2) {
        pvalue <- round(wilcox.test(result ~ sex, data = dataframe)$p.value, 2)
        test_used <- "Wilcoxon rank-sum test"
      }
      if (n_groups > 2) {
        pvalue <- round(kruskal.test(result ~ sex, data = dataframe)$p.value, 2)
        test_used <- "Kruskal-Wallis test"
      }
      mediantable <- dataframe %>% group_by(sex) %>% summarise(count = n(), median = median(result, na.rm = TRUE))
      mediandiff <- round((max(mediantable$median) - min(mediantable$median)) / min(mediantable$median) * 100, digits = 2)
      output$plot <- renderPlot({
        ggplot(dataframe, aes(sex, result)) + ylim(ylim.min, ylim.max) + labs(x = NULL) +
          theme_minimal() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
          geom_violin(fill = "#40668d", alpha = 0.2) +
          geom_boxplot(fill = "#40668d", outlier.shape = NA, width = 0.1)
      })
      output$placeholder = renderText({
        paste("<br>",
              kable(mediantable, "html", align = "c") %>%  kable_styling(),
              "<br><br><p style='font-size: 16px;'><b>Selected data:</b><br>",
              "n = ", casesfiltered, "<br>",
              if (dataremoved) "Groups with n<100 were removed<br>",
              if (agelimitsvalid) paste("age: from ", agell, " to ", ageul)
              else "age: no selection",
              "<br><br>",
              if (pvalue < 0) "More than one group needed for statistical evaluation."
              else if (pvalue < 0.01) paste("<b style='color: red;'>Significant difference </b> (p<0.01, ", test_used, ")<br>", sep = "")
              else if (pvalue < 0.05) paste("<b style='color: red;'>Significant difference</b> (p=", pvalue, ", ", test_used, ")<br>", sep = "")
              else paste("<b style='color: #10D010;'>No significant difference</b> (p=", pvalue, ", ", test_used, ")<br>", sep = ""),
              "<br><b>Max. deviation of the medians: </b> ",
              if (mediandiff > pu_percent) "<b style='color: red;'>"
              else if (mediandiff > (0.5 * pu_percent)) "<b style='color: orange;'>"
              else "<b style='color: #10D010;'>",
              mediandiff, "%</b><br>",
              "Estimated permissible uncertainty (pU): ", round(pu_percent,2), "%</p>"
        )
      })
      if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
      if (input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    }
  })
  
  # Visualize drift
  observeEvent(input$drift, {
    dataframe = hot_to_r(input$table)
    dataframe$result <- as.numeric(dataframe$result)
    dataframe <- dataframe[dataframe$result > 0, ]
    dataframe <- dataframe[!is.na(dataframe$result), ]
    sexradio <- isolate(input$sexradio)
    if (sexradio == 'M') dataframe <- dataframe[toupper(substr(dataframe$sex, 1, 1)) == 'M', ]
    else if (sexradio == 'F') dataframe <- dataframe[toupper(substr(dataframe$sex, 1, 1)) == 'F' | toupper(substr(dataframe$sex, 1, 1)) == 'W', ]
    else if (sexradio == 'D') dataframe <- dataframe[dataframe$sex == 'D', ]
    if ("trimester" %in% colnames(dataframe)) dataframe <- dataframe[dataframe$trimester == input$trimester, ]
    agell <- isolate(input$agell)
    ageul <- isolate(input$ageul)
    agelimitsvalid <- (ageul > 0 && (ageul > agell) && !is.na(agell) && !is.na(ageul))  
    if (agelimitsvalid) dataframe <- dataframe[(dataframe$age >= agell) & (dataframe$age <= ageul), ]
    cases <- sum(dataframe$result > 0 & !is.na(dataframe$result) & dataframe$age > 0 & !is.na(dataframe$age))
    
    # Subsampling and reduction of sample size for n>10000
    dataframe[sample(nrow(dataframe)), ]
    if (cases > 10000) dataframe <- dataframe[sample(1:length(dataframe$result), 10000), ]
    
    # Alert when n<500
    if (cases < 500) shinyalert("insufficient data...", paste("The selected dataset contains ", cases, " results with given age.\nPlease increase the sample size (minimum n=500)."), type = "error")
    else {
      q10 <- quantile(dataframe$result, probs = 0.1)
      q90 <- quantile(dataframe$result, probs = 0.9)
      ylim.min <- q10 - (q90 - q10) / 1.3
      ylim.max <- q90 + (q90 - q10) / 1.3
      
      pu_percent <- 2.39*(-0.25 + 100*(-1+exp(((log(reflim(dataframe$result)$limits[2])-log(reflim(dataframe$result)$limits[1]))/3.92)^2))^0.5)^0.5
      pu_absolute <- pu_percent * median(dataframe$result) / 100
      
      # Calculate number of digits for age with resulting 10^3 steps or groups.
      agedigits <- 3 - floor(log10(max(dataframe$age)))
      
      agegroups <- data.frame(from = seq(round(min(dataframe$age), agedigits), round(max(dataframe$age), agedigits), by=10^-agedigits)) %>%
        mutate(to = from + 10^-agedigits)
      
      # Predict medians and confidence intervals
      pred <- predict(qgam(result ~ s(age), data = dataframe, qu = 0.5), newdata = data.frame(age = (agegroups$from + agegroups$to) / 2), se=TRUE)
      
      agegroups <- agegroups %>% mutate(
        median = pred$fit,
        ci.low = pred$fit - 1.96 * pred$se.fit,
        ci.high = pred$fit + 1.96 * pred$se.fit) %>%
        rowwise() %>%
        mutate(count = sum(dataframe$age >= from & dataframe$age < to)) %>%
        ungroup()
      
      # Remove empty groups
      agegroups <- agegroups[agegroups$count > 0, ]
      
      stratplot <- ggplot() + ylim(ylim.min, ylim.max) +
        geom_point(data = dataframe, aes(x=age, y=result), color="grey", alpha=0.1) +
        geom_line(data = agegroups, aes(x = from, y = median), color = "#003B73", linewidth = 1) +
        geom_ribbon(data = agegroups, aes(x = from, ymin = ci.low, ymax = ci.high), fill = "#003B73", alpha = 0.2) +
        theme_minimal()
      
      agegroups$ci.low <- NULL
      agegroups$ci.high <- NULL
      
      while (TRUE) {
        if (length(agegroups$median) <= 1) break
        diffs <- diff(agegroups$median)
        idx <- which.min(diffs)
        if ((length(agegroups$median) > input$stratnum) && (diffs[idx] >= 0.5 * pu_absolute)) break
        if ((length(agegroups$median) <= input$stratnum) && (diffs[idx] >= 0.25 * pu_absolute)) break
        agegroups$to[idx] <- agegroups$to[idx + 1]
        agegroups$count[idx] <- agegroups$count[idx] + agegroups$count[idx + 1]
        agegroups$median[idx] <- median(dataframe[dataframe$age >= agegroups$from[idx] & dataframe$age <= agegroups$to[idx], ]$result)
        agegroups <- agegroups[-(idx + 1), ]
      }
      
      if (length(agegroups$median) > 1) output$showstratslider <- renderText({ '1' })
      else output$showstratslider <- renderText({ '0' })
      
      # Remove groups with n < 50
      agegroups <- agegroups[agegroups$count >= 50, ]
      
      agegroups$from <- round(agegroups$from, agedigits - 1)
      agegroups$to <- round(agegroups$to, agedigits - 1)
      if (max(agegroups$median) > 100) agegroups$median <- round(agegroups$median, 0)
      else if (max(agegroups$median) > 10) agegroups$median <- round(agegroups$median, 1)
      else agegroups$median <- round(agegroups$median, 2)
      agegroups$size <- round(100 * agegroups$count / sum(agegroups$count), 0)
      agegroups$count <- NULL
      agegroups$sizecol <- col_numeric(palette = c("red","red","#10D010","#10D010","#10D010","#10D010","#10D010","#10D010","#10D010","#10D010"), domain=c(0, 100))(agegroups$size)
      agegroups$size <- paste(agegroups$size, "%")
      mediandeviation <- round((100 * (max(agegroups$median) - min(agegroups$median)) / min(agegroups$median)), 2)
      
      if (length(agegroups$size) < 2) kabletable <- "No stratification necessary."
      else {
        kabletable <- agegroups[, -5] %>% kable("html",  align = "c") %>% kable_styling() %>%
          column_spec(4, bold = TRUE, color = agegroups$sizecol)
        
        stratplot <- stratplot +
          geom_rect(data = agegroups, aes(NULL, NULL, xmin = from, xmax = to, ymin = -Inf, ymax = Inf, fill = sizecol), alpha = 0.1) +
          scale_fill_identity() +
          theme(legend.position="none")
      }
      
      output$plot <- renderPlot({ stratplot })
      
      output$placeholder = renderText({
        paste(
          "<p style='font-size: 16px;'><b>Selected data:</b><br>",
          "n = ", cases,
          if (agelimitsvalid) paste("<br>age: from ", agell, " to ", ageul)
          else paste("<br>age: no selection"),
          "<br>sex: ",
          if (sexradio == 'M') "male"
          else if (sexradio == 'F') "female"
          else if (sexradio == 'D') "non-binary"
          else "no selection",
          "<br><br><b>Max. deviation of the medians: </b>",
          if (mediandeviation > pu_percent) "<b style='color: red;'>"
          else if (mediandeviation > (0.5 * pu_percent)) "<b style='color: orange;'>"
          else "<b style='color: #10D010;'>",
          mediandeviation, "%</b><br>",
          "Estimated permissible uncertainty (pU): ", round(pu_percent,2), "%<br><br>",
          "<b>Recommended age limits for stratification:</b>",
          kabletable,
          "</p>"
        )
      })
      if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
      if (input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    }
  })
  
  # Start calculation
  observeEvent(input$calc, {
    dataframe = hot_to_r(input$table)
    methodradio <- isolate(input$methodradio)
    dataframe$result <- gsub(",", ".", dataframe$result, fixed = TRUE)
    dataframe$age <- gsub(",", ".", dataframe$age, fixed = TRUE)
    if (methodradio == 'tmc') dataframe$result <- gsub("^<", "", dataframe$result)
    dataframe$result <- as.numeric(dataframe$result)
    dataframe <- dataframe[dataframe$result > 0, ]
    dataframe <- dataframe[!is.na(dataframe$result), ]
    sexradio <- isolate(input$sexradio)
    if (sexradio == 'M') dataframe <- dataframe[toupper(substr(dataframe$sex, 1, 1)) == 'M', ]
    else if (sexradio == 'F') dataframe <- dataframe[toupper(substr(dataframe$sex, 1, 1)) == 'F' | toupper(substr(dataframe$sex, 1, 1)) == 'W', ]
    else if (sexradio == 'D') dataframe <- dataframe[dataframe$sex == 'D', ]
    if ("trimester" %in% colnames(dataframe)) dataframe <- dataframe[dataframe$trimester == input$trimester, ]
    agell <- isolate(input$agell)
    ageul <- isolate(input$ageul)
    agelimitsvalid <- (ageul > 0 && (ageul > agell) && !is.na(agell) && !is.na(ageul))  
    if (agelimitsvalid) dataframe <- dataframe[(dataframe$age >= agell) & (dataframe$age <= ageul), ]
    referencelimits.low <- isolate(input$referencelimits.low)
    referencelimits.high <- isolate(input$referencelimits.high)
    referencelimitsvalid <- (referencelimits.high > 0 && (referencelimits.high > referencelimits.low) && !is.na(referencelimits.low) && !is.na(referencelimits.high))  
    estimatedlimits.low <- 0
    estimatedlimits.high <- 0
    if (length(na.omit(dataframe$result)) < 500) shinyalert("insufficient data...", paste("The selected dataset contains ", length(na.omit(dataframe$result)), " results.\nPlease increase the sample size (minimum n=500)."), type = "error")
    else {
      if (length(na.omit(dataframe$result)) < 2000) shinyalert("small sample size...", paste("The selected dataset contains only ", length(dataframe$result), " results.\nYou have been warned..."), type = "warning")
      if (length(na.omit(dataframe$result)) > 100000) shinyalert("large sample size...", paste("The selected dataset contains ", length(dataframe$result), " results.\nYou have been warned..."), type = "warning")
      
      if (methodradio == 'refiner') {
        resri <- findRI(na.omit(dataframe$result))
        estimatedlimits.low <- getRI(resri)[1, 2]
        estimatedlimits.high <- getRI(resri)[2, 2]
        output$plot <- renderPlot({
          plot(resri)
          if (referencelimitsvalid) plotcomparisonlimits(getRI(resri)[1, 2], getRI(resri)[2, 2], referencelimits.low, referencelimits.high)
        })
      }
      else if (methodradio == 'tmc') {
        output$plot <- renderPlot({
          temptmc <- tmc(na.omit(dataframe$result))
          estimatedlimits.low <<- temptmc$RL1
          estimatedlimits.high <<- temptmc$RL2
          replayPlot(temptmc$myplot)
          if (referencelimitsvalid)
            plotcomparisonlimits(temptmc$RL1, temptmc$RL2, referencelimits.low, referencelimits.high)
        })
      }
      else if (methodradio == 'tml') {
        if (quantile(dataframe$result, probs = 0.90) >= 1000) decimalcount <- 0
        else if (quantile(dataframe$result, probs = 0.90) >= 100) decimalcount <- 1
        else if (quantile(dataframe$result, probs = 0.90) >= 10) decimalcount <- 2
        else decimalcount <- 3
        if (input$fasttml == TRUE && decimalcount > 0) decimalcount <- decimalcount - 1
        dataframe$result <- round(dataframe$result, decimalcount)
        estimateleft <- length(dataframe$result[dataframe$result < reflim(dataframe$result)$limits[1]])
        estimateright <- length(dataframe$result[dataframe$result > reflim(dataframe$result)$limits[2]])
        if (estimateright > estimateleft) pathright <- TRUE
        else pathright <- FALSE
        output$plot <- renderPlot({
          temptml <- tml(na.omit(dataframe$result), pathright)
          estimatedlimits.low <<- temptml$DL25
          estimatedlimits.high <<- temptml$DL975
          replayPlot(temptml$myplot)
          if (referencelimitsvalid)
            plotcomparisonlimits(temptml$DL25, temptml$DL975, referencelimits.low, referencelimits.high)
        })
      }
      else if (methodradio == 'kosmic') {
        resri <- kosmic(na.omit(dataframe$result), decimals = 1)
        estimatedlimits.low <- summary(resri)[1]
        estimatedlimits.high <- summary(resri)[3]
        output$plot <- renderPlot({
          plot(resri)
        })
      }
      else if (methodradio == 'reflimr') {
        estimatedlimits.low <- reflim(dataframe$result)$limits[1]
        estimatedlimits.high <- reflim(dataframe$result)$limits[2]
        if (referencelimitsvalid)
          output$plot <- renderPlot(reflim(dataframe$result, targets = c(referencelimits.low, referencelimits.high)))
        else output$plot <- renderPlot(reflim(dataframe$result, targets = NULL))
      }
      
      if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
      if (input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
      output$placeholder = renderText({
        if (estimatedlimits.low > 100) decimalcount <- 0
        else if (estimatedlimits.low > 10) decimalcount <- 1
        else decimalcount <- 2
        if (!is.na(estimatedlimits.low)) estimatedlimits.lowtxt <- sprintf(paste("%.", decimalcount, "f", sep = ""), as.numeric(estimatedlimits.low))
        if (!is.na(estimatedlimits.high)) estimatedlimits.hightxt <- sprintf(paste("%.", decimalcount, "f", sep = ""), as.numeric(estimatedlimits.high))
        if (referencelimitsvalid) {
          referencelimits.low.txt <- sprintf(paste("%.", decimalcount, "f", sep = ""), referencelimits.low)
          referencelimits.high.txt <- sprintf(paste("%.", decimalcount, "f", sep = ""), referencelimits.high)
        }
        paste(
          "<p style='font-size: 16px;'><b>Selected data:</b><br>",
          "n = ", length(na.omit(dataframe$result)),
          if (agelimitsvalid) paste("<br>age: from ", agell, " to ", ageul)
          else paste("<br>age: no selection"),
          "<br>sex: ",
          if (sexradio == 'M') "male"
          else if (sexradio == 'F') "female"
          else if (sexradio == 'D') "non-binary"
          else "no selection",
          "<br>Method: ",
          if (methodradio == 'refiner') "refineR"
          else if (methodradio == 'TMC') "TMC"
          else if (methodradio == 'tml') "TML"
          else if (methodradio == 'kosmic') "kosmic"
          else if (methodradio == 'reflimr') "reflimR",
          "<br><br><b>Estimated reference limits: ", estimatedlimits.lowtxt, " - ", estimatedlimits.hightxt, "</b><br>",
          if (referencelimitsvalid) paste("<b style='color: #D0D0D0;'>Comparison reference limits: ", referencelimits.low.txt, " - ", referencelimits.high.txt, "</b><br>"),
          "</p>"
        )
      })
    }
  })
}

shinyApp(ui, server)
