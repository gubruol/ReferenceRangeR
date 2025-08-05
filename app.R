library(htmltools)
library(bslib)
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
buttoncolors2 <- "color: #000; background-color: #BBBBBB; border-color: #000;width: 90%;"
buttoncolors3 <- "color: #000; background-color: #4c779a; border-color: #000;width: 90%;"
buttoncolors4 <- "color: #000; background-color: #BBBBBB; border-color: #000;width: 90%; margin-left: 10px; margin-right: 10px;"


# Variables
malelist <- c('male', 'männlich', 'Mann', 'M', 'm')
femalelist <- c('female', 'weiblich', 'Frau', 'F', 'f', 'W', 'w')
diverselist = c("D", "d", "diverse","Diverse")
sexlist <- c(malelist, femalelist, diverselist)
tablesize <- 200000
input_sex=NULL
dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist))
fasttml=NULL


# Path for TMC and TML library
if (Sys.info()["sysname"] == "Windows") {
  source("Z:/R/kc.uol.de/referenceranger/TMC.R")
  source("Z:/R/kc.uol.de/referenceranger/TML.R")
} else if (Sys.info()["sysname"] == "Linux") {
  source("/srv/shiny-server/referenceranger/TMC.R")
  source("/srv/shiny-server/referenceranger/TML.R")
} else {
  source("./referenceranger/TMC.R")
  source("./referenceranger/TML.R")
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
    dashboardHeader(
      title = tags$span("ReferenceRangeR", class = "header-title")),
    
    dashboardSidebar(
      tags$head(
        tags$style(HTML("
        
      .header-title {
          font-size: 16px;
          font-style: italic;
          font-weight: bold;
          text-align: left;
      }
      .header-button {
          font-size: 14px;
          font-style: regular;
          text-align: center;
          padding-bottom: 50px;
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
      .radioButtons {
          line-hight: 10px;
          font-size: 12px;
        }

      .layout-sidebar1 {
          overflow-y: hidden !important;
          max-height: 700px !important;
        }
        
        .layout-sidebar2 {
           overflow-y: hidden !important;
           max-height: 750px !important;
        }
        .red-border input.form-control {
          border: 2px solid red !important;
          box-shadow: none !important;
        }
                        "))
        ),
      use_theme(uoltheme),
      tags$div(
        class = "sidebar-scroll",
        sidebarMenu(
          br(),
          div(img(
            src = 'rrr.webp', width = '150px'), style = 'text-align: center;'),
            menuItem(tags$div(
              style = "text-align:center; font-size: 0.8em; font-style: italic;",
              tags$a(
                href = "https://kc.uol.de/referenceranger/help.pdf", "help", target = "_blank"),
                "  |  ",
                tags$a(href = "https://github.com/gubruol/ReferenceRangeR/", "github", target = "_blank"),
                "  |  ",
                tags$a(href = "https://kc.uol.de/disclaimer/", "disclaimer", target = "_blank")
                 )
              ),
          br(),
          br(),
          actionButton("sexbox", HTML("Check sex differences"), style = buttoncolors1, width = '100%'),
          actionButton("drift", HTML("Check age drift"), style = buttoncolors1, width = '100%'),
          br(),
          br(),
          actionButton("calc", HTML("<b>Calculate</b>"), style = buttoncolors1, width = '100%'),
          br(),
          div(style = "display: flex; justify-content: space-between;",
              actionButton("clear", "Clear data", style = buttoncolors2, width = '100%'),
              actionButton("demo", "Demo data", style = buttoncolors2, style = 'margin-right: 10px;', width = '100%')
            )
          )
        ),
        conditionalPanel(
              condition = "output.advanced == 0",
              tags$div(class="header-button",
                       actionButton(inputId="advanced_button",HTML("Advanced mode"), style = buttoncolors3, width = '100%')
                  )),
            conditionalPanel(
              condition = "output.advanced == 1",
              tags$div(class="header-button",
              actionButton(inputId="standard_button",HTML("Standard mode"), style = buttoncolors3, width = '100%')
              )
              ),  
        tags$div(class = "sidebar-image",
               tags$img(src = "umo.svg"))
    ),
    
    dashboardBody(
      useShinyjs(),
      add_busy_spinner(spin = "folding-cube", color = "#003B73", position = "full-page", height = "150px", width = "150px", onstart = FALSE, timeout = 2000),
          as.card_item(
            box(
              title=actionLink("Data Input",inputId="DI"), 
              id = "boxtable",
              width = 8, 
              solidHeader = TRUE, 
              collapsible = TRUE, 
              status = "primary",
              sidebar= boxSidebar(
                title= "Settings",
                id= "initbox",
                width = 30,
                collapsible = TRUE,
                collapsed = TRUE,
                HTML("<div style='text-align:left;width:100%;font-size:100%'><i>Assign sex variables:</i></div>"),
                selectInput(
                  "Init_female",
                  "female:",
                  choices="F"),
                selectInput(
                  "Init_male",
                  "male:",
                  choices="M"),
                selectInput(
                  "Init_diverse",
                  "diverse:",
                  choices="D"),
                dropdownDivider(),
                conditionalPanel(
                  condition = "output.advanced==1",
                  conditionalPanel(
                    condition = "output.pregnancymode != 1", 
                    actionButton("pregnancy_button", HTML("add trimester column"), 
                                 style = buttoncolors4)
                ),
                conditionalPanel(
                  condition = "output.pregnancymode == 1", 
                  actionButton("remove_pregnancy", HTML("remove trimester column"), 
                               style = buttoncolors4)
                ))),
              card(
                max_height = "700px",
                rHandsontableOutput("table")
              )
            )
          ), 
        as.card_item(
          box(
            title= actionLink("Visualization and Stratification", inputId="SV"),
            id = "strat_boxplot", 
            width = 12, 
            solidHeader = TRUE, 
            collapsible = TRUE, 
            collapsed = TRUE, 
            status = "primary",
              card(
                max_height= 700,
                tags$div(class = "layout-sidebar1"
                         ,layout_sidebar( 
                sidebar=sidebar(
                  id = "sidebar_strat",
                  bg = "#003B73",
                  position = "left", 
                  open = "desktop",
                  HTML("<div style='text-align:center;width:100%;font-size:100%'><r>Visualization settings:</r></div>"),
                  tags$div(class = "radioButtons",
                  radioButtons("sexradio", "", c("all" = "A", "male" = "M", "female" = "F", "non-binary" = "D"))),
                  HTML("<div style='text-align:center;width:100%;font-size:80%'><i>Age limits from / to</i></div>"),
                  fluidRow(
                    column(
                      6, 
                      div(
                      id = "strat_agell_container",
                      numericInput("strat_agell", NULL, "0"))), 
                    column(
                      6,
                      div(
                        id = "strat_ageul_container",
                      numericInput("strat_ageul", NULL, "0")))),
                  conditionalPanel(
                    condition = "output.advanced == 1", 
                    sliderInput("stratnum", "max. no. of age groups", min = 2, max = 12, value = 3, step = 1,ticks = FALSE)
                    )
                  ),
                    card(
                      max_height= "320",
                      plotOutput("strat_plot"
                              )),
                card(
                  max_height= "330",
                    htmlOutput("strat_placeholder"
                               )))
                  ))
          )),
      as.card_item(
        box(
          title = actionLink("Reference Interval", inputId="RI"), 
          id = "boxplot", 
          width = 12, 
          solidHeader = TRUE, 
          collapsible = TRUE, 
          collapsed = TRUE, 
          status = "primary",
          card(
            max_height = "700px",
            tags$div(class = "layout-sidebar2"
                     ,layout_sidebar( 
            sidebar=sidebar( 
              id = "sidebar_ref",
              bg = "#003B73",
              fillable = T,
              position = "left", 
              open = "desktop",
              HTML("<div style='text-align:left;width:100%;line-height: 10px;font-size:90%'><r>Calculation settings:</r></div>"),
              tags$div(class = "radioButtons",
                       radioButtons("sexradioCalc", "",
                           choiceNames = list(
                             HTML("<r>all</r> "),
                             HTML("<r>male</r> "),
                             HTML("<r>female</r> "),
                             HTML("<r>non-binary</r> ")
                           ),
                           choiceValues = list("A","M","F","D")
                           )
                       ),
              conditionalPanel(
                condition = "output.pregnancymode == 1", 
                tags$div(class = "radioButtons",
                         radioButtons("trimesterCalc", "", 
                                      c("1. trimester" = 1, "2. trimester" = 2, "3. trimester" = 3)))
                ),
              HTML("<div style='text-align:center;width:100%;font-size:80%'><r>Select age: <i>from / to</i></div>"),
              fluidRow(
                column(
                  6,
                  div(
                    id = "agell_container",
                    numericInput("agell", NULL, "0"))), 
                column(
                  6,
                  div(
                    id = "ageul_container",
                    numericInput("ageul", NULL, "0")))),
              HTML("<div style='text-align:left;width:100%;line-height: 10px;font-size:90%'><r>Select method:</r></div>"),
              tags$div(class = "radioButtons",
                       radioButtons("methodradio", "", choiceNames = list(
                         HTML("<b>refineR</b> <a href = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8346497/pdf/41598_2021_Article_95301.pdf'>(Lit.)</a>"),
                         HTML("<b>TMC</b> <a href = 'https://www.degruyter.com/document/doi/10.1515/cclm-2018-1341/html'>(Lit.)</a>"),
                         HTML("<b>TML</b> <a href = 'https://www.degruyter.com/document/doi/10.1515/CCLM.2007.249/html'>(Lit.)</a>"),
                         HTML("<b>kosmic</b> <a href = 'https://www.nature.com/articles/s41598-020-58749-2'>(Lit.)</a>"),
                         HTML("<b>reflimR</b> <a href = 'https://doi.org/10.1515/labmed-2023-0042'>(Lit.)</a>")
                         ),
                         choiceValues = list("refiner", "tmc", "tml", "kosmic", "reflimr")
                         )
                       ),
              conditionalPanel(
                condition = "input.methodradio == 'tml'", 
                conditionalPanel(
                  condition = "output.advanced == 1",
                  checkboxInput("fasttml", 
                              HTML("Fast mode <i style = 'font-size:80%;'>(3 significant figures)</i>"),
                              value = TRUE
                              )
                )
                ),
              conditionalPanel(
                condition = "output.advanced == 1",
                conditionalPanel(
                  condition = "input.methodradio == 'refiner'",
                  checkboxInput("modboxcox", 
                                HTML("<i style = 'font-size:80%;line-height: 10px;'>use modified Box-Cox</i>"), 
                                value = FALSE),
                  sliderInput("nbootstrap", NULL, min=0, max=50, value=0, step=5, ticks=FALSE),
                HTML("<div style = 'text-align:center;line-height: 10px;font-size:90%'><i>bootstrap iterations</i></div>")
                )),
              checkboxInput(
                "compare", 
                "Compare RI limits", 
                value = FALSE),
              conditionalPanel(condition = "input.compare==1",
                               fluidRow(
                                 column(
                                   6,
                                   div(
                                     id = "ref_low_container",
                                   numericInput("referencelimits.low", NULL, "0.0", step = 0.1))),
                                 column(
                                   6,
                                   div(
                                     id = "ref_high_container",
                                     numericInput("referencelimits.high", NULL, "0.0", step = 0.1)))),
                               HTML("<div style='text-align:center;width:100%;font-size:80%'><i>lower limit / upper limit </i>")
                               )
              ),
            card(
              max_height = 400,
              plotOutput("plot")), 
            card(
              max_height = 300,
              htmlOutput("placeholder")))
            ))
          )
        )
      )
    )
  


# Server initialization
server <- function(input, output, session) {
  setwd(tempdir())
  dataframe1 = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist))
  dataframe1$result = as.character(dataframe$result)
  dataframe1$age = as.numeric(dataframe$age)
  output$table <- renderRHandsontable(rhandsontable(dataframe1, width = '400', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                        hot_col("result", validator = resultvalidator) %>%
                                        hot_col("sex", allowInvalid = TRUE)
  )
  
  
  raw_data <- reactive({
    dataframe <- hot_to_r(req(input$table))
    if (is.null(dataframe)) return(dataframe1)
    dataframe$result <- gsub("^<", "", dataframe$result) 
    dataframe$result <- gsub(",", ".", dataframe$result, fixed = TRUE)
    dataframe$result <- as.numeric(dataframe$result)
    dataframe$age <- gsub(",", ".", dataframe$result, fixed = TRUE)
    dataframe$age <- as.numeric(dataframe$age)
    dataframe <- dataframe[dataframe$result > 0 & !is.na(dataframe$result), ]
    femalelist = c(femalelist, input$Init_female)
    malelist = c(malelist,input$Init_male)
    diverselist = c(diverselist,input$Init_diverse)
    dataframe <- dataframe %>% mutate(sex = ifelse(sex %in% femalelist, 'F', ifelse(sex %in% malelist, 'M', ifelse(sex %in% diverselist, 'D', 'X'))))
    dataframe <- dataframe[!is.na(dataframe$sex), ]
    if ("trimester" %in% colnames(dataframe)) dataframe$trimester = plyr::mapvalues(dataframe$trimester,c(""),0)
    if (input$remove_pregnancy) dataframe = dataframe[,-4]
    return(dataframe)
  })
  
  
  raw_fasttml <- reactive({
    fasttml = input$fasttml
    if (is.null(fasttml)) fasttml == T
    else if (fasttml==F) fasttml==F
    else fasttml==T
    return(fasttml)
  })
  
  
  # Set initial Output conditions  
  
  output$advanced <- renderText({
    '0'
  })
  outputOptions(output, "advanced", suspendWhenHidden=FALSE)
  
  
  
  
  output$pregnancymode <- renderText({
    '0'
  })
  outputOptions(output, "pregnancymode", suspendWhenHidden=FALSE)

  output$fasttml = renderPrint({
    raw_fasttml()
  })
  outputOptions(output, "fasttml", suspendWhenHidden=FALSE)
  
  # Input advanced mode
  observeEvent(input$advanced_button, {
    updateBox("initbox", action = "toggle")
    output$advanced = renderText({
      '1'
      })
  })
  
  observeEvent(input$standard_button, {
    output$advanced = renderText({
      '0'
    })
  })
  
  
  # Input Observations oben/close boxes  
  observeEvent(input$DI, {
    if (input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    else (updateBox("boxtable", action = "toggle"))
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
  })
  
  observeEvent(input$SV, {
    if (input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    else (updateBox("strat_boxplot", action = "toggle"))
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
  })
  
  observeEvent(input$RI, {
    if (input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    else (updateBox("boxplot", action = "toggle"))
    if (!input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
  })
  
  
  
  # Input Table
  observeEvent(input$table, {
    dataframe = hot_to_r(input$table)
    input_sex= unique(dataframe$sex)
    updateSelectInput(
      session = getDefaultReactiveDomain(),
      inputId = "Init_female",
      label = "female:",
      choices= input_sex,
      selected = "F"
    )
    updateSelectInput(
      session = getDefaultReactiveDomain(),
      inputId = "Init_male",
      label = "male:",
      choices= input_sex,
      selected = "M"
    )
    updateSelectInput(
      session = getDefaultReactiveDomain(),
      inputId = "Init_diverse",
      label = "diverse:",
      choices= input_sex,
      selected = "D"
    )
    comparesexlist = c(femalelist, input$Init_female, malelist, input$Init_male, diverselist, input$Init_diverse, "", NA)
    if (all(input_sex %in% comparesexlist) == FALSE) shinyalert("Check sex column variables...", "Please assign the data manually") 
    if (all(input_sex %in% comparesexlist) == FALSE) updateBox("initbox", action = "toggle") 
  })
  
  
  
  # Clear data and table
  observeEvent(input$clear, {
    
    if (input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    
    dataframe = dataframe1
    updateNumericInput(session, "referencelimits.low", value = 0)
    updateNumericInput(session, "referencelimits.high", value = 0)
    updateNumericInput(session, "agell", value = "")
    updateNumericInput(session, "ageul", value = "")
    output$pregnancymode <- renderText({ '0' })
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '400', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                          hot_col("result", validator = resultvalidator) %>%
                                          hot_col("sex", allowInvalid = TRUE)
    )

    output$strat_plot <- renderPlot({
      ""
    })
    output$strat_placeholder <- renderText({
      ""
    })
    output$plot <- renderPlot({
      ""
    })
    output$placeholder <- renderText({
      ""
    })
    
  })
  
  # Observations of numericInputs: Input lower value > Input higher value => red marking added
    observe({
      strat_agell <- input$strat_agell
      strat_ageul <- input$strat_ageul
      
      if (!is.null(strat_agell) && !is.na(strat_agell) &&
          !is.null(strat_ageul) && !is.na(strat_ageul)) {
        
        if (strat_agell > strat_ageul) {
          shinyjs::addClass("strat_agell_container", "red-border")
          shinyjs::addClass("strat_ageul_container", "red-border")
        } else {
          shinyjs::removeClass("strat_agell_container", "red-border")
          shinyjs::removeClass("strat_ageul_container", "red-border")
        }
      } else {
        shinyjs::removeClass("strat_agell_container", "red-border")
        shinyjs::removeClass("strat_ageul_container", "red-border")
      }
    })
    
    observe({
      agell <- input$agell
      ageul <- input$ageul
      
      if (!is.null(agell) && !is.na(agell) &&
          !is.null(ageul) && !is.na(ageul)) {
        
        if (agell > ageul) {
          shinyjs::addClass("agell_container", "red-border")
          shinyjs::addClass("ageul_container", "red-border")
        } else {
          shinyjs::removeClass("agell_container", "red-border")
          shinyjs::removeClass("ageul_container", "red-border")
        }
      } else {
        shinyjs::removeClass("agell_container", "red-border")
        shinyjs::removeClass("ageul_container", "red-border")
      }
    })
    
    observe({
      ref_low <- input$referencelimits.low
      ref_high <- input$referencelimits.high
      
      if (!is.null(ref_low) && !is.na(ref_low) &&
          !is.null(ref_high) && !is.na(ref_high)) {
        
        if (ref_low > ref_high) {
          shinyjs::addClass("ref_low_container", "red-border")
          shinyjs::addClass("ref_high_container", "red-border")
        } else {
          shinyjs::removeClass("ref_low_container", "red-border")
          shinyjs::removeClass("ref_high_container", "red-border")
        }
      } else {
        shinyjs::removeClass("ref_low_container", "red-border")
        shinyjs::removeClass("ref_high_container", "red-border")
      }
    })
  
  # Generate demo data
  observeEvent(input$demo, {
    
    if (input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    
    dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist))
    dataframe$result = as.numeric(dataframe$result)
    dataframe$age = as.numeric(dataframe$age)
    demosamplesize <- tablesize * 0.25
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
    updateNumericInput(session, "agell", value = 18)
    updateNumericInput(session, "ageul", value = 45)
    updateNumericInput(session, "strat_agell", value = 18)
    updateNumericInput(session, "strat_ageul", value = 45)
    output$pregnancymode <- renderText({ '0' })
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '400', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                          hot_col("result", validator = resultvalidator) %>%
                                          hot_col("sex", allowInvalid = TRUE)
    )
  })
  
  # Add trimester column
  observeEvent(input$pregnancy_button, {
    dataframe = data.frame(result = rep(NA, tablesize), age = rep(NA, tablesize), sex = factor(rep(NA, tablesize), levels = sexlist), trimester = factor(rep(NA, tablesize), levels = c(0,1, 2, 3)))
    dataframe_raw=raw_data()
    dataframe$result <- c(dataframe_raw$result, rep(NA, tablesize)- length(dataframe_raw$result))[1:tablesize]
    dataframe$age <- c(dataframe_raw$age, rep(NA, tablesize)- length(dataframe_raw$age))[1:tablesize]
    dataframe$sex <- c(dataframe_raw$sex, rep(NA, tablesize)- length(dataframe_raw$sex))[1:tablesize]
    dataframe$trimester == 0
    updateNumericInput(session, "referencelimits.low", value = 0)
    updateNumericInput(session, "referencelimits.high", value = 0)
    updateNumericInput(session, "agell", value = "")
    updateNumericInput(session, "ageul", value = "")
    output$pregnancymode <- renderText({ '1' })
    updateBox("boxtable", action="update")
    updateBox("strat_boxplot", action="update")
    updateBox("boxplot", action="update")
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '450', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                          hot_col("result", validator = resultvalidator) %>%
                                          hot_col("sex", allowInvalid = TRUE)
    )
  })
  
  observeEvent(input$remove_pregnancy, {
    dataframe= raw_data()
    dataframe = dataframe[,-4]
    output$table <- renderRHandsontable(rhandsontable(dataframe, width = '400', height = 550, stretchH = "all", rowHeaderWidth = 65) %>%
                                          hot_col("result", validator = resultvalidator) %>%
                                          hot_col("sex", allowInvalid = TRUE)
    )
    output$pregnancymode <- renderText({ '0' })
  })
    
  # Visualize data for sex differences
  observeEvent(input$sexbox, {
    
    if (input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    
    dataframe=raw_data()
    
    updateRadioButtons(session,"trimesterCalc", selected= isolate(input$sexradio))
    
    if (length(nrow(dataframe$age))==0) {
      agell <- 0
      ageul <- 0
    }
    else {
      agell <- isolate(input$strat_agell)
      ageul <- isolate(input$strat_ageul)
    }
    agelimitsvalid <- (ageul > 0 && (ageul > agell) && !is.na(agell) && !is.na(ageul)) 
    
    updateNumericInput(session, "agell", value = isolate(input$strat_agell))
    updateNumericInput(session, "ageul", value = isolate(input$strat_ageul))
    
    sexradio<- isolate(input$sexradio)
    if (sexradio == 'M') dataframe <- dataframe[dataframe$sex == 'M', ]
    else if (sexradio == 'F') dataframe <- dataframe[dataframe$sex == 'F', ]
    else if (sexradio == 'D') dataframe <- dataframe[dataframe$sex == 'D', ]
    
    updateRadioButtons(session,"sexradioCalc", selected= isolate(input$sexradio))
  
    cases = length(dataframe$result)
    dataframe <- dataframe %>% group_by(sex) %>% filter(n() >= 100)
    casesfiltered = length(dataframe$result)
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
      output$strat_plot <- renderPlot({
        ggplot(dataframe, aes(sex, result)) + ylim(ylim.min, ylim.max) + labs(x = NULL) +
          theme_minimal() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
          geom_violin(fill = "#40668d", alpha = 0.2) +
          geom_boxplot(fill = "#40668d", outlier.shape = NA, width = 0.1)
      })
      output$strat_placeholder = renderText({
        paste("<br>",
              kable(mediantable, "html", align = "c") %>%  kable_styling(),
              "<br><p style='font-size: 14px;'><b>Selected data:</b>",
              "n = ", casesfiltered, "<br>",
              if (dataremoved) "Groups with n<100 were removed<br>",
              if (agelimitsvalid==T) paste("age: from ", agell, " to ", ageul)
              else "age: no selection",
              "<br>",
              if (pvalue < 0) "More than one group needed for statistical evaluation."
              else if (pvalue < 0.01) paste("<b style='color: red;'>Significant difference </b> (p<0.01, ", test_used, ")<br>", sep = "")
              else if (pvalue < 0.05) paste("<b style='color: red;'>Significant difference</b> (p=", pvalue, ", ", test_used, ")<br>", sep = "")
              else paste("<b style='color: #10D010;'>No significant difference</b> (p=", pvalue, ", ", test_used, ")<br>", sep = ""),
              "<br><b>Max. deviation of the medians: </b> ",
              if (mediandiff > pu_percent) "<b style='color: red;'>"
              else if (mediandiff > (0.5 * pu_percent)) "<b style='color: orange;'>"
              else "<b style='color: #10D010;'>",
              mediandiff, "%</b><br>",
              "Estimated permissible uncertainty (pU): ", round(pu_percent,2), "%</p>",
              br()
        )
      })
    }
  })
  

  # Visualize drift
  observeEvent(input$drift, {
    
    if (input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    if (!input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")

    dataframe=raw_data()
    
    sexradio<- isolate(input$sexradio)
    if (sexradio == 'M') dataframe <- dataframe[dataframe$sex == 'M', ]
    else if (sexradio == 'F') dataframe <- dataframe[dataframe$sex == 'F', ]
    else if (sexradio == 'D') dataframe <- dataframe[dataframe$sex == 'D', ]
    
    updateRadioButtons(session,"sexradioCalc", selected= isolate(input$sexradio))
    
    agell <- isolate(input$strat_agell)
    ageul <- isolate(input$strat_ageul)
    agelimitsvalid <- (ageul > 0 && (ageul > agell) && !is.na(agell) && !is.na(ageul))  
      
    updateNumericInput(session, "agell", value = isolate(input$strat_agell))
    updateNumericInput(session, "ageul", value = isolate(input$strat_ageul))
    
    cases = length(nrow(dataframe$age))
    
    # Alert when n<500 
    if (cases < 500) shinyalert("insufficient data...", paste("The selected dataset contains ", cases, " results with given age.\nPlease increase the sample size (minimum n=500)."), type = "error")
    else {
      # Subsampling and reduction of sample size for n>10000
      dataframe[sample(nrow(dataframe)), ]
      if (cases > 10000) dataframe <- dataframe[sample(1:length(dataframe$result), 10000), ]

        if (agelimitsvalid) dataframe <- dataframe[(dataframe$age >= agell) & (dataframe$age <= ageul), ]
        
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
      
      output$strat_plot <- renderPlot({ stratplot })
      
      output$strat_placeholder = renderText({
        paste(
          "<p style='font-size: 14px;'><b>Selected data:</b><br>",
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
      }
    
  })

    
  # Start calculation
  observeEvent(input$calc, {
    if (input$boxplot$collapsed) updateBox("boxplot", action = "toggle")
    if (!input$strat_boxplot$collapsed) updateBox("strat_boxplot", action = "toggle")
    if (!input$boxtable$collapsed) updateBox("boxtable", action = "toggle")
    
    dataframe = raw_data()
    
    sexradioCalc<- isolate(input$sexradioCalc)
    if (sexradioCalc == 'M') dataframe <- dataframe[dataframe$sex == 'M', ]
    else if (sexradioCalc == 'F') dataframe <- dataframe[dataframe$sex == 'F', ]
    else if (sexradioCalc == 'D') dataframe <- dataframe[dataframe$sex == 'D', ]
    
    if ("trimester" %in% colnames(dataframe)) trimesterCalc = isolate(input$trimesterCalc)
    if (all(dataframe$trimester==0,na.rm=T)) trimesterCalc = 0
    else dataframe = dataframe[dataframe$trimester == trimesterCalc ]
  
    if (length(nrow(dataframe$age))==0) {
      agell <- 0
      ageul <- 0
    }
    else {
      agell <- isolate(input$strat_agell)
      ageul <- isolate(input$strat_ageul)
    }
    agelimitsvalid <- (ageul > 0 && (ageul > agell) && !is.na(agell) && !is.na(ageul))  
    if (agelimitsvalid) dataframe <- dataframe[(dataframe$age >= agell) & (dataframe$age <= ageul), ]
    referencelimits.low <- isolate(input$referencelimits.low)
    referencelimits.high <- isolate(input$referencelimits.high)
    referencelimitsvalid <- (referencelimits.high > 0 && (referencelimits.high > referencelimits.low) && !is.na(referencelimits.low) && !is.na(referencelimits.high))  
    estimatedlimits.low <- 0
    estimatedlimits.high <- 0
    citext <- ""
    skewness <- ""
    
    if (length(na.omit(dataframe$result)) < 500) shinyalert("insufficient data...", paste("The selected dataset contains ", length(na.omit(dataframe$result)), " results.\nPlease increase the sample size (minimum n=500)."), type = "error")
    else {
      if (length(na.omit(dataframe$result)) < 2000) shinyalert("small sample size...", paste("The selected dataset contains only ", length(dataframe$result), " results.\nYou have been warned..."), type = "warning")
      if (length(na.omit(dataframe$result)) > 100000) shinyalert("large sample size...", paste("The selected dataset contains ", length(dataframe$result), " results.\nYou have been warned..."), type = "warning")
    
      methodradio <- isolate(input$methodradio)    
      if (methodradio == 'refiner') {
        nbootstrap <- isolate(input$nbootstrap)
        modboxcox <- isolate(input$modboxcox)
        if (modboxcox) boxcoxmode <- "modBoxCox"
        else boxcoxmode <- "BoxCox"
        resri <- findRI(Data = na.omit(dataframe$result), model=boxcoxmode, NBootstrap=nbootstrap)
        if (nbootstrap > 0) 
          citext <- paste("95% confidence intervals: Lower limit (", round(getRI(resri)[1, 3], 2), " - ", round(getRI(resri)[1, 4], 2),
                          ") - Upper limit ( ", round(getRI(resri)[2, 3], 2), " - ", round(getRI(resri)[2, 4], 2), ")<br>", sep="")
        skewness <- paste("</p><br><p style='font-size: 12px;'>Estimated skewness:", round(skewness(na.omit(dataframe$result)),2))
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
        fasttml=raw_fasttml()
     
        if (quantile(dataframe$result, probs = 0.90) >= 1000) decimalcount <- 0
        else if (quantile(dataframe$result, probs = 0.90) >= 100) decimalcount <- 1
        else if (quantile(dataframe$result, probs = 0.90) >= 10) decimalcount <- 2
        else decimalcount <- 3
        if (fasttml == TRUE && decimalcount > 0) decimalcount <- decimalcount - 1
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
          "<p style='font-size: 14px;'><b>Selected data:</b><br>",
          "n = ", length(na.omit(dataframe$result)),
          if (agelimitsvalid) paste("<br>age: from ", agell, " to ", ageul)
          else paste("<br>age: no selection"),
          "<br>sex: ",
          if (sexradioCalc == 'M') "male"
          else if (sexradioCalc == 'F') "female"
          else if (sexradioCalc == 'D') "non-binary"
          else "no selection",
          if ("trimester" %in% colnames(dataframe)) if (trimesterCalc== 0) "<br>no trimester in dataset " else if (trimesterCalc== 1)"<br>trimester: 1. " else if (trimesterCalc== 2)"<br>trimester: 2. " else if (trimesterCalc== 3)"<br>trimester: 3. ",
          "<br>method: ",
          if (methodradio == 'refiner') "refineR"
          else if (methodradio == 'tmc') "TMC"
          else if (methodradio == 'tml') "TML"
          else if (methodradio == 'kosmic') "kosmic"
          else if (methodradio == 'reflimr') "reflimR",
          "<br><br><b>Estimated reference limits: ", estimatedlimits.lowtxt, " - ", estimatedlimits.hightxt, "</b><br>",
          if (citext != "") citext,
          if (referencelimitsvalid) paste("<b style='color: #D0D0D0;'>Comparison reference limits: ", referencelimits.low.txt, " - ", referencelimits.high.txt, "</b>"),
          if (skewness != "") skewness,
          "</p>"
        )
      })
    }
  })
}

shinyApp(ui, server)
