# this script render the UI of analysis page of module: prognostic model
# cancer.datasets <- openxlsx::read.xlsx("./doc/datasetInfo/Bladder.xlsx",sheet = 1)

# set up with UI
ui.page_module_prognostic_model <- function() {
  fluidPage(
    style = "width:80%;",
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "./static/css/dataset.css")
    ),

    #some space on the head
    tags$br(),

    fluidRow(
      div(
        column(
          7,
          fluidRow(
            column(2, img(src = "model3.png", width = 125), style = "padding-top:10px;"),
            column(10,
                   h2("Customized Modelling", style = "color:#ff6b6b; font-weight: bold;"),
                   p("Users can fit a survival machine learning model based on a set of genes and explore its clinical value and biological implications.",
                     style="text-align:justify; "))
          ),
        ),
        column(
          5
        ),

        # background image
        style = "background-image:url(\"./module_background.png\");background-attachment: scroll; background-size: 1450px 150px"
      )
    ),
    tags$hr(),
    tags$br(),

    #### User input ####
    fluidRow(
      ##### event #####
      column(
        4,
        bs4Dash::bs4Card(
          fluidRow(
            p("Select an event outcome:"),

            shinyWidgets::awesomeRadio(
              inputId = "AnalysisPagePrognosticModelModuleEventOutcomeSelection",
              label = NULL,
              choices =
                c(
                  "Overall survival",
                  "Relapse-free survival",
                  "Disease-free survival",
                  "Progression-free survival",
                  "Disease-specific survival"
                ),
              selected = "",
              status = "info",
              checkbox = FALSE,
              width = "100%"
            )
          ),

          title = "1. Prognostic Outcomes",
          collapsible = FALSE,
          collapsed = FALSE,
          background = "gray",
          closable = FALSE,
          maximizable = FALSE,
          icon = icon("wheelchair"),
          width = 12,
          height = "300px",
          id = "AnalysisPagePrognosticModelModuleUserInput2"
        )
      ),

      ##### cohort #####
      column(
        4,
        bs4Dash::bs4Card(
          p("Determine datasets with the survival outcome you selected:"),

          # datasets
          div(
            style = "height:140px;overflow-y:scroll;overflow-x:scroll;",
            shinyWidgets::checkboxGroupButtons(
              inputId = "AnalysisPagePrognosticModelModuleDatasetGroup",
              label = NULL,
              choiceNames = "",
              choiceValues = "",
              individual = TRUE,
              direction = "horizontal",
              checkIcon = list(yes = icon("dot-circle-o"), no = icon("circle-o"))
            )
          ),

          tags$br(),

          ## select/disselect all btn
          fluidRow(
            # select all
            column(
              6,
              shiny::actionButton(
                inputId = "AnalysisPagePrognosticModelModuleDatasetGroupSelectAllBtn",
                label = "Select all",
                icon = icon("check"),
                width = "100%",
                style = "color:#ef5350; background-color:white;"
              )
            ),
            # disselect all
            column(
              6,
              shiny::actionButton(
                inputId = "AnalysisPagePrognosticModelModuleDatasetGroupDisSelectAllBtn",
                label = "Disselect all",
                icon = icon("close"),
                width = "100%",
                style = "color:#ef5350; background-color:white;"
              )
            )
          ),

          # card settings
          title = "2. Available Datasets",
          collapsible = FALSE,
          collapsed = FALSE,
          background = "gray",
          closable = FALSE,
          maximizable = FALSE,
          icon = icon("group"),
          width = 12,
          height = "300px",
          id = "AnalysisPagePrognosticModelModuleUserInput3"
        )
      ),

      ##### genelist  #####
      column(
        4,
        bs4Dash::bs4Card(
          fluidRow(
            p("Input your interest gene list:"),

            ## genelist
            shiny::textAreaInput(
              inputId = "AnalysisPagePrognosticModelModuleGenelistTextArea",
              label = NULL,
              placeholder = "CD8A\nCD8B\nGZMA\nGZMB\nPRF1\n",
              resize = "vertical",
              width = "100%",
              height = "80px"
            ),

            ## gene number summary
            htmlOutput(
              outputId = "AnalysisPagePrognosticModelModuleGeneNumberStat"
            ),

            ## genelist name
            shiny::textInput(
              inputId = "AnalysisPagePrognosticModelModuleGenelistNameInput",
              label = NULL,
              placeholder = "Enter your custom gene list name",
              width = "100%"
            ),

            # submit button
            fluidRow(
              column(1),
              column(
                10,

                shiny::actionButton(
                  inputId = "AnalysisPagePrognosticModelModuleGenelistSubmitBtn",
                  label = "Submit",
                  icon = icon("check"),
                  width = "100%",
                  style = "background-color:#40916C; color:#FFFFFF;"
                )
              ),
              column(1)
            )

          ),

          # card information
          title = "3. Gene List",
          collapsible = FALSE,
          collapsed = FALSE,
          background = "gray",
          closable = FALSE,
          maximizable = FALSE,
          icon = icon("google"),
          width = 12,
          height = "300px",
          id = "AnalysisPagePrognosticModelModuleUserInput1"
        ),
      ),
    ),

    tags$br(),

    div(
      fluidRow(
        ##### statistics #####
        column(
          4,

          bs4Dash::box(
            p("Summary of gene list across all selected datasets:"),

            # p("Number of genes in data sets"),

            div(
              fluidRow(
                column(
                  12,
                  div(
                    p("awaiting for user input...",
                      style = "font-size:80%; color:#3b3b3b;"),
                    style = "text-align:center;"
                  )
                )
              ),
              fluidRow(
                column(
                  12,
                  div(
                    img(src = "static/picture/waiting2.gif",
                        alt = "waiting for calculation",
                        style = "height:250px;"),

                    # div style to make img at center
                    style = "display:flex; align-items: center;justify-content: center;opacity:0.5;",
                  )
                )
              ),

              id = "AnalysisPagePrognosticModelModuleDatasetSummaryPlotDivWaiting"
            ),

            # barplot
            div(
              div(
                htmlOutput(
                  outputId = "AnalysisPagePrognosticModelModuleDatasetAndGeneSummaryBarplotHtml"
                ),
                style = "overflow-y:scroll; height:220px;"
              ),

              id = "AnalysisPagePrognosticModelModuleDatasetAndGeneSummaryBarplotDiv",
              class = "hidden-element"
            ),


            # stat information
            tags$br(),

            div(
              useShinyjs(),

              # p("Data sets with missing genes"),

              div(
                verbatimTextOutput(
                  outputId = "AnalysisPagePrognosticModelModuleDatasetAndGeneSummaryText"
                ),

                style = "overflow-y:scroll; height:140px;overflow-x:auto;"
              ),

              id = "AnalysisPagePrognosticModelModuleDatasetAndGeneSummaryTextDiv",
              class = "hidden-element"
            ),


            # card information
            title = "4. Gene Summary",
            collapsible = FALSE,
            collapsed = FALSE,
            background = "gray",
            closable = FALSE,
            maximizable = FALSE,
            icon = icon("calculator"),
            width = 12,
            height = "480px",
            id = "AnalysisPagePrognosticModelModuleUserInput4"
          )
        ),

        ##### dataset #####
        column(
          4,
          bs4Dash::bs4Card(

            ### datasets split to training and validation
            fluidRow(
              # training datasets
              fluidRow(
                p("Split into training or validation datasets:"),

                #tags$p("Traning set",style="font-size:120%;"),
                column(
                  12,

                  # use dropdown picker instead of text input
                  selectInput(
                    inputId = "AnalysisPagePrognosticModelModuleTrainDatasetSelection",
                    label = "Training dataset",
                    width = "100%",
                    choices = c("")
                  ),
                  tags$hr()
                )
              )
            ),


            p("Validation datasets"),
            # validation set
            # textAreaInput(
            #   inputId = "AnalysisPagePrognosticModelModuleValidationDatasetSelection",
            #   label = "Validation datasets",
            #   height = "150px",
            #   resize = "none",
            #   width = "100%"
            # ),

            div(
              verbatimTextOutput(
                outputId = "AnalysisPagePrognosticModelModuleValidationDataset"
              ),
              style = "height:180px; overflow-y:scroll;"
            ),
            tags$br(),
            p("* Once the training dataset is determined, the remaining datasets are automatically assigned to the validation datasets.",
              style = "color:#3b3b3b; font-size:80%; font-style:italic;"),

            # card information
            title = "5. Training and Validation",
            collapsible = FALSE,
            collapsed = FALSE,
            background = "gray",
            closable = FALSE,
            maximizable = FALSE,
            icon = icon("database"),
            width = 12,
            height = "480px",
            id = "AnalysisPagePrognosticModelModuleUserInput7"
          )
        ),

        ##### method #####
        column(
          4,

          useShinyjs(),

          bs4Dash::box(
            p("Choose a survival machine-learning algorithm:"),

            # pick a method
            selectInput(
              inputId = "AnalysisPagePrognosticModelModuleConstructionMethodSelection",
              label = NULL,
              choices = c(
                "StepCox",
                "RSF",
                "Lasso",
                "Ridge",
                "Enet",
                "GBM",
                "SVM",
                "plsRcox",
                "Coxboost",
                "SuperPC"
              ),
              width = "100%"
            ),

            tags$hr(),

            # algorithm settigns
            div(
              span(
                span("Tuning parameters",
                     style = "color:#000;font-size:110%;"),
                shiny::actionButton(
                  inputId = "PrognosticModelTuningParametersHelpInfoBtn",
                  label = NULL,
                  icon = icon("question-circle"),
                  style = "background-color:#fff;border-color:#fff;"
                )
              ),

              # shinyBS::bsTooltip(
              #   id = "PrognosticModelTuningParametersHelpInfoBtn",
              #   title = "help info",
              #   trigger = "hover",
              #   placement = "bottom",
              #   options = list(container = "body")
              # ),

              div(
                ## default show stepcox
                ###### stepcox #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),
                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_StepCox_direction",
                        label = NULL,
                        choices = c(
                          "forward",
                          "backward",
                          "both"
                        ),
                        width = "100%"
                      ),
                    ),
                    column(
                      6,
                      p("StepCox direction",
                        class = "rightInputTitle"),
                    )
                  ),

                  id = "PMStepCoxSettingDiv"
                ),

                ###### RSF ######
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      numericInput(
                        inputId = "PM_RSF_nodesize",
                        label = NULL,
                        min = 3,
                        max = 30,
                        value = 5,
                        step = 1,
                        width = "100%"
                      ),
                    ),
                    column(
                      6,
                      p("RSF nodesize",
                        class = "rightInputTitle")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      numericInput(
                        inputId = "PM_RSF_nsplit",
                        label = NULL,
                        min = 2,
                        max = 20,
                        value = 10,
                        step = 1,
                        width = "100%"
                      ),
                    ),
                    column(
                      6,
                      p("RSF nsplit",
                        class="rightInputTitle")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_RSF_splitrule",
                        label = NULL,
                        choices = c(
                          "logrank",
                          "bs.gradient",
                          "logrankscore"
                        ),
                        multiple = FALSE,
                        width = "100%"
                      )
                    ),
                    column(
                      6,
                      p("RSF splitrule",
                        class = "rightInputTitle")
                    )
                  ),


                  id = "PMRSFSettingDiv",
                  class = "hidden-element"
                ),

                ###### Lasso #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_Lasso_lamda_rule",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "lambda.min",
                          "lambda.1se"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("Lasso lamda rule",
                        class = "rightInputTitle")
                    )
                  ),

                  id = "PMLassoSettingDiv",
                  class = "hidden-element"
                ),

                ###### Ridge #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_Ridge_lamda_rule",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "lambda.min",
                          "lambda.1se"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("Ridge lamda rule",
                        class = "rightInputTitle")
                    )
                  ),


                  id = "PMRidgeSettingDiv",
                  class = "hidden-element"
                ),

                ###### Enet #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      numericInput(
                        inputId = "PM_Enet_alpha",
                        label = NULL,
                        min = 0.1,
                        max = 0.9,
                        value = 0.5,
                        width = "100%",
                        step = 0.1
                      ),
                    ),
                    column(
                      6,
                      p("Enet alpha",
                        class = "rightInputTitle")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_Enet_lamda_rule",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "lambda.min",
                          "lambda.1se"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("Enet lamda rule",
                        class = "rightInputTitle")
                    )
                  ),

                  id = "PMEnetSettingDiv",
                  class = "hidden-element"
                ),

                ###### GBM #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      numericInput(
                        inputId = "PM_GBM_nodesize",
                        label = NULL,
                        min = 3,
                        max = 30,
                        value = 5,
                        width = "100%",
                        step = 1
                      ),
                    ),
                    column(
                      6,
                      p("GBM nodesize",
                        class = "rightInputTitle")
                    )
                  ),


                  id= "PMGBMSettingDiv",
                  class = "hidden-element"
                ),

                ###### SVM #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_SVM_type",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "vanbelle1",
                          "vanbelle2",
                          "regression",
                          "hybrid"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("SVM type",
                        class = "rightInputTitle")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_SVM_diffmeth",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "makediff1",
                          "makediff2",
                          "makediff3"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("SVM diffmeth",
                        class = "rightInputTitle")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_SVM_optmeth",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "quadprog",
                          "ipop"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("SVM optmeth",
                        class = "rightInputTitle")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_SVM_kernel",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "add_kernel",
                          "lin_kernel",
                          "rbf_kernel",
                          "poly_kernel"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("SVM kernel",
                        class = "rightInputTitle")
                    )
                  ),


                  id = "PMSVMSettingDiv",
                  class = "hidden-element"
                ),

                ###### plsRcox #####
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_plsRcox_lambda_rule",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "lambda.min",
                          "lambda.1se"
                        ),
                        multiple = FALSE
                      )
                    ),
                    column(
                      6,
                      p("plsRcox lambda rule",
                        class = "rightInputTitle")
                    )
                  ),

                  id = "PMplsRcoxSettingDiv",
                  class = "hidden-element"
                ),

                ###### Coxboost ######
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      selectInput(
                        inputId = "PM_Coxboost_type",
                        label = NULL,
                        width = "100%",
                        choices = c(
                          "verweij",
                          "naive"
                        ),
                        multiple = FALSE
                      ),
                    ),
                    column(
                      6,
                      p("Coxboost type",
                        class = "rightInputTitle")
                    )
                  ),



                  id = "PMCoxboostSettingDiv",
                  class = "hidden-element"
                ),

                ####### SuperPC ######
                div(
                  # p("Algorithm settings"),
                  # tags$hr(),

                  fluidRow(
                    column(
                      6,
                      sliderInput(
                        inputId = "PM_SuperPC_ncomponents",
                        label = NULL,
                        min = 1,
                        max = 3,
                        step = 1,
                        value = 1,
                        animate = TRUE,
                        width = "100%"
                      ),
                    ),
                    column(
                      6,
                      p("SuperPC ncomponents",
                        class = "rightInputTitle")
                    )
                  ),



                  id = "PMSuperPCSettingDiv",
                  class = "hidden-element"
                ),

                style = "width:90%;"
              ),

              id = "PMMethodSettingDiv",
              style = "height:300px; overflow-y:scroll;"
            ),

            tags$br(),

            # card information
            title = "6. Machine-Learning",
            collapsible = FALSE,
            collapsed = FALSE,
            background = "gray",
            closable = FALSE,
            maximizable = FALSE,
            icon = icon("cogs"),
            width = 12,
            height = "480px",
            id = "AnalysisPagePrognosticModelModuleUserInput6"
          )
        )
      ),

      ##### modeling btn #####
      tags$br(),
      tags$br(),

      fluidRow(
        column(3),
        column(
          6,
          # analyze button
          shinyWidgets::actionBttn(
            inputId = "AnalysisPagePrognosticModelModuleStartAnalyzeBtn",
            label = "Everything is ready! Start Modeling",
            icon = icon("rocket"),
            color = "danger",
            style = "unite",
            size = "md",
            block = TRUE
          )
        ),
        column(3)
      ),
      tags$br(),

      id = "AnalysisPagePrognosticModelModuleUserInputRow2Div",
      class = "hidden-element"
    ),

    ##### Model result ####
    div(
      tags$hr(),

      ## results
      fluidRow(

        column(1),

        # score plot
        column(
          6,

          div(
            div(
              div(
                p("Gene Significance",
                  style = "text-align:center;font-weight:bold;padding-top:6px;padding-bottom:6px;color:#fff;"),
                style = "background-color:#0096c7;"
              ),


              # p("Risk score among all data sets"),
              div(
                htmlOutput(
                  outputId = "AnalysisPagePrognosticModelModuleRiskScorePlotHtml"
                ),
                style = "overflow-x:auto;margin-left:10px;margin-right:10px;"
              ),

              tags$br(),

              style = "height:450px;box-shadow:2px 2px 3px 2px #ccc;"
            ),
            style = "margin-right:4px;"
          )
        ),


        column(
          4,

          div(
            div(
              p("Performance Summary",
                style = "text-align:center;font-weight:bold;padding-top:6px;padding-bottom:6px;color:#fff;"),
              style = "background-color:#0096c7;"
            ),

            div(
              htmlOutput(
                outputId = "AnalysisPagePrognosticModelModulePlot1Html"
              ),
              style = "height:380px; overflow-y:auto;margin-left:10px;margin-right:10px;"
            ),

            style = "height:450px; box-shadow:2px 2px 3px 2px #ccc;"
          ),

        ),

        column(1)
      ),

      tags$br(),
      fluidRow(
        column(1),
        column(
          10,
          span(
            shiny::downloadButton(
              outputId = "AnalysisPagePrognosticModelModuleRiskScorePlotDownloadPdfBtn",
              label = "Download PDF",
              icon = icon("file-pdf-o"),
              style = "color:#FF6B6B;background-color:#ffffff;",
            ),
            shiny::downloadButton(
              outputId = "AnalysisPagePrognosticModelModuleRiskScorePlotDownloadPngBtn",
              label = "Download PNG",
              icon = icon("file-image-o"),
              style = "color:#FF6B6B;background-color:#ffffff;",
            ),
            shiny::downloadButton(
              outputId = "AnalysisPagePrognosticModelModuleRiskScorePlotDownloadTableBtn",
              label = "Download Table",
              icon = icon("file-image-o"),
              style = "color:#FF6B6B;background-color:#ffffff;",
            )
          )
        ),
        column(1)
      ),
      tags$br(),

      ###### plot 2-4 settings ######
      tags$br(),

      fluidRow(
        column(1),
        column(
          10,

          div(
            bs4Dash::box(
              #### user input, plot configuration

              fluidRow(
                column(
                  3,
                  # cohort
                  selectInput(
                    inputId = "AnalysisPagePrognosticModelModule4PlotsCohortSelection",
                    label = "Cohort",
                    width = "100%",
                    choices = ""
                  ),
                ),
                column(
                  4,
                  # plot color
                  selectInput(
                    inputId = "AnalysisPagePrognosticModelModule4PlotsKMPlotColorSelection",
                    label = "Color palette",
                    choices = c("npg","nejm","jco","d3","lancet","jama"),
                    width = "100%"
                  ),
                ),
                column(
                  5,
                  # alpha
                  sliderInput(
                    inputId = "AnalysisPagePrognosticModelModule4PlotsKMPlotTransparencySelection",
                    label = "Color transparency",
                    min = 0.2,
                    max = 1,
                    value = 0.8,
                    step = 0.1,
                    ticks = TRUE,
                    animate = TRUE,
                    width = "100%"
                  ),
                )
              ),

              fluidRow(
                column(
                  3,
                  # time point for dca
                  selectInput(
                    inputId = "AnalysisPagePrognosticModelModule4PlotsDCAPlotTimePoint",
                    label = "Timepoint for DCA",
                    choices = c(),
                    width = "100%"
                  )
                ),
                column(
                  4,
                  # km method
                  selectInput(
                    inputId = "AnalysisPagePrognosticModelModule4PlotsKMPlotMethodSelection",
                    label = "Group cutoff method for KM analysis",
                    width = "100%",
                    choices = c('median','mean','quantile','optimal','custom')
                  ),
                ),
                column(
                  5,
                  # cutoff
                  div(
                    sliderInput(
                      inputId = "AnalysisPagePrognosticModelModule4PlotsCutoffSelection",
                      label = "Cutoff value",
                      width = "100%",
                      min = 0.1,
                      max = 0.9,
                      value = 0.3,
                      step = 0.05,
                      ticks = T,
                      animate = T
                    ),

                    id = "AnalysisPagePrognosticModelModule4PlotsCutoffSelectionDiv",
                    class = "hidden-element"
                  ),
                )
              ),

              ###### model evlation btn ######
              # evaluation model btn
              tags$br(),

              fluidRow(
                column(
                  3
                ),
                column(
                  6,
                  # shinyWidgets::actionBttn(
                  #   inputId = "AnalysisPagePrognosticModelModuleModelEvaluationSubmitBtn",
                  #   size = "md",
                  #   label = "Model Evaluation",
                  #   icon = icon("shield"),
                  #   style = "unite",
                  #   color = "danger",
                  #   block = TRUE
                  # )
                  shiny::actionButton(
                    inputId = "AnalysisPagePrognosticModelModuleModelEvaluationSubmitBtn",
                    label = "Everything is ready! Start Model Evaluation",
                    icon = icon("shield"),
                    width = "100%",
                    style = "background-color:#bc8a5f;color:white;font-size:120%;border-radius:20px;"
                  )
                ),
                column(
                  3
                )
              ),

              tags$br(),

              id = "PMPageEvaluationSettingCard",
              title = "7. Evaluation settings for model evaluation",
              icon = icon("cogs"),
              collapsible = FALSE,
              width = 12,
              closable = FALSE,
              maximizable = FALSE
            ),

            # style = "box-shadow:2px 2px 5px 2px #ccc;border-color: #bc8a5f; border-width: 2px;"
          )
        ),
        column(1)
      ),



      id = "AnalysisPagePrognosticModelModuleRiskScoreDiv",
      class = "hidden-element"
    ),


    #### model evaluation plot 2-4 ####

    ## plots section
    div(
      id = "AnalysisPagePrognosticModelModule4PlotsDiv",

      ## space
      tags$hr(),
      tags$br(),

      # title of plots
      fluidRow(
        column(
          3,
          div(
            tags$p("Kaplan-Meier Analysis",style="color:#ffffff;"),
            style = "text-align:center;background-color:#4ea653;"
          )
        ),
        column(
          3,
          div(
            tags$p("Time-dependent ROC",style="color:#ffffff;"),
            style = "text-align:center;background-color:#0aa9a2;"
          )
        ),
        column(
          3,
          div(
            tags$p("Calibration Curve",style="color:#ffffff;"),
            style = "text-align:center;background-color:#ff6b6b;"
          )
        ),
        column(
          3,
          div(
            tags$p("Decision Curve Analysis",style="color:#ffffff;"),
            style = "text-align:center;background-color:#fd7e14;"
          )
        )
      ),

      # four plots
      fluidRow(
        column(
          3,
          div(
            plotOutput(
              outputId = "AnalysisPagePrognosticModelModuleKMPlot",
              height = "340px",
              width = "340px"
            ),
            style = "overflow-x:auto;"
          )

        ),
        column(
          3,
          div(
            plotOutput(
              outputId = "AnalysisPagePrognosticModelModuleTimeROCPlot",
              height = "340px",
              width = "340px"
            ),
            style = "overflow-x:auto;"
          )

        ),
        column(
          3,
          div(
            plotOutput(
              outputId = "AnalysisPagePrognosticModelModuleCalCurvePlot",
              height = "340px",
              width = "340px"
            ),
            style = "overflow-x:auto;"
          )

        ),
        column(
          3,
          div(
            plotOutput(
              outputId = "AnalysisPagePrognosticModelModuleDCAPlot",
              height = "340px",
              width = "340px"
            ),
            style = "overflow-x:auto;"
          )
        )
      ),

      # button group for saving plots
      # fluidRow(
      #   column(
      #     3,
      #     shinyWidgets::dropdown(
      #       # content of the dropdown
      #       tags$h6("Save the plots with different formats"),
      #       tags$br(),
      #
      #       # save pdf
      #       # shinyWidgets::actionBttn(
      #       #   inputId = "AnalysisPagePrognosticModelModuleSavePdfBtn",
      #       #   icon = icon("file-pdf-o"),
      #       #   label = "PDF",
      #       #   color = "danger",
      #       #   style = "pill",
      #       #   size = "sm"
      #       # ),
      #
      #       shiny::downloadButton(
      #         outputId = "AnalysisPagePrognosticModelModuleSavePdfDownloadBtn",
      #         label = "PDF",
      #         class = "btn-danger",
      #         icon = icon("file-pdf-o")
      #       ),
      #
      #       # save tiff
      #       # shinyWidgets::actionBttn(
      #       #   inputId = "AnalysisPagePrognosticModelModuleSaveImgBtn",
      #       #   icon = icon("file-image-o"),
      #       #   label = "TIF",
      #       #   color = "danger",
      #       #   style = "pill",
      #       #   size = "sm"
      #       # ),
      #
      #       shiny::downloadButton(
      #         outputId = "AnalysisPagePrognosticModelModuleSaveTiffDownloadBtn",
      #         label = "Tiff",
      #         class = "btn-danger",
      #         icon = icon("file-image-o")
      #       ),
      #
      #
      #       ## settings of the dropdown button
      #       icon = icon("save"),
      #       label = "Save plots",
      #       size = "sm",
      #       style = "unite",
      #       status = "danger",
      #       animate = animateOptions(
      #         enter = animations$fading_entrances$fadeInDown,
      #         exit = animations$fading_exits$fadeOutUp
      #       )
      #     )
      #   ),
      #   column(
      #     9
      #   )
      # ),

      tags$br(),
      span(
        shiny::downloadButton(
          outputId = "AnalysisPagePrognosticModelModulePlot234DownloadPdfBtn",
          label = "Download PDF",
          icon = icon("file-pdf-o"),
          style = "color:#FF6B6B;background-color:#ffffff;",
        ),
        shiny::downloadButton(
          outputId = "AnalysisPagePrognosticModelModulePlot234DownloadPngBtn",
          label = "Download PNG",
          icon = icon("file-image-o"),
          style = "color:#FF6B6B;background-color:#ffffff;",
        )
      ),

      class = "hidden-element"
    ),


    ## space
    tags$br(),
    tags$hr(),
    tags$br(),

    ## single gene analysis tab pages
    ui.page_module_prognostic_model_down_section()
  )
}
