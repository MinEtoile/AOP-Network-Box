###################################################
# Author: Minju Na
# Project: AOP Network Box
# File: ui.R
###################################################

## Load libraries ##
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(bslib)
library(DT)
library(visNetwork)

## Load Source Files ##
source("1_Protein.R", local = TRUE)
source("2_STRING.R", local = TRUE)
source("3_GO.R", local = TRUE)
source("4_Pathway.R", local = TRUE)
source("5_DO.R", local = TRUE)
source("6_AOP.R", local = TRUE)
source("tutorial.R", local = TRUE)
source("contact.R", local = TRUE)

## Set UI ##
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = NULL),
  title = "AOP Network Box",
  div(class = "main-header-banner",
    h1("AOP Network Box", id = "app_title", title = "Click to reset"),
      p("Adverse Outcome Pathway Analysis Platform")
  ),

  mainPanel(
    width = 12,
    tabsetPanel(
      id = "main_tabs",
      
      tabPanel("Prediction",
        div(class = "content-wrapper",
            div(class = "prediction-container",
                h3("Input Parameters"),
                hr(class="hr-styled"),
                fluidRow(
                  column(6,
          card(
                      card_header(h4("Chemical & Interaction")),
                      textAreaInput("input_smi", "Paste the SMILES", width = "100%", rows = 4,
                                    placeholder = "e.g., (fluoxetine) CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F"),
                      helpText("e.g., (fluoxetine) CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F"),
                      helpText("Each SMILES requires 1-2 minutes to analysis. Enter 1-3 SMILES for optimal performance."),
                      div(style = "margin-bottom: 1.5rem;",
                          h5("Select the minimum required interaction score", style = "font-size: 1.25rem; font-weight: 600; color: var(--text-primary); margin-bottom: 1rem;"),
                          div(style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 1rem;",
                              div(class = "score-card", style = "display: flex; align-items: center; gap: 0.8rem; padding: 1.2rem 1.5rem; border: 2px solid var(--border-color); border-radius: 12px; background: linear-gradient(135deg, #f8fafc 0%, #ffffff 100%); transition: all 0.3s ease; cursor: pointer; position: relative;",
                                  tags$input(type = "radio", name = "selectScore", value = "150", id = "score_low", class = "score-radio", checked = TRUE),
                                  tags$label("Low (0.15)", `for` = "score_low", style = "margin: 0; font-size: 1.25rem; font-weight: 500; color: var(--text-primary); cursor: pointer; flex: 1;")
                              ),
                              div(class = "score-card", style = "display: flex; align-items: center; gap: 0.8rem; padding: 1.2rem 1.5rem; border: 2px solid var(--border-color); border-radius: 12px; background: linear-gradient(135deg, #f8fafc 0%, #ffffff 100%); transition: all 0.3s ease; cursor: pointer; position: relative;",
                                  tags$input(type = "radio", name = "selectScore", value = "400", id = "score_medium", class = "score-radio"),
                                  tags$label("Medium (0.4)", `for` = "score_medium", style = "margin: 0; font-size: 1.25rem; font-weight: 500; color: var(--text-primary); cursor: pointer; flex: 1;")
                              ),
                              div(class = "score-card", style = "display: flex; align-items: center; gap: 0.8rem; padding: 1.2rem 1.5rem; border: 2px solid var(--border-color); border-radius: 12px; background: linear-gradient(135deg, #f8fafc 0%, #ffffff 100%); transition: all 0.3s ease; cursor: pointer; position: relative;",
                                  tags$input(type = "radio", name = "selectScore", value = "700", id = "score_high", class = "score-radio"),
                                  tags$label("High (0.7)", `for` = "score_high", style = "margin: 0; font-size: 1.25rem; font-weight: 500; color: var(--text-primary); cursor: pointer; flex: 1;")
                              ),
                              div(class = "score-card", style = "display: flex; align-items: center; gap: 0.8rem; padding: 1.2rem 1.5rem; border: 2px solid var(--border-color); border-radius: 12px; background: linear-gradient(135deg, #f8fafc 0%, #ffffff 100%); transition: all 0.3s ease; cursor: pointer; position: relative;",
                                  tags$input(type = "radio", name = "selectScore", value = "900", id = "score_highest", class = "score-radio"),
                                  tags$label("Highest (0.9)", `for` = "score_highest", style = "margin: 0; font-size: 1.25rem; font-weight: 500; color: var(--text-primary); cursor: pointer; flex: 1;")
                              )
                          )
                      ),
                      numericInput("selectInt", "Select the maximum number of interactors", width = "100%", value = 100, min = 10, max = 100000)
                    )
                  ),
                  column(6,
          card(
                      card_header(h4("Enrichment Analysis")),
                      pickerInput("selectOntology", "Select the Ontology", width = "100%",
                        choices = c("Molecular Function" = "MF", "Biological Process" = "BP", "Cellular Component" = "CC"),
                                  options = list('action-box' = TRUE),
                                  multiple = TRUE),
                      pickerInput("selectpval_method", "Select the p-value adjust method", width = "100%",
                        choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                  selected = "BH", options = list('action-box' = TRUE)),
                      div(style = "margin-top: 1rem;",
                          h5("Select the p-value cutoff", style = "font-size: 1.25rem; font-weight: 600; color: var(--text-primary); margin-bottom: 1rem;"),
                          div(style = "display: flex; gap: 1.5rem; align-items: center;",
                              div(class = "pval-card", style = "display: flex; align-items: center; gap: 0.8rem; padding: 1rem 1.5rem; border: 2px solid var(--border-color); border-radius: 12px; background: #f8fafc; transition: all 0.3s ease; cursor: pointer;",
                                  tags$input(type = "radio", name = "selectpval", value = "0.05", id = "pval_005", class = "pval-radio", checked = TRUE),
                                  tags$label("0.05", `for` = "pval_005", style = "margin: 0; font-size: 1.25rem; font-weight: 500; color: var(--text-primary); cursor: pointer;")
                              ),
                              div(class = "pval-card", style = "display: flex; align-items: center; gap: 0.8rem; padding: 1rem 1.5rem; border: 2px solid var(--border-color); border-radius: 12px; background: #f8fafc; transition: all 0.3s ease; cursor: pointer;",
                                  tags$input(type = "radio", name = "selectpval", value = "0.01", id = "pval_001", class = "pval-radio"),
                                  tags$label("0.01", `for` = "pval_001", style = "margin: 0; font-size: 1.25rem; font-weight: 500; color: var(--text-primary); cursor: pointer;")
                              )
                          )
                      )
                    )
                  )
                ),
                div(style = "text-align: center; margin-top: 2rem;",
                    useShinyjs(),
                    actionButton("InputSubmit", "Run Analysis", icon = icon("rocket"), class = "btn-submit-custom")
                )
            )
        )
      ),
      
      tabPanel("Result",
        div(class="content-wrapper",
            div(style = "margin-bottom: 2rem;",
                h3("Analysis Results", style = "text-align: center; color: var(--text-primary); font-weight: 700; margin-bottom: 1rem;"),
                p("Explore detailed analysis results across different biological domains", style = "text-align: center; color: var(--text-secondary); font-size: 1.6rem; font-weight: 500;")
            ),
            # SMILES-CID Mapping Display (compact) - only show when mapping exists
            uiOutput("smiles_cid_mapping_ui"),
            navset_card_tab(
              id = "result_tabs",
              nav_panel(
                value = "protein_tab",
                title = div(style = "display: flex; align-items: center; gap: 0.5rem;",
                           icon("dna", style = "color: #3b82f6;"),
                           span("Protein", style = "font-weight: 600;")
                ),
                proteinUI("1_Protein")
              ),
              nav_panel(
                value = "ppi_tab", 
                title = div(style = "display: flex; align-items: center; gap: 0.5rem;",
                           icon("project-diagram", style = "color: #10b981;"),
                           span("PPI", style = "font-weight: 600;")
                ),
                STRING_UI("2_STRING")
              ),
              nav_panel(
                value = "go_tab",
                title = div(style = "display: flex; align-items: center; gap: 0.5rem;",
                           icon("sitemap", style = "color: #8b5cf6;"),
                           span("GO", style = "font-weight: 600;")
                ),
                goUI("3_GO")
              ),
              nav_panel(
                value = "pathway_tab",
                title = div(style = "display: flex; align-items: center; gap: 0.5rem;",
                           icon("route", style = "color: #f59e0b;"),
                           span("Pathway", style = "font-weight: 600;")
                ),
                pathwayUI("4_Pathway_mod")
              ),
              nav_panel(
                value = "do_tab",
                title = div(style = "display: flex; align-items: center; gap: 0.5rem;",
                           icon("heart-pulse", style = "color: #ef4444;"),
                           span("DO", style = "font-weight: 600;")
                ),
                doUI("5_DO")
              ),
              nav_panel(
                value = "aop_tab",
                title = div(style = "display: flex; align-items: center; gap: 0.5rem;",
                           icon("network-wired", style = "color: #06b6d4;"),
                           span("AOP", style = "font-weight: 600;")
                ),
                AOPUI("6_AOP")
              )
            )
        )
      ),

      tabPanel("Tutorial", div(class="content-wrapper", tutorialUI("tutorial"))),
      tabPanel("Contact", div(class="content-wrapper", contactUI("contact")))
    )
  ),
  
  div(class = "footer-container",
    p("AOP Network Box", style = "font-size: 1.35rem; font-weight: 600; color: var(--text-primary); margin-bottom: 0.5rem;"),
    p("© 2025 BCBL Lab, KAIST", style = "font-size: 1.25rem; font-weight: 500; color: var(--text-secondary);"),
    p("Licensed under Apache 2.0", style = "font-size: 1.25rem; font-weight: 500; color: var(--text-secondary);")
  )
)
