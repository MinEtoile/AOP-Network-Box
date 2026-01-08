####################################################
# Project: AOP Network Box
# Author: Minju Na
# File: 2_STRING.R
####################################################

library(shiny)
library(dplyr)
library(DT)
library(httr)
library(png)
library(grid)
library(parallel)
library(shinyjs)
library(shinycssloaders)

# --- UI Module ---
STRING_UI <- function(id) {
  ns <- NS(id)
  
  tagList(
    navset_card_underline(
      nav_panel(
        "PPI Table",
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          selectInput(ns("select_cid"), "Filter by Compound (CID):", 
                      choices = NULL, multiple = TRUE, width = "100%")
        ),
        DTOutput(ns("PPI_Table")) %>% withSpinner(color = "#1e3a8a"),
        downloadButton(ns("downloadData_string"), "Download PPI Data", class = "mt-3")
      ),
      nav_panel(
        "PPI Plot",
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          fluidRow(
            column(6, selectInput(ns("select_cid_plot"), "Select CID for Plot:", choices = NULL)),
            column(6, numericInput(ns("ppi_max_interactions"), "Max Interactions:", value = 20, min = 1, max = 500))
          )
        ),
        plotOutput(ns("PPI_Plot"), height = "600px") %>% withSpinner(color = "#1e3a8a"),
        downloadButton(ns("downloadPlot_string"), "Download Plot (SVG)", class = "mt-3")
      )
    )
  )
}

# --- Server Module ---
STRING_Server <- function(input, output, session, selected_smiles, selected_score, selected_limit) {
  ns <- session$ns

  # 1. Initialize CID selection from shared mapping
  observe({
    mapping <- load_shared_data("smiles_cid_mapping")
    req(mapping)
    valid_map <- mapping[mapping$CID != "N/A", ]
    choices <- setNames(valid_map$CID, paste0(valid_map$SMILES, " (", valid_map$CID, ")"))
    updateSelectInput(session, "select_cid", choices = choices)
    updateSelectInput(session, "select_cid_plot", choices = choices)
  })

  # 2. Fetch PPI Results from STRING API
  ppi_results <- reactive({
    # Check cache
    cached <- load_shared_data("ppi_results")
    if (!is.null(cached)) return(cached)

    protein_data <- load_shared_data("table_protein")
    req(protein_data)

    # Parallel processing for multiple SMILES
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("selected_score", "selected_limit"))
    
    results_list <- parLapply(cl, unique(protein_data$SMILES), function(smi) {
      tryCatch({
        library(httr)
        # Construct API request
        url <- "https://version-12-0.string-db.org/api/tsv-no-header/interaction_partners"
        # Filter top proteins for the query
        query_prots <- protein_data[protein_data$SMILES == smi, ] %>% 
          arrange(desc(Score)) %>% head(5) %>% pull(Ensemble)
        
        res <- POST(url, body = list(
          identifiers = paste0("9606.", query_prots, collapse = "\n"),
          species = 9606,
          required_score = as.numeric(selected_score) * 1000,
          limit = 100
        ), encode = "form")

        if (status_code(res) == 200) {
          df <- read.table(text = content(res, "text"), sep = "\t")
          return(df[, c(1, 2, 3, 4, 6)]) # Standardize columns
        }
        return(NULL)
      }, error = function(e) NULL)
    })
    stopCluster(cl)

    final_df <- bind_rows(results_list)
    save_shared_data("ppi_results", final_df)
    
    # Signal completion to master UI
    shinyjs::runjs("if(window.updateProgress) { window.moduleCompletion.string = true; window.updateProgress(); }")
    return(final_df)
  })

  # 3. Render Table
  output$PPI_Table <- renderDT({
    df <- ppi_results()
    req(nrow(df) > 0)
    if (!is.null(input$select_cid)) df <- df[df$CID %in% input$select_cid, ]
    
    datatable(df, escape = FALSE, options = list(scrollX = TRUE, pageLength = 10))
  })

  # 4. Render Network Image
  output$PPI_Plot <- renderImage({
    req(input$select_cid_plot)
    outfile <- tempfile(fileext = '.png')
    
    # Fetch image from STRING static generator
    url <- "https://string-db.org/api/image/network"
    res <- POST(url, body = list(
      identifiers = input$select_cid_plot,
      species = 9606,
      required_score = as.numeric(selected_score) * 1000,
      add_white_nodes = input$ppi_max_interactions
    ), encode = "form")
    
    writeBin(content(res, "raw"), outfile)
    
    list(src = outfile, contentType = 'image/png', width = "100%", alt = "STRING PPI Network")
  }, deleteFile = TRUE)

  # 5. Export Handlers
  output$downloadData_string <- downloadHandler(
    filename = "PPI_Interactions.csv",
    content = function(file) write.csv(ppi_results(), file, row.names = FALSE)
  )
}
