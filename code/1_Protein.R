####################################################
# Project: AOP Network Box
# Author: Minju Na
# File: 1_Protein.R
####################################################

library(shiny)
library(dplyr)
library(tidyr)
library(DT)
library(httr)
library(jsonlite)
library(parallel)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- UI Module ---
proteinUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    navset_card_underline(
      nav_panel(
        "Analysis Table",
        # Filters Section
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          fluidRow(
            column(6, selectInput(ns("select_cid"), "Filter by Compound (CID):", choices = NULL, multiple = TRUE, width = "100%")),
            column(6, numericInput(ns("score_filter"), "Min Score (0.0 - 1.0):", value = 0.15, min = 0, max = 1, step = 0.01))
          )
        ),
        # Table Display
        div(style = "background: white; padding: 15px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.05);",
          DTOutput(ns("ProteinTable")) %>% withSpinner()
        ),
        downloadButton(ns("downloadData"), "Download Results (CSV)", class = "btn-block mt-3")
      )
    )
  )
}

# --- Server Module ---
proteinServer <- function(input, output, session, selected_smiles, selected_score, selected_limit) {
  ns <- session$ns

  # Update filter choices based on available data
  observe({
    data <- protein_data()
    req(nrow(data) > 0)
    mapping <- load_shared_data("smiles_cid_mapping")
    if (!is.null(mapping)) {
      choices <- setNames(mapping$CID, paste0(mapping$SMILES, " (", mapping$CID, ")"))
      updateSelectInput(session, "select_cid", choices = choices)
    }
  })

  # Helper function for parallel STRING API fetching
  fetch_string_data <- function(smi, threshold) {
    tryCatch({
      library(webchem)
      library(jsonlite)
      library(RCurl)
      
      cid <- get_cid(smi, from = "smiles")$cid
      # API returns Score in 0-1 range
      url <- paste0("https://cwtung.nhri.edu.tw/chemdis/api/chemdis/CID", 
                    sprintf("%09d", as.integer(cid)), "/5/", threshold, "/protein/1")
      
      res <- getURL(url, .opts = list(ssl.verifypeer = FALSE))
      json_data <- fromJSON(res)
      if (length(json_data) == 0) return(NULL)
      
      df <- as.data.frame(json_data)
      df$SMILES <- smi
      df$CID <- cid
      # Standardize Score to 0-1000 range
      df$Score <- round(as.numeric(df$Score) * 1000)
      return(df)
    }, error = function(e) NULL)
  }

  # Core Data Processing Reactive
  protein_data <- reactive({
    # 1. Check Cache
    cached <- load_shared_data("table_protein")
    if (!is.null(cached)) return(cached)

    # 2. Parallel API Execution
    withProgress(message = "Fetching Protein Interactions...", value = 0, {
      cl <- makeCluster(detectCores() - 1)
      clusterExport(cl, c("fetch_string_data", "selected_score"))
      
      api_results <- parLapply(cl, selected_smiles, function(smi) {
        fetch_string_data(smi, selected_score)
      })
      stopCluster(cl)
      
      combined_api <- bind_rows(api_results) %>% distinct()
      
      # 3. Merge with Python AI Predictions
      if (file.exists("./Data/aop_pred_prot_df.rda")) {
        load("./Data/aop_pred_prot_df.rda")
        # Standardize score and merge by Ensemble ID
        final_df <- full_join(combined_api, aop_pred_prot_df, by = "Ensemble")
      } else {
        final_df <- combined_api
      }

      # 4. AnnotationDbi Gene Mapping (UniProt/Symbol)
      try({
        ensp_ids <- gsub("^9606\\.", "", na.omit(final_df$Ensemble))
        mapping <- select(org.Hs.eg.db, keys = ensp_ids, 
                          columns = c("SYMBOL", "UNIPROT", "GENENAME"), 
                          keytype = "ENSEMBLPROT")
        
        final_df <- final_df %>%
          left_join(mapping, by = c("Ensemble" = "ENSEMBLPROT")) %>%
          mutate(Gene_symbol = coalesce(SYMBOL, Gene_symbol))
      }, silent = TRUE)

      # 5. Final Formatting & Saving
      final_df <- final_df %>%
        filter(!is.na(Ensemble)) %>%
        mutate(Score = as.numeric(Score))
      
      save_shared_data("table_protein", final_df)
      save(final_df, file = "./Data/table_protein.rda")
      
      # Signal global completion for UI progress bars
      shinyjs::runjs("if(window.updateProgress) { window.moduleCompletion.protein = true; window.updateProgress(); }")
      
      return(final_df)
    })
  })

  # Render Interactive DT with External Links
  output$ProteinTable <- renderDT({
    df <- protein_data()
    req(nrow(df) > 0)
    
    # Filter by user input
    if (!is.null(input$select_cid)) df <- df %>% filter(CID %in% input$select_cid)
    df <- df %>% filter(Score >= (input$score_filter * 1000))

    df_display <- df %>%
      mutate(
        CID = sprintf("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/%s' target='_blank' class='btn btn-xs btn-outline-primary'>%s</a>", CID, CID),
        Ensemble = sprintf("<a href='https://string-db.org/network/9606.%s' target='_blank'>%s</a>", Ensemble, Ensemble),
        Gene_symbol = sprintf("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank'>%s</a>", Gene_symbol, Gene_symbol),
        Score = round(Score / 1000, 3)
      ) %>%
      select(CID, Ensemble, Gene_symbol, Gene_name = GENENAME, Score)

    datatable(df_display, escape = FALSE, selection = 'none', 
              options = list(pageLength = 10, scrollX = TRUE))
  })

  # Download Handler
  output$downloadData <- downloadHandler(
    filename = "Protein_Binding_Results.csv",
    content = function(file) {
      write.csv(protein_data(), file, row.names = FALSE)
    }
  )
}
