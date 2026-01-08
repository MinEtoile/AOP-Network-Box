####################################################
# Project: AOP Network Box
# Author: Minju Na
# File: 5_DO.R
####################################################

library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(shinycssloaders)

# --- Helper: Format Gene IDs with Links ---
format_do_genes <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "") return("")
  # Remove potential HTML tags and split
  genes <- unlist(strsplit(gsub("<.*?>", "", gene_str), "[/|,]"))
  genes <- trimws(genes[genes != ""])
  
  if (length(genes) == 0) return("")
  
  # Generate clickable links for GeneCards
  links <- sapply(genes, function(g) {
    sprintf("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank' class='btn btn-xs btn-outline-success'>%s</a>", g, g)
  })
  
  # Truncate for UI performance; full list via toggle
  if (length(links) > 10) {
    truncated <- paste(links[1:10], collapse = " ")
    return(paste0(truncated, sprintf(" <span class='show-full-content' data-full-content='%s' style='color:#ef4444;cursor:pointer;'>...</span>", paste(genes, collapse = ", "))))
  }
  return(paste(links, collapse = " "))
}

# --- UI Module ---
doUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    navset_card_underline(
      nav_panel("Disease Table",
        div(style = "background: #fef2f2; padding: 20px; border-radius: 12px; margin-bottom: 20px; border: 1px solid #fee2e2;",
          selectInput(ns("select_cid_do"), "Filter by Compound (CID):", 
                      choices = NULL, multiple = TRUE, width = "100%")
        ),
        DTOutput(ns("DOTable")) %>% withSpinner(color = "#ef4444"),
        downloadButton(ns("downloadData_do"), "Download DO Results", class = "mt-3")
      ),
      nav_panel("Enrichment Plot",
        div(style = "background: #fef2f2; padding: 20px; border-radius: 12px; margin-bottom: 20px; border: 1px solid #fee2e2;",
          layout_columns(
            col_widths = c(8, 4),
            selectInput(ns("select_smi_plot"), "Select SMILES:", choices = NULL),
            sliderInput(ns("num_cat"), "Top Categories:", min = 5, max = 30, value = 15)
          ),
          actionButton(ns("update_plot"), "Update Visualization", class = "btn-danger w-100")
        ),
        plotOutput(ns("DOPlot"), height = "600px") %>% withSpinner(color = "#ef4444"),
        downloadButton(ns("downloadPlot_do"), "Download Plot (PNG)", class = "mt-3")
      )
    )
  )
}

# --- Server Module ---
doServer <- function(id, global_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Internal state for processed results
    display_state <- reactiveVal(list(table = NULL, plot_objects = NULL))

    # Sync with global enrichment results
    observeEvent(global_data()$do_results, {
      results <- global_data()$do_results
      req(results$table)
      
      # Update UI Filters based on SMILES-CID mapping
      mapping <- load_shared_data("smiles_cid_mapping")
      if (!is.null(mapping)) {
        choices <- setNames(mapping$CID, paste0(mapping$SMILES, " (", mapping$CID, ")"))
        updateSelectInput(session, "select_cid_do", choices = choices)
        updateSelectInput(session, "select_smi_plot", choices = unique(mapping$SMILES))
      }
      
      display_state(list(table = results$table, plot_objects = results$plot_objects))
    })

    # Render Main Data Table
    output$DOTable <- renderDT({
      df <- display_state()$table
      req(df)
      if ("Message" %in% colnames(df)) return(datatable(df))
      
      # Filter by user selection
      if (!is.null(input$select_cid_do)) {
        df <- df %>% filter(gsub("<.*?>", "", CID) %in% input$select_cid_do)
      }

      # Add DB links and format p-values
      df_display <- df %>%
        mutate(
          CID = sprintf("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/%s' target='_blank' class='btn btn-xs btn-outline-primary'>%s</a>", CID, CID),
          DO_ID = if("ID" %in% colnames(df)) ID else "",
          ID = sprintf("<a href='https://disease-ontology.org/term/%s' target='_blank'>%s</a>", DO_ID, DO_ID),
          geneID = sapply(geneID, format_do_genes),
          pvalue = formatC(pvalue, format = "e", digits = 2)
        )

      datatable(df_display, escape = FALSE, options = list(scrollX = TRUE, pageLength = 10))
    })

    # Render Dotplot Visualization
    output$DOPlot <- renderPlot({
      req(input$update_plot)
      results <- display_state()
      req(results$plot_objects)
      
      # Extract analysis object for selected SMILES
      edo_obj <- results$plot_objects[[input$select_smi_plot]]
      if (is.null(edo_obj)) edo_obj <- results$plot_objects[[1]] # Fallback
      
      req(edo_obj)
      dotplot(edo_obj, showCategory = input$num_cat) +
        ggtitle("Disease Ontology Enrichment") +
        theme_minimal() +
        scale_color_gradient(low = "#ef4444", high = "#fca5a5")
    })

    # Export Handlers
    output$downloadData_do <- downloadHandler(
      filename = paste0("DO_Results_", Sys.Date(), ".csv"),
      content = function(file) write.csv(display_state()$table, file, row.names = FALSE)
    )
  })
}

# --- Core Analysis Function ---
perform_do_analysis <- function(entrez_ids, pval_cutoff = 0.05) {
  # Standard DOSE enrichment call
  edo <- enrichDO(gene = entrez_ids, ont = "DO", 
                  pvalueCutoff = pval_cutoff, 
                  pAdjustMethod = "BH", 
                  readable = TRUE)
  return(edo)
}
