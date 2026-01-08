####################################################
# Project: AOP Network Box
# Author: Minju Na
# File: 3_GO.R
####################################################

library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# --- UI Module ---
goUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    navset_card_underline(
      nav_panel("Enrichment Table",
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          selectInput(ns("select_cid_go"), "Filter by Compound (CID):", 
                      choices = NULL, multiple = TRUE, width = "100%")
        ),
        DTOutput(ns("GOTable")) %>% withSpinner(color = "#8b5cf6"),
        downloadButton(ns("downloadData_go"), "Download GO Results (CSV)", class = "mt-3")
      ),
      nav_panel("Dot Plot",
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          layout_columns(
            col_widths = c(4, 4, 4),
            selectInput(ns("select_smi_plot"), "SMILES:", choices = NULL),
            selectInput(ns("select_ont"), "Ontology:", choices = c("BP", "CC", "MF")),
            sliderInput(ns("num_cat"), "Categories:", min = 5, max = 30, value = 15)
          ),
          actionButton(ns("update_plot"), "Update Visualization", class = "btn-primary w-100")
        ),
        plotOutput(ns("GOPlot"), height = "600px") %>% withSpinner(color = "#8b5cf6"),
        downloadButton(ns("downloadPlot_go"), "Download Plot (PNG)", class = "mt-3")
      )
    )
  )
}

# --- Server Module ---
goServer <- function(id, global_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 1. Internal state for processed results
    display_state <- reactiveVal(list(table = NULL, plot_data = NULL))

    # 2. Sync with global analysis results
    observeEvent(global_data()$go_results, {
      results <- global_data()$go_results
      req(results$table)
      
      # Populate selection filters
      mapping <- load_shared_data("smiles_cid_mapping")
      if (!is.null(mapping)) {
        choices <- setNames(mapping$CID, paste0(mapping$SMILES, " (", mapping$CID, ")"))
        updateSelectInput(session, "select_cid_go", choices = choices)
        updateSelectInput(session, "select_smi_plot", choices = unique(mapping$SMILES))
      }
      
      display_state(list(table = results$table, plot_objects = results$plot_objects))
    })

    # 3. Interactive Table Rendering
    output$GOTable <- renderDT({
      df <- display_state()$table
      req(df)
      if ("Message" %in% colnames(df)) return(datatable(df))
      
      # Filter by CID if selected
      if (!is.null(input$select_cid_go)) {
        df <- df %>% filter(gsub("<.*?>", "", CID) %in% input$select_cid_go)
      }

      # Format links and scientific notation
      df_display <- df %>%
        mutate(
          CID = sprintf("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/%s' target='_blank' class='btn btn-xs btn-outline-primary'>%s</a>", CID, CID),
          ID = sprintf("<a href='https://amigo.geneontology.org/amigo/term/%s' target='_blank'>%s</a>", ID, ID),
          pvalue = formatC(pvalue, format = "e", digits = 2)
        )

      datatable(df_display, escape = FALSE, options = list(scrollX = TRUE, pageLength = 10))
    })

    # 4. Dot Plot Visualization
    output$GOPlot <- renderPlot({
      req(input$update_plot)
      results <- display_state()
      req(results$plot_objects)
      
      # Extract specific ontology object (BP, CC, or MF)
      ego_obj <- results$plot_objects[[input$select_ont]]
      req(ego_obj)

      dotplot(ego_obj, showCategory = input$num_cat) +
        ggtitle(paste("GO Enrichment:", input$select_ont)) +
        theme_minimal() +
        scale_color_gradient(low = "#dc2626", high = "#3b82f6")
    })

    # 5. Export Handlers
    output$downloadData_go <- downloadHandler(
      filename = "GO_Enrichment_Results.csv",
      content = function(file) write.csv(display_state()$table, file, row.names = FALSE)
    )
  })
}

# --- Core Analysis Helper ---
perform_go_analysis <- function(entrez_ids, pval_cutoff = 0.05) {
  # Runs standard enrichment for all three sub-ontologies
  ontologies <- c("BP", "CC", "MF")
  results <- lapply(ontologies, function(ont) {
    enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = ont, 
             pAdjustMethod = "BH", pvalueCutoff = pval_cutoff, readable = TRUE)
  })
  names(results) <- ontologies
  return(results)
}
