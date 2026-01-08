####################################################
# Project: AOP Network Box
# Author: Minju Na
# File: 6_AOP.R
####################################################

library(shiny)
library(dplyr)
library(tidyr)
library(DT)
library(visNetwork)
library(bslib)

# --- UI Module ---
AOPUI <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      card_header("Network Visualization & Filtering"),
      card_body(
        layout_columns(
          col_widths = c(8, 4),
          selectInput(ns("select_cid_vis"), "Select Compound (CID):", choices = NULL),
          actionButton(ns("fullscreen_btn"), "Fullscreen Mode", icon = icon("expand"), class = "btn-outline-primary")
        ),
        visNetworkOutput(ns("network"), height = "600px") %>% withSpinner(color = "#ef4444")
      )
    ),
    navset_card_underline(
      nav_panel("AOP-GO Linkage",
        DTOutput(ns("AOP_GO_Table")),
        downloadButton(ns("download_go"), "Download AOP-GO Data", class = "mt-2")
      ),
      nav_panel("AOP-Protein Linkage",
        DTOutput(ns("AOP_Protein_Table")),
        downloadButton(ns("download_prot"), "Download AOP-Protein Data", class = "mt-2")
      )
    )
  )
}

# --- Server Module ---
AOPServer <- function(id, global_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 1. Reactive State for Linkage Data
    aop_state <- reactiveVal(list(go = NULL, prot = NULL))

    # 2. Perform Linkage Analysis
    observe({
      # Load internal AOP Wiki mappings (provided in /Data folder)
      req(file.exists("./Data/aop_go_df.rda"), file.exists("./Data/aop_prot_df.rda"))
      load("./Data/aop_go_df.rda")
      load("./Data/aop_prot_df.rda")
      load("./Data/df_all_aop_fin.rda")
      
      # Get results from previous modules
      go_results <- global_data()$go_results$table
      prot_results <- load_shared_data("table_protein")
      req(go_results, prot_results)

      # Link GO terms to AOPs
      go_linkage <- aop_go_df %>%
        inner_join(go_results, by = c("Process_Ontology_ID" = "ID")) %>%
        inner_join(df_all_aop_fin, by = c("AOP" = "AOP"))

      # Link Proteins to AOPs
      prot_linkage <- aop_prot_df %>%
        inner_join(prot_results, by = "Ensemble") %>%
        inner_join(df_all_aop_fin, by = c("AOP" = "AOP"))

      aop_state(list(go = go_linkage, prot = prot_linkage))
      
      # Update Filter Choices
      updateSelectInput(session, "select_cid_vis", choices = unique(prot_results$CID))
    })

    # 3. Network Visualization
    output$network <- renderVisNetwork({
      data <- aop_state()
      req(data$prot)
      
      # Filter nodes based on selected CID
      prot_sub <- data$prot %>% filter(CID == input$select_cid_vis)
      
      # Create Nodes: Compound -> Protein -> AOP
      nodes <- data.frame(
        id = c(input$select_cid_vis, unique(prot_sub$Ensemble), unique(prot_sub$AOP)),
        label = c(paste("CID:", input$select_cid_vis), unique(prot_sub$Ensemble), unique(prot_sub$AOP)),
        group = c("Compound", rep("Protein", length(unique(prot_sub$Ensemble))), rep("AOP", length(unique(prot_sub$AOP)))),
        color = c("#3b82f6", "#10b981", "#ef4444")
      )
      
      # Create Edges
      edges <- bind_rows(
        prot_sub %>% select(from = CID, to = Ensemble) %>% distinct(),
        prot_sub %>% select(from = Ensemble, to = AOP) %>% distinct()
      )

      visNetwork(nodes, edges) %>%
        visGroups(groupname = "AOP", shape = "diamond") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visPhysics(stabilization = TRUE)
    })

    # 4. Table Rendering (AOP-GO)
    output$AOP_GO_Table <- renderDT({
      df <- aop_state()$go
      req(df)
      df_display <- df %>%
        mutate(
          AOP = sprintf("<a href='https://aopwiki.org/aops/%s' target='_blank'>%s</a>", AOP, AOP),
          GO_Term = sprintf("<a href='https://amigo.geneontology.org/amigo/term/%s' target='_blank'>%s</a>", Process_Ontology_ID, Process_Ontology_ID)
        ) %>%
        select(CID, AOP, Event_Type, Description, GO_Term)
      
      datatable(df_display, escape = FALSE, options = list(scrollX = TRUE))
    })

    # 5. Table Rendering (AOP-Protein)
    output$AOP_Protein_Table <- renderDT({
      df <- aop_state()$prot
      req(df)
      df_display <- df %>%
        mutate(
          AOP = sprintf("<a href='https://aopwiki.org/aops/%s' target='_blank'>%s</a>", AOP, AOP),
          Ensemble = sprintf("<a href='https://string-db.org/network/9606.%s' target='_blank'>%s</a>", Ensemble, Ensemble)
        ) %>%
        select(CID, AOP, Event_Type, Ensemble, Score)
      
      datatable(df_display, escape = FALSE, options = list(scrollX = TRUE))
    })

    # 6. Download Handlers
    output$download_go <- downloadHandler(
      filename = "AOP_GO_Linkage.csv",
      content = function(file) write.csv(aop_state()$go, file)
    )
    
    output$download_prot <- downloadHandler(
      filename = "AOP_Protein_Linkage.csv",
      content = function(file) write.csv(aop_state()$prot, file)
    )
  })
}
