####################################################
# Project: AOP Network Box
# Author: Minju Na
# File: 4_Pathway_mod.R
####################################################

library(shiny)
library(dplyr)
library(DT)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(parallel)

# --- Helper: Format Gene IDs with Links ---
format_pathway_genes <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "") return("")
  genes <- unlist(strsplit(gsub("<.*?>", "", gene_str), "[/|,]"))
  genes <- trimws(genes[genes != ""])
  
  if (length(genes) == 0) return("")
  
  links <- sapply(genes, function(g) {
    sprintf("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank' class='btn btn-xs btn-outline-success'>%s</a>", g, g)
  })
  
  # Truncate for UI performance; full list available via tooltip/click
  if (length(links) > 10) {
    truncated <- paste(links[1:10], collapse = " ")
    return(paste0(truncated, sprintf(" <span class='show-full-content' data-full-content='%s' style='color:blue;cursor:pointer;'>...</span>", paste(genes, collapse = ", "))))
  }
  return(paste(links, collapse = " "))
}

# --- UI Module ---
pathwayUI <- function(id) {
  ns <- NS(id)
  tagList(
    navset_card_underline(
      nav_panel("KEGG",
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          selectInput(ns("select_cid_kegg"), "Filter by Compound (CID):", choices = NULL, multiple = TRUE, width = "100%")
        ),
        DTOutput(ns("KEGG_Table")) %>% withSpinner(color = "#1e3a8a"),
        downloadButton(ns("KEGG_Download"), "Download KEGG Results", class = "mt-3")
      ),
      nav_panel("Reactome",
        div(style = "background: #f8fafc; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
          selectInput(ns("select_cid_reactome"), "Filter by Compound (CID):", choices = NULL, multiple = TRUE, width = "100%")
        ),
        DTOutput(ns("Reactome_Table")) %>% withSpinner(color = "#e11d48"),
        downloadButton(ns("Reactome_Download"), "Download Reactome Results", class = "mt-3")
      )
    )
  )
}

# --- Server Module ---
pathwayServer <- function(input, output, session, selected_smiles, selected_score, selected_limit, pval) {
  ns <- session$ns
  pathway_data <- reactiveVal(NULL)

  # 1. Main Analysis Trigger
  observe({
    protein_df <- load_shared_data("table_protein")
    req(protein_df)
    
    # Extract unique genes and map to Entrez
    genes <- unique(protein_df$Gene_symbol[protein_df$Score >= (as.numeric(selected_score) * 1000)])
    genes <- genes[!is.na(genes) & genes != "" & genes != "Unknown"]
    
    entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    req(nrow(entrez_map) > 0)
    
    withProgress(message = "Running Pathway Enrichment...", {
      # Parallel execution of KEGG and Reactome
      results <- mclapply(c("KEGG", "Reactome"), function(type) {
        if (type == "KEGG") {
          res <- enrichKEGG(gene = entrez_map$ENTREZID, organism = 'hsa', pvalueCutoff = pval)
        } else {
          res <- enrichPathway(gene = entrez_map$ENTREZID, organism = 'human', pvalueCutoff = pval, readable = TRUE)
        }
        return(if(!is.null(res)) as.data.frame(res) else NULL)
      }, mc.cores = 2)
      
      pathway_data(list(kegg = results[[1]], reactome = results[[2]]))
      
      # Update Filter Choices
      mapping <- load_shared_data("smiles_cid_mapping")
      if (!is.null(mapping)) {
        updateSelectInput(session, "select_cid_kegg", choices = mapping$CID)
        updateSelectInput(session, "select_cid_reactome", choices = mapping$CID)
      }
      
      shinyjs::runjs("if(window.updateProgress) { window.moduleCompletion.pathway = true; window.updateProgress(); }")
    })
  })

  # 2. Table Rendering Helper
  render_pathway_dt <- function(df, type) {
    req(df)
    if (nrow(df) == 0) return(datatable(data.frame(Message = "No significant pathways found.")))
    
    df_display <- df %>%
      mutate(
        ID = if(type == "KEGG") 
               sprintf("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+%s' target='_blank'>%s</a>", ID, ID)
             else 
               sprintf("<a href='https://reactome.org/content/detail/%s' target='_blank'>%s</a>", ID, ID),
        geneID = sapply(geneID, format_pathway_genes)
      ) %>%
      select(ID, Description, p.adjust, GeneRatio, geneID)

    datatable(df_display, escape = FALSE, options = list(scrollX = TRUE, pageLength = 10))
  }

  output$KEGG_Table <- renderDT({ render_pathway_dt(pathway_data()$kegg, "KEGG") })
  output$Reactome_Table <- renderDT({ render_pathway_dt(pathway_data()$reactome, "Reactome") })

  # 3. Export Handlers
  output$KEGG_Download <- downloadHandler(
    filename = "KEGG_Pathways.csv",
    content = function(file) write.csv(pathway_data()$kegg, file, row.names = FALSE)
  )
  
  output$Reactome_Download <- downloadHandler(
    filename = "Reactome_Pathways.csv",
    content = function(file) write.csv(pathway_data()$reactome, file, row.names = FALSE)
  )
}
