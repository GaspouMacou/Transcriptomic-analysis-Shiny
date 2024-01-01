library(patchwork)
library(shiny)
library(Matrix)
library(shinydashboard)
library(ggplot2)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(ggrepel)
library(stringr)
library(dplyr)
library(Seurat)
library(hdf5r)
library(readr)
library(clustree)
library(tidyverse)
library(fontawesome)
library(cowplot)
library(Signac)
library(DT)
library(SeuratDisk)
library(data.table)
library(plotly)
library(colourpicker)
library(httr)
library(reticulate)
library(anndata)
library(SeuratObject)
#library(SeuratData)
library(jsonlite)
#library(monocle3)
#library(SeuratWrappers)
#Genomes
#library(EnsDb.Hsapiens.v86)
#library(biovizBase)

# trajectory_server.R

trajectory_server <-  function(input, output, session) {

    shinyjs::useShinyjs() #button deactivation


    # Tab 19
    seurat_monocle <- reactiveVal() # Seurat object loaded by the user
    monocle_object <- reactiveVal() # Monocle object converted by the user from a seurat object
    selected_cell_id <- reactiveVal(NULL) # Cell id to connect the pointer and the id of cell when the user is clicking on the umap for the trajectory


    # Désactiver initialement les boutons
    shinyjs::disable("convertToMonocle")
    shinyjs::disable("constructGraph")
    shinyjs::disable("startRootSelection")

    # For Seurat loading
    observeEvent(input$load_seurat_file_monocle, {
      message("Attempting to read file at: ", input$load_seurat_file_monocle$datapath)

      tryCatch({
        temp_seurat <- readRDS(input$load_seurat_file_monocle$datapath)
        message("File successfully read.")
        message("Log: Seurat object loaded.")

        seurat_monocle(temp_seurat)
        showNotification("Seurat object has been successfully loaded!")

        shinyjs::enable("convertToMonocle")
      }, error = function(e) {
        message("Log: Error loading Seurat object: ", e$message)
        showNotification(paste("An error occurred: ", e$message), type = "error")
      })
    })

    # For conversion
    observeEvent(input$convertToMonocle, {
      req(seurat_monocle())

      tryCatch({
        # Convert Seurat object to Monocle CellDataSet object using SeuratWrappers
        sce <- as.cell_data_set(seurat_monocle())
        monocle_temp <- as(sce, "cell_data_set")

        # Estimate size factors using Monocle
        monocle_temp <- estimate_size_factors(monocle_temp)

        # Add gene names to the Monocle object
        rowData(monocle_temp)$gene_short_name <- rownames(seurat_monocle()[["RNA"]])

        monocle_object(monocle_temp)

        showNotification("Converted to Monocle object successfully!")
        message("Log: Conversion to Monocle successful.")

        # Enable the buttons after conversion
        shinyjs::enable("constructGraph")
      }, error = function(e) {
        showNotification(paste("An error occurred during conversion: ", e$message), type = "error")
        message("Log: Error during conversion:", e$message)
      })
    })



    cell_data <- reactive({
      data <- as.data.frame(reducedDims(monocle_object())$UMAP)
      data$cell_id <- as.character(rownames(data))
      return(data)
    })

    get_cell_id_from_coordinates <- function(x, y, cell_data_func) {
      data <- cell_data_func()
      distances <- rowSums((data[,1:2] - c(x, y))^2)
      closest_cell <- which.min(distances)
      return(data$cell_id[closest_cell])
    }


    observeEvent(input$constructGraph, {
      req(monocle_object())

      k_value <- as.integer(input$cluster_k)
      monocle_temp <- cluster_cells(monocle_object(), k = k_value)
      monocle_temp <- learn_graph(monocle_temp)

      monocle_object(monocle_temp)

      # Créez le graphique sans utiliser geom_text_repel()
      p <- plot_cells(monocle_temp, color_cells_by = "cluster", cell_size = 0.5, group_label_size = 3, label_groups_by_cluster = F, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 5)
      p <- p + theme_void() # Pour retirer les éléments inutiles

      # Convertir en plotly et enregistrer l'événement
      p_plotly <- ggplotly(p, source = "trajectoryPlot") %>% event_register("plotly_click")

      output$trajectoryPlot <- renderPlotly({
        p_plotly
      })

      showNotification("Graph construction completed!")
      cat("Log: Graph construction completed.\n")
      shinyjs::enable("startRootSelection")

    })

    root_selection_started <- reactiveVal(FALSE)

    observeEvent(input$startRootSelection, {
      showNotification("Please click on the graphic to select the root cell.")
      root_selection_started(TRUE)
    })
    observeEvent(event_data("plotly_click", source = "trajectoryPlot"), {
      if(root_selection_started()) {
        clicked_point <- event_data("plotly_click", source = "trajectoryPlot")
        # Ajout des instructions de débogage ici
        cat("Clicked coordinates: x = ", clicked_point$x, ", y = ", clicked_point$y, "\n")
        print(head(cell_data()))

        # Utiliser la fonction pour obtenir l'ID de cellule
        sel_id <- get_cell_id_from_coordinates(clicked_point$x, clicked_point$y, cell_data)
        cat("Selected ID from coordinates: ", sel_id, "\n")

        selected_cell_id(get_cell_id_from_coordinates(clicked_point$x, clicked_point$y, cell_data))
        cat("Log: Selected root cell ID: ", sel_id, "\n")

        if (is.null(sel_id)) {
          cat("Log: No root cell selected.\n")
        } else {
          showModal(modalDialog(
            title = "Confirm root selection",
            paste("Have you selected the right point with the cell ID:", sel_id, "?"),
            footer = tagList(
              actionButton("confirmRoot", "Confirm"),
              modalButton("Annuler")
            )
          ))
        }
      }
    })

    observeEvent(input$confirmRoot, {
      removeModal()
      root_selection_started(FALSE)

      req(monocle_object())
      if (is.null(selected_cell_id())) {
        showNotification("No root cell selected!", type = "error")
        cat("Log: No root cell selected.\n")
        return()
      }

      monocle_temp <- order_cells(monocle_object(), reduction_method = "UMAP", root_cells = selected_cell_id())
      output$pseudotimePlot <- renderPlot({
        plot_cells(monocle_temp, color_cells_by = "pseudotime", cell_size = 0.5, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 4)
      })

      showNotification("Cell ordering completed!")
      cat("Log: Cell ordering completed.\n")
    })

    #Onglet 20

    pseudotime_diff_gene_results <- reactiveVal()
    observeEvent(input$diffGeneTestButton, {
      req(monocle_object())
      
      # Effectuez le test des gènes différentiels
      diff_gene_results_temp <- graph_test(monocle_object(), neighbor_graph = "principal_graph", cores = 4)
      
      # Filtrez les résultats pour obtenir ceux avec q_value < 0.05
      pr_deg_ids <- row.names(subset(diff_gene_results_temp, q_value < 0.05))
      
      # Mettez à jour les résultats des tests de gènes différentiels filtrés
      pseudotime_diff_gene_results(subset(diff_gene_results_temp, row.names(diff_gene_results_temp) %in% pr_deg_ids))
      
      # Mettez à jour l'affichage du tableau dans l'UI
      output$diffGeneTable <- DT::renderDataTable({
        DT::datatable(pseudotime_diff_gene_results())
      })
    })
    
    
    output$download_pseudotime_diff_genes <- downloadHandler(
      filename = function() {
        paste("pseudotime-diff-genes-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        tryCatch({
          diff_gene_results <- pseudotime_diff_gene_results()
          req(!is.null(diff_gene_results))
          
          # Formattez les colonnes comme vous le souhaitez
          diff_gene_results$p_val <- formatC(diff_gene_results$p_val, format = "e", digits = 2)
          diff_gene_results$p_val_adj <- formatC(diff_gene_results$p_val_adj, format = "e", digits = 2)
          
          # Trier le tableau et sauvegarder
          diff_gene_results_sorted <- diff_gene_results[order(diff_gene_results$p_val), ]
          write.csv(diff_gene_results_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Erreur lors du téléchargement du tableau de gènes différentiels de pseudotime: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )
    
    observeEvent(input$geneSelectButton, {
      req(monocle_object())

      # Récupérer le gène saisi
      entered_gene <- trimws(input$geneTextInput)

      # Afficher pour le debug
      print(paste("Entered gene:", entered_gene))

      # Vérifiez si le gène saisi est non vide et existe dans vos données
      if (entered_gene == "" || !(entered_gene %in% rownames(rowData(monocle_object())))) {
        showNotification("Entered gene is not valid or not found in the data!", type = "error")
        return(NULL)
      }

      # Générer le graphique
      output$geneTrajectoryPlot <- renderPlot({
        plot_cells(monocle_object(), genes = entered_gene, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)
      })
    })


}
