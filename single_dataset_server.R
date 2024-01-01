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

# single_dataset_server.R

single_dataset_server <-  function(input, output, session) {
  shinyjs::useShinyjs() 
  shinyjs::disable("pca_button")
  shinyjs::disable("scale_button")

    # Reactives objects for the single dataset part
    single_dataset_object <- reactiveVal(NULL)  # This variable will stock the seurat object for all the analysis
    gene_list <- reactiveValues(features = NULL) # list of all the genes sequenced



    # Tab 1 : Loading Data

    # Loading raw data
    observeEvent(input$file, {
      showNotification("Uploading and processing data...", type = "message")
      mt_pattern <- ifelse(input$species_choice == "mouse", "^mt-", "^MT-")
      if (dir.exists("unzipped")) {
        unlink("unzipped", recursive = TRUE)
      }
      zipPath <- input$file$datapath
      dataDir <- "unzipped"
      unzip(zipPath, exdir = dataDir)
      single_dataset_object.data <- Read10X(dataDir)
      if (input$dataset_type == "snRNA") {
        seuratObj <- CreateSeuratObject(counts = single_dataset_object.data, project = "single_dataset", min.cells = 3, min.features = 200)
        seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
        single_dataset_object(seuratObj)
        showNotification("snRNA-seq data processed successfully.", type = "message")
      } else if (input$dataset_type == "multiome") {
        rna.data <- single_dataset_object.data$`Gene Expression`
        atac.data <- single_dataset_object.data$`Peaks`
        rna.seurat <- CreateSeuratObject(counts = rna.data, project = "RNA", min.cells = 3, min.features = 200)
        atac.seurat <- CreateSeuratObject(counts = atac.data, project = "ATAC", min.cells = 3, min.features = 200)
        rna.seurat[["percent.mt"]] <- PercentageFeatureSet(rna.seurat, pattern = mt_pattern)
        single_dataset_object(rna.seurat)
        showNotification("Multiome data processed successfully.", type = "message")
      }
    })

    # Loading seurat object
    observeEvent(input$load_seurat_file, {
      message("Attempting to read file at: ", input$load_seurat_file$datapath)
      tryCatch({
        loaded_seurat <- readRDS(input$load_seurat_file$datapath)
        message("File successfully read.")
        single_dataset_object(loaded_seurat)
        showNotification("L'objet Seurat a été chargé avec succès!")

      }, error = function(e) {
        showNotification(paste("An error occurred: ", e), type = "error")
      })
    })


    # Tab 2 : QC metrics and normalization

    #number of nuclei:
    nuclei_count <- reactive({
      req(single_dataset_object())
      seuratRNA_subset <-subset(single_dataset_object(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
      if (exists("seuratRNA_subset")) {
        return(dim(seuratRNA_subset)[2])
      } else {
        return(0)
      }
    })
    output$nuclei_count <- renderInfoBox({
      # S'assurer que la fonction nuclei_count est exécutée et renvoie une valeur
      count <- tryCatch({
        count_value <- nuclei_count()  # Appel de la fonction pour obtenir la valeur
        infoBox(
          title = HTML("Number<br>of Nuclei"),  # Utiliser HTML avec br() pour le saut de ligne
          value = count_value,  # Utiliser la valeur retournée ici
          color = "blue",
          icon = icon("dna"),
          width = 12  # Augmenter la largeur ici
        )
      }, error = function(e) {
        # Gérer l'erreur si nécessaire
        infoBox(
          title = HTML("Number<br>of Nuclei"),
          value = "Error",
          icon = icon("dna"),
          color = "blue",
          width = 12
        )
      })
    })
    
    # Display QC metrics on a VlnPlot chart
    output$vlnplot <- renderPlot({
      tryCatch({
        if (input$QCmetrics) {
            VlnPlot(single_dataset_object(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        }
      }, error = function(e) {
        showNotification(paste0("Error displaying VlnPlot: ", e$message), type = "error")
      })
    })

    # Display the scatter plot 1
    output$scatter_plot1 <- renderPlot({
      tryCatch({
        if (input$show_plots) {
            seuratRNA_subset <- subset(single_dataset_object(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
            FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")

        }
      }, error = function(e) {
        showNotification(paste0("Erreur lors de l'affichage du scatter plot 1 : ", e$message), type = "error")
      })
    })

    # Display the scatter plot 2
    output$scatter_plot2 <- renderPlot({
      tryCatch({
        if (input$show_plots) {
            seuratRNA_subset <- subset(single_dataset_object(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
            FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        }
      }, error = function(e) {
        showNotification(paste0("Error displaying scatter plot 2: ", e$message), type = "error")
      })
    })

    # Apllying filtering with the qc aprameter
    observeEvent(input$apply_qc, {
      seurat_obj <-
        single_dataset_object()
      seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > input$nFeature_range[1] &
                             nFeature_RNA < input$nFeature_range[2] &
                             percent.mt < input$percent.mt_max)
      single_dataset_object(seurat_obj)
      num_nuclei_after_qc <- dim(seurat_obj)[2]
      msg <- paste("QC filters applied. Number of nuclei retained:", num_nuclei_after_qc)
      showNotification(msg, type = "message")
    })

    # Function for normalizing data and finding variable features
    observeEvent(input$normalize_data, {
      req(single_dataset_object())
      tryCatch({
        withProgress(message = "Normalizing data...", value = 0, {
          # Normalisation des données
          normalized_seurat <- NormalizeData(single_dataset_object(), normalization.method = "LogNormalize", scale.factor = input$scale_factor)
          incProgress(0.5)
          print("Data normalization completed")

          # Identification of variable features
          normalized_seurat <- FindVariableFeatures(normalized_seurat, selection.method = "vst", nfeatures = 2000)
          incProgress(0.5)
          print("Variable features identification completed")
          single_dataset_object(normalized_seurat)
          shinyjs::enable("scale_button")
        })
        showNotification("Data normalization and variable features identification completed successfully.", type = "message")
      }, error = function(e) {
        print(paste0("Error during data normalization or variable features identification: ", e$message))
        shinyjs::disable("scale_button")
        showNotification(paste0("Error during data normalization or variable features identification: ", e$message), type = "error")
      })
    })



    # Plot rendering of variable features
    output$variable_feature_plot <- renderPlot({
      req(input$show_plots, single_dataset_object())
      if (!is.null(single_dataset_object()) && length(VariableFeatures(single_dataset_object())) > 0) {
        plot1 <- VariableFeaturePlot(single_dataset_object())
        top10 <- head(VariableFeatures(single_dataset_object()), 10)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        plot2
      } else {
        return(NULL)
      }
    })




    #  Third tab: Scaling and PCA

    shinyjs::disable("pca_button")

    # Scaling data
    observeEvent(input$scale_button, {
      shinyjs::disable("pca_button")
      showNotification("Scaling data...", type = "message")
      tryCatch({
        all_genes <- rownames(single_dataset_object())
        print("Scaling data for all genes")
        scaled_seurat <-ScaleData(single_dataset_object(), features = all_genes)
        single_dataset_object(scaled_seurat)
        print("Data scaling completed")
        shinyjs::enable("pca_button")
        showNotification("Data scaling completed successfully.", type = "message")
      }, error = function(e) {
        print(paste("Error during data scaling:", e$message))
        shinyjs::disable("pca_button")
        showNotification(paste("Error during data scaling:", e$message), type = "error")
      })
    })

    # PCA calculation
    observeEvent(input$pca_button, {
      req(single_dataset_object())
      showNotification("Running PCA...", type = "message")
      tryCatch({
        pca_result <- RunPCA(single_dataset_object(), features = VariableFeatures(single_dataset_object()))
        single_dataset_object(pca_result)
      }, error = function(e) {
        showNotification(paste("Error running PCA:", e$message), type = "error")
      })
    })

    #Printing PCA results
    output$pca_results <- renderPrint({
      req(!is.null(single_dataset_object()))
      if ("pca" %in% names(single_dataset_object())) {
        print(single_dataset_object()$pca, dims = 1:5, nfeatures = 5)
      }
    })

    #Vizdimloading plot
    output$loading_plot <- renderPlot({
      tryCatch(VizDimLoadings(single_dataset_object(), dims = 1:2, reduction = "pca"), error = function(e){})
    })

    # Dimplot
    output$dim_plot <- renderPlot({
      tryCatch(DimPlot(single_dataset_object(), reduction = "pca"), error = function(e){})
    })


# 4 tab: Heatmaps

    observe({
      # Vérifiez si single_dataset_object() est non-null et un objet Seurat
      if (!is.null(single_dataset_object()) && inherits(single_dataset_object(), "Seurat")) {
        # Mettre à jour les choix de gènes pour la heatmap
        updateSelectInput(session, "gene_select_heatmap",
                          choices = c("", rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts'))))
      }
    })

    observe({
      # Vérifiez si single_dataset_object() est non-null et un objet Seurat
      if (!is.null(single_dataset_object()) && inherits(single_dataset_object(), "Seurat")) {
        # Mettre à jour le champ de texte avec les gènes sélectionnés
        selected_genes <- paste(input$gene_select_heatmap, collapse = ", ")
        updateTextInput(session, "selected_genes_heatmap", value = selected_genes)
      }
    })

    observeEvent(input$runHeatmap, {
      # Vérifiez si single_dataset_object() est non-null et un objet Seurat
      if (!is.null(single_dataset_object()) && inherits(single_dataset_object(), "Seurat")) {
        gene_list <- unlist(strsplit(input$selected_genes_heatmap, ","))
        gene_list <- unique(trimws(gene_list)) # Éliminer les espaces et doublons
        if (length(gene_list) > 0) {
          output$heatmap <- renderPlot({
            # Utiliser DoHeatmap pour générer la heatmap
            DoHeatmap(single_dataset_object(), features = gene_list)
          })
        } else {
          showNotification("Please select at least one gene.", type = "warning")
        }
      } else {
        showNotification("Seurat object not loaded or invalid.", type = "error")
      }
    })



# 5th tab: elbowplot an jacktraw


    observeEvent(input$run_jackstraw, {
      showNotification("Starting JackStraw analysis...", type = "message", duration = 3)
      withProgress(message = "Running JackStraw analysis...", value = 0, {
        seurat_tmp <- single_dataset_object()
        seurat_tmp <- JackStraw(seurat_tmp, num.replicate = 100)
        seurat_tmp <- ScoreJackStraw(seurat_tmp, dims = 1:input$num_dims)
        single_dataset_object(seurat_tmp)
        incProgress(1)
      })

      # Display the JackStraw plot
      output$jackstraw_plot <- renderPlot({
        req(single_dataset_object())
        JackStrawPlot(single_dataset_object(), dims = 1:input$num_dims)
      })
      showNotification("JackStraw analysis completed.", type = "message")
    })

    # Show the ElbowPlot
    observeEvent(input$run_elbow, {
      output$elbow_plot <- renderPlot({
        req(single_dataset_object())
        ElbowPlot(single_dataset_object(), ndims = 50)
      })
    })

    # 6th tab: Neighbors calculation and clustering

    # Finding neighbors
    observeEvent(input$run_neighbors, {
      tryCatch({
        req(!is.null(single_dataset_object()))
        showNotification("Looking for neighbors", type = "message", duration = 5)
        seurat_tmp <- single_dataset_object()
        seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:input$dimension_1)
        seurat_tmp <- RunUMAP(seurat_tmp, dims = 1:input$dimension_1)
        single_dataset_object(seurat_tmp)
        incProgress(1)
        print("Neighbors found and UMAP calculated.")  # Debug line
        output$neighbors_plot <- renderPlot({
          DimPlot(single_dataset_object(), group.by = "ident")
        })
      }, error = function(e) {
        showNotification(paste0("Neighbor search error : ", e$message), type = "error")
      })
    })

    # Clustering
    observeEvent(input$run_clustering, {
      tryCatch({
        req(!is.null(single_dataset_object()))
        showNotification("Clustering process started...", type = "message", duration = 5)
        seurat_tmp <- single_dataset_object()
        seurat_tmp <- FindClusters(seurat_tmp, resolution = input$resolution_step1)
        single_dataset_object(seurat_tmp)
        print("Clustering performed.")  # Debug line
      }, error = function(e) {
        showNotification(paste0("Clustering error: ", e$message), type = "error")
      })
    })


    # 7th tab: Calculate differentially expressed genes for each cluster


    # Variables réactives pour cet onglet
    gene_tables <- reactiveValues()
    all_markers <- reactiveVal()

    # Function to calculate all differential markers once only
    calculate_all_markers <- function() {
      if (is.null(single_dataset_object())) {
        showNotification("single_dataset_object is NULL, please check the previous steps.", type = "error")
        return(NULL)
      }
      tryCatch({
        showNotification("Calculation of differential markers in progress...", type = "message", duration = 10)
        markers <- FindAllMarkers(single_dataset_object(),
                                  min.pct = input$min_pct_all_single,
                                  logfc.threshold = input$logfc_threshold_all_single,)
        markers <- as.data.frame(markers)
        markers$gene <- paste0('<span class="gene-name">', markers$gene, '</span>')
        markers$p_val <- format(markers$p_val, scientific = TRUE, digits = )
        markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
        markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(markers), '">', rownames(markers), '</a>')
        all_markers(markers)
      }, error = function(e) {
        showNotification(paste0("Error when calculating differential markers: ", e$message), type = "error")
      })
    }

    # Function to update gene tables according to the number selected
    update_gene_tables_display <- function() {
      markers <- all_markers()
      num_genes_to_display <- input$number_genes
      for (cluster in unique(markers$cluster)) {
        gene_tables[[paste0("table_", cluster)]] <- head(markers[markers$cluster == cluster, ], n = num_genes_to_display)
      }
      output$diff_genes_tables <- renderUI({
        tagList(
          lapply(names(gene_tables), function(name) {
            tags$div(style = "width: 100%; font-size: 75%;",
                     tagList(
                       h3(paste0("Cluster ", stringr::str_replace(name, "table_", ""))),
                       DTOutput(name),
                       hr()
                     )
            )
          })
        )
      })
    }

    # Observer for run_DE button to calculate all markers
    observeEvent(input$run_DE, {
      tryCatch({
        calculate_all_markers()
        update_gene_tables_display()      }, error = function(e) {
        showNotification(paste0("Error executing run_DE:", e$message), type = "error")
      })
    })

    # Observer to update display only when gene number changes
    observeEvent(input$number_genes, {
      if (!is.null(input$number_genes) && !is.null(all_markers())) {
        update_gene_tables_display()
      }
    })

    #update the output genes table
    observe({
      lapply(names(gene_tables), function(name) {
        output[[name]] <- renderDataTable({
          datatable(gene_tables[[name]], escape = FALSE)
        })
      })
    })



    # Download seurat object
    output$save_seurat <- downloadHandler(
      filename = function() {
        showModal(modalDialog(
          title = "Please wait",
          "Wait, saving of seurat object in process...",
          footer = NULL,
          easyClose = FALSE
        ))
        return("seurat_object.rds")
      },
      content = function(file) {
        saveRDS(single_dataset_object(), file)
        removeModal()
      }
    )


    # Download markers as a CSV
    output$download_DE <- downloadHandler(
      filename = function() {
        paste("DE_genes_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        markers <- all_markers()
        write.csv(markers, file)
      }
    )

    # Tab 8 Visualization of expressed genes::

    #VLN plot for each gene selected
    observeEvent(input$show_vln, {
      print("show_vln triggered")
      req(!is.null(single_dataset_object()))
      gene_list <- strsplit(input$gene_list_vln, ",")[[1]]
      gene_list <- unique(trimws(gene_list))  # Éliminer les doublons
      print(paste("gene_list:", paste(gene_list, collapse = ", ")))
      withProgress(message = "Generating VlnPlot...", value = 0, {
        output$vln_plot <- renderPlot({
          if(input$hide_vln_points) {
            VlnPlot(single_dataset_object(), features = gene_list, log = TRUE, pt.size = 0)
          } else {
            VlnPlot(single_dataset_object(), features = gene_list, log = TRUE)
          }
        })
        incProgress(1)
      })
    })

    observe({
      req(!is.null(single_dataset_object()))
      unique_idents <- unique(Idents(single_dataset_object()))
      updateSelectInput(session, "ident_1", choices = unique_idents)
      updateCheckboxGroupInput(session, "ident_2", choices = unique_idents, selected = "")
    })

    #FeaturePlot plot for each gene selected
    observeEvent(input$show_feature, {
      req(input$gene_list_feature, single_dataset_object())
      idents <- strsplit(input$gene_list_feature, ",")[[1]]
      idents <- unique(trimws(idents))  # Éliminer les doublons
      if(all(idents %in% rownames(single_dataset_object()))) {
        output$feature_plot <- renderPlot({
          FeaturePlot(single_dataset_object(), features = idents)
        })
      } else {
        print("The requested genes are not present in the dataset.")
      }
    })

    #Gene selector
    observeEvent(single_dataset_object(), {
      gene_list <- c("", rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts')))
      updateSelectInput(session,
                        "gene_select",
                        "Select a gene:",
                        choices = gene_list,
                        selected = "")
    })

    observeEvent(input$gene_select, {
      # Get the selected gene
      selected_gene <- input$gene_select
      update_genes <- function(existing_genes, new_gene) {
        if (existing_genes == "") {
          return(new_gene)
        } else {
          genes <- strsplit(existing_genes, ",")[[1]]
          genes <- unique(c(trimws(genes), new_gene))  # Add new gene and eliminate duplicates
          return(paste(genes, collapse = ", "))
        }
      }

      updateTextInput(session, "gene_list_vln", value = update_genes(input$gene_list_vln, selected_gene))
      updateTextInput(session, "gene_list_feature", value = update_genes(input$gene_list_feature, selected_gene))
      updateTextInput(session, "gene_list_dotplot", value = update_genes(input$gene_list_dotplot, selected_gene))
    })

    #DotPlot for each gene selected
    observeEvent(input$show_dot, {
      req(!is.null(single_dataset_object()))
      gene_list <- strsplit(input$gene_list_dotplot, ",")[[1]]
      gene_list <- unique(trimws(gene_list))  # Éliminer les doublons

      # Generate the Dot plot
      output$dot_plot <- renderPlot({
        features <- gene_list
        DotPlot(single_dataset_object(), features = features) + RotatedAxis()
      })
    })

    # Tab 9: Final UMAP

    #Reactive variable for that tab
    cluster_colours <- reactiveVal()  # Initialize a list to store cluster colors

    observe({
      if (!is.null(single_dataset_object())) {
        updateSelectInput(session, "select_cluster", choices = unique(Idents(single_dataset_object())))
      }
    })

    #function for renaming each cluster
    observeEvent(input$rename_single_cluster_button, {
      req(input$select_cluster, input$rename_single_cluster, single_dataset_object())

      updated_seurat <- single_dataset_object()

      # Si le nouveau nom de cluster existe déjà, fusionnez les clusters
      if (input$rename_single_cluster %in% unique(Idents(updated_seurat))) {
        cells_to_merge <- which(Idents(updated_seurat) %in% c(input$select_cluster, input$rename_single_cluster))
        Idents(updated_seurat, cells = cells_to_merge) <- input$rename_single_cluster
        showNotification("Clusters merged under the name: ", input$rename_single_cluster, type = "message")
      } else {
        # Sinon, renommez simplement le cluster sélectionné
        Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster)) <- input$rename_single_cluster
        showNotification("Cluster renamed to: ", input$rename_single_cluster, type = "message")
      }

      single_dataset_object(updated_seurat)  # Mettez à jour l'objet Seurat renommé
      # Mettez à jour les options de sélection des clusters avec les nouveaux noms de clusters
      updateSelectInput(session, "select_cluster", choices = unique(Idents(single_dataset_object())))
      updateSelectInput(session, "cluster_select", choices = unique(Idents(single_dataset_object())))
    })


    observe({
      if (!is.null(single_dataset_object())) {
        message("Update cluster colors")

        current_colours <- scales::hue_pal()(length(unique(Idents(single_dataset_object()))))
        names(current_colours) <- sort(unique(Idents(single_dataset_object())))
        cluster_colours(current_colours)
      }
    })

    observeEvent(input$update_colour, {
      message("Update the color of the selected cluster")

      current_colours <- cluster_colours()
      current_colours[input$cluster_select] <- input$cluster_colour
      cluster_colours(current_colours)
    })

    output$umap_finale <- renderPlotly({
      req(single_dataset_object())
      message("Generating the final UMAP")


      plot_data <- DimPlot(single_dataset_object(), group.by = "ident", pt.size = input$pt_size, label = TRUE, label.size = input$label_font_size) +
        theme(axis.line = element_line(size = 0.5)) +
        scale_color_manual(values = cluster_colours())

      plot_data <- plot_data + NoLegend()
      interactive_plot <- ggplotly(plot_data, tooltip = "text")

      interactive_plot <- interactive_plot %>%
        layout(
          title = list(text = input$plot_title, font = list(size = 24)),
          hovermode = "closest"
        )

      message("UMAP final généré")
      return(interactive_plot)
    })




    # Onglet 10:   Find markers for a specific cluster
    # Variables réactives pour cet onglet
    gene_tables_10 <- reactiveValues()
    all_markers_10 <- reactiveVal()

    # Ajout de l'UMAP en haut de l'onglet
    output$umap_cluster_10 <- renderPlot({
      req(single_dataset_object())
      message("UMAP generation for a specific cluster")

      plot_data <- DimPlot(single_dataset_object(), group.by = "ident", pt.size = input$pt_size, label = TRUE, label.size = input$label_font_size) +
        theme(axis.line = element_line(size = 0.5)) +
        scale_color_manual(values = cluster_colours())

      plot_data <- plot_data + NoLegend()

      message("UMAP for a specific cluster generated")
      return(plot_data)
    })

    # Mettre à jour les choix pour le selectInput cluster
    observe({
      req(single_dataset_object())
      updateSelectInput(session, "cluster", choices = levels(single_dataset_object()))
    })



    # Mettre à jour les choix pour le selectInput cluster
    observe({
      req(single_dataset_object())
      updateSelectInput(session, "cluster", choices = levels(single_dataset_object()))
    })

    # Fonction pour calculer tous les marqueurs différentiels pour le cluster choisi
    calculate_markers_for_cluster <- function() {
      # Vérifier que les objets et les entrées nécessaires sont disponibles
      req(single_dataset_object(), input$cluster, input$min_pct_single, input$logfc_threshold_single)

      tryCatch({
        # Calcul des marqueurs différentiels
        markers <- FindMarkers(single_dataset_object(),
                               ident.1 = input$cluster,
                               min.pct = input$min_pct_single,
                               logfc.threshold = input$logfc_threshold_single)

        # Formatage des valeurs p et ajustement des valeurs p pour l'affichage
        markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
        markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
        markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(markers), '">', rownames(markers), '</a>')

        # Stocker les résultats dans une valeur réactive
        all_markers_10(markers)

        # Retourner TRUE en cas de succès
        return(TRUE)
      }, error = function(e) {
        # Gestion des erreurs, par exemple, affichage d'une notification dans l'application
        showNotification(paste0("Error in calculate_markers_for_cluster: ", e$message), type = "error")
        # Retourner FALSE en cas d'erreur
        return(FALSE)
      })
    }


    # Observeur pour le bouton find_markers
    observeEvent(input$find_markers, {
      calculate_markers_for_cluster()
      update_gene_tables_display_10()
    })



    # Fonction pour mettre à jour le tableau de gènes pour l'onglet 10
    update_gene_tables_display_10 <- function() {
      markers <- all_markers_10()
      num_genes_to_display_10 <- input$number_genes_10

      gene_tables_10$table <- head(markers, n = num_genes_to_display_10)

      output$gene_tables_10 <- renderUI({
        tags$div(style = "width: 100%; font-size: 75%;",
                 tagList(
                   DTOutput("table_10")
                 )
        )
      })
    }
    # Observeur pour la mise à jour de l'affichage uniquement quand le nombre de gènes change
    observeEvent(input$number_genes_10, {
      if (!is.null(input$number_genes_10) && !is.null(all_markers_10())) {
        update_gene_tables_display_10()
      }
    })

    # Ce bloc s'exécute pour mettre à jour le tableau de sortie pour l'onglet 10
    observe({
      output$table_10 <- renderDataTable({
        datatable(gene_tables_10$table, escape = FALSE)
      })
    })

    observe({
      req(single_dataset_object())
      updateSelectInput(session, "cluster1", choices = levels(single_dataset_object()))
    })

    # Cluster selection
    observeEvent(input$cluster1, {
      updateSelectInput(session, "cluster2", choices = setdiff(levels(single_dataset_object()), input$cluster1))
    })

    # Variables réactives pour cet onglet
    gene_tables_new <- reactiveValues()
    all_markers_new <- reactiveVal()

    # Fonction pour calculer tous les marqueurs différentiels pour la comparaison de clusters
    calculate_markers_for_comparison <- function() {
      # Vérifier que les objets et les entrées nécessaires sont disponibles
      req(single_dataset_object(), input$cluster1, input$cluster2, input$min_pct_comparison, input$logfc_threshold_comparison)

      tryCatch({
        # Vérifier que les clusters choisis sont différents
        if (input$cluster1 == input$cluster2) {
          showNotification("Les clusters sélectionnés ne doivent pas être identiques.", type = "error")
          return(FALSE)
        }

        # Calcul des marqueurs différentiels
        markers <- FindMarkers(single_dataset_object(),
                               ident.1 = input$cluster1,
                               ident.2 = input$cluster2,
                               min.pct = input$min_pct_comparison,
                               logfc.threshold = input$logfc_threshold_comparison)

        # Formatage des valeurs p et ajustement des valeurs p pour l'affichage
        markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
        markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
        markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(markers), '">', rownames(markers), '</a>')

        # Stocker les résultats dans une valeur réactive
        all_markers_new(markers)

        # Retourner TRUE en cas de succès
        return(TRUE)
      }, error = function(e) {
        # Gestion des erreurs, par exemple, affichage d'une notification dans l'application
        showNotification(paste0("Error in calculate_markers_for_comparison: ", e$message), type = "error")
        # Retourner FALSE en cas d'erreur
        return(FALSE)
      })
    }


    # Fonction pour mettre à jour le tableau de gènes pour cet onglet
    update_gene_tables_display_new <- function() {
      markers <- all_markers_new()
      num_genes_to_display_new <- input$n_diff_markers

      gene_tables_new$table <- head(markers, n = num_genes_to_display_new)

      output$gene_tables_new <- renderUI({
        tags$div(style = "width: 100%; font-size: 75%;",
                 tagList(
                   DTOutput("table_new")
                 )
        )
      })
    }

    # Observer for compare_markers button
    observeEvent(input$compare_markers, {
      calculate_markers_for_comparison()
      update_gene_tables_display_new()
    })

    # Observer to update display only when gene count changes
    observeEvent(input$n_diff_markers, {
      if (!is.null(input$n_diff_markers) && !is.null(all_markers_new())) {
        update_gene_tables_display_new()
      }
    })

    # Observer to update the output table for this tab
    observe({
      output$table_new <- renderDataTable({
        datatable(gene_tables_new$table, escape = FALSE)
      })
    })



    # Download gene comparison table
    output$download_markers_single_cluster <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster, "-VS-", "all_others_clusters", "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        tryCatch({
          req(!is.null(all_markers_10))
          all_markers_10_temp <- all_markers_10()
          all_markers_10_temp_sorted <- all_markers_10_temp[order(all_markers_10_temp$p_val), ]
          write.csv(all_markers_10_temp_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )

    # Download gene comparison table
    output$download_markers_multiple_clusters <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster1, "-VS-", input$cluster2, "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        tryCatch({
          req(!is.null(all_markers_new))

          all_markers_new_temp <- all_markers_new()
          all_markers_new_temp_sorted <- all_markers_new_temp[order(all_markers_new_temp$p_val), ]
          write.csv(all_markers_new_temp_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )
    # Onglet 11 : Subseting

    # Variables wich stock the subset of single_dataset_object()
    subset_seurat <- reactiveVal()

    observe({
{
        subset_seurat(single_dataset_object())
      }
    })

    observe({
      if (!is.null(single_dataset_object())) {
        updateSelectInput(session, "select_ident_subset", choices = unique(Idents(single_dataset_object())))
      }
    })

    #Run subset
    observeEvent(input$apply_subset, {
      req(single_dataset_object())
      subset_seurat_temp <- single_dataset_object()

      if (length(input$select_ident_subset) > 0) {
        subset_seurat_temp <- subset(x = subset_seurat_temp, idents = input$select_ident_subset)

        if (nrow(subset_seurat_temp@meta.data) == 0) {
          showNotification("No cells found with the selected identities.", type = "error")
          return()
        }
      }

      subset_seurat(subset_seurat_temp)
    })


    # Sous-ensemble basé sur l'expression génique (au moins un gène exprimé)
    observeEvent(input$apply_gene_subset, {
      tryCatch({
        req(single_dataset_object())
        gene_list <- unlist(strsplit(input$gene_list, ","))
        gene_list <- trimws(gene_list)

        if (length(gene_list) == 0) {
          showNotification("Please enter valid gene names.", type = "error")
          return()
        }

        # Créer une matrice logique pour l'expression génique
        expression_matrix <- sapply(gene_list, function(gene) {
          FetchData(single_dataset_object(), vars = gene) >= input$expression_threshold
        })

        # Sélection des cellules exprimant au moins un des gènes spécifiés
        cells_to_keep <- which(rowSums(expression_matrix) > 0)

        # Vérifier si des cellules sont sélectionnées
        if (length(cells_to_keep) == 0) {
          showNotification("No cells found with the specified gene expression criteria.", type = "error")
          return()
        }

        subset_seurat_temp <- subset(single_dataset_object(), cells = cells_to_keep)
        subset_seurat(subset_seurat_temp)
        showNotification(paste("Subsetting applied. Number of cells retained:", length(cells_to_keep)), type = "message")
      }, error = function(e) {
        showNotification(paste0("Error during gene-based subset: ", e$message), type = "error")
      })
    })




    # Downloading Seurat subset object
    output$download_subset_seurat <- downloadHandler(
      filename = function() {
        paste("seurat_subset_", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        # Display the modal box
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the seurat object for download...",
          easyClose = FALSE,
          footer = NULL
        ))

        tryCatch({
          req(subset_seurat())
          saveRDS(subset_seurat(), file)

          # Remove the modal box after saving the Seurat subset
          removeModal()

        }, error = function(e) {
          removeModal()
          showNotification(paste0("Error while downloading Seurat subset: ", e$message), type = "error")
        })
      }
    )


    # Umap plot with all clusters
    reactivePlotAll <- reactive({
      req(single_dataset_object())
      plot_data_all <- DimPlot(single_dataset_object(), group.by = "ident", label = TRUE) +
        theme(axis.line = element_line(size = 0.5)) +
        NoLegend()
      return(plot_data_all)
    })

    # Umap plot with filtered data
    reactivePlotSubset <- reactive({
      req(subset_seurat())
      plot_data_subset <- DimPlot(subset_seurat(), group.by = "ident", label = TRUE) +
        theme(axis.line = element_line(size = 0.5)) +
        NoLegend()
      return(plot_data_subset)
    })

    # Rendering Umap plot with all clusters
    output$global_umap <- renderPlot({
      plot_data_all <- reactivePlotAll()
      print(plot_data_all)
    })

    # Rendering Umap plot with filtered data
    output$subset_umap <- renderPlot({
      plot_data_subset <- reactivePlotSubset()
      print(plot_data_subset)
    })





}
