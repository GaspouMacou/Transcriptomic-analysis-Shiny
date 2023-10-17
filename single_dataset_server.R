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
library(data.table)
library(plotly)
library(colourpicker)
library(httr)
library(jsonlite)
library(monocle3)
library(SingleCellExperiment)
library(SeuratWrappers)

# single_dataset_server.R

single_dataset_server <-  function(input, output, session) {
  shinyjs::useShinyjs() #button deactivation


    # Create reactives objects for the single dataset part
    normalized_seuratRNA <- reactiveVal(NULL)
    scaled_seuratRNA <- reactiveVal(NULL)
    seuratRNA <- reactiveVal(NULL)
    variable_features <- reactiveValues(var_features = NULL)
    cluster_mapping <- reactiveValues()
    seurat_tmp <- reactiveVal(NULL)
    renamed_seurat <- reactiveVal()
    rv <- reactiveValues()
    gene_tables <- reactiveValues()
    integrated_seurat <- reactiveVal(NULL)
    gene_list <- reactiveValues(features = NULL)
    seurat_integrated_tmp <- reactiveVal(NULL)
    seurat_integrated <- reactiveVal(NULL)

    observe({
      if (!is.null(scaled_seuratRNA())) {
        seurat_tmp(scaled_seuratRNA())
      }
    })

    shinyjs::disable("pca_button")
    shinyjs::disable("scale_button")

    # Onglet 1 : Loading Data

    # Loading raw data
    observeEvent(input$file, {
      showNotification("Uploading and processing data...", type = "message")

      # Définissez le motif mitochondrial basé sur le choix de l'espèce.
      mt_pattern <- ifelse(input$species_choice == "mouse", "^mt-", "^MT-")

      if (dir.exists("unzipped")) {
        unlink("unzipped", recursive = TRUE)
      }
      zipPath <- input$file$datapath
      dataDir <- "unzipped"
      unzip(zipPath, exdir = dataDir)
      seuratRNA.data <- Read10X(dataDir)

      if (input$dataset_type == "snRNA") {
        seuratObj <- CreateSeuratObject(counts = seuratRNA.data, project = "seuratRNA", min.cells = 3, min.features = 200)
        seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
        seuratRNA(seuratObj)
        showNotification("snRNA-seq data processed successfully.", type = "message")
      } else if (input$dataset_type == "multiome") {
        rna.data <- seuratRNA.data$`Gene Expression`
        atac.data <- seuratRNA.data$`Peaks`
        rna.seurat <- CreateSeuratObject(counts = rna.data, project = "RNA", min.cells = 3, min.features = 200)
        atac.seurat <- CreateSeuratObject(counts = atac.data, project = "ATAC", min.cells = 3, min.features = 200)
        rna.seurat[["percent.mt"]] <- PercentageFeatureSet(rna.seurat, pattern = mt_pattern)
        seuratRNA(rna.seurat)
        showNotification("Multiome data processed successfully.", type = "message")
      }
    })

    # Loading seurat object
    observeEvent(input$load_seurat_file, {
      message("Attempting to read file at: ", input$load_seurat_file$datapath)
      tryCatch({
        loaded_seurat <- readRDS(input$load_seurat_file$datapath)
        message("File successfully read.")
        seuratRNA(loaded_seurat)
        scaled_seuratRNA(loaded_seurat)
        showNotification("L'objet Seurat a été chargé avec succès!")
      }, error = function(e) {
        showNotification(paste("An error occurred: ", e), type = "error")
      })
    })


    # 2ème onglet :
    # Display QC metrics on a VlnPlot chart
    output$vlnplot <- renderPlot({
      if (input$QCmetrics) {
        if (!is.null(normalized_seuratRNA())) {
          VlnPlot(normalized_seuratRNA(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        } else {
          VlnPlot(seuratRNA(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        }
      }
    })

    # Display the scatter plot 1 only when the user clicks on the button
    output$scatter_plot1 <- renderPlot({
      if (input$show_plots) {
        if (!is.null(normalized_seuratRNA())) {
          seuratRNA_subset <- subset(normalized_seuratRNA(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
          FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
        } else {
          seuratRNA_subset <- subset(seuratRNA(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
          FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
        }
      }
    })

    # Display the scatter plot 2 only when the user clicks on the button
    output$scatter_plot2 <- renderPlot({
      if (input$show_plots) {
        if (!is.null(normalized_seuratRNA())) {
          seuratRNA_subset <- subset(normalized_seuratRNA(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
          FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        } else {
          seuratRNA_subset <- subset(seuratRNA(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
          FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        }
      }
    })



    # Function to normalize data
    observeEvent(input$normalize_data, {
      # Désactiver le scale pendant la normalisation

      withProgress(message = "Normalizing data...", value = 0, {
        normalized_seuratRNA(NormalizeData(seuratRNA(), normalization.method = "LogNormalize", scale.factor = input$scale_factor))
        incProgress(1)
        print(normalized_seuratRNA)

        # Activez le scale une fois la normalisation terminée
        shinyjs::enable("scale_button")
      })
    })


    # Display the plot of variable features with the 10 most variable genes for normalized data
    output$variable_feature_plot <- renderPlot({
      if (input$show_plots && !is.null(normalized_seuratRNA())) {
        seuratRNA_norm <- FindVariableFeatures(normalized_seuratRNA(), selection.method = "vst", nfeatures = 2000)
        variable_features$var_features <- VariableFeatures(seuratRNA_norm)
        top10 <- head(VariableFeatures(seuratRNA_norm), 10)
        plot1 <- VariableFeaturePlot(seuratRNA_norm)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
        plot2
      }
    })

    #  Third tab: Scaling and PCA

    shinyjs::disable("pca_button")

    # Scale function:
    observeEvent(input$scale_button, {
      showNotification("Starting data scaling...", type = "message")
      all_genes <- rownames(normalized_seuratRNA())
      withProgress(message = "Scaling in progress...", value = 0, {
        tryCatch({
          scaled_seuratRNA(ScaleData(normalized_seuratRNA(), features = all_genes))
          shinyjs::enable("pca_button")  # Activate the PCA button only if scaling succeeds
          showNotification("Data scaling completed successfully.", type = "message")
        }, error = function(e) {
          shinyjs::disable("pca_button")
        })

      })
    })

    # PCA calculation:
    observeEvent(input$pca_button, {
      req(!is.null(scaled_seuratRNA()))
      runPCA()
      showNotification("PCA analysis completed successfully.", type = "message")
    }, ignoreNULL = FALSE, ignoreInit = TRUE)


    runPCA <- function() {
      req(!is.null(scaled_seuratRNA()))

      # Utilisez les features variables stockées, sinon retournez une erreur.
      if (is.null(variable_features$var_features)) {
        shinyjs::disable("pca_button")
        return(NULL)
      }

      tryCatch({
        scaled_seuratRNA(RunPCA(scaled_seuratRNA(), features = variable_features$var_features))
      }, error = function(e) {
        shinyjs::disable("pca_button")
      })
    }

    #Printing PCA results
    output$pca_results <- renderPrint({
      req(!is.null(scaled_seuratRNA()))

      ("pca" %in% names(scaled_seuratRNA()))
      print(scaled_seuratRNA()$pca, dims = 1:5, nfeatures = 5)

    })

    #Vizdimloading plot
    output$loading_plot <- renderPlot({
      tryCatch(VizDimLoadings(scaled_seuratRNA(), dims = 1:2, reduction = "pca"), error = function(e){})
    })

    # Dimplot
    output$dim_plot <- renderPlot({
      tryCatch(DimPlot(scaled_seuratRNA(), reduction = "pca"), error = function(e){})
    })



    # Fourth tab: Heatmaps
    observeEvent(input$runHeatmap, {
      output$heatmap <- renderPlot({
        req(!is.null(scaled_seuratRNA()))
        DimHeatmap(scaled_seuratRNA(), dims = input$dims, cells = input$cells, balanced = TRUE)
      })
    })


    # 5th tab: elbowplot an jacktraw

    observeEvent(input$run_jackstraw, {
      showNotification("Starting JackStraw analysis...", type = "message", duration = 3)
      withProgress(message = "Running JackStraw analysis...", value = 0, {
        seurat_tmp <- scaled_seuratRNA()
        seurat_tmp <- JackStraw(seurat_tmp, num.replicate = 100)
        seurat_tmp <- ScoreJackStraw(seurat_tmp, dims = 1:input$num_dims)  # Use the number of dimensions chosen by the user
        scaled_seuratRNA(seurat_tmp)
        incProgress(1)
      })

      # Display the JackStraw plot
      output$jackstraw_plot <- renderPlot({
        req(scaled_seuratRNA())
        JackStrawPlot(scaled_seuratRNA(), dims = 1:input$num_dims) # Use the number of dimensions chosen by the user
      })
      showNotification("JackStraw analysis completed.", type = "message")
    })

    # Show the ElbowPlot
    observeEvent(input$run_elbow, {
      output$elbow_plot <- renderPlot({
        req(scaled_seuratRNA())
        ElbowPlot(scaled_seuratRNA(), ndims = 50)
      })
    })

    # 6th tab: Neighbors calculation and clustering

    #Finding neighbors
    observeEvent(input$run_neighbors, {
      req(!is.null(scaled_seuratRNA()))
      showNotification("Looking for neighbors", type = "message", duration = 5)
      seurat_tmp <- scaled_seuratRNA()
      seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:input$dimension_1)
      seurat_tmp <- RunUMAP(seurat_tmp, dims = 1:input$dimension_1)
      scaled_seuratRNA(seurat_tmp)
      incProgress(1)
      print("Neighbors found and UMAP calculated.")  # Debug line
      output$neighbors_plot <- renderPlot({
        DimPlot(scaled_seuratRNA(), group.by = "ident")
      })
    })

    #Clustering
    observeEvent(input$run_clustering, {
      req(!is.null(scaled_seuratRNA()))
      showNotification("Clustering process started...", type = "message", duration = 5)
      seurat_tmp <- scaled_seuratRNA()
      seurat_tmp <- FindClusters(seurat_tmp, resolution = input$resolution_step1)
      scaled_seuratRNA(seurat_tmp)
      print("Clustering performed.")  # Debug line
    })


    # 7th tab: Calculate differentially expressed genes for each cluster


    # Variables réactives pour cet onglet
    gene_tables <- reactiveValues()
    all_markers <- reactiveVal()



    # Fonction pour calculer tous les marqueurs différentiels une seule fois
    calculate_all_markers <- function() {

      # Vérification si scaled_seuratRNA() n'est pas NULL
      if (is.null(scaled_seuratRNA())) {
        return(NULL)
      }

      markers <- FindAllMarkers(scaled_seuratRNA(), min.pct = 0.25)
      markers <- as.data.frame(markers)
      markers$gene <- paste0('<span class="gene-name">', markers$gene, '</span>')
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 2)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 2)
      all_markers(markers)
    }

    # Fonction pour mettre à jour les tableaux de gènes en fonction du nombre choisi
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

    # Observeur pour le bouton run_DE pour calculer tous les marqueurs
    observeEvent(input$run_DE, {
      calculate_all_markers()
      update_gene_tables_display()
    })

    # Observeur pour la mise à jour de l'affichage uniquement quand le nombre de gènes change
    observeEvent(input$number_genes, {
      if (!is.null(input$number_genes) && !is.null(all_markers())) {
        update_gene_tables_display()
      }
    })

    # Ce bloc s'exécute pour mettre à jour le tableau de sortie
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
        return("seurat_object.rds")
      },
      content = function(file) {
        saveRDS(scaled_seuratRNA(), file)
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
      req(!is.null(seurat_tmp()))
      gene_list <- strsplit(input$gene_list_vln, ",")[[1]]
      gene_list <- trimws(gene_list)
      print(paste("gene_list:", paste(gene_list, collapse = ", ")))
      withProgress(message = "Generating VlnPlot...", value = 0, {
        output$vln_plot <- renderPlot({
          VlnPlot(seurat_tmp(), features = gene_list, log = TRUE)
        })
        incProgress(1)
      })
    })

    observe({
      req(!is.null(seurat_tmp()))
      unique_idents <- unique(Idents(seurat_tmp()))
      updateSelectInput(session, "ident_1", choices = unique_idents)
      updateCheckboxGroupInput(session, "ident_2", choices = unique_idents, selected = "")
    })

    #  #FeaturePlot plot for each gene selected
    observeEvent(input$show_feature, {
      req(input$gene_list_feature, seurat_tmp())
      idents <- strsplit(input$gene_list_feature, ",")[[1]]
      idents <- trimws(idents)
      if(all(idents %in% rownames(seurat_tmp()))) {
        output$feature_plot <- renderPlot({
          FeaturePlot(seurat_tmp(), features = idents)
        })
      } else {
        print("The requested genes are not present in the dataset.")
      }
    })

    #Gene selector
    observeEvent(seuratRNA(), {
      gene_list <- c("", rownames(seuratRNA()@assays$RNA@counts))
      updateSelectInput(session,
                        "gene_select",
                        "Select a gene:",
                        choices = gene_list,
                        selected = "")
    })

    observeEvent(input$gene_select, {
      # Get the selected gene
      selected_gene <- input$gene_select

      if (input$gene_list_vln == "") {
        updateTextInput(session, "gene_list_vln", value = selected_gene)
      } else {
        updateTextInput(session, "gene_list_vln", value = paste(input$gene_list_vln, selected_gene, sep = ", "))
      }

      if (input$gene_list_feature == "") {
        updateTextInput(session, "gene_list_feature", value = selected_gene)
      } else {
        updateTextInput(session, "gene_list_feature", value = paste(input$gene_list_feature, selected_gene, sep = ", "))
      }
    })


    # Tab 9: Final UMAP

    #Reactive variable for that tab
    cluster_colours <- reactiveVal()  # Initialize a list to store cluster colors


    observe({
      if (!is.null(seurat_tmp())) {
        renamed_seurat(seurat_tmp())
      }
    })

    observe({
      if (!is.null(renamed_seurat())) {
        updateSelectInput(session, "select_cluster", choices = unique(Idents(renamed_seurat())))
      }
    })

    #function for renaming each cluster
    observeEvent(input$rename_single_cluster_button, {
      req(input$select_cluster, input$rename_single_cluster, renamed_seurat())

      updated_seurat <- renamed_seurat()

      if (input$rename_single_cluster %in% unique(Idents(updated_seurat))) {
        showNotification("This cluster name already exists!", type = "error")
        return()
      }

      Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster)) <- input$rename_single_cluster
      print(unique(Idents(updated_seurat)))

      renamed_seurat(updated_seurat)
      print(unique(Idents(updated_seurat)))

    })


    observe({
      req(renamed_seurat())
      message("Mise à jour des options de sélection du cluster")

      updateSelectInput(session, "select_cluster", choices = unique(Idents(renamed_seurat())))
      updateSelectInput(session, "cluster_select", choices = unique(Idents(renamed_seurat())))
    })

    observe({
      if (!is.null(renamed_seurat())) {
        message("Mise à jour des couleurs des clusters")

        current_colours <- scales::hue_pal()(length(unique(Idents(renamed_seurat()))))
        names(current_colours) <- sort(unique(Idents(renamed_seurat())))
        cluster_colours(current_colours)
      }
    })

    observeEvent(input$update_colour, {
      message("Mise à jour de la couleur du cluster sélectionné")

      current_colours <- cluster_colours()
      current_colours[input$cluster_select] <- input$cluster_colour
      cluster_colours(current_colours)
    })

    output$umap_finale <- renderPlotly({
      req(renamed_seurat())
      message("Génération du UMAP final")


      plot_data <- DimPlot(renamed_seurat(), group.by = "ident", pt.size = input$pt_size, label = TRUE, label.size = input$label_font_size) +
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

    # Fonction pour calculer tous les marqueurs différentiels pour le cluster choisi
    calculate_markers_for_cluster <- function() {
      req(renamed_seurat())
      markers <- FindMarkers(renamed_seurat(), ident.1 = input$cluster, min.pct = 0.25)
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 2)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 2)

      # Convertir les noms des gènes en liens
      markers$gene <- paste0('<span class="gene-name">', rownames(markers), '</span>')

      all_markers_10(markers)
    }


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

    # Observeur pour le bouton find_markers
    observeEvent(input$find_markers, {
      calculate_markers_for_cluster()
      update_gene_tables_display_10()
    })

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

    # Mettre à jour les choix pour le selectInput cluster
    observe({
      req(renamed_seurat())
      updateSelectInput(session, "cluster", choices = levels(renamed_seurat()))
    })

    observe({
      req(renamed_seurat())
      updateSelectInput(session, "cluster1", choices = levels(renamed_seurat()))
    })

    # Cluster selection
    observeEvent(input$cluster1, {
      updateSelectInput(session, "cluster2", choices = setdiff(levels(renamed_seurat()), input$cluster1))
    })

    # Variables réactives pour cet onglet
    gene_tables_new <- reactiveValues()
    all_markers_new <- reactiveVal()

    # Fonction pour calculer tous les marqueurs différentiels pour la comparaison de clusters
    calculate_markers_for_comparison <- function() {
      req(renamed_seurat())
      markers <- FindMarkers(renamed_seurat(), ident.1 = input$cluster1, ident.2 = input$cluster2, min.pct = 0.25)
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 2)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 2)

      # Convertir les noms des gènes en liens
      markers$gene <- paste0('<span class="gene-name">', rownames(markers), '</span>')

      all_markers_new(markers)
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

    # Observeur pour le bouton compare_markers
    observeEvent(input$compare_markers, {
      calculate_markers_for_comparison()
      update_gene_tables_display_new()
    })

    # Observeur pour la mise à jour de l'affichage uniquement quand le nombre de gènes change
    observeEvent(input$n_diff_markers, {
      if (!is.null(input$n_diff_markers) && !is.null(all_markers_new())) {
        update_gene_tables_display_new()
      }
    })

    # Ce bloc s'exécute pour mettre à jour le tableau de sortie pour cet onglet
    observe({
      output$table_new <- renderDataTable({
        datatable(gene_tables_new$table, escape = FALSE)
      })
    })



    #Download data for a cluster
    output$download_markers_single_cluster <- downloadHandler(
      filename = function() {
        paste("DE_markers_", input$cluster, "_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(rv$markers)
        write.csv(rv$markers, file)
      },
      contentType = "text/csv"
    )

    #Download multi-cluster comparison data
    output$download_markers_multiple_clusters <- downloadHandler(
      filename = function() {
        paste("diff_markers_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(rv$diff_markers, file)
      }
    )

    # Onglet 11 : Subseting

    #Variables for that tab
    subset_seurat <- reactiveVal()

    observe({
      if (!is.null(seurat_tmp())) {
        if (is.null(renamed_seurat())) { # Ajout de cette condition pour éviter d'écraser des renommages précédents
          renamed_seurat(seurat_tmp())
        }
        subset_seurat(renamed_seurat()) # Initialize subset_seurat with renamed_seurat
      }
    })

    observe({
      if (!is.null(renamed_seurat())) {
        updateSelectInput(session, "select_ident_subset", choices = unique(Idents(renamed_seurat())))
      }
    })

    #Run subset
    observeEvent(input$apply_subset, {
      req(renamed_seurat())
      subset_seurat_temp <- renamed_seurat()
      if (length(input$select_ident_subset) > 0) {
        subset_seurat_temp <- subset(x = subset_seurat_temp, idents = input$select_ident_subset)
      }
      subset_seurat(subset_seurat_temp)
    })

    #Downloading Seurat subset object
    output$download_subset_seurat <- downloadHandler(
      filename = function() {
        paste("seurat_subset_", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        req(subset_seurat())
        saveRDS(subset_seurat(), file)
      }
    )

    # Reactive to create a UMAP plot
    reactivePlot <- reactive({
      req(subset_seurat())
      plot_data <- DimPlot(subset_seurat(), group.by = "ident", label = TRUE) +
        theme(axis.line = element_line(size = 0.5)) +
        NoLegend()
      return(plot_data)
    })

    # Rendering the UMAP plot
    output$final_umap_subset <- renderPlot({
      plot_data <- reactivePlot()
      print(plot_data)  # Assurez-vous que plot_data n'est pas NULL avant de l'imprimer
    })


}
