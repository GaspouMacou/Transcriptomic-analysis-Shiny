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

# multiple_datasets_server.R

multiple_datasets_server <-  function(input, output, session) {




    #Variables for the merge dataset part
    seurat_objects <- reactiveValues()
    data_loaded <- reactiveValues()
    rv_merge <- reactiveValues(seurat_integrated = NULL, added_metadata = character(0))
    renamed_seurat_merge <- reactiveVal(NULL)
    merged_gene_tables <- reactiveValues()
    renamed_datasets <- NULL
    rv_metadata <- reactiveValues(num_fields = 1)

    # Créez une liste réactive pour stocker les noms des champs de métadonnées
    reactive_metadata_fields <- reactive({
      req(rv_merge$seurat_integrated)
      # Obtenez toutes les colonnes de métadonnées
      all_metadata_fields <- colnames(rv_merge$seurat_integrated@meta.data)
      # Liste des champs à exclure
      exclude_fields <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "integrated_snn_res.0.5")
      # Obtenez la liste des champs à inclure en excluant les champs spécifiés de la liste complète
      include_fields <- setdiff(all_metadata_fields, exclude_fields)
      return(include_fields)
    })


    shinyjs::disable("add_field")
    shinyjs::disable("add_metadata")

    observeEvent(input$load_seurat_file_merge, {
      message("Attempting to read file at: ", input$load_seurat_file_merge$datapath)
      tryCatch({
        loaded_seurat <- readRDS(input$load_seurat_file_merge$datapath)
        message("File successfully read.")
        rv_merge$seurat_integrated <- loaded_seurat
        showNotification("The Seurat object has been successfully loaded!")

        # Activer les éléments après le chargement réussi
        shinyjs::enable("add_field")
        shinyjs::enable("add_metadata")
        shinyjs::enable("runScalePCA")
      }, error = function(e) {
        showNotification(paste("An error occurred: ", e), type = "error")
      })
    })


    #Processing of loaded datas
    observe_file_input <- function(index) {
      observeEvent(input[[paste0("merge", index)]], {
        if (is.null(input[[paste0("merge", index)]])) {
          return(NULL)
        }

        if (!is.null(data_loaded[[paste0("loaded", index)]]) && data_loaded[[paste0("loaded", index)]] == TRUE) {
          return(NULL)
        }

        seurat_object <- NULL
        mt_pattern_merge <- ifelse(input$species_choice_merge == "mouse", "^mt-", "^MT-")
        dataset_type <- input[[paste0("dataset_type_merge", index)]]

        if (is.null(dataset_type)) {
          return(NULL)
        }

        if(dataset_type == "seurat_object_merge") {
          seurat_object <- readRDS(input[[paste0("merge", index)]]$datapath)
          message("Chargement de l'objet Seurat terminé")

        }
        else {
          if (dir.exists(paste0("unzipped", index))) {
            unlink(paste0("unzipped", index), recursive = TRUE)
          }
          unzip(input[[paste0("merge", index)]]$datapath, exdir = paste0("unzipped", index))
          message("Décompression terminée")
          data <- Read10X(paste0("unzipped", index))
          message("Lecture des données terminée")

          if (dataset_type == "snRNA_merge") {
            seurat_object <- CreateSeuratObject(counts = data, project = input[[paste0("dataset_name", index)]], min.cells = 3, min.features = 200)
            seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,pattern = mt_pattern_merge)
            message("Création de l'objet Seurat pour snRNA-seq terminée")
          }
          else if (dataset_type == "multiome_merge") {
            rna.data <- data$`Gene Expression`
            atac.data <- data$`Peaks`

            rna.seurat <- CreateSeuratObject(counts = rna.data, project = input[[paste0("dataset_name", index)]], min.cells = 3, min.features = 200)
            atac.seurat <- CreateSeuratObject(counts = atac.data, project = "ATAC", min.cells = 3, min.features = 200)

            rna.seurat[["percent.mt"]] <- PercentageFeatureSet(rna.seurat, pattern = mt_pattern_merge)
            seurat_object <- rna.seurat
            message("Création de l'objet Seurat pour Multiome terminée")
          }

          seurat_object <- NormalizeData(seurat_object)
          seurat_object <- FindVariableFeatures(seurat_object, nfeatures = 4000)
          seurat_object <- ScaleData(seurat_object)
          seurat_object <- RunPCA(seurat_object)
          message("Traitement de l'objet Seurat terminé")
        }

        seurat_object@meta.data$dataset <- input[[paste0("dataset_name", index)]]
        seurat_object$orig.ident <- seurat_object@meta.data$dataset
        seurat_objects[[paste0("seurat_object", index)]] <- seurat_object
        data_loaded[[paste0("loaded", index)]] <- TRUE
        showNotification("Data processed successfully!", type = "message")

      })
    }

    output$fileInputs <- renderUI({
      lapply(1:input$num_datasets, function(i) {
        fluidRow(
          column(3, textInput(paste0('dataset_name', i), paste0('Nom du Dataset ', i), value = paste0("Dataset ", i))),
          column(3, selectInput(paste0("dataset_type_merge", i), "Type de données :",
                                choices = list("snRNA-seq" = "snRNA_merge",
                                               "Multiome" = "multiome_merge",
                                               "Objet Seurat" = "seurat_object_merge"))),
          column(6, fileInput(paste0('merge', i), paste0('Choisir un fichier Seurat ou .gz ', i), accept=c('rds', 'gz')))
        )
      })
    })



    # Observe the change in the number of datasets and create the necessary entries
    observeEvent(input$num_datasets, {
      lapply(1:input$num_datasets, function(i) {
        observe_file_input(i)
      })
    })

    # Count the number of datasets loaded
    observe({
      if(is.numeric(input$num_datasets) && !is.null(input$num_datasets)) {
        num_loaded_datasets <- sum(sapply(1:input$num_datasets, function(i) {
          !is.null(data_loaded[[paste0("loaded", i)]]) && data_loaded[[paste0("loaded", i)]] == TRUE
        }))
        if (num_loaded_datasets == input$num_datasets) {
          shinyjs::enable("integrate")
        } else {
          shinyjs::disable("integrate")
        }
      } else {
        shinyjs::disable("integrate")
        warning("input$num_datasets is not valid: ", input$num_datasets)
      }
    })


    # Integration function
    integrate_data <- function(seurat_list) {
      print("Starting integration")

      # Sélection des features pour l'intégration
      print("Selecting integration features")
      features <- SelectIntegrationFeatures(
        object.list = seurat_list,
        nfeatures = 2000 # Augmentez cette valeur pour inclure plus de gènes
      )

      # Trouver les ancres d'intégration
      print("Finding Integration Anchors")
      anchors <- FindIntegrationAnchors(
        object.list = seurat_list,
        dims = 1:30
      )

      # Intégration des données
      print("Integrating Data")
      seurat_integrated_temp <- IntegrateData(
        anchorset = anchors,
        dims = 1:30,
        features.to.integrate = features # Utilisez les features sélectionnés ici
      )
      DefaultAssay(seurat_integrated_temp) <- "integrated"
      new_metadata <- seurat_integrated_temp@meta.data
      colnames(new_metadata)[colnames(new_metadata) == "orig.ident"] <- "dataset"
      seurat_integrated_temp <- AddMetaData(seurat_integrated_temp, metadata = new_metadata)

      print("Finished integration")
      print(paste("Integrated object class: ", class(seurat_integrated_temp)))
      print(paste("Integrated object dimensions: ", dim(seurat_integrated_temp)))

      return(seurat_integrated_temp)
    }



    # Integration function
    observeEvent(input$integrate, {
      print("Integrate button pressed")
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_list[[i]] <- seurat_objects[[paste0("seurat_object", i)]]
      }
      print("Copied Seurat objects from reactiveValues to a standard list")
      print(paste("seurat_list length:", length(seurat_list)))
      if (length(seurat_list) < 2) {
        stop("Please upload at least two datasets for integration")
      }
      rv_merge$seurat_integrated <- integrate_data(seurat_list)
      print("Finished data integration")
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      shinyjs::enable("runScalePCA")

    })



    #Add metadata field
    observeEvent(input$add_field, {
      rv_metadata$num_fields <- rv_metadata$num_fields + 1
    })

    output$metadata_inputs <- renderUI({
      req(rv_merge$seurat_integrated)
      datasets <- unique(rv_merge$seurat_integrated@meta.data$dataset)

      lapply(1:rv_metadata$num_fields, function(j) {
        fluidRow(
          column(6, textInput(paste0("metadata_name_", j), paste0("Metadata Field Name ", j), value = "")),
          column(6, lapply(datasets, function(dataset_name) {
            textInput(paste0("metadata_value_", dataset_name, "_", j), paste0("Value for ", dataset_name, " (", input[[paste0("dataset_name", which(datasets == dataset_name))]], ") ", j), value = "")
          }))
        )
      })
    })

    # Add metadata by the user
    observeEvent(input$add_metadata, {
      req(rv_merge$seurat_integrated)
      datasets <- unique(rv_merge$seurat_integrated@meta.data$dataset)
      for (j in 1:rv_metadata$num_fields) {
        metadata_field_name <- input[[paste0("metadata_name_", j)]]
        if (metadata_field_name == "") {
          next
        }
        if (!metadata_field_name %in% colnames(rv_merge$seurat_integrated@meta.data)) {
          rv_merge$seurat_integrated@meta.data[[metadata_field_name]] <- NA
          rv_merge$added_metadata <- c(rv_merge$added_metadata, metadata_field_name)
        }
        for (dataset in datasets) {
          rows <- which(!is.na(rv_merge$seurat_integrated@meta.data$dataset) & rv_merge$seurat_integrated@meta.data$dataset == dataset)
          metadata_field_value <- input[[paste0("metadata_value_", dataset, "_", j)]]
          rv_merge$seurat_integrated@meta.data[[metadata_field_name]][rows] <- metadata_field_value
        }
      }
      print(head(rv_merge$seurat_integrated@meta.data))
      new_metadata_fields <- colnames(rv_merge$seurat_integrated@meta.data)
      reactive_metadata_fields <- reactive({ new_metadata_fields })
    })

    # Identifiez et extrayez le premier objet Seurat chargé
    # pour obtenir la liste de gènes de l'assay RNA
    reactive_gene_list_merge <- reactive({
      # Vérifiez si l'objet intégré est disponible
      req(rv_merge$seurat_integrated)

      # Obtenez le premier dataset à partir des métadonnées
      first_dataset <- unique(rv_merge$seurat_integrated@meta.data$dataset)[1]

      # Obtenez les cellules de ce premier dataset
      cells_in_first_dataset <- which(rv_merge$seurat_integrated@meta.data$dataset == first_dataset)

      # Obtenez la liste unique de gènes à partir de l'assay RNA de ces cellules
      unique_genes <- unique(rownames(rv_merge$seurat_integrated[cells_in_first_dataset, ]@assays$RNA@counts))

      # Ajoutez une option vide au début de la liste de gènes
      c("", unique_genes)
    })


    # Tab 13: Scaling and PCA reduction

    shinyjs::disable("runScalePCA")
    shinyjs::disable("runFindNeighbors")
    shinyjs::disable("runFindClusters")


    # Scaling, PCA and display elbow plot
    observeEvent(input$runScalePCA, {
      req(rv_merge$seurat_integrated)
      rv_merge$seurat_integrated <- scale_PCA(rv_merge$seurat_integrated, input)
      output$elbow_plot2 <- renderPlot({
        ElbowPlot(rv_merge$seurat_integrated)
      })
    })

    scale_PCA <- function(seurat_integrated, input) {
      showNotification("Scaling and PCA started...", type = "message")
      if (is.null(input$resolution_step2) || is.null(input$dimension_2)) {
        stop("Missing input.")
      }
      all_genes <- rownames(seurat_integrated)
      seurat_integrated_temp <- seurat_integrated
      seurat_integrated_temp <- FindVariableFeatures(seurat_integrated_temp, selection.method = "vst", nfeatures = 3000)
      seurat_integrated_temp <- ScaleData(seurat_integrated_temp, features = all_genes)
      seurat_integrated_temp <- RunPCA(seurat_integrated_temp, npcs = 50)
      print(head(seurat_integrated_temp@meta.data))
      shinyjs::enable("runFindNeighbors")

      return(seurat_integrated_temp)


    }

    # For the "Find neighbors and display UMAP" button
    observeEvent(input$runFindNeighbors, {
      showNotification("Finding neighbors started...", type = "message")
      req(rv_merge$seurat_integrated)
      rv_merge$seurat_integrated <- find_neighbors(rv_merge$seurat_integrated, input)
      output$UMAPPlot <- renderPlot({
        UMAPPlot(rv_merge$seurat_integrated, group.by="orig.ident")
      })
    })
    find_neighbors <- function(seurat_integrated, input) {
      seurat_integrated_temp <- seurat_integrated
      seurat_integrated_temp <- FindNeighbors(seurat_integrated_temp, dims = 1:input$dimension_2)
      seurat_integrated_temp <- RunUMAP(seurat_integrated_temp, dims = 1:input$dimension_2)
      shinyjs::enable("runFindClusters")
      return(seurat_integrated_temp)
    }

    # Clustering
    observeEvent(input$runFindClusters, {
      req(rv_merge$seurat_integrated)
      rv_merge$seurat_integrated <- find_clusters(rv_merge$seurat_integrated, input)
      output$UMAPPlot <- renderPlot({
        p1 <- DimPlot(rv_merge$seurat_integrated, reduction = "umap", group.by = "orig.ident")
        p2 <- DimPlot(rv_merge$seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE)
        p1 + p2
      })
    })
    find_clusters <- function(seurat_integrated, input) {
      showNotification("Clustering started...", type = "message")
      seurat_integrated_temp <- seurat_integrated
      seurat_integrated_temp <- FindClusters(seurat_integrated_temp, resolution = input$resolution_step2)
      print(head(seurat_integrated_temp@meta.data))
      return(seurat_integrated_temp)
    }


    # Tab 14: Visualize gene expression

    observe({
      updateSelectInput(session, "group_by_select", choices = reactive_metadata_fields())
    })

    observe({
      updateSelectInput(session, "geneInput_merge", choices = reactive_gene_list_merge())
    })

    #Vlnplot of the selected gene
    observeEvent(input$runVlnPlot, {
      req(input$geneInput_merge, rv_merge$seurat_integrated)

      # Préparez la liste des gènes en se basant sur l'input de l'utilisateur
      gene_list <- strsplit(input$geneInput_merge, ",")[[1]]
      gene_list <- trimws(gene_list)

      # Récupérez le paramètre sélectionné pour le groupement
      selected_group_by <- input$group_by_select

      # Vérifiez que tous les gènes demandés sont présents dans l'assay "RNA"
      if(all(gene_list %in% rownames(rv_merge$seurat_integrated[["RNA"]]@counts))) {
        output$VlnPlot2 <- renderPlot({
          # Créez le VlnPlot en utilisant les paramètres sélectionnés
          VlnPlot(rv_merge$seurat_integrated, features = gene_list, group.by = selected_group_by, assay = "RNA")
        })
      } else {
        print("Requested genes are not present in the dataset.")
      }
    })

    # FeaturePlot du gène sélectionné
    observeEvent(input$runFeaturePlot, {
      req(input$geneInput_merge, rv_merge$seurat_integrated)

      gene_list <- strsplit(input$geneInput_merge, ",")[[1]]
      gene_list <- trimws(gene_list)

      # Vérifiez que tous les gènes demandés sont présents dans l'assay "RNA"
      if(all(gene_list %in% rownames(rv_merge$seurat_integrated[["RNA"]]@counts))) {
        # Définissez l'assay par défaut à "RNA"
        DefaultAssay(rv_merge$seurat_integrated) <- "RNA"

        output$FeaturePlot2 <- renderPlot({
          # Créez le FeaturePlot en utilisant les paramètres sélectionnés
          FeaturePlot(rv_merge$seurat_integrated, features = gene_list)
        })
      } else {
        print("Requested genes are not present in the dataset.")
      }
    })



    # Onglet 15: Calculer les gènes différentiellement exprimés pour chaque cluster dans l'ensemble de données fusionné

    # Variables pour l'onglet 15
    gene_tables_merge <- reactiveValues()
    all_markers_merge <- reactiveVal()

    shinyjs::disable("download_DE_merged")


    # Fonction pour calculer tous les marqueurs différentiels pour le dataset fusionné
    calculate_merged_markers <- function() {
      req(rv_merge$seurat_integrated)
      markers <- FindAllMarkers(rv_merge$seurat_integrated, min.pct = 0.25)
      markers <- as.data.frame(markers)
      markers$gene <- paste0('<span class="gene-name">', rownames(markers), '</span>')

      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 2)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 2)

      all_markers_merge(markers)
    }



    update_gene_tables_display_merge <- function() {
      markers <- all_markers_merge()
      num_genes_to_display <- input$number_genes_merge
      for (cluster in unique(markers$cluster)) {
        gene_tables_merge[[paste0("table_", cluster)]] <- head(markers[markers$cluster == cluster, ], n = num_genes_to_display)
      }

      output$diff_genes_tables_merge <- renderUI({
        tagList(
          lapply(names(gene_tables_merge), function(name) {
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


    # Observeur pour le bouton run_DE_merged
    observeEvent(input$run_DE_merged, {
      calculate_merged_markers()
      update_gene_tables_display_merge()
      shinyjs::enable("download_DE_merged")

    })



    # Observeur pour la mise à jour de l'affichage uniquement quand le nombre de gènes change
    observeEvent(input$number_genes_merge, {
      if (!is.null(input$number_genes_merge) && !is.null(all_markers_merge())) {
        update_gene_tables_display_merge()
      }
    })


    observe({
      lapply(names(gene_tables_merge), function(name) {
        output[[name]] <- renderDT({
          datatable(gene_tables_merge[[name]], escape = FALSE)
        })
      })
    })


    # Notification for users
    output$previous_tab_notification <- renderUI({
      if (!is.null(session$userData$previous_tab_notification_msg)) {
        list(
          tags$hr(),
          tags$strong(session$userData$previous_tab_notification_msg),
          tags$hr()
        )
      }
    })





    # Download seurat  object
    output$save_seurat_merge <- downloadHandler(
      filename = function() {
        return("seurat_object.rds")
      },
      content = function(file) {
        saveRDS(rv_merge$seurat_integrated, file)
      }
    )

    # Download markers as a CSV
    output$download_DE_merged <- downloadHandler(
      filename = function() {
        paste("DE_genes_merged_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        markers <- all_markers_merge()
        write.csv(markers, file)
      }
    )

    # Tab 16: Assigning cell identity merge

    # Reactive variable for that tab
    renamed_seurat_merge <- reactiveVal()
    cluster_colours_merge <- reactiveVal()

    observe({
      if (!is.null(rv_merge$seurat_integrated)) {
        renamed_seurat_merge(rv_merge$seurat_integrated)
      }
    })

    observe({
      if (!is.null(renamed_seurat_merge())) {
        updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(renamed_seurat_merge())))
      }
    })

    # Renaming each cluster
    observeEvent(input$rename_single_cluster_merge_button, {
      req(input$select_cluster_merge, input$rename_single_cluster_merge, renamed_seurat_merge())
      updated_seurat <- renamed_seurat_merge()
      if (input$rename_single_cluster_merge %in% unique(Idents(updated_seurat))) {
        showNotification("This cluster name already exists!", type = "error")
        return()
      }
      Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster_merge)) <- input$rename_single_cluster_merge
      renamed_seurat_merge(updated_seurat)
    })


    observe({
      if (!is.null(renamed_seurat_merge())) {
        updateSelectInput(session, "select_color_merge", choices = unique(Idents(renamed_seurat_merge())))
        current_colours <- scales::hue_pal()(length(unique(Idents(renamed_seurat_merge()))))
        names(current_colours) <- sort(unique(Idents(renamed_seurat_merge())))
        cluster_colours_merge(current_colours)
      }
    })

    # Color selection
    observeEvent(input$update_colour_merge_button, {
      current_colours <- cluster_colours_merge()
      current_colours[input$select_color_merge] <- input$select_cluster_merge_color
      cluster_colours_merge(current_colours)
    })

    # Finale UMAP
    output$umap_finale_merge <- renderPlotly({
      req(renamed_seurat_merge())
      plot_data <- DimPlot(renamed_seurat_merge(), group.by = "ident", pt.size = input$pt_size_merge, label = TRUE, label.size = input$label_font_size_merge) +
        theme(axis.line = element_line(size = 0.5)) +
        scale_color_manual(values = cluster_colours_merge()) +
        NoLegend()
      interactive_plot <- ggplotly(plot_data, tooltip = "text")
      interactive_plot <- interactive_plot %>%
        layout(
          title = list(text = input$plot_title_merge, font = list(size = 24)),
          hovermode = "closest"
        )
      return(interactive_plot)
    })




    # Tab 17: Calculation of differentially expressed genes

    # Reactive variable for that tab
    markers_df_reactive <- reactiveVal() # Create a new reactive value to store markers_df
    diff_genes_compare <- reactiveVal()   # Define a new reactive value to store comparison tables

    shinyjs::disable("download_markers_single_cluster_merge")
    shinyjs::disable("download_markers_multiple_clusters_merge")
    shinyjs::disable("download_diff_dataset_cluster")

    # Drop-down menu to filter by Dataset
    output$dataset_filter_ui <- renderUI({
      req(renamed_seurat_merge())
      selectInput("dataset_filter",
                  label = "Filtrer par Dataset",
                  choices = c( unique(renamed_seurat_merge()@meta.data$dataset)),
                  selected = unique(renamed_seurat_merge()@meta.data$dataset),
                  multiple = TRUE)
    })

    # Display of filtered UMAP graph
    output$filtered_umap_plot <- renderPlot({
      req(renamed_seurat_merge())
      req(input$dataset_filter)  # assurez-vous que dataset_filter est défini

      # Vérifiez si "dataset" est une colonne valide des métadonnées
      if (!"dataset" %in% colnames(renamed_seurat_merge()@meta.data)) {
        stop("The 'dataset' column was not found in the Seurat object's metadata.")
      }

      # Si "Tous les Datasets" est sélectionné, ou si la sélection est NULL, affichez tous les datasets
      if ("Tous les Datasets" %in% input$dataset_filter || is.null(input$dataset_filter)) {
        return(DimPlot(renamed_seurat_merge(), group.by = "ident", label = TRUE, label.size = 5) + NoLegend())
      } else {
        # Assurez-vous que les datasets choisis existent réellement dans les données
        valid_datasets <- input$dataset_filter %in% unique(renamed_seurat_merge()@meta.data$dataset)
        if (any(!valid_datasets)) {
          warning("Some selected datasets were not found in the data:",
                  paste(input$dataset_filter[!valid_datasets], collapse = ", "))
        }

        # Si aucune correspondance n'est trouvée, arrêtez de rendre le plot pour éviter les erreurs
        if (all(!valid_datasets)) {
          stop("None of the selected datasets were found in the data.")
        }

        # Effectuez la sous-ensemble et tracez
        subset_seurat <- subset(renamed_seurat_merge(), subset = dataset %in% input$dataset_filter[valid_datasets])
        return(DimPlot(subset_seurat, group.by = "ident", label = TRUE, label.size = 5) + NoLegend())
      }
    })


    # Reactively update the cluster drop-down menu
    observe({
      req(renamed_seurat_merge())
      cluster_choices <- unique(Idents(renamed_seurat_merge()))
      updateSelectInput(session, "selected_cluster", choices = cluster_choices, selected = cluster_choices[1])
    })

    # Calculation of differentially expressed genes
    observeEvent(input$calculate_DE, {
      showNotification("Finding differentially expressed features...", type = "message")
      req(renamed_seurat_merge())
      subset_seurat <- renamed_seurat_merge()
      if (input$dataset_filter != "Tous les Datasets") {
        subset_seurat <- subset(renamed_seurat_merge(), subset = dataset == input$dataset_filter)
      }
      markers <- FindMarkers(subset_seurat,
                             ident.1 = input$selected_cluster,
                             min.pct = 0.25)
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 2)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 2)
      markers_df <- as.data.frame(markers)
      markers_df$gene <- rownames(markers_df)
      markers_df_reactive(markers_df)
      shinyjs::enable("download_markers_single_cluster_merge")
      markers_df <- head(markers_df[order(markers_df$p_val),], n = input$num_genes)
      output$DE_genes_table <- renderTable({
        markers_df
      })
    })


    # Download gene table
    output$download_markers_single_cluster_merge <- downloadHandler(
      filename = function() {
        paste("genes-differents-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        markers_df <- markers_df_reactive()
        if (!is.null(markers_df)) {
          write.csv(markers_df, file)
        } else {
          stop("No data to download. Perform an analysis first.")
        }
      },
      contentType = "text/csv"
    )

    # Drop-down menus to select clusters for comparison
    output$cluster1_compare_ui <- renderUI({
      req(renamed_seurat_merge())
      selectInput("cluster1_compare",
                  label = "Selec first cluster",
                  choices = unique(Idents(renamed_seurat_merge())),
                  selected = NULL)
    })

    output$cluster2_compare_ui <- renderUI({
      req(renamed_seurat_merge())
      selectInput("cluster2_compare",
                  label = "Select second cluster",
                  choices = unique(Idents(renamed_seurat_merge())),
                  selected = NULL)
    })

    # Cluster comparison
    observeEvent(input$compare_clusters_button, {
      showNotification("Finding differentially expressed features...", type = "message")
      req(input$cluster1_compare, input$cluster2_compare, renamed_seurat_merge())
      if(input$cluster1_compare == input$cluster2_compare) {
        showNotification("Please select two different clusters for comparison!", type = "error")
        return()
      }
      temp_res <- FindMarkers(renamed_seurat_merge(), ident.1 = input$cluster1_compare, ident.2 = input$cluster2_compare) # Store results in a temporary variable
      temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 2)
      temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 2)
      temp_res$gene <- rownames(temp_res)
      shinyjs::enable("download_markers_multiple_clusters_merge")
      diff_genes_compare(temp_res)
    })

    output$diff_genes_table_compare <- renderTable({
      diff_genes_compare_df <- diff_genes_compare()
      req(!is.null(diff_genes_compare_df))
      head(diff_genes_compare_df[order(diff_genes_compare_df$p_val), ], n = input$num_genes_compare)
    })

    # Download gene comparison table
    output$download_markers_multiple_clusters_merge <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster1_compare, "-VS-", input$cluster2_compare, "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        diff_genes_data <- diff_genes_compare()
        req(!is.null(diff_genes_data))
        diff_genes_data_sorted <- diff_genes_data[order(diff_genes_data$p_val), ]
        write.csv(diff_genes_data_sorted, file)
      },
      contentType = "text/csv"
    )



    output$dataset1_compare_ui <- renderUI({
      req(renamed_seurat_merge())
      selectInput("dataset1_compare",
                  label = "Select the first dataset",
                  choices = unique(renamed_seurat_merge()$dataset),
                  selected = NULL)
    })

    output$dataset2_compare_ui <- renderUI({
      req(renamed_seurat_merge())
      selectInput("dataset2_compare",
                  label = "Select the second Dataset",
                  choices = unique(renamed_seurat_merge()$dataset),
                  selected = NULL)
    })

    output$cluster_compare_ui <- renderUI({
      req(renamed_seurat_merge())
      selectInput("cluster_compare",
                  label = "Select cluster",
                  choices = unique(Idents(renamed_seurat_merge())),
                  selected = NULL)
    })

    observeEvent(input$compare_datasets_button, {
      showNotification("Finding differentially expressed features...", type = "message")
      req(input$dataset1_compare, input$dataset2_compare, input$cluster_compare, renamed_seurat_merge())

      if(input$dataset1_compare == input$dataset2_compare) {
        showNotification("Please select two different datasets for comparison!", type = "error")
        return()
      }

      subset_seurat1 <- subset(renamed_seurat_merge(), subset = dataset == input$dataset1_compare & seurat_clusters == input$cluster_compare)
      subset_seurat2 <- subset(renamed_seurat_merge(), subset = dataset == input$dataset2_compare & seurat_clusters == input$cluster_compare)

      if(dim(subset_seurat1)[1] == 0 || dim(subset_seurat2)[1] == 0) {
        showNotification("The selected cluster is not present in any of the selected datasets. Please choose another cluster.", type = "error")
        return()
      }

      idents_subset1 <- unique(Idents(subset_seurat1))
      idents_subset2 <- unique(Idents(subset_seurat2))

      diff_genes <- FindMarkers(object = renamed_seurat_merge(),
                                ident.1 = idents_subset1[1],
                                ident.2 = idents_subset2[1])

      diff_genes$p_val <- formatC(diff_genes$p_val, format = "e", digits = 2)
      diff_genes$p_val_adj <- formatC(diff_genes$p_val_adj, format = "e", digits = 2)
      shinyjs::enable("download_diff_dataset_cluster")
      diff_genes_compare(diff_genes)
    })

    output$download_diff_dataset_cluster <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster_compare, "-DS1-", input$dataset1_compare, "-DS2-", input$dataset2_compare, "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        diff_genes_data <- diff_genes_compare()
        req(!is.null(diff_genes_data))

        diff_genes_data$p_val <- formatC(diff_genes_data$p_val, format = "e", digits = 2)
        diff_genes_data$p_val_adj <- formatC(diff_genes_data$p_val_adj, format = "e", digits = 2)

        diff_genes_data_sorted <- diff_genes_data[order(diff_genes_data$p_val), ]
        write.csv(diff_genes_data_sorted, file)
      },
      contentType = "text/csv"
    )



    # Tab 18: Subseting seurat object

    # Reactive variable for that tab
    subset_seurat_merge <- reactiveVal(NULL)
    shinyjs::disable("download_subset_merge")

    observe({
      if (!is.null(rv_merge$seurat_integrated)) {
        renamed_seurat_merge(rv_merge$seurat_integrated)
      }
    })

    observe({
      if (!is.null(renamed_seurat_merge())) {
        updateSelectInput(session, "select_ident_subset_merge", choices = unique(Idents(renamed_seurat_merge())))
      }
    })


    observe({
      if (!is.null(rv_merge$seurat_integrated)) {
        renamed_seurat_merge(rv_merge$seurat_integrated)
        subset_seurat_merge(rv_merge$seurat_integrated)
      }
    })

    # Apply the subset selection
    observeEvent(input$apply_subset_merge, {
      req(input$select_ident_subset_merge, renamed_seurat_merge())
      if (length(input$select_ident_subset_merge) > 0) {
        subset_seurat_merge(subset(x = renamed_seurat_merge(), idents = input$select_ident_subset_merge))
      } else {
        subset_seurat_merge(renamed_seurat_merge())
      }
      shinyjs::enable("download_subset_merge")

    })

    # Download seurat subset object
    output$download_subset_merge <- downloadHandler(
      filename = function() {
        paste("seurat_subset-", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        req(subset_seurat_merge())
        saveRDS(subset_seurat_merge(), file)
      }
    )

    # Umap subset
    output$final_umap_subset_merge <- renderPlot({
      req(subset_seurat_merge())
      plot_data <- DimPlot(subset_seurat_merge(), group.by = "ident", label = TRUE) +
        theme(axis.line = element_line(size = 0.5)) +
        NoLegend()
      print(plot_data)
    })


}
