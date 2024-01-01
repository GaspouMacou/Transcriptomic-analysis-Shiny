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


# multiple_datasets_server.R

multiple_datasets_server <-  function(input, output, session) {




    #Variables for the merge dataset part
    multiple_datasets_object <- reactiveVal ()
    seurat_objects <- reactiveValues()
    data_loaded <- reactiveValues()
    merged_gene_tables <- reactiveValues()
    rv_metadata <- reactiveValues(num_fields = 1)

    shinyjs::disable("add_field")
    shinyjs::disable("add_metadata")
    
    # Load a seurat object
    observeEvent(input$load_seurat_file_merge, {
      message("Attempting to read file at: ", input$load_seurat_file_merge$datapath)
      tryCatch({
        loaded_seurat <- readRDS(input$load_seurat_file_merge$datapath)
        message("File successfully read.")
        multiple_datasets_object(loaded_seurat) # Notez l'utilisation des parenthèses ici
        showNotification("The Seurat object has been successfully loaded!")
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
                                               "Seurat Object" = "seurat_object_merge"))),
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

    # Count the number of dataset loaded
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
      print("Selecting integration features")
      features <- SelectIntegrationFeatures(
        object.list = seurat_list,
        nfeatures = 2000 
      )
      print("Finding Integration Anchors")
      anchors <- FindIntegrationAnchors(
        object.list = seurat_list,
        dims = 1:30
      )
      print("Integrating Data")
      seurat_integrated_temp <- IntegrateData(
        anchorset = anchors,
        dims = 1:30,
        features.to.integrate = features
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



    # Integration button
    observeEvent(input$integrate, {
      tryCatch({
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
        
        integrated_object <- integrate_data(seurat_list)
        multiple_datasets_object(integrated_object)  # Mise à jour de l'objet intégré
        print("Finished data integration")
        
        # Vérifier que l'objet a été mis à jour
        if (!is.null(multiple_datasets_object())) {
          print("The integrated object is updated in multiple_datasets_object.")
          shinyjs::enable("add_field")
          shinyjs::enable("add_metadata")
          shinyjs::enable("runScalePCA")
        } else {
          print("Failed to update the integrated object in multiple_datasets_object.")
        }
      }, error = function(e) {
        showNotification(paste0("Error during integration: ", e$message), type = "error")
      })
    })
    
    

    #Add metadata field
    observeEvent(input$add_field, {
      rv_metadata$num_fields <- rv_metadata$num_fields + 1
    })

    # Create a responsive list to store metadata field names
    reactive_metadata_fields <- reactive({
      req(multiple_datasets_object())
      all_metadata_fields <- colnames(multiple_datasets_object()@meta.data)
      exclude_fields <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "integrated_snn_res.0.5", "RNA_snn_res.0.5","RNA_snn_res.0.4","RNA_snn_res.0.3","RNA_snn_res.0.2")
      include_fields <- setdiff(all_metadata_fields, exclude_fields)
      return(include_fields)
    })
    
    
    output$metadata_inputs <- renderUI({
      req(multiple_datasets_object())
      datasets <- unique(multiple_datasets_object()@meta.data$dataset)
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
      tryCatch({
        req(multiple_datasets_object()) # Assurez-vous que l'objet Seurat est chargé
        datasets <- unique(multiple_datasets_object()@meta.data$dataset)
        req(datasets) # Assurez-vous que les datasets sont présents
        
        # Créer une copie temporaire pour la modification
        seurat_temp <- multiple_datasets_object()
        
        for (j in 1:rv_metadata$num_fields) {
          metadata_field_name <- input[[paste0("metadata_name_", j)]]
          req(metadata_field_name) # Assurez-vous que le nom du champ de métadonnées n'est pas nul
          
          if (!metadata_field_name %in% colnames(seurat_temp@meta.data)) {
            seurat_temp@meta.data[[metadata_field_name]] <- NA
          }
          
          for (dataset in datasets) {
            metadata_field_value <- input[[paste0("metadata_value_", dataset, "_", j)]]
            req(metadata_field_value) # Assurez-vous que la valeur du champ de métadonnées n'est pas nulle
            
            rows <- which(seurat_temp@meta.data$dataset == dataset)
            seurat_temp@meta.data[rows, metadata_field_name] <- metadata_field_value
          }
        }
        
        # Mettre à jour l'objet Seurat intégré avec les nouvelles métadonnées
        multiple_datasets_object(seurat_temp)
        
        # Après mise à jour
        print(head(multiple_datasets_object()@meta.data))
      }, error = function(e) {
        showNotification(paste0("Error adding metadata: ", e$message), type = "error")
      })
    })
    
    
    
    # Identify and extract the first loaded Seurat object to obtain the RNA assay gene list
    reactive_gene_list_merge <- reactive({
      req(multiple_datasets_object())
      unique_genes <- rownames(LayerData(multiple_datasets_object(), assay = "RNA", layer = 'counts'))
      c("", unique_genes)
    })


    # Tab 13: Scaling and PCA reduction

    shinyjs::disable("runScalePCA")
    shinyjs::disable("runFindNeighbors")
    shinyjs::disable("runFindClusters")
    
    # Scaling, PCA and Elbowplot
    observeEvent(input$runScalePCA, {
      req(multiple_datasets_object())
      showNotification("Le Scaling et PCA ont commencé...", type = "message")
      tryCatch({
        all_genes <- rownames(multiple_datasets_object())
        seurat_object_temp <- multiple_datasets_object()
        seurat_object_temp <- FindVariableFeatures(seurat_object_temp, selection.method = "vst", nfeatures = 3000)
        seurat_object_temp <- ScaleData(seurat_object_temp, features = all_genes)
        seurat_object_temp <- RunPCA(seurat_object_temp, npcs = 50)
        multiple_datasets_object(seurat_object_temp)
        req(multiple_datasets_object()[["pca"]])
        output$elbow_plot2 <- renderPlot({
          ElbowPlot(multiple_datasets_object())
        })
        shinyjs::enable("runFindNeighbors")
        showNotification("Scaling and PCA are complete.", type = "message")
      }, error = function(e) {
        showNotification(paste0("Scaling and PCA errors:", e$message), type = "error")
      })
    })
    
    
  # Neigbors calculation 
    observeEvent(input$runFindNeighbors, {
      showNotification("Finding neighbors started...", type = "message")
      req(multiple_datasets_object())
      tryCatch({
        seurat_object_temp <- multiple_datasets_object()
        seurat_object_temp <- FindNeighbors(seurat_object_temp, dims = 1:input$dimension_2)
        seurat_object_temp <- RunUMAP(seurat_object_temp, dims = 1:input$dimension_2)
        multiple_datasets_object(seurat_object_temp)
        req(multiple_datasets_object()[["umap"]])
        output$UMAPPlot <- renderPlot({
          UMAPPlot(multiple_datasets_object(), group.by="orig.ident")
        })
        shinyjs::enable("runFindClusters")
        showNotification("Finding neighbors and UMAP completed.", type = "message")
      }, error = function(e) {
        showNotification(paste0("Error during find_neighbors: ", e$message), type = "error")
      })
    })
    
    
    # Observer forfind Clusters and display UMAP
    observeEvent(input$runFindClusters, {
      req(multiple_datasets_object())
      showNotification("Clustering begins...", type = "message")
      tryCatch({
        seurat_integrated_temp <- FindClusters(multiple_datasets_object(), resolution = input$resolution_step2)
        multiple_datasets_object(seurat_integrated_temp)
        req(multiple_datasets_object())
        output$UMAPPlot_cluster <- renderPlot({
          p2 <- DimPlot(multiple_datasets_object(), reduction = "umap", label = TRUE, repel = TRUE)
          p2
        })
        showNotification("Clustering and UMAP rendering completed.", type = "message")
      }, error = function(e) {
        showNotification(paste0("Error during clustering: ", e$message), type = "error")
      })
    })
    
    
    # Tab 14: Visualize gene expression
    
    # Updated selectInput to choose how to group data
    observe({
      updateSelectInput(session, "group_by_select", choices = reactive_metadata_fields())
    })
    
    
    
    # Updated selectInput to choose the genes to visualize  
    observe({
      updateSelectInput(session, "geneInput_merge", choices = reactive_gene_list_merge())
    })
    
    observeEvent(input$geneInput_merge, {
      selected_genes <- input$geneInput_merge
      
      # Update text fields with selected genes
      updateTextInput(session, "gene_list_vln_merge", value = paste(selected_genes, collapse = ", "))
      updateTextInput(session, "gene_list_feature_merge", value = paste(selected_genes, collapse = ", "))
      updateTextInput(session, "gene_list_dot_merge", value = paste(selected_genes, collapse = ", "))
    })
    
    # VlnPlot of the selected gene
    observeEvent(input$runVlnPlot, {
      tryCatch({
        req(input$gene_list_vln_merge)
        gene_list <- unique(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]))
        
        # Assurez-vous que multiple_datasets_object n'est pas NULL
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        
        selected_group_by <- input$group_by_select
        if(all(gene_list %in% rownames(LayerData(seurat_object, assay="RNA", layer='counts')))) {
          output$VlnPlot2 <- renderPlot({
            DefaultAssay(seurat_object) <- "RNA" # S'assurer que l'assay par défaut est "RNA"
            if(input$hide_vln_points_merge) {
              VlnPlot(seurat_object, features = gene_list, group.by = selected_group_by, pt.size = 0)
            } else {
              VlnPlot(seurat_object, features = gene_list, group.by = selected_group_by)
            }
          })
        } else {
          print("Requested genes are not present in the dataset.")
        }
      }, error = function(e) {
        showNotification(paste0("Error during VlnPlot generation: ", e$message), type = "error")
      })
    })
    
    
    
    #FeaturePlot of selected gene
    observeEvent(input$runFeaturePlot, {
      tryCatch({
        req(input$gene_list_feature_merge)
        gene_list <- unique(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]]))
        
        # Assurez-vous que multiple_datasets_object n'est pas NULL
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        
        selected_group_by <- input$group_by_select
        if(all(gene_list %in% rownames(LayerData(seurat_object, assay="RNA", layer='counts')))) {
          DefaultAssay(seurat_object) <- "RNA"
          output$FeaturePlot2 <- renderPlot({
            FeaturePlot(seurat_object, features = gene_list)
          })
        } else {
          print("Requested genes are not present in the dataset.")
        }
      }, error = function(e) {
        showNotification(paste0("Error during FeaturePlot generation: ", e$message), type = "error")
      })
    })
    
    
    # DotPlot
   
    observeEvent(input$runDotPlot, {
      tryCatch({
        req(input$gene_list_dot_merge, multiple_datasets_object())
        gene_list <- unique(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]))
        
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        
        selected_group_by <- input$group_by_select
        if(selected_group_by %in% colnames(seurat_object@meta.data)) {
          if(all(gene_list %in% rownames(LayerData(seurat_object, assay="RNA", layer='counts')))) {
            DefaultAssay(seurat_object) <- "RNA"
            output$DotPlot2 <- renderPlot({
              dot_plot <- DotPlot(seurat_object, features = gene_list, group.by = selected_group_by) + RotatedAxis()
              dot_plot
            })
          } else {
            message("Some of the requested genes are not present in the RNA assay of the dataset.")
          }
        } else {
          message("Selected group by option is not available in Seurat object metadata.")
        }
      }, error = function(e) {
        showNotification(paste0("Error during DotPlot generation: ", e$message), type = "error")
        message(paste0("Error during DotPlot generation: ", e$message))
      })
    })
    
    
    

    # Tab 15: Calculate differentially expressed genes for each cluster in the merged dataset

    # Variables for tab 15
    gene_tables_merge <- reactiveValues()
    all_markers_merge <- reactiveVal()

    shinyjs::disable("download_DE_merged")


    # Function to calculate all differential markers for the merged dataset
    calculate_merged_markers <- function() {
      tryCatch({
        req(multiple_datasets_object())
        markers <- FindAllMarkers(multiple_datasets_object(), min.pct = input$min_pct_all_multiple,
                                  logfc.threshold = input$logfc_threshold_all_multiple,)
        markers <- as.data.frame(markers)
        markers$gene <- paste0('<span class="gene-name">', rownames(markers), '</span>')
        
        markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
        markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
        markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(markers), '">', rownames(markers), '</a>')
        all_markers_merge(markers)
      }, error = function(e) {
        showNotification(paste0("Error calculating markers:", e$message), type = "error")
        NULL
      })
    }
    
    update_gene_tables_display_merge <- function() {
      tryCatch({
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
      }, error = function(e) {
        showNotification(paste0("Error updating gene table display:", e$message), type = "error")
      })
    }
    
    # Observer for run_DE_merged button
    observeEvent(input$run_DE_merged, {
      tryCatch({
        calculate_merged_markers()
        update_gene_tables_display_merge()
        shinyjs::enable("download_DE_merged")
      }, error = function(e) {
        showNotification(paste0("Error running run_DE_merged:", e$message), type = "error")
      })
    })
    
    # Observer for updating the display only when the number of genes changes
    observeEvent(input$number_genes_merge, {
      tryCatch({
        if (!is.null(input$number_genes_merge) && !is.null(all_markers_merge())) {
          update_gene_tables_display_merge()
        }
      }, error = function(e) {
        showNotification(paste0("Erreur lors de la mise à jour des gènes: ", e$message), type = "error")
      })
    })
    
    
    # Observer for gene tables output
    observe({
      tryCatch({
        lapply(names(gene_tables_merge), function(name) {
          output[[name]] <- renderDT({
            datatable(gene_tables_merge[[name]], escape = FALSE)
          })
        })
      }, error = function(e) {
        showNotification(paste0("Error rendering gene tables:", e$message), type = "error")
      })
    })
    
    # Notification for users
    output$previous_tab_notification <- renderUI({
      tryCatch({
        if (!is.null(session$userData$previous_tab_notification_msg)) {
          list(
            tags$hr(),
            tags$strong(session$userData$previous_tab_notification_msg),
            tags$hr()
          )
        }
      }, error = function(e) {
        showNotification(paste0("Erreur lors de la création de la notification: ", e$message), type = "error")
      })
    })
    
    # Download merged seurat object
    output$save_seurat_merge <- downloadHandler(
      filename = function() {
        return("seurat_object.rds")
      },
      content = function(file) {
        # Affichez la boîte modale
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the seurat object for download...",
          easyClose = FALSE,
          footer = NULL
        ))
        tryCatch({
          saveRDS(multiple_datasets_object(), file)
          removeModal()
        }, error = function(e) {
          removeModal()
          showNotification(paste0("Error saving Seurat object: ", e$message), type = "error")
        })
      }
    )
    
    
    # Download markers as a CSV
    output$download_DE_merged <- downloadHandler(
      filename = function() {
        paste("DE_genes_merged_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        tryCatch({
          markers <- all_markers_merge()
          write.csv(markers, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading markers:", e$message), type = "error")
        })
      }
    )
    

    # Tab 16: Assigning cell identity merge
    
    # Reactive variable for that tab
    cluster_colours_merge <- reactiveVal()
  
    observe({
      if (!is.null(multiple_datasets_object())) {
        updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(multiple_datasets_object())))
        update_cluster_colours(multiple_datasets_object())
      }
    })
    
    # Function to update cluster colors
    update_cluster_colours <- function(seurat_object) {
      unique_idents <- unique(Idents(seurat_object))
      current_colours <- scales::hue_pal()(length(unique_idents))
      names(current_colours) <- sort(unique_idents)
      cluster_colours_merge(current_colours)
    }
    
    # Renaming clusters
    observeEvent(input$rename_single_cluster_merge_button, {
      req(input$select_cluster_merge, input$rename_single_cluster_merge, multiple_datasets_object())
      updated_seurat <- multiple_datasets_object()
      if (input$rename_single_cluster_merge %in% unique(Idents(updated_seurat))) {
        cells_to_merge <- which(Idents(updated_seurat) %in% c(input$select_cluster_merge, input$rename_single_cluster_merge))
        Idents(updated_seurat, cells = cells_to_merge) <- input$rename_single_cluster_merge
        showNotification("Clusters merged under the name: ", input$rename_single_cluster_merge, type = "message")
      } else {
        Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster_merge)) <- input$rename_single_cluster_merge
        showNotification("Cluster renamed to: ", input$rename_single_cluster_merge, type = "message")
      }
      
      multiple_datasets_object(updated_seurat)  
      update_cluster_colours(updated_seurat)  
      updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(updated_seurat)))
    })
    
    observeEvent(input$update_colour_merge_button, {
      req(input$select_color_merge, input$select_cluster_merge_color, cluster_colours_merge())
      current_colours <- cluster_colours_merge()
      current_colours[input$select_color_merge] <- input$select_cluster_merge_color
      cluster_colours_merge(current_colours)  
    })
    
    
    # Finale UMAP
    output$umap_finale_merge <- renderPlotly({
      req(multiple_datasets_object(), cluster_colours_merge())
      plot_data <- DimPlot(multiple_datasets_object(), group.by = "ident", pt.size = input$pt_size_merge, label = TRUE, label.size = input$label_font_size_merge) +
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
    diff_genes_compare <- reactiveVal() # Create a new reactive value to store markers_df
    diff_genes_compare_cluster <- reactiveVal()   # Define a new reactive value to store comparison tables
    diff_genes_compare_datasets <- reactiveVal()   # Define a new reactive value to store comparison tables    
    
    shinyjs::disable("download_markers_single_cluster_merge")
    shinyjs::disable("download_markers_multiple_clusters_merge")
    shinyjs::disable("download_diff_dataset_cluster")
    
    # Drop-down menu to filter by Dataset
    output$dataset_filter_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("dataset_filter",
                  label = "Filter by Dataset",
                  choices = c(unique(multiple_datasets_object()@meta.data$dataset)),
                  selected = unique(multiple_datasets_object()@meta.data$dataset),
                  multiple = TRUE)
    })
    
    # Display of filtered UMAP graph
    output$filtered_umap_plot <- renderPlot({
      req(multiple_datasets_object(), input$dataset_filter)
      
      if (!"dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
        stop("The dataset column was not found in the Seurat object's metadata.")
      }
      
      # If all datasets are selected, display them all
      if ("Tous les Datasets" %in% input$dataset_filter || is.null(input$dataset_filter)) {
        return(DimPlot(multiple_datasets_object(), group.by = "ident", label = TRUE, label.size = 5) + NoLegend())
      } else {
        valid_datasets <- input$dataset_filter %in% unique(multiple_datasets_object()@meta.data$dataset)
        subset_seurat <- subset(multiple_datasets_object(), subset = dataset %in% input$dataset_filter[valid_datasets])
        return(DimPlot(subset_seurat, group.by = "ident", label = TRUE, label.size = 5) + NoLegend())
      }
    })
    
    # Reactively update the cluster drop-down menu
    observe({
      req(multiple_datasets_object())
      cluster_choices <- unique(Idents(multiple_datasets_object()))
      updateSelectInput(session, "selected_cluster", choices = cluster_choices, selected = cluster_choices[1])
    })
    
    # Reactive function for markers
    observeEvent(input$calculate_DE, {
      tryCatch({
        showNotification("Calculating differentially expressed genes...", type = "message")
        
        subset_seurat <- multiple_datasets_object()
        
        markers <- FindMarkers(subset_seurat,
                               ident.1 = input$selected_cluster, 
                               min.pct = input$min_pct_merge,
                               logfc.threshold = input$logfc_threshold_merge)
        
        # Vérifiez si markers contient des données
        if (nrow(markers) > 0) {
          markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
          markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
          markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(markers), '">', rownames(markers), '</a>')
          markers_df <- as.data.frame(markers)
          diff_genes_compare(markers_df) 
          shinyjs::enable("download_markers_single_cluster_merge")
        } else {
          showNotification("No differentially expressed genes found for the selected cluster.", type = "message")
        }
      }, error = function(e) {
        showNotification(paste0("Error when calculating differentially expressed genes: ", e$message), type = "error")
      })
    })
    
    output$DE_genes_table <- renderDataTable({
      tryCatch({
        datatable(diff_genes_compare(), escape = FALSE) 
      }, error = function(e) {
        showNotification(paste0("Error when rendering the gene table: ", e$message), type = "error")
      })
    })
    

    # Download gene comparison table
    output$download_markers_single_cluster_merge <- downloadHandler(
      filename = function() {
        paste("genes-differents-", Sys.Date(), ".csv", sep="")      },
      content = function(file) {
        tryCatch({
          req(!is.null(diff_genes_compare))
          diff_genes_compare_DL <- diff_genes_compare()
          diff_genes_compare_DL_sorted <-diff_genes_compare_DL[order(diff_genes_compare_DL$p_val), ]
          write.csv(diff_genes_compare_DL_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )
    
    
    # Drop-down menus to select clusters for comparison
    output$cluster1_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster1_compare",
                  label = "Selec first cluster",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL)
    })

    output$cluster2_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster2_compare",
                  label = "Select second cluster",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL)
    })

    # Cluster comparison
    observeEvent(input$compare_clusters_button, {
      tryCatch({
        showNotification("Finding differentially expressed features...", type = "message")
        req(input$cluster1_compare, input$cluster2_compare, multiple_datasets_object())
        if(input$cluster1_compare == input$cluster2_compare) {
          showNotification("Please select two different clusters for comparison!", type = "error")
          return()
        }
        temp_res <- FindMarkers(multiple_datasets_object(), ident.1 = input$cluster1_compare, ident.2 = input$cluster2_compare, 
                                min.pct = input$min_pct_compare_merge,
                                logfc.threshold = input$logfc_threshold_compare_merge)  
        temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
        temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
        temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')
        shinyjs::enable("download_markers_multiple_clusters_merge")
        diff_genes_compare_cluster(temp_res)
      }, error = function(e) {
        showNotification(paste0("Cluster comparison error:", e$message), type = "error")
      })
    })
    

    output$diff_genes_table_compare <- renderDataTable({
      req(diff_genes_compare_cluster()) 
      tryCatch({
        datatable(diff_genes_compare_cluster(), escape = FALSE) 
      }, error = function(e) {
        showNotification(paste0("Error when rendering the gene table: ", e$message), type = "error")
      })
    })
    
    # Download gene comparison table
    output$download_markers_multiple_clusters_merge <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster1_compare, "-VS-", input$cluster2_compare, "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        tryCatch({
          diff_genes_compare_cluster_DL <- diff_genes_compare_cluster()
          req(!is.null(diff_genes_compare_cluster_DL))
          diff_genes_compare_cluster_DL_sorted <- diff_genes_compare_cluster_DL[order(diff_genes_compare_cluster_DL$p_val), ]
          write.csv(diff_genes_compare_cluster_DL_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )
    

    output$dataset1_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("dataset1_compare",
                  label = "Select the first dataset",
                  choices = unique(multiple_datasets_object()@meta.data$dataset),
                  selected = NULL)
    })
    
    output$dataset2_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("dataset2_compare",
                  label = "Select the second dataset",
                  choices = unique(multiple_datasets_object()@meta.data$dataset),
                  selected = NULL)
    })
    
    output$cluster_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster_compare",
                  label = "Select a cluster",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL)
    })
    
    # Comparez les gènes différentiellement exprimés entre datasets pour un cluster spécifique
    observeEvent(input$compare_datasets_button, {
      tryCatch({
        req(input$dataset1_compare, input$dataset2_compare, multiple_datasets_object())
        if(input$dataset1_compare == input$dataset2_compare) {
          showNotification("Please select two different datasets for comparison!", type = "error")
          return()
        }
        if(input$all_clusters) {
          cluster_data <- multiple_datasets_object()
        } else {
          req(input$cluster_compare)
          cluster_data <- subset(multiple_datasets_object(), seurat_clusters == input$cluster_compare)
        }
        diff_genes <- FindMarkers(cluster_data, 
                                  group.by = "dataset", 
                                  ident.1 = input$dataset1_compare, 
                                  ident.2 = input$dataset2_compare, 
                                  min.pct = input$min_pct_compare_dataset_merge,
                                  logfc.threshold = input$logfc_threshold_datasets,
                                  )
        
        if (nrow(diff_genes) == 0) {
          showNotification("No differentially expressed genes found!", type = "error")
          return()
        }
        diff_genes$p_val <- formatC(diff_genes$p_val, format = "e", digits = 3)
        diff_genes$p_val_adj <- formatC(diff_genes$p_val_adj, format = "e", digits = 3)
        
        diff_genes$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(diff_genes), '">', rownames(diff_genes), '</a>')
        
        shinyjs::enable("download_diff_dataset_cluster")
        diff_genes_compare_datasets(diff_genes)
        
        output$diff_dataset_cluster <- renderDataTable({
          datatable(diff_genes, escape = FALSE)
        })
      }, error = function(e) {
        showNotification(paste0("Error when comparing datasets for a specific cluster: ", e$message), type = "error")
      })
    })
    
    # To download the table of diff genes
    output$download_diff_dataset_cluster <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster_compare, "-DS1-", input$dataset1_compare, "-DS2-", input$dataset2_compare, "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        tryCatch({
          diff_genes_data <- diff_genes_compare_datasets()
          req(!is.null(diff_genes_data))
          
          diff_genes_data$p_val <- formatC(diff_genes_data$p_val, format = "e", digits = 2)
          diff_genes_data$p_val_adj <- formatC(diff_genes_data$p_val_adj, format = "e", digits = 2)
          
          diff_genes_data_sorted <- diff_genes_data[order(diff_genes_data$p_val), ]
          write.csv(diff_genes_data_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table between datasets: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )
    

    # Tab 18: Subseting seurat object
    
    # Variable réactive pour cet onglet
    subset_seurat_merge <- reactiveVal(NULL)
    shinyjs::disable("download_subset_merge")
    
    # Mise à jour des choix d'identité de cellules pour le sous-ensemble
    observe({
      if (!is.null(multiple_datasets_object())) {
        updateSelectInput(session, "select_ident_subset_merge", choices = unique(Idents(multiple_datasets_object())))
      }
    })
    
    # Réinitialisation de l'objet subset_seurat_merge lors du changement de multiple_datasets_object
    observe({
      if (!is.null(multiple_datasets_object())) {
        subset_seurat_merge(multiple_datasets_object())
      }
    })
    
    # Affichage de l'UMAP global
    output$global_umap_merge <- renderPlot({
      req(multiple_datasets_object())
      plot_data <- DimPlot(multiple_datasets_object(), group.by = "ident", label = TRUE) +
        theme(axis.line = element_line(size = 0.5)) +
        NoLegend()
      print(plot_data)
    })
    
    # Affichage de l'UMAP du sous-ensemble
    output$subset_umap_merge <- renderPlot({
      req(subset_seurat_merge())
      plot_data <- DimPlot(subset_seurat_merge(), group.by = "ident", label = TRUE) +
        theme(axis.line = element_line(size = 0.5)) +
        NoLegend()
      print(plot_data)
    })
    
    # Téléchargement de l'objet Seurat sous-ensemble
    output$download_subset_merge <- downloadHandler(
      filename = function() {
        paste("seurat_subset-", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the seurat object for download...",
          easyClose = FALSE,
          footer = NULL
        ))
        
        tryCatch({
          req(subset_seurat_merge())
          saveRDS(subset_seurat_merge(), file)
          removeModal()
        }, error = function(e) {
          removeModal()
          showNotification(paste0("Error while downloading Seurat subset: ", e$message), type = "error")
        })
      }
    )
    
    # Mise à jour de la sélection de sous-ensemble
    observeEvent(input$apply_subset_merge, {
      tryCatch({
        req(multiple_datasets_object())
        
        # Si aucun identifiant n'est sélectionné, retourne l'objet original
        if (length(input$select_ident_subset_merge) == 0) {
          subset_seurat_merge(multiple_datasets_object())
          shinyjs::enable("download_subset_merge")
          return()
        }
        
        subsetted_seurat <- subset(multiple_datasets_object(), idents = input$select_ident_subset_merge)
        
        # Vérifiez s'il y a des cellules dans le sous-ensemble
        if (nrow(subsetted_seurat@meta.data) == 0) {
          showNotification("No cells found with the selected identities.", type = "error")
          return()
        }
        
        subset_seurat_merge(subsetted_seurat)
        shinyjs::enable("download_subset_merge")
      }, error = function(e) {
        showNotification(paste0("Error updating subset selection: ", e$message), type = "error")
      })
    })
    # Sous-ensemble basé sur l'expression génique (au moins un gène exprimé)
    observeEvent(input$apply_gene_subset_merge, {
      tryCatch({
        req(multiple_datasets_object())
        gene_list <- unlist(strsplit(input$gene_list_merge, ","))
        gene_list <- trimws(gene_list)
        
        if (length(gene_list) == 0) {
          showNotification("Please enter valid gene names.", type = "error")
          return()
        }
        
        # Créer une matrice logique pour l'expression génique
        expression_matrix <- sapply(gene_list, function(gene) {
          FetchData(multiple_datasets_object(), vars = gene) >= input$expression_threshold_merge
        })
        
        # Sélection des cellules exprimant au moins un des gènes spécifiés
        cells_to_keep <- which(rowSums(expression_matrix) > 0)
        
        # Vérifier si des cellules sont sélectionnées
        if (length(cells_to_keep) == 0) {
          showNotification("No cells found with the specified gene expression criteria.", type = "error")
          return()
        }
        
        subsetted_seurat <- subset(multiple_datasets_object(), cells = cells_to_keep)
        subset_seurat_merge(subsetted_seurat)
        shinyjs::enable("download_subset_merge")
      }, error = function(e) {
        showNotification(paste0("Error during gene-based subset: ", e$message), type = "error")
      })
    })
    
}
