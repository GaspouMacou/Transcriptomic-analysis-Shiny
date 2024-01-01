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

atac_server <- function(input, output, session) {
  
  shinyjs::useShinyjs()
  seurat_object_atac <- reactiveVal(NULL)
  
  
  
  
# Tab 1 loading ATAC object:
  
  output$fileInputs_atac <- renderUI({
    lapply(1:input$num_datasets_atac, function(i) {
      fluidRow(
        column(4, textInput(paste0('dataset_name_atac', i), paste0('ATAC-seq Dataset Name ', i))),
        column(4, selectInput(paste0('data_type_atac', i), 'Choose Data Type:',
                              choices = c('ATAC-seq h5' = 'atac_seq_h5', 'Multiome zip' = 'multiome_zip'))),
        column(4, fileInput(paste0('file_atac', i), paste0('Choose File ', i), accept = c(".h5", ".zip")))
      )
    })
  })
  
  
  # Function to process Multiome zip data
  process_multiome_zip_data <- function(zip_path, dataset_name, species_choice_atac) {
    tryCatch({
      showNotification("Uploading and processing data...", type = "message")
      mt_pattern <- ifelse(species_choice_atac == "mouse", "^mt-", "^MT-")
      
      if (!file.exists(zip_path)) {
        stop("Specified zip file does not exist.")
      }
      if (dir.exists("unzipped")) {
        unlink("unzipped", recursive = TRUE)
      }
      
      unzip(zip_path, exdir = "unzipped")
      single_dataset_data <- Read10X("unzipped")
      rna_data <- single_dataset_data$`Gene Expression`
      atac_data <- single_dataset_data$`Peaks`
      
      rna_seurat <- CreateSeuratObject(counts = rna_data, project = paste0(dataset_name, "_RNA"), min.cells = 3, min.features = 200)
      atac_seurat <- CreateSeuratObject(counts = atac_data, project = paste0(dataset_name, "_ATAC"), min.cells = 3, min.features = 200)
      
      rna_seurat[["percent.mt"]] <- PercentageFeatureSet(rna_seurat, pattern = mt_pattern)
      
      showNotification("Multiome data processed successfully for dataset: ", dataset_name, type = "message")
      return(list(rna = rna_seurat, atac = atac_seurat))
    }, error = function(e) {
      showNotification("Error processing Multiome data: ", e$message, type = "error")
      return(NULL)
    })
  }
  
  
  # Intégration des datasets ATAC-seq si nécessaire
  integrate_atac_datasets <- function(atac_datasets) {
    # Votre code d'intégration ici
    # Retourne un objet Seurat intégré
  }
  
  # Observe file inputs and update reactive Seurat object
  observe({
    lapply(1:input$num_datasets_atac, function(i) {
      observeEvent(input[[paste0('file_atac', i)]], {
        if (!is.null(input[[paste0('file_atac', i)]])) {
          zip_path <- input[[paste0('file_atac', i)]]$datapath
          dataset_name <- input[[paste0('dataset_name_atac', i)]]
          
          seurat_object <- process_atac_zip_data(zip_path, dataset_name)
          
          if (!is.null(seurat_object)) {
            current_data <- seurat_object_atac()
            if (is.null(current_data)) {
              current_data <- list()
            }
            current_data[[dataset_name]] <- seurat_object
            # Vérifier si l'intégration est nécessaire
            if (length(current_data) > 1) {
              integrated_seurat <- integrate_atac_datasets(current_data)
              seurat_object_atac(integrated_seurat)
            } else {
              seurat_object_atac(seurat_object) # Stocker un seul objet Seurat si un seul dataset est chargé
            }
          }
        }
      })
    })
  })
  
  # Loading seurat object
  observeEvent(input$load_seurat_file_atac, {
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
  
  # Fonction pour traiter les données ATAC-seq à partir d'un fichier zip
  process_atac_zip_data <- function(zip_path, dataset_name) {
    tryCatch({
      showNotification("Uploading and processing ATAC-seq data...", type = "message")
      print(paste("Processing ATAC-seq data for dataset:", dataset_name))
      
      # Décompression du fichier zip
      print("Checking if zip file exists...")
      if (!file.exists(zip_path)) {
        stop("Specified zip file does not exist.")
      }
      print("Unzipping the file...")
      if (dir.exists("unzipped_atac")) {
        unlink("unzipped_atac", recursive = TRUE)
      }
      unzip(zip_path, exdir = "unzipped_atac")
      
      # Lecture des fichiers décompressés
      print("Reading unzipped files...")
      h5_file <- list.files("unzipped_atac", pattern = "\\.h5$", full.names = TRUE)
      csv_file <- list.files("unzipped_atac", pattern = "\\.csv$", full.names = TRUE)
      tsv_file <- list.files("unzipped_atac", pattern = "\\.tsv\\.gz$", full.names = TRUE)
      
      print("Checking the existence of required files...")
      if (length(h5_file) != 1 || length(tsv_file) != 1) {
        stop("Required files (.h5, .csv, .tsv.gz)not found in the zip archive.")
      }
      
      print("Reading counts from h5 file...")
      counts <- Read10X_h5(filename = h5_file)
      metadata <- read.csv(file = csv_file, header = TRUE, row.names = 1)
      print("Creating ChromatinAssay...")
      chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        fragments = tsv_file,
        min.cells = 10,
        min.features = 200
      )
      
      print("Creating Seurat object...")
      seurat_object <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = metadata
      )
      
      # Extraction et ajout des annotations génomiques
      print("Extracting and adding genomic annotations...")
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
      genome(annotations) <- "hg19"
      Annotation(seurat_object) <- annotations
      
      print("ATAC-seq data processed successfully.")
      print(granges(seurat_object)) # Affichage des informations GRanges
      
      showNotification("ATAC-seq data processed successfully for dataset: ", dataset_name, type = "message")
  
      return(seurat_object)
    }, error = function(e) {
      print(paste("Error processing ATAC-seq data:", e$message))
      showNotification("Error processing ATAC-seq data: ", e$message, type = "error")
      return(NULL)
    })
  }
  
 
  # Tab 2 QC for ATAC-seq
  
  # Fonction pour calculer les métriques de QC avec gestion d'erreurs
  calculate_qc_metrics <- function(seurat_object) {
    tryCatch({
      message("Calculating QC metrics for ATAC-seq data...")
      
      # Ajout de l'argument 'slot'
      # Remplacer 'counts' par le slot approprié si nécessaire
     
      seurat_object <- NucleosomeSignal(object = seurat_object, assay="peaks")
      seurat_object <- TSSEnrichment(object = seurat_object, fast = FALSE, assay= "peaks")
      
      # Calcul et ajout des métriques QC
      seurat_object$pct_reads_in_peaks <- seurat_object$peak_region_fragments / seurat_object$passed_filters * 100
      seurat_object$blacklist_ratio <- seurat_object$blacklist_region_fragments / seurat_object$peak_region_fragments
      
      print("QC metrics calculated successfully.")
      print(head(seurat_object@meta.data)) # Affichage des premières lignes des métadonnées
      
      return(seurat_object)
    }, error = function(e) {
      print(paste("Error during QC metric calculation:", e$message))
      stop("Error calculating QC metrics: ", e$message)
    })
  }
  
  # Observateur pour le bouton de calcul des métriques de QC
  observeEvent(input$calculate_qc_btn, {
    if (!is.null(seurat_object_atac())) {
      message("Starting QC metric calculation...")
      updated_object <- calculate_qc_metrics(seurat_object_atac())
      
      # Mettre à jour l'objet Seurat réactif avec les nouvelles métriques de QC
      seurat_object_atac(updated_object)
      
      message("QC metrics calculated successfully.")
      print("Updated Seurat object after QC metrics calculation:")
      print(head(updated_object@meta.data)) # Affichage des premières lignes des métadonnées mises à jour
      showNotification("QC metrics calculated successfully.", type = "message")
    } else {
      message("No ATAC-seq data available for QC calculation.")
      showNotification("No ATAC-seq data available for QC calculation.", type = "error")
    }
  })
  
  
  
  # Rendu pour DensityScatter
  output$density_scatter_plot <- renderPlot({
    req(seurat_object_atac())
    if ("TSS.enrichment" %in% names(seurat_object_atac()@meta.data)) {
      DensityScatter(seurat_object_atac(), x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
    }
  })
  
  # Rendu pour TSSPlot
  output$tss_plot <- renderPlot({
    req(seurat_object_atac())
    if ("TSS.enrichment" %in% names(seurat_object_atac()@meta.data)) {
      seurat_object <- seurat_object_atac()
      seurat_object$high.tss <- ifelse(seurat_object$TSS.enrichment > 3, 'High', 'Low')
      TSSPlot(seurat_object, group.by = 'high.tss') + NoLegend()
    }
  })
  
  # Rendu pour l'histogramme des fragments
  output$fragment_histogram_pYeslot <- renderPlot({
    req(seurat_object_atac())
    if ("nucleosome_signal" %in% names(seurat_object_atac()@meta.data)) {
      seurat_object <- seurat_object_atac()
      seurat_object$nucleosome_group <- ifelse(seurat_object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
      FragmentHistogram(object = seurat_object, group.by = 'nucleosome_group')
    }
  })
  

}
