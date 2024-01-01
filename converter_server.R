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

  
converter_server <- function(input, output, session) {
  
  shinyjs::useShinyjs()
  
  seurat_object_loaded <- reactiveVal(NULL) # Seurat object loaded by the user 
  path_to_converted_file <- reactiveVal(NULL)
  path_to_monocle_file <- reactiveVal(NULL)  # Pour Monocle
  
  # Chargement de l'objet Seurat
  observeEvent(input$file_seurat_conversion, {
    message("Attempt to read the file: ", input$file_seurat_conversion$datapath)
    tryCatch({
      loaded_seurat <- readRDS(input$file_seurat_conversion$datapath)
      message("File read successfully.")
      seurat_object_loaded(loaded_seurat)
      showNotification("The Seurat object has been successfully loaded!")
    }, error = function(e) {
      showNotification(paste("An error has occurred: ", e), type = "error")
    })
  })
  
  observeEvent(input$convert_button_anndata, {
    showModal(modalDialog(
      title = "Please wait",
      "Conversion en cours...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    message("Start of conversion...")
    
    tryCatch({
      seurat_for_conversion <- seurat_object_loaded()
      temp_seurat_file <- tempfile(fileext = ".h5Seurat")
      message("Save Seurat object in HDF5 format...")
      
      SaveH5Seurat(seurat_for_conversion, filename = temp_seurat_file, overwrite = TRUE)
      
      converted_file_path <- tempfile(fileext = ".h5ad")
      message("Convert to anndata format...")
      Convert(source = temp_seurat_file, dest = converted_file_path, assay = "RNA")
      
      file.remove(temp_seurat_file)
      path_to_converted_file(converted_file_path)
      
      removeModal()
      
      showNotification("Successful conversion!")
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Conversion error : ", e$message), type = "error")
      print(e)
      traceback()
    })
  })
  
  # Téléchargement du fichier converti (AnnData)
  output$download_anndata_converted_file <- downloadHandler(
    filename = function() {
      paste("anndata_object_", Sys.Date(), ".h5ad", sep="")
    },
    content = function(file) {
      showModal(modalDialog(
        title = "Please wait",
        "Download in progress...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      file_to_download <- path_to_converted_file()
      if(!is.null(file_to_download) && file.exists(file_to_download)) {
        file.copy(file_to_download, file)
      }
      
      removeModal()
    }
  )
  
  # Conversion à Monocle
  observeEvent(input$convert_button_monocle, {
    showModal(modalDialog(
      title = "Please wait",
      "Conversion to Monocle in progress...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      seurat_for_monocle <- seurat_object_loaded()
      sce <- as.cell_data_set(seurat_for_monocle)
      monocle_temp <- as(sce, "cell_data_set")
      
      monocle_temp <- estimate_size_factors(monocle_temp)
      
      rowData(monocle_temp)$gene_short_name <- rownames(seurat_object_loaded()[["RNA"]])
      
      monocle_file_path <- tempfile(fileext = ".RDS")
      saveRDS(monocle_temp, file=monocle_file_path)
      path_to_monocle_file(monocle_file_path)
      
      removeModal()
      
      showNotification("Successful conversion to Monocle!")
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error converting to Monocle : ", e$message), type = "error")
      print(e)
      traceback()
    })
  })
  
  # Téléchargement de l'objet Monocle
  output$download_monocle_converted_file <- downloadHandler(
    filename = function() {
      paste("monocle_object_", Sys.Date(), ".RDS", sep="")
    },
    content = function(file) {
      showModal(modalDialog(
        title = "Please wait",
        "Download Monocle object in progress...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      monocle_to_download <- path_to_monocle_file()
      if(!is.null(monocle_to_download) && file.exists(monocle_to_download)) {
        file.copy(monocle_to_download, file)
      }
      
      removeModal()
    }
  )
}
