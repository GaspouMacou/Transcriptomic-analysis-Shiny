

# Créer l'interface utilisateur

ui <- dashboardPage(
             # Ajouter un en-tête de page
                  dashboardHeader(title = tags$img(src = "header.png",   height = "53px", width = "218px")),
                  # Ajouter une barre latérale
                  dashboardSidebar(
                    tags$head(
                      tags$style(HTML("
                        /* Change header color */
                        .main-header .navbar {
                          background-color: rgb(122, 207, 176) !important;
                        }
                        /* Aligner le contenu de la sidebar avec le contenu principal */
                        .sidebar .sidebar-menu {
                          margin-top: 5px;
                        }
                  
                        /* Augmenter la largeur de la bulle d'information */
                        .popover {
                          width: 300px;
                        }
                  
                        /* Style pour ajouter l'image au header */
                        .main-header {
                          background-image: url('www/header.png');
                          background-repeat: no-repeat;
                          background-position: center center;
                        }
                        /* Style pour info-box */
                        .info-box { 
                          width: 100%; 
                        }
                      ")),
                      tags$script(src = "my_script.js")
                    ),
                    useShinyjs(),
                    # Ajouter des onglets à la barre latérale
                    sidebarMenu(
                      menuItem("Introduction", tabName = "introduction", icon = icon("mug-hot")),
                      menuItem(
                        "Single Dataset Analysis", icon = icon("virus-covid"),
                        menuItem("Load dataset", tabName = "intro", icon = icon("file-arrow-down")),
                        menuItem("Data cleanup & Variable features", tabName = "qc", icon = icon("bug")),
                        menuItem("Dimensional reduction", tabName = "perform_linear_dimensional_reduction", icon = icon("home")),
                        menuItem("Heat maps", tabName = "heat_maps", icon = icon("thermometer-half")),
                        menuItem("Dimensionality of the dataset", tabName = "dimensionality_of_the_dataset", icon = icon("th-large")),
                        menuItem("Clustering", tabName = "cluster_the_cells", icon = icon("brain")),
                        menuItem("Differentially expressed genes", tabName = "gene_cluster", icon = icon("react")),
                        menuItem("Genes expressions", tabName = "visualizing_marker_expression", icon = icon("flask")),
                        menuItem("Assigning cell type identity", tabName = "assigning_cell_type_identity", icon = icon("id-card")),
                        menuItem("Cluster biomarkers", tabName = "cluster_biomarkers", icon = icon("earth-americas")),
                        menuItem("Subset", tabName = "subset", icon = icon("project-diagram"))
                      ),
                      menuItem(
                        "Multiple Datasets Analysis", icon = icon("viruses"),
                        menuItem("Load datasets", tabName = "merge_dataset", icon = icon("file-arrow-down")),
                        menuItem("Clustering", tabName = "plot_merge", icon = icon("brain")),
                        menuItem("Differentially expressed genes", tabName = "DE_merged_dataset", icon = icon("react")),
                        menuItem("Genes expressions", tabName = "visualization_merge", icon = icon("flask")),
                        menuItem("Assigning cell type identity", tabName = "assigning_cell_type_identity_merge", icon = icon("id-card")),
                        menuItem("Cluster biomarkers", tabName = "combined_analysis", icon = icon("earth-americas")),
                        menuItem("Subset", tabName = "subset_merge", icon = icon("project-diagram"))
                      ),
                      menuItem("ATAC seq",  icon = icon("viruses"),
                               menuItem("Load ATAC-seq data", tabName = "atac_load", icon = icon("file-arrow-down")),
                               menuItem("ATAC prepocessing", tabName = "atac_preprocessing", icon = icon("file-arrow-down")),
                               menuItem("Load ATAC-seq data", tabName = "atac_load", icon = icon("file-arrow-down")),
                               menuItem("Load ATAC-seq data", tabName = "atac_load", icon = icon("file-arrow-down"))

                     ),
                      menuItem(
                        "Trajectory Analysis", icon = icon("project-diagram"),
                        menuItem("Trajectory Analysis", tabName = "trajectory", icon = icon("clock")),
                        menuItem("Genes pseudotime", tabName = "genes_pseudotime", icon = icon("fa-solid fa-hourglass"))
                      ),
                      menuItem("The converter", tabName = "converter",icon = icon("file")),
                      menuItem("Acknowledgement & Licence", tabName = "acknowledgement", icon = icon("heart"))
                    )
                  ),

                    # Ajouter un corps principal avec des onglets
                    dashboardBody(
                      tabItems(
                        tabItem(
                          tabName = "introduction",
                          h1("Introduction"),
                          p("Welcome to our Single-Cell and Single-Nucleus RNA sequencing (scRNA-seq and snRNA-seq) data analysis application.
                  The data you will be exploring are produced using 10x Genomics technology, a cutting-edge method for genome-wide profiling of single cells and nuclei. This allows us to analyze gene expression at the level of a single cell or nucleus, offering unprecedented resolution for understanding biological mechanisms.
                  This application utilizes the Seurat library, a popular R platform for scRNA-seq and snRNA-seq analysis. Seurat offers a suite of tools for quality, analysis, exploration, and visualization of such data.
                  A key principle of this application is the clustering of cells based on their gene expression profiles. This means we group cells into subsets (clusters) based on the similarity of their gene expression profiles. These clusters can often correspond to different cell types or distinct cellular states.
                  Additionally, we focus on differentially expressed genes, that is, genes that are significantly more or less expressed in one cluster compared to others. These differentially expressed genes provide us with valuable clues about the unique characteristics of cells in each cluster.
                  We hope you will find this application useful and informative. Happy exploring!"),
                          div(style = "text-align: center; display: block;",
                              tags$img(src = 'pipeline.png', style="max-height: 100vh; max-width: 140vw; height: auto; width: auto;")
                          )
                        )
                      ,

                        # Ajouter le contenu de l'onglet Introduction
                        tabItem(tabName = "intro",
                                fluidPage(
                                  # Ajouter un titre
                                  h2("Single Dataset Analysis"),
                                  p("To load your transcriptomics data, please compress the three 10X files (barcodes.tsv.gz, matrix.mtx.gz and features.tsv.gz) into a single .zip file."),

                                  # Adding a row to place the radio buttons side by side
                                  fluidRow(
                                    column(6,
                                           radioButtons("species_choice", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "mouse")
                                    ),
                                    column(6,
                                           radioButtons("dataset_type", "Choose Dataset Type", choices = list("snRNA-seq" = "snRNA", "Multiome" = "multiome"))
                                    )
                                  ),

                                  fileInput('file', 'Choose a ZIP file for the experiment', accept = c('.zip')),
                                  fileInput("load_seurat_file", "Choose Seurat Object", accept = ".rds"),

                                  tableOutput('contents'),

                                  p("This section provides a comprehensive guide for uploading and analyzing single-cell transcriptomic data using the advanced 10X Genomics technology. Follow the steps below to navigate through the data processing pipeline."),

                                  h3(style = "font-weight: bold;", "Uploading Data"),
                                  p("Begin by compressing and uploading your 10X Genomics files (barcodes.tsv.gz, matrix.mtx.gz, features.tsv.gz) into a single .zip file. For pre-processed data, you can directly load a Seurat object (.rds file) for analysis."),

                                  h3(style = "font-weight: bold;", "Quality Control Metrics"),
                                  p("Quality control is crucial in single-cell data analysis. Set thresholds for gene detection and mitochondrial DNA contamination. Apply quality control filters to remove low-quality cells and normalize the data to prepare for further analysis."),

                                  h3(style = "font-weight: bold;", "Dimensionality Reduction"),
                                  p("Reduce the complexity of your dataset with Principal Component Analysis (PCA). Use Elbow plots and JackStraw analysis to determine the optimal number of principal components for capturing the most significant data variations."),

                                  h3(style = "font-weight: bold;", "Clustering Cells"),
                                  p("Cluster your cells based on gene expression profiles to identify distinct cellular populations or states. Adjust the resolution parameter to fine-tune the granularity of these clusters, revealing subtle cellular subtypes or broader cell populations."),

                                  h3(style = "font-weight: bold;", "Differential Expression Analysis"),
                                  p("Explore the differentially expressed genes across various clusters. This step is pivotal for identifying potential biomarkers and understanding the biological functions and pathways involved in each cellular subset."),

                                  h3(style = "font-weight: bold;", "Visualization"),
                                  p("Visualize your analysis results with various graphical representations. Utilize UMAP plots for cluster visualization, gene expression heatmaps for in-depth gene analysis, and FeaturePlots for specific marker exploration."),

                                  h3(style = "font-weight: bold;", "Assigning Cell Type Identity"),
                                  p("Annotate each cluster with cell type identities using known genetic markers. This step helps in interpreting the biological context of your data, associating clusters with specific cell types or states."),

                                  h3(style = "font-weight: bold;", "Subset Analysis"),
                                  p("Perform targeted analysis by focusing on specific clusters or genes. This allows for a more detailed exploration of particular aspects of your data, such as examining the behavior of certain cell types or the expression of key genes.")
                                )
                        ),

                        # Ajouter le contenu de l'onglet QC Metrics
                        tabItem(
                          tabName = "qc",

                          sidebarLayout(

                            # Sidebar panel for inputs ----
                            sidebarPanel(

                              div(style = "display: inline-block; width: 80%;",
                                  sliderInput("nFeature_range", "Unique genes detected in each cell", min = 0, max = 5000, value = c(200, 2500))
                              ),
                              div(style = "display: inline-block; width: 18%;",actionButton("exclamation1",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "The number of unique genes detected in each cell. <ul><li>Low-quality cells or empty droplets will often have very few genes</li><li>Cell doublets or multiplets may exhibit an aberrantly high gene count.</li><li>Be careful, in most cases it's better to stay between 200 and 2500</li></ul>")
                              )
                              ,
                              div(style = "display: inline-block; width: 80%;",sliderInput("percent.mt_max", "Maximum value for  mitochondrial genome:", value = 5, min = 0, max = 100)
                              ),
                              div(style = "display: inline-block; width: 18%;",actionButton("exclamation3", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-content` = "It is normal to have mitochondrial DNA contamination for single cell, but be careful with this value for single nuclei")
                              ),
                              div(style = "display: inline-block; width: 80%;",sliderInput(inputId = "scale_factor", label = "Scale factor", value = 10000, min = 500, max = 50000, step = 500)
                              ),
                              div(style = "display: inline-block; width: 18%;",actionButton("exclamation3", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-content` = "Feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result")
                              ),
                              actionButton("QCmetrics", "QC metrics plot"),tags$br(),
                              actionButton("show_plots", "Feature Scatter Plots"),tags$br(),
                              actionButton("apply_qc", "Apply QC Filters"),
                              actionButton("normalize_data", "Normalize Data")
                            ),
                            mainPanel(
                              infoBoxOutput("nuclei_count"),

                              plotOutput("vlnplot"),
                              fluidRow(
                                column(6, plotOutput("scatter_plot1")),
                                column(6, plotOutput("scatter_plot2"))
                              ),
                              fluidRow(
                                column(8, plotOutput("variable_feature_plot"))
                              )))),


                        tabItem(
                          tabName = "perform_linear_dimensional_reduction",
                          div(
                            class = "flex-container",
                            div(
                              class = "button-group",
                              actionButton("scale_button", "Scale Data"),
                              actionButton("exclamation4",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",
                                           `data-content` = "Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:<ul><li>Shifts the expression of each gene, so that the mean expression across cells is 0</li><li>Scales the expression of each gene, so that the variance across cells is 1</li><li>This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate") ),
                            div(
                              class = "button-group",
                              actionButton("pca_button", "Run PCA", disabled = TRUE), actionButton("exclamation5",label = icon("exclamation-triangle"),
                                                                                                   `data-toggle` = "popover", `data-html` = "true",
                                                                                                   `data-content` = "PCA is a very powerful technique and can have much utility if you can get a grasp of it and what it means. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed (for example, a data-set consisting of a large number of samples, each with their own data points/variables). It extracts the fundamental structure of the data without the need to build any model to represent it. This ‘summary’ of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the ‘principal components'), whilst at the same time being capable of easy interpretation on the original data.")) ),
                          verbatimTextOutput("pca_results"),
                          plotOutput("loading_plot"),
                          plotOutput("dim_plot")
                        ),


                        # Ajouter le contenu de l'onglet Heat
                        tabItem(
                          tabName = "heat_maps",
                          h2("Heatmaps"),
                          sliderInput("dims", "Number of dimensions:", min = 1, max = 15, value = 1),
                          sliderInput("cells", "Number of cells :", min = 100, max = 1000, value = 500),
                          selectInput("gene_select_heatmap",
                                      "Select genes for Heatmap:",
                                      choices = NULL, # Les choix seront mis à jour depuis le serveur
                                      selected = "",
                                      multiple = TRUE) ,

                          textInput("selected_genes_heatmap", "Selected genes for Heatmap:", ""),

                          actionButton("runHeatmap", "Generate Heatmap"),
                          plotOutput("heatmap")
                        )
                        ,

                        # Ajouter le contenu de l'onglet dimensionality_of_the_dataset
                        tabItem(
                          tabName="dimensionality_of_the_dataset",
                          div(
                            class = "flex-container",
                            div(
                              class = "button-group",
                              actionButton("run_elbow", "Show ElbowPlot"),
                              actionButton("exclamation6",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one") ),
                            plotOutput("elbow_plot"),
                            div(
                              class = "button-group",
                              actionButton("run_jackstraw", "Launch JackStraw analysis"), actionButton("exclamation7",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "The JackStrawPlot function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).")) ),
                          numericInput("num_dims", "Number of dimensions for JackStraw", value = 20, min = 1),
                          plotOutput("jackstraw_plot"),
                        ),

                        # Ajouter le contenu de l'onglet cluster the cells
                        tabItem(
                          tabName = "cluster_the_cells",
                          div(
                            class = "flex-container",
                            div(
                              class = "button-group",
                              numericInput("dimension_1", "Number of dimensions :", value = 10, min = 1),
                              actionButton("run_neighbors", "Find neighbors"),
                              actionButton("exclamation8",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "In single-cell analyses, calculating nearest neighbors involves determining which cells are most similar to each other based on their gene expression profiles. This is vital for several downstream analyses, including clustering and defining cellular trajectories. The choice of dimensions influences how similarities are computed, directly impacting the precision of identified neighbors. Selecting an optimal number of dimensions—often via methods like PCA—ensures that the variability explained is maximized while avoiding overfitting, thus providing a robust basis for identifying true biological relationships among cells.") ),
                            plotOutput("neighbors_plot"),
                            div(
                              class = "button-group",
                              sliderInput("resolution_step1", "Resolution :", min = 0.1, max = 2, step = 0.1, value = 0.5),
                              actionButton("run_clustering", "Find clusters"), actionButton("exclamation9",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "Clustering in single-cell analysis groups cells together based on the similarity of their expression profiles, aiming to identify distinct cellular populations or states. The resolution parameter in clustering algorithms, like those used in Seurat, controls the granularity of these identified clusters: a higher resolution often yields more, smaller clusters, whereas a lower resolution generates fewer, larger clusters. Selecting an appropriate resolution is pivotal, as it influences the biological interpretations by determining the discernibility of subtle cellular subtypes and the comprehensiveness of identified cell populations.")) ),
                          plotOutput("clustering_plot"),
                        )
                        ,

                        # Ajouter le contenu de l'onglet umap
                        tabItem(
                          tabName = "gene_cluster",
                          column(width = 4,
                            actionButton("run_DE", "Finding differentially expressed genes"),
                            numericInput("logfc_threshold_all_single", label = "Log2 Fold Change threshold", value = 0.25),
                            ),
                        column(width = 4,
                             downloadButton("save_seurat", "Save Seurat Object"),
                             numericInput("min_pct_all_single", "Percentage threshold ", value = 0.01, min = 0, max = 1, step = 0.01)
                            ),
                        column(width = 4,
                               downloadButton('download_DE', 'Download differentially expressed genes'),
                               sliderInput("number_genes","Number of genes:",min = 10,max = 1000, value = 10,step = 50)
                          ),
                          uiOutput("diff_genes_tables")
                        )
                        ,

                        # visualizing_marker_expression
                        tabItem(
                          tabName = "visualizing_marker_expression",
                          selectInput("gene_select", "Select a gene:", choices = NULL),
                          textInput("gene_list_feature", "Selected genes for FeaturePlot:", value = ""),
                          actionButton("show_feature", "Display FeaturePlot"),
                          plotOutput("feature_plot"),
                          textInput("gene_list_vln", "Selected genes for VlnPlot:", value = ""),
                          actionButton("show_vln", "Display VlnPlot"),
                          checkboxInput("hide_vln_points", "Hide points on VlnPlot", value = FALSE),
                          plotOutput("vln_plot"),
                          textInput("gene_list_dotplot", "Selected genes for DotPlot:", value = ""),
                          actionButton("show_dot", "Display DotPlot"),
                          plotOutput("dot_plot")
                        )
                        ,

                        # Ajouter le contenu de l'onglet assigning_cell_type_identity
                        tabItem(
                          tabName = "assigning_cell_type_identity",
                          fluidRow(
                            box(  title = "Compares one cluster with all others", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,

                            column(width = 4,
                                   selectInput("select_cluster", "Select the cluster to rename:", choices = NULL),
                                   textInput("rename_single_cluster", "New name for the cluster:"),
                                   actionButton("rename_single_cluster_button", "Rename the selected cluster")
                            ),
                            column(width = 4,
                                   textInput("plot_title", "Plot title:", "UMAP Final"),
                                   numericInput("label_font_size", "Label font size", 5, min = 1, max = 20, step = 0.5),
                                   numericInput("pt_size", "Point size:", min = 0.1, max = 3, value = 0.3, step = 0.1)
                            ),
                            column(width = 4,
                                   selectInput("cluster_select", "Select Cluster:",  choices = NULL),
                                   colourInput("cluster_colour", "Select Colour:", value = "red"),
                                   actionButton("update_colour", "Update Colour")
                            )
                          )),
                          plotlyOutput("umap_finale")
                        ),


                        # Ajouter le contenu de l'onglet cluster_biomarkers
                        tabItem(
                          tabName = "cluster_biomarkers",
                          plotOutput("umap_cluster_10"),

                          box(title = "Compares one cluster with all others", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,

                          column(width = 4,
                          selectInput("cluster", "Cluster identifier:", choices = NULL),
                          sliderInput("number_genes_10", "Number of differential genes to display:",
                                      min = 10, max = 1000, value = 10, step = 50)
                          ),

                          column(width = 4,
                          numericInput("logfc_threshold_single","Log2 Fold Change threshold", value = 0.1, min = 0, max = 5, step = 0.1),
                          actionButton("find_markers", "Start analysis")
                          ),
                          column(width = 4,
                          numericInput("min_pct_single", "Minimum percentage threshold", value = 0.01, min = 0, max = 1, step = 0.01),
                          downloadButton("download_markers_single_cluster", "Download table (csv)")
                          ),
                          uiOutput("gene_tables_10")),
                          box(title = "Compares one cluster with one other cluster", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                          column(width = 4,

                          selectInput("cluster1", "Identifier of the first cluster:", choices = NULL),
                          selectInput("cluster2", "Second cluster identifier:", choices = NULL),
                          sliderInput("n_diff_markers", "Number of differential genes to display:",
                                      min = 10, max = 1000, value = 10, step = 50)
                          ),
                          column(width = 4,
                          numericInput("logfc_threshold_comparison", "Log2 Fold Change threshold:",value = 0.1, min = 0, max = 5, step = 0.1),
                          actionButton("compare_markers", "Start analysis")
                          ),

                          column(width = 4,
                          numericInput("min_pct_comparison",
                                       "Minimum percentage threshold",
                                       value = 0.01, min = 0, max = 1, step = 0.01),

                          downloadButton('download_markers_multiple_clusters', 'Download table (csv)')
                          ),
                          uiOutput("gene_tables_new")
                          ))

                        ,


                        # Ajouter le contenu de l'onglet subset
                        tabItem(
                          tabName = "subset",
                          plotOutput("global_umap"),
                          selectInput("select_ident_subset", "Select the clusters to include:", choices = NULL, multiple = TRUE),
                          textInput("gene_list", "Enter genes separated by comma (,)", value = ""),
                          numericInput("expression_threshold", "Expression threshold (default 0.1)", value = 0.1),
                          actionButton("apply_gene_subset", "Apply Gene Subsetting"),
                          numericInput("num_genes_to_express", "Number of genes to be expressed from the list:", value = 1, min = 1),
                          actionButton("apply_subset", "Apply the subset"),
                          plotOutput("subset_umap"),
                          downloadButton("download_subset_seurat", "Save subset as .RDS")


                        ),

                        #Merge de dataset
                        tabItem(
                          tabName = "merge_dataset",
                          p("To load your transcriptomics data, please compress the three 10X files (barcodes.tsv.gz, matrix.mtx.gz and features.tsv.gz) into a single .zip file. "),
                          # Adding a row to place the radio buttons side by side
                          fluidRow(
                            column(6,numericInput('num_datasets', 'Enter number of datasets to upload:', value = 2, min = 2)

                            ),
                            column(6,
                                   radioButtons("species_choice_merge", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "mouse")
                            )
                          ),
                          uiOutput("fileInputs"),
                          actionButton('integrate', 'Integrate'),
                          actionButton('add_field', 'Add Metadata Field'),
                          uiOutput("metadata_inputs"),
                          actionButton('add_metadata', 'Add Metadata'),
                          fileInput("load_seurat_file_merge", "Load Seurat object", accept = ".rds")

                        ),

                        #processing
                        tabItem(
                          tabName = "plot_merge",
                          actionButton("runScalePCA", "Scaling, PCA and et Elbow Plot"),
                          plotOutput("elbow_plot2"),
                          numericInput("dimension_2", "Number of dimensions:", value = 15, min = 1),
                          actionButton("runFindNeighbors", "Find Neighbors eand run UMAP"),
                          plotOutput("UMAPPlot"),
                          sliderInput("resolution_step2", "Resolution for clustering:", min = 0.01, max = 2, step = 0.01, value = 0.5),
                          actionButton("runFindClusters", "Find clusters"),
                          plotOutput("UMAPPlot_cluster")

                        )

                        ,



                        # Onglet Visualisation des gènes exprimés
                        tabItem(
                          tabName = "visualization_merge",

                          # Sélection du groupement
                          selectInput("group_by_select", "Group By:", choices = c("dataset", "cluster")),

                          # Sélection des gènes à visualiser
                          selectInput("geneInput_merge", "Select genes:", choices = NULL, multiple = TRUE, selectize = TRUE),

                          # TextInputs pour afficher les gènes sélectionnés
                          textInput("gene_list_vln_merge", "Selected genes for VlnPlot:", value = ""),

                          # Boutons pour exécuter les plots
                          actionButton("runFeaturePlot", "Run Feature Plot"),
                          plotOutput("FeaturePlot2"),
                          textInput("gene_list_feature_merge", "Selected genes for FeaturePlot:", value = ""),

                          actionButton("runVlnPlot", "Run Vln Plot"),
                          checkboxInput("hide_vln_points_merge", "Hide points on VlnPlot", value = FALSE),
                          plotOutput("VlnPlot2"),
                          textInput("gene_list_dot_merge", "Selected genes for DotPlot:", value = ""),

                          actionButton("runDotPlot", "Run Dot Plot"),
                          plotOutput("DotPlot2")
                        )
                        ,

                        tabItem(
                          tabName = "DE_merged_dataset",
                          actionButton("run_DE_merged", "Differentially expressed genes"),
                          downloadButton("save_seurat_merge", "Save Seurat Object"),
                          downloadButton('download_DE_merged', 'Download differentially expressed genes'),
                          numericInput("logfc_threshold_all_multiple", label = "Log2FC threshold", value = 0.25),
                          numericInput("min_pct_all_multiple", "Pct threshold ", value = 0.01, min = 0, max = 1, step = 0.01),
                          sliderInput("number_genes_merge", "Number of genes to display:",min = 1, max = 100, value = 10),
                          uiOutput("diff_genes_tables_merge"),
                          uiOutput("previous_tab_notification")
                        ),

                        # Onglet "assigning_cell_type_identity_merge"
                        tabItem(
                          tabName = "assigning_cell_type_identity_merge",
                          fluidRow(
                            column(width = 4,
                                   selectInput("select_cluster_merge", "Select Cluster:", choices = NULL),
                                   textInput("rename_single_cluster_merge", "New name for selected cluster:"),
                                   actionButton("rename_single_cluster_merge_button", "Rename the selected cluster")
                            ),
                            column(width = 4,
                                   textInput("plot_title_merge", "Plot title:", value = "UMAP Finale Merge"),
                                   numericInput("label_font_size_merge", "Label font size", value = 5, min = 1, max = 20, step = 0.5),
                                   numericInput("pt_size_merge", "Points size:", value = 0.1, min = 0.01, max = 5, step = 0.1)
                            ),
                            column(width = 4,
                                   selectInput("select_color_merge", "Select the cluster to be modified:", choices = NULL),
                                   colourInput("select_cluster_merge_color", "Choose a new color for the cluster:", value = "red"),
                                   actionButton("update_colour_merge_button", "Update cluster color")
                            )
                          ),
                          plotlyOutput("umap_finale_merge")
                        ),


                        # Onglet compare multiple clusters
                        tabItem(
                          tabName = "combined_analysis",
                          plotOutput("filtered_umap_plot"),
                          uiOutput("dataset_filter_ui"),
                              box(  title = "Compares one cluster with all others", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                            column(width = 4,
                                  selectInput("selected_cluster", label = "Select a cluster for comparison", choices = NULL),
                                  actionButton("calculate_DE", "Start analysis")
                                  ),
                            column(width = 4,
                                   numericInput("logfc_threshold_merge", label = "Log2 Fold Change threshold:", value = 0.25),
                                   downloadButton('download_markers_single_cluster_merge', 'Download differentially expressed genes')
                                   ),
                            column(width = 4,
                                  numericInput("min_pct_merge", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01)
                                   ),
                            DTOutput('DE_genes_table')
                            ),

                            box(  title = "Compares one cluster with one other cluster", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                              column(width = 4,
                                     uiOutput("cluster1_compare_ui"),
                                     numericInput("min_pct_compare_merge", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01),
                              ),
                              column(width = 4,
                                     uiOutput("cluster2_compare_ui"),

                                     actionButton("compare_clusters_button", "Start analysis between those Clusters")
                                     ),
                              column(width = 4,
                                     numericInput("logfc_threshold_compare_merge", label = "Log2 Fold Change threshold:", value = 0.25),

                                     downloadButton('download_markers_multiple_clusters_merge', 'Download differentially expressed genes')
                              ),
                            DTOutput("diff_genes_table_compare")
                            ),

                            box(title = "Compares a cluster between two datasets ", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                              column(width = 4,
                                     uiOutput("dataset1_compare_ui"),
                                     numericInput("logfc_threshold_datasets", label = "Log2 Fold Change threshold:", value = 0.25),
                                     actionButton("compare_datasets_button", "start analysis  between those datasets"),
                                     downloadButton("download_diff_dataset_cluster", "Download differentially expressed genes")

                                      ),
                              column(width = 4,
                                     uiOutput("dataset2_compare_ui"),
                                     numericInput("min_pct_compare_dataset_merge", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01),
                                     checkboxInput("all_clusters", "Compare all clusters", value = FALSE)
                              ),
                              column(width = 4,
                                     uiOutput("cluster_compare_ui")
                              ),
                              DTOutput("diff_dataset_cluster")
                            ),
                        )
                        ,

                        tabItem(
                          tabName = "subset_merge",
                          plotOutput("global_umap_merge"),
                          selectInput("select_ident_subset_merge", "Select the clusters to include:", choices = NULL, multiple = TRUE),
                          actionButton("apply_subset_merge", "Apply the subset"),
                          textInput("gene_list_merge", "Enter Gene Names (comma-separated)", value = ""),
                          numericInput("expression_threshold_merge", "Expression Threshold", value = 1, min = 0),
                          actionButton("apply_gene_subset_merge", "Apply Gene-based Subset"),
                          plotOutput("subset_umap_merge"),
                          downloadButton("download_subset_merge", "Save subset as .RDS")
                        ),

                      tabItem(
                        tabName = "atac_load",
                        fluidRow(
                          column(6, numericInput("num_datasets_atac", "Number of ATAC-seq datasets to load:", value = 1, min = 1)),
                          column(6, radioButtons("species_choice_atac", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "human"))
                        ),
                        uiOutput("fileInputs_atac"),
                        fileInput("load_seurat_file_atac", "Choose Seurat Object", accept = ".rds")

                      ),


                      tabItem(
                        tabName = "atac_preprocessing",
                        h2("QC Metrics for ATAC-seq Data"),
                        actionButton("calculate_qc_btn", "Calculate QC Metrics"),

                        fluidRow(
                          column(6, plotOutput("density_scatter_plot")),
                          column(6, plotOutput("tss_plot"))
                        ),
                        plotOutput("fragment_histogram_plot")
                      )


                      ,



                        tabItem(
                          tabName = "trajectory",
                          fluidRow(
                            # Colonne de gauche
                            column(width = 4,
                                   fileInput('load_seurat_file_monocle', 'Choose Seurat file')),

                            column(width = 4,
                                   actionButton("convertToMonocle", "Convert to Monocle")),

                            column(width = 4,
                                   actionButton("constructGraph", "Construct Graph"),
                                   verbatimTextOutput("selected_root_cell"),
                                   actionButton("startRootSelection", "Select root cell")
                            )),
                          numericInput("cluster_k", "Number of clusters (k):", value = 5, min = 1, step = 1),
                          plotlyOutput("trajectoryPlot"),
                          plotOutput("pseudotimePlot")
                        )
                        ,

                        tabItem(
                          tabName = "genes_pseudotime",
                          actionButton("diffGeneTestButton", "Run Differential Gene Test"),
                          DTOutput("diffGeneTable"),
                          downloadButton('download_pseudotime_diff_genes', 'Download pseudotime differential genes'),
                          textInput("geneTextInput", "Enter gene name"),
                          actionButton("geneSelectButton", "Visualize Gene Trajectory"),
                          plotOutput("geneTrajectoryPlot")
                        )
                        ,
                      tabItem(
                        tabName = "converter",
                        fileInput("file_seurat_conversion", "Choose Seurat RDS File", accept = c(".rds")),
                        actionButton("convert_button_anndata", "Convert to AnnData"),
                        downloadButton("download_anndata_converted_file", " Download AnnData object"),
                        actionButton("convert_button_monocle", "Convert to Monocle"),
                        downloadButton("download_monocle_converted_file", " Download Monocle object"),
                      )

                        ,

                        tabItem(
                          tabName = "acknowledgement",
                          fluidPage( tags$div(style = "position:relative; width:630px; height:800px;",
                                              tags$img(src = "muscle.png",
                                                       style = "position:absolute; top:0; left:0; width:100%; height:100%;"),
                                              tags$div("Acknowledgement & License",
                                                       style = "position:absolute; top:5%; left:30%; color:white; font-size:30px;"),
                                              tags$div("This application was developed by Gaspard Macaux and is the property of the Neuromuscular Development, Genetics and Physiopathology laboratory directed by Dr. Pascal Maire.",
                                                       style = "position:absolute; top:20%; left:4%; color:white; font-size:20px;"),
                                              tags$div("Many thanks to Edgar Jauliac, Léa Delivry and Hugues Escoffier for their advices and expertise in transcriptomic analysis",
                                                       style = "position:absolute; top:35%; left:4%; color:white; font-size:20px;"),
                                              tags$div("Many thanks to Seurat, who has built a very useful and well-documented library.",
                                                       style = "position:absolute; top:45%; left:4%; color:white; font-size:20px;"),
                                              tags$div("Many thankd to Shiny, which makes it so easy to develop graphical user interfaces.",
                                                       style = "position:absolute; top:55%; left:4%; color:white; font-size:20px;"),
                                              tags$div("Many thanks to Rstudio.",
                                                       style = "position:absolute; top:65%; left:4%; color:white; font-size:20px;"),
                                              tags$div("This application is licensed under the GPL3.",
                                                       style = "position:absolute; top:75%; left:4%; color:white; font-size:20px;"),
                          ),



                          ))
                      )



                      , tags$script(HTML('$(function () { $("[data-toggle=\'popover\']").popover(); });')))



)
