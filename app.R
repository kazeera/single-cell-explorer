# Load libraries
library(shiny)
library(Seurat)
library(viridis) #color palettes
library(cowplot) #plot_grid


# Color palettes
ui <- navbarPage("Khokha Lab MEC Explorer",
                 
                 # Tab 1
                 tabPanel("Human RM scRNA-seq",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            # Sidebar panel for inputs ----
                            sidebarPanel( 
                              # CSS
                              id="hmPanel", 
                              style="height:800px;",
                              helpText("Query human mammary epithelial cell (MEC) single cell RNA seq data (first published in Mahendralingam et al., Nature Metabolism, 2021)"),
                              
                              # Make drop down menu with discrete meta data variables
                              selectInput("sc_group.by", label = h5("Group by"), selected = "CellType", choices = c("CellType","CellStates","Seurat_clusters","Seurat_cellCyclePhase","Source","Patient")),# colnames(sc_data@meta.data)),
                              
                              # Make drop down menu with all gene names in single cell object
                              selectizeInput("sc_geneIDs", multiple = T, label = h5('Select genes'), choices = NULL, options = list(create = FALSE)), # if TRUE, allows newly created inputs),
                              
                              # Point size
                              numericInput("sc_pt.size", h5("Point size"), value = "0.1", min = 0.0001, max=10, step = 0.1, width = 100),
                              
                              # # Gradient color picker
                              # selectInput("sc_color_low", label = "Select color (low)", choices = names(sc_colors), selected = "Gray"),
                              # selectInput("sc_color_high", label = "Select color (high)", choices = names(sc_colors), selected = "Vermillion"),
                              
                              # Save data
                              
                              downloadButton('downloadData', 'Download Data'),
                              h5(""),
                              # Save plot
                              downloadButton(outputId = 'downloadPlot_sc', label = 'Download Plots'),
                              h5(""),
                              div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("file_height", h5("Plot height"), value = 10, min = 1, max=20, step = 1, width = 90)),
                              div(style="display: inline-block;vertical-align:top; width: 100px;",numericInput("file_width", h5("Plot width"), value = 10, min = 1, max=20, step = 1, width = 90)),
                              div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput("file_format_sc", label = h5("Save plots as"), choices = c("png", "jpeg", "svg", "pdf", "tiff"), selected = "jpeg"))
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              plotOutput(outputId = "sc", height = "800px") # sc = single cell plot
                            )
                          ))
)


# Define server logic ----
server <- function(input, output, session) {
  ## scRNA ----------------------------------------
  
  # Import data
  sc_data <- readRDS("data/sc_DietSeurat.rds")
  # For drop down menu with gene names, update and filter server-side to speed up program
  updateSelectizeInput(session, 'sc_geneIDs', choices = rownames(sc_data), selected = c("KRT18","VIM"), server = TRUE)

  # Umap for annotations
  umap1 <- reactive({
    DimPlot(sc_data, pt.size = input$sc_pt.size, reduction = "umapharmony", cols= turbo(length(unique(unlist(sc_data[[input$sc_group.by]])))), group.by = input$sc_group.by)
  })

  # Umap for features
  umap2 <- reactive({
    FeaturePlot(sc_data, features = input$sc_geneIDs, pt.size = input$sc_pt.size, reduction = "umapharmony", 
                cols = viridis(3))#, ncol = ifelse(length(input$sc_geneIDs) < 3, 2, 3)) #ifelse((length(input$sc_geneIDs)/4)<=1, length(input$sc_geneIDs), length(input$sc_geneIDs)/2))#c(sc_colors[[input$sc_color_low]], sc_colors[[input$sc_color_high]]),
  })

  # Stacked violin plot for features
  vln <- reactive({
    VlnPlot(sc_data, features = input$sc_geneIDs, stack = T, group.by = input$sc_group.by)
  })
  
  # Render output
  output$sc <- renderPlot({
    plot_grid(plot_grid(umap1(), vln(), ncol=2), umap2(), ncol=1)
  })
  
  # Download plot
  output$downloadPlot_sc <- downloadHandler(
    filename = function() { paste0('MEC_single_cell.', input$file_format_sc)},
    content = function(file) {
      # ggsave(filename, plot = plot_grid(umap1(), vln(), umap2(), sc_hm()), device = input$file_format_sc)#,  width = units(8.5, "in"), height = units(5, "in"), device = input$file_format_sc)
      ggsave(file, plot = plot_grid(plot_grid(umap1(), vln(), ncol=2), umap2(), ncol=1), 
             device = input$file_format_sc, height = input$file_width, width = input$file_width, units = "in")
    }
  )
  
  # Download csv of processed dataset
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("MEC_single_cell.csv")
    },
    
    content = function(file) {
      # Combine annotation with features
      df <- cbind(sc_data[[input$sc_group.by]], FetchData(sc_data, input$sc_geneIDs))
      # Save to file
      write.csv(df, file, row.names = TRUE)
    }
  )
}

# Run app 
shinyApp(ui = ui, server = server)
