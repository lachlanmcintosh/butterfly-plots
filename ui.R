rm(list=ls())
# update.packages(c("shiny"))
# install.packages("shiny")
library(shinydashboard)
library(shiny)
library("shinyBS")

# only need to edt this line in the whole script:
# dir = "/Users/lmcintosh/chris_CBSfits_example_copy/"

# full_files <- sapply(list.files(path=dir, pattern=".RData$"), function(x) paste(dir,x[[1]],sep=''))

CN_track_description <- "A copy number track shows the major (in colour) and minor (in black) copy numbers on the vertical axis against location on the horizontal axis. This plot is interactive with the butterfly plot in that you can highlight points (by dragging the mouse) on one plot and they will also be highlighted on the other."
bar_plot_description <- ""
data_table_description <- ""
density_plots_description <- ""
density_change_plots_description <- ""
gam_plots_description <- ""
subclone_barchart_description <- ""
Mahanalobis_control_plotui_description <- ""
GridIT_dialPlotui_description <- ""
butterfly_plot_description <- "A butterfly plot shows the major (vertical axis) and minor (horizontal axis) copy numbers for each location in the genome. This plot is interactive with the copy number track plot in that you can highlight points (by dragging the mouse) on one plot and they will also be highlighted on the other."


G=10
spread = pi/3
# Define UI for application
shinyUI(
  fluidPage(
  
    #  Application title
    titlePanel("Butterfly plots"),
    
    sidebarLayout(
      sidebarPanel(id = "tPanel",style = "overflow-y:scroll; height: 90vh;", width=4,
        bsCollapse(id = "collapseExample", open = "Input",
          ############################### INPUT ################################
          bsCollapsePanel("Input", "You have three options to select data for input to the application.",
                          style = "info",
            # INPUT PSCBS              
            div( 
              div(
                fileInput(inputId='file1', label='Option 1: Choose a segmented PSCBS File'),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("file1_q", label = "", icon = icon("question"), 
                         # always make the bsButton name correlated with the inputId
                  style = "info",size="extra-small"),
                bsPopover(id = "file1_q", title = "Load a PSCBS file",
                  content = paste0("Select a file which contains the entire PSCBS segmentation. For details on how to segment with PSCBS see the ", a("paper", href="https://www.ncbi.nlm.nih.gov/pubmed/21666266"),"."),
                  placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            # INPUT GRIDITH style
            div(
              fileInput(inputId='file2', label='Option 2: Choose gridith input CSV File'),
              style="width:92%;display:inline-block; vertical-align: middle;"
            ),
            div(
              style="display:inline-block; vertical-align: middle;",
              bsButton("file2_q", label = "", icon = icon("question"),
                style = "info",size="extra-small"),
              bsPopover(id = "file2_q", title = "GRIDITH file format.",
                content = "Select a file which contains a segmentation in the format as specified by GRIDITH. The file format is a csv delimited file with a header (,,,,,).",
                placement = "right", trigger = "click", options = list(container = "body")
              )
            ),
            # INPUT by folder
            div(
              textInput(inputId="folder", placeholder = "path to a folder containing segmentations", label = "Option 3: Select a folder of segmentations", value = NULL),
              checkboxInput("PSCBS_input", label = "Segmentations in this folder are of PSCBS format", value = TRUE),
              actionButton("next_file", "Go to the next file in the folder"),
              style="width:92%;display:inline-block; vertical-align: middle;"
            ),
            div(
              style="display:inline-block; vertical-align: middle;",
              bsButton("folder_q", label = "", icon = icon("question"),
                       style = "info",size="extra-small"),
              bsPopover(id = "folder_q", title = "GRIDITH file format.",
                        content = "Select a folder which only contains segmentations of interest.",
                        placement = "right", trigger = "click", options = list(container = "body")
              )
            )
          ),
          ############################### SELCT PLOTS ################################
          bsCollapsePanel("Select plots to display",style = "info",
            div( 
              div(
                checkboxInput("butterfly_plot", label = "Display the butterfly plot", value = TRUE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("butterfly_plot_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "butterfly_plot_q", title = "Butterfly plot",
                          content = butterfly_plot_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            div( 
              div(
                checkboxInput("track_plot", label = "Genome CN track", value = TRUE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("track_plot_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "track_plot_q", title = "Track plot",
                          content = CN_track_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            sliderInput("num_rows","Number of rows to split the CN track plot over", min=1,max=24,value=1,step=1),
            div( 
              div(
                checkboxInput("bar_plots", label = "CN bar plots", value = FALSE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("bar_plot_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "bar_plot_q", title = "CN bar plots",
                          content = bar_plot_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            div( 
              div(
                checkboxInput("density_plots", label = "CN density plots", value = FALSE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("density_plots_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "density_plots_q", title = "CN density plots",
                          content = density_plots_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            div( 
              div(
                checkboxInput("density_change_plots", label = "Adjacent CN change density plots", value = FALSE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("density_change_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "density_change_q", title = "Adjacent CN change density plots",
                          content = density_change_plots_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            div( 
              div(
                checkboxInput("gam_plots", label = "GAM Plots", value = FALSE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("gam_plots_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "gam_plots_q", title = "GAM plots",
                          content = gam_plots_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            div( 
              div(
                checkboxInput("subclone_barchart", label = "Subclonal clustering", value = FALSE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("subclone_barchart_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "subclone_barchart_q", title = "Subclonal clustering",
                          content = subclone_barchart_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            ),
            div( 
              div(
                checkboxInput("data_table", label = "Data table", value = FALSE),
                style="width:92%;display:inline-block; vertical-align: middle;"
              ),
              div(
                style="display:inline-block; vertical-align: middle;",
                bsButton("data_table_q", label = "", icon = icon("question"),
                         style = "info",size="extra-small"),
                bsPopover(id = "data_table_q", title = "Data table",
                          content = data_table_description,
                          placement = "right", trigger = "click", options = list(container = "body")
                )
              )
            )
          ),
          ############################### BUTTERFLY PLOTS ################################
          bsCollapsePanel("Butterfly plot options",style = "info",
                          checkboxInput("points", label = "Points", value = TRUE),
                          checkboxInput("error_bars", label = "Error bars", value = FALSE),
                          checkboxInput("lines", label = "Lines", value = TRUE),
                          checkboxInput("mirror_on", label = "Show mirrored points and lines", value = TRUE),
                          checkboxInput("selected_points", label = "Show the points that were selected with the shiny brush", value = TRUE),
                          checkboxInput("plot_clusters", label = "Plot clusters. Only do this after the dirichelet process clustering.", value = FALSE),
                          checkboxInput("save_plot", label = "Automatically save all butterfly plots as they are produced", value = FALSE)
          ),
          bsCollapsePanel("Butterfly plot contour options",style = "info",
            checkboxInput("contours", label = "Plot contours of a kernal density estimate (KDE)", value = FALSE),
            checkboxInput("contours_shade", label = "Shade KDE contours", value = FALSE),
            sliderInput("num_grid_points4contours", label = "Select the number of grid points to use when calculating the KDE", min = 50,max=1000,value = 200,step=1),
            sliderInput("bandwidth4contours", label = "Select the bandwidth of the kernal when calculating the KDE", min=0,max=1,value = 0.3,step=0.01)
          ),
          ############################### LIMIT WHAT YOU CAN SEE ################################
          bsCollapsePanel("Limit what you can see",style = "info",
            sliderInput("length", "segment length quantiles:",min = 0.0, max = 100, value =c(0,100),step=1),
            sliderInput("number", "segment snp number quantiles:",min = 0.0, max = 100, value =c(0,100),step=1),
            checkboxInput("equal_scales", label = "Equal scales", value = TRUE),
            actionButton("update_scale_max","zoom into 99% of the data"),
            sliderInput("xlimits", "xlimits:",min = -0.5, max = 10, value = c(0,10),step=0.01),
            sliderInput("ylimits", "ylimits:",min = -0.5, max = 10, value = c(0,10),step=0.01)
          ),
          ############################### PLOTTING AESTHETICS ################################
          bsCollapsePanel("Plotting aesthetics",style = "info",
            checkboxInput("same_sizes", label = "Different sized points", value = FALSE),
            sliderInput("size_points", "Size points:",
                        min = 0.0, max = 1, value = 1,step=0.01),
            sliderInput("size_lines", "Size lines:",
                        min = 0.0, max = 1, value = 1,step=0.01),
            sliderInput("size_density", "Size density:",
                        min = 0.0, max = 10, value = 1,step=0.01),
            sliderInput("shade_points", "Shade points:",
                        min = 0.0, max = 100, value = 10,step=1),      
            sliderInput("shade_selected_points", "Shade selected points:",
                        min = 0.0, max = 100, value = 50,step=1),
            sliderInput("shade_lines", "Shade lines:",
                        min = 0.0, max = 100, value = 10,step=1),
            sliderInput("user_shade", "Shade of sleected area:",
                        min = 0.0, max = 100, value = 20,step=1),
            sliderInput("grid_spacing", "Grid spacing:",
                        min = 0.0, max = 2.0, value = 0.25,step=0.01)   
          ),
          ############################### MAIN NORMALISATION PROCEDURES ################################
          bsCollapsePanel("Initial normalisation procedures",style = "info",
                          actionButton("shear_on_changes","Shear on 'change' density"),
                          actionButton("shear_on_density","Shear on density"),
                          sliderInput("shear_adjust", label = "Shear adjust", min = 0,max=2,value = 0.1,step=0.1),
                          actionButton("Do_IMBA","Do IMBA"),
                          actionButton("resize_by_change_ratio","Resize by change ratio"),
                          actionButton("resize_by_ab_ratio","Resize by AB ratio"),
                          #actionButton("slide_on_density","Slide on density"), # DEPRECATED
                          checkboxInput("draw_rect","Draw 11 to 20")
          ),
          ############################### KDE labelling of PARAMETERS ################################
          bsCollapsePanel("Kernal labelling paramters",style = "info",
                          checkboxInput("KDE", label = "Label clusters", value = FALSE),
                          sliderInput("num_grid_points", label = "KDE number of grid points", min = 50,max=1000,value = 200,step=1),
                          sliderInput("bandwidth", label = "KDE Bandwith", min=0,max=1,value = 0.3,step=0.01),
                          textInput("max_min_density_height", label = "Maximum minimum KDE height", value =0.1),
                          sliderInput("min_density_height", label = "Minimum KDE height", min=0,max=0.1,value = 0,step=0.01)
          ),
          bsCollapsePanel("Secondary normalisation procedures",style = "info",
                          actionButton("normalise_via_gam","nonlinear normalisation - GAM"),
                          actionButton("cluster","dirichelet process guasian mixture clustering"),
                          textInput("alpha", label = "alpha", value = 0.02),
                          sliderInput("sameness", label = "Prop adj same", min = 0,max=1,value = 0,step=0.01),
                          sliderInput("adjsame", label = "Prop adj same", min = 0,max=1,value = 0,step=0.01),
                          sliderInput("oneadjsame", label = "Prop one adj same", min = 0,max=1,value = 0,step=0.01),
                          sliderInput("allelicbal", label = "Prop allelic bal", min = 0,max=1,value = 0,step=0.01)
          ),
          ############################### MAHALANOBIS ################################
          bsCollapsePanel("Mahalanobis",style = "info",
            actionButton("Mahalanobis","Mahalanobis"),
            checkboxInput("Mahanalobis_control_plot", label = "Mahanalobis control plots", value = FALSE),
            h4("Mahalanobis parameters"),
            textInput("Mahalanobis_cut_top", label = "cut_top", value = 3.5),
            textInput("Mahalanobis_cut_bottom", label = "cut_bottom", value = 0.8),
            textInput("Mahalanobis_cut_right", label = "cut_right", value = 3.5),
            textInput("Mahalanobis_cut", label = "cut_distance", value = 7)
          ),
          ############################### OUTPUT ################################
          bsCollapsePanel("Output",style = "info",
            downloadButton('downloadImage', 'Download image'),
            downloadButton('downloadInfo', 'Download info'),
            downloadButton('downloadData', 'Download data')
          )
        ),
        ############################### OTHER ################################
        actionButton("redo","redo"),
        actionButton("undo","undo"),
        checkboxInput("automatically_update", label = "Auto update", value = TRUE),
        actionButton("accept","accept"),
        
        h3("Not sure where this stuff goes..."),
        actionButton("GridIt","modal Grid It"),
        #actionButton("GridIt","Grid It"), # inactive?
        #h3('Redundant'),
        #sliderInput("strechx", "strech x:",
        #            min = 0.0, max = 3.0, value = 1.0,step=0.01),
        #sliderInput("strechy", "strech y:",
        #            min = 0.0, max = 3.0, value = 1.0,step=0.01),
        #textInput("max_angle", label = "Max angle", value = pi/3),
        #textInput("max_ratio", label = "Max ratio", value = 3.0),
        #textInput("max_strech", label = "Max Stech", value = 3.0),
        #sliderInput("slider_bonus", "slider_bonus:",
        #            min=0, max=1, value=0.1,step=0.01),
        #sliderInput("density_adjust", "density adjust:",
        #            min=0, max=5, value=1,step=0.01),
        #sliderInput("angle1", "angle1:",
        #            min=-round(pi/spread*G,2), max=round(pi/spread*G,2), value=0,step=0.01),
        #sliderInput("angle2", "angle2:",
        #            min=-round(pi/spread*G,2), max=round(pi/spread*G,2), value=0,step=0.01),
        # h3("Focus on a genomic region"),
        textInput("genome_positions", label = "Genome positions, format: chr:(start,end),chr:pos,...", value = "")
        ############################### END SIDEBAR PANEL ################################
      ),
      
      mainPanel(
        ############################### MAIN PANEL START ################################
        # column(6,uiOutput("butterflyplotui"),uiOutput("trackplotui"),uiOutput("trackplot_histogramui"),uiOutput("densityplotui"),uiOutput("densitychangeplotui"),uiOutput("subclone_barchartui"),
        #        verbatimTextOutput("info"),
        #        dataTableOutput("contents")
               
        div(
          uiOutput("butterflyplotui"),
          # bsPopover(id = "butterflyplotui", title = "Load a PSCBS file",
          #           content = butterfly_plot_description,
          #           placement = "right", trigger = "click", options = list(container = "body")
          # )
          bsPopover(id = "butterflyplotui", title = "Load a PSCBS file",
                    content = butterfly_plot_description,
                    placement = "right", trigger = "click", options = list(container = "body")
          )
        ),
        # bsPopover(id="butterflyplotui_q", title="Butterfly plot", content=butterfly_plot_description, placement = "left", trigger = "hover", options = NULL),
        
        uiOutput("trackplotui"),
        #uiOutput("trackplot_histogramui"),
        uiOutput("densityplotui"),
        uiOutput("densitychangeplotui"),
        uiOutput("subclone_barchartui"),
        #uiOutput("Mahanalobis_control_plotui"),
        htmlOutput("ITH_res_text"),
        verbatimTextOutput("info"),
        textOutput("purity_out"),
        dataTableOutput("contents"),
        bsModal("modalExample", "Data Table", "GridIt", size = "large",
                actionButton("reset", "Reset"),
                uiOutput("GridIT_dialPlotui",
                         click = "plot_click",
                         brush = "plot_brush")
                )
        #      bsModal(id = "GridIt", title = "Grid It", trigger = "GridIt", size = "large",dataTableOutput("distTable"))      
        # bsModal("GridIt", "Grid It", "GridIt", size = "large",
        #         fluidRow(
        #           column(3, wellPanel(
        #             sliderInput('sampleSize', 'Sample Size', 
        #                         min=1, max=2, value=1, 
        #                         step=500, round=0)
        #             #uiOutput("GridIt_plot")
        #           )))
        # )      
    
        # modal_GridIt=modalDialog(
        #   fluidPage(
        #     h3(strong("GridIt"),align="center"),
        #     verbatimTextOutput("info")
        #     #uiOutput("GridIT_plotui")
        #   ),
        #   size="l"
        # )
        ############################### MAIN PANEL END ################################
        
      )
    )      # tab
  )
)