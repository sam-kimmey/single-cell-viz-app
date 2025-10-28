#
# Author: Sam Kimmey, PhD
# Date: July 15, 2024
#
# Purpose: This app was developed in order to viz data from general segmented
# data with QuPath.
# Instructions:

# Source file (Cmnd + Shift + Enter)
# Run app by simply entering "runApp('.')" into the console

# info on shiny plots: https://shiny.posit.co/r/articles/build/selecting-rows-of-data/

library(shiny)
library(data.table)
library(ggplot2)
library(plotly)

app.colors <- c(
  "light blue" = "#0f7d1cff",
  "vivid blue" = "#5fdf90ff",
  "dark blue" = "#0c7c53ff",
  "midnight blue" = "#0c462cff"
)
app.Ex.colors <- c(
  "light green 4" = "#5bc076ff",
  "light white warm" = "#ffffffff"
)

# Define UI for application that visualizes single-cell dataset generated from MIBI segmented data
# UI --------------
ui <- fluidPage(

  tags$head(
    # Note the wrapping of the string in HTML()
    # below tags$sytle enable some modifying of the shiny styling
    # below I updated the font
    tags$style(HTML("
        @import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,100..1000;1,9..40,100..1000&display=swap');
        body {
        background-color: white;
        color: #0c7c53ff;
        }
        
        h2 {
        font-family: 'DM Sans', sans-serif;
        font-weight: 700;
        }
        
        .shiny-input-container {
        color: #474747;
        }"))
  ),
    ## Application title ----
    titlePanel("MIBIviz2: An interactive single-cell segmentation map and expression visualization tool"),
    
    # Sidebar with a slider input for number of bins  ----
    sidebarLayout(
        sidebarPanel(
          ## logo ----          
          # img(src = "logo.png"),  # sidebar logo

          shinyjs::useShinyjs(), # to gray out color box with density
          
          ## input file -----
          fileInput("file", "Select .csv file", accept = ".csv"),
          
          ## Col entries & buttons -----
          uiOutput("columnSelectUI_X"),
            
          uiOutput("columnSelectUI_Y"),
            
          uiOutput("colorOverlaySelectUI"),
          
          radioButtons("density", "Exp or Density:",
                       c("Expression" = "exp",
                         "Density" = "dens")),
          
          actionButton("gate_info_toggle", "Gate information"), # TODO - Right now this button is non-functional - make into something useful or remove
          actionButton("exclude_reset", "Reset color"),
          actionButton("gate_name", "Gate name"), # added for naming the gate
          actionButton("clear_all_anno", "Clear Annotations"),
          
          ### Choose facet -----
          # this parameter is what will be displayed as multiple facets below the main plot
          # this is hard coded for any experimental groupings to be viz'd
          # first in concat list is default display
          selectInput("group", "Group:", 
                        c("slide view" = "slide.type",
                          "tile ROI view" = "tile.ROI",
                          "MIBIscope view" = "MIBIscope",
                          "Slide with ROI view" = "slide.with.ROI"
                          )),
          ## text output -----
          verbatimTextOutput("gateCoords"),
          downloadButton("download", class = "btn-block", label = "Download gate-annotated dataset as .csv")
          ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("biAxial1", brush = "plot1_brush", 
                      dblclick = "plot1_dblclick"), # brush object
           plotOutput("biAxial2", width = "800px", height = "1000px"), # can modify width and height depending on facet to group on lower plot
           # set to: plotOutput("biAxial2", width = "800px", height = "1000px") - in order to better view lower plot with multiple rows
           # tableOutput("data")
           )
        )
    )

#Server ----- 
# Define server logic required to draw a biaxial plot
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2) # allows larger file size import
  
  # reactivate value to store data.table
  data <- reactiveVal(NULL)

  # read in the data.table object
  ### input file ---------------
  observeEvent(input$file, {
    req(input$file)
    dt <- fread(input$file$datapath)
    
    # check for gateAnnotation col, and add new col for if necessary
    if(any(colnames(dt) %like% "gateAnnotation")){
      print("existing gateAnnotation column.")
    }else{
      dt[,gateAnnotation:="none"] # set to "none" instead of NA
    }
    data(dt)
    # selected <- reactiveVal(rep(FALSE, nrow(data()))) # added to try and iteritively add to selected cells
  })
  
  # X axis - default to centroid axis
  ### choose X axis ---------------
  output$columnSelectUI_X <- renderUI({
    req(data())
    selectInput("column_X", "Select X axis", choices = colnames(data()), 
                # selected = "HLA_1A_B_C") # can use for quick viz - will need to comment out below line if using this one
                selected = "Centroid.X.um")
  })# X axis
    
  ### choose Y axis ---------------
  output$columnSelectUI_Y <- renderUI({
    req(data())
    selectInput("column_Y", "Select Y axis", choices = colnames(data()), 
                # selected = "CD45") # can use for quick viz - will need to comment out below line if using this one
                selected = "Centroid.Y.um") # switch back to centroid for default
  })# Y axis
  
  ### choose color axis ---------------
  output$colorOverlaySelectUI <- renderUI({
    req(data())
    selectInput("column_color", "Select color axis", 
                choices = colnames(data()), 
                selected = "cell_label") # defaults to CD45
  })# Color axis
  
  ## Observe for density -----
  observeEvent(input$density, {
    if(input$density == "dens") {
      shinyjs::disable('colorOverlaySelectUI') 
    } else {
      shinyjs::enable('colorOverlaySelectUI')
    }
  }, ignoreNULL = T)
  
  ## Zoomable plot ---------
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  ### biaxial1 ggplot (top pane) ---------------
  output$biAxial1 <- renderPlot({
    req(data(), input$column_X, input$column_Y, input$column_color, input$density) # req data and coordinates to be loaded before plot appears
    rows.rand <- sample(nrow(data())) # randomized rows used for plotting
    selected_color_col <- input$column_color # string of selected col for color
    eval_color_data_type <- data()[[selected_color_col]][1] # extract first value for selected col to evaluate
    
    ### color scale evaluation, if statement to select color option for number of factor/string
    #### Logic for ggplot color -----
    if(is.numeric(eval_color_data_type)){
      # color for numeric (protein expression)
      expression_color_scale = scale_color_viridis_c(option = "magma", direction = 1)
    }else{
      # color for non-numerical column - i.e. factor colum (like slide type, MIBIscope, etc)
      expression_color_scale = scale_color_brewer(palette = "Dark2", direction = 1)
      # expression_color_scale = paletteer::scale_colour_paletteer_d("nationalparkcolors::Acadia") # alt pallet
      # this particular paletteer (above) includes some white - but this line is functional and can be the basis for different paletteer pallets
    }
    
    ### ggplot obj -----
    g <- ggplot(data()[rows.rand,], # data()[rows.rand,] - removing [rows.rand,] to check if that is leading to additional annotated cells in gate
           aes_string(x = input$column_X, # X and Y entered in by drop down
                      y = input$column_Y
                      )) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + # line for dynamic view/zoom of plot
      theme_minimal() + # theme
      labs(title = paste("Biaxial of", input$column_X, "x",input$column_Y)) # biaxial plot
    
    ### Exp/Density switch ----
    colors <- switch(
      input$density,
      "exp" = T,
      "dens" = F
    )
    
    # reverse centroid Y ordering to match expected orientation
    if(input$column_Y == "Centroid.Y.um"){g = g + scale_y_reverse()}
    
    if(isTRUE(colors)){ # if NOT plotting density
      g + 
        geom_point(alpha= 0.5, aes_string(color = input$column_color)) + # initial alpha
        geom_point( # geom point objects for those highlighted with the brush
          data = brushedPoints(data(), brush), # brush object created below
          alpha= 0.75, 
          color = app.colors["vivid blue"]) + # new color of cells that are highlighted
        labs(title = "Overlaid single-cell expression") + 
        expression_color_scale
    }else{ # else plotting density of data points
      g + 
        geom_point(alpha= 0.25, color = "black", ) +
        geom_density2d() +
        geom_density_2d_filled(alpha = 0.25, contour_var = "count") + # need to rename legend
        # geom_hex(bins = 70) +
        geom_point( # geom point objects for those highlighted with the brush
          data = brushedPoints(data(), brush), # brush object created below
          alpha= 0.75, 
          color = app.colors["vivid blue"]) + # new color of cells that are highlighted
        labs(title = "Cell density plot", fill = "Cells per contour")
    }
  })# biaxial ggplot end
  
  # Zoomable plot observation -------
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  ### biaxial1 brush ---------------
  brush <- NULL # brush object
  makeReactiveBinding("brush")
  
  observeEvent(input$plot1_brush, { # create reactive brush object
    brush <<- input$plot1_brush
  })
  
  ### Gate info button ---------------
  observeEvent(input$gate_info_toggle, {
    res <- brushedPoints(data(), input$plot1_brush, allRows = T)$selected_
  })# reset data to all true
  
  ### reset colored points button ---------------
  # eliminates colored points
  observeEvent(input$exclude_reset, {
    
    session$resetBrush("plot1_brush")
    # vals$keeprows <- rep(TRUE, nrow(data()))
    brush <<- NULL
  })# Reset all points
  
  ### reset Anno button ---------------
  # eliminates colored points
  observeEvent(input$clear_all_anno, {
    data()[,gateAnnotation:="NA"]
    # https://mastering-shiny.org/action-graphics.html - visit this link to help with coding this button to
    # remove annotation and update plot once that is done
    brush <<- NULL
  })# Reset all points
  
  ### Gate Name button ---------------
  l <- reactiveValues()
  observeEvent(input$gate_name, {
    # display a modal dialog with a header, textinput and action buttons
    showModal(modalDialog(
      tags$h2('Please enter gate name'),
      textInput('gatename', 'Gate Name'),
      footer=tagList(
        actionButton('submit', 'Submit name'),
        modalButton('cancel')
      )
    ))
  })# gate name

  # only store the information if the user clicks submit
  observeEvent(input$submit, {
    removeModal()
    toname <- brushedPoints(data(), input$plot1_brush, allRows = T)
    print(table(toname$selected_))
    l$name <- input$gatename
    print("length of named cells")
    print(nrow(toname))
    if(is.null(brush)){
      print("no cells selected")
    }else{
      print("naming cells...")
      print(l$name)
      ObjIDs <- toname[selected_ == TRUE,]$'cell_label'
      print("length of selected cells:")
      print(length(ObjIDs))
      data()[data()$'cell_label' %in% ObjIDs, gateAnnotation:= l$name]
      # session$resetBrush("plot1_brush")
      # vals$keeprows <- rep(TRUE, nrow(data()))
      brush <<- NULL
    }
  })
  
  ### displays selected cell info -----------
  output$gateCoords <- renderPrint({
    req(data(), input$column_X, input$column_Y)
    df <- brushedPoints(data(), input$plot1_brush, allRows = F)
    print(noquote(paste( "Total cells in gate:", nrow(df))))

    print(noquote("Count of cells in gate, for each condition:"))
    print(table(df$sample.Group))
    
  })# displays data table
  
  ### biaxial2 ggplot (bottom) ---------------
  output$biAxial2 <- renderPlot({
    # req data and coordinates to be loaded before plot appears
    req(data(), input$column_X, input$column_Y, input$column_color)
    rows.rand2 <- sample(nrow(data())) # randomized rows used for plotting
    ggplot(data()[rows.rand2,], 
           aes_string(
             x = "Centroid.X.um", # X and Y entered in by drop down
             y = "Centroid.Y.um")) + 
      geom_point(# mainplot style
        alpha= 0.5, 
        color = "black") + 
      geom_point( # geom point objects for those highlighted with the brush
        data = brushedPoints(data(), brush), # brush object created below
        size = 2,
        shape = 21,
        alpha= 0.75, 
        color = app.Ex.colors["light white warm"],
        fill = app.Ex.colors["light green 4"]) + # new color of cells that are highlighted
      theme_minimal() + # theme
      facet_wrap(~get(input$group), ncol = 2, ) +
      scale_y_reverse() + # reverse Y axis so the indexing matches default (counts from 0 at top left for Y axis)
      labs(title = paste("Cell centroid biaxial")) # biaxial plot
  })
  
  # Download -------------------------------------------------------
  # the Download .csv button will download the dataset with added annotations to the "/Downloads/" directory
  output$download <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$file$name),"gate_Annotated" ,".csv")
    },
    content = function(file) {
      fwrite(data(), file)
      # vroom::vroom_write(tidied(), file)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
