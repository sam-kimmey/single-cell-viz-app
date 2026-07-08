# Author: Sam Kimmey, PhD, Josh Kramer
# Create Date: July 15, 2024

# Purpose: This app was developed in order to viz single-cell datasets that are derived from segmented imaging data 
# and subsequent count extraction, followed by any necessary data transformation, normalization, cluster analysis, etc.

# This tool expects single-cell data to be saved in a .CSV format and in a file path that is on the computer used to 
# access this tool. 
# 
# Each columm in the loaded file can be visualized on the top biaxial plot. The user can select density, or to overlay a numeric or
# grouping column. Data can be gated to select, and those selected data points will be highlighted in cell centroid biaxial plot below.
# Qualities of the selected data points are displayed on the left. Double-clicking within the selected cell square will zoom into that region.
# Selected cells can be further annotated by clicking "Gate name" and typing in a name for the gate. A new column will be added to the loaded 
# dataset that can be saved by clicking "Download gate-annotated dataset as .csv" and chosing a name and file save location.
# 
# The Group setting can be used to change the facet-parameter for the centroid biaxial, allowing multiple ROIs or other grouping 
# variables to display multiple spatial plots.

# Source file (Cmnd + Shift + Enter)
# Run app by simply entering "runApp('.')" into the console

# info on shiny plots: https://shiny.posit.co/r/articles/build/selecting-rows-of-data/

# rsconnect::writeManifest(appPrimaryDoc = "app.R") # run this single line of code to update the manifest.json file after new lib added, otherwise comment out.

library(shiny)
library(data.table)
library(ggplot2)
library(plotly)
library(shinymanager)
library(tidyverse)
library(pals)
library(paletteer)
library(Polychrome)

app.colors = c(
  "light blue" = "#0f7d1cff",
  "vivid blue" = "#5fdf90ff",
  "dark blue" = "#0c7c53ff",
  "midnight blue" = "#0c462cff"
)
app.Ex.colors = c(
  "light green 4" = "#5bc076ff",
  "light white warm" = "#ffffffff"
)

deploy_msg = paste0("Last update: July 2026. Developed by Josh Kramer and Sam Kimmey, PhD")
# Define UI for application that visualizes single-cell dataset generated from MIBI segmented data
# UI --------------
ui = fluidPage(

  # absolutePanel(
  #   fixed = TRUE,
  #   bottom = 5,
  #   right = 10,
  #   width = 500,
  #   height = 20,
  #   tags$div(
  #     p(deploy_msg),
  #     style = "color: #aeaea0ff;font-size: 10px;"
  #   )
  # ),

  tags$head(
    # Note the wrapping of the string in HTML()
    # below tags$sytle enable some modifying of the shiny styling
    # below I updated the font
    tags$style(HTML("
        @import url('https://fonts.googleapis.com/css2?family=Montserrat:ital,wght@0,100..900;1,100..900&family=Oswald:wght@200..700&family=Rammetto+One&display=swap');        body {
        background-color: white;
        color: #0c7c53ff;
        }
        
        h2 {
        font-family: 'Rammetto One', sans-serif;
        font-weight: 700;
        }
        h3 {
        font-family: 'Montserrat', sans-serif;
        font-weight: 300;
        }        
        h5 {
        color: #000000ff;
        }
        .shiny-input-container {
        color: #474747;
        }"))
    ),

    tags$style(HTML("
    #top_right_image {
      position: absolute;
      top: 10px;
      right: 10px;
      z-index: 1000;
    }
    ")),

    tags$style(HTML("
    #top_right_image_left {
      position: absolute;
      top: 10px;
      right: 120px;
      z-index: 1000;
    }
    ")),

  # commented out below TO RESTORE IMG in upper right of plot page
  div(
    id = "top_right_image",
    tags$a(
      href = "https://mibiscope.com", # The destination URL
      target = "_blank",
      tags$img(
        # src = "logo.png", 
        src = "mibiscope_logo_dark.svg",
        height = "40px", 
        width = "100px",
        alt = "MIBIscope Logo Link" #https://www.oregon-physics.com
      )
    )
  ),

  div(
    id = "top_right_image_left",
    tags$a(
      href = "https://www.oregon-physics.com", # The destination URL
      target = "_blank",
      tags$img(
        src = "logo.png", 
        # src = "mibiscope_logo_dark.svg",
        height = "40px", 
        width = "100px",
        alt = "Oregon Physics Logo Link" #https://www.oregon-physics.com
      )
    )
  ),

    ## Application title ----
    titlePanel("CELLviz"),

    # subtitle using h3
    h3("an interactive spatial single-cell data visualization tool"),

    # Authors
    h6("Last update: July 2026. Developed by Josh Kramer and Sam Kimmey, PhD"),

    # Sidebar with a inputs for plot ----
    sidebarLayout(
        sidebarPanel(
          ## logo ----          
          # img(src = "logo.png"),  # sidebar logo

          shinyjs::useShinyjs(), # to gray out color box with density
          
          ## input file -----
          fileInput("file", "Select .csv file", accept = ".csv"),

          
          ## Col entries & buttons -----
          uiOutput("roi_selector"),

          uiOutput("facetWrap"),

          uiOutput("columnSelectUI_X"),
            
          uiOutput("columnSelectUI_Y"),
            
          uiOutput("colorOverlaySelectUItop"),

          uiOutput("subsetSelection"),
        
          uiOutput("colorOverlaySelectUIbottom"),
          
          radioButtons("density", "Expression or Density:",
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
          # selectInput("group", "Group:", 
          #               c("Slide Type" = "slide_type",
          #                 "ROI View" = "roi_id",
          #                 "MIBIscope view" = "mibi_instr",
          #                 "Phenotype View" = "phenotype",
          #                 "Cell Neighborhood View" = "neigh_kmeans"
          #                 )),
          downloadButton("download", class = "btn-block", label = "Download gate-annotated dataset as .csv"),
          ## text output -----
          h5("Selected data information:"), # eval labeling printed section

          verbatimTextOutput("gateCoords"),
          verbatimTextOutput("clickCoords")
          ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("biAxial1", 
           brush = "plot1_brush", 
           dblclick = "plot1_dblclick",
           click = "biax1_click"
          ), # brush object
           plotOutput("biAxial2", width = "1000px", height = "1000px") # can modify width and height depending on facet to group on lower plot
           # set to: plotOutput("biAxial2", width = "800px", height = "1000px") - in order to better view lower plot with multiple rows
           # tableOutput("data")
           )
        )
    )

make_creds = function(){
  users_raw = Sys.getenv("APP_ALLOWED_USERS", "")
  shared_pw = Sys.getenv("APP_SHARED_PWD", "")

  users = trimws(strsplit(users_raw, ",")[[1]])
  users = users[nzchar(users)]

  if (length(users) == 0 || shared_pw == ""){
    stop("Auth Vars not set")
  }

  ### User Cred ---------------
  data.frame(
    user = users,
    password = rep(shared_pw, length(users)),
    stringsAsFactors = FALSE
  )
}

# Server ----- 
# Define server logic required to draw a biaxial plot
server = function(input, output, session) {

  options("shinymanager.pwd_failure_limit" = 5) # allows larger file size import

  # user_creds = make_creds() # *** COMMENT TO RESTORE LOGIN ***

  ### Login check ---------------
  res_auth = secure_server(

    # check_credentials = check_credentials(user_creds) # *** COMMENT TO RESTORE LOGIN ***

  )
  
  options(shiny.maxRequestSize = 100*1024^2) # allows larger file size import - do not adjust
  
  # reactivate value to store data.table
  data = reactiveVal(NULL)

  # read in the data.table object
  ### input file ---------------
  observeEvent(input$file, {
    req(input$file)
    print(input$file$datapath)
    dt = fread(input$file$datapath)
    
    # check for gateAnnotation col, and add new col for if necessary
    if(any(colnames(dt) %like% "gateAnnotation")){
      print("existing gateAnnotation column.")
    }else{
      dt[,gateAnnotation:="none"] # set to "none" instead of NA
    }
    data(dt)
  })

  ### If there are mutliple ROIs, select ROI to view  ---------------
  output$roi_selector = renderUI({
    req(data())

    roi_choices = na.omit(unique(data()$roi_id))

    selectInput(
      "roi",
      "Select ROI to Visualize",
      choices = if (length(roi_choices) > 1) {
        c("All", roi_choices)
      } else {
        roi_choices
      },
      selected = (if (length(roi_choices) > 1) "All" else roi_choices)
    )
  })

  # If there are multiple ROIs in the dataset, allow the user to facet wrap the top plot by ROI
  output$facetWrap = renderUI({
    req(data())
    req(input$roi)

    if (length(unique(data()$roi_id)) > 1 && input$roi == "All") {
      radioButtons("facetWrap", "Facet Wrap Top Plot by ROI?", 
                  c("Yes" = "yes", 
                    "No" = "no"))
    } else {
      NULL
    }
  })
  
  # X axis - default to centroid axis
  ### choose X axis ---------------
  output$columnSelectUI_X = renderUI({
    req(data())
    selectInput("column_X", 
                "Select X axis", 
                choices = colnames(data()), 
                selected = "centroid_X_um")
  })# X axis
    
  ### choose Y axis ---------------
  output$columnSelectUI_Y = renderUI({
    req(data())
    selectInput("column_Y", 
                "Select Y axis", 
                choices = colnames(data()), 
                selected = "centroid_Y_um") # switch back to centroid for default
  })# Y axis
  # TESTING
  ### choose color axis for top plot ---------------
  output$colorOverlaySelectUItop = renderUI({
    req(data())
    selectInput("column_color_top", 
                "Select colors - Top Plot", 
                choices = colnames(data()), 
                selected = "cell_label") # defaults to CD45
  })# Color axis

  ### If top plot is colored by phenotype or neighborhood, allow user to subset to a specific one
  output$subsetSelection = renderUI({
    req(input$column_color_top)
    req(data())

    # If phenoype or neigh_kmeans, gives all options, or is NULL and does not appear
    choices = switch(
      input$column_color_top,
      "phenotype" = c(
        "All",
        unique(data()$phenotype)
      ),
      "neigh_kmeans" = c(
        "All",
        unique(data()$neigh_kmeans)
      ),
      NULL
    )

    req(choices)
    selectInput(
      "overlay_option",
      "Subset by phenotype or cell neighborhood",
      choices = choices
    )
  })
    
  ### choose color axis for bottom plot ---------------
  output$colorOverlaySelectUIbottom = renderUI({
    req(data())
    selectInput("column_color_bottom", 
                "Select colors - Bottom Plot", 
                choices = c("phenotype", "neigh_kmeans", "default"), 
                selected = "default")
  })# Color axis
  
  ## Observe for density -----
  observeEvent(input$density, {
    if(input$density == "dens") {
      shinyjs::disable('colorOverlaySelectUItop') 
    } else {
      shinyjs::enable('colorOverlaySelectUItop')
    }
  }, ignoreNULL = T)

  # Filter data by ROI if selected option is not "All"
  data_roi_filter = reactive({
    req(input$roi)

    if (input$roi != "All") {
      data() |>
        filter(roi_id == input$roi)
    } else {
      data()
    }
  })
  
  ## Zoomable plot ---------
  ranges = reactiveValues(x = NULL, y = NULL)
  
  ### biaxial1 ggplot (top pane) ---------------
  output$biAxial1 = renderPlot({
    ### Subset the if required
    data_filtered = data_roi_filter()

    if (!is.null(input$overlay_option) &&
        input$overlay_option != "All") {
      data_filtered = data_filtered |>
        filter(
          .data[[input$column_color_top]] == input$overlay_option
        )
    }

    req(data(), input$column_X, input$column_Y, input$column_color_top, input$density) # req data and coordinates to be loaded before plot appears
    rows.rand = sample(nrow(data_filtered)) # randomized rows used for plotting
    selected_color_col = input$column_color_top # string of selected col for color
    eval_color_data_type = data_filtered[[selected_color_col]][1] # extract first value for selected col to evaluate
    
    ### color scale evaluation, if statement to select color option for number of factor/string
    #### Logic for ggplot color -----
    if(is.numeric(eval_color_data_type)){
      # color for numeric (protein expression)
      expression_color_scale = scale_color_viridis_c(option = "magma", direction = 1)
    }else{
      # color for non-numerical column - i.e. factor colum (like slide type, MIBIscope, etc)
      # expression_color_scale = scale_color_brewer(palette = "Dark2", direction = 1) # only has 8 colors

      expression_color_scale = scale_color_paletteer_d("Polychrome::palette36") # may want to update to better pallete, OK for now.
      
    }
    
    ### ggplot obj -----
    g = ggplot(data_filtered[rows.rand,], # data()[rows.rand,] - removing [rows.rand,] to check if that is leading to additional annotated cells in gate
            do.call(aes, list(x = as.name(input$column_X), # X and Y entered in by drop down
                                          y = as.name(input$column_Y)
                                  ))) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) + # line for dynamic view/zoom of plot
      theme_bw()

    # Facet Wrap the top graph by ROI if user requested                   
    if (isTruthy(input$facetWrap) && input$facetWrap == "yes") {
      g = g + facet_wrap(~ roi_id, scales = "free")
    } else {
      g
    }

    # Title the top plot
    if (input$roi == "All") {
      g = g  + labs(title = paste("Biaxial of All ROIs with", input$column_X, "by", input$column_Y, "-", input$column_color_top))
    } else {
      g = g +labs(title = paste("Biaxial of", input$roi, "with", input$column_X, "by", input$column_Y, "-", input$column_color_top)) 
    }
    
    ### Exp/Density switch ----
    colors = switch(
      input$density,
      "exp" = T,
      "dens" = F
    )
    
    # reverse centroid Y ordering to match expected orientation
    # if(input$column_Y == "centroid_Y_um"){g = g + scale_y_reverse()}
    
    if(isTRUE(colors)){ # if NOT plotting density
      g + 
        geom_point(alpha= 0.5, do.call(aes, list(color = as.name(input$column_color_top)))) + # initial alpha
        geom_point( # geom point objects for those highlighted with the brush
          data = brushedPoints(data_filtered, brush), # brush object created below
          alpha= 0.75, 
          color = app.colors["vivid blue"]) + # new color of cells that are highlighted
          expression_color_scale
    }else{ # else plotting density of data points
      g + 
        geom_point(alpha= 0.25, color = "black") +
        geom_density2d() +
        geom_density_2d_filled(alpha = 0.25, contour_var = "count") + # need to rename legend
        geom_point( # geom point objects for those highlighted with the brush
          data = brushedPoints(data(), brush), # brush object created below
          alpha= 0.75, 
          color = app.colors["vivid blue"]) + # new color of cells that are highlighted
        labs(title = paste("Density plot", "-", input$column_color_top), fill = "Cells per contour")
    }
  })# biaxial ggplot end
  
  # Zoomable plot observation -------
  observeEvent(input$plot1_dblclick, {
    brush = input$plot1_brush
    if (!is.null(brush)) {
      ranges$x = c(brush$xmin, brush$xmax)
      ranges$y = c(brush$ymin, brush$ymax)
    } else {
      ranges$x = NULL
      ranges$y = NULL
    }
  })
  
  ### biaxial1 brush ---------------
  brush = NULL # brush object
  makeReactiveBinding("brush")
  
  observeEvent(input$plot1_brush, { # create reactive brush object
    brush <<- input$plot1_brush
  })
  
  ### Gate info button ---------------
  observeEvent(input$gate_info_toggle, {
    res <<- brushedPoints(data(), input$plot1_brush, allRows = T)$selected_
  })# reset data to all true
  
  ### reset colored points button ---------------
  # eliminates colored points
  observeEvent(input$exclude_reset, {
    session$resetBrush("plot1_brush")
    # vals$keeprows = rep(TRUE, nrow(data()))
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
  l = reactiveValues()
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
    toname = brushedPoints(data(), input$plot1_brush, allRows = T)
    print(table(toname$selected_))
    l$name = input$gatename
    print("length of named cells")
    print(nrow(toname))
    if(is.null(brush)){
      print("no cells selected")
    }else{
      print("naming cells...")
      print(l$name)
      ObjIDs = toname[selected_ == TRUE,]$'cell_label'
      print("length of selected cells:")
      print(length(ObjIDs))
      data()[data()$'cell_label' %in% ObjIDs, gateAnnotation:= l$name]
      brush <= NULL
    }
  })
  
  ### displays selected cell info -----------
  output$gateCoords = renderPrint({
    req(data(), input$column_X, input$column_Y)
    df = brushedPoints(data(), input$plot1_brush, allRows = FALSE)

    cat("Total cells in gate:", nrow(df), "\n\n")

    # column name -> display label
    summary_cols = c(
      roi_id         = "Count of cells per ROI",
      phenotype      = "Count of cell phenotypes",
      neigh_kmeans   = "Count of cell neighborhoods",
      sample_group1  = "Count of cells per Sample Grouping (1)",
      sample_group2  = "Count of cells per Sample Grouping (2)"
    )

    for (col in names(summary_cols)) {
      if (col %in% colnames(data())) {
        cat(summary_cols[[col]], ":\n", sep = "")
        print(table(df[[col]]))
        cat("\n")
      }
    }
  }, width = 50)

  ### displays click data info -----------
  output$clickCoords = renderPrint({
    req(input$biax1_click, input$column_X, input$column_Y)

    x_name = input$column_X
    y_name = input$column_Y

    cat("Click coordinates for:\n")
    cat(sprintf("  X (%s): %.2f\n", x_name, input$biax1_click$x))
    cat(sprintf("  Y (%s): %.2f\n", y_name, input$biax1_click$y))

  }, width = 50)
    
    output$biAxial2 = renderPlot({
      req(data(), input$column_color_bottom)
      
      color_col = input$column_color_bottom  
      req(color_col, color_col != "")          

      data_filtered2 = data_roi_filter()

      rows.rand2 = sample(nrow(data_filtered2))

      # If selection is "default" - everything is black, if not, it is scaled by color
      if (color_col == "default") {
        color_layer = geom_point(alpha = 0.5, color = "black")
        scale_layer = NULL
      } else {
        color_layer = geom_point(
          alpha = 0.5,
          aes(color = .data[[color_col]])
        )
        scale_layer = scale_color_paletteer_d("Polychrome::palette36")
      }

      ggplot(data_filtered2[rows.rand2, ],
            aes(x = .data[["centroid_X_um"]],
                y = .data[["centroid_Y_um"]])) +
        geom_point( # geom point objects for those highlighted with the brush
          data = brushedPoints(data_filtered2, brush), # brush object created below
          shape = "square") + # new shape of cells that are highlighted
        color_layer +
        scale_layer +
        theme_bw() +
        # theme(
        #   axis.text.x = element_blank(),
        #   axis.text.y = element_blank()
        # ) +
        # scale_y_reverse() +
        facet_wrap(~ roi_id, scales = "free") +
        labs(title = paste("ROI Visualization -", input$column_color_bottom), x= "Centroid (um)", y = "Centroid (um)")
    })
  
  # Download -------------------------------------------------------
  # the Download .csv button will download the dataset with added annotations to the "/Downloads/" directory
  output$download = downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$file$name),"gate_Annotated" ,".csv")
    },
    content = function(file) {
      fwrite(data(), file)
    }
  )
  }

# Run App -------------------------------------------------------

# Run the secure shiny app 
shinyApp( 
  ui = secure_app(ui,
  tags_top = tags$div(
    tags$h3("CELLviz"),
    tags$h6("Developed by Oregon Physics")
    
    # COMMENT TO RESTORE IMG - including photo seemed to slow down/lead to crash when loading data
    # tags$a(
    #   href = "https://www.oregon-physics.com", # The destination URL
    #   target = "_blank",
    #   tags$img(
    #     src = "logo.png", 
    #     height = "40px", 
    #     width = "100px",
    #     alt = "Oregon Physics Logo Link" #https://www.oregon-physics.com
    #   )
    # )

  )
), server = server)

# Run the application with no secure login
shinyApp(ui = ui, server = server) # *** COMMENT out TO RESTORE LOGIN ***
