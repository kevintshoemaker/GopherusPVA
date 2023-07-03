# SHINY APP: 'MITCH-clim 0.1'
#   Authors: Kevin Shoemaker (UNR), Elizabeth Hunter (USGS), Kevin Loope (VT), Michael Lam (UNR)
#   Funding provided by SERDP: RC18-1103 (Critical Habitat Breadth for Gopherus Tortoises)


# Step 1: load required packages (and install if not already installed)

# Load packages ---------------------

library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyFeedback)
library(popbio)
library(tidyverse)
library(raster)
library(terra)
library(leaflet)
library(sf)
library(fontawesome)
library(rgdal)

# Load functions ---------------------
make_compdf <- function(parr,ndx){
  temp <- apply(parr,c(1,2),sum)
  tempdf <- as.data.frame(t(apply(temp,2,function(t) quantile(t,c(0.05,0.5,0.95)))))
  tempdf[["scen"]] <- ndx
  names(tempdf) <- c("lb","median","ub","scenario")
  tempdf <- tempdf[,c("scenario","lb","median","ub")]
  tempdf[['year']] <- 0:nyears
  return(tempdf)
}

# Prepare the Workspace -----------------------------

nreps=20
nyears=89
nages=3

r <- raster("www/habmap.tif")  # should be raster to interface well with leaflet- don't change to terra yet
# v <- st_read("www/DESERT_FOCALAREA.shp")


# Read in list of scens and reps 

allscenarios<-readRDS("allscenarios_urls_A_and_B.rds") #use later for downloading data bundle
allscenarios$binscen <- sapply(1:nrow(allscenarios),function(t) paste(as.numeric(allscenarios[t,1:6]),collapse="") )
allscenarios$index <- as.integer(allscenarios$index) 
# allscenarios$index2 <- as.character(allscenarios$index)
allscenarios_short <- allscenarios[,1:8]

# Prepare workspace  

choicesnames <- names(allscenarios)[grepl("_clim",names(allscenarios))]
clim_activations <- data.frame(code=choicesnames)
clim_activations$fullnames <- c("Adult survival","Egg production","Probability of reproduction",
                                "Hatching success", "Sex ratio", "Somatic growth")
choiceslist <- 1:nrow(clim_activations)
names(choiceslist) = clim_activations$fullnames
choiceslist=as.list(choiceslist)

#LOADING DOD SHAPEFILE AND ACCESSING THE OBJECT IDS INSIDE

dod_sites = st_read("www/DOD_SITESINRANGE_F.shp")
# head(dod_sites)
dodSiteList = dod_sites$SITE_NAME

# history_list <- list() 

#Text strings 
instruction <- "To select a region of interest (ROI), you can \n (1) Draw a polygon to define the ROI on the map window to the right (one vertex generated per click), \n (2) Select a DoD installation from a preset list, or \n (3) upload a polygon from a shapefile (.shp)"
introtext <- "Conserving threatened and endangered species (TES) on Department of Defense (DoD) installations, without sacrificing vital military objectives, is necessarily a complex balancing act – especially given the uncertainties and risks associated with climate change."
objectivetext <- "This project advances the concept of critical habitat breadth as the foundation for rigorous TES conservation planning and vulnerability assessment. Gopherus tortoises - particularly Desert Tortoise (G. agassizii) and Gopher Tortoise (G. polyphemus) -- provide an excellent model system for applying the concept of critical habitat breadth because their populations have been extensively studied, and many prior translocations have been conducted which enable investigation into the inherent ability to acclimate to novel environments)."
methodtext <- "We compiled and re-analyzed previously collected data from across the range of both species to assess how demographic vital rates (e.g., survival, fecundity, age-at-maturity) respond to spatiotemporal environmental and climatic gradients. In addition, we conducted nesting surveys in the field to investigate how hatching success and hatchling sex ratios respond to climate gradients. We integrated these statistical models into comprehensive, spatially explicit predictive models capable of forecasting annual range-wide population dynamics for both species. We used this simulation model to forecast and visualize population dynamics through the year 2099, and to predict when and where populations are likely to be self-sustaining. We used machine-learning methods to assess which environmental and climatic drivers were most influential for determining population sustainability."
resulttext <- "Our simulation results for Desert Tortoise indicated that range-wide population abundance appears to be relatively stable. However, the southeastern portion of the Desert Tortoise range appeared to be the most heavily (and negatively) impacted by climate change, whereas habitat regions in the northern part of the species range appeared to respond more favorably to anticipated changes in climate.. Precipitation fluctuations led to strong temporal variance in projected reproductive output and population abundance. 
Our simulations indicated that most Gopher Tortoise populations appear to be declining, but that these declines may be driven more by low habitat quality rather than climate change per se. In particular, the observed decline in survival rate with years since fire strongly influenced population dynamics. Climate change, which is expected to reduce the frequency of prescribed fires across much of the southeastern United States, only exacerbated this effect in our simulations. Projected local extinction events in the Gopher Tortoise model tended to be most severe in the western and northern parts of the range, underscoring concerns about the status of Gopher Tortoise populations in this region.
"
benefittext <- "The results of our data collection linking climate to growth, reproduction, hatching success and sex ratio in both species will be directly beneficial to all future efforts to model climate effects on Gopherus populations.  Our range-wide viability analyses provide a tool to help resource managers to prioritize conservation efforts in regions with projected climate resiliency, to identify potential future translocation sites outside of DoD lands that will remain suitable as habitat for many decades into the future, and to direct resources toward ameliorating the most damaging effects of climate change, including the reduced potential for critical habitat management (e.g. prescribed fire). A decision support tool is currently in development to enable more effective information transfer to land managers. However, harnessing the value of this tool for management will require additional consultations and demonstrations with key stakeholders."


# UI ------------------
ui <- fluidPage( #theme=shinytheme("superhero"),title="Page title",
  #shinyFeedback::useShinyFeedback(), #future error code development
  titlePanel(h2("MITCH-Clim v0.1",
             tags$i( h4("Model-based Insights into Tortoise Critical Habitat under CLIMate change")))),
  navbarPage("",id="mynavbar",collapsible=TRUE, fluid=TRUE,
             #tabsetPanel(#"",
             #######navlistPanel(
             
             # Overview tab --------------
             
             tabPanel("Overview", #starts tabpanel1
                      #"THIS IS THE OVERVIEW", #for debugging 
                      div(style="color:darkgreen",
                          p(strong(h4("Purpose"))),), 
                      p(p("The purpose of this application is to assist resource managers in identifying mitigation and conservation strategies capable of sustaining viable metapopulations of Desert Tortoises under changing and uncertain future conditions. Defining the critical habitat breadth allows us to understand and predict the full range of habitat and environmental conditions capable of harboring viable populations now and in the future. This application allows the user to define a region of interest from which to display simulated abundance and critical habitat over time."),""),
                      p(p("We used these range-wide, spatially explicit population simulation models to produce projections of population growth rates and critical habitat over the next century, and to assess metapopulation viability under multiple plausible future climate scenarios. In this study, each cell represents a 5 X 5 km grid area.")),
                      p(p("The sections below provide detailed instructions on how to utilize this application.")),
                      div(p(""),br()),
                      
                      div(style="color:darkgreen",
                          p(strong(h4("General Application Flow"))),
                      ),
                      fluidRow(column(1),column(11,
                                                "1. At top of navigation panel, click on the 'Application' tab",br(), 
                                                "2. Select vital rates to activate",br(),
                                                "3. Select RCP climate scenario",br(), 
                                                "4. Load the scenario (vital rate x RCP combination)",br(),
                                                "5. Choose an ROI (region of interest)",br(),
                                                "6. Click on the 'Run Analysis' button (takes two minutes to process)",br(), 
                                                "7. Click on the new 'Output Viewer' tab to see results",br(), 
                                                "8. Load more scenarios (Application tab) to compare by rerunning analyses"
                      )),
                      div(p(""),br()),
                      
                      div(style="color:darkgreen",
                          p(strong(h4("Selecting Inputs (Step 1)"))),), 
                      p(p("The user must <check> which climate effect variables to include in the display of results. All CE variables are checked 'on' as the default setting. Once inputs are selected, click on the 'Load Scenario' button to begin loading. The loading process takes 60-90 seconds."),""),br(),
                      div(style="color:blue",
                          p(em(h4("Tortoise Climate Effect Variables: "))),
                      ),
                      p(strong("Six key vital rates modeled for fecudity and survival.")),
                      fluidRow(
                        column(1),
                        column(11,
                      p(strong("Adult survival"),": This parameter represents the annual survival rate for all reproductive adults (individuals above the age-at-maturity). Modeled using data from existing capture-mark-recapture and telemetry data, and published literature."),
                      p(strong("Egg production"),": This parameter was measured as the total annual per-capita egg production, also known as clutch size. Modeled using data from diagnostic imaging technology, and telemetry."),
                      p(strong("Probabability of reproduction"),": This parameter was drawn directly from our rangewide analysis of Desert Tortoise reproductive output. Modeled using data from diagnostic imaging technology, and telemetry."),
                      p(strong("Hatching success"),": This parameter represents the fraction of eggs within each nest that hatch successfully, assuming the nest does not fall victim to predation. Modeled using data from nest site monitoring."),
                      p(strong("Sex ratio"),": The sex ratio parameter (proportion of female hatchlings within each nest; PF) reflects the fraction of successful hatchlings from each nest that are female. Since this species has temperature-dependent sex determination, this parameter was modeled as a function of temperature, using data from nest site monitoring and hormone assays."),
                      p(strong("Somatic growth"),": This parameter represents growth by age as modeled using the von Bertalanffi growth function (VBGF), involving a K factor growth coefficient and an A factor asymptotic body size. Modeled using data from diagnostic imaging technology, and morphometrics."),
                        )#ends column
                      ),#ends fluidRow 
                      div(p(""),br()),
                      
                      div(style="color:blue",
                          p(em(h4("Climate Scenarios: "))),
                      ),
                      p("The Representative Concentration Pathways (RCPs) describe four different 21st century pathways of greenhouse gas (GHG) emissions and atmospheric concentrations, air pollutant emissions and land use."),
                      fluidRow(
                        column(1),
                        column(11,
                      p(strong("RCP2.6"),": This scenario aims to keep global warming likely below 2°C above pre-industrial temperatures. This is the most stringent of the scenarios."),
                      p(strong("RCP4.5"),": This intermediate scenario assumes a consistent decrease in emissions as a consequence of assumed air pollution control and GHG mitigation policy. More likely than not to exceed 2°C temperature increase by end of the 21st century."),
                      p(strong("RCP6.0"),": This intermediate scenario assumes a consistent decrease in emissions as a consequence of assumed air pollution control and GHG mitigation policy. Likely to exceed 2°C temperature increase by end of the 21st century."),
                      p(strong("RCP8.5"),": This 'baseline' scenario follows higher GHG emissions compared to the other three scenarios."),
                        )#ends column 
                      ),#ends fluidRow 
                      div(p(""),br()),
                      
                      
                      div(style="color:darkgreen",
                          p(strong(h4("Choosing a Region of Interest (Step 2)"))),
                      ),
                      p(em(strong("Drawing a polygon:"))),
                      fluidRow(column(1),column(11,
                      p("1. Click inside the Map box to initiate first vertex of polygon. If you make a wrong vertex at any point, click on 'Clear ROI' to start over.",br(),
                        "2. A vertex will be created where 'clicked' and a message a will appear.",br(),
                        "3. Click on a second location in the Map box to draw a second vertex. A blue line connecting these two points will now display.",br(),
                        "4. Continue to add vertices to your polygon until you have a least four vertices.",br(),
                        "5. Your polygon should now be a closed figure with the last line/vertex connected back to your first vertex."),)), 
                      p(em(strong("Selecting a polygon from dropdown menu:"))),
                      fluidRow(column(1),column(11,
                      p("1. Click on dropdown menu",br(),
                        "2. Select the desired preset DoD installation site (Note that this menu lists only the DoD sites that overlap with the desert tortoise range).",br(),
                        "3. The chosen DoD site polygon will automatically generate on the map.",br(), 
                        "4. Proceed to 'Run Analysis' (if scenario already downloaded) immediately. Clicking on Map will start a new polygon of interest."), 
                      )), 
                      p(em(strong("Uploading a polygon from file:"))),
                      fluidRow(column(1),column(11,
                      p("1. This must be an associated set of polygon files in the form of a shape files, with file extensions .dbf, .prj, .shp, and .shx.",br(),
                        "2. Select all four file components together and then click ok.",br(),
                        "3. A popup should then appear in a second to confirm that the file selected has been successfully downloaded.",br(),
                        "4. The polygon should automatically render onto the map.",br(),
                        "5. User will be able to download any current polygon to their PC by clicking the green 'Download Polygon' button under the Map box."), 
                      )),#end fluidRow 
                      #),
                      
                      div(p(""),br()),
                      
                      div(style="color:darkgreen",
                          p(strong(h4("Analysis and Output (Step 3)"))),
                      ),
                      p(p("Upon clicking the 'Run Analysis' button, a popup message on the bottom right should appear saying 'Processing scenario data...'. Please wait for the duration for the analysis to run (this takes two to three minutes). A second message will appear once the analysis is complete. A new tab called 'Output Viewer' will appear next to the 'Application' tab. Click on this tab. This tab contains several results panels which display the following: ")),
                      div(#style="color:black",
                      p(tags$ul(
                        tags$li("Current Run -- Plot of change in population abundance over time for ROI vs Full Range for the recent scenario run"),
                        tags$li("Age Structure -- Plot of population over time broken down by age structure class (hatchling, juvenile, and adult stages) of ROI vs Full Range of user-selected scenarios for comparison"),
                        tags$li("Comparison -- Plot output of ROI vs Full Range of user-selected scenarios for comparison"),
                        tags$li("Map -- Map output of population density ROI vs Full Range of user-selected scenarios for comparison"),
                      ),),
                      ),
                      #p("Interactive map gif (as shown in image) compiled from the 10-year interval maps."),
                      #p("Map displaying the net change in abundance from initial stage to end of simulation"),
                      #p("Map displaying the percent change in abundance from initial stage to end of simulation"),
                      div(p(""),br()),
                      
                      div(style="color:darkgreen",
                          p(strong(h4("Downloading Results"))),
                      ),
                      p(p("Users will be able to download the maps and/or plots individually to their own device as a Tif file. Users will also be able to download any current selected polygon ROI as a shapefile, which can be reuploaded in future app sessions to run more analyses."),""),
                      
                      div(p(""),br()),
                      
                      div(style="color:darkgreen",
                          p(strong(h4("Starting New Analyses"))),
                      ),
                      p(p(tags$ul(
                      tags$li("To begin a new analysis, click back to the 'Application' tab"),
                      tags$li("Click on the red 'Clear ROI' button located under the Map panel"),
                      tags$li("This will get rid of the polygon(s) as well as clear any plots, tables and maps displayed in the Output Viewer tab"),
                      ""))),
                      #p(em(tags$ul(
                      #  tags$li("To begin a new analysis, click back to the 'Application' tab"),
                      #  tags$li("Click on the red 'Clear ROI' button located under the Map panel"),
                      #  tags$li("This will get rid of the polygon(s) as well as clear any tables and maps displayed in the results area (Output Viewer tab)"),
                      #  ""))),
                      div(p(""),br()),
                      
             ), #ends tabpanel1
             
             # Application tab----------------
             
             tabPanel(id="tabpanel1", "Application", #starts tabpanel1   
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          h3("STEP 1: Select Scenario"),
                          # Copy the chunk below to make a group of checkboxes
                          checkboxGroupInput("climParams", label = h4("Select climate effects to be activated:"), 
                                             choices = choiceslist,
                                             selected = 1:length(choiceslist) ),
                          # textOutput("climParaTxt"), 
                          
                          selectInput("climScen", label = h4("Climate scenario"), c("rcp2.6","rcp4.5","rcp6.0","rcp8.5")),
                          # actionButton("generateROI", "Generate ROI Polygon"),
                          
                          actionButton("go_download", "Load Scenario"), # action variable called 'go_download' internally 
                          p(em("typical download time: 60-90 seconds")),
                          h3("________________",br(),br(),"STEP 2: Select ROI"),
                          # textOutput("instruction"),
                          p(em(instruction),""),
                          # select shapefile option 
                          selectizeInput(
                            'dodPolygon', h4("Choose DoD Installation"), choices = dodSiteList,
                            options = list(
                              placeholder = '(Optional) select installation',
                              onInitialize = I('function() { this.setValue(""); }')
                            )
                          ),
                          # selectInput("dodPolygon", label = h4("Preset polygons"), choices=dodSiteList, selected=NULL),
                          
                          # load shapefile 
                          fileInput(
                            inputId = "user_poly",
                            label = h4("Upload a polygon",
                                       tags$style(type = "text/css", "#q1 {vertical-align: top;}"),
                                       bsButton("q1", label = "", icon = icon("question"), style = "info", size = "extra-small")
                            ),
                            multiple = TRUE,
                            accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj")
                          ),
                          
                          bsPopover(id = "q1", title = "Detailed Instructions",
                                    content = paste0("<p>Click the browse button below and locate your shapefile, which should ",
                                                     "include one or more polygon(s) that represent your region of interest. ",
                                                     "Make sure to select all four objects of your ",
                                                     "shapefile (.dbf, .prj, .shp, and .shx). Then click Open at the ",
                                                     "bottom right on the window.</p>"),
                                    placement = "right", 
                                    trigger = "hover", 
                                    options = list(container = "body")),
                          h3("________________",br(),br(),"STEP 3: Run it!"),
                          actionButton("go_analysis", "Run Analysis"),br(), # action variable called 'go' internally
                          p(em("typical analysis runtime: 120-180 seconds")),
                          
                          textOutput("summaryText"), #,
                          # textOutput("Click_text") #debug
                          tags$div(class = "header", checked = NA,
                                   tags$br(),
                                   #tags$p("Michael was here >.< Kevinx2 was here -.-"),
                                   #tags$a(href = "https://external-preview.redd.it/SJ5Q4wfxekQPExJUCi4oF_NI2aWr0LR0nU6yLUmjz30.jpg?width=640&crop=smart&auto=webp&s=5aa466748bd3a5c5d6bf0bfd2dad9715dc7c6579", "Click Here!")
                          )
                          
                          
                        ),#end sidebarpanel
                        
                        #put main panel here 
                        mainPanel(id="mainpanel1",
                          # "THIS IS THE MAIN PANEL",#for debugging
                          #textOutput("Outputting text"),
                          leafletOutput("map"),
                          
                          actionButton("clearROI", "Clear ROI" ,
                                       icon("close"),
                                       style = "color: #fff; background-color: red; border-color: #fff;width:130;padding: 5px 5px 5px 5px;margin: 5px 5px 5px 5px; "),
                          
                          downloadButton('downloadPolygon',
                                         'Download Polygon',
                                         style = "color: #fff; background-color: #27ae60; border-color: #fff;padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; "),
                          
                         
                          # actionButton("clearROI", "Clear ROI"), # action variable called 'clearROI' internally 
                          # downloadButton('downloadPolygon', 'Download Polygon', style="display: block; margin: 0 auto; width: 230px;color: black;"),
                          
                          div(p(""),br()),
                          # div(p(""),br()),
                          
                          # "THIS IS THE MAIN PANEL",#for debugging
                          #hr(),
                          # fluidRow(column(3, verbatimTextOutput("value")))
                          tabsetPanel(id="ts1",type="tabs", # Outputs put into tabs
                                     tabPanel("Messages", 
                                        textOutput("main_message"))
                                      # tabPanel("Empty for now")
                          )
                        ),#end main panel
                      
                        #tabsetPanel(id="ts1",type="tabs", # Outputs put into tabs
                        #            tabPanel("Messages", 
                        #                     textOutput("main_message"))
                        #),
                          
                      ),#end sidebarlayout
             ),#end tabpanel1
             
             # Output Viewer tab ---------------
             
             tabPanel("Output Viewer",id="tabpanel2",#"Output Viewer", #starts tabpanel2
                      #"SECOND PART", #for debugging 
                      p(),
                      strong(),
                      em(),
                      br(),
                      
                      #div(style="color:blue",
                      #   p(strong(h4("Main Results: "))),),

                      #p(strong("Adult survival"),": This variable represents..."),
                      
                      tabsetPanel(id="tabsetPanel2", type="tabs", # Outputs put into tabs
                                  tabPanel("Current Run",fluid=TRUE,br(),
                                           splitLayout(strong(h4("Region of Interest"),style="color:indigo"),strong(h4("Full Range"),style="color:indigo")),
                                           splitLayout(
                                             plotOutput("trajec_plot_main"),
                                             plotOutput("trajec_plot_full") 
                                           )
                                  ),  # end 'current run' panel
                                  tabPanel("Age Structure",fluid=TRUE,br(),
                                           splitLayout(strong(h4("Region of Interest"),style="color:indigo"),strong(h4("Full Range"),style="color:indigo")),
                                           splitLayout(
                                             plotOutput("trajec_plot_main_age"),
                                             plotOutput("trajec_plot_full_age") 
                                           )
                                  ),  # end 'age structure' panel
                                  #tabPanel("Full Range",fluid=TRUE,),
                                  tabPanel("Comparisons",fluid=TRUE,br(),
                                           fluidRow(
                                             column(4,p("TRUE = climate effect on parameter turned on",br(),"FALSE = no climate effect on parameter")),
                                             #column(1,checkboxInput('ones', 'Ones'),
                                             #       checkboxInput('twos', 'Twos')),
                                             column(8,tableOutput("history_table")),
                                           ),
                                           #splitLayout(cellWidths = c("25%", "75%"),
                                          #             p("TRUE = climate effect on parameter turned on",br(),"FALSE = no climate effect on parameter"),
                                          #             tableOutput("history_table") ), 
                                           splitLayout(
                                             strong(h4("Region of Interest"),style="color:indigo"),strong(h4("Full Range"),style="color:indigo")
                                           ),
                                           splitLayout(
                                             checkboxGroupInput("compare_scen_clipped", label = h5("Select scenario(s) to display:"), 
                                                                choices = NULL,
                                                                ),
                                             checkboxGroupInput("compare_scen", label = h5("Select scenario(s) to display:"), 
                                                              choices = NULL,)
                                           ),
                                           splitLayout(
                                             plotOutput("trajec_plot_compare_clipped"),
                                             plotOutput("trajec_plot_compare"),
                                           ),
                                           # fluidRow(
                                           #   div(style="color:purple", #displays on right side of screen 
                                           #       p(strong(h4("Scenarios: "))),
                                           #   ),
                                           #   p(strong("RCP26"),": This scenario represents..."),
                                           #   #),
                                           # )   
                                           ),  #end tab panel
                                  # tabPanel("Generated Scenarios",fluid=TRUE,
                                  #          # fluidRow(
                                  #          # "Display text here for now",
                                  #          # tableOutput("history_table"),),
                                  #          # fluidRow(
                                  #          #   div(style="color:purple", #displays on right side of screen 
                                  #          #       p(strong(h4("Scenarios: "))),
                                  #          #   ),
                                  #          #   p(strong("RCP26"),": This scenario represents..."),
                                  #          #   p(strong("RCP45"),": This scenario represents..."),
                                  #          #   p(strong("RCP60"),": This scenario represents..."),
                                  #          #   p(strong("RCP85"),": This scenario represents..."),
                                  #          # #),
                                  #          # )   
                                  #          ),#end tabpanel
                              
                                    tabPanel("Maps",fluid=TRUE,
                                             sliderInput("map_year", "Simulation Year", min=0, 
                                                         max=nyears, value=0, 
                                                         step = 1, round = TRUE, ticks = TRUE, animate = animationOptions(interval=4000)),
                                             sliderInput("map_percentile", "Percentile to Display", min=0, 
                                                         max=100, value=50, 
                                                         step = 1, round = TRUE, ticks = TRUE, animate = FALSE),
                                             splitLayout(strong(h4("Region of Interest"),style="color:indigo"),strong(h4("Full Range"),style="color:indigo")),
                                             splitLayout(
                                               plotOutput("map_plot_main"),
                                               plotOutput("map_plot_full") 
                                             )
                                    )#end tabpanel 'Maps'
                                  # tabPanel("Vital Rates",fluid=TRUE,
                                  #          
                                  #          leafletOutput("map2"),
                                  #          p("Click on a pixel to view change in vital rates over time"),
                                  #          fluidRow(
                                  #            div(style="color:purple", #displays on right side of screen 
                                  #                p(strong(h4("Scenarios: "))),
                                  #            ),
                                  #            p(strong("RCP26"),": This scenario represents..."),
                                  #            p(strong("RCP45"),": This scenario represents..."),
                                  #            p(strong("RCP60"),": This scenario represents..."),
                                  #            p(strong("RCP85"),": This scenario represents..."),
                                  #            #),
                                  #          )   
                                  # )#end tabpanel
                      ),
                      #),
                      div(p(""),br()),
                      
                      # fluidRow(
                      # div(style="color:purple", #displays on right side of screen 
                      #     p(strong(h4("Climate Scenarios: "))),
                      # ),
                      # p(strong("RCP26"),": This scenario represents..."),
                      # p(strong("RCP45"),": This scenario represents..."),
                      # p(strong("RCP60"),": This scenario represents..."),
                      # p(strong("RCP85"),": This scenario represents..."),
                      # ),
                      
                      div(p(""),br()),
                      #textOutput("Hello"),
             ), #ends tabpanel3
             
             
             tabPanel("Figures", #starts tabpanel4
                      #"FOURTH PART", #for debugging 
                      div(#style="color:blue",
                        p(strong(h4("Below are a selection of figures highlighting results of the general desert tortoise range simulation model.")),br()),
                      ),
                      
                      div(align="center",
                          p(strong("Figure 1. "),"Schematic of population modeling for determining population resilience and critical habitat."),img(src='FIG_ES_1.png',height="50%",width="50%"),br(),br(),br(),
                          p(strong("Figure 2. "),"Comparison of Final Abundances between RCP Climate Scenarios"),img(src='FA_COMP2.jpg',height="50%",width="50%"),br(),br(),br(),
                          #p(strong("Figure 2. "),"Difference in Change in Abundance between the Two Most Likely Climate Scenarios"),img(src='ABUNDANCECHANGE_DT_60MINUS45.png',height="50%",width="50%"),br(),br(),br(), #renderImage(all_im[[1]]),#imageOutput(all_im[[1]]), 
                          p(strong("Figure 3. "),"Comparison of Critical Habitat Breadth between RCP Climate Scenarios"),img(src='LAMBDA_COMP2.jpg',height="50%",width="50%"),br(),br(),br(),#renderPlot(all_im[[1]]),#img(src=all_im[[1]],height="25%",width="25%"),
                          #p(strong("Figure 4. "),"Comparison of Critical Habitat Breadth between Best and Worst Case Scenarios"),img(src='60v85_LAMBDA2.jpg',height="75%",width="75%"),br(),br(),br(),
                          p(strong("Figure 5. "),"Critical Habitat Breadth for Averaged Likely Scenario 'RCP4.5, RCP6.0, RCP8.5' Overlay with DoD Sites"),img(src='AVG456085LAMBDA3DoD.png',height="50%",width="50%"),br(),br(),br(),#img(src=img_files[[1]]),
                          p(strong("Figure 6. "),"Comparison of Clim Effect vs No Clim Effect between RCP Climate Scenarios"),img(src='C_O_NC_COMPPINK.jpg',height="50%",width="50%"),br(),br(),br(), 
                          p(strong("Figure 7. "),"Clim Effect for Averaged Likely Scenario 'RCP4.5, RCP6.0, RCP8.5'"),img(src='AVG456085C_O_NCPINK.png',height="50%",width="50%"),br(),br(),br(),
                      ),
                      
                      div(p(""),br()),
             ), #ends tabpanel4
             
             tabPanel("About", #starts tabpanel5
                      #"THIRD PART", #for debugging 
                      div(style="color:darkgreen",
                          p(strong(h4("Introduction"))),
                      ),
                      p(em(introtext),""),
                      div(p(""),br()),
                      div(style="color:darkgreen",
                          p(strong(h4("Objectives"))),
                      ),
                      p(em(objectivetext),""),
                      div(p(""),br()),
                      div(style="color:darkgreen",
                          p(strong(h4("Technical Approach"))),
                      ),
                      p(em(methodtext),""),
                      div(p(""),br()),
                      div(style="color:darkgreen",
                          p(strong(h4("Results"))),
                      ),
                      p(em(resulttext),""),
                      div(p(""),br()),
                      div(style="color:darkgreen",
                          p(strong(h4("Benefits"))),
                      ),
                      p(em(benefittext),""),
                      div(p(""),br()),
             ), #ends tabpanel5
             
             tabPanel("More Info", #starts tabpanel6
                      #"SIxTH PART", #for debugging 
                      div(#style="color:blue",
                        p(strong(h4("For more information, please see: "))),
                      ),
                      p("Paper 1"),
                      p("Paper 2"),
                      p("Paper 3"),
                      p("Paper 4"),
                      p("Paper 5"),
                      p("Paper 6"),br(),
                      p("The code for this application as well as the underlying models are available on GitHub",tags$a(href="https://github.com/kevintshoemaker/GopherusPVA", " HERE!")),
                      div(p(""),br()),
                      div(style="font-size:12px",
                          p("This research is supported by the Strategic Environmental Research and Development Program (grant award #RC18-1103). Research was authorized under the following permits: UNR IACUC #00703, USFWS Recovery Permits TE-030659 and TE-50049D, NV-NDOW permit S36421, and GSU IACUC #19007, GA DNR permit 1000838720, FL FWC permit LSSC-18-00023B."),br()),
             ) , #ends tabpanel6
             #)
             
  ), #end narbarpage
  
  navbarMenu(Title="Project",
             tabPanel("End"),
             "(c) 2023, Kevin T Shoemaker & Team. All rights reserved.")

) #ends fluidpage


# SERVER-------------------
server <- function(input, output, session) {
  
  hideTab(inputId="mynavbar", target="Output Viewer") #"Compare"   # hide tab on start
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # upload and plot the DOD installation shapefile 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  dod_polygon <- reactive({
    req(input$dodPolygon)
    
    poly <- dod_sites[dod_sites$SITE_NAME==input$dodPolygon,]
    
    poly <- st_transform(poly, crs=4326 )   # CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
    
    # poly <- st_cast(poly,"POLYGON")
    
    # col <- SpatialPolygonsDataFrame(col, data.frame(IDs=1:length(col@polygons)))
    # col@data$col <- sample(c('purple','orange','pink','yellow','green',"blue","red"), nrow(col@data),1)
    return(poly)
    
  })
  
  # observer that adds user uploaded polygon if available and saves it to the polygon object
  observe({
    map_proxy = leafletProxy("map",data=dod_polygon()) %>%
      clearPopups() %>% 
      clearShapes() %>% 
      addPolygons()
    
    myReactives$polygon <- dod_polygon()
  })
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # upload and plot the user defined ROI shapefile 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  user_polygon <- reactive({
    req(input$user_poly)
    
    # shpdf is a data.frame with the name, size, type and
    # datapath of the uploaded files
    shpdf <- input$user_poly
    # Name of the temporary directory where files are uploaded
    tempdirname <- dirname(shpdf$datapath[1])
    # Rename files
    for (i in 1:nrow(shpdf)) {
      file.rename(
        shpdf$datapath[i],
        paste0(tempdirname, "/", shpdf$name[i])
      )
    }
    
    # read-in user defined shapefile
    poly <- st_read(paste(tempdirname,
                         shpdf$name[grep(pattern = "*.shp$", shpdf$name)],
                         sep = "/"
    ))
    
    poly <- st_transform(poly, crs=4326 )   # CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
    
    # poly <- st_cast(poly,"POLYGON")
    
    # col <- SpatialPolygonsDataFrame(col, data.frame(IDs=1:length(col@polygons)))
    # col@data$col <- sample(c('purple','orange','pink','yellow','green',"blue","red"), nrow(col@data),1)
    return(poly)
    
  })
  
  # observer that adds user uploaded polygon if available and saves it to the polygon object
  observe({
    map_proxy = leafletProxy("map",data=user_polygon()) %>%
      clearPopups() %>% 
      clearShapes() %>% 
      addPolygons()
    
    myReactives$polygon <- user_polygon()
  })

  
  addPopover(session, id="map", title="Interpretation of Map", content = paste0("<p>This map reflects the range of the desert tortoise. ",
                                                                       "Colors represent habitat suitability (blue is more suitable, red is less suitable.",
                                                                       " Polygons may be uploaded, selected from a dropdown list (DoD sites), ",
                                                                       "or drawn directly on the map. Click the 'Run Analysis' button",
                                                                       " when you are satisfied with the shapefile. If there is no ",
                                                                       "shapefile, then the app will run analyses for the ",
                                                                       "entire range.</p>"), trigger = 'hover')
  
  # END SHAPEFILE LOADING 
  
  
  
  
  
  
  
  output$climParaTxt <- renderText({input$climParams})
  
  myReactives <- reactiveValues(numClicks=0)  # make a list for storing reactive values... 
  # output$value <- renderPrint({ input$checkGroup })
  myReactives$scensRun <- NULL   # may not be needed, but can't hurt
  
  # static map elements
  output$map <- renderLeaflet({
    map <- leaflet() %>% 
      addTiles %>% 
      addRasterImage(r,opacity=0.4) %>% 
      # addPolygons() %>% 
      setView(lng=-116.84725141459298, lat=34.96726660170186, zoom=6)
  })
  
  #output$map2 <- renderLeaflet({
  #   map2 <- leaflet() %>% 
  #    addTiles %>% 
  #    addRasterImage(r,opacity=0.4) %>% 
  #    # addPolygons() %>% 
  #    setView(lng=-116.84725141459298, lat=34.96726660170186, zoom=6)
  #})
  
  observeEvent(input$map_click, {
    click = input$map_click
    myReactives$numClicks <- myReactives$numClicks + 1
    text<-paste("Latitude: ", click$lat, ", Longtitude: ", click$lng)
    text2<-paste("You've selected point ", text)
    text3 <- "Make as many polygon vertices as you'd like. When you're ready, click *run analysis*"
    myReactives$thisROI <- rbind(myReactives$thisROI,click) #creates the object region of interest ROI
    df <- as.data.frame(myReactives$thisROI[,c(1,2)])
    temp <- lapply(1:ncol(df),function(t) df[[t]] <<- as.numeric(df[[t]])) #was producing error
    output$ROItable <- renderTable({myReactives$thisROI})
    
    map_proxy = leafletProxy("map") %>%
      clearPopups() 
    
    # output$Click_text<-renderText({
    #   text2
    # })
    
    # write.csv(myReactives$thisROI,"thisROI.csv", row.names = F)
    
    if(nrow(df)>3){
      polygon <- df %>%
        st_as_sf(coords = c("lng", "lat"), crs = 4326) %>%
        summarise(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON")
      myReactives$polygon <- polygon
      map_proxy = leafletProxy("map",data=polygon) %>%
          clearShapes() %>% 
          addPolygons()
    }else if(nrow(df)>1){
      line <- df[, c("lng","lat")] %>%
        as.matrix() %>%
        st_linestring() %>% #converts to line segments
        st_sfc() %>% #formats above
        st_sf() #formats above to a "polyline" object 
      map_proxy = leafletProxy("map",data=line) %>%
        addPolylines () #draws ACTUAL line 
    }else if(nrow(df)==1){
      map_proxy = leafletProxy("map",data=line) %>%
        addPopups(click$lng, click$lat, text3)
    }
    
  })
  
  observeEvent(input$clearROI, {
    myReactives$thisROI <- NULL
    myReactives$polygon <- NULL
    myReactives$numClicks <- 0
    leafletProxy("map",session) %>%
      clearShapes() 
  })
  
  observeEvent(input$tabpanel2_click, {
    if(length(input$compare_scen)==0){
      showNotification("No results will display since no analysis has been run yet.")
    }
  })
  
  # observeEvent(input$map_draw_new_feature, {
  #   # leafletProxy("map") %>%
  #   #   removeDrawToolbar(clearFeatures = T)
  #   
  #   feat <- input$map_draw_new_feature
  #   coords <- unlist(feat$geometry$coordinates)
  #   coords <- matrix(coords, ncol = 2, byrow = T)
  #   poly <- st_sf(st_sfc(st_polygon(list(coords))), crs = 4326)
  #   print(st_bbox(poly))
  # })
  
  # output$instruction <- renderText({
  #   "Draw a polygon to define the region of interest (ROI). \nEach click on the map will generate one vertex."
  # })
  
  thisclim <- reactive({
    gsub("\\.","", input$climScen)
  })
  
  scenndx <- reactive({
    climParmsSelected <- clim_activations$code[as.numeric(input$climParams)]
    binps <- paste(as.numeric(names(allscenarios)[1:6]%in%climParmsSelected),collapse="") #converts selection to T/F list
    # myReactives$thisclim <- gsub("\\.","", input$climScen)
    which(allscenarios$binscen==binps&allscenarios$clim_scenario==thisclim() ) # scenario index
    # scenndx = 1
  })
  
  my_dl_list <- reactive({
    file.path(tempdir(),allscenarios$name[[scenndx()]])
  })  # #tempdir()
  
     # myrasts_clipped: this is a recipe for creating a clipped version of myrasts... 
  myrasts_clipped <- reactive({
    req(myReactives$polygon)
    theserasts <- myrasts()
    thisclip <- lapply(theserasts, function(t) terra::mask(t, thispoly()) )
    return(thisclip)
  })
  
  # myrasts: recipe for reading downloaded rasters using terra::rast
  myrasts<- reactive({
    if(all(file.exists(my_dl_list()))){
      lapply(my_dl_list(),rast)
    }else{
      NULL
    }
    
  })  #convert list to rasters (read in the downloaded tifs)
  
     # poparray: this is a recipe for creating an unclipped version of poparray
  poparray <- reactive({
    # if(file.exists(sprintf("%s/poparray_%s.rds",tempdir(),scenndx()))){
    #   thisarray <- readRDS(sprintf("%s/poparray_%s.rds",tempdir(),scenndx()))
    #   return(thisarray)
    # }else 
    if(!is.null(myrasts())){
      # make the poparray object...  [here includes all replicates, years, ages]
      thisarray <- array(0,dim=c(nreps,nyears+1,nages))  # hardcode with 3 stages
      rep=1
      for(rep in 1:nreps){  
        thisrep <- myrasts()[[rep]]
        y=1
        for(y in 1:(nyears+1)){
          st <- (y-1)*nages+1
          end <- y*nages
          thisyear <- thisrep[[st:end]]  # extract this year
          abunds <- terra::global(thisyear,"sum",na.rm=T)[,1]
          thisarray[rep,y,] <- abunds
        }
      }
      return(thisarray)
    }else{
      NULL
    }
    
  })
  
  # poparray_clipped: this is a recipe for creating a clipped version of poparray
  poparray_clipped <- reactive({
    if(!is.null(myReactives$polygon)){
      if(!is.null(myrasts())){
        # make the poparray object...  [here includes all replicates, years, ages]
        thisarray <- array(0,dim=c(nreps,nyears+1,nages))  # hardcode with 3 stages
        rep=1
        for(rep in 1:nreps){   
          thisrep <- myrasts_clipped()[[rep]]
          y=1
          for(y in 1:(nyears+1)){
            st <- (y-1)*nages+1
            end <- y*nages
            thisyear <- thisrep[[st:end]]  # extract this year
            abunds <- terra::global(thisyear,"sum",na.rm=T)[,1]
            thisarray[rep,y,] <- abunds
          }
        }
        return(thisarray)
      }else{
        NULL
      }
    }else{
      NULL
    }
  
  })
  
  observeEvent(input$go_download,{
    # download data bundle 
    #make a tempfile with the name we want, saving in the tempdir()
    withProgress(message = "Downloading in progress...",{
      # if the files don't already exist in the temp dir, then download them!
      if(!all(file.exists(my_dl_list()))){
        #download the files to those locations
        #these are the paths to where the files should be saved
        tar1<-file.path(tempdir(),allscenarios$nameA[[scenndx()]])
        tar2<-file.path(tempdir(),allscenarios$nameB[[scenndx()]])
        
        #download files A and B (each has 10 reps in it; or can just do the first one and only have 10 reps and a single download per scenario)
        download.file(allscenarios$urlA[[scenndx()]],file.path(tempdir(),allscenarios$nameA[[scenndx()]]),mode="wb")
        incProgress(1/2)
        download.file(allscenarios$urlB[[scenndx()]],file.path(tempdir(),allscenarios$nameB[[scenndx()]]),mode="wb")
        incProgress(1/2)
        #untar both files into scenpath (which is "sc_thisscen#" inside of tempdir)
        untar(tar1,exdir=tempdir())
        untar(tar2,exdir=tempdir())
        # list.files(tempdir()) #it works!  
        
        # for(i in 1:length(allscenarios$url[[scenndx()]])){ #input the scenario selected by the user 
        #   download.file(allscenarios$url[[scenndx()]][i],my_dl_list()[i],mode = "wb")
        #   incProgress(1/nreps)
        #   Sys.sleep(1) #don't bombard the server... 
        #   # increment progress bar here??
        # }
      }
      # save the scenario to the myreactives object to indicate it has been run  [not used rn]
      myReactives$scensRun <- unique(c(myReactives$scensRun,scenndx()))   # could be used in GUI to allow user to refer back to already-run scenarios...
      showNotification("Download complete!")
      output$main_message <- renderText("Your selected scenario has been successfully downloaded!")
    })


    #plot(myrasts[[1]]) #test to make sure they look ok
    
    
    # myReactives$scensRun_df <- rbind(myReactives$scensRun_df,allscenarios[scenndx(),])
    
    # temp <- poparray()  # making poparray takes quite a bit of time the first time.. 
    # 
    
  })   # end download files from GDrive
  
  output$history_table <- renderTable(allscenarios_short[as.numeric(myReactives$scensRun),] )
  output$history_list <- renderText(as.numeric(myReactives$scensRun))
  
  
  thispoly <- reactive({
    req(myReactives$polygon)
    thisvect <- terra::vect(myReactives$polygon)
    terra::project(thisvect,r)
  })
  
  
  comp_pick <- reactive({
    sort(unique(myReactives$compdf$scenario))
  })
  
  comp_pick_clipped <- reactive({
    sort(unique(myReactives$compdf_clipped$scenario))
  })
  
  # observe({
  #   req(poparray())
  #   thispa <- poparray()
  #   saveRDS(thispa,file=sprintf("%s/poparray_%s.rds",tempdir(),scenndx()))
  # })
  
  # output$main_message <- renderText("") %>% 
  #   bindEvent(input$analysis)
  # 
  # output$main_message <- renderText("") %>% 
  #   bindEvent(input$go_download)
  # 
  observeEvent(input$go_analysis,{  
    
    # temp = dod_polygon()   # FOR DEBUG
    # temp=user_polygon()   # FOR DEBUG!! 
    output$main_message <- renderText("")  # grey out the main_message text
    
    showNotification("Processing scenario data (this can take over a minute)...",duration=150,id="startgo")
    
    if(!all(file.exists(my_dl_list()))){
      showModal(modalDialog( #shows popup window for below message
        title = "Careful!",
        "Hey, you need to click the 'download' button to load the scenario files before running analyses!",
        easyClose = TRUE,
        footer = NULL
      ))
    }
    req(all(file.exists(my_dl_list())))
    
    # use validate messages?
    
    if((myReactives$numClicks>0) & (myReactives$numClicks<4)){
      showModal(modalDialog( #shows popup window for below message
        title = "Careful!",
        "you need at least 4 vertices to make a valid ROI polygon. \n
        If you want to run the analysis for a specific ROI, \n
        please keep clicking on the map \n
        until you have created a valid polygon object (at least 4 vertices). \n 
        If you intend to run a rangewide analysis \n
        clear the ROI and then proceed. Clear the ROI if you don't want to \n 
        see this message.
      ",
        easyClose = TRUE,
        footer = NULL
      ))
      removeNotification(id="startgo")
    }
    
    # check to make sure polygon overlaps with DT range... 
    
    if(!is.null(myReactives$polygon)){
      
      check <- any(!is.na(terra::extract(myrasts()[[1]][[1]],thispoly())[,2]))
      if(!check){
        showModal(modalDialog( #shows popup window for below message
          title = "Careful!",
          "Your polygon does not intersect the Desert Tortoise range. Please \n 
          draw or load another polygon and try again!",
          easyClose = TRUE,
          footer = NULL
        ))
      }
      req(check)
    }
    
    # run poparray and update the comparison df
    
    temp <- poparray()   # run poparray
    
    ## make reactive object that stores plotting info for all scenarios that have been run
    
    myReactives$compdf <- rbind(myReactives$compdf,make_compdf(poparray(),scenndx()))
    
    updateCheckboxGroupInput(session, "compare_scen",
                             choices = comp_pick(), selected = comp_pick()) 
  
    if(is.null(myReactives$polygon)){    # myReactives$numClicks==0
     
      removeNotification(id="startgo")
      showNotification("Processing complete!")
      output$main_message <- renderText("Processing complete!")
    }else{
      # output$main_message <- renderText("Processing started...")
      temp <- poparray_clipped()    # need to run both clipped and unclipped poparrays- can take a while
      # temp2 <- poparray()
      
      myReactives$compdf_clipped <- rbind(myReactives$compdf_clipped,make_compdf(poparray_clipped(),scenndx()))
      
      updateCheckboxGroupInput(session, "compare_scen_clipped",
                               choices = comp_pick_clipped(), selected = comp_pick_clipped())

      removeNotification(id="startgo")
      showNotification("Processing complete!")
      output$main_message <- renderText("Processing complete!")
    }
    
    showTab(inputId="mynavbar", target="Output Viewer") #"Compare"   # show tab when analyses are complete

  })  # end 'go button
  
  
  # render the trajectory plots!
  
  output$trajec_plot_full <- renderPlot({
      req(poparray())
      maxabund <- max(sapply(1:nreps,function(t) max(apply(poparray()[t,,],1,sum)/1000) ))*1.2
      
      plot(1,1,pch="",xlim=c(1,nyears+1),
           ylim=c(0,maxabund),
           xaxt="n",xlab="Years",ylab="Abundance (thousands)", main="")  #paste0("Scenario ",thisscen," ",thisclim())
      
      temp <- lapply(1:nreps,function(t) lines(1:(nyears+1),apply(poparray()[t,,],1,sum)/1000,col=gray(0.6)) ) 
      
      lines(1:(nyears+1),apply(sapply(1:nreps,function(t) apply(poparray()[t,,],1,sum)),1,median)/1000,lwd=2)
      
      axis(1,at=seq(1,nyears+1,length=5),labels = round(seq(2010,2010+nyears,length=5))) 
      
  })
  
  output$trajec_plot_main <- renderPlot({
    
    if(is.null(myReactives$polygon)){
      req(poparray())
      maxabund <- max(sapply(1:nreps,function(t) max(apply(poparray()[t,,],1,sum)/1000) ))*1.2
      
      plot(1,1,pch="",xlim=c(1,nyears+1),
           ylim=c(0,maxabund),
           xaxt="n",xlab="Years",ylab="Abundance (thousands)", main="")  #paste0("Scenario ",thisscen," ",thisclim())
      
      temp <- lapply(1:nreps,function(t) lines(1:(nyears+1),apply(poparray()[t,,],1,sum)/1000,col=gray(0.6)) ) 
      
      lines(1:(nyears+1),apply(sapply(1:nreps,function(t) apply(poparray()[t,,],1,sum)),1,median)/1000,lwd=2)
      
      axis(1,at=seq(1,nyears+1,length=5),labels = round(seq(2010,2010+nyears,length=5))) 
      
    }else{   # if full ROI polygon is available
      # display population trajectory over time 
      req(poparray_clipped())
      maxabund <- max(sapply(1:nreps,function(t) max(apply(poparray_clipped()[t,,],1,sum)/1000) ))*1.2
      
      plot(1,1,pch="",xlim=c(1,nyears+1),
           ylim=c(0,maxabund),
           xaxt="n",xlab="Years",ylab="Abundance (thousands)", main="")  #paste0("Scenario ",thisscen," ",thisclim())
      
      temp <- lapply(1:nreps,function(t) lines(1:(nyears+1),apply(poparray_clipped()[t,,],1,sum)/1000,col=gray(0.6)) ) 
      
      lines(1:(nyears+1),apply(sapply(1:nreps,function(t) apply(poparray_clipped()[t,,],1,sum)),1,median)/1000,lwd=2)
      
      axis(1,at=seq(1,nyears+1,length=5),labels = round(seq(2010,2010+nyears,length=5)))
      
    }
  })
  
  makeagedf <- function(parr){
    df <- NULL
    # y=2
    for(y in 1:(nyears+1)){
      dat <- parr[,y,]/1000
      this = data.frame(
        year= rep(2010+y,nages),
        age=c("Hatchling","Juvenile","Adult"),
        abund = apply(dat,2,median),
        lb = apply(dat,2,function(t) quantile(t,0.05)),
        ub = apply(dat,2,function(t) quantile(t,0.95))
      )
      df <- rbind(df,this)
    }
    return(df)
  }
  
  agedf <- reactive({
    req(poparray())
    df <- makeagedf(poparray())
    return(df)
  })
  
  agedf_clipped <- reactive({
    req(poparray_clipped())
    df <- makeagedf(poparray_clipped())
    return(df)
  })
  
  output$trajec_plot_full_age <- renderPlot({
    req(poparray())
    # maxabund <- max(sapply(1:nreps,function(t) max(poparray()[t,,]/1000) ))*1.2
    
    df <- agedf()
    
    ageplot <- ggplot(df,aes(year,abund,col=age)) +
      geom_ribbon(aes(ymin=lb,ymax=ub,fill=age),alpha=0.4) +
      geom_path(lwd=1.5) +
      labs(y="Total abundance (thousands)",x="Year",title = "Full Range")
    print(ageplot)
    
  })
  
  output$trajec_plot_main_age <- renderPlot({
    
    if(is.null(myReactives$polygon)){
      req(poparray())
      df <- agedf()
      title <- "Full Range"
    }else{
      req(poparray_clipped())
      df <- agedf_clipped()
      title = "Selected ROI"
    }
    ageplot <- ggplot(df,aes(year,abund,col=age)) +
      geom_ribbon(aes(ymin=lb,ymax=ub,fill=age),alpha=0.4) +
      geom_path(lwd=1.5) +
      labs(y="Total abundance (thousands)",x="Year",title = title)
    print(ageplot)
  })
  
  output$trajec_plot_compare_clipped <- renderPlot({
    
    if(is.null(myReactives$polygon)){
      # display reg plot
      these_scens <- as.character(input$compare_scen)
      thisdf <- myReactives$compdf
      thisdf[,2:4] <- thisdf[,2:4]/1000
      thisdf$year2 <- 2010+thisdf$year
      thisdf$scenario <- as.character(thisdf$scenario)
      thisdf <- thisdf %>% filter(scenario%in%these_scens)
      title="Full range"
    }else{
      these_scens <- as.character(input$compare_scen_clipped)
      thisdf <- myReactives$compdf_clipped
      thisdf[,2:4] <- thisdf[,2:4]/1000
      thisdf$year2 <- 2010+thisdf$year
      thisdf$scenario <- as.character(thisdf$scenario)
      thisdf <- thisdf %>% filter(scenario%in%these_scens)
      title="Selected ROI"
    }
    
    if(nrow(thisdf)>0){
      compplot <- ggplot(thisdf,aes(year2,median,col=scenario)) +
        geom_ribbon(aes(ymin=lb,ymax=ub,fill=scenario),alpha=0.4) +
        geom_path(lwd=2) +
        labs(y="Total abundance (thousands)",x="Year",title = title)
      print(compplot)
    }else{
      NULL
    }
    
  })
  
  output$trajec_plot_compare <- renderPlot({
    these_scens <- as.character(input$compare_scen)
    thisdf <- myReactives$compdf
    thisdf[,2:4] <- thisdf[,2:4]/1000
    thisdf$year2 <- 2010+thisdf$year
    thisdf$scenario <- as.character(thisdf$scenario)
    thisdf <- thisdf %>% filter(scenario%in%these_scens)
    if(nrow(thisdf)>0){
      compplot <- ggplot(thisdf,aes(year2,median,col=scenario)) +
        geom_ribbon(aes(ymin=lb,ymax=ub,fill=scenario),alpha=0.4) +
        geom_path(lwd=2) +
        labs(y="Total abundance (thousands)",x="Year",title = "Full range")
      print(compplot)
    }else{
      NULL
    }
    
      
  })
  
  #NEW#download ROI polygon 
  #output$downloadPolygon <- downloadHandler({
  #  filename = function() {
  #    paste('data-', Sys.Date(), '.csv', sep='')
  #  }
  #  content = function(file) {
  #    #csv_write<-array(0,dim=c(length(GHI_D),15))
  #    csv_write<-cbind("hello")
  #    write.csv(csv_write,row.names=FALSE, na="")
  #  }
  #})
  
  output$downloadPolygon <- downloadHandler(
      # check to make sure polygon is available??
    filename = function() {
      sprintf("scen_%s_polygonid_%s.zip", scenndx(), round(runif(1,1000,9999)) )
    },
    content = function(file) {
      data = myReactives$polygon # get polygon object
      temp_shp1 <- sprintf("%s/shp",tempdir())  # use tempdir to create space for shapefile download
      if(!dir.exists(temp_shp1)) dir.create(temp_shp1)
      temp_shp2 <- sprintf("%s/dl.shp",temp_shp1)
      st_write(data,temp_shp2,append=F)   # write the shapefile, overwrite if needed
      # zip all the shp files
      zip_file <- file.path(temp_shp1, "myROI.zip")  # 
      # zip_file2 <- "myROI.zip"
      shp_files <- list.files(temp_shp1,
                              "dl",
                              full.names = TRUE)
      zip(zip_file, files = shp_files, extras="-j")  # zip the files  [NOTE: this adds the whole filepath to the zip... ] [not sure this works on other platforms...]
      # copy the zip file to the file argument
      file.copy(zip_file, file)
      # remove all the files created
      file.remove(zip_file, shp_files)
    }
  )   # end download handler
  
  output$map_plot_main <- renderPlot({

    thisyear <- as.numeric(input$map_year)
    thisquant <- as.numeric(input$map_percentile)/100
    ndx <- (thisyear*nages+1):((thisyear+1)*nages)
    
    if(is.null(myReactives$polygon)){
      req(poparray())
      thismap <- do.call("c",lapply(1:nreps, function(t) sum(myrasts()[[t]][[ndx]]) ) )
    }else{
      req(poparray_clipped())
      thismap <- do.call("c",lapply(1:nreps, function(t) sum(myrasts_clipped()[[t]][[ndx]]) ) )
    }
    
    globalmax <- max(global(thismap,"max",na.rm=T))
    qmap <- quantile(thismap,thisquant)
    qmap <- trim(qmap)
    plot(qmap,range=c(0,globalmax) )
  })
  
  output$map_plot_full <- renderPlot({
  
    thisyear <- as.numeric(input$map_year)
    thisquant <- as.numeric(input$map_percentile)/100
    ndx <- (thisyear*nages+1):((thisyear+1)*nages)
    
    req(poparray())
    
    
    thismap <- do.call("c",lapply(1:nreps, function(t) sum(myrasts()[[t]][[ndx]]) ) )
    globalmax <- max(global(thismap,"max",na.rm=T))
    qmap <- quantile(thismap,thisquant)
    qmap <- trim(qmap)
    plot(qmap,range=c(0,globalmax) )
  })
 
}  # end server

# RUN APP -----------------

# Run the app
shinyApp(ui, server)


# 
# ## scratch
# df <- read.csv("thisROI.csv")
# polygon <- df %>%
#   st_as_sf(coords = c("lng", "lat"), crs = 4326) %>%
#   summarise(geometry = st_combine(geometry)) %>%
#   st_cast("POLYGON")
# polygon
# plot(polygon)
# leaflet(polygon) %>% addTiles() %>% addPolygons()
# 
# line <- df[, c("lng","lat")] %>%
#   as.matrix() %>%
#   st_linestring() %>%
#   st_sfc() %>%
#   st_sf()
# line
# plot(line)
# 
# leaflet(line) %>% addTiles() %>% addPolylines()
