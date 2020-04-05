  #'
  #' #' Outbreaks Predictor and Policy Analyzer ("FeverFighter")
  #' (c) 2016-2020 By Authors 
  #' You are free to use this software under CC-BY 4.0 https://creativecommons.org/licenses/by/4.0/
  #'
  #' Shiny app.  Some code from rstudio.com
  #'
  library(shiny)
  library(lubridate)
  library(markdown)
  
  print('sourcing the model code...')
  source('schools_solver.R', local = T)
  print('source completed successfully!')
  
  # Define UI for the app ----
  ui <- fluidPage(
    titlePanel("FeverFighter - the outbreak predictor and policy analyzer"),
    tags$head(includeScript('google_analytics.js')),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        actionButton("runButton", "Run model ...", style = "fill",
                     color = "blue", autofocus="autofocus"), 
        br(),# element to introduce extra vertical spacing ----
        h4('Enter outbreak setting'),
        #students per grade
        selectInput("symptoms_data", label=NULL, c("COVID-19" = "coronavirus_symptoms.csv", "Influenza" = "influenza_symptoms.csv"), multiple = F, selected = 'coronavirus_symptoms.csv'),
        dateRangeInput('outbreak_dates', label='Outbreak dates', start = params_def_elementary$day_start, 
                       end = params_def_elementary$day_end, 
                       min = NULL, max = NULL, format = "yyyy-mm-dd", startview = "month", weekstart = 0,
                       language = "en", separator = " to ", width = NULL),
        textInput('cohort_sizes', label='Persons per cohort (persons, students, etc)', placeholder = '70,70,70,70,70,70', value = cat(params_def_elementary$cohort_sizes, sep = ',')),
        textInput('initial_school', label='Initially infected not isolated',       placeholder = '1,0,0,0,0,0', value = cat(params_def_elementary$initially_infected_school, sep = ',')),
        textInput('initial_home',   label='Initially infected isolated at home',   placeholder = '0,0,0,0,0,0', value = cat(params_def_elementary$initially_infected_home, sep = ',')),
        selectInput("focal_symptom", "Symptom being monitored", c("Systemic (fever)" = "fever", "Respiratory (coughing)" = "cough", "Nasal (sneezing)" = "sneezing"), multiple = F, selected = 'fever'),
        numericInput('symptom_propensity', label = 'Fraction with symptoms', min=0, max=1.0, step=0.05, value = params_def_elementary$symptom_propensity),
        h4('Policy for isolation'),
        numericInput('exclusion_days', label  = 'Home isolation after last symptom (days)', min=0, max=9, step=1, value = params_def_elementary$exclusion_days),
        numericInput('workweek_days', label   = 'Workweek (days)', min=0, max=7, step=1, value = params_def_elementary$workweek_days),
        numericInput('compliance', label      = 'Policy compliance (0=none, 1=total)', min=0, max=1.0, step=0.01, value = params_def_elementary$compliance),
        numericInput('symptom_attention', label = 'Symptom attention (0=none, 1=total)', min=0, max=1.0, step=0.01, value = params_def_elementary$symptom_attention),
        dateRangeInput('wb_object', label='Winter Holiday Range', start = params_def_elementary$wb_start, end = params_def_elementary$wb_end, min = NULL, 
                       max = NULL, format = "yyyy-mm-dd", startview = "month", weekstart = 0,
                       language = "en", separator = " to ", width = NULL),
        dateRangeInput('sb_object', label='Spring Holiday Range', start = params_def_elementary$sb_start, end = params_def_elementary$sb_end, min = NULL, 
                        max = NULL, format = "yyyy-mm-dd", startview = "month", weekstart = 0,
                        language = "en", separator = " to ", width = NULL),
        dateRangeInput('sv_object', label='Summer Holiday Range', start = params_def_elementary$sv_start, end = params_def_elementary$sv_end, min = NULL, 
                       max = NULL, format = "yyyy-mm-dd", startview = "month", weekstart = 0,
                       language = "en", separator = " to ", width = NULL),
        
        h4('Set or calibrate parameters'),
        numericInput('vac_rate', label = 'Vaccination rate (0=none, 1=total)', min=0, max=1.0, step=0.01, value = params_def_elementary$vac_rate),
        numericInput('vac_efficacy', label = 'Vaccination efficacy (0=none, 1=sterilizing)', min=0, max=1.0, step=0.01, value = params_def_elementary$vac_efficacy),
        numericInput('transmissibility', label = 'Transmissibility', min=0, max=10.0, step=0.001, value = params_def_elementary$transmissibility),
        numericInput('transmissibility_weekend_ratio', label = 'Relative transmissibility on weekends (0=no transmission, 1=unmodified)', min=0, max=2.0, step=0.05, value = params_def_elementary$transmissibility_weekend_ratio),
        numericInput('transmissibility_closure_ratio', label = 'Relative transmissibility holiday and closure (0=no transmission, 1=unmodified)', min=0, max=2.0, step=0.05, value = params_def_elementary$transmissibility_closure_ratio),
        numericInput('seasonality', label = 'Transmission seasonality (0=none, 1=complete)', min=0, max=1.0, step=0.05, value = params_def_elementary$seasonality),
        numericInput('contacts_cross_grade', label = 'Cross-cohort contact rate (0=none, 1=unlimited)', min=0, max=1.0, step=0.01, value = -1),
        fileInput("contact_data_new", label = NULL, buttonLabel = "Contact matrix CSV File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
        fileInput("symptoms_data_new", label = NULL, buttonLabel = "Symptoms by day CSV File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
        fileInput("loaded_parameters", label = NULL, buttonLabel = "Load parameters",multiple = FALSE,accept=c("text/rds",".rds")),
        downloadButton('save_inputs', 'Download parameters (sug. use browser)')
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        tabsetPanel(type = "tabs", id = "maintabs",
                    tabPanel("Introduction", includeMarkdown('introduction.md')),
                    tabPanel("Summary",  
                             #h2('Outbreak summary'),
                             plotOutput("finalResultsPlot"),
                             textOutput("finalResultsDetailed")
                             #h2('All metrics'),
                             #tableOutput("finalResultsTable")
                             #verbatimTextOutput("finalResultsPlot2")
                    ),
                    tabPanel("Curves", h3("Epidemic curves"), plotOutput("timePlot"), 
                             h3("Calibration to data"), 
                             selectInput("calibration_dataset_name", label=NULL, c("Smith pH1H1 2009", "Marchbanks pH1N1 2009", "McCann Influenza B 2012", "Custom calibration"), multiple = F, selected = "Custom calibration"),
                             fileInput("calibration_data", label = NULL, buttonLabel = "Calibration CSV File", multiple = FALSE, accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                             h5(textOutput("calibScore")), plotOutput("calibPlot")), 
                    tabPanel("Daily Data", tableOutput("outbreak_summary_table")),
                    tabPanel("Help", includeMarkdown('help.md'))
                    
        )
        
      )
    )
  )
  
  run_outbreak <- function(params, input, output) {
    new_params <- params
    new_params$seasonality      = input$seasonality
    new_params$transmissibility = input$transmissibility
    new_params$transmissibility_weekend_ratio = input$transmissibility_weekend_ratio
    new_params$transmissibility_closure_ratio = input$transmissibility_closure_ratio
    new_params$focal_symptom      = input$focal_symptom
    
    if(input$focal_symptom != "fever" && input$symptoms_data == "coronavirus_symptoms.csv") {
      showNotification("Focal symptom must be 'fever' for COVID-19", type="error")
      return(NULL)
    }
    new_params$symptom_propensity = input$symptom_propensity
    
    new_params$vac_rate         = input$vac_rate
    new_params$vac_efficacy     = input$vac_efficacy
    new_params$exclusion_days   = input$exclusion_days
    new_params$workweek_days    = input$workweek_days
    new_params$compliance       = input$compliance
    new_params$symptom_attention  = input$symptom_attention
  
    #wishlist: the start and end day appear off by one day in the plot
    new_params$day_start     = ymd(input$outbreak_dates[1], quiet = FALSE, tz = NULL)
    new_params$day_end       = ymd(input$outbreak_dates[2], quiet = FALSE, tz = NULL)
    new_params$wb_start = ymd(input$wb_object[1], quiet = FALSE, tz = NULL)
    new_params$wb_end   = ymd(input$wb_object[2], quiet = FALSE, tz = NULL)
    new_params$sb_start = ymd(input$sb_object[1], quiet = FALSE, tz = NULL)
    new_params$sb_end   = ymd(input$sb_object[2], quiet = FALSE, tz = NULL)
    new_params$sv_start = ymd(input$sv_object[1], quiet = FALSE, tz = NULL)
    new_params$sv_end   = ymd(input$sv_object[2], quiet = FALSE, tz = NULL)
    
    if (new_params$day_start >= new_params$day_end) {
      showNotification("Please ensure that outbreak end date > start date", type="error")
      return(NULL)
    }
    if ((new_params$wb_start >= new_params$wb_end) || 
        (new_params$sb_start >= new_params$sb_end) ||
        (new_params$sv_start >= new_params$sv_end)){
      showNotification("Please ensure that for all holidays, end date > start date", type="error")
      return(NULL)
    }
    contact_data_new <- input$contact_data_new
    if(is.null(contact_data_new)) {
      contact_data_new <- default_elementary_contact_data
    } else {
      contact_data_new <- read.csv(input$contact_data_new$datapath, header=T, row.names=1) 
      if(! all(diag(as.matrix(contact_data_new)) == 1)) {
        showNotification("Contact rate matrix should have 1.0 on the diagnoal", type="error")
        #wishlist: in the future we could remove this restriction
      }
    }
    if (input$contacts_cross_grade != -1) {
      contact_data_new <- update_off_diagonal_values(default_elementary_contact_data, input$contacts_cross_grade)
    }
    if(!is.null(input$symptoms_data_new)) {
      symptoms_data_new <- read.csv(input$symptoms_data_new$datapath) 
    } else {
      symptoms_data_new <- read.csv(input$symptoms_data)
    }
    #'' is the value in these fields by default, if the user doesn't enter anything
    if(input$cohort_sizes != '') {
      students_per_grade <- as.numeric(unlist(strsplit(input$cohort_sizes,",")))
    } else {
      students_per_grade <- new_params$cohort_sizes
    }
    if(input$initial_school != '') {
      initial_school <- as.numeric(unlist(strsplit(input$initial_school,",")))
    } else {
      initial_school <- new_params$initial_infected_school
    }
    if(input$initial_home != '') {
      initial_home <- as.numeric(unlist(strsplit(input$initial_home,",")))
    } else {
      initial_home <- new_params$initial_infected_home
    }
    
    if(length(students_per_grade) == 1 && 
       length(initial_school) == 1 && 
       length(initial_home) == 1) {
      students_per_grade <- c(students_per_grade, 0)
      initial_school     <- c(initial_school, 0)
      initial_home       <- c(initial_home, 0)
      #wishlist: this is a workaround to a bug with handling a single cohort. fix properly.
    }
      
    new_params <- tryCatch(
      params_update_school(new_params, 
                                         students_per_cohort=students_per_grade, 
                                         initially_infected_home=initial_home, 
                                         initially_infected_school=initial_school, 
                                         contact_data=contact_data_new,
                                         symptoms_data=symptoms_data_new)
      ,error = function(e) {
        traceback(e)
        showNotification("Error with parameters values.  Update values or restart the app", type="error")
        return(NULL)
      }
    )
    if(is.null(new_params)) {
      return(NULL)
    }

    #TODO: ensure that there are no trailing empty columns in initial_infected_home
    showNotification("Running ..", type="message")
    
    outbreak_result <- compute_school_outbreak(params = new_params)
    
    return(outbreak_result)
  }
  
  run_outbreak_show_results <- function(input, output, sessionData, session) {
    if(! (sessionData$ui_ready)) {
      return(NULL)
    }
    res <- tryCatch(run_outbreak(params = params_def_elementary, input=input, output=output)
                    ,error = function(e) {
                      traceback(e)
                      showNotification("Error running the epidemic. Please review parameter values or restart the app", type="error")
                      return(NULL)
                    })
    
    if(!is.null(res)) {
      sessionData$outbreakResult <- res
      if(input$maintabs == 'Introduction') {
        updateNavbarPage(session, inputId='maintabs', selected = 'Summary')
      }
    }
  }
  
  # Define server logic for random distribution app ----
  server <- function(input, output, session) {
    sessionData     <- reactiveValues()
    sessionData$ui_ready <- FALSE
     
    output$finalResultsPlot <- renderPlot({
      if(is.null(sessionData$outbreakResult)) {return(NULL)}
      final_results <- summarize_epidemic_metrics(sessionData$outbreakResult)
      metrics <- c('total_students', 'vaccinated', 'infected',
                   'persondays_isolated', 'persondays_unisolated')
      metrics <- rev(metrics)
      final_results <- data.frame(metric=metrics,
                                  val=unlist(final_results[metrics]))
      
      p <- ggplot(data=final_results, 
                  aes(x=metric, y=val)) + 
        geom_bar(stat="identity") + coord_flip() + 
        scale_x_discrete(limits=metrics, title('Cases')) +
        scale_y_continuous(title(NULL)) +
        #geom_text(aes(x=metric,y=val,label=val),vjust=0, hjust=0.2) +
        geom_text(aes(label=val), position=position_dodge(width=0.5), 
                  hjust=1.0, color='white') + 
        theme_minimal()
  
        #http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
      return(p)  
    })
    output$finalResultsTable <- renderTable({
      if(is.null(sessionData$outbreakResult)) {return(NULL)}
      outbreak_result <- sessionData$outbreakResult 
      ret <- summarize_epidemic_metrics(outbreak_result)
  
      ret$peak_date_new_ILI <- strftime(ret$peak_date_new_ILI, '%Y-%m-%d')
      return(as.data.frame(ret))  
    })
    output$finalResultsDetailed <- renderText({
      if(is.null(sessionData$outbreakResult)) {return(NULL)}
      outbreak_result <- sessionData$outbreakResult 
      ret <- summarize_epidemic_metrics(outbreak_result)
      
      return(sprintf('Attack rate: %.1f%%. Peak infected persons: %d persons. Peak duration: %d days', 100*ret$attack_rate, ret$peak_infected, ret$peak_duration))
    })
    output$timePlot <- renderPlot({
      if(is.null(sessionData$outbreakResult)) {return(NULL)}
      outbreak_result <- sessionData$outbreakResult
      return(plot_curves_app(outbreak_result))  
    })
    
    output$calibPlot <- renderPlot({
      if(is.null(sessionData$outbreakResult)) {
        return(NULL)
      }
      outbreak_result <- sessionData$outbreakResult
      
      if(!is.null(input$calibration_data) && (input$calibration_dataset_name == "Custom calibration")) {
        calibration_data <- fread(input$calibration_data$datapath) 
      }else{
        calibration_data <- get_calibration_dataset(input$calibration_dataset_name) 
        if(is.null(calibration_data)) {
          return(NULL)
        }
      }
      return(plot_calibration(outbreak_result, calibration_data))
    })
    
    output$calibScore <- renderText({
      outbreak_result <- sessionData$outbreakResult 
      if(is.null(outbreak_result)){
        return(NULL)
      }
      final_metrics   <- summarize_epidemic_metrics(outbreak_result)
      
      final_calibration_metrics <- get_calibration_dataset_results(input$calibration_dataset_name)
      if(is.null(final_calibration_metrics) || (!is.null(input$calibration_data))) {
        return(NULL)
        #wishlist: implement final_calibration_metrics   <- summarize_epidemic_metrics(calibration_data)
      } 
      
      result_difference <- calibration_result_difference(final_calibration_metrics, final_metrics)
      return(sprintf(" Calibration error=%.3f", result_difference))
    })
    
    output$outbreak_summary_table <- renderTable({
      if(is.null(sessionData$outbreakResult)) {return(NULL)}
      outbreak_result <- sessionData$outbreakResult
      outbreak_result2           <- outbreak_result[,c('dates', 'Infected', 'InfectedIsolated', 'InfectedUnisolated', 'New_ILI', 'ILI')]
      outbreak_result2$dates     <- strftime(as.Date(outbreak_result$dates), format="%Y-%m-%d")
      outbreak_result2$Recovered <- outbreak_result$R
      
      return(outbreak_result2)
    })
  
    observeEvent(input$calibration_data, {
      updateSelectInput(session, inputId="calibration_dataset_name", selected="Custom calibration")
    })
    
    #observeEvent(input$calibration_dataset_name, {
    #  if(input$calibration_dataset_name != 'Custom calibration') {
    #    #updateTextInput(session, inputId="calibration_data", value=NULL)
    #  }
    #})
    
    
    observeEvent(input$runButton, {
      sessionData$ui_ready <- TRUE #user must hit this once
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$vac_rate, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$vac_efficacy, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$focal_symptom, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$symptom_propensity, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$exclusion_days, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$workweek_days, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$compliance, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$symptom_attention, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$transmissibility, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    
    
    observeEvent(input$outbreak_dates, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$wb_object, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$sb_object, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$sv_object, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    
    
    observeEvent(input$transmissibility_weekend_ratio, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$seasonality, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    observeEvent(input$contacts_cross_grade, {
      run_outbreak_show_results(input, output, sessionData, session)
    })
    
   
    observeEvent(input$loaded_parameters, {
      sessionData$ui_ready <- FALSE #temporarily disable triggering on input change
      
      inputs_file <- input$loaded_parameters$datapath
      savedInputs <- readRDS(inputs_file)

      if(!is.null(savedInputs$fever_attention)) { #wishlist: remove this is to support old parameter names in old files
        savedInputs$symptom_attention = savedInputs$fever_attention
      }
      if(is.null(savedInputs$workweek_days)) { #wishlist: some calibrations don't have this parameter
        savedInputs$workweek_days = 5
      }
      
      updateTextInput(session, inputId ="cohort_sizes", value=savedInputs$cohort_sizes)
      updateTextInput(session, inputId ="initial_school", value=savedInputs$initial_school)
      updateTextInput(session, inputId ="initial_home", value=savedInputs$initial_home)
      updateTextInput(session, inputId ="vac_rate", value=savedInputs$vac_rate)
      updateTextInput(session, inputId ="vac_efficacy", value=savedInputs$vac_efficacy)
      updateTextInput(session, inputId ="contacts_cross_grade", value=savedInputs$contacts_cross_grade)
      updateDateRangeInput(session, "outbreak_dates", start = savedInputs$outbreak_dates[1], end = savedInputs$outbreak_dates[2])
      updateTextInput(session, inputId ="exclusion_days", value=savedInputs$exclusion_days)
      updateTextInput(session, inputId ="workweek_days", value=savedInputs$workweek_days)
      updateTextInput(session, inputId ="compliance", value=savedInputs$compliance)
      updateTextInput(session, inputId ="symptom_attention", value=savedInputs$symptom_attention)
      updateDateRangeInput(session, "wb_object", start = savedInputs$wb_object[1], end = savedInputs$wb_object[2])
      updateDateRangeInput(session, "sb_object", start = savedInputs$sb_object[1], end = savedInputs$sb_object[2])
      updateDateRangeInput(session, "sv_object", start = savedInputs$sv_object[1], end = savedInputs$sv_object[2])
      updateTextInput(session, inputId ="transmissibility", value=savedInputs$transmissibility)
      updateTextInput(session, inputId ="transmissibility_weekend_ratio", value=savedInputs$transmissibility_weekend_ratio)
      updateTextInput(session, inputId ="transmissibility_closure_ratio", value=savedInputs$transmissibility_closure_ratio)
      updateTextInput(session, inputId ="seasonality", value=savedInputs$seasonality)
      
      updateSelectInput(session, inputId="calibration_dataset_name", selected=savedInputs$calibration_dataset_name)
      updateSelectInput(session, inputId="symptoms_data", selected=savedInputs$symptoms_data)
      
      #sessionData$ui_ready <- TRUE, #this cannot be done, sadly, since the above observer events will be triggered after this function ends
    }) 
    
    output$save_inputs <- downloadHandler(
      filename = function() {
        paste("outbreak_parameters_", Sys.Date(), ".rds", sep = "")
      },
      content = function(conn) {
        saveRDS( reactiveValuesToList(input), file = conn)
      },
      contentType = "application/octet-stream"
    )
  }
  
  # Create Shiny app ----
  shinyApp(ui, server)