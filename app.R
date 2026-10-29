library(shiny)
library(shinyjs) # Required for enabling/disabling inputs
library(GLM3r)
library(glmtools)
library(ggplot2)
library(scales)

# Determine app directory cleanly
rm(list=ls())	
glm_dir <- getwd()

# Master baseline file paths (read-only references)
nml_file_base <- file.path(glm_dir, "glm3_1.nml")
met_file_base <- file.path(glm_dir, "bcs", "to_glm.csv")

# Dynamic runtime file paths (generated fresh per simulation run)
nml_file_run     <- file.path(glm_dir, "glm3.nml")
met_file_run     <- file.path(glm_dir, "bcs", "to_glm_run.csv")
inflow_file_run  <- file.path(glm_dir, "bcs", "inflow_run.csv")
outflow_file_run <- file.path(glm_dir, "bcs", "outflow_run.csv")

# Ensure template exists, then initialize glm3.nml from glm3_1.nml
if (file.exists(nml_file_base)) {
  file.copy(nml_file_base, nml_file_run, overwrite = TRUE)
} else if (!file.exists(nml_file_run)) {
  stop("Neither 'glm3_1.nml' nor 'glm3.nml' was found in the working directory.")
}

# Fix potential 'incomplete final line' warnings by re-writing line endings
if (file.exists(nml_file_run)) {
  lines <- readLines(nml_file_run, warn = FALSE)
  writeLines(lines, nml_file_run)
}

# Helper to read baseline values from glm3_1.nml safely
get_initial_val <- function(nml, param, default_val) {
  val <- tryCatch(get_nml_value(nml, param), error = function(e) NULL)
  if (is.null(val) || length(val) == 0) return(default_val)
  return(val)
}

# Fetch initial unscaled base values from your master template (glm3_1.nml)
initial_nml <- read_nml(nml_file_base)
base_kw     <- get_initial_val(initial_nml, "Kw", 0.86)
base_wind   <- get_initial_val(initial_nml, "wind_factor", 1.96)
base_sw     <- get_initial_val(initial_nml, "sw_factor", 1.02)
default_kja <- 8000   # Baseline total KJA plots 17786, assuming only 50% operational 

# Custom color palettes
temp_palette <- c("#00007F", "blue", "#007FFF", "cyan", 
                  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")

# Conversion factors
mmol_to_mgL_O2    <- 0.0319988
mmol_to_mgL_N     <- 0.0140067
mmol_to_mgL_P     <- 0.0309738
mmolC_to_ugL_Chla <- 0.30

# ==========================================================
# Meteorology Driver Modifier Function (Seasonal Wind Forcing)
# ==========================================================
apply_climate_forcing <- function(base_file, run_file, temp_offset = 0.0, rain_factor = 1.0, wind_factor = 1.0) {
  if (!file.exists(base_file)) {
    stop(paste("Base meteorology file not found at:", base_file))
  }
  
  met <- read.csv(base_file, stringsAsFactors = FALSE)
  
  # Air Temperature Modification
  if ("AirTemp" %in% colnames(met)) {
    met$AirTemp <- met$AirTemp + temp_offset
  }
  
  # Rainfall Modification
  if ("Rain" %in% colnames(met)) {
    met$Rain <- met$Rain * rain_factor
  }
  
  # Seasonal Wind Speed Modification (End of Dec - Feb & Jul - Aug)
  if ("WindSpeed" %in% colnames(met) && "time" %in% colnames(met)) {
    dates <- as.POSIXct(met$time, format = "%Y-%m-%d %H:%M:%S")
    months <- as.numeric(format(dates, "%m"))
    days   <- as.numeric(format(dates, "%d"))
    
    # Identify windy season timestamps:
    # 1. Late Dec (Dec 15 onwards) through February (Months 12 [>=15], 1, 2)
    # 2. July and August (Months 7, 8)
    is_windy_season <- (months == 12 & days >= 15) | (months %in% c(1, 2)) | (months %in% c(7, 8))
    
    # Multiply wind speed ONLY during the designated windy periods
    met$WindSpeed[is_windy_season] <- met$WindSpeed[is_windy_season] * wind_factor
  }
  
  write.csv(met, run_file, row.names = FALSE, quote = FALSE)
}

# ==========================================================
# Inflow & Outflow Generator Helper Function (Non-Destructive)
# ==========================================================
generate_inflow_outflow <- function(kja_count, driver_file, out_inflow, out_outflow) {
  
  if (!file.exists(driver_file)) {
    stop(paste("Driver file not found at:", driver_file))
  }
  
  met <- read.csv(driver_file, sep = ",")
  met$time <- as.POSIXct(strptime(met$time, "%Y-%m-%d %H:%M:%S"))
  
  KJA.num <- kja_count
  
  mg_L_to_mmol_m3 <- function(value_mg_L, molar_mass_g_mol) {
    value_mg_L * 1000 / molar_mass_g_mol
  }
  
  M_P   <- 30.97
  M_PO4 <- 94.97
  M_N   <- 14.01
  
  DON_mg_L <- KJA.num * 0.098 * 4.5E-5 * 114.16
  PON_mg_L <- KJA.num * 0.042 * 4.5E-5 * 114.16
  DOP_mg_L <- KJA.num * 0.01  * 4.8E-5 * 114.16
  POP_mg_L <- KJA.num * 0.02  * 4.8E-5 * 114.16
  PO4_mg_L <- KJA.num * 0.7   * 1e-3    * 114.16
  
  DON <- mg_L_to_mmol_m3(DON_mg_L, M_N) 
  PON <- mg_L_to_mmol_m3(PON_mg_L, M_N) 
  DOP <- mg_L_to_mmol_m3(DOP_mg_L, M_P) 
  POP <- mg_L_to_mmol_m3(POP_mg_L, M_P)
  PO4 <- mg_L_to_mmol_m3(PO4_mg_L, M_PO4)
  
  inflow <- data.frame(
    time    = met$time, 
    temp    = 24.3, 
    flow    = 1, 
    OGM_pop = POP, 
    OGM_dop = DOP, 
    PHS_frp = PO4,
    OGM_pon = PON, 
    OGM_don = DON
  )
  
  outflow <- data.frame(
    time    = inflow$time, 
    flow    = 1
  )
  
  ensure_full_timestamp <- function(df) {
    time_cols <- sapply(df, function(x) inherits(x, "POSIXt"))
    df[time_cols] <- lapply(df[time_cols], function(x) format(x, "%Y-%m-%d %H:%M:%S"))
    return(df)
  }
  
  inflow  <- ensure_full_timestamp(inflow)
  outflow <- ensure_full_timestamp(outflow)
  
  if (!dir.exists("bcs")) dir.create("bcs")
  write.csv(inflow,  out_inflow,  quote = FALSE, row.names = FALSE)
  write.csv(outflow, out_outflow, quote = FALSE, row.names = FALSE)
}

# ==========================================================
# UI
# ==========================================================
ui <- fluidPage(
  useShinyjs(), # Enables shinyjs functions (disable/enable)
  
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 0px; border-bottom: 2px solid #eee; margin-bottom: 20px;",
    uiOutput("app_header_title"),
    div(
      style = "display: flex; align-items: center; gap: 15px;",
      selectInput("lang", NULL, choices = c("English" = "en", "Bahasa Indonesia" = "id"), selected = "en", width = "160px"),
      img(src = "brin_logo.png", style = "height: 55px; width: auto;")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      uiOutput("sidebar_header_climate"),
      
      selectInput(
        inputId  = "climate_preset",
        label    = "IPCC Climate Scenario Presets",
        choices  = c(
          "Custom Settings"                          = "custom",
          "Baseline Present-Day"                     = "baseline",
          "SSP2-4.5 Near-Term (+0.8°C)"               = "ssp245_near",
          "SSP3-7.0 Mid-Century (+1.5°C, -10% Rain)"  = "ssp370_mid",
          "SSP5-8.5 Late-Century (+2.5°C, +15% Rain)" = "ssp585_late"
        ),
        selected = "baseline"
      ),
      
      sliderInput(
        inputId = "temp_offset",
        label   = "Air Temperature Offset (ΔT, °C)",
        min     = -1.0, max = 4.0, value = 0.0, step = 0.1
      ),
      sliderInput(
        inputId = "rain_factor",
        label   = "Rainfall Multiplier Factor",
        min     = 0.5, max = 2.0, value = 1.0, step = 0.05
      ),
      sliderInput(
        inputId = "wind_factor", 
        label   = "Windy Season Multiplier Factor (Dec-Feb & Jul-Aug)",
        min     = 0.5, max = 2.0, value = 1.0, step = 0.05
      ),
      sliderInput(
        inputId = "sw_factor", 
        label   = "Solar Radiation Multiplier Factor",
        min     = 0.5, max = 2.0, value = 1.0, step = 0.05
      ),
      sliderInput(
        inputId = "Kw", 
        label   = "Light Attenuation Coefficient (Kw, m⁻¹)",
        min     = 0.05, max = 2.5, value = base_kw, step = 0.01
      ),
      
      hr(),
      uiOutput("sidebar_header_kja"),
      sliderInput(
        inputId = "kja_count", 
        label   = "Number of KJA Plots (Floating Net Cages)",
        min     = 0, max = 53358, value = default_kja, step = 500
      ),
      
      hr(),
      uiOutput("sidebar_buttons"),
      br(), br(),
      verbatimTextOutput("status")
    ),
    
    mainPanel(
      uiOutput("main_tabs_content")
    )
  ),
  
  hr(style = "margin-top: 30px; margin-bottom: 15px;"),
  div(
    style = "text-align: center; padding-bottom: 20px; color: #666; font-size: 13px;",
    p("© 2026 Tropical Lake Ecosystem Dynamics Research Group - BRIN. All rights reserved.")
  )
)

# ==========================================================
# SERVER
# ==========================================================
server <- function(input, output, session) {
  
  # Reactive Header Title
  output$app_header_title <- renderUI({
    title_txt <- "BESt : Batur Ecosystem Simulation V.1.0"
    h2(title_txt, style = "margin: 0; font-weight: 600;")
  })
  
  # Reactive Sidebar Headers & Action Buttons
  output$sidebar_header_climate <- renderUI({
    h4(if (input$lang == "id") "Skenario Iklim & Fisik" else "Climate & Physical Scenarios")
  })
  
  output$sidebar_header_kja <- renderUI({
    h4(if (input$lang == "id") "Tekanan Antropogenik: Beban Budidaya (KJA)" else "Anthropogenic Pressures: Aquaculture Loading (KJA)")
  })
  
  output$sidebar_buttons <- renderUI({
    is_id <- input$lang == "id"
    tagList(
      actionButton("run_btn", if (is_id) "Jalankan Simulasi" else "Run Model Simulation", icon = icon("play"), class = "btn-primary"),
      actionButton("reset_btn", if (is_id) "Reset ke Baseline" else "Reset to Baseline", icon = icon("undo"), class = "btn-default")
    )
  })
  
  # Dynamically update labels when switching language
  observeEvent(input$lang, {
    is_id <- input$lang == "id"
    
    updateSelectInput(
      session, "climate_preset",
      label = if (is_id) "Preset Skenario Iklim IPCC" else "IPCC Climate Scenario Presets",
      choices = if (is_id) c(
        "Pengaturan Kustom"                            = "custom",
        "Baseline Saat Ini"                            = "baseline",
        "SSP2-4.5 Jangka Pendek (+0.8°C)"              = "ssp245_near",
        "SSP3-7.0 Pertengahan Abad (+1.5°C, -10% Hujan)" = "ssp370_mid",
        "SSP5-8.5 Akhir Abad (+2.5°C, +15% Hujan)"     = "ssp585_late"
      ) else c(
        "Custom Settings"                          = "custom",
        "Baseline Present-Day"                     = "baseline",
        "SSP2-4.5 Near-Term (+0.8°C)"               = "ssp245_near",
        "SSP3-7.0 Mid-Century (+1.5°C, -10% Rain)"  = "ssp370_mid",
        "SSP5-8.5 Late-Century (+2.5°C, +15% Rain)" = "ssp585_late"
      ),
      selected = input$climate_preset
    )
    
    updateSliderInput(session, "temp_offset", label = if (is_id) "Deviasi Suhu Udara (ΔT, °C)" else "Air Temperature Offset (ΔT, °C)")
    updateSliderInput(session, "rain_factor", label = if (is_id) "Faktor Pengali Curah Hujan" else "Rainfall Multiplier Factor")
    updateSliderInput(session, "wind_factor", label = if (is_id) "Faktor Pengali Angin Musim Berangin (Des-Feb & Jul-Agt)" else "Windy Season Multiplier Factor (Dec-Feb & Jul-Aug)")
    updateSliderInput(session, "sw_factor",   label = if (is_id) "Faktor Pengali Radiasi Matahari" else "Solar Radiation Multiplier Factor")
    updateSliderInput(session, "Kw",          label = if (is_id) "Koefisien Attenuasi Cahaya (Kw, m⁻¹)" else "Light Attenuation Coefficient (Kw, m⁻¹)")
    updateSliderInput(session, "kja_count",   label = if (is_id) "Jumlah Petak KJA (Keramba Jaring Apung)" else "Number of KJA Plots (Floating Net Cages)")
  })
  
  # Dynamic Main Panel Tabs
  output$main_tabs_content <- renderUI({
    is_id <- input$lang == "id"
    
    tabsetPanel(
      type = "tabs",
      
      tabPanel(
        if (is_id) "Bagaimana Proses di Danau?" else "How Lake Behaves?",
        br(),
        wellPanel(
          h3(if (is_id) "Proses Fisik & Biogeokimia" else "Physical & Biogeochemical Processes"),
          p(if (is_id) "Diagram ini mengilustrasikan hubungan timbal balik antara hidrodinamika fisik, stratifikasi termal, kualitas air, dan fluks hara di danau." 
            else "This diagram illustrates typical coupled physical hydrodynamics, thermal stratification, water quality, and nutrient flux feedback loops in lakes."),
          hr(),
          div(
            style = "text-align: center; padding: 10px;",
            img(src = "conceptualization.png", style = "max-width: 90%; height: auto; border: 1px solid #ddd; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);")
          ),
          br(),
          h4(if (is_id) "Proses Utama" else "Key Process Highlights"),
          tags$ul(
            tags$li(strong(if (is_id) "Stratifikasi Termal & Hidrodinamika: " else "Thermal Stratification & Hydrodynamics: "), 
                    if (is_id) "Angin dan fluks panas permukaan mengendalikan pencampuran lapisan batas, pembentukan lapisan termal, dan stabilitas termoklin." 
                    else "Wind forcing and surface heat fluxes drive boundary-layer mixing, thermal layer formation, and thermocline stability."),
            tags$li(strong(if (is_id) "Beban Nutrisi: " else "Nutrient Loads: "), 
                    if (is_id) "Masukan dari daerah tangkapan air serta sisa pakan dan kotoran ikan dari Keramba Jaring Apung (KJA) melepaskan bahan organik (karbon, nitrogen, dan fosfor)." 
                    else "Both catchment and unconsumed feed and fish waste from floating net cages (KJA) release organic matter and nutrients, i.e. carbon, nitrogen, and phosphorus."),
            tags$li(strong(if (is_id) "Transformasi Biogeokimia: " else "Biogeochemical Transformations: "), 
                    if (is_id) "Mineralisasi bahan organik mengonsumsi oksigen, memicu anoksia di dasar danau dan daur ulang nutrisi dari sedimen." 
                    else "Organic matter mineralization consumes oxygen, driving hypolimnetic anoxia and benthic nutrient recycling.")
          )
        )
      ),
      
      tabPanel(
        if (is_id) "Tentang Dashboard" else "Dashboard Overview",
        br(),
        wellPanel(
          h3(if (is_id) "Memahami Dinamika Ekosistem Danau" else "Better Understand Lake Ecosystem Dynamics"),
          p(if (is_id) "Dashboard interaktif ini memodelkan limnologi fisik, kualitas air, dan interaksi biogeokimia di " else "This interactive dashboard models the physical limnology, water quality, and biogeochemical interactions in ", 
            strong("Lake Batur (Bali, Indonesia)"), 
            if (is_id) " menggunakan " else " using the ", 
            strong("General Lake Model coupled with Aquatic EcoDynamics (GLM-AED).")),
          hr(),
          
          h4(if (is_id) "Danau Batur & Konteks Lingkungan" else "Lake Batur & Environmental Context"),
          p(if (is_id) "Danau Batur adalah danau kaldera tropis tertutup (endorheic) yang terletak pada ketinggian ~1.050 meter di atas permukaan laut di Bali timur laut. Tanpa outlet sungai permukaan yang signifikan, tinggi muka air dan siklus biogeokimianya dikendalikan oleh curah hujan langsung, limpasan permukaan, pertukaran air tanah, dan evaporasi." 
            else "Lake Batur is a tropical, endorheic (closed-basin) caldera lake situated at ~1,050 meters above sea level in northeastern Bali. Lacking prominent surface river outlets, its hydrologic level and biogeochemical cycles are controlled by direct rainfall, surface runoff, groundwater exchanges, and surface evaporation."),
          
          h4(if (is_id) "Proyeksi Iklim & Sensitivitas Tropis Dataran Tinggi (IPCC AR6)" else "Climate Projections & Highland Tropical Sensitivity (IPCC AR6)"),
          p(if (is_id) "Berdasarkan Laporan Penilaian Ke-6 IPCC (AR6) dan proyeksi iklim regional CMIP6 untuk Indonesia:" 
            else "According to the IPCC 6th Assessment Report (AR6) and regional CMIP6 climate projections for Indonesia:"),
          tags$ul(
            tags$li(strong(if (is_id) "Tren Pemanasan Termal: " else "Thermal Warming Trends: "), 
                    if (is_id) "Suhu rata-rata lingkungan di dataran tinggi Indonesia diproyeksikan meningkat antara +0,5°C hingga +0,8°C pada jangka pendek (2021–2040) dan hingga +1,2°C hingga +2,5°C pada akhir abad di bawah skenario emisi tinggi (SSP3-7.0 / SSP5-8.5)." 
                    else "Mean ambient temperatures in Indonesian highlands are projected to rise between +0.5°C to +0.8°C in the near-term (2021–2040) and up to +1.2°C to +2.5°C by late century under higher-emission pathways (SSP3-7.0 / SSP5-8.5)."),
            tags$li(strong(if (is_id) "Stratifikasi & Pendinginan Malam Hari: " else "Stratification & Nocturnal Cooling: "), 
                    if (is_id) "Di danau tropis dataran tinggi, pemanasan air permukaan melemahkan pencampuran konvektif malam hari selama musim basah, meningkatkan stabilitas stratifikasi termal dan membatasi suplai oksigen ke lapisan hipolimnion. Sebaliknya, udara dingin dan angin kencang saat musim kering meningkatkan pencampuran saat malam hari" 
                    else "In tropical high-elevation lakes, surface water warming weakens nocturnal convective mixing during wet seasons, enhancing thermal stratification stability and restricting oxygen delivery to the hypolimnion. On the other hand, the cold air temperature, together with strong wind, during the dry enhances nocturnal mixing"),
            tags$li(strong(if (is_id) "Curah Hujan Ekstrem: " else "Precipitation Extremes: "), 
                    if (is_id) "Proyeksi menunjukkan intensifikasi kejadian hujan monsun bersamaan dengan periode kekeringan kemarau yang lebih panjang, memperburuk penurunan muka air dan mengubah laju pembilasan nutrisi." 
                    else "Projections indicate intensified monsoonal rainfall events alongside prolonged dry season droughts, exacerbating lake level drawdown and shifting nutrient flushing rates.")
          ),
          
          # -------------------------------------------------------------
          # RESTORED SECTION 1: Anthropogenic Pressures (KJA)
          # -------------------------------------------------------------
          h4(if (is_id) "Tekanan Antropogenik: Keramba Jaring Apung (KJA)" else "Anthropogenic Pressures: Floating Fish Cages (KJA)"),
          p(if (is_id) "Budidaya perikanan intensif menggunakan keramba jaring apung" else "Intensive aquaculture using floating net cages", 
            # em("KJA"), 
            if (is_id) "merupakan pendorong utama beban nutrisi antropogenik di Danau Batur. Sisa pakan yang tidak terkelola dan kotoran metabolik ikan memasok nitrogen organik, fosfor organik, serta ortofosfat dalam jumlah tinggi ke kolom air." 
            else "represents the primary anthropogenic nutrient driver in Lake Batur. Uneaten feed pellets and metabolic fish waste load high quantities of organic nitrogen, organic phosphorus, and orthophosphate into the water column."),
          p(if (is_id) "Pengguna dapat mengubah total " else "Users can modify the total ", 
            strong(if (is_id) "Jumlah petak KJA" else "Number of KJA plots"), 
            if (is_id) " bersama dengan parameter perubahan iklim untuk mensimulasikan dampak lingkungan sinergis terhadap hipoksia air dasar, daur ulang nutrisi internal, dan intensitas blooming fitoplankton." 
            else " alongside climate change levers to simulate synergistic environmental impacts on bottom-water hypoxia, internal nutrient recycling, and phytoplankton bloom intensity."),
          
          # -------------------------------------------------------------
          # RESTORED SECTION 2: Model Architecture (GLM-AED)
          # -------------------------------------------------------------
          h4(if (is_id) "Arsitektur Model: GLM-AED" else "Model Architecture: GLM-AED"),
          tags$ul(
            tags$li(strong(if (is_id) "General Lake Model (GLM): " else "General Lake Model (GLM): "), 
                    if (is_id) "Mensimulasikan lapisan termal dinamik 1D, kedalaman lapisan tercampur, fluks panas, dan tegangan geser angin." 
                    else "Simulates 1D dynamic thermal layering, mixed-layer depth, heat fluxes, and wind-shear stress."),
            tags$li(strong(if (is_id) "Aquatic EcoDynamics (AED): " else "Aquatic EcoDynamics (AED): "), 
                    if (is_id) "Memodelkan proses biogeokimia termasuk penyerapan oksigen terlarut, siklus nitrogen/fosfor, dan kinetika pertumbuhan fitoplankton." 
                    else "Models biogeochemical processes including dissolved oxygen uptake, nitrogen/phosphorus cycling, and phytoplankton growth kinetics.")
          ),
          
          h4(if (is_id) "Cara Menggunakan Alat Ini" else "How to Use This Tool"),
          tags$ol(
            tags$li(if (is_id) "Pilih " else "Select an ", strong(if (is_id) "Preset Skenario Iklim IPCC" else "IPCC Climate Scenario Preset"), if (is_id) " atau pilih Pengaturan Kustom untuk mengatur parameter secara bebas." else " or choose Custom Settings to freely adjust climate parameters."),
            tags$li(if (is_id) "Atur slider " else "Adjust the ", strong(if (is_id) "Jumlah Petak KJA" else "KJA Plot Count"), if (is_id) " untuk menguji berbagai skenario pengelolaan dan kepadatan tebar." else " slider to test different management and stocking density scenarios."),
            tags$li(if (is_id) "Klik " else "Click ", strong(if (is_id) "Jalankan Simulasi" else "Run Model Simulation"), if (is_id) " untuk mengeksekusi GLM-AED." else " to execute GLM-AED."),
            tags$li(strong(if (is_id) "Navigasi" else "Navigate"), if (is_id) " melalui tab hasil dan periksa grafik heatmap kedalaman-waktu dari " else " through the result tabs and inspect the depth-time heatmaps of the ", strong(if (is_id) "Output Model" else "Model Outputs"), ".")
          )
        )
      ),
      
      tabPanel(if (is_id) "Profil Suhu" else "Temperature Profile", br(), plotOutput("tempPlot", height = "600px")),
      tabPanel(if (is_id) "Oksigen Terlarut" else "Dissolved Oxygen", br(), plotOutput("doPlot", height = "600px")),
      tabPanel(if (is_id) "Nitrogen Total" else "Total Nitrogen", br(), plotOutput("tnPlot", height = "600px")),
      tabPanel(if (is_id) "Fosfor Total" else "Total Phosphorus", br(), plotOutput("tpPlot", height = "600px")),
      tabPanel(if (is_id) "Klorofil a" else "Chlorophyll a", br(), plotOutput("chlaPlot", height = "600px")),
      
      tabPanel(
        if (is_id) "Galeri" else "Gallery",
        br(),
        wellPanel(
          h3(if (is_id) "Galeri Sampling Lapangan & Data Observasi" else "Field Sampling & Observational Data Gallery"),
          p(if (is_id) "Dokumentasi visual kegiatan lapangan dan dataset observasi empiris yang dikumpulkan dari Danau Batur." 
            else "Visual documentation of field sampling activities and empirical observation datasets collected from Lake Batur."),
          hr(),
          fluidRow(
            column(
              width = 6,
              div(
                style = "background: #fff; padding: 15px; border: 1px solid #e3e3e3; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); margin-bottom: 20px;",
                h4(if (is_id) "Kegiatan Sampling Lapangan" else "Fieldwork & Sampling Operations", style = "margin-top: 0; text-align: center; font-weight: 600;"),
                hr(style = "margin: 10px 0 15px 0;"),
                div(style = "text-align: center;", img(src = "fieldwork.png", style = "max-width: 100%; height: auto; border-radius: 6px; border: 1px solid #ccc;")),
                p(style = "margin-top: 12px; font-size: 13px; color: #555; text-align: justify;",
                  if (is_id) "Kegiatan pemantauan lapangan di Danau Batur, termasuk pemprofilan CTD, pengambilan sampel kualitas air, dan pemeliharaan sensor suhu dan meteorologi frekuensi tinggi." 
                  else "Field monitoring activities at Lake Batur, including CTD profiling, water quality sampling, and maintenance of high-frequency water column temperature and meteorology sensors.")
              )
            ),
            column(
              width = 6,
              div(
                style = "background: #fff; padding: 15px; border: 1px solid #e3e3e3; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); margin-bottom: 20px;",
                h4(if (is_id) "Plot Data Observasi" else "Observation Data Plot", style = "margin-top: 0; text-align: center; font-weight: 600;"),
                hr(style = "margin: 10px 0 15px 0;"),
                div(style = "text-align: center;", img(src = "observation.png", style = "max-width: 100%; height: auto; border-radius: 6px; border: 1px solid #ccc;")),
                p(style = "margin-top: 12px; font-size: 13px; color: #555; text-align: justify;",
                  if (is_id) "Data terukur in-situ dari suhu kolom air dan pengamatan meteorologi yang digunakan untuk kalibrasi dan validasi model." 
                  else "In-situ measured water column temperature and meteorological observations used for model calibration and validation.")
              )
            )
          )
        )
      ),
      
      tabPanel(
        if (is_id) "Referensi" else "References",
        br(),
        wellPanel(
          h3(if (is_id) "Referensi Utama & Literatur" else "Key References & Literature"),
          p(if (is_id) "Berikut adalah studi ilmiah dasar, dokumentasi manual model, laporan penilaian IPCC, dan literatur proyeksi iklim regional yang mendukung pengembangan dan parameterisasi alat simulasi ini:" 
            else "Below are the foundational studies, model manual documentation, IPCC assessment reports, and regional climate projection literature supporting the development and parameterization of this simulation tool:"),
          hr(),
          tags$ol(
            style = "padding-left: 20px; line-height: 1.8;",
            tags$li(p("Santoso, A. B., Hamilton, D. P., Nomosatryo, S., Sunaryani, A., Rustini, H. A., Rahmadya, A., … Imananda, A. (2026). ", strong("Stratification and mixing dynamics in tropical polymictic Lake Batur, Indonesia."), em(" Inland Waters,"), " 16(1). ", a(href = "https://doi.org/10.1080/20442041.2025.2598457", target = "_blank", "https://doi.org/10.1080/20442041.2025.2598457"))),
            tags$li(p("Sunaryani, A., Soewondo, P., Santoso, A. B., Suharyanto, Wulan, D. R., Nomosatryo, S., & Rahmadya, A. (2026). ", strong("Persistent Eutrophication in a Tropical Endorheic Lake Driven by Sediment–Water Interactions."), em(" Limnological Review,"), " 26(3), 42. ", a(href = "https://doi.org/10.3390/limnolrev26030042", target = "_blank", "https://doi.org/10.3390/limnolrev26030042"))),
            tags$li(p("IPCC (2021). ", strong("Climate Change 2021: The Physical Science Basis. Contribution of Working Group I to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change."), " Chapter 12, pp. 1767–1926. Cambridge University Press. ", a(href = "https://doi.org/10.1017/9781009157896.014", target = "_blank", "https://doi.org/10.1017/9781009157896.014"))),
            tags$li(p("Tangang, F., Chung, J. X., Juneng, L., Supari, S., Salimun, E., Ngai, S. T., Cruz, F., & Aldrian, E. (2020). ", strong("Projected future changes in rainfall in Southeast Asia based on CORDEX–SEA multi-model simulations."), em(" Climate Dynamics,"), " 55, 1247–1267. ", a(href = "https://doi.org/10.1007/s00382-020-05322-2", target = "_blank", "https://doi.org/10.1007/s00382-020-05322-2"))),
            tags$li(p("Supari, S., Tangang, F., Juneng, L., Cruz, F., Chung, J. X., Ngai, S. T., & Aldrian, E. (2020). ", strong("Multi-model projections of precipitation extremes in Southeast Asia based on CORDEX-Southeast Asia simulations."), em(" Environmental Research,"), " 184, 109350. ", a(href = "https://doi.org/10.1016/j.envres.2020.109350", target = "_blank", "https://doi.org/10.1016/j.envres.2020.109350"))),
            tags$li(p("Hipsey, M. R., Bruce, L. C., Boon, C., Busch, B., Carey, C. C., Hamilton, D. P., Hanson, P. C., Read, J. S., de Sousa, E., Weber, M., and Winslow, L. A. (2019). ", strong("A General Lake Model (GLM 3.0) for linking with high-frequency sensor data from the Global Lake Ecological Observatory Network (GLEON)."), em(" Geoscientific Model Development,"), " 12, 473–523. ", a(href = "https://doi.org/10.5194/gmd-12-473-2019", target = "_blank", "https://doi.org/10.5194/gmd-12-473-2019"))),
            tags$li(p("Hipsey, M. R., Boon, C., Paraska, D., Bruce, L. C., & Huang, P. (2019). ", strong("AquaticEcoDynamics/libaed2: v1.3.0-rc2."), " In ", em("Aquatic EcoDynamics (AED) Model Library & Science Manual"), " (Version v1.3.0, pp. 1–34). The University of Western Australia Technical Report. ", a(href = "https://doi.org/10.5281/zenodo.2538495", target = "_blank", "https://doi.org/10.5281/zenodo.2538495")))
          )
        )
      )
    )
  })
  
  # Lock or unlock sliders based on Preset Selection
  climate_sliders <- c("temp_offset", "rain_factor", "wind_factor", "sw_factor", "Kw")
  
  observeEvent(input$climate_preset, {
    preset <- input$climate_preset
    
    if (preset == "custom") {
      lapply(climate_sliders, shinyjs::enable)
    } else {
      lapply(climate_sliders, shinyjs::disable)
      
      if (preset == "baseline") {
        updateSliderInput(session, "temp_offset", value = 0.0)
        updateSliderInput(session, "rain_factor", value = 1.0)
        updateSliderInput(session, "wind_factor", value = 1.0)
        updateSliderInput(session, "sw_factor",   value = 1.0)
        updateSliderInput(session, "Kw",          value = base_kw)
      } else if (preset == "ssp245_near") {
        updateSliderInput(session, "temp_offset", value = 0.8)
        updateSliderInput(session, "rain_factor", value = 1.05)
        updateSliderInput(session, "wind_factor", value = 1.0)
        updateSliderInput(session, "sw_factor",   value = 1.02)
        updateSliderInput(session, "Kw",          value = base_kw)
      } else if (preset == "ssp370_mid") {
        updateSliderInput(session, "temp_offset", value = 1.5)
        updateSliderInput(session, "rain_factor", value = 0.90)
        updateSliderInput(session, "wind_factor", value = 0.95)
        updateSliderInput(session, "sw_factor",   value = 1.05)
        updateSliderInput(session, "Kw",          value = base_kw)
      } else if (preset == "ssp585_late") {
        updateSliderInput(session, "temp_offset", value = 2.5)
        updateSliderInput(session, "rain_factor", value = 1.15)
        updateSliderInput(session, "wind_factor", value = 1.0)
        updateSliderInput(session, "sw_factor",   value = 1.05)
        updateSliderInput(session, "Kw",          value = base_kw)
      }
    }
  }, ignoreInit = FALSE)
  
  # Reset button action
  observeEvent(input$reset_btn, {
    updateSelectInput(session, "climate_preset", selected = "baseline")
    updateSliderInput(session, "temp_offset",  value = 0.0)
    updateSliderInput(session, "rain_factor",  value = 1.0)
    updateSliderInput(session, "Kw",           value = base_kw)
    updateSliderInput(session, "wind_factor",  value = 1.0)
    updateSliderInput(session, "sw_factor",    value = 1.0)
    updateSliderInput(session, "kja_count",    value = default_kja)
  })
  
  # Reactive execution triggered ONLY when "Run Model Simulation" button is clicked
  run_result <- eventReactive(input$run_btn, {
    
    res <- tryCatch({
      withProgress(message = if (input$lang == "id") "Menjalankan Simulasi..." else "Running Simulation...", value = 0.1, {
        
        incProgress(0.15, detail = "Applying climate perturbations & seasonal wind forcing...")
        apply_climate_forcing(
          base_file   = met_file_base,
          run_file    = met_file_run,
          temp_offset = input$temp_offset,
          rain_factor = input$rain_factor,
          wind_factor = input$wind_factor
        )
        
        incProgress(0.2, detail = "Updating aquaculture nutrient loadings...")
        generate_inflow_outflow(
          kja_count   = input$kja_count,
          driver_file = met_file_run,
          out_inflow  = inflow_file_run,
          out_outflow = outflow_file_run
        )
        
        incProgress(0.15, detail = "Updating NML parameters...")
        scaled_sw <- base_sw * input$sw_factor
        
        glm_nml <- read_nml(nml_file_base)
        glm_nml <- set_nml(glm_nml, arg_list = list(
          'Kw'           = input$Kw,
          'wind_factor'  = base_wind,
          'sw_factor'    = scaled_sw,
          'meteo_fl'     = "bcs/to_glm_run.csv",
          'inflow_fl'    = "bcs/inflow_run.csv",
          'outflow_fl'   = "bcs/outflow_run.csv"
        ))
        write_nml(glm_nml, file = nml_file_run)
        
        incProgress(0.3, detail = "Executing GLM-AED2...")
        run_glm(glm_dir)
        
        incProgress(0.1, detail = "Rendering heatmaps...")
      })
      
      list(
        success     = TRUE,
        time        = Sys.time(),
        temp_offset = input$temp_offset,
        rain_factor = input$rain_factor,
        Kw          = input$Kw,
        wind_mult   = input$wind_factor,
        sw_mult     = input$sw_factor,
        sw_actual   = base_sw * input$sw_factor,
        kja         = input$kja_count
      )
      
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })
    
    res
  })
  
  # Status log output
  output$status <- renderText({
    res <- run_result()
    is_id <- input$lang == "id"
    if (!res$success) {
      paste("Error executing GLM:", res$error)
    } else {
      paste0(
        if (is_id) "Simulasi selesai pada: " else "Run completed at: ", format(res$time, "%H:%M:%S"), "\n",
        "───────────────────────────\n",
        if (is_id) "Parameter Diterapkan:\n" else "Applied Parameters:\n",
        " - Temp Offset (ΔT): +", res$temp_offset, " °C\n",
        " - Rain Factor     : ", res$rain_factor, "x\n",
        " - Kw              : ", res$Kw, "\n",
        " - Windy Season Factor: ", res$wind_mult, "x (Dec 15-Feb 28 & Jul-Aug)\n",
        " - SW Multiplier   : ", res$sw_mult, "x (nml sw_factor = ", round(res$sw_actual, 3), ")\n",
        " - Total KJA        : ", format(res$kja, big.mark = ","), " plots"
      )
    }
  })
  
  get_nc_path <- function() {
    nc_path <- file.path(glm_dir, "output", "output.nc")
    req(file.exists(nc_path))
    return(nc_path)
  }
  
  # Plot outputs
  output$tempPlot <- renderPlot({
    res <- run_result(); req(res$success)
    plot_temp(nc_file = get_nc_path(), reference = "surface", legend.title = "Temp (°C)", interval = 0.5, text.size = 12, show.legend = TRUE, legend.position = "right") +
      scale_fill_gradientn(colors = temp_palette, limits = c(22, 32), guide = guide_colorbar(title = "Temp (°C)"))
  })
  
  output$doPlot <- renderPlot({
    res <- run_result(); req(res$success)
    plot_var_nc(nc_file = get_nc_path(), var_name = "OXY_oxy", reference = "surface", xlab = "Date", ylab = "Depth (m)", legend_title = "DO (mg/L)", interval = 0.5, text_size = 12, show_legend = TRUE, legend_position = "right") +
      scale_fill_gradientn(colors = temp_palette, limits = c(0, 296.875), labels = function(x) round(x * mmol_to_mgL_O2, 1), guide = guide_colorbar(title = "DO (mg/L)"))
  })
  
  output$tnPlot <- renderPlot({
    res <- run_result(); req(res$success)
    plot_var_nc(nc_file = get_nc_path(), var_name = "TOT_tn", reference = "surface", xlab = "Date", ylab = "Depth (m)", legend_title = "TN (mg/L)", interval = 0.5, text_size = 12, show_legend = TRUE, legend_position = "right") +
      scale_fill_gradientn(colors = temp_palette, limits = c(0, 356.8879), oob = scales::squish, labels = function(x) round(x * mmol_to_mgL_N, 2), guide = guide_colorbar(title = "TN (mg/L)"))
  })
  
  output$tpPlot <- renderPlot({
    res <- run_result(); req(res$success)
    plot_var_nc(nc_file = get_nc_path(), var_name = "TOT_tp", reference = "surface", xlab = "Date", ylab = "Depth (m)", legend_title = "TP (mg/L)", interval = 0.5, text_size = 12, show_legend = TRUE, legend_position = "right") +
      scale_fill_gradientn(colors = temp_palette, limits = c(0, 48.43397), oob = scales::squish, labels = function(x) round(x * mmol_to_mgL_P, 3), guide = guide_colorbar(title = "TP (mg/L)"))
  })
  
  output$chlaPlot <- renderPlot({
    res <- run_result(); req(res$success)
    plot_var_nc(nc_file = get_nc_path(), var_name = "PHY_TCHLA", reference = "surface", xlab = "Date", ylab = "Depth (m)", legend_title = "Chl-a (µg/L)", interval = 0.5, text.size = 12, show_legend = TRUE, legend_position = "right") +
      scale_fill_gradientn(colors = temp_palette, limits = c(0, 20), guide = guide_colorbar(title = "Chl-a (µg/L)"))
  })
}

shinyApp(ui = ui, server = server)