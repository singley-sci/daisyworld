# ------------------------------------------------------------------------------
# Daisyworld: an interactive model of planetary self-regulation
#
# A teaching Shiny app built on the Watson and Lovelock (1983) Daisyworld model,
# following Dave Bice's STELLA implementation (Penn State) and the SERC/Carleton
# activity. Written in R with tidyverse conventions and the native pipe.
#
# Sources for the equations:
#   Watson, A. J., and J. E. Lovelock (1983), Biological homeostasis of the
#     global environment: the parable of Daisyworld, Tellus B, 35(4), 284-289.
#   Bice, D. (n.d.), Modeling Daisyworld, Penn State,
#     https://personal.ems.psu.edu/~dmb53/DaveSTELLA/Daisyworld/daisyworld_model.htm
#   SERC/Carleton, https://serc.carleton.edu/introgeo/mathstatmodels/examples/DBDaisyworld.html
#
# Run with:  shiny::runApp("app.R")
# ------------------------------------------------------------------------------

# ---- Packages ----------------------------------------------------------------
needed <- c("shiny", "bslib", "ggplot2", "dplyr", "tidyr", "tibble",
            "readr", "scales", "DT")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install these packages first:\n  install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
    call. = FALSE
  )
}

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(scales)

# ---- Visual identity ---------------------------------------------------------
# A cool "observatory" palette. Data colors are colorblind-safe and are always
# paired with a distinct linetype so the plots never rely on color alone.
pal <- list(
  ink    = "#22333B",  # deep slate teal (headings, structure)
  paper  = "#F4F6F5",  # cool off-white surface
  sun    = "#E4A11B",  # amber solar accent (primary actions, the brightening sun)
  teal   = "#3E7C8B",  # secondary accent, planetary albedo line
  black_d = "#2B2B2B", # black daisies
  white_d = "#4F86C6", # white daisies (cool blue stands in for reflective/cool)
  bare   = "#D3A85C",  # bare ground (warm tan)
  living = "#2E6E4E",  # living-planet temperature
  dead   = "#C1442E",  # lifeless-planet temperature
  band   = "#CFE3D6"   # shaded habitable-temperature band
)

app_theme <- bs_theme(
  version = 5,
  bg = pal$paper,
  fg = pal$ink,
  primary = pal$sun,
  secondary = pal$teal,
  base_font = font_collection(font_google("Public Sans", local = FALSE),
                              "system-ui", "sans-serif"),
  heading_font = font_collection(font_google("Fraunces", local = FALSE),
                                 "Georgia", "serif"),
  code_font = font_collection(font_google("IBM Plex Mono", local = FALSE),
                              "monospace")
)

# A shared ggplot theme.
theme_daisy <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(family = "Fraunces", face = "bold",
                                colour = pal$ink, size = rel(1.05)),
      plot.subtitle = element_text(colour = "#5A6A72", size = rel(0.85)),
      plot.caption = element_text(colour = "#5A6A72", size = rel(0.75),
                                  hjust = 0, margin = margin(t = 8)),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E4E9E8"),
      axis.title = element_text(colour = pal$ink),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# ---- Model -------------------------------------------------------------------

#' Daisy growth factor: a downward parabola in local temperature.
#' Peaks at 1 at the optimum temperature and is clamped at zero outside the
#' viable range (a negative growth factor is not physically meaningful in the
#' canonical model, so we hold it at zero, following Watson and Lovelock, 1983).
growth_factor <- function(temp_local, opt_temp, curvature, clamp = TRUE) {
  g <- 1 - curvature * (opt_temp - temp_local)^2
  if (clamp) pmax(g, 0) else g
}

#' Effective death rate at time t. Optionally injects periodic "plagues"
#' (short windows of elevated mortality) as in the Bice experiment set.
effective_death_rate <- function(t, base_rate, plague, plague_times,
                                  plague_rate, plague_width) {
  if (!plague || length(plague_times) == 0) return(base_rate)
  in_plague <- any(abs(t - plague_times) <= plague_width / 2)
  if (in_plague) plague_rate else base_rate
}

#' Simulate Daisyworld with explicit Euler steps (matching the STELLA default),
#' returning one tidy row per time step.
#'
#' Luminosity modes:
#'   "ramp"  luminosity increases linearly from lum_start to lum_end over the run
#'   "fixed" luminosity held constant at lum_fixed (approach to equilibrium)
simulate_daisyworld <- function(t_max = 200, dt = 0.25,
                                init_black = 0, init_white = 0,
                                albedo_black = 0.25, albedo_white = 0.75,
                                albedo_ground = 0.5, death_rate = 0.3,
                                opt_temp = 22.5, curvature = 0.003265,
                                heat_absorp = 20, solar_flux = 917,
                                sb_constant = 5.669e-8,
                                lum_mode = "ramp", lum_start = 0.6,
                                lum_end = 1.8, lum_fixed = 1.0,
                                seed_nudge = 0.001, clamp_growth = TRUE,
                                plague = FALSE, plague_times = c(50, 100, 150),
                                plague_rate = 1.0, plague_width = 4) {

  n_steps <- max(1, round(t_max / dt))
  black <- init_black
  white <- init_white

  rows <- vector("list", n_steps + 1)

  for (i in 0:n_steps) {
    t <- i * dt
    bare <- max(0, 1 - black - white)

    lum <- if (lum_mode == "fixed") {
      lum_fixed
    } else {
      lum_start + t * (lum_end - lum_start) / t_max
    }

    # Planetary albedo is the area-weighted average of the three surfaces.
    albedo <- bare * albedo_ground + black * albedo_black + white * albedo_white

    # Energy balance (Stefan-Boltzmann) gives planetary temperature in degrees C.
    temp_planet <- (lum * solar_flux * (1 - albedo) / sb_constant)^0.25 - 273
    temp_dead   <- (lum * solar_flux * (1 - albedo_ground) / sb_constant)^0.25 - 273

    # Local land temperatures differ from the planetary mean by albedo contrast.
    temp_black_land <- heat_absorp * (albedo - albedo_black) + temp_planet
    temp_white_land <- heat_absorp * (albedo - albedo_white) + temp_planet

    dr <- effective_death_rate(t, death_rate, plague, plague_times,
                               plague_rate, plague_width)

    g_black <- growth_factor(temp_black_land, opt_temp, curvature, clamp_growth)
    g_white <- growth_factor(temp_white_land, opt_temp, curvature, clamp_growth)

    rows[[i + 1]] <- tibble(
      time = t, luminosity = lum,
      area_black = black, area_white = white, area_bare = bare,
      albedo = albedo,
      temp_planet = temp_planet, temp_dead = temp_dead,
      temp_black_land = temp_black_land, temp_white_land = temp_white_land,
      growth_black = g_black, growth_white = g_white,
      death_rate = dr
    )

    # Biflow area change (seed_nudge keeps a trace population so growth can start).
    d_black <- black * (bare * g_black - dr) + seed_nudge
    d_white <- white * (bare * g_white - dr) + seed_nudge

    black <- min(1, max(0, black + dt * d_black))
    white <- min(1, max(0, white + dt * d_white))
  }

  bind_rows(rows)
}

# ---- Built-in scenarios ------------------------------------------------------
# Each scenario is a set of control overrides applied on top of the standard case.
scenario_defaults <- list(
  lum_mode = "ramp", lum_start = 0.6, lum_end = 1.8, lum_fixed = 1.0,
  t_max = 200, dt = 0.25,
  albedo_black = 0.25, albedo_white = 0.75, albedo_ground = 0.5,
  death_rate = 0.3, opt_temp = 22.5, curvature = 0.003265,
  heat_absorp = 20, solar_flux = 917,
  init_black = 0, init_white = 0,
  clamp_growth = TRUE, seed_nudge = 0.001,
  plague = FALSE, plague_rate = 1.0, plague_width = 4,
  plague_times = "50, 100, 150"
)

scenarios <- list(
  "Standard case" = list(),
  "Extreme albedos (0.05 / 0.95)" = list(albedo_black = 0.05, albedo_white = 0.95),
  "Moderate albedos (0.40 / 0.60)" = list(albedo_black = 0.40, albedo_white = 0.60),
  "Cooler-loving daisies (optimum 15 C)" = list(opt_temp = 15),
  "Narrow temperature tolerance" = list(curvature = 0.05),
  "Plagues at t = 50, 100, 150" = list(plague = TRUE),
  "Fixed sun (approach to equilibrium)" =
    list(lum_mode = "fixed", lum_fixed = 1.0, t_max = 120)
)

apply_scenario <- function(name) {
  modifyList(scenario_defaults, scenarios[[name]])
}

# ---- Parameter guide (also reused for the parameter export) ------------------
param_guide <- tribble(
  ~parameter,              ~symbol,   ~default,     ~range,          ~units,        ~meaning,
  "Black daisy albedo",    "A_black", "0.25",       "0 to 1",        "fraction",    "Reflectivity of ground covered by black daisies. Lower means more absorbed energy and local warming.",
  "White daisy albedo",    "A_white", "0.75",       "0 to 1",        "fraction",    "Reflectivity of ground covered by white daisies. Higher means more reflected energy and local cooling.",
  "Bare ground albedo",    "A_ground","0.50",       "0 to 1",        "fraction",    "Reflectivity of uncovered ground, midway between the two daisy types.",
  "Death rate",            "gamma",   "0.30",       "0 to 1",        "per time",    "Fraction of daisy area lost per unit time to natural mortality.",
  "Optimum temperature",   "T_opt",   "22.5",       "5 to 40",       "deg C",       "Local temperature at which the growth factor peaks at 1.",
  "Growth curve width",    "k",       "0.003265",   "0.001 to 0.06", "per (deg C)^2","Curvature of the growth parabola. 0.003265 gives growth from 5 to 40 C; larger values narrow the viable range.",
  "Heat absorption factor","q",       "20",         "0 to 40",       "deg C",       "Strength of local temperature contrast between daisy patches and the planetary mean.",
  "Solar flux constant",   "S",       "917",        "500 to 1500",   "W m-2",       "Baseline solar energy at Daisyworld's surface, before the luminosity factor.",
  "Luminosity (fixed)",    "L",       "1.0",        "0.6 to 1.8",    "relative",    "Relative brightness of the sun when held constant.",
  "Luminosity start / end","L",       "0.6 / 1.8",  "0.6 to 1.8",    "relative",    "Start and end brightness when the sun ramps over the run.",
  "Simulation length",     "t_max",   "200",        "20 to 400",     "time units",  "Total run length. One time unit is roughly 10 million years in the parable.",
  "Time step",             "dt",      "0.25",       "0.05 to 1",     "time units",  "Euler integration step. Smaller steps are more accurate but slower.",
  "Seed nudge",            "epsilon", "0.001",      "0 to 0.01",     "area per time","Trace seeding that lets a population start from zero. Without it the daisies never begin.",
  "Initial black cover",   "-",       "0",          "0 to 1",        "fraction",    "Starting area of black daisies.",
  "Initial white cover",   "-",       "0",          "0 to 1",        "fraction",    "Starting area of white daisies."
)

# ---- Plot builders -----------------------------------------------------------

plot_coverage <- function(df) {
  long <- df |>
    select(time, area_bare, area_black, area_white) |>
    pivot_longer(-time, names_to = "surface", values_to = "fraction") |>
    mutate(surface = recode(surface,
                            area_bare = "Bare ground",
                            area_black = "Black daisies",
                            area_white = "White daisies"),
           surface = factor(surface,
                            levels = c("Bare ground", "Black daisies", "White daisies")))

  ggplot(long, aes(time, fraction, fill = surface)) +
    geom_area(alpha = 0.9, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("Bare ground" = pal$bare,
                                 "Black daisies" = pal$black_d,
                                 "White daisies" = pal$white_d)) +
    scale_y_continuous(labels = percent_format(), expand = expansion()) +
    scale_x_continuous(expand = expansion()) +
    labs(title = "Surface coverage over time",
         subtitle = "Fraction of the planet covered by each surface type",
         x = "Time", y = "Coverage",
         caption = "Stacked areas sum to the whole planet.") +
    theme_daisy()
}

plot_temperature <- function(df, opt_temp, curvature) {
  half_width <- sqrt(1 / curvature)
  t_lo <- opt_temp - half_width
  t_hi <- opt_temp + half_width

  long <- df |>
    select(time, temp_planet, temp_dead) |>
    pivot_longer(-time, names_to = "series", values_to = "temp") |>
    mutate(series = recode(series,
                           temp_planet = "Living planet",
                           temp_dead = "Lifeless planet"),
           series = factor(series, levels = c("Living planet", "Lifeless planet")))

  ggplot(long, aes(time, temp)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = t_lo, ymax = t_hi,
             fill = pal$band, alpha = 0.55) +
    annotate("text", x = min(df$time), y = t_hi, hjust = 0, vjust = -0.4,
             label = "Habitable range for daisies", size = 3, colour = "#4B6B58") +
    geom_line(aes(colour = series, linetype = series), linewidth = 1) +
    scale_colour_manual(values = c("Living planet" = pal$living,
                                   "Lifeless planet" = pal$dead)) +
    scale_linetype_manual(values = c("Living planet" = "solid",
                                     "Lifeless planet" = "dashed")) +
    labs(title = "Planetary temperature over time",
         subtitle = "The lifeless planet tracks changes in the sun.",
         x = "Time", y = "Temperature (deg C)",
         caption = "Shaded band shows where daisy growth is possible.") +
    theme_daisy()
}

plot_albedo <- function(df, albedo_black, albedo_white, albedo_ground) {
  ggplot(df, aes(time, albedo)) +
    geom_hline(yintercept = c(albedo_black, albedo_ground, albedo_white),
               colour = "#C3CDCB", linetype = "dotted") +
    annotate("text", x = max(df$time), y = albedo_black, label = "black",
             hjust = 1.05, vjust = -0.4, size = 3, colour = pal$black_d) +
    annotate("text", x = max(df$time), y = albedo_ground, label = "bare",
             hjust = 1.05, vjust = -0.4, size = 3, colour = pal$bare) +
    annotate("text", x = max(df$time), y = albedo_white, label = "white",
             hjust = 1.05, vjust = -0.4, size = 3, colour = pal$white_d) +
    geom_line(colour = pal$teal, linewidth = 1.1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Planetary albedo over time",
         subtitle = "Whole-planet reflectivity as the daisy mix shifts",
         x = "Time", y = "Albedo (fraction reflected)",
         caption = "Dotted lines mark the three surface albedos. The planetary value moves between them as coverage changes.") +
    theme_daisy()
}

plot_growth_curve <- function(opt_temp, curvature, clamp_growth) {
  half_width <- sqrt(1 / curvature)
  temps <- seq(opt_temp - half_width - 5, opt_temp + half_width + 5,
               length.out = 400)
  curve_df <- tibble(
    temp = temps,
    growth = growth_factor(temps, opt_temp, curvature, clamp_growth)
  )
  ggplot(curve_df, aes(temp, growth)) +
    geom_hline(yintercept = 0, colour = "#C3CDCB") +
    geom_vline(xintercept = opt_temp, colour = pal$sun, linetype = "dashed") +
    geom_line(colour = pal$living, linewidth = 1.1) +
    annotate("text", x = opt_temp, y = 1.02, label = "optimum",
             colour = pal$sun, size = 3, vjust = 0) +
    scale_y_continuous(limits = c(NA, 1.08)) +
    labs(title = "Growth response",
         subtitle = "How local temperature sets the growth factor",
         x = "Local temperature (deg C)", y = "Growth factor") +
    theme_daisy()
}

# ---- UI ----------------------------------------------------------------------

controls_card <- card(
  card_header("Controls"),
  selectInput("scenario", "Load a scenario", choices = names(scenarios)),
  actionButton("apply_scenario", "Apply scenario",
               class = "btn-primary btn-sm mb-2"),
  helpText("Applying a scenario overwrites the settings below. Then adjust freely."),

  radioButtons("lum_mode", "The sun",
               choices = c("Brightens over the run" = "ramp",
                           "Held constant" = "fixed"),
               selected = "ramp"),
  conditionalPanel(
    "input.lum_mode == 'ramp'",
    sliderInput("lum_start", "Luminosity at start", 0.6, 1.8, 0.6, step = 0.05),
    sliderInput("lum_end", "Luminosity at end", 0.6, 1.8, 1.8, step = 0.05)
  ),
  conditionalPanel(
    "input.lum_mode == 'fixed'",
    sliderInput("lum_fixed", "Luminosity (constant)", 0.6, 1.8, 1.0, step = 0.05)
  ),

  sliderInput("albedo_black", "Black daisy albedo", 0, 1, 0.25, step = 0.01),
  sliderInput("albedo_white", "White daisy albedo", 0, 1, 0.75, step = 0.01),
  sliderInput("albedo_ground", "Bare ground albedo", 0, 1, 0.50, step = 0.01),
  sliderInput("death_rate", "Death rate", 0, 1, 0.30, step = 0.01),
  sliderInput("opt_temp", "Optimum growth temperature (deg C)", 5, 40, 22.5,
              step = 0.5),
  sliderInput("curvature", "Growth curve width (smaller = wider tolerance)",
              0.001, 0.06, 0.003265, step = 0.001),

  accordion(
    open = FALSE,
    accordion_panel(
      "Advanced settings",
      sliderInput("heat_absorp", "Heat absorption factor", 0, 40, 20, step = 1),
      sliderInput("solar_flux", "Solar flux constant (W m-2)", 500, 1500, 917,
                  step = 1),
      sliderInput("t_max", "Simulation length (time units)", 20, 400, 200,
                  step = 10),
      sliderInput("dt", "Time step", 0.05, 1, 0.25, step = 0.05),
      sliderInput("init_black", "Initial black cover", 0, 1, 0, step = 0.01),
      sliderInput("init_white", "Initial white cover", 0, 1, 0, step = 0.01),
      sliderInput("seed_nudge", "Seed nudge", 0, 0.01, 0.001, step = 0.001),
      checkboxInput("clamp_growth", "Hold growth factor at zero below viable range",
                    value = TRUE)
    ),
    accordion_panel(
      "Plagues",
      checkboxInput("plague", "Enable periodic plagues", value = FALSE),
      textInput("plague_times", "Plague times (comma separated)", "50, 100, 150"),
      sliderInput("plague_rate", "Death rate during a plague", 0, 1, 1.0,
                  step = 0.05),
      sliderInput("plague_width", "Plague duration (time units)", 1, 20, 4,
                  step = 1)
    )
  )
)

downloads_card <- card(
  card_header("Save your run"),
  layout_columns(
    col_widths = c(6, 6),
    downloadButton("dl_timeseries", "Time-step data (CSV)", class = "btn-sm"),
    downloadButton("dl_params", "Parameter settings (CSV)", class = "btn-sm"),
    downloadButton("dl_coverage", "Coverage plot (PNG)", class = "btn-sm"),
    downloadButton("dl_temperature", "Temperature plot (PNG)", class = "btn-sm"),
    downloadButton("dl_albedo", "Albedo plot (PNG)", class = "btn-sm")
  )
)

simulate_tab <- layout_sidebar(
  sidebar = sidebar(width = 340, controls_card, downloads_card),
  layout_columns(
    col_widths = c(4, 4, 4),
    value_box("Final planetary temperature", textOutput("vb_temp"),
              theme = value_box_theme(bg = pal$living, fg = "white")),
    value_box("Final living cover", textOutput("vb_cover"),
              theme = value_box_theme(bg = pal$teal, fg = "white")),
    value_box("Temperature regulation", textOutput("vb_reg"),
              theme = value_box_theme(bg = pal$ink, fg = "white"))
  ),
  layout_columns(
    col_widths = c(6, 6),
    card(full_screen = TRUE, plotOutput("plot_coverage", height = 340)),
    card(full_screen = TRUE, plotOutput("plot_temperature", height = 340))
  ),
  layout_columns(
    col_widths = c(6, 6),
    card(full_screen = TRUE, plotOutput("plot_albedo", height = 320)),
    card(full_screen = TRUE, plotOutput("plot_growth", height = 320))
  ),
  card(
    card_header("Time-step data"),
    DT::DTOutput("table")
  )
)

about_tab <- card(
  card_body(
    class = "mx-auto", style = "max-width: 820px;",
    htmltools::HTML(
      "<h2 style='font-family:Fraunces'>What Daisyworld is</h2>
      <p>Daisyworld is a thought-model of a planet whose only life is two kinds of
      daisies, black and white, growing on bare ground. It was built by Andrew
      Watson and James Lovelock (1983) to answer a fair challenge to the Gaia idea:
      how could life regulate a planet's climate without any shared plan or
      intention? Daisyworld shows that it can happen through ordinary feedback, with
      no coordination and no goal.</p>

      <h3 style='font-family:Fraunces'>The idea in one loop</h3>
      <p>The daisies differ only in color, so they differ in albedo (reflectivity).
      Black daisies absorb sunlight and warm the ground around them; white daisies
      reflect it and cool their surroundings. Each type grows best near an optimum
      temperature and stops growing when it is too cold or too hot. That single
      difference sets up a self-correcting loop:</p>
      <ul>
        <li>When the planet is cool, black daisies warm their patches and spread,
        which warms the whole planet.</li>
        <li>When the planet is warm, white daisies cool their patches and spread,
        which cools the whole planet.</li>
        <li>Either way, the mix of daisies shifts in the direction that pushes the
        temperature back toward what the daisies prefer.</li>
      </ul>
      <p>This is negative feedback, and it is emergent. No daisy senses the planet
      or acts for the common good. Regulation is a byproduct of many local
      responses to local temperature.</p>

      <h3 style='font-family:Fraunces'>How temperature is computed</h3>
      <p>The model follows an energy chain. Energy received depends on the sun's
      brightness (the luminosity factor) times a solar constant. Some is reflected
      according to the planet's albedo, and the rest is absorbed. At steady state,
      energy emitted equals energy absorbed, and the Stefan-Boltzmann law turns that
      balance into a temperature. Each daisy patch then feels a local temperature
      that is warmer or cooler than the planetary average depending on how its
      albedo compares with the whole planet.</p>

      <h3 style='font-family:Fraunces'>What to look for</h3>
      <p>In the standard case the sun slowly brightens. A lifeless planet would heat
      up steadily. The living planet instead holds close to the daisies' comfort
      range across a wide span of brightness, because the coverage keeps shifting
      from black toward white. Watch the green line stay inside the shaded band on
      the temperature plot while the dashed line climbs away. Near the end the sun
      grows too bright for even all-white cover to compensate, the daisies collapse,
      and regulation fails.</p>

      <h3 style='font-family:Fraunces'>Reading the model well</h3>
      <p>Daisyworld is deliberately simple. It assumes a well-watered planet, a fixed
      greenhouse, and no clouds, and it compresses a stellar lifetime into a short
      run. It is a way to reason about feedback and emergence, not a forecast for any
      real planet. That simplicity is the point: it makes the systems logic visible
      and gives everyone an accessible way into quantitative systems thinking,
      whatever their background. The visuals here use colorblind-safe colors backed
      by distinct line styles so the patterns read without depending on color.</p>

      <h3 style='font-family:Fraunces'>References</h3>
      <p style='font-size:0.9em'>Watson, A. J., and J. E. Lovelock (1983),
      Biological homeostasis of the global environment: the parable of Daisyworld,
      <em>Tellus B</em>, 35(4), 284-289.</p>
      <p style='font-size:0.9em'>Bice, D. (n.d.), Modeling Daisyworld, Penn State.
      Available at
      https://personal.ems.psu.edu/~dmb53/DaveSTELLA/Daisyworld/daisyworld_model.htm</p>
      <p style='font-size:0.9em'>SERC/Carleton (2005), Daisyworld: Stella Mac or PC.
      Available at
      https://serc.carleton.edu/introgeo/mathstatmodels/examples/DBDaisyworld.html</p>"
    )
  )
)

parameters_tab <- card(
  card_header("Parameter guide"),
  card_body(
    p("Every control on the Simulate page is listed below with its default, a
       typical range, units, and what it does. Two model choices are worth noting.
       The seed nudge is a trace amount of new daisy area added each step so a
       population can start from zero (without it the daisies never begin). Holding
       the growth factor at zero below its viable range keeps mortality from turning
       into unphysical negative growth outside the temperature window."),
    DT::DTOutput("param_table")
  )
)

ui <- page_navbar(
  title = "Daisyworld",
  theme = app_theme,
  fillable = FALSE,
  nav_panel("Simulate", simulate_tab),
  nav_panel("About the model", about_tab),
  nav_panel("Parameter guide", parameters_tab),
  nav_spacer(),
  nav_item(tags$span(style = "color:#5A6A72;font-size:0.85em;",
                     "J.G. Singley (Roger Williams University)"))
)

# ---- Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # Apply a scenario by updating every control.
  observeEvent(input$apply_scenario, {
    s <- apply_scenario(input$scenario)
    updateRadioButtons(session, "lum_mode", selected = s$lum_mode)
    updateSliderInput(session, "lum_start", value = s$lum_start)
    updateSliderInput(session, "lum_end", value = s$lum_end)
    updateSliderInput(session, "lum_fixed", value = s$lum_fixed)
    updateSliderInput(session, "t_max", value = s$t_max)
    updateSliderInput(session, "dt", value = s$dt)
    updateSliderInput(session, "albedo_black", value = s$albedo_black)
    updateSliderInput(session, "albedo_white", value = s$albedo_white)
    updateSliderInput(session, "albedo_ground", value = s$albedo_ground)
    updateSliderInput(session, "death_rate", value = s$death_rate)
    updateSliderInput(session, "opt_temp", value = s$opt_temp)
    updateSliderInput(session, "curvature", value = s$curvature)
    updateSliderInput(session, "heat_absorp", value = s$heat_absorp)
    updateSliderInput(session, "solar_flux", value = s$solar_flux)
    updateSliderInput(session, "init_black", value = s$init_black)
    updateSliderInput(session, "init_white", value = s$init_white)
    updateSliderInput(session, "seed_nudge", value = s$seed_nudge)
    updateCheckboxInput(session, "clamp_growth", value = s$clamp_growth)
    updateCheckboxInput(session, "plague", value = s$plague)
    updateTextInput(session, "plague_times", value = s$plague_times)
    updateSliderInput(session, "plague_rate", value = s$plague_rate)
    updateSliderInput(session, "plague_width", value = s$plague_width)
  })

  parse_times <- function(txt) {
    vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",")[[1]])))
    vals[!is.na(vals)]
  }

  # Collect current controls, debounced so dragging a slider settles first.
  current_params <- reactive({
    list(
      t_max = input$t_max, dt = input$dt,
      init_black = input$init_black, init_white = input$init_white,
      albedo_black = input$albedo_black, albedo_white = input$albedo_white,
      albedo_ground = input$albedo_ground, death_rate = input$death_rate,
      opt_temp = input$opt_temp, curvature = input$curvature,
      heat_absorp = input$heat_absorp, solar_flux = input$solar_flux,
      lum_mode = input$lum_mode, lum_start = input$lum_start,
      lum_end = input$lum_end, lum_fixed = input$lum_fixed,
      seed_nudge = input$seed_nudge, clamp_growth = input$clamp_growth,
      plague = input$plague, plague_times = parse_times(input$plague_times),
      plague_rate = input$plague_rate, plague_width = input$plague_width
    )
  }) |>
    debounce(400)

  sim_data <- reactive({
    do.call(simulate_daisyworld, current_params())
  })

  # Reactive plot objects, shared by the display and the download handlers.
  coverage_plot <- reactive(plot_coverage(sim_data()))
  temperature_plot <- reactive(
    plot_temperature(sim_data(), input$opt_temp, input$curvature)
  )
  albedo_plot <- reactive(
    plot_albedo(sim_data(), input$albedo_black, input$albedo_white,
                input$albedo_ground)
  )
  growth_plot <- reactive(
    plot_growth_curve(input$opt_temp, input$curvature, input$clamp_growth)
  )

  output$plot_coverage <- renderPlot(coverage_plot(), res = 96)
  output$plot_temperature <- renderPlot(temperature_plot(), res = 96)
  output$plot_albedo <- renderPlot(albedo_plot(), res = 96)
  output$plot_growth <- renderPlot(growth_plot(), res = 96)

  # Value boxes.
  output$vb_temp <- renderText({
    last <- tail(sim_data(), 1)
    sprintf("%.1f deg C", last$temp_planet)
  })
  output$vb_cover <- renderText({
    last <- tail(sim_data(), 1)
    percent(last$area_black + last$area_white, accuracy = 1)
  })
  output$vb_reg <- renderText({
    df <- sim_data()
    established <- df |> filter(area_black + area_white > 0.2)
    if (nrow(established) < 2) return("not established")
    dead_spread <- diff(range(established$temp_dead))
    live_spread <- diff(range(established$temp_planet))
    if (live_spread <= 0) return("fully buffered")
    sprintf("%.1fx steadier than a lifeless planet", dead_spread / live_spread)
  })

  # Data table.
  output$table <- DT::renderDT({
    sim_data() |>
      mutate(across(where(is.numeric), \(x) round(x, 3))) |>
      DT::datatable(options = list(pageLength = 10, scrollX = TRUE),
                    rownames = FALSE)
  })

  # Parameter guide table.
  output$param_table <- DT::renderDT({
    param_guide |>
      rename(Parameter = parameter, Symbol = symbol, Default = default,
             Range = range, Units = units, Meaning = meaning) |>
      DT::datatable(options = list(pageLength = 15, dom = "t"),
                    rownames = FALSE)
  })

  # A tidy snapshot of the current settings for export.
  params_snapshot <- reactive({
    p <- current_params()
    tibble(
      parameter = c("luminosity_mode", "luminosity_start", "luminosity_end",
                    "luminosity_fixed", "simulation_length", "time_step",
                    "initial_black_cover", "initial_white_cover",
                    "albedo_black", "albedo_white", "albedo_ground",
                    "death_rate", "optimum_temperature_C", "growth_curve_width",
                    "heat_absorption_factor", "solar_flux_constant_W_m2",
                    "seed_nudge", "clamp_growth_factor", "plague_enabled",
                    "plague_times", "plague_death_rate", "plague_width"),
      value = c(p$lum_mode, p$lum_start, p$lum_end, p$lum_fixed, p$t_max, p$dt,
                p$init_black, p$init_white, p$albedo_black, p$albedo_white,
                p$albedo_ground, p$death_rate, p$opt_temp, p$curvature,
                p$heat_absorp, p$solar_flux, p$seed_nudge, p$clamp_growth,
                p$plague, paste(p$plague_times, collapse = "; "),
                p$plague_rate, p$plague_width)
    )
  })

  stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

  # Downloads.
  output$dl_timeseries <- downloadHandler(
    filename = \() sprintf("daisyworld_timeseries_%s.csv", stamp()),
    content = \(file) write_csv(sim_data(), file)
  )
  output$dl_params <- downloadHandler(
    filename = \() sprintf("daisyworld_parameters_%s.csv", stamp()),
    content = \(file) write_csv(params_snapshot(), file)
  )
  output$dl_coverage <- downloadHandler(
    filename = \() sprintf("daisyworld_coverage_%s.png", stamp()),
    content = \(file) ggsave(file, coverage_plot(), width = 8, height = 5,
                             dpi = 150, bg = "white")
  )
  output$dl_temperature <- downloadHandler(
    filename = \() sprintf("daisyworld_temperature_%s.png", stamp()),
    content = \(file) ggsave(file, temperature_plot(), width = 8, height = 5,
                             dpi = 150, bg = "white")
  )
  output$dl_albedo <- downloadHandler(
    filename = \() sprintf("daisyworld_albedo_%s.png", stamp()),
    content = \(file) ggsave(file, albedo_plot(), width = 8, height = 5,
                             dpi = 150, bg = "white")
  )
}

shinyApp(ui, server)
