
# ----- Panel TWFE Staggered -----

TWFE_staggered <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 25 largest fixed effects to
  # treatment, but randomizes when exactly those treated units will get treated.
  if (periods == 2){
    g <- c(rep(0, 25), rep(2, 25))
  }else{
    g <- c(rep(0, 25), replicate(25, sample(2:periods, 1)))
  }
  
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$b <- effect_size
  
  df$y <- df$d * df$b + df$unit_fe + df$time_fe + df$e
  
  return(df)
  
}

# ----- Panel TWFE Staggered: Cohort Variation -----

TWFE_staggered_cohort <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = t,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this code assigns the units with the 25 largest fixed effects to
  # treatment. For eventually-treated units, it randomizes which periods
  # (beginning with period 2) will be periods in which units first get treated
  # (the "cohorts"), then uniformly at random assigns treated units to cohorts.
  if (periods >= 6){
    treat_periods <- sample(seq(from = 2, to = periods, by = floor(periods/5)))[1:5]
    g <- c(rep(0, 25), sample(rep(treat_periods, 5)))
  }
  
  if (periods < 6){
    a <- 2:periods
    m <- length(a)
    treat_periods_all <-  append(rep(a[1:(25 %% m)],
                                     each = ceiling(25 / m)),
                                 rep(a[((25 %% m) + 1):m],
                                     each = floor(25 / m)))
    g <- c(rep(0, 25), sample(treat_periods_all))
  }
  
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$b <- effect_size
  
  df$y <- df$d * df$b + df$unit_fe + df$time_fe + df$e
  
  return(df)
  
}

# ----- Staggered Panel - Violated Parallel Trends -----


staggered_violated_PT <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 25 largest fixed effects to
  # treatment, but randomizes when exactly those treated units will get treated.
  if (periods == 2){
    g <- c(rep(0, 25), rep(2, 25))
  }else{
    g <- c(rep(0, 25), replicate(25, sample(2:periods, 1)))
  }
    
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  # ADDING IN A TIME-VARYING COVARIATE $\nu$ TO VIOLATE PARALLEL TRENDS
  
  df$nu <- 0
  
  df$nu <- ifelse(test = (df$g == 0), yes = (-0.1 * df$t), no = df$nu)
  
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$b <- effect_size
  
  df$y <- df$d * df$b + df$unit_fe + df$time_fe + df$nu + df$e
  
  return(df)
  
}

# ----- Staggered Panel - HTE v.1 -----

staggered_hte1 <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 25 largest fixed effects to
  # treatment, but randomizes when exactly those treated units will get treated.
  g <- c(rep(0, 25), replicate(25, sample(2:periods, 1)))
  
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  # MODELING THE TREATMENT EFFECTS
  n_obs <- periods * 50
  
  df$b <- effect_size * replicate(n_obs, sample(c(0.9, 1.1), 1, replace = TRUE)) + rnorm(n_obs, 0, 0.04)
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$y <- df$d * df$b + df$unit_fe + df$time_fe + df$nu + df$e
  
  return(df)
  
}

# ----- Staggered Panel - HTE v.2 -----

staggered_hte2 <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 25 largest fixed effects to
  # treatment, but randomizes when exactly those treated units will get treated.
  g <- c(rep(0, 25), replicate(25, sample(2:periods, 1)))
  
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  # MODELING THE TREATMENT EFFECTS
  
  n_obs <- periods * 50
  
  df$b <- effect_size - 0.01 * df$unit_fe  + rnorm(n_obs, 0, 0.001)
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$y <- df$d * df$b + df$unit_fe + df$time_fe + df$nu + df$e
  
  return(df)
  
}

# ----- Staggered Panel - HTE v.3 -----

staggered_hte3 <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 25 largest fixed effects to
  # treatment, then assigns g in a way that is correlated with the fixed effects
  g <- assign_g_cor_i_FEs(periods = periods)

  # MODELING THE TREATMENT EFFECTS
  # Tau here is monotonically increasing in g and, for treated units, ranges
  # from 0.02 to 0.5
  tau <-  ifelse(g > 0, i - 25, 0) * 0.02 
  
  # PUTTING TOGETHER THE DATA FRAME
  df <- merge(data.frame(i = i, g = g, tau = tau), data.frame(t = t), by = NULL)
  
  

  # IDENTIFYING IN-TREATMENT UNITS
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  # ADDING IN RAW AVERAGE TAU
  df$mean_tau_raw <- mean(tau[g > 0])
  
  # ADDING IN TIME-WEIGHTED AVERAGE TAU
  df$mean_tau_t_weighted <- mean(df$tau[df$d])
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$y <- df$d * df$tau + df$unit_fe + df$time_fe + df$e
  
  return(df)
  
}


# ----- Staggered Panel - Dynamic TEs -----

staggered_dte <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 25 largest fixed effects to
  # treatment, then assigns g in a way that is correlated with the fixed effects
  g <- assign_g_cor_i_FEs(periods = periods)
  
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  # MODELING THE TREATMENT EFFECTS
  df$tau <-  ifelse(i > 25, t - 1, 0) * 0.02
  
  # ADDING IN RAW AVERAGE TAU
  df$mean_tau_raw <- mean(df$tau[g > 0])
  
  # ADDING IN TIME-WEIGHTED AVERAGE TAU
  df$mean_tau_t_weighted <- mean(df$tau[df$d])
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$y <- df$d * df$tau + df$unit_fe + df$time_fe + df$e
  
  return(df)
  
}


# ----- Panel TWFE Staggered: Hate Crime Law Implementation -----

TWFE_staggered_hate_crime_law <- function(periods, effect_size){
  
  # INITIALIZING VECTORS
  
  i <- 1:50
  t <- 1:periods
  
  # UNIT/TIME FACTORS 
  
  # Year fixed effects are a simulated random walk. 
  time <- data.frame(
    t = 1:periods,
    time_fe = random_walk_simulated(n = periods)
  )
  
  # Unit fixed effects are literally just index values 1:50 times 0.1. This is
  # quite arbitrary.
  unit <- data.frame(
    i = 1:50,
    unit_fe = 1:50 * 0.1
  )
  
  # TREATMENT ASSIGNMENT 
  
  # The unit fixed effects increase with the index of the unit;
  # this line of code assigns the units with the 48 largest fixed effects to
  # treatment, but randomizes when exactly those treated units will get treated.
  g <- c(rep(0, 2), replicate(48, sample(2:periods, 1)))
  
  df <- merge(data.frame(i = i, g = g), data.frame(t = t), by = NULL)
  
  df$d <- 0
  
  df$d <- ifelse(test = (df$g != 0 & df$t >= df$g), yes = 1, no = df$d)
  
  # ADDING IN THE FACTORS 
  
  df <- merge(df, unit, by = "i")
  df <- merge(df, time, by = "t")
  
  
  # OUTCOMES 
  
  df$e <- rnorm(nrow(df))
  
  df$b <- effect_size
  
  df$y <- df$d * df$b + df$unit_fe + df$time_fe + df$e
  
  return(df)
  
}


