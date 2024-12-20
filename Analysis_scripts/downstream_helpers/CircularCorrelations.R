Fischer_Circular_Correlations <- function(rphases, sphases) {
  n1 <- length(rphases)
  n2 <- length(sphases)
  
  # Ensure lengths are equal
  if (n1 != n2) {
    stop("Lengths of rphases and sphases must be equal")
  }
  num <- n1
  
  # Wrap angles to [0, 2*pi]
  rphases <- (rphases + 2 * pi) %% (2 * pi)
  sphases <- (sphases + 2 * pi) %% (2 * pi)
  
  # Initialize totals
  numtot <- 0
  d1tot <- 0
  d2tot <- 0
  
  # Compute correlation components
  for (i in 1:(num - 1)) {
    for (j in (i + 1):num) {
      numeratorterm <- sin(sphases[i] - sphases[j]) * sin(rphases[i] - rphases[j])
      denomterm1 <- sin(sphases[i] - sphases[j])^2
      denomterm2 <- sin(rphases[i] - rphases[j])^2
      numtot <- numtot + numeratorterm
      d1tot <- d1tot + denomterm1
      d2tot <- d2tot + denomterm2
    }
  }
  
  # Compute Fisher correlation
  fischercor <- numtot / (sqrt(d1tot) * sqrt(d2tot))
  return(fischercor)
}

# Function to compute the circular mean
Circular_Mean <- function(phases) {
  sinterm <- sum(sin(phases))
  costerm <- sum(cos(phases))
  atan2(sinterm, costerm) %% (2 * pi) # atan2 ensures correct quadrant handling
}

# Function to compute the Jammalamadaka Circular Correlation
Jammalamadaka_Circular_Correlations <- function(rphases, sphases) {
  # Compute circular means
  rbar <- (2 * pi + Circular_Mean(rphases)) %% (2 * pi)
  sbar <- (2 * pi + Circular_Mean(sphases)) %% (2 * pi)
  
  # Compute numerator and denominators
  numtot <- sum(sin(rphases - rbar) * sin(sphases - sbar))
  d1tot <- sqrt(sum(sin(rphases - rbar)^2))
  d2tot <- sqrt(sum(sin(sphases - sbar)^2))
  
  # Jammalamadaka Circular Correlation
  Jammalamadaka <- numtot / (d1tot * d2tot)
  return(Jammalamadaka)
}

Circular_Rank_Phases <- function(rphases) {
  n <- length(rphases)
  rphases <- (rphases + 2 * pi) %% (2 * pi)
  rranks <- rank(rphases) # Using R's rank function
  rranks * 2 * pi / n
}

Jammalamadaka_Rank_Circular_Correlations <- function(rphases, sphases) {
  rphases <- Circular_Rank_Phases(rphases)
  sphases <- Circular_Rank_Phases(sphases)
  
  r_minus_s_bar <- atan2(sum(sin(rphases - sphases)), sum(cos(rphases - sphases))) %% (2 * pi)
  r_plus_s_bar <- atan2(sum(sin(rphases + sphases)), sum(cos(rphases + sphases))) %% (2 * pi)
  
  Ntot <- length(rphases)
  
  term1 <- cos(rphases - sphases - r_minus_s_bar)
  term2 <- cos(rphases + sphases - r_plus_s_bar)
  
  (1 / Ntot) * sum(term1) - (1 / Ntot) * sum(term2)
}

FindComponentAngles <- function(angle_sum, angle_diff) {
  rang <- (angle_sum + angle_diff) / 2
  sang <- (angle_sum - angle_diff) / 2
  c(rang, sang)
}

Jammalamadaka_Uniform_Circular_Correlations <- function(rphases, sphases) {
  # This is a modification of he Jammalamadka Circular Correlation
  # Described in Topics in Circular Statistics
  # It is required beacuse the circular average is not well defined with circular uniform data
  # This measure should only be used when 1 or both of the data sets being compared are uniform
  rphases <- rphases %% (2 * pi)
  sphases <- sphases %% (2 * pi)
  
  r_minus_s_bar <- atan2(sum(sin(rphases - sphases)), sum(cos(rphases - sphases))) %% (2 * pi)
  r_plus_s_bar <- atan2(sum(sin(rphases + sphases)), sum(cos(rphases + sphases))) %% (2 * pi)
  
  bars <- FindComponentAngles(r_plus_s_bar, r_minus_s_bar)
  rbar <- bars[1]
  sbar <- bars[2]
  
  numtot <- sum(sin(rphases - rbar) * sin(sphases - sbar))
  d1tot <- sqrt(sum(sin(rphases - rbar)^2))
  d2tot <- sqrt(sum(sin(sphases - sbar)^2))
  
  numtot / (d1tot * d2tot)
}

Fischer_Circular_CorrelationMeasures <- function(rphases, sphases) {
  rrankphases <- Circular_Rank_Phases(rphases)
  srankphases <- Circular_Rank_Phases(sphases)
  
  F <- Fischer_Circular_Correlations(rphases, sphases)
  FR <- Fischer_Circular_Correlations(rrankphases, srankphases)
  
  list(F = F, FR = FR)
}

Jammalamadaka_Circular_CorrelationMeasures <- function(rphases, sphases) {
  J <- Jammalamadaka_Circular_Correlations(rphases, sphases)
  JU <- Jammalamadaka_Uniform_Circular_Correlations(rphases, sphases)
  JR <- Jammalamadaka_Rank_Circular_Correlations(rphases, sphases)
  
  list(J = J, JU = JU, JR = JR)
}


# SANITY CHECKS
# # Example usage:
# rphases <- runif(n = 1000, 0, 1)
# sphases <- rnorm(n = 1000, 0, 1)
# 
# print("R Fischer Circular Correlation Measures:")
# Fischer <- Fischer_Circular_CorrelationMeasures(rphases, sphases)
# print(paste(Fischer$F, Fischer$FR))
# 
# print("R Jammalamadka Circular Correlation Measures")
# Jamma <- Jammalamadaka_Circular_CorrelationMeasures(rphases, sphases)
# print(paste(Jamma$J, Jamma$JU, Jamma$JR))
# 
# 
# #Calls julia script with trusted functions:
# rphases_julia <- paste(rphases, collapse = ",")
# sphases_julia <- paste(sphases, collapse = ',')
# 
# command <- paste("julia", "~/Box\\ Sync/Henry_stuff/AD_project/scROSMAP/scrosmap_newAtlas/sanity_checks/Circular_stats.jl", rphases_julia, sphases_julia)
# system(command)
