library(mclust)
library(tidyverse)

messagef <- function(...) message(sprintf(...))
printf <- function(...) print(sprintf(...))

BUR_areas <- 
  tribble(~xmin, ~xmax, ~ymin, ~ymax, ~type,
        -2.0,  -.5, 0, Inf, "snap",
        -.5,   .5, 0, Inf, "even",
        .5,  2.0, 0, Inf, "swung")

sample_by_id <-function(data, id_var = "id", frac = NULL, n = NULL, ...){
  ids <- unique(data[[id_var]])
  l <- length(ids)
  if(is.null(frac)){
    if(is.null(n)){
      stop("Either frac and n must be set")
    }
  }
  else{
    if(is.null(frac)){
      stop("Either frac and n must be set")
    }
    n <- round(length(ids) * frac)
  }
  n <- min(n, l)
  sample_ids <- sample(unique(data[[id_var]]), size = n, ...)  
  data %>% filter(!!sym(id_var) %in% sample_ids) 
}

limit_rnorm <- function(n, mean, sd, lower = NULL, upper = NULL, exact_n = T){
  stopifnot(n > 0)
  x <- rnorm(n, mean, sd)
  if(!is.null(lower)){
    x <- x[x>= lower]
  }
  if(!is.null(upper)){
    x <- x[x<= upper]
  }
  if(!exact_n){
    return(x)
  }
  l <- length(x)
  while(l < n ){
    #messagef("Filling up vector, %d elements of %d missing (m = %.2f, s = %.2f)", n - l, n, mean, sd)
    fill <- limit_rnorm(n, mean, sd, lower, upper, exact_n = F)
    if(length(fill) == 0){
      stop("Check parameters")
    }
    x <- c(x, fill)
    l <- length(x)
    
  }
  #messagef("*** DONE ***")
  x[1:n]
}
simulate_swing_eighths <- function(beat_dur = .5, 
                                   onset1_sigma = .03, 
                                   onset2_sigma = .03, 
                                   beat_sigma = .03, 
                                   sr = 2, n = 1000){
  beat1_onset <- rnorm(n, mean = 0, beat_sigma)
  beat2_onset <- rnorm(n, mean = beat_dur, beat_sigma)
  onset1 <- map_dbl(rnorm(n, 0, onset1_sigma), function(x) max(-.03, x))
  onset2 <- map_dbl(rnorm(n, beat_dur/(1 + 1/sr), onset2_sigma), function(x) min(x, beat_dur-.02))
  ioi <- onset2-onset1
  beat_dur <- beat2_onset - beat1_onset
  ooi <- beat_dur - ioi
  
  tibble(nom_sr = sr, 
         beat_dur = beat_dur, 
         onset1 = onset1,  
         onset2 = onset2, 
         ioi = ioi, 
         ooi = ooi,
         sr = ioi/ooi)
}

simulate_swing_eighths2 <- function(sr = 2, beat_dur = .5, sigma1 = .03, n = 10000, min_sr = 1, max_sr = 5, type = "interval"){
  m_ioi1 <- beat_dur*sr/(sr + 1)
  m_ioi2 <- beat_dur/(sr + 1)
  if(type == "point"){
    onset0 <- rnorm(n = n, mean = 0, sd =  sigma1) 
    onset1 <- rnorm(n = n, mean = m_ioi1, sd =  sigma1) 
    onset2 <- rnorm(n = n, mean = m_ioi1 + m_ioi2, sd =  sigma1)
    
  }
  else{
    onset0 <- limit_rnorm(n = n, mean = 0, sd =  sigma1, lower = -sigma1, upper = sigma1) 
    onset1 <- limit_rnorm(n = n, mean = m_ioi1, sd =  sigma1, lower = -sigma1 + m_ioi1, upper = sigma1 + m_ioi1) 
    onset2 <- limit_rnorm(n = n, mean = m_ioi1 + m_ioi2, sd =  sigma1, lower = -sigma1 +m_ioi1 + m_ioi2, upper = sigma1 + m_ioi1 + m_ioi2)
    
  }
  ioi1 <- onset1 - onset0
  ioi2 <- onset2 - onset1
  
  #if((ioi1/ioi2) <0 ){
  #  browser()
  #}
  tibble(nom_sr = sr, 
         beat_dur = beat_dur,
         nom_longer = m_ioi1,
         nom_shorter = m_ioi2,
         onset0 = onset0, 
         onset1 = onset1, 
         onset2 = onset2,
         ioi1 = ioi1,  
         ioi2 = ioi2,
         sr = ioi1/ioi2, 
         log2_sr = log2(sr)) %>% filter(sr >= min_sr, sr <= max_sr)
}

estimate_sr_bias <- function(sr = 2, beat_dur = .5, sigma = .03,  n = 100){
  #browser()
  map_dfr(sr, function(s){
    map_dfr(beat_dur, function(b){
      map_dfr(sigma, function(sig){
          messagef("Estimating sr = %.2f with beat_dur = %.2f, sigma = %.3f, n = %d", s, b, sig, n)
          srs <- (limit_rnorm(n, s*b/(1+s), sig, 0)/limit_rnorm(n, b/(1 + s), sig, 0))
          srs <- srs[srs >0 & srs < 4 ]
          tibble(beat_dur = b, nom_sr  = s, sr = mean(srs, na.rm = T), sigma = sig)
      })
    })
  })
}

estimate_sr_bias2 <- function(sr = 2, beat_dur = .5, sigma = .03,  n = 1000){
  #browser()
  map_dfr(sr, function(s){
    map_dfr(beat_dur, function(b){
      map_dfr(sigma, function(sig){
        messagef("Estimating sr = %.2f with beat_dur = %.2f, sigma = %.3f, n = %d", s, b, sig, n)
        #srs <- (limit_rnorm(n, s*b/(1+s), sig, 0)/limit_rnorm(n, b/(1 + s), sig, 0))
        srs <- simulate_swing_eighths2(s, b, sig, min_sr = 0, max_sr = 4, n = n) %>% 
          summarise(sr = mean(sr), log2_sr = mean(log2(sr))) 
        #browser()
        tibble(beat_dur = b, nom_sr  = s, sr = srs$sr, log2_sr = srs$log2_sr, d_sr = sr - nom_sr, sigma = sig)
      })
    })
  })
}

interpolate_sr <- function(data, sigma = .025, precision = .1, max_sr = 3, n = 500){
  data <- data %>% 
    filter(!is.na(sr), sr > 1, sr < max_sr) %>% 
    group_by(id, avgtempo) %>% 
    summarise(sr = mean(sr, na.rm = T),
              sr_median = median(sr, na.rm = T),
              sd_sr = sd(sr, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(beat_dur = 60/avgtempo)
  beat_durs <- 60/unique(data$avgtempo)
  esr <- estimate_sr_bias(sr = seq(1, max_sr, precision), 
                          beat_dur = beat_durs, 
                          sigma = sigma, 
                          n = n) %>% mutate(avgtempo = 60/beat_dur)
  sr_model <- esr %>% lm(nom_sr ~ sr + beat_dur, data = .)
  predict_df <- tibble(id = data$id, predict_sr = predict(sr_model, newdata = data))
  browser()
  esr %>% 
    rename(sr_est = sr) %>%  
    left_join(data %>% select(-beat_dur), by = "avgtempo") %>% 
    filter(!is.na(id)) %>% 
    mutate(d_sr = sr_est - sr) %>% 
    group_by(id) %>% 
    filter(abs(d_sr) == min(abs(d_sr))) %>% 
    ungroup() %>% 
    left_join(predict_df, by = "id") %>% 
    rename(measured_sr = sr)
}

read_swing_ratios_raw <- function(fname = "data/swing_ratios_raw.csv"){
  srw <- read.csv(fname, header = T, stringsAsFactors = F, sep = ";") %>% as_tibble()
  srw <- srw %>% mutate(melid = as.integer(factor(id)))
  srw <- srw %>% mutate(beat_id = sprintf("%s_%s_%s", melid, bar, beat))
  srw <- srw %>% group_by(beat_id) %>% mutate(events_per_beat = n()) %>% ungroup()
  srw <- srw %>% mutate(diff_tatum = c(NA, diff(tatum))) 
  srw <- srw %>% 
    group_by(id) %>% 
    mutate(ioi = c(diff(onset), NA),
           abs_mcm = (bar - 1)*48 + mcm_48,
           next_ioi = lead(ioi), 
           next_abs_mcm = lead(abs_mcm),
           ioi_mcm = c(diff(abs_mcm), NA),
           next_tatum = lead(tatum)) %>% 
    ungroup()
  #%>% filter(division %in% c(2,3,4), events_per_beat == 2, tatum == 1 | diff_tatum == division -1) %>% select(bar, beat, tatum, division, events_per_beat, beat_id, onset)  
  srw
}

read_swing_ratios_sum <- function(fname = "c:/MeloSpyGUI/analysis/feature+viz/swing_ratios_sum.csv"){
  srs <- read.csv(fname, header = T, stringsAsFactors = F, sep = ";") %>% as_tibble()
}

read_single_beattrack <- function(bt_file){
  id <- gsub("_beattrack.csv", ".sv", basename(bt_file), fixed = T) 
  read.csv(bt_file, sep = ";", header = T, stringsAsFactors = F) %>% 
    as_tibble() %>% 
    mutate(id = id) %>% select(onset, duration, bar, beat, id)
  
}

read_beattracks <- function(bt_dir = "./data/beattracks"){
  bts <- list.files(bt_dir, pattern = "*.csv", full.names = T)
  #browser()
  map_dfr(bts, read_single_beattrack) %>% 
    mutate(melid = as.integer(factor(id))) %>% 
    mutate(beat_id = sprintf("%s_%s_%s", melid, bar, beat, id))
}

get_swing_ratios <- function(data, beat_divisions = c(2, 3, 4)){
  if(!is.null(beat_divisions)){
    data <- data %>% filter(division %in% beat_divisions)
  }
  data %>% 
    filter(events_per_beat == 2, tatum == 1 |tatum == division ) %>%  
    select(bar, beat, tatum, division, events_per_beat, beat_id, onset, ioi, next_ioi, beat_onset, beat_dur)  %>% 
    group_by(beat_id) %>% 
    mutate(sr = ioi/pmin(beat_dur-ioi, next_ioi),
           sr2 = ioi/(beat_dur - ioi), 
           n = n()) %>% 
    ungroup() %>% 
    filter(!is.na(sr), n == 2, tatum == 1)  %>% 
    select(beat_id, sr, sr2)
}

get_swing_ratios2 <- function(data, beat_divisions = c(2, 3, 4)){
  if(!is.null(beat_divisions)){
    data <- data %>% filter(division %in% beat_divisions)
  }
  #browser()
  data %>% 
    filter(events_per_beat == 2, 
           #ioi <= beat_dur, next_ioi <= beat_dur,
           ioi + next_ioi <= beat_dur + .03,
           tatum == 1 |(tatum == division)) %>%  
    group_by(beat_id) %>% 
    mutate(sr = ioi/next_ioi,
           #sr2 =  ioi/pmin(beat_dur-ioi, next_ioi),
           #sr3 = ioi/(beat_dur - ioi),
           n = n()) %>% 
    ungroup() %>% 
    filter(!is.na(sr), n == 2, tatum == 1)  %>% 
    select(beat_id, sr)
}

name_sr_cluster <- function(gmm_data){
  gmm_data %>%
    mutate(sr_type = case_when(
      log2_sr < -.5 ~ "snap", 
      log2_sr > .5 ~ "swung", 
      TRUE ~ "even"))
}

fit_all_GMMS <- function(data, min_n_sr = 15){
  solo_gmm <- fit_GMM(data, min_n_sr = min_n_sr) %>% 
    left_join(data %>% distinct(id, avgtempo, style, tempoclass, recordingyear, decade, swing_feel), by = "id")
  performer_gmm <- fit_GMM(data, group_var = "performer", min_n_sr = min_n_sr)
  decade_gmm <- fit_GMM(data, "decade", min_n_sr = min_n_sr)
  swing_gmm <- fit_GMM(data, "swing_feel", min_n_sr = min_n_sr)
  style_gmm <- fit_GMM(data, "style", min_n_sr = min_n_sr) 
  tempo_gmm <- fit_GMM(data, "tempoclass", min_n_sr = min_n_sr) 
  assign("solo_gmm", solo_gmm, globalenv())
  assign("performer_gmm", performer_gmm, globalenv())
  assign("decade_gmm", decade_gmm, globalenv())
  assign("swing_gmm", swing_gmm, globalenv())
  assign("style_gmm", style_gmm, globalenv())
  assign("tempo_gmm", tempo_gmm, globalenv())
}

fit_GMM <- function(data, group_var = "id", max_sr = 4, min_n_sr = 15, use_bootstrap = F){
  if(use_bootstrap){
    set.seed(666)
  }
  if(is.null(group_var) || is.na(group_var) || !is.character(group_var) || nchar(group_var) == 0){
    group_var <- "dummy"
    data$dummy <- 1
  }
  data <- data %>% 
    group_by(!!sym(group_var)) %>% 
    mutate(n_sr = sum(!is.na(sr) & tatum == 1)) %>% 
    ungroup()
  data <- data %>% filter(!is.na(sr), tatum == 1, abs(sr) <= max_sr, n_sr >= min_n_sr)
  groups <- unique(data[[group_var]])
  purrr::map_dfr(groups, function(g){
    #browser()
    log_sr <- data %>% filter(!!sym(group_var) == g) %>% pull(sr) %>% log2()
    a <- mclust::Mclust(log_sr)
    cis <- list()
    if(use_bootstrap) {
        cis <- tryCatch({
        aa <- MclustBootstrap(a)
        log2_sr_ci <- summary(aa, what = "ci")$mean %>% 
          as_tibble() %>% 
          t() %>% 
          as_tibble() %>% 
          rename(log2_sr_lo = 1, log2_sr_hi = 2)
        sr_ci <- 2^summary(aa, what = "ci")$mean %>% 
          as_tibble() %>% 
          t() %>% 
          as_tibble() %>% 
          rename(sr_lo = 1, sr_hi = 2)  
        list(log2_sr_ci = log2_sr_ci, sr_ci = sr_ci )
      }, 
      error = function(e){
        if(a$G == 1){
          m <- summary(aa, what = "ave")$mean[1,1]
          se <- summary(aa, what = "se")$mean[1,1]
          return(list(log2_sr_ci  = tibble(log2_sr_lo = m - se, log2_sr_hi = m + se),
                      sr_ci =  tibble(sr_lo = 2^(m - 1.96*se), sr_hi = 2^(m + 1.96*se))))
        }
        else{
          return(list(log2_sr_ci =  tibble(log2_sr_lo = rep(NA, a$G), log2_sr_hi = rep(NA, a$G)),
          sr_ci = tibble(sr_lo = rep(NA, a$G), sr_hi = rep(NA, a$G))))
        }
        browser()
      })

    }
    #if(a$G > 3){
    #  print(g)
    #  plot(a)
    #}
    ret <- tibble(!!sym(group_var) := g, 
           num_clust = a$G, 
           n = length(log_sr), 
           a$parameters$mean %>% 
             as_tibble() %>% 
             rename(log2_sr = value), 
           a$parameters$pro %>% 
             as_tibble() %>% 
             rename(prob = value) ) %>%
      mutate(sr = 2^log2_sr) %>% 
      arrange(sr) %>% 
      mutate(cluster_id = 1:nrow(.))
    if(use_bootstrap){
      ret <- ret %>% bind_cols(cis[["log2_sr_ci"]], cis[["sr_ci"]])
    }
    ret
  }) %>% name_sr_cluster()
}

setup_workspace <- function(recalc_GMM = F){
  srw <- read_swing_ratios_raw()
  #srs <- read_swing_ratios_sum()
  bts <- read_beattracks()
  srw <- srw %>% left_join(bts %>% select(beat_id, beat_onset = onset, beat_dur = duration), by = "beat_id") 
  sr <- get_swing_ratios2(srw, beat_divisions = NULL) 
  srw <- srw %>% left_join(sr, by = "beat_id") %>% 
    group_by(id) %>%
    mutate(n_sr = sum(!is.na(sr) & tatum == 1),
           n_sr_ls = sum(!is.na(sr) & tatum == 1 & sr > 1)) %>% 
    ungroup() 
  srw <- srw %>% 
    mutate(swing_feel = factor(rhythmfeel %in% c("SWING", "TWOBEAT"), labels = c("Non-Swing", "Swing Feel"))) %>%
    mutate(decade = factor(floor((recordingyear - 1900)/10)*10)) %>% 
    group_by(id) %>%   
    mutate(num_notes = n(), nom_first = ioi/beat_dur, nom_second = next_ioi/beat_dur) %>% 
    ungroup()
  srw_f <- srw %>% filter(!is.na(sr), abs(sr) < 4, tatum == 1) %>%   
    group_by(id) %>%
    mutate(n_sr = n(),
           n_sr_ls = sum(sr > 1),
           n_sr_sl = sum(sr < 1)) %>% 
    ungroup() 
  srw_f <- srw_f %>%  mutate(log2_sr = log2(sr)) %>% name_sr_cluster()
  srw_f_sum <-  srw_f %>% 
    group_by(id, sr_type) %>% 
    summarise(mean_log2_sr = mean(log2_sr)) %>% 
    ungroup() %>% 
    left_join(srw_f %>% distinct(id, .keep_all = T) %>% select(-sr_type), by ="id") 
  assign("srw", srw, globalenv())
  assign("srw_f", srw_f, globalenv())
  assign("srw_f_sum", srw_f_sum, globalenv())
  if(recalc_GMM){
    fit_all_GMMS(srw_f)
  }
}

sr_density_plot <- function(data, min_n_sr = 50){
  q <- 
    data %>% 
    filter(n_sr >= min_n_sr) %>% 
    ggplot() 
  q <- q + geom_rect(data = BUR_areas,
                     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type), 
                     alpha = .1)
  q <- q + geom_density(aes(x = log2_sr, y = ..density..), 
                        fill = "lightblue4", colour ="black", alpha = .2) 
  q <- q + theme_bw() 
  #q <- q + facet_wrap(~melid, scale = "free_y") 
  q <- q + geom_segment(aes(x = log2_sr, xend = log2_sr, y = 0, yend = prob),
                        colour ="blue", linetype = "solid")
  q <- q + geom_vline(xintercept = 0, linetype = "dashed") 
  q <- q + geom_vline(xintercept = 1, linetype = "dashed") 
  q <- q + labs(x ="Log2 BUR", y = "Density")
  q <- q + theme(legend.title = element_blank())
  q 
}

global_GMM_analysis <- function(solo_gmm = NULL, swingers = NULL, min_n_sr = 15){
  if(is.null(solo_gmm)){
    solo_gmm <- fit_GMM(srw_f, min_n_sr = min_n_sr)
    swingers <- solo_gmm %>% filter(abs(2-sr) < .5)  %>%  pull(id) %>% unique()
    swingers_melids <- srw %>% filter(id %in% swingers) %>% distinct(id, melid)

  }
  browser()
  srw_f <- srw_f %>% left_join(solo_gmm %>% select(cluster_log2_sr = log2_sr, prob, id), by = "id") 
  q <- 
    srw_f %>% 
      filter(id %in% swingers, n_sr >= min_n_sr) %>% 
      ggplot() 
  q <- q + geom_rect(data = BUR_areas,
                     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type), 
                     alpha = .1)
  q <- q + geom_density(aes(x = log2_sr, y = ..density..), 
                        fill = "lightblue4", colour ="black", alpha = .2) 
  q <- q + theme_bw() 
  q <- q + facet_wrap(~melid, scale = "free_y") 
  q <- q + geom_segment(aes(x = cluster_log2_sr, xend = cluster_log2_sr, y = 0, yend = prob),
                        colour ="blue", linetype = "solid")
  #q <- q + geom_rect(xmin = -2, xmax = -.5, ymin = 0, ymax = Inf, fill = "blue", alpha = .1)
  #q <- q + geom_vline(xintercept = GMM$log_sr[2], colour ="blue", linetype = "dashed", size = 1) 
  #q <- q + geom_vline(xintercept = GMM$log_sr[3], colour ="blue", linetype = "dashed", size = 1)
  q <- q + geom_vline(xintercept = 0, linetype = "dashed") 
  q <- q + geom_vline(xintercept = 1, linetype = "dashed") 
  #q <- q + geom_vline(aes(xintercept = m), colour ="blue")
  q <- q + labs(x ="Log2 BUR", y = "Density")
  q <- q + theme(legend.title = element_blank())
  ggsave("figs/sr_distribution_150.png", dpi = 600)
  #sprint(q)
  #return()
  a <- mclust::Mclust(srw_f %>%  pull(log2_sr))
  #browser()
  GMM <- 
    tibble(a$parameters$mean %>% as_tibble() %>% rename(log2_sr = value), 
         a$parameters$pro %>% as_tibble() %>% rename(prob = value) ) %>% 
    mutate(sr = 2^log2_sr)
  #browser()
  q <- 
    srw_f %>% 
    ggplot() 
  q <- q + geom_histogram(aes(x = log2_sr, y = ..density..), 
                          fill = "lightblue4", colour ="black", binwidth = .05, alpha = .5)
  q <- q + geom_density(aes(x = log2_sr, y = ..density..),
                        fill = "lightblue4", colour ="black", alpha = .2) 
  q <- q + theme_bw() 
  q <- q + geom_rect(data = BUR_areas,
                     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax, fill = type), 
                     alpha = .1)
  
  q <- q + geom_vline(xintercept = GMM$log2_sr[1], colour ="black", linetype = "dashed", size = 1) 
  q <- q + geom_vline(xintercept = GMM$log2_sr[2], colour ="black", linetype = "dashed", size = 1) 
  #q <- q + geom_vline(xintercept = GMM$log2_sr[3], colour ="blue", linetype = "dashed", size = 1)
  q <- q + theme(legend.title = element_blank())
  q <- q + labs(x ="Log2 BUR", y = "Density")
  print(q)
  ggsave("figs/sr_distribution_all.png", dpi = 600)
  #solo_gmm <- fit_GMM(srw_f, group_var = "id", min_n_sr = 30) 
  
}

tempo_analysis_GMM <- function(solo_gmm){
  deg <- 1
  solo_gmm <- solo_gmm %>% mutate(tempo_var = avgtempo)  
  q <- 
    solo_gmm %>% 
    ggplot(aes(x = tempo_var, y = log2_sr, colour = sr_type)) 
  q <- q + geom_point() 
  q <- q + geom_smooth(method = "lm", formula = "y~poly(x,1)") 
  q <- q + theme_bw() 
  q <- q + theme(legend.title = element_blank()) 
  q <- q + labs(x = "Avg. Tempo (bpm)", y = "log2-BUR") 
  print(q)
  ggsave("figs/sr_components_by_tempo.png", dpi = 600)

  fit <- lmer(log2_sr ~ poly(tempo_var, deg) + (1|sr_type), data = solo_gmm)

  print(broom.mixed::glance(fit))  
  print(broom.mixed::tidy(fit)) 
  print(anova(fit))
  fit <- lm(log2_sr ~ poly(tempo_var,deg), data = solo_gmm %>% filter(sr_type == "swung"))
  print(broom::glance(fit))  
  print(broom::tidy(fit)) 
  print(anova(fit))
  srw_f_sum <- srw_f_sum %>%  
    distinct(id, .keep_all = T) %>% 
    left_join(srw_f %>% 
                distinct(id, avgtempo), by = "id") %>% 
    mutate(tempo_var = avgtempo)  
  #browser()
  p <- 
    srw_f_sum %>% 
    ggplot(aes(x = tempo_var, y = mean_log2_sr, colour = sr_type)) 
  p <- p + geom_point() 
  p <- p + geom_smooth(method = "lm", formula = "y~poly(x,1)") 
  p <- p + theme_bw() 
  p <- p + theme(legend.title = element_blank()) 
  p <- p + labs(x = "Avg. Tempo (bpm)", y = "Mean log2-BUR")
  print(p)
  ggsave("figs/sr_group_by_tempo.png", dpi = 600)
  fit <- lmer(mean_log2_sr ~ poly(tempo_var, deg) + (1 +  poly(tempo_var, deg)|sr_type), data = srw_f_sum)
  print(broom::glance(fit))  
  print(broom::tidy(fit)) 
  print(anova(fit))
  #browser()
  fit <- lm(mean_log2_sr ~ poly(tempo_var, deg), data = srw_f_sum %>% filter(sr_type == "swung")) 
  print(broom::glance(fit))  
  print(broom::tidy(fit)) 
  print(anova(fit))
  deg <- 2
  srw_f_sum2 <- srw_f %>% 
    group_by(id) %>% mutate(mean_log2_sr = mean(log2_sr)) %>% 
    ungroup() %>% 
    distinct(id, mean_log2_sr, avgtempo, sr_type) %>% 
    mutate(tempo_var = avgtempo)  
  browser()
  p <- 
    srw_f_sum2 %>% 
    ggplot(aes(x = tempo_var, y = mean_log2_sr)) 
  p <- p + geom_point() 
  p <- p + geom_smooth(method = "lm", formula = "y~poly(x,1)") 
  p <- p + theme_bw() 
  p <- p + theme(legend.title = element_blank()) 
  p <- p + labs(x = "Avg. Tempo (bpm)", y = "Mean log2-BUR")
  print(p)
  ggsave("figs/sr_by_tempo.png", dpi = 600) 
  fit <- lm(mean_log2_sr ~ poly(avgtempo, deg), data = srw_f_sum2) 
  print(broom::glance(fit))  
  print(broom::tidy(fit)) 
  print(anova(fit))
}

get_counts <- function(data){
  l <- list()
  l[["all"]] <- nrow(data)
  l[["two_events"]] <- data %>% filter(events_per_beat == 2) %>% nrow()
  data <- data %>% filter(!is.na(sr), tatum == 1)
  l[["abs_sr_lt_4"]] <- data %>% filter(abs(sr) <= 4) %>% nrow()
  l[["min_gt_sext"]] <- data %>% filter(nom_second > 1/6 | nom_first > 1/6) %>% nrow()
  data <- data %>% filter(!is.na(sr), tatum == 1, abs(sr) <= 4)
  l[["nsr_15_all"]] <- data %>% filter(n_sr > 15) %>% nrow()
  l[["nsr_30_all"]] <- data %>% filter(n_sr > 30) %>% nrow()
  l[["nsr_15_sl"]] <- data %>% filter(sr >= 1, n_sr >= 15) %>% nrow()
  l[["nsr_30_sl"]] <- data %>% filter(sr >= 1, n_sr >= 30) %>% nrow()
  l[["nsr_15_ls"]] <- data %>% filter(sr < 1, n_sr >= 15) %>% nrow()
  l[["nsr_30_ls"]] <- data %>% filter(sr < 1, n_sr >= 30) %>% nrow()
  l[["nsr_15_sl_solos"]] <- data %>% filter(sr >= 1, n_sr >= 15) %>% pull(id) %>% n_distinct()
  l[["nsr_30_sl_solos"]] <- data %>% filter(sr >= 1, n_sr >= 30) %>% pull(id) %>% n_distinct()
  l[["nsr_100_all_solos"]] <- data %>% filter(n_sr >= 100) %>% pull(id) %>% n_distinct()
  l[["nsr_tempoclass"]] <- data %>% count(tempoclass, n_sr > 15)
  #l %>% as_tibble()
  l
}
se <- function(x, na.rm = T, ...){
  if(length(x) <= 1){
    return(NA)
  }
  sd(x, na.rm = na.rm, ...)/sqrt(length(x))
}

swing_ratio_analysis <- function(data, n_select = 20, max_sr = 1000){
  #t <- data %>% 
  #  group_by(id) %>% 
  #  mutate(num_notes = n()) %>% 
  #  distinct(id, .keep_all = T) 
  tmp <- data %>% 
    group_by(id) %>% 
    mutate(num_notes = n(), nom_shorter = next_ioi/beat_dur) %>% 
    ungroup() %>% 
    filter(!is.na(sr), sr>= 1, sr <= max_sr, nom_shorter >= 1/6, tatum == 1) %>% 
    group_by(id) %>%
    mutate(n_sr = n()) %>% 
    ungroup() %>% 
    filter(n_sr >= n_select)   %>% 
    select(id, avgtempo, 
           instrument, performer, 
           recordingyear, rhythmfeel, decade,
           style, tempoclass, swing_feel,
           ioi, next_ioi,
           beat_dur, sr, n_sr, 
           num_notes)
  
  tmp_sum <- tmp %>% 
    group_by(id) %>% 
    summarise(mean_sr = mean(sr, na.rm = T), 
              median_sr = median(sr, na.rm = T), 
              sd_sr = sd(sr),
              se_sr = se(sr),
              IQR_sr = IQR(sr),
              hdi_lower = HDInterval::hdi(sr)["lower"],
              hdi_upper = HDInterval::hdi(sr)["upper"],
              n_sr = mean(n_sr)) %>% 
    ungroup()
  master <- tmp %>% left_join(tmp_sum %>% select(-n_sr), by = "id")
  master_sum <- tmp_sum %>% 
    left_join(tmp %>%  select(-n_sr) %>% distinct(id, .keep_all = T), by = "id")
  #Regression models
  fit_lm   <- lm(median_sr ~ poly(avgtempo/120, 1) + swing_feel, data = master_sum) 
  fit_lm_w <- lm(median_sr ~ poly(avgtempo/120, 1) + swing_feel, weights = 1/sd_sr^2,  data = master_sum) 
  fit_lmer <- lmer(median_sr ~ poly(avgtempo, 1) + swing_feel + (1|performer), data = master_sum) 
  fit_lmer_w <- lmer(median_sr ~ poly(avgtempo, 1) + swing_feel+ (1|performer), weights = 1/sd_sr^2, data = master_sum) 
  #fit_monster <- lmer(sr ~ poly(log2(avgtempo/120), 2) + swing_feel+ (1|performer) + (1|id), data = master)
  fit_lm_l   <- lm(median_sr ~ poly(log2(avgtempo/120), 1) + swing_feel, data = master_sum) 
  fit_lm_w_l <- lm(median_sr ~ poly(log2(avgtempo/120), 1) + swing_feel, weights = 1/sd_sr^2,  data = master_sum) 
  fit_lmer_l <- lmer(median_sr ~ poly(log2(avgtempo/120), 1) + swing_feel + (1|performer), data = master_sum) 
  fit_lmer_w_l <- lmer(median_sr ~ poly(log2(avgtempo/120), 1) + swing_feel+ (1|performer), weights = 1/sd_sr^2, data = master_sum) 
  #fit_monster <- lmer(sr ~ poly(log2(avgtempo/120), 2) + swing_feel+ (1|performer) + (1|id), data = master)
  print(BIC(fit_lm, fit_lm_l, fit_lm_w, fit_lm_w_l, fit_lmer, fit_lmer_l, fit_lmer_w, fit_lmer_w_l))
  print(fit_lmer_w_l %>% broom::tidy())
  print(fit_lmer_w_l %>% broom::glance())
  #Tempo vs. Swing Ratio 
  mean_tempo <- median(master_sum$avgtempo)
  mean_sr <- mean(master_sum$mean_sr, na.rm = T)
  median_sr <- median(master_sum$mean_sr, na.rm = T)
  #q <- master_sum %>% ggplot(aes(x =  log2(avgtempo/60), y = median_sr, colour = swing_feel)) 
  q <- master_sum %>% ggplot(aes(x =  avgtempo, y = mean_sr, colour = swing_feel)) 
  q <- q + geom_point() 
  q <- q + geom_errorbar(aes(ymin = mean_sr - 1.96*se_sr, ymax = mean_sr +  1.96*se_sr), alpha = .5) 
  q <- q + geom_smooth(method = "lm", formula = "y~poly(x,2)")
  q <- q + geom_hline(yintercept = mean_sr)
  q <- q + geom_hline(yintercept = median_sr, colour = "yellow")
  q <- q + theme_bw()
  #q <- q + labs(x = "Log Mean Tempo/120 bpm ", y = "Median BUR") 
  q <- q + labs(x = "Mean Tempo (bpm) ", y = "BUR", title = sprintf("Mean/Median BUR = %.2f/%.2f", mean_sr, median_sr))
  q <- q + theme(legend.title = element_blank())
  q <- q + scale_x_continuous(breaks = c(60, 120, 180, 240, 300))
  #print(fit_lmer_w %>% broom::glance())
  #print(fit_lmer_w %>% broom::tidy())
  browser()
  assign("master_sum", master_sum, globalenv())
  assign("master", master, globalenv())
  q
}

historic_trends <- function(){
  fit1 <- lmer(log2_sr ~ recordingyear * sr_type + (1|id) + swing_feel + tempoclass, srw_f )
  fit2 <- lmer(log2_sr ~ recordingyear * sr_type + (1|id) +              tempoclass, srw_f )
  fit3 <- lmer(log2_sr ~ recordingyear * sr_type + (1|id) + swing_feel             , srw_f )
  fit4 <- lmer(log2_sr ~ recordingyear * sr_type + (1|id)                          , srw_f )
  fit5 <- lmer(log2_sr ~ recordingyear * sr_type + (1|id) +                avgtempo, srw_f )
  AIC(fit1, fit2, fit3, fit4, fit5)
  print(fit2 %>% anova())
  srw_f_sum <- srw_f_sum %>% left_join(srw_f %>% select(id, swing_feel, tempoclass), by = "id") %>% distinct(id, .keep_all = T)
  fit1 <- lm(mean_log2_sr ~ recordingyear * sr_type + swing_feel + tempoclass, srw_f_sum )
  
}
performer_analysis <- function(min_n_sr = 100, min_n_solo = 5){
 
  good_performer <- srw_f %>% 
    distinct(performer, id) %>% 
    count(performer) %>% 
    filter(n >= min_n_solo ) %>% 
    pull(performer)
  
  performer_d2 <- srw_f %>% 
    group_by(sr_type) %>% 
    mutate(mean_log2_sr = mean(log2_sr)) %>%  
    ungroup() %>% 
    group_by(performer, sr_type) %>% 
    summarise(d_m = mean(log2_sr) - mean_log2_sr[1], m = mean_log2_sr[1]) %>% 
    mutate(d_sr = 2^d_m, sr = 2^(d_m + m)) 
  
  performer_d <- performer_gmm %>% 
    group_by(performer) %>% 
    mutate(min_sr = min(sr)) %>% 
    ungroup() %>% 
    mutate(performer_label = sprintf("%s (%s)", performer, n))
  browser()
  q <- performer_d %>% filter(performer %in% good_performer) %>% 
    ggplot(aes(x = fct_reorder(performer_label, sr, max), y = sr)) 
  q <- q + geom_linerange(aes(xmin = performer_label, xmax = performer_label, ymin = min_sr, ymax = sr)) 
  q <- q + geom_point(aes(colour = factor(sr_type, 
                                          levels = c("snap", "even", "swung")), 
                          size = sqrt(prob/2)), 
                      show.legend = F) 
  q <- q + theme_bw() 
  
  q <- q + coord_flip() 
  q <- q + theme(legend.title = element_blank(), 
                 legend.position = c(.9,.75), 
                 legend.background = element_rect(colour = "black"))
  q <- q  + geom_hline(yintercept = c(.5, 1, 2), linetype = "dashed")
  q <- q + labs(x = "Performer",  y= "Mean BUR")
  
  print(q)
  ggsave("figs/bur_comp_by_performer.png", dpi = 600, height  = 8.98, width =  5.67) 
}
raster_plot <- function(data, max_ioi = 1.0, alpha = .02){
  data <- data %>% filter(ioi < max_ioi, next_ioi < max_ioi)
  data <- data %>% 
    mutate(tot_ioi_rank = rank(tot_ioi), short_first  = ioi <= next_ioi) %>% 
    filter(tot_ioi_rank > 0)
  
  q <- data %>% 
    #filter(short_first) %>% 
    ggplot()
  #q <- q + geom_point(aes(x = ioi, y = tot_ioi_rank), alpha  = alpha, color = "blue")
  #q <- q + geom_point(aes(x = tot_ioi, y = tot_ioi_rank), alpha = alpha, colour = "blue") 
  q <- q + geom_point(aes(x = tot_ioi, y = ioi/tot_ioi), alpha  = alpha, color = "blue") + geom_smooth(aes(x = tot_ioi, y = ioi/tot_ioi))
  #q <- q + geom_point(data = data %>% filter(!short_first), aes(x = -next_ioi, y = tot_ioi_rank), alpha = alpha, colour= "red") 
  #q <- q + geom_point(data = data %>% filter(!short_first), aes(x = ioi, y = tot_ioi_rank), alpha = alpha, colour = "red")
  q <- q + theme_bw()
  q
}

group_gmm_analysis <- function(data, min_n_sr = 15){
  solo_gmm_bt_15 <- fit_GMM(data, "id", min_n_sr = min_n_sr, use_bootstrap = T) %>% 
    left_join(data %>% distinct(id, style, tempoclass, recordingyear, decade, swing_feel), by = "id")
  #basic stats
  knitr::kable(
    solo_gmm %>% distinct(id, num_clust) %>% count(num_clust))
  knitr::kable(
    solo_gmm %>% filter(num_clust == 1) %>% summarise(mean_bur = 2^mean(log2_sr), 
                                                          sr_range = 2^range(log2_sr), 
                                                          s = 2^sd(log2_sr))
    )
  knitr::kable(
    solo_gmm %>% 
      filter(num_clust == 2) %>%
      group_by(cluster_id) %>% 
      summarise(mean_bur = 2^mean(log2_sr), 
                sr_range = 2^range(log2_sr), 
                s = 2^sd(log2_sr))
  )
  solos_with_swung_triples =  solo_gmm %>% distinct(id, sr_type) %>% filter(sr_type == "swung") %>% nrow()
  mean_swung_triples <- srw_f %>% count(id, sr_type) %>% filter(sr_type =="swung") %>% summarise(m = mean(n), s = sd(n)) 
  mean_even_triples <- srw_f %>% count(id, sr_type) %>% g(sr_type =="swung") %>% summarise(m = mean(n), s = sd(n)) 
  mean_triple_numbers_by_type <- srw_f %>% 
    count(id, sr_type) %>% 
    ungroup() %>% 
    group_by(sr_type) %>% 
    summarise(m = mean(n), s = sd(n), span = sprintf("%s-%s", min(n), max(n))) 
  
  knitr::kable(
    tibble(solos_with_swung_cluster = solos_with_swung_triples, 
              total_selection = solos_with_swung_triples/n_distinct(solo_gmm),
              total_wjd = solos_with_swung_triples/456)) 
  n_solos_with_swung_triples <- srw_f %>% count(id, sr_type) %>% filter(sr_type =="swung", n >= 1) %>% nrow()
  #triple types
  knitr::kable(
    srw_f %>% filter(!is.na(sr) & tatum == 1) %>% count(sr_type) %>% mutate(freq = n/sum(n)) 
  )
  solo_gmm  %>% 
    group_by(style, sr_type) %>% 
    summarise(m_sr = 2^mean(log2_sr), n = n(), solos = n_distinct(id)) %>% 
    mutate(min_n_sr = min_n_sr)
}
style_gmm_analysis <- function(data, min_n_sr = 15){
  style_gmm <- fit_GMM(srw, "style", use_bootstrap = T, min_n_sr = min_n_sr) %>% 
    select(-num_clust, -log2_sr, -log2_sr_lo, -log2_sr_hi, -cluster_id) %>%  
    mutate_if(is.numeric, round, 2) %>% 
    mutate(ci_95 = sprintf("[%s, %s]", sr_lo, sr_hi)) %>% 
    select(-sr_lo, -sr_hi) %>% 
    mutate(style = factor(style, levels = c("TRADITIONAL", "SWING", "BEBOP", "COOL", "HARDBOP", "FREE", "FUSION", "POSTBOP"))) %>% 
    arrange(style) %>% 
    select(style, sr_type, everything()) %>% 
    write.table(file = "style_gmm.csv", sep = ";", dec = ",", row.names = F)  
  
}
boots_trap_gmm <- function(data, group_var = "id", size = 10000, iteration = 10, min_n_sr = 15){
  data <- data %>% filter(!is.na(sr))
  map_dfr(min_n_sr, function(y){
    map_dfr(1:iteration, function(x){
      tmp <- data %>% group_by(!!sym(group_var)) %>% sample_n(n(), replace = T) %>% ungroup()
      #browser()
      fgmm <- 
        tryCatch({
          fit_GMM(tmp, min_n_sr = y) %>% mutate(iter = x, min_n_sr = y)},
          error = function(e){return(NULL)}  
      )
      fgmm
    })  
  })
}

get_triple_stats <- function(data = srw_f){
  mean_swung_triples <- data %>% 
    count(id, sr_type) %>% 
    group_by(sr_type) %>% 
    summarise(m = mean(n), s = sd(n))
  id_stats <- data %>% 
    group_by(sr_type) %>% 
    summarise(n= n_distinct(id), freq = n/456)
  type_stats <- data %>% 
    mutate(total = nrow(.)) %>% 
    count(sr_type, total) %>% 
    mutate(freq = n/total )
  type_stats %>% 
    left_join(mean_swung_triples, by = "sr_type") %>% 
    left_join(id_stats, by = "sr_type") %>% 
    select(-total, 
           n_triples = n.x, 
           f_triples = freq.x, 
           mean_triples = m, 
           sd_triples = s, 
           n_solos = n.y, 
           f_solos =freq.y)
}