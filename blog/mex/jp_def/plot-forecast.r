library(ggplot2)
library(tidyr)
library(dplyr)
library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)

source("utils/geom-stepribbon.r")
#---------------------------------------------------------------------------

load(paste0('results/',"base",'-',Sys.Date(),'-stanfit.Rdata'))
  
data_interventions <- read.csv("interventions.csv", 
                                 stringsAsFactors = FALSE)
  


    N <- length(dates[[1]])
    N2 <- N + 7
    country <- "MEX"
    
    predicted_cases <- colQuantiles(prediction[,1:N,1], probs=.5)
    predicted_cases_li <- colQuantiles(prediction[,1:N,1], probs=.025)
    predicted_cases_ui <- colQuantiles(prediction[,1:N,1], probs=.975)
    
    estimated_deaths <- colQuantiles(estimated.deaths[,1:N,1], probs=.5)
    estimated_deaths_li <- colQuantiles(estimated.deaths[,1:N,1], probs=.025)
    estimated_deaths_ui <- colQuantiles(estimated.deaths[,1:N,1], probs=.975)
    
    estimated_deaths_forecast <- cumsum(colQuantiles(estimated.deaths[,1:N2,1], probs=.5)[N:N2])
    estimated_deaths_li_forecast <- cumsum(colQuantiles(estimated.deaths[,1:N2,1], probs=.025)[N:N2])
    estimated_deaths_ui_forecast <- cumsum(colQuantiles(estimated.deaths[,1:N2,1], probs=.975)[N:N2])
    
    estimated_deaths_cf <- cumsum(colQuantiles(estimated.deaths.cf[,1:N2,1], probs=.5))
    estimated_deaths_li_cf <- cumsum(colQuantiles(estimated.deaths.cf[,1:N2,1], probs=.025))
    estimated_deaths_ui_cf <- cumsum(colQuantiles(estimated.deaths.cf[,1:N2,1], probs=.975))
    
    rt <- colMeans(out$Rt[,1:N,1])
    rt_li <- colQuantiles(out$Rt[,1:N,1],probs=.025)
    rt_ui <- colQuantiles(out$Rt[,1:N,1],probs=.975)
    
    data_country <- data.frame("time" = as_date(as.character(dates[[1]])),
                               "country" = rep(country, length(dates[[1]])),
                               #"country_population" = rep(country_population, length(dates[[i]])),
                               "reported_cases" = reported_cases[[1]], 
                               "reported_cases_c" = cumsum(reported_cases[[1]]), 
                               "predicted_cases_c" = cumsum(predicted_cases),
                               "predicted_min_c" = cumsum(predicted_cases_li),
                               "predicted_max_c" = cumsum(predicted_cases_ui),
                               "predicted_cases" = predicted_cases,
                               "predicted_min" = predicted_cases_li,
                               "predicted_max" = predicted_cases_ui,
                               "deaths" = deaths_by_country[[1]],
                               "deaths_c" = cumsum(deaths_by_country[[1]]),
                               "estimated_deaths_c" =  cumsum(estimated_deaths),
                               "death_min_c" = cumsum(estimated_deaths_li),
                               "death_max_c"= cumsum(estimated_deaths_ui),
                               "estimated_deaths" = estimated_deaths,
                               "death_min" = estimated_deaths_li,
                               "death_max"= estimated_deaths_ui,
                               "rt" = rt,
                               "rt_min" = rt_li,
                               "rt_max" = rt_ui)
    
    times <- as_date(as.character(dates[[1]]))
    times_forecast <- times[length(times)] + 0:7
    data_country_forecast <- data.frame("time" = times_forecast,
                                        "country" = rep(country, 8),
                                        "estimated_deaths_forecast" = estimated_deaths_forecast+head(tail(data_country$estimated_deaths_c, n=2), n=1),
                                        "death_min_forecast" = estimated_deaths_li_forecast+head(tail(data_country$death_min_c, n=2), n=1),
                                        "death_max_forecast"= estimated_deaths_ui_forecast+head(tail(data_country$death_max_c, n=2), n=1))
    
   
    times_cf <- times[1] + 0:(N2-1)
    data_country_cf <- data.frame("time" = times_cf,
                                        "country" = rep(country, (N2)),
                                        "estimated_deaths_cf" = estimated_deaths_cf,
                                        "death_min_cf" = estimated_deaths_li_cf,
                                        "death_max_cf"= estimated_deaths_ui_cf)
    

  
  data_deaths <- data_country %>%
    select(time, deaths, estimated_deaths) %>%
    gather("key" = key, "value" = value, -time)
  
  data_deaths_forecast <- data_country_forecast %>%
    select(time, estimated_deaths_forecast) %>%
    gather("key" = key, "value" = value, -time)
  
  # Force less than 1 case to zero
  data_deaths$value[data_deaths$value < 1] <- NA
  data_deaths_forecast$value[data_deaths_forecast$value < 1] <- NA
  data_deaths_all <- rbind(data_deaths, data_deaths_forecast)
  
  p <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = deaths_c), 
             fill = "coral4", stat='identity', alpha=0.5) + 
    geom_line(data = data_country, aes(x = time, y = estimated_deaths_c), 
              col = "deepskyblue4") + 
    geom_line(data = data_country_forecast, 
              aes(x = time, y = estimated_deaths_forecast), 
              col = "deepskyblue4", alpha = 0.5) + 
    geom_line(data = data_country_cf, 
              aes(x = time, y = estimated_deaths_cf), 
              col = "red", alpha = 0.5)+
    geom_ribbon(data = data_country, aes(x = time, 
                                         ymin = death_min_c, 
                                         ymax = death_max_c),
                fill="deepskyblue4", alpha=0.3) +
    geom_ribbon(data = data_country_forecast, 
                aes(x = time, 
                    ymin = death_min_forecast, 
                    ymax = death_max_forecast),
                fill = "deepskyblue4", alpha=0.35) +
    geom_ribbon(data = data_country_cf, 
                aes(x = time, 
                    ymin = death_min_cf, 
                    ymax = death_max_cf),
                fill = "red", alpha=0.35)+
    geom_vline(xintercept = data_deaths$time[length(data_deaths$time)], 
               col = "black", linetype = "dashed", alpha = 0.5) + 
    #scale_fill_manual(name = "", 
    #                 labels = c("Confirmed deaths", "Predicted deaths"),
    #                 values = c("coral4", "deepskyblue4")) + 
    xlab("Date") +
    ylab("Daily number of deaths\n") + 
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    scale_y_continuous(trans='log10', labels=comma) + 
    coord_cartesian(ylim = c(1, 100000), expand = FALSE) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill=guide_legend(ncol=1, reverse = TRUE)) + 
    annotate(geom="text", x=data_country$time[length(data_country$time)]+8, 
             y=10000, label="Forecast",
             color="black")
  print(p)
  
  ggsave(file= paste0("figures/", country, "_forecast_", filename, ".pdf"), 
         p, width = 10)

#-----------------------------------------------------------------------------------------------
make_forecast_plot()

