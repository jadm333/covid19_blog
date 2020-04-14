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

  
load(paste0('results/',"base",'-',Sys.Date(),'.Rdata'))
  
data_interventions <- read.csv("interventions.csv", 
                                 stringsAsFactors = FALSE)
covariates <- data_interventions[1, c(1,2,3,4,5,6, 7, 8)]
  
N <- length(dates[[1]])
country <- "MEX"
    
predicted_cases <- colQuantiles(prediction[,1:N,1], probs=.5)
predicted_cases_li <- colQuantiles(prediction[,1:N,1], probs=.025)
predicted_cases_ui <- colQuantiles(prediction[,1:N,1], probs=.975)
predicted_cases_li2 <- colQuantiles(prediction[,1:N,1], probs=.25)
predicted_cases_ui2 <- colQuantiles(prediction[,1:N,1], probs=.75)
    
    
estimated_deaths <- colQuantiles(estimated.deaths[,1:N,1], probs=.5)
estimated_deaths_li <- colQuantiles(estimated.deaths[,1:N,1], probs=.025)
estimated_deaths_ui <- colQuantiles(estimated.deaths[,1:N,1], probs=.975)
estimated_deaths_li2 <- colQuantiles(estimated.deaths[,1:N,1], probs=.25)
estimated_deaths_ui2 <- colQuantiles(estimated.deaths[,1:N,1], probs=.75)
    
rt <- colQuantiles(out$Rt[,1:N,1],probs=.5)
rt_li <- colQuantiles(out$Rt[,1:N,1],probs=.025)
rt_ui <- colQuantiles(out$Rt[,1:N,1],probs=.975)
rt_li2 <- colQuantiles(out$Rt[,1:N,1],probs=.25)
rt_ui2 <- colQuantiles(out$Rt[,1:N,1],probs=.75)
    

covariates_country <- covariates[which(covariates$ï..Country == country), 2:8]   
    
covariates_country$lockdown = NULL
covariates_country$sport = NULL 
covariates_country$travel_restrictions = NULL 
covariates_country_long <- gather(covariates_country[], key = "key", 
                                      value = "value")
covariates_country_long$x <- rep(NULL, length(covariates_country_long$key))
un_dates <- unique(covariates_country_long$value)
    
    for (k in 1:length(un_dates)){
      idxs <- which(covariates_country_long$value == un_dates[k])
      max_val <- round(max(rt_ui)) + 0.3
      for (j in idxs){
        covariates_country_long$x[j] <- max_val
        max_val <- max_val - 0.3
      }
    }
    
    
    covariates_country_long$value <- mdy(covariates_country_long$value) 
    covariates_country_long$country <- rep(country, 
                                           length(covariates_country_long$value))
    
    data_country <- data.frame("time" = as_date(as.character(dates[[1]])),
                               "country" = rep(country, length(dates[[1]])),
                               "reported_cases" = reported_cases[[1]], 
                               "reported_cases_c" = cumsum(reported_cases[[1]]), 
                               "predicted_cases_c" = cumsum(predicted_cases),
                               "predicted_min_c" = cumsum(predicted_cases_li),
                               "predicted_max_c" = cumsum(predicted_cases_ui),
                               "predicted_min_c2" = cumsum(predicted_cases_li2),
                               "predicted_max_c2" = cumsum(predicted_cases_ui2),
                               "predicted_cases" = predicted_cases,
                               "predicted_min" = predicted_cases_li,
                               "predicted_max" = predicted_cases_ui,
                               "predicted_min2" = predicted_cases_li2,
                               "predicted_max2" = predicted_cases_ui2,
                               "deaths" = deaths_by_country[[1]],
                               "deaths_c" = cumsum(deaths_by_country[[1]]),
                               "estimated_deaths_c" =  cumsum(estimated_deaths),
                               "death_min_c" = cumsum(estimated_deaths_li),
                               "death_max_c"= cumsum(estimated_deaths_ui),
                               "death_min_c2" = cumsum(estimated_deaths_li2),
                               "death_max_c2"= cumsum(estimated_deaths_ui2),
                               "estimated_deaths" = estimated_deaths,
                               "death_min" = estimated_deaths_li,
                               "death_max"= estimated_deaths_ui,
                               "death_min2" = estimated_deaths_li2,
                               "death_max2"= estimated_deaths_ui2,
                               "rt" = rt,
                               "rt_min" = rt_li,
                               "rt_max" = rt_ui,
                               "rt_min2" = rt_li2,
                               "rt_max2" = rt_ui2)


  data_cases_95 <- data.frame(data_country$time, data_country$predicted_min_c, 
                              data_country$predicted_max_c)
  names(data_cases_95) <- c("time", "cases_min", "cases_max")
  data_cases_95$key <- rep("nintyfive", length(data_cases_95$time))
  data_cases_50 <- data.frame(data_country$time, data_country$predicted_min_c2, 
                              data_country$predicted_max_c2)
  names(data_cases_50) <- c("time", "cases_min", "cases_max")
  data_cases_50$key <- rep("fifty", length(data_cases_50$time))
  data_cases <- rbind(data_cases_95, data_cases_50)
  levels(data_cases$key) <- c("ninetyfive", "fifty")
  
  p1 <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = reported_cases_c), 
             fill = "coral4", stat='identity', alpha=0.5) + 
    geom_ribbon(data = data_cases_95, 
                aes(x = time, ymin = cases_min, ymax = cases_max, fill = key)) +
    geom_line(aes(x=time,y=predicted_cases_c,color="Predicción"))+
    xlab("") +
    ylab("Número acumulado de casos confirmados") +
    scale_y_continuous(labels = scales::comma_format(accuracy = 1),n.breaks = 10)+
    scale_x_date(date_breaks = "3 days", labels = date_format("%e %b")) + 
    scale_fill_manual(name = "", labels = c("95% CrI"),
                      values = c(alpha("deepskyblue4", 0.35))) + 
    scale_color_manual(name = "",
                      values ="black") +
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill=guide_legend(ncol=1))
  
  data_deaths_95 <- data.frame(data_country$time, data_country$death_min_c, 
                               data_country$death_max_c)
  names(data_deaths_95) <- c("time", "death_min", "death_max")
  data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
  data_deaths_50 <- data.frame(data_country$time, data_country$death_min_c2, 
                               data_country$death_max_c2)
  names(data_deaths_50) <- c("time", "death_min", "death_max")
  data_deaths_50$key <- rep("fifty", length(data_deaths_50$time))
  data_deaths <- rbind(data_deaths_95, data_deaths_50)
  levels(data_deaths$key) <- c("ninetyfive", "fifty")
  
  
  p2 <-   ggplot(data_country, aes(x = time)) +
    geom_bar(data = data_country, aes(y = deaths_c, fill = "reported"),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_line(aes(y=estimated_deaths_c,color="Predicción"))+
    geom_ribbon(
      data = data_deaths_95,
      aes(ymin = death_min, ymax = death_max, fill = key)) +
    scale_x_date(date_breaks = "3 days", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("95% CrI"),
                      values = c(alpha("deepskyblue4", 0.35))) + 
    scale_color_manual(name = "",
                       values ="black")+
    scale_y_continuous(labels = scales::comma_format(accuracy = 1),
                       n.breaks = 10)+
    xlab("")+
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill=guide_legend(ncol=1))+
    ylab("Número acumulado de muertes")
  
  
  plot_labels <- c("Cancelación de eventos",
                   "Cierre de escuelas",
                   "Aislamiento después de presentar síntomas",
                   "Distanciamiento Social")
  
  # Plotting interventions
  data_rt_95 <- data.frame(data_country$time, 
                           data_country$rt_min, data_country$rt_max)
  names(data_rt_95) <- c("time", "rt_min", "rt_max")
  data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
  data_rt_50 <- data.frame(data_country$time, data_country$rt_min2, 
                           data_country$rt_max2)
  names(data_rt_50) <- c("time", "rt_min", "rt_max")
  data_rt_50$key <- rep("fifty", length(data_rt_50$time))
  data_rt <- rbind(data_rt_95, data_rt_50)
  levels(data_rt$key) <- c("ninetyfive", "fifth")
  
  p3 <- ggplot(data_country) +
    geom_stepribbon(data = data_rt_95, aes(x = time, ymin = rt_min, ymax = rt_max, 
                                        group = key,
                                        fill = key)) +
    geom_hline(yintercept = 1, color = 'black', size = 0.1) + 
    geom_step(aes(x=time,y=rt,linetype="Predicción"))+
    scale_linetype_manual(name = "",
                       values ="solid")+
    geom_segment(data = covariates_country_long,
                 aes(x = value, y = 0, xend = value, yend = max(x)), 
                 linetype = "dashed", colour = "black", alpha = 0.75) +
    geom_point(data = covariates_country_long, aes(x = value, 
                                                   y = x, 
                                                   group = key, 
                                                   shape = key, 
                                                   col = key), size = 2) +
    xlab("") +
    scale_y_continuous(n.breaks = 10)+
    ylab("Número efectivo de reproducción") +
    scale_fill_manual(name = "", labels = c("95% CrI"),
                      values = c(alpha("deepskyblue4", 0.55))) + 
    scale_shape_manual(name = "Intervenciones", labels = plot_labels,
                       values = c(21, 22, 23, 24, 25, 12)) + 
    scale_colour_discrete(name = "Intervenciones", labels = plot_labels) + 
    scale_x_date(date_breaks = "3 days", labels = date_format("%e %b"), 
                 limits = c(data_country$time[1], 
                            data_country$time[length(data_country$time)])) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="right")
  
  p <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 2))
