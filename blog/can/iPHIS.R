library(tidyverse)
library(lubridate)
library(countrycode)
library(ggrepel)
library(prismatic)
library(ggsci)
library(paletteer)
library(incidence)
library(EpiEstim)
library(projections)
library(epitrix)
library(distcrete)
library(magrittr)


countries<-c("CAN","ITA","ESP","USA","MEX","BRA","ETH","COL","MYS")

confirmed_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", col_types = cols())
deaths_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv", col_types = cols())

ontario2 <- read_csv("https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv", col_types = cols())

ontario2<-ontario2 %>% rename(date=ACCURATE_EPISODE_DATE) %>% arrange(date) %>%
  group_by(date) %>% tally(name = "new_confirmed") %>% mutate(iso3c="Ontario",data_type="iPHIS",confirmed=cumsum(new_confirmed)) %>% 
  select(iso3c,date,confirmed,new_confirmed,data_type)

confirmed_data<-confirmed_raw %>%
  select(-Lat, -Long) %>%
  rename(country = `Country/Region`,Province_Re=`Province/State`) %>%
  mutate(iso3c = countrycode(country,
                             origin = "country.name",
                             destination = "iso3c")) %>%
  select(-country) %>%
  filter(iso3c %in% countries) %>% 
  pivot_longer(
    c(-Province_Re,-iso3c), 
    names_to = "date_str", 
    values_to = "confirmed"
  ) %>%
  mutate(date = mdy(date_str)) %>%
  group_by(Province_Re)  %>% 
  mutate(new_confirmed=c(0,diff(confirmed)))  %>% ungroup() %>%
  select(iso3c,Province_Re, date,confirmed,new_confirmed) 

deaths_data<-deaths_raw %>%
  select(-Lat, -Long) %>%
  rename(country = `Country/Region`,Province_Re=`Province/State`) %>%
  mutate(iso3c = countrycode(country,
                             origin = "country.name",
                             destination = "iso3c")) %>%
  select(-country) %>%
  filter(iso3c %in% countries) %>% 
  pivot_longer(
    c(-Province_Re,-iso3c), 
    names_to = "date_str", 
    values_to = "deaths"
  ) %>%
  mutate(date = mdy(date_str)) %>%
  group_by(Province_Re)  %>% 
  mutate(new_deaths=c(0,diff(deaths)))  %>% ungroup() %>%
  select(iso3c,Province_Re, date,deaths,new_deaths) 


data_covid <- confirmed_data %>%
  full_join(deaths_data, by = c("iso3c","Province_Re", "date")) 


ontario<-data_covid %>% filter(Province_Re=="Ontario") %>%
  mutate(iso3c="Ontario",data_type="JH") %>% select(iso3c,date,confirmed,new_confirmed,data_type)


plot_project<-function(country="Ontario",conf=20,method=c("Estimate_R","Simulate_R"),
                       pred_days=21,n_sim=1000,mu_R=2.5,sigma_R=0.5,cumul=F,
                       config_si=list(mean_si = 7.5, std_mean_si = 2, 
                                      min_mean_si = 1, max_mean_si = 8.4, std_si = 3.4, std_std_si = 1, 
                                      min_std_si = 0.5, max_std_si = 4, n1 = 1000, n2 = 1000)) {
  
  df<-data_covid %>% bind_rows(ontario) %>% bind_rows(ontario2) %>%
    group_by(iso3c,date,data_type) %>%
    summarise_at(vars(-Province_Re),funs(sum)) %>% filter(confirmed>=conf,iso3c==country) %>%
    rename(I=new_confirmed,dates=date) 
  
  
  
  max_cumul<-df %>% filter(dates=="2020-03-23") %>% summarise(max_cumul=max(confirmed))
  max_cumul<-as.numeric(max_cumul[3])
  
  incidence_data<- df %>% filter(data_type=="iPHIS",dates<="2020-03-23") %>% uncount(I)
  
  incidence_object <- incidence(dates = incidence_data$dates) 
  param <- gamma_mucv2shapescale(config_si$mean_si, config_si$std_si / config_si$mean_si)
  w <- distcrete("gamma", interval = 1,
                 shape = param$shape,
                 scale = param$scale, w = 0.5)
  date_range <- 1:(length(get_dates(incidence_object)))
  
  df2<-df %>% filter(data_type=="iPHIS",dates<="2020-03-23") %>% select(-deaths,-new_deaths)
  

  R_sample <- sample_posterior_R(re_est, n = 1000,window = 1:nrow(re_est$R))

  re_est<- estimate_R(df2, method = "uncertain_si",
                      config = make_config(config_si))
  
  R_sample <- sample_posterior_R(re_est, n = 1000,window = nrow(re_est$R))
  pred_growth_var <- project(incidence_object[date_range],
                             R = list(R_sample, R_sample_Korea_3),
                             R_fix_within = F,
                             si = w,
                             n_days = 60, time_change = c(15),
                             n_sim = n_sim)

  if (cumul ==T) {
    pred_growth_var<-cumulate(pred_growth_var)
    df<-df %>% mutate(I=confirmed)
  }
  
  ### functions
  
  quantile_pal <- grDevices::colorRampPalette(
    c("#b3c6ff", "#d147a3", "#993366"), bias = 2)
  
  color_quantiles <- function(x, palette = quantile_pal) {
    labels <- as.character(unique(x))
    dist_from_median <- 50 - abs(50-as.numeric(sub("%", "", labels)))
    out <- palette(51)[dist_from_median + 1]
    names(out) <- labels
    out
  }
  
  
  
  transp <- function(col, alpha = .5){
    res <- apply(grDevices::col2rgb(col), 2,
                 function(c)
                   grDevices::rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
  }
  
  ###################
  
  dates <- attr(pred_growth_var, "dates")
  stats <- t(apply(pred_growth_var, 1, stats::quantile, c(0.01, .99)))
  if (cumul ==T) {
    stats<-stats+max_cumul
  }
  
  temp <- cbind.data.frame(dates, stats)
  names(temp) <- c("dates", "ymin", "ymax")
  
  df<- df %>% filter(data_type=="iPHIS")
  
  ribbon_color <- color_quantiles(c(0.01, .99), quantile_pal)[1]
  
  ribbon_color <- transp(ribbon_color, 0.2)
  out<-ggplot()+
    geom_ribbon(
      data = temp,
      aes(x = dates, ymin = ymin, ymax = ymax),
      fill = ribbon_color)+geom_point(aes(x=dates,y=I),data = df,size=0.9,color="darkblue")+geom_line(aes(x=dates,y=I),color="darkblue",
                                                                                                      size=.6,data = df)+
    scale_y_continuous(labels = scales::comma_format(accuracy = 1),trans = "log2",
                       n.breaks = 10)
  stats <- t(apply(pred_growth_var, 1, stats::quantile, c(.01,0.25,.5,.75,.99)))
  if (cumul ==T) {
    stats<-stats+max_cumul
  }
  quantiles <- rep(colnames(stats), each = nrow(stats))
  quantiles <- factor(quantiles, levels = unique(quantiles))
  temp <- cbind.data.frame(dates = rep(dates, ncol(stats)),
                           quantile = quantiles,
                           value = as.vector(stats),
                           stringsAsFactors = FALSE)
  
  
  
  colors <- color_quantiles(temp$quantile, quantile_pal)
  colors <- transp(colors, 1)
  ylab <- ifelse(cumul==T,
                 "Cumulative incidence",
                 "Daily incidence")
  titlab <- ifelse(cumul==T,
                   "Cumulative incidence",
                   "Daily incidence")
  
  out <- suppressMessages(
    out +
      ggplot2::geom_line(
        data = temp,
        aes(x = dates, y = value, color = quantile),
        linetype = 1,
        size = 1
      ) + scale_color_manual(values = colors))+
    labs(title= paste(titlab,"of observed and predicted confirmed cases, ",df$iso3c),
         x = "Date", color="Quantile",
         y = paste(ylab,"(log2 scale)"),
         caption = "iPHIS data",
         subtitle = paste("Data as of", format(max(data_covid$date), "%A, %B %e, %Y")))+
    scale_x_date(date_breaks = "1 week")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),plot.title = element_text(size = rel(2), face = "bold"),
          plot.subtitle = element_text(size = rel(1.5)),
          axis.text.y = element_text(size = rel(2)),
          axis.title.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.5)),
          panel.grid.major = element_line(colour = "grey90", size = .1),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey50", size = .9))
  return<-list(re_est=re_est,out=out,temp=temp)
  return(return)
}

plot_Ri <- function(estimate_R_obj) {
  p_I <- plot(estimate_R_obj, "incid") # plots the incidence
  p_SI <- plot(estimate_R_obj, "SI") # plots the serial interval distribution
  p_Ri <- plot(estimate_R_obj, "R")
  return(gridExtra::grid.arrange(p_I, p_SI, p_Ri, ncol = 1))
}

result<-plot_project(method = "Estimate_R",cumul = F)

