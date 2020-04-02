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

setwd("~/tool")
countries<-c("KOR","ITA","MEX")

date_cut<-"03-30-2020"
date_cut<-mdy(date_cut)


confirmed_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", col_types = cols())
deaths_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv", col_types = cols())

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



  
  df<-data_covid %>% 
    group_by(iso3c,date) %>%
    summarise_at(vars(-Province_Re),funs(sum)) %>% filter(confirmed>=50,iso3c=="KOR",date<=date_cut) %>%
    mutate(date = ymd(date),data_type="Observed") %>%
    select(date, new_confirmed,iso3c,data_type,confirmed) %>% rename(I=new_confirmed,dates=date) 
  
  
  incidence_data<- df %>% uncount(I)
  mean_si = 7.5
  std_si = 3.4
  incidence_object <- incidence(dates = incidence_data$dates) 
  param <- gamma_mucv2shapescale(mean_si, std_si / mean_si)
  w <- distcrete("gamma", interval = 1,
                 shape = param$shape,
                 scale = param$scale, w = 0.5)
  date_range <- 1:(length(get_dates(incidence_object)))
    t_start <- c(2,11,20) # starting at 2 as conditional on the past observations
    t_end <- c(10,19,40)
    re_est<- estimate_R(df, method = "uncertain_si",
                        config = make_config(list(t_start=t_start,t_end=t_end,mean_si = 3.96, std_mean_si = 2, 
                                                  min_mean_si = 1, max_mean_si = 8.4, std_si = 4.75, std_std_si = 1, 
                                                  min_std_si = 0.5, max_std_si = 6, n1 = 1000, n2 = 1000)))
    
    str(tt$data)
    plot(re_est, "R") + geom_hline(aes(yintercept = 1), color = "red", lty = 2)+
      scale_x_date(date_breaks = "1 week",date_labels = "%b %d")+
      labs(title= expression(Effective~reproductive~number~R[e]),
           x = "Date",
           y = expression(R[0]),
           subtitle = "Based on daily incidence from South Korea")+scale_y_continuous(n.breaks = 7)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),plot.title = element_text(size = rel(2), face = "bold"),
                                                                       plot.subtitle = element_text(size = rel(1.5)),
                                                                       axis.text.y = element_text(size = rel(2)),
                                                                       axis.title.x = element_text(size = rel(1.5)),
                                                                       axis.title.y = element_text(size = rel(1.5)),
                                                                       panel.grid.major = element_line(colour = "grey90", size = .1),
                                                                       panel.background = element_blank(),
                                                                       axis.line = element_line(colour = "grey50", size = .9))
    
    plot(re_est,"incid")+geom_vline(aes(xintercept = as.numeric(df[t_end[1],1]-1)), color = "red", lty = 2)+
      geom_vline(aes(xintercept = as.numeric(df[t_end[2]-1,1])), color = "red", lty = 2)+
    scale_x_date(date_breaks = "1 week",date_labels = "%b %d")+
      labs(title= "Daily incidence of confirmed cases",
           x = "Date",
           y = "Daily incidence",
           subtitle = "South Korea (Data as of March 30)")+scale_y_continuous(n.breaks = 7)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),plot.title = element_text(size = rel(2), face = "bold"),
            plot.subtitle = element_text(size = rel(1.5)),
            axis.text.y = element_text(size = rel(2)),
            axis.title.x = element_text(size = rel(1.5)),
            axis.title.y = element_text(size = rel(1.5)),
            panel.grid.major = element_line(colour = "grey90", size = .1),
            panel.background = element_blank(),
            axis.line = element_line(colour = "grey50", size = .9))
    
     R_sample_Korea_1 <-sample_posterior_R(re_est, n = 1000,window = 1)
    hist(R_sample_Korea_1)
    R_sample_Korea_2 <-sample_posterior_R(re_est, n = 1000,window = 2)
    hist(R_sample_Korea_2)
    R_sample_Korea_3 <-sample_posterior_R(re_est, n = 1000,window = 3)
    hist(R_sample_Korea_3)
    
    
    
    
    
    

    
    
    
    
    
    
    
    

 