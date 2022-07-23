rm(list = ls())

3+3

        
        weight_kg <- 55
        weight_kg
        weight_kg * 10
        weight_kg <- 100
        weight_kg

#Challenge
        
        mass <- 47.5            # mass?
        age  <- 122             # age?
        mass <- mass * 2.0      # mass?
        age  <- age - 20        # age?
        mass_index <- mass/age  # mass_index?

# now calculate the mass index
        
        mass_index <- mass/age
        mass_index
        round(3.1415926, digits = 3)

#Vectors and data types 

        weight_g <- c(50, 60, 65, 82)
        weight_g        
        animals <- c("mouse", "rat", "dog")
        months <- c("1", "2", "3")
        months
        as.character(animals)
        
        #dir.create("data_raw")
        download.file(url = "https://ndownloader.figshare.com/files/2292169",
                      destfile = "data_raw/portal_data_joined.csv")
        
#installing packages
        
        #install.packages("tidyverse")
        
#loading packages
        library(tidyverse)
        
#laoding data
        surveys <- read_csv("data_raw/portal_data_joined.csv")
        str(surveys)
        view(surveys)        
        glimpse(surveys)    
        head(surveys)

        #Q: what is the class of the object surveys? 
        class(surveys) # it is a dataframe
        # how many lows and howm any columns in this object?
        str(surveys) # 34786 rows and 13 columns
        

        
        
######Latter session (after 7:40 PM)
        
#date and time
        
        library(lubridate)
        my_date <- ymd("2015-01-01")        
        str(my_date)

#Selecting data
                
        # select only columns plot_id, species_id and weight
        
        select(surveys, plot_id, species_id, weight)
        
        #Everything but plot_id, species_id and weight
        select(surveys, -plot_id, -species_id, -weight)
        
        # just get the row wehre the eyar equals 1995
        
        filter(surveys, year == 1995)
        
        survey1995 <- filter(surveys, year == 1995)
        survey1995

        
#Pipes
        #Creating a data frame where
        #weight < 5
        #colums of species_id, sew, weight
        new_survey <- surveys %>%
                        filter(weight < 5) %>%
                        select(species_id, sex, weight)
        new_survey

#task
        #Using pipes, subset the surveys data to include animals collcted
        #before 1995
        #columns only yea, sex and weight
        surveys %>%
                filter( year < 1995) %>%
                select(year, sex, weight)
        
        
