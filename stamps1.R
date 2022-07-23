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
        
        dir.create("data_raw")
        download.file(url = "https://ndownloader.figshare.com/files/2292169",
                      destfile = "data_raw/portal_data_joined.csv")
        list.files()
        data <- read.csv("data_raw/portal_data_joined.csv")
        
#installing packages
        
        install.packages()