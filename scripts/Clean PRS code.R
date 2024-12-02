#upload data
combined_data <- read.csv("combined_data.csv")
shark_data <- read.csv("shark_data.csv")
#load packages
library(dplyr)
library(readr)
library(rnaturalearth)
library(geosphere)
library(rnaturalearthdata)
library(sf)
library(lme4) 
library(lmerTest) 
library(DHARMa) 
library(rfishbase) 
library(ggOceanMaps)
library(ggplot2)
library(MASS)
library(emmeans)
library(car)
library(sjPlot)
library(ggpubr)
library(broom)
library(ggeffects)
library(stringr)
library(patchwork)


#Section 1 -------------calculate slope values and add column---------------------

#change value column name
colnames(combined_data)<- c("year","z_value","species","unit","ocean","lat","long","order","superorder","ID")
colnames(shark_data)<- c("year","z_value","species","unit","ocean","lat","long","order","superorder","ID","slope","length","vulnerability")

#basic test plot looking at standardization

plot(data=combined_data, z_value~year)

#mean z_value each year use dplyr to group by

mean_z_value_all <- combined_data %>% group_by(year) %>% summarise(mean_z=mean(z_value, na.rm = TRUE))

#plot mean z vs year

plot(data=mean_z_value_all, mean_z~year)

model <- lm(mean_z ~ year, data = mean_z_value_all)

# Add the line of best fit
abline(model, col = "red")
summary(model)
#shows sig decrease in ocean life over the last century

#validate species names
species_data <- as.data.frame(combined_data$species) #obtain list of species in filtered dataframe
species_data <- unique(species_data) #remove duplicates

species_data$species <- validate_names(species_data$`combined_data$species`)
species_input <- as.data.frame(unique(species_data$species)) #removes repeated records

#trend column - regression slope
combined_data <- combined_data %>%
  group_by(ID) %>%
  summarise(slope = coef(lm(z_value ~ year, data = pick(cur_data())))[2], .groups = "drop")

combined_data1 <- combined_data %>%
  group_by(ID) %>%
  summarise(
    slope = coef(lm(z_value ~ year))[2],  
    .groups = "drop")
shark_data <- left_join(combined_data, combined_data1, by = "ID")

#write csv to save og shark data
write.csv(shark_data,"shark_data.csv",row.names = FALSE)

#add a time series length column
year_range <- shark_data %>%
  group_by(ID) %>%
  summarise(min_year = min(year), max_year = max(year))
year_range <- year_range %>% group_by(ID) %>% summarise(max_year-min_year)
colnames(year_range) <- c("ID","tsrange")

shark_data <- left_join(shark_data, year_range, by = "ID")

#Section 2 ---------add fishbase and IUCN data columns and make linear models-------

#find col names of fish base data for traits to use
sample_data <- species(species_input$'unique(species_data$species)')
colnames(sample_data)

estimate_data <- estimate(species_input$'unique(species_data$species)')
colnames(estimate_data)


#adding fishbase data
length_data <- species(species_input$`unique(species_data$species)`, fields = c("Species","Length"))#extracts data
colnames(length_data)<- c("species","length")
shark_data <- left_join(shark_data,length_data,by="species")

troph_data <- estimate(species_input$`unique(species_data$species)`, fields = c("Species","Troph"))
colnames(troph_data)<- c("species","troph")
shark_data <- left_join(shark_data,troph_data,by="species")

vuln_data <- species(species_input$`unique(species_data$species)`, fields = c("Species","Vulnerability"))
colnames(vuln_data)<- c("species","vulnerability")
shark_data <- left_join(shark_data,vuln_data,by="species")

depthmax_data <- estimate(species_input$`unique(species_data$species)`, fields = c("Species","DepthMax"))
colnames(depthmax_data)<- c("species","depthmax")
shark_data <- left_join(shark_data,depthmax_data,by="species")

feedingpath_data <- estimate(species_input$`unique(species_data$species)`, fields = c("Species","FeedingPath"))
colnames(feedingpath_data)<- c("species","feedingpath")
shark_data <- left_join(shark_data,feedingpath_data,by="species")


#attaching IUCN data
IUCN_csv <- read.csv("assessments.csv", header = TRUE, sep = ",")
IUCN_csv <- IUCN_csv[, c(3, 4, 14)]
shark_data <- left_join(shark_data,IUCN_csv, by="species")

#make models of different traits using lme4
length_model <- lmer(slope~length + (1 |species), data=shark_data)
summary(length_model)

ocean_model <- lmer(slope~ocean + (1 |species), data=shark_data)
summary(ocean_model) #all sig diff from atlantiv

depthmaxlength_model <- lmer(slope~depthmax*length + (1 |species), data=shark_data)
summary(depthmax_model)
summary(depthmaxlength_model)

#then use DHARMa to test assumptions
length_modelx <- simulateResiduals(length_model)
plot(length_modelx) #don't fit

ocean_modelx <- simulateResiduals(ocean_model)
plot(ocean_modelx) #don't fit

#assumptions not met - try transformations of response to fit lmer
#transformation
shark_data$log_slope <- log(shark_data$slope+1)

#make a binomial slope column
shark_data$binomial_slope <- ifelse(shark_data$slope>0, 1, 0)

#rerun lmer with transformations
length_model_log <- lmer(log_slope~length + (1 |species), data=shark_data)
summary(length_model_log)

#use DHARMa again to test transformed residuals
length_model_logx <- simulateResiduals(length_model_log)
plot(length_model_logx)

#Section 3 ----------try to use general linear models-------------

#if residuals don't fit - need to use general linear model
g_length_model <- glmer(slope~length + (1 |species), family=gaussian("identity"), data=shark_data)

#DHARMa
g_length_modelx <- simulateResiduals(g_length_model)
plot(g_length_modelx)


#Section 4 -------try binomial glmer tests using binomial slope------------

#glmer of binomial slope
#for length (significant)
glength_model <- glmer(binomial_slope~length+(1|species), data=shark_data, family = binomial("logit") )
summary(glength_model) #more likely to be increasing in abundance if you have a smaller body size
glength_modelx <- simulateResiduals(glength_model)
plot(glength_modelx)  #much better graph

#for vulnerability  (significant)
gvuln_model <- glmer(binomial_slope~vulnerability+(1|species), data=shark_data, family = binomial("logit"))
summary(gvuln_model)
gvuln_modelx <- simulateResiduals(gvuln_model)
plot(gvuln_modelx)

#for IUCN status (requires post hoc)
gIUCN_model <- glmer(binomial_slope~redlist_cat+(1|species), data=shark_data, family = binomial("logit"))
summary(gIUCN_model)  #strange - might not be useable
gIUCN_modelx <- simulateResiduals(gIUCN_model)
plot(gIUCN_modelx)

#for IUCN population trend
gtrend_model <- glmer(binomial_slope~populationTrend+(1|species), data=shark_data, family = binomial("logit"))
summary(gtrend_model)
#for feeding path (not significant)
gfeeding_model <- glmer(binomial_slope~feedingpath+(1|species), data=shark_data, family = binomial("logit"))
summary(gfeeding_model) #not significant

#for ocean (requires post hoc)
gocean_model <- glmer(binomial_slope~ocean+(1|species), data=shark_data, family = binomial("logit"))
summary(gocean_model)

#for superorder
gsuperorder_model <- glmer(binomial_slope~superorder+(1|species), data=shark_data, family = binomial("logit"))
summary(gsuperorder_model)

#for depthmax and length
glengthdepth_model <- glmer(binomial_slope~depthmax*length+(1|species), data=shark_data, family = binomial("logit") )
summary(glengthdepth_model) 
glengthdepth_modelx <- simulateResiduals(glengthdepth_model)
plot(glength_modelx) 

#trying glmer with 2 fixed variables in the model
gsuper_ocean_model <- glmer(binomial_slope~superorder+ocean+(1|species), data=shark_data, family = binomial("logit"))
summary(gsuper_ocean_model)
gsuper_ocean_model_results <- emmeans(gsuper_ocean_model, pairwise ~ superorder+ocean, adjust = "BH")
Anova(gsuper_ocean_model) #both vairiables significant
gsuper_ocean_model_results <- as.data.frame(gsuper_ocean_model_results$contrasts)

#posthoc of multi-category glmers with emmeans

#for IUCN
emmeans(gIUCN_model, pairwise ~ redlist_cat, adjust = "BH")

#for oceans
emmeans(gocean_model, pairwise ~ ocean, adjust = "BH")

#for superorder (all sig diff)
emmeans(gsuperorder_model, pairwise ~ superorder, adjust = "BH")

#Section 5 ------------visualise data with ggplot-------------

#remove rows with na values
shark_data_cut <- shark_data %>% filter(!is.na(redlist_cat)&(!is.na(binomial_slope)))
shark_data_cut$predictedbislope <- predict(gIUCN_model, type="response")

#use factor to reorder the boxplot
shark_data_cut$redlist_cat <- factor(shark_data_cut$redlist_cat,
                                     levels = c("Data Deficient", "Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered"))

#plot the predicted probabilities by redlist category
ggplot(shark_data_cut, aes(x = redlist_cat, y = predictedbislope)) +
  geom_boxplot() +
  labs(title = "",
       x = "Redlist Category",
       y = "Predicted Category (Binomial Slope)") +
  theme_minimal()

#remove rows with na values again for superorder
shark_data_cut <- shark_data %>% filter(!is.na(superorder)&(!is.na(binomial_slope)))
shark_data_cut$predictedbislope <- predict(gsuperorder_model, type="response")

#plot superorder (above 0.5 = numbers are increase, below = decreasing, over the whole time series)
ggplot(shark_data_cut, aes(x = superorder, y = predictedbislope)) +
  geom_boxplot() +
  labs(title = "",
       x = "Superorder",
       y = "Predicted Probability (Binomial Slope)") +
  theme_minimal()

#repeat for ocean
shark_data_cut <- shark_data %>% filter(!is.na(ocean)&(!is.na(binomial_slope)))
shark_data_cut$predictedbislope <- predict(gocean_model, type="response")

#plot
ggplot(shark_data_cut, aes(x = ocean, y = predictedbislope)) +
  geom_boxplot() +
  labs(title = "",
       x = "Ocean",
       y = "Predicted Probability (Binomial Slope)") +
  theme_minimal()

#Section 6 ---------make plots using sjplot-------------

#sjplot example with gg plot and ggpubr
plot <- plot_model(glength_model, type = "pred", terms = c("length"))
plot #50% - bigger (less than 250ish) = increasing in abundance, less than 50% (bigger than 250) = decresing in abundance.

#make new data frame
new <- data.frame(binomial_slope)
matrix_data <- as.matrix(new)
binomial_slope <- c(25,50,75)

#use pred function to find exact value for 50%
v <- c(226,227,228,229,230)
predict_response(glength_model, terms = "length [v]") #cutoff is 229

#test for interactions between groups
glengthsuper_model <- glmer(binomial_slope~length*superorder+(1|species), data=shark_data, family = binomial("logit") )
summary(glengthsuper_model)#no interaction - all doing the same thing except sharks are slightly more vulnerable at the same size
Anova(glength_model)
plot3 <- plot_model(glengthsuper_model, type = "pred", terms = c("length [all]", "superorder"))
plot3 <- plot3 + theme_classic()
plot3

#Section 7 ----using sjplot to plot length models for superorder subsets------

#sjplot glength_model with all superorders
glengthplot<- plot_model(glength_model, type= "pred", terms = c("length[all]"), color= "black")+  
  labs(
    title = "",
    x = "",
    y = "Probabiltiy of Abundance Trend"
  )+ 
  theme_classic()
glengthplot

#make superorder subsets
shark_subset <- subset(shark_data, superorder == "Selachimorpha")
ray_subset <-subset(shark_data, superorder == "Batoidea")
chimera_subset <- subset(shark_data, superorder == "Holocephalimorpha")

#make binomial glm's for all three superorders testing length of species as a predictor
glength_shark_model <- glmer(binomial_slope~length+(1|species), data=shark_subset, family = binomial("logit") )
summary(glength_shark_model)
glength_ray_model <- glmer(binomial_slope~length+(1|species), data=ray_subset, family = binomial("logit") )
summary(glength_ray_model)
glength_chimera_model <- glmer(binomial_slope~length+(1|species), data=chimera_subset, family = binomial("logit") )
summary(glength_chimera_model)

#sj plot the models
glength_shark_plot<- plot_model(glength_shark_model, type= "pred", terms = c("length[all]"), color = "darkblue") +
  labs(
    title = "",
    x = "",
    y = "Probabiltiy of Abundance Trend"
  )+ 
  theme_classic()
glength_ray_plot<- plot_model(glength_ray_model, type= "pred", terms = c("length[all]"), color = "red")+
  labs(
    title = "",
    x = "",
    y = "Probabiltiy of Abundance Trend"
  )+ 
  theme_classic()
glength_chimera_plot<- plot_model(glength_chimera_model, type= "pred", terms = c("length[all]"), color= "darkgreen")+
  labs(
    title = "Chimeras",
    x = "Length (cm)",
    y = "Predicted Probability (binomial slope)"
  )+ 
  theme_classic()
glength_shark_plot
glength_ray_plot
glength_chimera_plot


#remove chimeras because of insufficient data and make a histogram to show the distribution of the legnth data to replace it
#dist of length
lengthdist_plot <- ggplot(shark_data, aes(x=length))+ 
  geom_histogram(aes(y=..density..),binwidth = 10, center = 5, colour= "black",fill = "white") +
  geom_vline(aes(xintercept = mean(length)), color= "red", linetype= "dashed", size=1)+
  theme_classic()+
  labs(x="", y="Proportion of Data")
lengthdist_plot

#make a figure of these 4 graphs
#use ggpubr package to arrange and add labels
figure1 <- ggarrange(glengthplot, glength_shark_plot, lengthdist_plot, glength_ray_plot, labels= c("A", "B","","C"), ncol = 2, nrow = 2)
figure1 <- annotate_figure(
  figure1,
  bottom = text_grob("Maximum Total Length (cm)", size = 14)
)
figure1 

#Section 8 -----test interaction between length and other variables

#length and tsrange interaction
glengthtsrange_model <- glmer(binomial_slope~length*tsrange+(1|species), data=shark_data, family = binomial("logit") )
summary(glengthtsrange_model)
Anova(glengthtsrange_model) #each significant individually but not their interaction

glengthtsrange_plot<- plot_model(glengthtsrange_model, type= "pred", terms = c("length[all], tsrange"), color = "darkblue") +
  labs(
    title = "",
    x = "Length / cm",
    y = "Predicted Probability (binomial slope)"
  )+ 
  theme_classic()
glengthtsrange_plot

#length and ocean
control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
glengthocean_model <- glmer(binomial_slope~length*ocean+(1|species), data=shark_data, family = binomial("logit"), control = control )
summary(glengthocean_model)
Anova(glengthocean_model)

#Section 9 ------create a barplot showing frequency of timeseries in dataset by order--------

#summarise dataframe - complete to provide a count of species and timeseries, then add the order information
order_ts_data <- shark_data %>% 
  group_by(species, ID, order, superorder) %>% 
  summarise(n = n())


order_ts_data2 <- order_ts_data %>% 
  group_by(order, superorder) %>% 
  summarise(n = n())


ordertsbarplot <- ggplot(data = order_ts_data2, aes(x = reorder(order, -n), y=n))+ 
  geom_bar(stat = "identity", width=0.5,fill='red4') +
  coord_flip() + theme_classic() +
  xlab("Taxonomic Order")+ylab("Number of population abundance timeseries")

ordertsbarplot

#clean order_ts_data to remove na values
shark_data <- shark_data %>%
  mutate(order = str_trim(order))

order_ts_data2 <- order_ts_data %>% 
  group_by(order, superorder) %>% 
  summarise()
#make acorrection
order_ts_data2$order[order_ts_data2$order == "Echinorhiniformes"] <- "Squaliformes"

#create a horizontal bar chart
order_ts_plot <- ggplot(data = order_ts_data2, aes(x = reorder(order, n), y = n, fill = superorder)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # This makes the bar chart horizontal
  scale_fill_brewer(palette = "Set2") +  # Assign colors to each superorder
  labs(x = "Order", y = "Number of Time Series", title = "") +  #Time Series Distribution by Order
  theme_classic()+
  scale_y_continuous(limits=c(0, 300), breaks = seq(from=0, to=300, by=50), expand = c(0, 0))+
  labs(fill="Superorder") +
  theme(
    legend.position = c(0.65, 0.3), # Set normalized x, y coordinates for legend position
    legend.justification = c("left", "center"), # Align the legend anchor
    legend.text = element_text(size = 10),  # Adjust the size of the legend text
    legend.title = element_text(size = 11),  # Adjust the size of the legend title
    legend.key.size = unit(1, "cm")
  )
order_ts_plot 

#add order silhouettes after from phylopic

#make other histograms to see if they may be useful
#dist of ts length
tsdist_plot <- ggplot(shark_data, aes(x=tsrange))+ 
  geom_histogram(aes(y=..density..),binwidth = 10, center = 5, colour= "black",fill = "white") +
  geom_vline(aes(xintercept = mean(tsrange)), color= "red", linetype= "dashed", size=1)+
  theme_classic()+
  labs(x="Length of time series (Years)", y="Proportion of Timeseries")+
  scale_x_continuous(breaks = seq(from=0, to=120, by= 20))
tsdist_plot

#dist of year
yeardist_plot <- ggplot(shark_data, aes(x=year))+ 
  geom_histogram(aes(y=..density..),binwidth = 10, center = 5, colour= "black",fill = "white") +
  geom_vline(aes(xintercept = mean(year)), color= "red", linetype= "dashed", size=1)+
  labs(x="Years", y="Proportion of Timeseries")+
  theme_classic()
yeardist_plot

#use ggoceanmaps to create a plot showing the global distribution of the data

#group data by latitude and longitude and count occurrences
#define map limits based on data extent
lat_min <- -90
lat_max <- 90
lon_min <- min(location_data$long, na.rm = TRUE)
lon_max <- max(location_data$long, na.rm = TRUE)

#plot the data on a global map with defined limits
ggoceanmap <- basemap(limits = c(lon_min, lon_max, lat_min, lat_max)) +
  geom_point(data = location_data, aes(x = long, y = lat, size = Count), 
             color = "blue", alpha = 0.6) +
  scale_size_continuous(range = c(1, 10)) +  # Adjust circle size range as needed
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue"), # Set sea color
    panel.border = element_blank(),                     # Remove panel border
    plot.background = element_rect(fill = "white", color = NA),  # Remove outer border
    legend.background = element_rect(fill = "white"),   # Legend background
    axis.line = element_blank(),                        # Remove axis lines
    panel.grid = element_blank()                        # Remove gridlines
  ) +
  labs(title = "", size = "Number")  # Global Distribution of Elasmobranch Data
ggoceanmap

#create a figure that describes the distribution of data in the dataset being used
figure2 <- ggarrange(order_ts_plot, tsdist_plot, yeardist_plot, ggoceanmap, labels= c("A", "B", "C", "D"), ncol = 2, nrow = 2)
figure2 

#Section 10 -----test if distance to land can be a predictor of slope decrease------

#add dist 2 land column with ggoceanmaps
dist2land_data <- dist2land(shark_data, lon = ncol(9), lat = ncol(8))
shark_data <- dist2land_data

#make a model
gldist_model <- glmer(binomial_slope~ldist+(1|species), data=shark_data, family = binomial("logit"))
summary(gldist_model)
gldist_modelx <- simulateResiduals(gldist_model)
plot(gldist_modelx)

gldistplot<- plot_model(gldist_model, type= "pred", terms = c("ldist[all]"), color= "black")+  
  labs(
    title = "",
    x = "Distance From Shore (km)",
    y = "Slope of Trend in Abundance"
  )+ 
  theme_classic()
gldistplot

#make a histogram of the distribution of the distances
ldist_plot <- ggplot(shark_data, aes(x=ldist))+
  geom_histogram(aes(y=..density..),binwidth = 20, center = 5, colour= "black",fill = "white") +
  geom_vline(aes(xintercept = mean(ldist)), color= "red", linetype= "dashed", size=1)+
  theme_classic()+
  geom_density(alpha=.2, fill="#334bed")+
  labs(x="Ldist (km)", y="Density")
ldist_plot

#Section 11 -------split data into sections to see if there is any change over the years---------

#make a loop to subset and test this for 20 year slots moving along one year at a time

#initialize an empty list to store the results
model_results <- list()

#set the initial starting and ending years
start_year <- 1950
end_year <- 2020

#loop through overlapping 20-year windows
for (current_start in start_year:(end_year - 19)) {
  
  #define the end of the 20-year window
  current_end <- current_start + 19
  
  #filter the data for this 20-year window
  decade_data <- shark_data %>%
    filter(year >= current_start & year <= current_end)
  
  #check if there's enough data for the model
  if (nrow(decade_data) > 1) {
    
    #fit a linear model (e.g., predicting z_value from year)
    model <- lm(z_value ~ year, data = decade_data)
    
    #store relevant results
    model_results[[length(model_results) + 1]] <- data.frame(
      start_year = current_start,
      end_year = current_end,
      slope = coef(model)["year"],  #slope of the model
      intercept = coef(model)["(Intercept)"],  #intercept of the model
      p_value = summary(model)$coefficients[2, 4],  #p-value for the slope
      r_squared = summary(model)$r.squared,  #R-squared value
      std_error = summary(model)$coefficients[2,2]  #standard error
    )
    
  } else {
    #store NA if there's not enough data
    model_results[[length(model_results) + 1]] <- data.frame(
      start_year = current_start,
      end_year = current_end,
      slope = NA,
      intercept = NA,
      p_value = NA,
      r_squared = NA,
      std_error = NA
    )
  }
}

#combine the list into a data frame
model_results_df <- do.call(rbind, model_results)

#view the results
print(model_results_df)

#plot this
abundancetime_plot <- ggplot(model_results_df, aes(x = start_year, y = slope)) +
  geom_ribbon(aes(ymin=(slope-std_error), ymax=(slope+std_error)), fill = "grey") +
  geom_line(colour = "blue") +
  geom_hline(yintercept = 0, linetype = "solid") +  #add a solid line for y-axis at y=0
  theme_classic() +
  labs(title = "", #sliding Window Analysis of Elasmobranchs Abundance Trends
       x = "Start Year",
       y = "Slope of Trend in Abundance") 
abundancetime_plot
 

#create 2 more ggoceanmaps for each half of the data according to year

#subset the data for 1905-1965 and 1965-2020
data_1905_1972 <- shark_data %>% filter(year >= 1905 & year <= 1972)
data_1973_2020 <- shark_data %>% filter(year > 1973 & year <= 2020)


plot_slope_map <- function(shark_data, title) {
  lat_min <- -90
  lat_max <- 90
  lon_min <- min(location_data$long, na.rm = TRUE)
  lon_max <- max(location_data$long, na.rm = TRUE)
  
  basemap(limits = c(lon_min, lon_max, lat_min, lat_max)) +
    geom_point(data = shark_data, 
               aes(x = long, y = lat, size = abs(slope), 
                   color = slope > 0),  #logical value: TRUE for positive slope, FALSE for negative
               alpha = 0.7) +
    scale_size_continuous(range = c(1, 10)) +  #adjust size range as needed
    scale_color_manual(
      values = c("red", "green"),  #red for negative, green for positive
      labels = c("Negative", "Positive"),  #custom labels
      na.translate = FALSE  #remove NA from the legend
    ) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "lightblue"),  #set sea color
          plot.background = element_rect(fill = "white"),  #background for the whole plot
          legend.background = element_rect(fill = "white")) +
    labs(title = title, size = "Slope Magnitude", color = "Slope Sign")
}

#now name the function for each dataset
#plot the map for 1905-1965
ggoceanmap2 <- plot_slope_map(data_1905_1965, "1905-1965") +
  labs(x = "", y = "") +
  theme(
    panel.background = element_rect(fill = "lightblue"), # Set sea color
    panel.border = element_blank(),                     # Remove panel border
    plot.background = element_rect(fill = "white", color = NA),  # Remove outer border
    legend.background = element_rect(fill = "white"),   # Legend background
    axis.line = element_blank(),                        # Remove axis lines
    panel.grid = element_blank(),                        # Remove gridlines
    legend.position = "none"                             #remove legends
  )
ggoceanmap2
#plot the map for 1965-2020
ggoceanmap3 <- plot_slope_map(data_1965_2020, "1965-2020") +
  scale_color_manual(values = c("red", "blue")) + # Replace "green" with "blue"
  labs(x = "", y = "") +
  theme(
    panel.background = element_rect(fill = "lightblue"), # Set sea color
    panel.border = element_blank(),                     # Remove panel border
    plot.background = element_rect(fill = "white", color = NA),  # Remove outer border
    legend.background = element_rect(fill = "white"),   # Legend background
    axis.line = element_blank(),                        # Remove axis lines
    panel.grid = element_blank(),                        # Remove gridlines
    legend.position = "none"
  )
ggoceanmap3

#Section 12 -----create a figure showing how abundance has changed across time and that the distance from shore is a predictor-------
figure3 <- ggarrange(abundancetime_plot, gldistplot, ggoceanmap2, ggoceanmap3, labels= c("A", "B", "C", "D"), ncol = 2, nrow = 2)
figure3 <- annotate_figure(figure3,
  bottom = text_grob("Longitude (degrees)", size = 12, vjust = 0.01),
  left = text_grob("Latitude (degrees)", size = 12, , rot = 90, vjust = 0.5))

figure3 
