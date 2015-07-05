
# Created by:    Emily S. Darling
# Created:       29 June 2015
# Last modified: 29 June 2015
# Purpose:       reassign new WA surveys to coral life histories, FD variables

# List of variables to calculate: 
# *** check new Mouillot code for FD from December
# Competitive  StressTolerant	Weedy	Generalist	BleachSuscept	no_genera
# FRic  FEve	FDis	FDiv	RaoQ	
# CWMBranching	CWMPlating	CWMDomed	CWMMaxsize	CWMGrowthrate	CWMBrooding	CWMSpawning	CWMFecundity


# load new surveys 
setwd("/Users/emilydarling/Dropbox/1-On the go/Workshops/WA Corals/IAS workshop May 2014/Data/Corals to update_June2015")
d <- read.csv("rf_data_use_workshop_w_AB-depth_JG2.csv", 
              header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)   
head(d)
nrow(d) 
names(d)

# pull ID and cover by genus
d <- d[,c(1,16:80)]
head(d)

d2 <- melt(d,id.var = 1, variable.name = "genus", value.name = "cover")
head(d2)
nrow(d2)

porites <- d2 %>%
  filter(genus == "Porites")
hist(porites$cover)

##################################
# calculate sumcover and no_genera
d3 <- d2 %>%
  group_by(ID) %>%
  mutate(sumcover = sum(cover)) %>%
  filter(cover != 0) %>%
  mutate(no_genera = n()) %>%
  summarise(sumcover = mean(sumcover), no_genera = mean(no_genera))

hist(d3$no_genera)
hist(d3$sumcover)

####################################
# calculate bleaching susceptibility
# load bleaching responses for each taxa
br <- read.csv("WIO taxa bleaching responses_30May2012.csv", 
               header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE) 
br
head(br)
head(d2)
levels(as.factor(d2$genus))

# join WA cover data to taxa bleaching responses, check join
d4 <- inner_join(d2, br[,1:2])
head(d4)
test <- d4 %>%
  ungroup() %>%
  group_by(genus) %>%
  summarise(br = mean(taxa_br))

# multiple by relative abundance and then sum to site
ssi <- d4 %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(total_cover = sum(cover)) %>%
  mutate(rel_cover = round((cover / total_cover), 3)) %>%
  mutate(suscept = rel_cover * taxa_br) %>%
  summarize(site_suscept = sum(suscept, na.rm = TRUE))

hist(ssi$site_suscept)
names(ssi)[2] <- "BleachSuscept"

##########################
# calculate life histories
# all corals are not grouped to genus - no lifeforms
# dammit cover combined for Porites too - how common is Porites? 
# check genus splits





# create combined file
coral_vars <- join(d3,ssi)
plot(coral_vars$sumcover ~ coral_vars$BleachSuscept)
text(coral_vars$BleachSuscept,coral_vars$sumcover, labels = coral_vars$ID )
