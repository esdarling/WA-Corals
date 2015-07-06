
# Created by:    Emily S. Darling
# Created:       29 June 2015
# Last modified: 29 June 2015
# Purpose:       reassign new WA surveys to coral life histories, FD variables

# List of variables to calculate: 
# *** check new Mouillot code for FD from December
# no_genera sumcover
# BleachSuscept	
# Competitive  StressTolerant  Weedy  Generalist	
# FRic  RaoQ	**from new Mouillot code
# CWMBranching	CWMPlating	CWMDomed	CWMMaxsize	CWMGrowthrate	CWMBrooding	CWMSpawning	CWMFecundity **from package FD


# load new surveys 
setwd("/Users/emilydarling/Dropbox/1-On the go/Workshops/WA Corals/IAS workshop May 2014/Data/Corals to update_June2015")
d <- read.csv("rf_data_use_workshop_w_AB-depth_JG2_JZ_PorBrMass_v2.csv", 
              header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)   
head(d)


# pull ID and cover by genus
d <- d[,c(1,16:81)]
head(d)

d2 <- melt(d,id.var = 1, variable.name = "genus", value.name = "cover")
head(d2)
nrow(d2)
hist(d2$cover)
min(d2$cover); max(d2$cover)


porites <- d2 %>%
  filter(genus == "Porites.massive")
min(porites$cover); max(porites$cover)
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

# upload traits, LH data for genus
setwd("/Users/emilydarling/Dropbox/1-On the go/Workshops/WA Corals/IAS workshop May 2014/Data/Coral Traits info")
lh_info <- read.csv("IAS Genera with traits_LH.csv", 
              header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)   

d5 <- left_join(d2, lh_info[,1:2])
#write.csv(d5, "test.csv", row.names = FALSE)

lh <- d5 %>%
  ungroup() %>%
  group_by(ID, LH) %>%
  summarize(cover = sum(cover))
levels(as.factor(lh$LH))
head(lh)


#Competitive
lh$Competitive <- ifelse(lh$LH == "Competitive", lh$cover, 
                         ifelse(lh$LH == "1C_1Gen",lh$cover * 1/2,
                                ifelse(lh$LH == "5C_10G_2ST", lh$cover * 5/17,0)))
head(lh)
hist(lh$Competitive)
min(lh$Competitive, na.rm = TRUE); max(lh$Competitive, na.rm = TRUE) 

#Stress-Tolerant
lh$StressTolerant <- ifelse(lh$LH == "Stress-tolerant", lh$cover, 
                            ifelse(lh$LH == "4ST_2W",lh$cover * 4/6,
                                   ifelse(lh$LH == "5C_10G_2ST", lh$cover * 2/17,0)))

hist(lh$StressTolerant)
min(lh$StressTolerant, na.rm = TRUE); max(lh$StressTolerant, na.rm = TRUE)   				

#Weedy			
lh$Weedy <- ifelse(lh$LH == "Weedy", lh$cover, 
                   ifelse(lh$LH == "4G_1W", lh$cover * 1/5,
                          ifelse(lh$LH == "4ST_2W", lh$cover * 2/6, 0)))

hist(lh$Weedy)
min(lh$Weedy, na.rm = TRUE); max(lh$Weedy, na.rm = TRUE)    

#Generalist			
lh$Generalist <- ifelse(lh$LH == "Generalist", lh$cover, 
                        ifelse(lh$LH == "1C_1Gen" | lh$LH == "2C_3Gen_1ST", lh$cover * 1/2,
                               ifelse(lh$LH == "4G_1W", lh$cover * 1/3, 
                                      ifelse(lh$LH == "5C_10G_2ST", lh$cover * 10/17, 0)))) 						

hist(lh$Generalist)
min(lh$Generalist, na.rm = TRUE); max(lh$Generalist, na.rm = TRUE)            

head(lh)
levels(as.factor(lh$LH))

#_NA - no life history classification
lh$no_data <- ifelse(is.na(lh$LH) == TRUE, lh$cover, 0)
hist(lh$no_data)
min(lh$no_data, na.rm = TRUE); max(lh$no_data, na.rm = TRUE) 

# sum to life history
head(lh)
lh %>%
  filter(ID == 252)

lh_melt <- melt(lh[,-c(2:3)], id.vars = 1, variable.name = "LH")
head(lh_melt)

lh2 <- lh_melt %>%
  ungroup() %>%
  group_by(ID, LH) %>%
  summarize(cover = sum(value, na.rm = TRUE))
head(lh2)

lh2 %>%
  filter(ID == 252)

lh3 <- dcast(lh2, ID~LH)
head(lh3)
hist(lh3$Competitive)
hist(lh3$no_data)

test <- lh3 %>%
  mutate(sumLH = rowSums(.[2:6])) %>%
  mutate(prop_no_data = no_data / sumLH) %>%
  arrange(desc(prop_no_data))

#site w/ 50% no data has 1.5% cover.. 
min(test$prop_no_data); max(test$prop_no_data)
mean(test$prop_no_data)  
hist(test$prop_no_data)      

head(lh3)

##########################
# calculate FD, traits
# FD has two parts: 'x' trait table + 'a' abundance table
# 'x' has genera as rows, traits as columns

setwd("/Users/emilydarling/Dropbox/1-On the go/Workshops/WA Corals/IAS workshop May 2014/Data/Coral Traits info")
x <- read.csv("IAS Genera with traits_LH_matched 5July2015.csv", 
                    header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)   
head(x)
names(x)
x$genus

x2 <- x
x2 <- x[,-c(2,11)]
row.names(x2) <- x2$genus   
x2 <- x2[,-1]
head(x2)
names(x2)
row.names(x2)


# load 'a' abundance file
setwd("/Users/emilydarling/Dropbox/1-On the go/Workshops/WA Corals/IAS workshop May 2014/Data/Corals to update_June2015")
d <- read.csv("rf_data_use_workshop_w_AB-depth_JG2_JZ_PorBrMass_v2.csv", 
              header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)   
head(d)
names(d)
a <- d[,c(1,16:81)]
head(a)
names(a)

a2 <- a
a2 <- a2[order(a2$ID),]
row.names(a2) <- a2$ID  
a2 <- a2[,-1]
head(a2)

# dbFD error: at least one species does not occur in community (zero total abundance)
# find species in WA dataset with 0 abundance.. 
head(a)
spp.test <- melt(a, id.var = 1, variable.name = "genus")
spp.test2 <- spp.test  %>% 
  ungroup() %>%
  group_by(genus) %>%
  summarize(colSums = sum(value)) %>%
  filter(colSums == 0)
spp.test2

head(spp.test)
a3 <- left_join(spp.test, spp.test2)
a3 <- a3[-which(a3$colSums == 0),]

a3_cast <- dcast(a3[,-4], ID ~ genus)
row.names(a3_cast) <- a3_cast$ID  
a3_cast <- a3_cast[,-1]
head(a3_cast)
a4 <- a3_cast


# combined x2, a2 in FD code (and cross fingers)
head(x2)
row.names(x2)

head(a4)
names(a4)

ias_fd <- dbFD(x2, a4, w = c(0.33,0.33,0.33,1,1,0.5,0.5,1), corr="cailliez", calc.CWM = TRUE) 
levels(as.factor(ias_fd$sing.sp))
head(ias_fd$CWM) 

CWM <- data.frame(ias_fd$CWM)
names(CWM) <- c("CWM.Branching","CWM.Plating","CWM.Domed","CWM.Maxsize","CWM.Growthrate",
                           "CWM.Brooding","CWM.Spawning","CWM.Fecundity")                
CWM$ID <- row.names(CWM)
CWM <- CWM[,c(9,1:8)]
head(CWM) 
names(CWM)

head(ias_fd$FRic)
head(ias_fd$RaoQ)

CWM$FRic <- ias_fd$FRic
CWM$RaoQ <- ias_fd$RaoQ
head(CWM)

plot(CWM$FRic ~ CWM$RaoQ)

CWM$ID <- as.numeric(CWM$ID)

# check w/new Mouillot code (did I get code?)
# CODE for FD
# Hard to understand - problem in the dist[sp,sp] matrix multiplicatin
# Using the code from dbFD in package "FD"

#Reading species data file

#data=read.csv(file="UniqueCorals + Traits_2029 sites_26Aug2014.csv")

#lab.com=xtabs(Perc_cover2~Site+Coral,data=data)
head(d)
names(d)
d2 <- d[,c(1,16:81)]
d2 <- d2[order(d2$ID),]
head(d2)

d3 <- melt(d2, id.var = 1, variable.name = "coral", value.name = "cover")
head(d3)

lab.com=xtabs(cover ~ ID + coral,data=d3)
head(lab.com)

# Traits
head(x2)
#trait=data[,c(5,7:14)]
#trait=unique(trait)
#rownames(trait)=trait$Coral
#trait=trait[,-1]
trait <- x2
dist.coral <- daisy(trait, metric = "gower", stand = TRUE)

# Rao et FRIC
# Based on relative cover
x <- as.matrix(lab.com)
S <- length(x[1, ])
N <- length(x[, 1])
total <- apply(x, 1, sum)
relab=x/total
head(relab)

res=matrix(0,nrow=N,ncol=2)

rownames(res)=rownames(lab.com)
colnames(res)=c("nb_taxa","RAO")

dist2<- as.matrix(dist.coral)


# creation d'axes quanti avec PCoA
# PCOA of coral axes based on traits 
coral.pcoa=cmdscale(dist.coral, k = 4)

plot(coral.pcoa)

head(relab)

#N is no. of sites (here, 444 in WA)
for (i in 1:2)
{
  
  #i=1
  #sp collects species with >0% abundance at each site
  sp<-which(relab[i,]!=0)
  #sp<-which(relab[1,]!=0)
  #sp
  
  #relab[1,"Acropora"] # relative percent cover of Acropora
  #dist2["Acropora","Acropora"]
  #head(dist2)
  
  # Rao's quadratic entropy (Q): sum of abundance-weighted distances
  # matrix multiplication
  d2ab <- relab[i,sp] %*% dist2[sp,sp] %*% relab[i,sp]
  
  
  # Q expressed as an equivalent number of species
  Q<- 1/(1-d2ab)
  
  res[i,2]=Q
  
  res[i,1]=length(sp)
  
}


write.table(res,file="FD-Corals.txt"



##########################
# create combined file by ID, merge back in after
# d3 = sumCover, no_Genera
# ssi = BleachSuscept
# lh3 = life histories
# CWM = traits, FD metrics

coral_vars <- left_join(d3,ssi)
coral_vars <- left_join(coral_vars, lh3[,-6])
coral_vars <- left_join(coral_vars, CWM)

#incorporate into final file
head(d)
names(d)
d_clean <- d[,c(1:9,16:81,95:141)]

coral_vars2 <- left_join(d_clean, coral_vars)
names(coral_vars2)

#plot(coral_vars2$sum_cover ~ coral_vars2$sumcover)

coral_vars2 <- coral_vars2[,c(1:9,124:139,10:122)]

plot(coral_vars2$BleachSuscept ~ coral_vars2$Competitive)


setwd("/Users/emilydarling/Dropbox/1-On the go/Workshops/WA Corals/IAS workshop May 2014/Data/Corals to update_June2015")
write.csv(coral_vars2, 
          "rf_data_use_workshop_coralupdate_esd_6Jul2015.csv", row.names = FALSE)

