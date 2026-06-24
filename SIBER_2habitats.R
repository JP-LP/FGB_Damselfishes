#######################
## June 08 -2026
## SIBER damselfish
## FGB, Texas
#######################

# loading libraries.
graphics.off()
rm(list=ls())
library(geometry)
library(ape)
library(rcdd)
library(bayestestR)
library(dplyr)
library(ggplot2)
library(ggdist)
library(ggtext)
library(nichetools)
library(purrr)
library(SIBER)
library(tidyr)
library(viridis)
library(skimr)
library(tidybayes)
library(hdrcde)
library(dplyr)
library(nicheROVER)

##Set working directory
setwd("/Users/jl82726/Documents/Casey_Damselfish/Data folder")

#Read dataset containing all metadata
damsel.metadata <- read.csv("Mihalek_RegalDamsel.csv") %>%
  filter(!is.na(Fish.ID))

#Read dataset containing isotopic data
damsel.isodata <- read.csv("Damselfish_Casey_April2026.csv") %>%
  rename(Fish.ID = Damsel.Fish.ID )

#Merge both datasets. Select useful columns and rename those with weird names
isotopic.dataset <-  damsel.isodata %>%
  left_join(damsel.metadata, by = "Fish.ID") %>%
  select(Fish.ID, Species, d15N,d13C, C.percentage, N.percentage, `Location.Type`,
         `Location.Specific`, , `Weight..g.`, `Total.Length..cm.`, `Standard.Length..cm.`) %>%
  rename(weight = `Weight..g.`) %>%
  rename(tl = `Total.Length..cm.`) %>%
  rename(sl = `Standard.Length..cm.`) %>%
  mutate(species = case_when(Species == "NC" ~ "N. cyanomos",
                             Species == "AM" ~ "A. multilineata",
                             Species == "SP" ~ "S. partitus"))

skim(isotopic.dataset)
levels(as.factor(isotopic.dataset$habitat))

#Merge FGB and Stetson bank as "Natural reefs"
#Delete N.cyanomos from HIA-474A due to unusual high 
#d15N values (n=13) + No other species were collected at this site.

#Also, prepare data to be read by the SIBER package.
#Create group and community indexes
damsel.data.sia <- isotopic.dataset %>% 
  mutate(habitat = case_when(Location.Type == "Stetson bank" ~ "Natural",
                             Location.Type == "FGB" ~ "Natural",
                             Location.Type == "Artificial reef" ~ "Artificial")) %>% 
  filter(Location.Specific != "HIA-474A") %>% 
  mutate(community = case_when(habitat == "Artificial" ~ 1,
                               habitat == "Natural" ~ 2,
                               TRUE ~ NA_real_)) %>%
  mutate(group = case_when(species == "A. multilineata" ~ 1,
                           species == "N. cyanomos" ~ 2,
                           species == "S. partitus" ~ 3,
                           TRUE ~ NA_real_)) %>%
  mutate(group = as.factor(group)) %>%
  mutate(community = as.factor(community)) %>%
  rename(iso1 = d13C) %>%
  rename(iso2 = d15N) %>%
  rename(group_name = species) %>%
  rename(community_name = habitat)

#compile community and group information
cg_names <- damsel.data.sia %>%
  distinct(community,
           group,
           community_name,
           group_name) %>%
  arrange(community, group)

#Filter data to be read by SIBER
damsel.sia <- damsel.data.sia %>%
  select (iso1, iso2, group, community) %>%
  drop_na() %>%  #drop all rows containing na's
  arrange(community, group) # reorder community groups in ascending order as SIBER required it that way
str(damsel.sia)

#Read data with SIBER
damsel.SIBER.data <- createSiberObject(damsel.sia)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(damsel.SIBER.data)
print(group.ML)

parms <- list()
parms$n.iter <- 1000000   # number of iterations to run the model for
parms$n.burnin <- 500000 # discard the first set of values
parms$n.thin <- 500    # thin the posterior by this many
parms$n.chains <- 3        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the
# means. Fitting is via the JAGS method.
#This function loops over each community and then loops over each
#group member, fitting a Bayesian multivariate (bivariate in this case)
#normal distribution to each group of data.
set.seed(27)
ellipses.posterior <- siberMVN(damsel.SIBER.data, parms, priors)

#This function loops over each group within each community and 
#calculates the posterior distribution describing the corresponding 
#Standard Ellipse Area.
SEA.B <- siberEllipses(ellipses.posterior)
head(SEA.B)

#Extract SEAb Credible intervals
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)
print(SEA.B.credibles)
# calculate bayesian mode of each group
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T) 
print(SEA.B.modes)

print(group.ML) #Glimpse of niche width

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('permille' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group")

#----- Reshaping SEA.B(Posterior of SEA) to print boxplot of niche width-----
SEA.Data <- as.data.frame(SEA.B)
head(SEA.Data)
SEA.Data <- rename(SEA.Data, #Based on print(group.ML) (community.group)
                   "1.1" = V1, 
                   "1.2" = V2,
                   "1.3" = V3, 
                   "2.1" = V4,
                   "2.2" = V5, 
                   "2.3" = V6)
head(SEA.Data)

# Organized dataset of allSEA posterior for each species on both habitats
SEA.b_df <- SEA.Data %>% 
  gather(key = community.group, value= "sea") %>% 
  mutate(Habitat = case_when(
    grepl("^1", community.group) ~ "Artificial",
    grepl("^2", community.group) ~ "Natural",
    TRUE ~ NA_character_)) %>% 
  mutate(Species = case_when(
    grepl("\\.1$", community.group) ~ "A. multilineata",
    grepl("\\.2$", community.group) ~ "N. cyanomos",
    grepl("\\.3$", community.group) ~ "S. partitus"))

head(SEA.b_df)

#Extract summary of results
seab.summary <- SEA.b_df %>%
  group_by(Habitat, Species) %>%
  summarize(SEAb = median(sea),
            ci.50.low =  hdr(sea, prob = 50)$hdr[,1],
            ci.50.high = hdr(sea, prob = 50)$hdr[,2],
            ci.95.low =  hdr(sea, prob = 95)$hdr[,1],
            ci.95.high = hdr(sea, prob = 95)$hdr[,2])
seab.summary

#Plot Ellipses for each group on each habitat

damsels.siber.plot <- ggplot(damsel.data.sia, aes(x=iso1, y=iso2,
                                                  colour=factor(group_name),
                                                  shape=factor(group_name), fill=factor(group_name)))+
  facet_wrap( ~ community_name)+
  geom_point(data = subset(damsel.data.sia, group_name == 'A. multilineata'), aes(x = iso1, y = iso2),
             shape = 21, colour="#4DBBD5",
             alpha=0.5, size =2) +
  geom_point(data = subset(damsel.data.sia, group_name == 'N. cyanomos'), 
             aes(x = iso1, y = iso2),
             shape = 21, colour="#00A087",
             alpha=0.5, size=2) +
  geom_point(data = subset(damsel.data.sia, group_name == 'S. partitus'), 
             aes(x = iso1, y = iso2),
             shape = 21, colour="#3C5488",
             alpha=0.5, size=2) +
  stat_ellipse(position="identity", level= pchisq(1,df = 2),
               segments=1000, linewidth = 1,
               type = "norm",
               geom = "polygon",
               alpha=0.7)+
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030'))+
  theme_bw()+
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.025, 0.375),
    legend.justification = c(0, 1),
    legend.text = element_text(face = "italic", size = 7),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    text = element_text(family = "sans"),
    axis.title = element_text(size= 10))+
  scale_colour_manual(values =c("#4DBBD5", "#00A087", "#3C5488"))+
  scale_fill_manual(values =c("#4DBBD5", "#00A087", "#3C5488"))

damsels.siber.plot

ggsave(filename = "Ellipse_damsels_2habss.png", width = 5, height = 3, dpi = 999, device = "png")
damsels.siber.plot
dev.off()

#### 
#Estimate overlap using nicheROVER
# Split communities dataset (Artificial reef, flower garden banks and stetson bank)
ar.damsel.sia <- damsel.sia %>% 
  filter(community == "1")

nat.damsel.sia <- damsel.sia %>% 
  filter(community == "2")

# Estimate mu and sigma (centroid and dispersion) for each species within each community

#Take all isotope measurements for one species, fit a Bayesian multivariate 
#normal niche model to those observations, and save 1000 posterior estimates 
#of that species' niche centroid and niche covariance matrix."
#tapply() then repeats that procedure for every species in the 
#community and stores all of the results in niche.par

# niw.post() is not fitting an MCMC model. Instead, it uses the fact that for 
#a multivariate normal likelihood with a Normal-Inverse-Wishart prior
#the function Calculates the posterior parameters of the NIW distribution from your data.
#Directly samples 1000 realizations from that posterior distribution.

#Artificial reefs
niche.par.ar <- tapply( #Applies a function for each index of a dataset (in this case for each group)
  1:nrow(ar.damsel.sia), #takes total number of row numbers
  ar.damsel.sia$group, # and splits them based on 
  function(ii) #extracts the isotopic value for a single species
    niw.post( #Fit a Normal Inverse-Wishart model ()
      nsamples = 10000, #Fit a Bayesian multivariate normal model
      X = as.matrix(ar.damsel.sia[ii, c("iso1", "iso2")]))) 

#Natural reefs
niche.par.fgb <- tapply( #Applies a function for each index of a dataset (in this case for each group)
  1:nrow(nat.damsel.sia), #takes total number of row numbers
  nat.damsel.sia$group, # and splits them based on 
  function(ii) #extracts the isotopic value for a single species
    niw.post( #Fit a Normal Inverse-Wishart model ()
      nsamples = 10000, #Fit a Bayesian multivariate normal model
      X = as.matrix(nat.damsel.sia[ii, c("iso1", "iso2")]))) 


##calculate the mean overlap metric between each species. 

#This function (overlap()) Calculates the distribution of a niche region 
#overlap metric for each pairwise species combination and user-specified niche region sizes (40%).
# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

ar.over.stat <- overlap(niche.par.ar, nreps = 10000, nprob = 10000, alpha = 0.95)
fgb.over.stat <- overlap(niche.par.fgb, nreps = 10000, nprob = 10000, alpha = 0.95)
                         
 #diagnostic plot                        
clrs <- c("black", "red", "blue", "orange") 
overlap.plot(ar.over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE)
overlap.plot(fgb.over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE)                       

#Convert Overlap estimates as dataframes to do a better plot
#Remember that the overlap metric is directional, such that it 
#represents the probability that an individual 
#from Species 𝐴 (rows) will be found in the niche of Species 𝐵(columns in the output from the overlap()).
df.ar.over <- as.data.frame(ar.over.stat) %>% 
  tibble::rownames_to_column(var = "Species A") %>% 
  gather(key= "Species B", value = "Prop.overlap", 2:30001) %>% 
  mutate(spp.a = case_when(`Species A` == "1" ~ "A. multilineata",
                           `Species A` == "2" ~ "N. cyanomos",
                           `Species A` == "3" ~ "S. partitus",
                           TRUE ~ NA_character_)) %>% 
  mutate(spp.b = case_when(grepl("^1", `Species B`) ~ "A. multilineata", #INDEX =group.community (extract starts with ^1)
                           grepl("^2", `Species B`) ~ "N. cyanomos",
                           grepl("^3", `Species B`) ~ "S. partitus",
                           TRUE ~ NA_character_)) %>% 
  filter(!is.na(Prop.overlap)) %>% 
  mutate(habitat = "Artificial")

df.fgb.over <- as.data.frame(fgb.over.stat) %>% 
  tibble::rownames_to_column(var = "Species A") %>% 
  gather(key= "Species B", value = "Prop.overlap", 2:30001) %>% 
  mutate(spp.a = case_when(`Species A` == "1" ~ "A. multilineata",
                           `Species A` == "2" ~ "N. cyanomos",
                           `Species A` == "3" ~ "S. partitus",
                           TRUE ~ NA_character_)) %>% 
  mutate(spp.b = case_when(grepl("^1", `Species B`) ~ "A. multilineata",
                           grepl("^2", `Species B`) ~ "N. cyanomos",
                           grepl("^3", `Species B`) ~ "S. partitus",
                           TRUE ~ NA_character_)) %>% 
  filter(!is.na(Prop.overlap)) %>% 
  mutate(habitat = "Natural")


#Combine all probalibilities into a single dataset and filter for the probabilities that
# represent the likelihood for an individual of Neopomacentrus cyanomos to be found in the
#niche space of one of the native species

overlap.probability <- rbind(df.ar.over, df.fgb.over) %>% 
  filter(spp.a == "N. cyanomos") %>% 
  mutate(Prop.overlap = as.numeric(Prop.overlap))

#Extract summary of posterior probability
Overlap.summary <- overlap.probability %>%
  group_by(habitat, spp.b) %>%
  summarize(median = median(Prop.overlap),
            ci.50.low =  hdr(Prop.overlap, prob = 50)$hdr[,1],
            ci.50.high = hdr(Prop.overlap, prob = 50)$hdr[,2],
            ci.95.low =  hdr(Prop.overlap, prob = 95)$hdr[,1],
            ci.95.high = hdr(Prop.overlap, prob = 95)$hdr[,2])
Overlap.summary

#Plot probability of overlap
Overlap.plot <- overlap.probability %>%
  ggplot(aes(x= spp.b, y = Prop.overlap,
             color = factor(spp.b),
             fill = factor(spp.b))) + 
  stat_dotsinterval(aes(slab_color = factor(spp.b)),
                    shape = 21, alpha =0.6,
                    .width = c(0.5, 0.95),
                    scale = 0.6, justification = -0.2)+
  facet_wrap(~ habitat)+
  theme_bw()+
  theme(
    strip.text = element_text(face = "bold"), #remove title text from panels
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    #legend.justification = c(0, 1),
    #legend.text = element_text(face = "italic", size = 8),
    #legend.title = element_text(face = "bold", size = 8),
    axis.text.x = element_text(size =8, face = "italic"),
    axis.ticks.length.x = unit(-0.1, "cm"),
    axis.title = element_text(size= 10))+
  scale_color_manual(values =c("#4DBBD5",  "#3C5488"),
                     aesthetics = c("color", "fill", "slab_color"))+ #very important line!
  scale_y_continuous(lim = c(0,1),
                     breaks= seq(0,1,0.20))+
  xlab("")+
  ylab(expression("Probability of overlap"))

Overlap.plot

#Extract plot
ggsave(filename = "2habs_Overlap_damsels.png", width = 5, height = 3, dpi = 999, device = "png")
Overlap.plot
dev.off()

##########
## Extract Layman's metrics at the assemblage level
###

# extract the posterior means
mu.post <- extractPosteriorMeans(damsel.SIBER.data, ellipses.posterior)
head(mu.post)

# calculate the corresponding distribution of layman metrics
layman.Bayes <- bayesianLayman(mu.post)

#Extract Layman metrics (using niche tools)
Layman.data <- extract_layman(layman.Bayes, community_df = cg_names[,1:2]) %>% 
  mutate(habitat = case_when(community == "1" ~ "Artificial",
                             community == "2" ~ "Natural",
                             TRUE ~ NA_character_)) %>%
  mutate(species = case_when(group == "1" ~ "A. multilineata",
                             group == "2" ~ "N. cyanomos",
                             group == "3" ~ "S. partitus",
                             TRUE ~ NA_character_)) 


#### Extract summary of posterior results of the layman metrics
Layman.summary <- Layman.data %>%
  group_by(metric, habitat) %>%
  summarize(median = median(post_est),
            ci.50.low =  hdr(post_est, prob = 50)$hdr[,1],
            ci.50.high = hdr(post_est, prob = 50)$hdr[,2],
            ci.95.low =  hdr(post_est, prob = 95)$hdr[,1],
            ci.95.high = hdr(post_est, prob = 95)$hdr[,2])
Layman.summary


Layman.cd.plot <-  ggplot(data = subset(Layman.data, metric == 'CD'), 
                          aes(x= habitat, y = post_est, color = factor(habitat),
                              fill = factor(habitat)))+
  stat_halfeye(aes(y = post_est), .width = c(.50, .95),
               alpha = 0.75, justification = -0.2)+
  # position = position_nudge(y = -0.2))+
  xlab("Reef type")+
  ylab(expression("MDC " ('\u2030')))+
  theme_bw()+
  theme(
    axis.text = element_text(face = "plain"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position ="none",
    legend.title = element_blank(),
    text = element_text(family = "sans"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
    axis.title = element_text(size= 10))+
  scale_y_continuous(lim = c(0,2.2),
                     breaks= seq(0,2.2,0.5))+
  scale_colour_manual(values =c("#CC79A7", "#EFBF04"))+
  scale_fill_manual(values =c("#CC79A7", "#EFBF04"))


Layman.cd.plot

ggsave(filename = "2habs_MDCplot.png", width = 5, height = 3, dpi = 999, device = "png")
Layman.cd.plot
dev.off()


#Now we are going to perform a comparison for each species among habitats
library(geometry)
library(ape)
library(rcdd)

#scaleSI_range01: function to standardize Stable Isotope values of a dataset given a reference dataset of organisms (individuals or group of individuals)
#
# INPUTS:- 'raw_data': a dataframe or matrix with stable isotope values for at least 2 elements (columns) for several groups (rows).
#						Columns' names should be a subset of ('d13C','d15N','dD','d34S').
#		 - 'all_data': a dataframe or matrix with stable isotope values for the same elements than 'raw_data' (columns) with an equal or larger set of samples (rows).
#				By default, 'all_data' is identical to 'raw_data'.
#
# OUTPUT: a matrix similar to 'raw_data' with scaled stable isotope values (each ranging from 0 to 1).
#         If standard deviation values are provided, they are scaled based on the range of mean values.

raw_data <- damsel.data.sia %>%
  select(iso1, iso2, group, community, group_name, community_name) %>% 
  rename(d13C = iso1,
         d15N = iso2)

str(raw_data)

####
#Since we  aims at comparing stable isotope diversity spatially, it is frequent
#that basal resources differ and/or have different stable isotope values. 
#In these cases, is better to pool the stable isotope values of ALL organisms
#studied before scaling them to give each isotope the same weight (instead of doing 
#the scaling independently for each ecosystem). Such a procedure will guarantee
#than the diversity of resources is accounted for in the computation of isotopic diversity

#So we are scaling isotopic values of all speices on both "Artifical" and "Natural" reef
#at the same time
# codes for stable isotope ratios
nm_si<-c("d13C","d15N","dD","d34S")
scaled.damsel.sia <- scaleSI_range01(raw_data, all_data=raw_data) %>% 
rename(iso1 = d13C,
       iso2 = d15N) 

#Preparefor SIBER
sc.damsel.sia <- scaled.damsel.sia %>% 
  select(iso1, iso2, group, community) %>% 
  arrange(community, group)

#Now, lets estimate SIBER ovelap estimates
damsel.SIBER.data.sc <- createSiberObject(sc.damsel.sia)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(damsel.SIBER.data.sc)
print(group.ML)

#Set parameter to fit ellipses
parms <- list()
parms$n.iter <- 1000000   # number of iterations to run the model for
parms$n.burnin <- 500000 # discard the first set of values
parms$n.thin <- 500    # thin the posterior by this many
parms$n.chains <- 3        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

#Fit bayesian ellipses for scaled data
set.seed(270)
ellipses.posterior.sc <- siberMVN(damsel.SIBER.data.sc, parms, priors)

#This function loops over each group within each community and 
#calculates the posterior distribution describing the corresponding 
#Standard Ellipse Area.
SEA.B.sc <- siberEllipses(ellipses.posterior.sc)
head(SEA.B.sc)

SEA.Data.ellipses <- rename(SEA.Data.scaled, #Based on print(group.ML) 
                   "1.1" = V1, 
                   "1.2" = V2,
                   "1.3" = V3, 
                   "2.1" = V4,
                   "2.2" = V5, 
                   "2.3" = V6) %>% as.matrix() 
head(SEA.Data.ellipses)



#diagnostic plot
siberDensityPlot(SEA.B.sc, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('permille' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group")

#Plot ellipses for each habitat on panels for each species
scaled.damsels.siber.plot <- ggplot(scaled.damsel.sia, aes(x=iso1, y=iso2,
                            colour=factor(community_name),shape=factor(community_name), 
                                                  fill=factor(community_name)))+
  facet_wrap( ~ group_name)+
  geom_point(data = scaled.damsel.sia, aes(x = iso1, y = iso2,
             colour=factor(community_name),fill=factor(community_name)),
             shape = 21,
             alpha=0.5, size =2)+
  stat_ellipse(position="identity", level= pchisq(1,df = 2),
               segments=1000, linewidth = 1,
               type = "norm",
               geom = "polygon",
               alpha=0.7)+
  xlab("Scaled"~ {delta}^13*C)+
  ylab("Scaled"~ {delta}^15*N)+
  theme_bw()+
  theme(
    strip.text = element_text(face = "bold.italic"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.05, 0.25),
    legend.justification = c(0, 1),
    legend.text = element_text(face = "plain", size = 7),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    text = element_text(family = "sans", size= 8),
    axis.title = element_text(size= 10))+
  scale_colour_manual(values =c("#CC79A7", "#EFBF04"))+
  scale_fill_manual(values =c("#CC79A7", "#EFBF04"))

scaled.damsels.siber.plot

ggsave(filename = "2habs_scaledSIBER.png", width = 5, height = 3, dpi = 999, device = "png")
scaled.damsels.siber.plot
dev.off()

# Define ellipses index (community.group) #Based on: SEA.Data.ellipses
str(SEA.Data.ellipses)
ellipse1 <- "1.1" # AR.AM
ellipse2 <- "1.2" # AR.NC
ellipse3 <- "1.3" #AR.SP
ellipse4 <- "2.1" #Nat.AM
ellipse5 <- "2.2" #Nat.NC
ellipse6 <- "2.3" #Nat.SP

# Calculate the relative overlap among scaled ellipses
AM.overlap <- bayesianOverlap(ellipse1, ellipse4, ellipses.posterior.sc,
                            draws = 1000, p.interval = 0.95, n = 100)

NC.overlap <- bayesianOverlap(ellipse2, ellipse5, ellipses.posterior.sc,
                              draws = 1000, p.interval = 0.95, n = 100)

SP.overlap <-  bayesianOverlap(ellipse3, ellipse6, ellipses.posterior.sc,
                               draws = 1000, p.interval = 0.95, n = 100)


#Calculate the proportion of overlap relative to the non overlapping area
AM.overlap.prop <- ((AM.overlap$overlap / (AM.overlap$area1 + AM.overlap$area2 - AM.overlap$overlap) 
                     * 100)) %>% as.data.frame() %>% mutate(comparison = "A. multilineata")

NC.overlap.prop <- ((NC.overlap$overlap/ (NC.overlap$area1+NC.overlap$area2-NC.overlap$overlap) 
                     * 100)) %>% as.data.frame() %>% mutate(comparison = "N. cyanomos")

SP.overlap.prop <- ((SP.overlap$overlap / (SP.overlap$area1+SP.overlap$area2-SP.overlap$overlap) 
                     * 100)) %>% as.data.frame() %>% mutate(comparison = "S. partitus")

#Convert into a tidy dataset 
sc.overlap.props <- bind_rows(
  AM.overlap.prop,
  NC.overlap.prop,
  SP.overlap.prop) %>% 
  rename(overlap.sc = ".")

#Extract summary of posterior estimates of overlap

sc.overlap.summary <- sc.overlap.props %>% 
  group_by(comparison) %>%
  mutate(overlap.scaled = as.numeric(overlap.sc)) %>% 
  summarize(median = median(overlap.sc),
            ci50.low = quantile(overlap.sc, 0.25),
            ci50.high = quantile(overlap.sc, 0.75),
            ci95.low = quantile(overlap.sc, 0.025),
            ci95.high = quantile(overlap.sc, 0.975))
sc.overlap.summary

####Estimate overlap with MaxLikOverlap function

AM.overlap.max <- maxLikOverlap("1.1", "2.1",
  damsel.SIBER.data.sc, p.interval = 0.95, n = 100)

NC.overlap.max <- maxLikOverlap("1.2","2.2", damsel.SIBER.data.sc,
  p.interval = 0.95, n = 100)

SP.overlap.max <- maxLikOverlap("1.3","2.3", damsel.SIBER.data.sc,
  p.interval = 0.95, n = 100)

#Calculate Jaccard overlap
AM.overlap.prop <- (AM.overlap.max[3] /
    (AM.overlap.max[1] + AM.overlap.max[2] - AM.overlap.max[3]) * 100)


#Calculate Jaccard overlap
AM.overlap.max <- (AM.overlap.max[3] /
                      (AM.overlap.max[1] + AM.overlap.max[2] - AM.overlap.max[3]) * 100)

NC.overlap.max <- (NC.overlap.max[3] /
                      (NC.overlap.max[1] + NC.overlap.max[2] - NC.overlap.max[3]) * 100)

SP.overlap.max <- (SP.overlap.max[3] /
                     (SP.overlap.max[1] + SP.overlap.max[2] - SP.overlap.max[3]) * 100)


AM.overlap.max
NC.overlap.max
SP.overlap.max

#########Perform a two way anova
N.mod <- lm(iso2 ~ group_name * community_name, data = damsel.data.sia)
anova(N.mod)
shapiro.test(residuals(N.mod))
car::leveneTest (iso2 ~ group_name * community_name, data = damsel.data.sia)


C.mod <- lm(iso1 ~ group_name * community_name, data = damsel.data.sia)
anova(C.mod)
shapiro.test(residuals(C.mod))
car::leveneTest(iso1 ~ group_name * community_name, data = damsel.data.sia)

#Compare differences among species for NITROGEN
library(emmeans)
# Compare species within each habitat
emmeans(N.mod, pairwise ~ group_name | community_name,
        adjust = "tukey")

# Compare habitats within each species
emmeans(N.mod, pairwise ~ community_name | group_name,
        adjust = "tukey")

#Compare differences among species for Carbon

# Compare species within each habitat
emmeans(C.mod, pairwise ~ group_name | community_name,
        adjust = "tukey")

# Compare habitats within each species
emmeans(C.mod, pairwise ~ community_name | group_name,
        adjust = "tukey")
