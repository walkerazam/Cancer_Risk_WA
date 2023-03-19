library(knitr)
# Setting up RMarkdown
opts_chunk$set(collapse = TRUE, tidy = TRUE, fig.align='center',
               tidy.opts = list(blank = TRUE, strip.white = TRUE), warning = FALSE, 
               message = FALSE, cache = FALSE, echo=FALSE)

# Loading Libraries
library(RColorBrewer) # creates nice color schemes
library(classInt)     # finds class intervals for continuous variables
library(ggplot2)
library(rgdal)
library(SUMMER)
library(maps)
library(maptools)
library(dplyr)
library(sf)
library(INLA)
library(prioritizr)
library(truncnorm)

# Reading in WA state County Shapefiles
# Data from: https://geo.wa.gov/datasets/wadnr::wa-county-boundaries/explore?location=47.182649%2C-120.817600%2C7.57
washington <- rgdal::readOGR(dsn="WA_County_Boundaries",layer="WA_County_Boundaries")

# Converting WA census tract to a polygon DF
wmap <- SpatialPolygons(washington@polygons) # getting Geometries
washington$geometry <- st_as_sfc(wmap) # adding geometry column
wmap <- st_as_sf(washington)
# selecting important columns
wmap <- wmap[,c("JURISDIC_2", "geometry")]
colnames(wmap) <- c('County', "geometry")

# Reading cancer data (pre-cleaned)
cancer <- read.csv("Lung_Cancer_Deaths_2015-2019.csv", header = TRUE)
# Selecting Columns
cancer <- cancer[, c("County", "Count", "Population", "Age.Adjusted.Rate.per.100.000")]

# Dropping the 3 counties with missing data:
cancer$Age.Adjusted.Rate.per.100.000[is.na(cancer$Age.Adjusted.Rate.per.100.000)] <- 36.4
cancer$Count[is.na(cancer$Count)] <- 0
#cancer <- na.omit(cancer)

# Merging with spatial object
merge <- merge(cancer, wmap, by.x='County', by.y='County')
wa <- st_as_sf(merge)

# Expected Morbidity:
wa$Expected <- wa$Population / 100000 * wa$Age.Adjusted.Rate.per.100.000

# Mapping The Populations per County
pal = function(n) brewer.pal(n, "Oranges")
plot(wa["Population"], pal = pal, nbreaks = 8, breaks = "equal",
     main=NULL)

# Mapping Mortality per County
pal = function(n) brewer.pal(n, "Purples")
# Assigning it to be Numeric
wa$Count <- as.numeric(wa$Count)
plot(wa["Count"], pal = pal, nbreaks = 8, breaks = "equal",
     main=NULL)

table <- cbind('Male Population Total' = sum(wa$Population), 'Mortality Count Count' = sum(wa$Count))
knitr::kable(table, caption = "Counts")

options(scipen=999) # no sci.notation
# Adding SMR column
wa$SMR <- wa$Count/wa$Expected
pal = function(n) brewer.pal(n, "BuPu")
plot(wa["SMR"], pal = pal, nbreaks = 8, breaks = "equal", main=NULL)

# Adding SE column
pal = function(n) brewer.pal(n, "BuPu")
# Calculating Standard Errors
wa$se <- sqrt(wa$SMR/wa$Expected)
plot(wa["se"], pal = pal, nbreaks = 8, breaks = "equal", main=NULL)

# Viewing a table of the counties with highest SMRs
top10 <- head(arrange(wa,desc(SMR)), n = 5)
top10 <- subset(top10, select = c("County", "SMR", "se"))
top10 <- st_drop_geometry(top10)
colnames(top10) <- c('County', 'SMR', 'Standard Error')
knitr::kable(top10, caption = "5 Counties with Highest SMRs")

# Fit Poisson-Lognormal model in INLA:
model.fit0 <- inla(Count ~ 1 + f(County, model="iid"),
                   data=wa, family="cenpoisson", 
                   E=Expected, 
                   control.family = list(cenpoisson.I = c(0, 5)),
                   control.predictor = list(compute = TRUE), 
                   control.compute= list(return.marginals=TRUE))

### beta: model.fit0$summary.fixed
### sigmˆ2a_e: model.fit0$summary.hyper

# Showing results as a table
iid_results <- data.frame("Posterior Median"=c(model.fit0$summary.fixed$`0.5quant`, 1/sqrt(model.fit0$summary.hyper$`0.5quant`)),
                     "Lower CI"=c(model.fit0$summary.fixed$`0.025quant`, 1/sqrt(model.fit0$summary.hyper$`0.975quant`)),
                     "Upper CI"=c(model.fit0$summary.fixed$`0.975quant`, 1/sqrt(model.fit0$summary.hyper$`0.025quant`)), check.names = F)
rownames(iid_results) = c("Beta","Sigma")
knitr::kable(round(iid_results,2), "simple", caption = "Poission Lognormal Non-Spatial")

options(scipen=999) # no sci.notation
# Extracting the posterior medians using our fitted model
wa$fit0fitted <- model.fit0$summary.fitted.values$`0.5quant`
wa$fit0se <- model.fit0$summary.fitted.values$`sd`

# Mapping our results [SUPPLEMENTARY]

library(gridExtra)
# Plotting Risk Estimates vs SMR
p1 <- ggplot(data.frame(iid=wa$fit0fitted, smr=wa$SMR),
       aes(y=iid, x=smr)) + geom_point() + labs(y="Poisson-Lognormal Estimates", x="SMR") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + xlim(0,3) + ylim(0,3) + theme_classic()

# Plotting standard deviations vs standard errors
p2 <- ggplot(data.frame(psd=model.fit0$summary.fitted.values$`sd`,se=wa$se),
       aes(y=psd,x=se)) + geom_point() + labs(y="IID Standard Deviation",x="SMR Standard Error") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + xlim(0,0.6) + ylim(0,0.6) + theme_classic()

grid.arrange(grobs = list(p1, p2), ncol = 2)

# compute adjacency matrix
#sf_use_s2(FALSE)
mat <- adjacency_matrix(wa)  
# Setting row and col name to HRA names in our neighbor matrix
colnames(mat) <- rownames(mat) <- wa$County
mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])

# Assigning a region column with a numeric identifier
wa$region <- 1:nrow(wa)
# Running INLA with given specifications (no covariates) Default Priors
formula <- Count ~ 1 + f(region, model="bym2", graph=mat)

model.fit1 <- inla(formula, data=wa, family="cenpoisson",E=Expected, control.family = list(cenpoisson.I = c(0, 5)), control.predictor=list(compute=TRUE),
                   control.compute=list(return.marginals.predictor=TRUE, config = TRUE))

spatial_result <- data.frame("Posterior Median"=c(model.fit1$summary.fixed$`0.5quant`, 
                                           1/(model.fit1$summary.hyper$`0.5quant`[1]), 
                                           model.fit1$summary.hyper$`0.5quant`[2]),
                      "Lower CI"=c(model.fit1$summary.fixed$`0.025quant`,
                                   1/(model.fit1$summary.hyper$`0.975quant`[1]),
                                   model.fit1$summary.hyper$`0.025quant`[2]),
                      "Upper CI"=c(model.fit1$summary.fixed$`0.975quant`,
                                   1/(model.fit1$summary.hyper$`0.025quant`[1]),
                                   model.fit1$summary.hyper$`0.975quant`[2]),
                      check.names = F)

rownames(spatial_result) = c("Beta", "Total Variance","Proportion Spatial")
knitr::kable(round(spatial_result,2), "simple", caption = " Spatial Poission Lognormal")

pal = function(n) brewer.pal(n,"BuPu")
# Extracting posterior median values and sd
wa$model1fitted <- model.fit1$summary.fitted.values$`0.5quant`
wa$fit1se <- model.fit1$summary.fitted.values$`sd`

# plotting
plot(wa["model1fitted"], pal = pal, nbreaks=8, breaks = "equal", main=NULL)
plot(wa["fit1se"], pal = pal, nbreaks = 8, breaks = "equal", main=NULL)

# Plotting Posterior Median vs SMR
p1 <- ggplot(data.frame(pmedian=model.fit1$summary.fitted.values$`0.5quant`,SMR=wa$SMR),
       aes(y=pmedian,x=SMR, size = wa$se)) + geom_point(alpha = 0.5) + 
  labs(y="Spatial Posterior Medians",x="SMR", size='SMR error') +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + xlim(0, 2.5) + ylim(0, 2.5) + theme_classic()

p1

p2 <- ggplot(data.frame(psd=wa$fit1se, se=wa$se),
       aes(y=psd, x=se)) + geom_point() + labs(y="Spatial Standard Deviation",x="SMR Standard Error") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + xlim(0, 0.6) + ylim(0, 0.6) + theme_classic()

p3 <- ggplot(data.frame(psd=wa$fit1se, se=wa$fit0se),
       aes(y=psd, x=se)) + geom_point() + labs(y="Spatial Standard Deviation",x="Non-Spatial Standard Deviation") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + xlim(0, 0.4) + ylim(0, 0.4) + theme_classic()

grid.arrange(grobs = list(p2, p3), ncol = 2)

# Reading in access to personal health care provider, 2015-2019
provider <- read.csv("healthcare_provider.csv", header = TRUE)
provider <- provider[, c('County', 'Count')]
colnames(provider) <- c('County', 'Provider.Count')

income <- read.csv("Median_Household_Income.csv", header = TRUE)
income <- income[, c('County', 'Median.Household.Income')]
# using 1000s of the income to let INLA get good enough initial values
income$Median.Household.Income <- income$Median.Household.Income / 1000

# merging 
merge <- merge(wa, provider)
merge <- merge(merge, income)

# Adding all three covariates
formula <- Count ~ 1 + Provider.Count + Median.Household.Income + f(region, model="bym2", graph=mat)

model.fit2 <- inla(formula, data=merge, family="cenpoisson", E=Expected, control.family = list(cenpoisson.I = c(0, 5)), control.predictor=list(compute=TRUE),
                   control.compute=list(return.marginals.predictor=TRUE, config = TRUE))

#model.fit2$summary.fixed[,1:5]
#model.fit2$summary.hyper[,1:5]

covariates_result <- data.frame("Posterior Median"=c(exp(model.fit2$summary.fixed$`0.5quant`), 
                                           1/(model.fit2$summary.hyper$`0.5quant`[1]), 
                                           model.fit2$summary.hyper$`0.5quant`[2]),
                      "Lower CI"=c(exp(model.fit2$summary.fixed$`0.025quant`),
                                   1/(model.fit2$summary.hyper$`0.975quant`[1]),
                                   model.fit2$summary.hyper$`0.025quant`[2]),
                      "Upper CI"=c(exp(model.fit2$summary.fixed$`0.975quant`),
                                   1/(model.fit2$summary.hyper$`0.025quant`[1]),
                                   model.fit2$summary.hyper$`0.975quant`[2]),
                      check.names = F)

rownames(covariates_result) = c("Beta", "Provider Access", "Median Household Income", "Total Variance","Proportion Spatial")
knitr::kable(round(covariates_result,4), "simple", caption = 'Spatial Model with Covariates')

# Extracting posterior median values and sd
wa$model2fitted <- model.fit2$summary.fitted.values$`0.5quant`
wa$fit2se <- model.fit2$summary.fitted.values$`sd`

p1 <- ggplot(data.frame(m1=wa$model1fitted, m2=wa$model2fitted),
       aes(y=m2, x=m1)) + geom_point() + labs(y="Spatial RR with Covariates",x="Spatial RR") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + theme_classic()

p2 <- ggplot(data.frame(psd=wa$fit1se, se=wa$fit2se),
       aes(y=se, x=psd)) + geom_point() + labs(y="Standard Deviations (with Covariates)",x="Standard Deviations (no Covariates)") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_abline(intercept=0,slope=1,color="red") + theme_classic()

grid.arrange(grobs = list(p1, p2), ncol = 2)

pal = function(n) brewer.pal(n,"BuPu")
# Extracting posterior median values and sd
wa$model2fitted <- model.fit2$summary.fitted.values$`0.5quant`
wa$fit2se <- model.fit2$summary.fitted.values$`sd`

# plotting
plot(wa["model2fitted"], pal = pal, nbreaks=8, breaks = "equal", main=NULL)
plot(wa["fit2se"], pal = pal, nbreaks = 8, breaks = "equal", main=NULL)
