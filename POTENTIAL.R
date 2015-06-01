### Basic model selection and multilevel model fitting for POTENTIAL submodel 
#############################################################################

# 1. Input only variables that are not highly correlated
# 2. Will fit model with only varying intercepts by level (will not vary slopes)
# 3. Need to add output for warnings or errors
# 4. Models will not be automatically checked


# GENERIC CODE ####

# Import data
potential.data = read.csv('C:/PATH/DATA.csv')

# Load required libraries
library(MuMIn)
library(lme4)

# Create global model with all variables
global.model<-glmer(change~var0+var1+var2+var3+var4+var5+var6+var7+var8+var9+(1|county), family=binomial, data=potential.data)

# Create all possible models with a minimum of 3 and maximum of 7 variables, always include county as the level
select.model<-dredge(global.model, evaluate=TRUE, rank="AIC", fixed=~(1|county), m.min=3, m.max=7, trace=FALSE)

# Save the best model
model.best<-get.models(select.model, 1)
print(model.best)

# Rerun the best model and extract coefficients (will need to extract the correct variables from the best model and insert them here)
model<-glmer(change~var0+var1+var2+var3+var4+(1|county), family=binomial, data=potential.data)
summary(model)
coef(model)



# WORKING EXAMPLE ####

# Import Data
potential.data = read.csv('C:/<INSERT PATH HERE>/data_network.csv')

# Load required libraries
library(MuMIn)
library(lme4)

# Create global model with all variables
global.model = glmer(convert ~ slope + d2urban + d2water + d2protect + d2inter + d2rds + gdp10_2 + (1|county), family = binomial, data=potential.data)

# Create all possible models with a minimum of 3 and maximum of 7 variables, always include county as the level
select.model<-dredge(global.model, evaluate=TRUE, rank="AIC", fixed=~(1|county), m.min=3, m.max=7, trace=FALSE)

# Save the best model
model.best<-get.models(select.model, 1)
print(model.best)

# Rerun the best model and extract coefficients
model = glmer(convert ~ slope + d2urban + d2inter + d2rds + gdp10_2 + (1|county), family = binomial, data=potential.data)
summary(model) 
coef(model) 
