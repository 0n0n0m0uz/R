
install.packages("BiocManager")
BiocManager::install("Rgraphviz")
install.packages('bnlearn')

library(bnlearn)
library(emdbook)
library(Rgraphviz)




##############################
# What is a Bayesian Network?
# Lizards!!
# Courtest of the bnlearn package
##############################

# get data
data("lizards")
summary(lizards)

# Create and plot the expected network structure.
lizmod <- model2network("[Species][Diameter|Species][Height|Species]") 
graphviz.plot(lizmod)
lizmod


# Fit the network! Find the best parameter values

# Maximum liklihood estimate
lizfit <- bn.fit(lizmod, lizards, method = 'mle')
lizfit

# Bayesian parameter estimation
lizfit2 <- bn.fit(lizmod, lizards, method = 'bayes')
lizfit2

### A little different because they're two different ways for estimation!

##############################
# Lizards!! (cont.)
# Model Learning
##############################

##### BNs must satisfy the local Markov property
#### Essentially, this asks about conditional independence

#### Hill-climbing and incremental association Markov blankets measure this in different ways

### Hill-climbing uses d-seperation
## A kinda complicated idea with links below
dsep(lizfit, 'Height', 'Diameter', 'Species')
## This example staistifies d-seperation
# Two nodes are not connected given a parent node 


#### But, iamb uses conditional independence tests, such as mutual information before finding the rest of the network using d-seperation
### Mutual information quantifies how much information you can get from one variable 
#    from observation of another variable

# Test if variables are conditionally independent
# i.e. Given the species, are height and diameter then independent?
ci.test("Height", "Diameter", "Species", test = "mi", data = lizards) 
### No!



##### Let's go back to our lizard example....

graphviz.plot(lizmod)

#### What do machine learning algorithms think compared to our model?

# Find our model structure!

# hc() - a score-based algorithm
# hill-climbing algorithm
liznet <- hc(lizards)
liznet

graphviz.plot(liznet)
### The same!


# iamb() - a constraint-based algorithm 
# Incremental association Markov blanket
# These may be better for smaller, discrete datasets (there's some debate and emerging research!)
liznet2 <- iamb(lizards)
liznet2

graphviz.plot(liznet2)
### There's no directionality here....

graphviz.compare(liznet,liznet2)

##### Why are these so different?



# #### How can we decide what model to use?
# ### We can try and scoring it
# 
# # By default, score uses BIC
# bnlearn::score(liznet, data = lizards)
# bnlearn::score(liznet2, data = lizards)
# 
# ### Ack! One isn't fully directed, so we can use score!
# ### What else could we do?
# 
# ## Model validation methods may be helpful...
# 
# # Here, we use k-fold validation to measure expected loss
# # Expected loss: how far is the estimate from the true value?
# 
# bn.cv(lizards, liznet)
# bn.cv(lizards, liznet2)
# # The hill-climbing algorithm has a lower expected loss. Barely.
# 
# # Run it with more group splitting?
# bn.cv(lizards, liznet, k = 100)
# bn.cv(lizards, liznet2, k = 100)
# # Still close!
# 
# # More runs?
# bn.cv(lizards, liznet, k = 100, runs = 100)
# bn.cv(lizards, liznet2, k = 100, runs = 100)
# 
# ### The expected loss is essential identical
# # Q1: What does the closness of the expected loss imply?




##############################
# Let's try and fit a slightly more complicated example!
# Damselfish Settlement
# Courtesy of Bolker book!
##############################

# get data and make easier to work with
data("DamselSettlement")
summary(DamselSettlement)

damsel <- DamselSettlement
# all values must be numeric or integer for learning algorithms to learn
damsel[,c(2,3)] <- as.numeric(unlist(DamselSettlement[,c(2,3)]))


# Find structure

# hill-climbing algorithm
damnet <- hc(damsel)
graphviz.plot(damnet)

# iamb() - a constraint-based algorithm 
damnet2 <- iamb(damsel)
graphviz.plot(damnet2)

graphviz.compare(damnet,damnet2)

#### You get the same result!
### Q2: Why are these structures the same? 
#    hint: compare hc() and iamb() using the dsep and ci.test functions to test independence of children from grandparents given parents
### Challenge Q3: Why is dsep(damfit, 'site', 'obs', 'density') false?
#    hint: check out 'Bayesian Networks without Tears'!


### Sweet! We got agreement between our models!
## Let's optimize some conditional probability tables!

# Maximum likelihood method (Bayesian parameter estimation only for discrete networks (for now!))
damfit <- bn.fit(damnet, damsel, method = 'mle')
damfit

#### Check it!




##############################
# MyxoTiter Data!
# Our favorite! Courtesy of Bolker book
##############################

data("MyxoTiter_sum")
summary(MyxoTiter_sum)

titer <- MyxoTiter_sum
titerclass <- unlist(lapply(titer, is.integer))  
titer[,which(titerclass == TRUE)] <- as.numeric(unlist(titer[,which(titerclass == TRUE)]))


# Find structure

titerRelations <- hc(titer) 
graphviz.plot(titerRelations)

titerRelations2 <- iamb(titer) 
graphviz.plot(titerRelations2)

graphviz.compare(titerRelations,titerRelations2)

### They're different! Let's figure out which one is better...

# We can use BIC for this one!
bnlearn::score(titerRelations, data = titer)
bnlearn::score(titerRelations2, data = titer)

## The IAMB algorithm did a better job! Let's use that one to fit our network

# Fit the Bayesian network
titerfit <- bn.fit(titerRelations2, titer, method = 'mle')
titerfit


#### We can now find the probability of things happening, given evidence!

# What's the probability that the viral titer has a grade over 3, given it's been less than two days of infection?
cpquery(titerfit, grade > 3, day < 2)

## What about the other way?
# What's the probability that the individual has been infected for less than two days, given the grade is more than 3?
cpquery(titerfit, day < 2, grade > 3)




##############################
# Even bigger example!
# Seed Predation
# Courtesy of Bolker book
##############################

# Get data and clean it!
data("SeedPred")
summary(SeedPred)

pred <- na.omit(SeedPred)
predclass <- unlist(lapply(pred, is.integer))  
pred[,which(predclass == TRUE)] <- as.numeric(unlist(pred[,which(predclass == TRUE)]))
pred <- subset(pred, select = -(c(date, tint, tcum))) # network fitting has issues with dates
### Challenge Q4: How can yo make it so you can still include date?


# Find structure
prednet <- hc(pred)
graphviz.plot(prednet)

prednet2 <- iamb(pred)
graphviz.plot(prednet2)

### Way different!
## Why not?

# Hill climbing is less sensitive to individual error and uses different test for conditional independence

graphviz.compare(prednet,prednet2)

## Can we score them?
# Use a hybrid scoring algorithm
bnlearn::score(prednet, data = pred, type = "bic-cg")
bnlearn::score(prednet2, data = pred, type = "bic-cg")
# Again, partially directed networks can't be scored

bn.cv(pred, prednet)
bn.cv(pred, prednet2)
# Hill-climbing wins!


# Fit the Bayesian network
predfit <- bn.fit(prednet, data = pred, method = 'mle')
predfit

graphviz.plot(predfit)


### Relevant questions?

# Here's one:
#  What's the probability that there are more than 1 seeds taken, given the tree species is D. abyssinica?
cpquery(predfit, taken > 1, species == 'abz')


#### What if we have prior info about a dataset we want to include?

### Let's say we know there is an arc between seeds and taken?
prednet3 <- hc(pred, whitelist = c('seeds','taken'))
graphviz.plot(prednet3)

#### It updates!!!!!


#######################################
# custom.fit() - expert opinion
# Generally not ideal, TBH

### Creating our own example!
### Polar bear presence in relation to disease prevalence and habitat quality
### Based off USFWS BN!!
#######################################

disease <- matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("LOW", "HIGH")))
disease

habitatQuality <- matrix(c(0.8, 0.2), ncol = 2, dimnames = list(NULL, c("GOOD", "BAD")))
habitatQuality

conditionalProb <- c(0.5, 0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8)
dim(conditionalProb) <- c(2, 2, 2)
dimnames(conditionalProb) <- list("Presence" = c("TRUE", "FALSE"), "Disease" =  c("LOW", "HIGH"), 
                                  "Habitat_Quality" = c("GOOD", "BAD"))
conditionalProb

net <- model2network("[Disease][Habitat_Quality][Presence|Disease:Habitat_Quality]")
dfit <- custom.fit(net, dist = list(Disease = disease, Habitat_Quality = habitatQuality, 
                                    Presence = conditionalProb))
dfit
graphviz.plot(dfit)

bn.fit.barchart(dfit$Presence, ylab = 'Presence', xlab = 'Probability of Presence')

# How likely is it polar bears are present, given it's a poor quality habitat?
cpquery(dfit, Presence == 'TRUE', Habitat_Quality == 'BAD')

















