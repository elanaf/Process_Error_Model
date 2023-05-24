load("main_dfs.RData")

library(ggplot2)
library(magrittr)
library(dplyr)
library(gtools)


## Make a new dataframe that is in a better format ####
#only select columns that I need for the analysis
native.dat <- greenhouse %>%
  filter(Species == "HENU" | Species == "RUMA" | Species == "EPCI" | 
           Species == "EUMA" | Species == "SYCI" | Species == "SOCA" |
           Species == "MUAS" | Species == "PUNU" | Species == "DISP" | Species == "EUOC") %>%
  filter(Density == "L" & Phrag_Presence == "W")%>%
  select(Species, Block, Density, Phrag_Presence, Date_Cleaned, Cover.Native)  %>%
  arrange(Species) #put same species next to each other

#save this for graphing later
native.dat2 <- greenhouse %>%
  filter(Species == "HENU" | Species == "RUMA" | Species == "EPCI" | 
           Species == "EUMA" | Species == "SYCI" | Species == "SOCA" |
           Species == "MUAS" | Species == "PUNU" | Species == "DISP"| Species == "EUOC") %>%
  filter(Density == "L" & Phrag_Presence == "W")%>%
  select(Species, Block, Density, Phrag_Presence, Date_Cleaned, Cover.Native)  %>%
  arrange(Species) #put same species next to each other

native.dat$ID_Col <- with(native.dat, paste0(Species, Block))

#make vectors to keep track of the species
#each species gets a number 1 - 13, so I repeat those numbers 3times for the 3 blocks
species.vec <- rep(1:10, each = 3)

#rearrange so date is column
native.dat <- reshape2::dcast(native.dat, ID_Col ~ Date_Cleaned, value.var = "Cover.Native")

#now make a matrix
native.mat <- as.matrix(native.dat[,-1]) #make it a matrix, without the ID_Col
native.mat[is.na(native.mat)] <- 0 #make all NAs 0

native.mat[native.mat == 0] <- 0.025 #get rid of 0s

#add the extra days so that it becomes a daily timestep
first <- matrix(NA, nrow = 30, ncol=17)
second <- matrix(native.mat[,1])
third <- matrix(NA, nrow = 30, ncol = 6)
fourth <- matrix(native.mat[,2])
fifth <- matrix(NA, nrow = 30, ncol = 6)
sixth <- matrix(native.mat[,3])
seventh <- matrix(NA, nrow = 30, ncol = 6)
eighth <- matrix(native.mat[,4])

native.mat <- cbind(first, second, third, fourth, fifth, sixth, seventh, eighth)


##State Prediction Funcions: ####
#Model with single and multi r and process vs observation error


#predict the mean (expected value) of every timestep when you don't have an observation- see below
#Multi process error
multi.func.p<-function(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, obs,n0, species.vec){
  rvec <- c(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
  
  dims <- dim(obs)
  ntubs <- dims[1]
  ts <- dims[2]
  
  Nout <- matrix(0, nrow = ntubs, ncol = ts)
  for(i in 1:ntubs){
    Nout[i,1]<-n0*(1+rvec[species.vec[i]]*(1-n0/.995))
  }
  
  for(i in 2:ts){
    for(j in 1:ntubs){
      if(!is.na(obs[j, i-1])) {
        Nout[j, i]<-obs[j, i-1]*(1+rvec[species.vec[j]]*(1-obs[j, i-1]/.995))
      }
      if(is.na(obs[j, i-1])){ #if it is an NA, do off the last predicted
        Nout[j, i] <- Nout[j, i-1]*(1+rvec[species.vec[j]]*(1-Nout[j, i-1]/.995))
      }
    }
  }
  return(Nout)
}

#test
multi.func.p(r1 = .7, r2 = .9, r3 = .9, r4 = .7, r5 = .9, r6 = .2, r7 = .7, r8 = .6,
             r9 = .9, r10 = .1,species.vec = species.vec,
             obs = native.mat, n0 = .1)


##NLL Function ####


#Multi process error - for the process error model, I just need to predict a number based on the previous timstep
#I would calculate the timestep between recorded values similar to how we did with the spatial models
#then the mean is your last observed value
#variance is the variance of the last timestepx the number of time steps because the variance increases linearly with the number of timesteps out 
#tbefore putting it into dnorm, you calculate the sd from the variance (sqrt(timesteps x variance))
#qlogis of the predicted N 
nll.multi.func.p<-function(lr1, lr2, lr3, lr4, lr5, lr6, lr7, lr8, lr9,
                           lr10, species.vec,
                           obs,ln0, lsd){
  r1<-exp(lr1)
  r2<-exp(lr2)
  r3<-exp(lr3)
  r4<-exp(lr4)
  r5 <- exp(lr5)
  r6 <- exp(lr6)
  r7 <- exp(lr7)
  r8 <- exp(lr8)
  r9 <- exp(lr9)
  r10 <- exp(lr10)
  
  s <-exp(lsd)
  n0 <- exp(ln0)
  
  predN<-multi.func.p(r1=r1,r2=r2,r3=r3,r4=r4,r5=r5,r6=r6, r7=r7,r8=r8, r9=r9, r10=r10,
                      obs=obs,species.vec = species.vec, n0 = n0)
  
  # obs2 <- obs[!is.na(obs)]
  # predN2 <- predN[!is.na(obs)]
  
  # param <- dampack::beta_params(mean = predN2, sigma = sd)
  # alpha <- param[1]
  # alpha <- unlist(alpha)
  # beta <- param[2]
  # beta <- unlist(beta)
  
  predN[predN==0]<-.01
  predN[predN==1]<-.99
  # print(obs)
  # print(predN)
  liks<-0
  
  for(j in 1:nrow(obs)){
    lastobs <- 0
    for(i in 1:ncol(obs)){
      if(!is.na(obs[j, i])){
        tbtwn<-i-lastobs
        # print(tbtwn)
        liks<-c(liks, dnorm(x=qlogis(obs[j, i]),mean=qlogis(predN[j, i]),sd=sqrt(tbtwn*s^2)))
        lastobs<-i
        # print(liks)
      }
    }
  }
  
  
  
  # likes <- dbeta(x=obs2,shape1 = alpha, shape2 = beta)
  # likes <- dnorm(x=qlogis(obs2), mean=plogis(predN2), sd = sd)
  
  nll<--1*sum(log(liks[-1]))  
  return(nll)
}

#test
nll.multi.func.p(lr1 = log(.4), lr2 = log(.4),lr3 = log(.4),
                 lr4 = log(.4),lr5 = log(.4),lr6 = log(.4),
                 lr7 = log(.4),lr8 = log(.4),lr9 = log(.4), lr10 = log(.4),
                 ln0 = log(.01), obs = native.mat, lsd = log(.05), species.vec = species.vec)



##Find MLE parameters####


library(bbmle)

# Create list of starting guesses for the multi model
start.list<-list(lr1=log(.1),
                 lr2 = log(.2),
                 lr3 = log(.1),
                 lr4 = log(.1),
                 lr5=log(.2),
                 lr6 = log(.1),
                 lr7 = log(.1),
                 lr8 = log(.1),
                 lr9 = log(.1),
                 lr10 = log(.1),
                 lsd = log(.05),
                 ln0 = log(0.001))
# Create list of observed data for model
data.list<-list(obs=native.mat, species.vec = species.vec, method = "SANN")

# # Find MLE parameter estimates
fit_mp<-mle2(minuslogl=nll.multi.func.p,start=start.list,data=data.list, method = "SANN")
# # store MLE parameter estimates
cfs_mp<-exp(coef(fit_mp)) 
cfs_mp


## Predict historical dynamics from MLE parameter estimates ####


pred_mp<-multi.func.p(r1 = cfs_mp[1], r2 = cfs_mp[2], r3 = cfs_mp[3],
                      r4 = cfs_mp[4], r5 = cfs_mp[5], r6 = cfs_mp[6],
                      r7 = cfs_mp[7], r8 = cfs_mp[8], r9 = cfs_mp[9],
                      r10 = cfs_mp[10], n0 = cfs_mp[11], obs = native.mat, species.vec = species.vec)
plot(native.mat, pred_mp, xlab = "Observed", ylab = "Predicted",pch = 16, las = 1, ylim = c(0,1))
abline(0, 1) # Add 1:1 line on figure indicating perfect fit
summary(lm(as.vector(native.mat) ~ as.vector(pred_mp)))
#mtext("R2 = 0.88", side=3)


