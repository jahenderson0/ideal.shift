############################################
# John Henderson
# - scaling congress issue-by-issue
# - explatory analysis
# - April 18, 2019 [updated May 11, 2022]
############################################
#
# GOAL: produce a model of ideal points that can shift isssue-by-issue; reflects variation in how MCs represent districts
#   depending on the issue areas at stake
#
# APPROACH: scale cosponsor choices with top-coded issue areas using the Policy Agendas Project codings;
#   see here: http://congressionalbills.org/index.html
#    - to replicate topic mixtures, the idea is use words in bill titles to predict PAP issue areas, this gives topic proportion
#   estimates for each bill
#   - these topic proportions are used as weights in an unsupervised IRT model, so that ideal points are allowed to float left or right
#   on cosponsorship choices that have more topic weight given to topic
#   - the cutpoint model then is Pr(Cosponsor=1) = f{(alpha+t(z)%*%phi)*beta + gamma + delta};
#   - where alpha are ideal points over all choices; z reflect topic-specific shifts in ideal points
#   - beta are bill discrimination parameters; delta are bill location parameters; gamma are legislator connectedness paramaters
#   - a threshold model is also explored: Pr(Cosponsor=1) = f{(alpha+t(z)%*%phi - beta)^2 + gamma + delta};

rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load('glmnet','tm','stringr','doParallel','parallel','Rcpp')
pacman::p_load('RcppArmadillo','inline','RcppHNSW','emIRT','cmdstanr','posterior')

#require(glmnet)
#require(tm)
#require(stringr)
#require(doParallel)
#require(parallel)
#require(Rcpp)
#require(RcppArmadillo)
#require(inline)
#require(RcppHNSW)
#require(emIRT)
#require(cmdstanr)
#require(posterior)

setwd('~/Documents/CV/NonAC/yougov/datascience')

############################################
# source all the functions
############################################
source('ideal.shift/funs/legisIDEAL.R')
source('ideal.shift/funs/smoteFun.R')
source('ideal.shift/funs/cosineDistF.R')

############################################
# load the cosponsorship data :: to predict PAP topic
############################################

load("ideal.shift/data/houseTitlesCosponsors.Rdata")
#congress=110
#titles=titleBills[[congress]]
#cosponsor=cosponsorBills[[congress]]

dim(cosponsor)
#[1]  449 7347

dim(titles)
#[1] 7340   12

txts=titles[,3]
txts=gsub(txts,pattern="\\W",replace=' ')
txts=gsub(txts,pattern="[0-9]",replace=' ')
#txts=stemDocument(txts)
txts=tolower(txts)
txts=removeWords(txts,stopwords("english"))
txts=stripWhitespace(txts)
txts=str_squish(txts)
tdm=DocumentTermMatrix(VCorpus(VectorSource(txts)))
tdm=as.matrix(tdm)
y=titles$Major
y[which(y==99)]=0

############################################
# focus where data has meaningful density
############################################
iDs=cosponsor[-c(1),1:7]
cosD=cosponsor[-c(1),-c(1:7)]
for(j in 1:ncol(cosD)){
  cosD[which(is.na(cosD[,j])),j]=0
}
iDs=iDs[which(rowSums(cosD)>25),]
cosD=cosD[which(rowSums(cosD)>25),]

cx=colSums(cosD)>15 & y>0

y=y[which(cx)]
tdm=tdm[which(cx),]
cosD=cosD[,which(cx)]

party=as.numeric(iDs$party==200)
party[which(party==0)]=-1

############################################
# outcome indicator for PAP topics
############################################

lvls=sort(unique(y))
len_lvls=length(lvls)
ymat=matrix(0,length(y),len_lvls)
for(i in 1:len_lvls){
  ymat[,i]=as.numeric(y==lvls[i])
}
Nk=ncol(ymat)

############################################
# smote => ridge on full data
############################################
# loop is pretty fast all-in-all

############################################
# do ridge in parallel loops || -- decent timing
############################################
# set tuning range :: pre-tested to fit better w/ less shrinkage generally
lambda=seq(from=.5,to=0.025,by=-0.025)
#alpha=0#seq(from=0,to=1,by=0.1)

# set up parallel cores
n.cores = parallel::detectCores() - 1
my.cluster = parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
#source('ideal.shift/funs/smoteFun.R')

preds=foreach(i = 1:Nk, .packages=c("glmnet","RcppHNSW")) %dopar% {
  dts=smoteFun(y=ymat[,i],x=tdm,k=20)
  xtmp=dts$x
  ytmp=dts$y

  ############################################
  # set training and testing data on augmented data
  ############################################
  set.seed(1005)
  ix.tr=sort(sample(which(!is.na(ytmp)),size=nrow(xtmp)/2,replace=F))

  x.train=as.matrix(xtmp)[ix.tr,]
  x.test=as.matrix(xtmp)[-c(ix.tr),]

  y.test=ytmp[-c(ix.tr)]
  y.train=ytmp[ix.tr]
  ytr=(y.train-.5)/.5
  oo=cv.glmnet(nfolds=10,x=x.train, y=ytr, family=c("gaussian"),nlambda=length(lambda),alpha=0,lambda=lambda)
  lam_min=oo$lambda.min
  op=glmnet(x=x.train, y=ytr, family=c("gaussian"),nlambda=length(lambda),alpha=0,lambda=lam_min)

  pred=as.numeric(predict(op,x.test)>0)
  #pred_tm=as.numeric(predict(op,tdm)>.5)
  pred_tm=1/(1+exp(-predict(op,tdm)))
  list(pred_tm,table(pred,y.test),op,lam_min)
}

parallel::stopCluster(cl = my.cluster)

### analogy to topic proportions ###
beta_mat=matrix(NA,Nk,ncol(tdm))
pred_mat=matrix(NA,Nk,nrow(tdm))

for(jj in 1:Nk){
  beta_mat[jj,]=preds[[jj]][[3]]$beta[,1]
  pred_mat[jj,]=preds[[jj]][[1]][,1]
}
topic_prop=matrix(NA,nrow(tdm),Nk)
for(i in 1:nrow(tdm)){
  Rs=rowSums(beta_mat[,which(tdm[i,]>0)]>0)
  topic_prop[i,]=Rs/sum(Rs)
}

colnames(beta_mat)=colnames(tdm)

############################################
# plotting to validate models by coarse inspection
############################################

#table(y)
 #1   2   3   4   5   6   7   8  10  12  13  14  15  16  17  18  19  20  21
#93  44 374  27 104 139  77 100  38 117  44  39 104 170  37  24  87 217  64

topcodes=c(
  '1. Macroeconomics','2. Civil Rights, Minority Issues, and Civil Liberties',
  '3. Health','4. Agriculture','5. Labor and Employment',
  '6. Education','7. Environment','8. Energy','10. Transportation',
  '12. Law, Crime and Family Issues','13. Social Welfare',
  '14. Community Development and Housing Issues','15. Banking, Finance, and Domestic Commerce',
  '16. Defense','17. Space, Science, Technology, and Communications',
  '18. Foreign Trade','19. International Affairs and Foreign Aid',
  '20. Government Operations','21. Public Lands and Water Management'
)

for(kk in 1:nrow(beta_mat)){
  mm=sort(beta_mat[kk,],decreasing=F)
  mm=c(mm[(length(mm)-50):length(mm)])

  pdf(height=5,width=6,file=paste('ideal.shift/figures/',gsub(topcodes[kk],pattern=' ',replace='_'),'_pap.pdf',sep=''))
  plot(cex=.8,x=mm,y=1:length(mm),col='white',xlim=c(-0,1.2),ylim=c(-6,length(mm))+6,xlab='Coefficient',axes=F,ylab='Top-50 Predictors')
  box(col='grey')
  text(x=mm,y=1:length(mm),label=names(mm),cex=log(1+7*abs(mm))/2,pos=4)

  points(x=c(mm),y=c(1:length(mm)),col='red',lty=2,cex=.25)
  axis(1,)
  legend('topleft',legend=topcodes[kk],cex=.75)
  dev.off()
}

###### remove some structure from cosD that slow things down
y=as.matrix(cosD)

###### define phi topic proportions for models below
phi=topic_prop

############################################
# legisIDEAL -- produces supervised IRT estimates of legislative choices, w/ issue/topic-specific ideological shifts
############################################
#source('ideal.shift/funs/legisIDEAL.R')

mod_outs=legisIDEAL(y,party=party,eps=.1,cutpoint=FALSE,phi=phi)

############################################
# fast IRT ideal point model for checking (appproximate) model fit
############################################
#library(emIRT)

rc = convertRC(rollcall(cosD,yea=1,nay=-1,missing=0))
p = makePriors(rc$n, rc$m, 1)
s  = getStarts(rc$n, rc$m, 1)

outs=binIRT(.rc=rc,
  .starts = s,
  .priors = p,
  .control = {
    list(threads = 1,
    verbose = FALSE,
    thresh = 1e-6)
  },
  .anchor_subject = 1
)

alpha0=outs$means[[1]][,1]
beta0=outs$means[[2]][,2]
delta0=outs$means[[2]][,1]

############################################
# cmdscale for checking (appproximate) model fit
############################################
#source('ideal.shift/funs/cosineDistF.R')

dS=cosineDistF(y)[[1]]
cmd_alpha=cmdscale(dS,k=1)
if(mean(cmd_alpha[which(party==-1)])>mean(cmd_alpha[which(party==1)])){
  cmd_alpha=-cmd_alpha
}

############################################
############################################
# STAN IDEAL POINT MODEL w/ SHIFT :: alternative way to estimate the shift model
############################################
############################################
#library(cmdstanr)
#library(posterior)

############################################
# shift model estimating unsupervised IRT parameters as data
############################################

scode <- "data {
    int<lower=2> T;              // num topics
    int<lower=1> J;              // number of legislators
    int<lower=1> K;              // number of bills
    int<lower=1> N;              // number of bills-legislator pairs
    array[N] int<lower=1,upper=J> jj;  // legislator for observation n
    array[N] int<lower=1,upper=K> kk;  // bill for observation n
    array[N] int<lower=0,upper=1> y;   // cosponsored for observation n
    matrix[T,K] phi;  	       	 // topic bill probabilities bill x topic
}
parameters {
    real mu_beta;                // mean bill discrimination
    real mu_delta;               // mean bill difficulty
    real mu_gamma;               // mean legislator connectedness
    vector[J] alpha;  	       	 // legislator location
    vector[J] gamma;  	       	 // legislator connectedness
    vector[K] beta;              // bill discrimination
    vector[K] delta;             // bill difficulty
    real<lower=0> sigma_beta;    // scale of discrimination
    real<lower=0> sigma_delta;   // scale of difficulty
    real<lower=0> sigma_gamma;   // scale of connectedness
    real mu_z;									 // mean document x topic shift
    matrix[J,T] z;							 // document shift on T topic
    real<lower=0> sigma_z;			 // scale of topic location shifts
}
model {
    to_vector(z) ~ normal(mu_z, sigma_z);
    mu_z ~ cauchy(0, 5);
    sigma_z ~ cauchy(0, 5);
    alpha ~ std_normal();
    beta ~ normal(mu_beta, sigma_beta);
    delta ~ normal(mu_delta, sigma_delta);
    gamma ~ normal(mu_gamma, sigma_gamma);
    mu_beta ~ cauchy(0, 5);
    sigma_beta ~ cauchy(0, 5);
    mu_delta ~ cauchy(0, 5);
    sigma_delta ~ cauchy(0, 5);
    mu_gamma ~ cauchy(0, 5);
    sigma_gamma ~ cauchy(0, 5);
    {
      matrix[J,K] zphi = z * phi;
      vector[N] a_z;
      for (n in 1:N) {
        a_z[n] = (alpha[jj[n]] + zphi[jj[n],kk[n]]);
      }
      y ~ bernoulli_logit(beta[kk] .* a_z + delta[kk] + gamma[jj]);
    }
}
"

writeLines(scode, "~/.cmdstan/cmdstan-2.29.1/examples/shift/ideal.shift.stan")

# cmdstan must be installed and linked using correct path to makefile
set_cmdstan_path("~/.cmdstan/cmdstan-2.29.1")
file = file.path(cmdstan_path(), "examples", "shift", "ideal.shift.stan")
mod = cmdstan_model(file)

############################################
# vectorizes data for use in stan models
############################################

vecFun=function(mat,w=NULL){
	if(is.null(w)){
		w=rep(1,ncol(mat))
	}
	jj=rep(1,ncol(mat))
	for(i in 2:nrow(mat)){
	  jj=c(jj,rep(i,ncol(mat)))
	}
	kk=rep(1:ncol(mat),nrow(mat))
	mat=mat*w
	yy=c(t(mat))
	return(list("yy"=yy,"jj"=jj,"kk"=kk))
}

vout=vecFun(as.matrix(cosD))
jj=vout$jj
kk=vout$kk
yy=vout$yy
J=nrow(cosD)
K=ncol(cosD)
N=J*K

list2 = list("jj"=jj,"kk"=kk,"y"=yy,"J"=J,"K"=K,"N"=N,"T"=ncol(topic_prop),"phi"=t(topic_prop))

do.stan=FALSE # can take a very long time to get good chain convergence (7 hrs+!)
if(do.stan==TRUE){
  fit = mod$sample(
    data = list2,iter_sampling=250,iter_warmup=250,
    seed = 1005,
    chains = 4,
    parallel_chains = 4,
    refresh = 5 # print update every 500 iters
  )
  save(fit,file='ideal.shift/stan/fit.Rdata')
  outs1=summarise_draws(fit)
  outs1=as.matrix(outs1)
  alpha1=as.numeric(outs1[,2][grep(outs1[,1],pattern='alpha\\[')])
  gamma1=as.numeric(outs1[,2][grep(outs1[,1],pattern='gamma\\[')])
  beta1=as.numeric(outs1[,2][grep(outs1[,1],pattern='beta\\[')])
  delta1=as.numeric(outs1[,2][grep(outs1[,1],pattern='delta\\[')])
  z=as.numeric(outs1[,2][grep(outs1[,1],pattern='z\\[')])

  if(mean(alpha1[which(party==-1)])>mean(alpha1[which(party==1)])){
    alpha1=-alpha1
  }
  save(alpha1,file='ideal.shift/stan/stan.alpha1.Rdata')
} else {
  load(file='ideal.shift/stan/fit.Rdata')
  load('ideal.shift/stan/stan.alpha1.Rdata')
}



############################################
# plot all the alpha estimates
############################################
alphas=cbind(mod_outs$alpha/sd(mod_outs$alpha),alpha0/sd(alpha0),cmd_alpha/sd(cmd_alpha),alpha1/sd(alpha1))
nms=c('legisIDEAL','emIRT','cmdscale','stan')


pdf(file='ideal.shift/figures/alpha_cor.pdf')
vv=c()
for(i in 1:4){
  for(j in 1:4){
    if(i < j){
      v=round(cor.test(alphas[,i],alphas[,j])$est,digit=2)
      vv=c(vv,v)
    }
  }
}


par(mfrow = c(3,3), oma = c(2,3,2,0) + 0.8, mar = c(0,0,0,0) + 0.0)
for(i in 1:4){
  for(j in 1:4){

    if(i==4){
      break()
    }

    if(i==j){

    } else if(i>j){
        plot(x=-1000,y=-1000,xlim=c(-2.25,2.25),ylim=c(-2.25,2.25),ylab='',xlab='',axes=F)
    } else if(i<j){
      plot(x=alphas[,i],y=alphas[,j],xlim=c(-2.25,2.25),ylim=c(-2.25,2.25),ylab='',xlab='',axes=F)
      box(col='grey',lty=2)
      v=round(cor.test(alphas[,i],alphas[,j])$est,digit=2)
      text(x=-1.75,y=2.,label=paste('p = ',v,sep=''),col="black",cex=0.95)
      lines(lowess(x=alphas[which(!is.na(alphas[,i]) & !is.na(alphas[,j])),i],y=alphas[which(!is.na(alphas[,i]) & !is.na(alphas[,j])),j]),col="#FF000099",lwd=2)
      if(j-i==1){
          mtext(text=nms[i],side=2,line=1.,cex=.95,col='red')
      }
      if(i==1){
        mtext(text=nms[j],side=3,line=1.,cex=.95,col='red')
      }
    }
  }
}
dev.off()

############################################
# plot the shifts from the mod_outs model
############################################

zshift=mod_outs$z
zalpha=mod_outs$alpha
zshift=zshift[order(zalpha),]
zparty=party[order(zalpha)]
zalpha=sort(zalpha)


for(kk in 1:nrow(beta_mat)){
  pdf(height=3,width=6,file=paste('ideal.shift/figures/',gsub(topcodes[kk],pattern=' ',replace='_'),'_shifts.pdf',sep=''))
  plot(y=rep(0,length(zalpha)),x=zalpha,cex=.5,col='white',axes=F,ylab=topcodes[kk],xlab='alpha + shift estimates',ylim=c(0,1),xlim=c(-5,5))
  axis(1,)
  points(y=rep(0,length(zalpha[which(zparty==-1)])),x=zalpha[which(zparty==-1)],cex=.5,col='blue')
  points(y=rep(0,length(zalpha[which(zparty==1)])),x=zalpha[which(zparty==1)],cex=.5,col='red')
  points(y=rep(1,length((zshift[,kk]+zalpha)[which(zparty==-1)])),x=(zshift[,kk]+zalpha)[which(zparty==-1)],cex=.5,col='blue')
  points(y=rep(1,length((zshift[,kk]+zalpha)[which(zparty==1)])),x=(zshift[,kk]+zalpha)[which(zparty==1)],cex=.5,col='red')
  for(i in 1:length(zalpha)){
    if(zparty[i]==-1){
      lines(x=c(zalpha[i],(zshift[,kk]+zalpha)[i]),y=c(0,1),col='blue')
    } else {
      lines(x=c(zalpha[i],(zshift[,kk]+zalpha)[i]),y=c(0,1),col='red')
    }
  }
  dev.off()
}

#END shift.ideal.cosponsorship.R
