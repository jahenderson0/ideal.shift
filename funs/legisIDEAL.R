############################################
# John Henderson
# - legisIDEAL: ideal point estimator incorporating issue-by-issue shifts
# - June 5, 2018
############################################

######################################################################
# legisIDEAL -- produces supervised IRT estimates of legislative choices, w/ issue/topic-specific ideological shifts
#    cutpoint ::
# 		== TRUE is cutpoint implementation :: alpha*beta
# 		== FALSE is threshold distance implementation :: (alpha-beta)^2
#
#    party :: length of MCs being scaled
# 		== -1 for Democrats or Liberal
# 		== +1 for Republicans or Conservative
#
#    phi :: matrix of bill topic proportions = K x T
#
#    y :: matrix of MC x bill cosponsor choices = N x K
#
######################################################################

legisIDEAL=function(y,party=party,eps=.1,cutpoint=FALSE,phi=phi){ # better=FALSE; only sweep start ...

	control=list('alpha'=FALSE,'beta'=FALSE,'delta'=FALSE,'gamma'=FALSE)
  alpha=NULL;gamma=NULL;beta=NULL;delta=NULL;

	# unsupervised IRT
	if(is.null(c(alpha,gamma,beta,delta))){
		cat("##########################\n  Running fully unsupervised IRT model to estimate alpha, gamma, beta, delta and z-shift parameters using\n    --",nrow(y),"MCs and",ncol(y),"bills")
	 	control$alpha=TRUE
		control$beta= TRUE
		control$delta=TRUE
		control$gamma=TRUE
	}

	if(cutpoint==TRUE){
		cat("\n    -- with cutpoints")
	} else if(cutpoint==FALSE){
		cat("\n    -- with distances")
	}
	cat("\n    -- with convergence threshold",eps,"\n")
	cat("##########################\n")
	cat("Iteration   Convergence\n")
  # end print frontmatter

	if(cutpoint==TRUE){
		# liklihood function in dimension i=1,2,...,n
		llik.i=function(pars,y,beta,delta,phi,lam_a,lam_g,lam_z,methods=0){
			alpha=pars[1]
			gamma=pars[2]
      zi=pars[-c(1,2)]
      zphi=phi%*%zi
			if(methods==0){
				mu=(1/(1+exp(- ((((alpha+zphi)*beta))+gamma+delta))))
				llk=(sum(y*log(mu)+(1-y)*log(1-mu))
					-lam_a*sum(alpha^2)
					-lam_g*sum(gamma^2)
          -lam_z*sum(zi^2)
				)
			}
			return(-llk)
		}

    # liklihood functions in dimension j=1,2,...,m
		llik.j=function(pars,y,alpha,gamma,zphi,lam_b,lam_d,methods=0){
			beta=pars[1]
			delta=pars[2]
			if(methods==0){
				mu=(1/(1+exp(- ((((alpha+zphi)*beta))+gamma+delta))))
				llk=(sum(y*log(mu)+(1-y)*log(1-mu))
				 	-lam_b*sum(beta^2)
				 	-lam_d*sum(delta^2)
				)
			}
			return(-llk)
		}

    # sum liklihood function
    llik.sum=function(y,alpha,gamma,delta,beta,zphi,n,m,z){
      lsum=0
      for(i in 1:n){
	      mu=(1/(1+exp(- ((((alpha[i]+zphi[i,])*beta))+gamma[i]+delta))))
        lsum=lsum+sum(y[i,]*log(mu)+(1-y[i,])*log(1-mu))
      }
      for(j in 1:m){
        mu=(1/(1+exp(- ((((alpha+zphi[,j])*beta[j]))+gamma+delta[j]))))
        lsum=lsum+sum(y[,j]*log(mu)+(1-y[,j])*log(1-mu))
      }

      lsum=lsum/2-lam_b*sum(beta^2)-lam_d*sum(delta^2)-lam_a*sum(alpha^2)-lam_g*sum(gamma^2)-lam_z*sum(z[i,]^2)
      return(lsum)
    }

	} else if(cutpoint==FALSE){
    # liklihood function in dimension i=1,2,...,n
		llik.i=function(pars,y,beta,delta,phi,lam_a,lam_g,lam_z,methods=0){
			alpha=pars[1]
			gamma=pars[2]
      zi=pars[-c(1,2)]
      zphi=phi%*%zi
			if(methods==0){
				mu=(1/(1+exp(- ((-((alpha+zphi)-beta)^2)-gamma-delta))))
				llk=(sum(y*log(mu)+(1-y)*log(1-mu))
					-lam_a*sum(alpha^2)
					-lam_g*sum(gamma^2)
          -lam_z*sum(zi^2)
				)
			}
			return(-llk)
		}

		# liklihood functions in dimension j=1,2,...,m
		llik.j=function(pars,y,alpha,gamma,zphi,lam_b,lam_d,methods=0){
			beta=pars[1]
			delta=pars[2]
			if(methods==0){
				mu=(1/(1+exp(- ((-((alpha+zphi)-beta)^2)-gamma-delta) )))
				llk=(sum(y*log(mu)+(1-y)*log(1-mu))
			 		-lam_b*sum(beta^2)
			 		-lam_d*sum(delta^2)
				)
			}
			return(-llk)
		}

    # sum liklihood function
		llik.sum=function(y,alpha,gamma,delta,beta,zphi,n,m,z){
      lsum=0
      for(i in 1:n){
        mu=(1/(1+exp(- ((-((alpha[i]+zphi[i,])-beta)^2)-gamma[i]-delta))))
        lsum=lsum+sum(y[i,]*log(mu)+(1-y[i,])*log(1-mu))
      }
      for(j in 1:m){
        mu=(1/(1+exp(- ((-((alpha+zphi[,j])-beta[j])^2)-gamma-delta[j]))))
        lsum=lsum+sum(y[,j]*log(mu)+(1-y[,j])*log(1-mu))
      }

		  lsum=lsum/2-lam_b*sum(beta^2)-lam_d*sum(delta^2)-lam_a*sum(alpha^2)-lam_g*sum(gamma^2)-lam_z*sum(z[i,]^2)
			return(lsum)
		}
	}

	# function to create starting values sweeping over naive model estimates
	sweepStart=function(y=y,rest=F,party=party,alpha=alpha,gamma=gamma,beta=beta,delta=delta){
    N=nrow(y)
    M=ncol(y)

		if(is.null(beta)){
			# party to be in the direction of -1 for Dems or Left party; +1 for Reps or Right party
			beta=foreach(j = 1:M,.combine=c) %dopar% {
				mean(party[which(y[,j]>0)])/(1/sum(y[,j]))
			}
      beta[which(is.na(beta))]=median(beta,na.rm=T)
		}

		if(is.null(alpha)){
			alpha=foreach(i = 1:N,.combine=c) %dopar% {
				mean(beta[which(y[i,]>0)])
			}
			alpha[which(is.na(alpha))]=median(alpha,na.rm=T)
		}

		# converge the sweep
		set.seed(1005)
		alpha=(alpha-mean(beta))/sd(beta)
		beta=jitter((beta-mean(beta))/sd(beta),factor=2)  # small jitter to climb in probability

		difs=100
		while(difs>1){
			beta=foreach(j = 1:M,.combine=c) %dopar% {
				mean(alpha[which(y[,j]>0)])
			}

			beta[is.na(beta)]=mean(beta[!is.na(beta)])
			beta=jitter((beta-mean(beta))/sd(beta),factor=2)

			alpha=foreach(i = 1:N,.combine=c) %dopar% {
				mean(beta[which(y[i,]>0)])
			}
			alpha[is.na(alpha)]=mean(alpha[!is.na(alpha)])
			alpha_old=alpha

			difs=sum((alpha_old-alpha)^2,na.rm=T)
		}

		if(rest==T){
			if(is.null(gamma)){
				gamma=-(rowMeans(y)/rowVars(y))
				gamma[which(is.na(gamma))]=mean(gamma[which(!is.na(gamma))])
				gamma=(gamma-mean(gamma))/sd(gamma)
			}
			if(is.null(delta)){
				delta=-(colMeans(y)/colVars(y))
				delta=(delta-mean(delta))/sd(delta)
			}
		} else{
			set.seed(1005)
			if(is.null(gamma)){
				gamma=rnorm(n=nrow(y))
			}
			if(is.null(delta)){
				delta=rnorm(n=ncol(y))
			}
		}

		return(list("alpha"=alpha,"beta"=beta,"gamma"=gamma,"delta"=delta))
	}


  ##########################################
	# 0. initialize starting values
  ##########################################

  # set up parallel cores
  n.cores = parallel::detectCores() - 1
  my.cluster = parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)

	n=nrow(y); m=ncol(y)
  # normal priors on spatial terms; floating/diffuse priors on valence terms
	lam_b=lam_a=.5
	lam_d=lam_g=.5
	lam_z=.5

  ##########################################
  # 1. initialize optimization
  ##########################################

  z=matrix(rnorm(0,1,n=nrow(y)*ncol(phi)),nrow(y),ncol(phi))
	starts=sweepStart(y=y,rest=F,party=party,alpha=alpha,gamma=gamma,beta=beta,delta=delta)
	alpha=starts$alpha
	beta=starts$beta
	gamma=starts$gamma
	delta=starts$delta

  ##########################################
  # 2. running optimization while loop until convergence
  ##########################################

  llk0=100
  llk1=1000
  difs=llk1-llk0
  counts=0
  while(difs>eps){
    if(counts>0){
      cat(counts,"   ",difs,"\n")
    }

    alpha_old=alpha
    gamma_old=gamma
    beta_old=beta
    delta_old=delta

    counts=counts+1

    pars1=foreach(i = 1:n,.combine=cbind) %dopar% {
      pars=c(alpha[i],gamma[i],z[i,])
      optim(pars,llik.i,y=y[i,],beta=beta,delta=delta,phi=phi,lam_z=lam_z,lam_a=lam_a,lam_g=lam_g,methods=0,method='BFGS')$par
    }
    alpha=pars1[1,]
    gamma=pars1[2,]
    z=t(pars1[-c(1,2),])

    if(!is.null(party)){
      if(mean(na.rm=T,alpha[which(party==-1)])>mean(na.rm=T,alpha[which(party==1)])){
        alpha=alpha*-1
      }
    }

    zphi=z%*%t(phi)
    pars1=foreach(j = 1:m,.combine=cbind) %dopar% {
      pars=c(beta[j],delta[j])
      optim(pars,llik.j,y=y[,j],alpha=alpha,gamma=gamma,zphi=zphi[,j],lam_b=lam_b,lam_d=lam_d,methods=0,method='BFGS')$par
    }
    beta=pars1[1,]
    delta=pars1[2,]

    #difs=sum(c((alpha_old-alpha)^2))

    llk1=llik.sum(y,alpha,gamma,delta,beta,zphi,n,m,z)
    difs=abs(llk1-llk0)
    llk0=llk1
    if(counts>50){
      warning("Did not fully converge, breaking at ",difs)
      break
    }
  }

  ##########################################
  # 3. return estimates
  ##########################################
  parallel::stopCluster(cl = my.cluster)

	return(list('alpha'=alpha,'beta'=beta,'delta'=delta,'gamma'=gamma,"z"=z,'llk.converged'=difs))
}

#END legisIDEAL
