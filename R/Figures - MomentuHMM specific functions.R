## Assign HMM state to dive batches and test covariates effect 
##Author: Brett McClintock - https://github.com/bmcclintock/momentuHMM
##Purpose: This code contains some functions developed in the momentuHMM package that needs to be loaded in the environment
#          to recreate the stationary and transition probabilities figures.

#Reference
#McClintock, B. T. and T. Michelot (2018). "momentuHMM : R package for generalized hidden Markov models of animal movement." Methods in Ecology and Evolution 9(6): 1518-1530.

#FORMATRECHARGE FUNCTION
formatRecharge <- function(nbStates,formula,data,covs=NULL,par=NULL){
  
  nbAnimals <- length(unique(data$ID))
  dataNames <- colnames(data)
  
  newForm <- newFormulas(formula,nbStates)
  newformula <- newForm$newformula
  formulaStates <- newForm$formulaStates
  recharge <- hierRecharge <- newForm$recharge
  
  aInd <- NULL
  for(i in 1:nbAnimals){
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  }
  
  if(!is.null(recharge)){
    
    recharge <- expandRechargeFormulas(hierRecharge)
    
    if(inherits(data,"hierarchical")){
      
      recLevels <- length(hierRecharge)
      recLevelNames <- names(hierRecharge)
      
      tmpg0covs <- model.matrix(recharge$g0,data)
      g0covs <- tmpg0covs[rep(aInd,each=recLevels)+rep(0:(recLevels-1),nbAnimals),,drop=FALSE]
      for(i in 1:nbAnimals){
        for(iLevel in 1:recLevels){
          g0covs[(i-1)*recLevels+iLevel,] <- tmpg0covs[which(data$level==iLevel & data$ID==unique(data$ID)[i])[1],,drop=FALSE]
        }
      }
      recovs <- model.matrix(recharge$theta,data)
      rechargeNames <- paste0("recharge",gsub("level","",recLevelNames))
      colInd <- lapply(recLevelNames,function(x) which(grepl(paste0("I((level == \"",gsub("level","",x),"\")"),colnames(recovs),fixed=TRUE)))
    } else {
      g0covs <- model.matrix(recharge$g0,data[aInd,])
      recovs <- model.matrix(recharge$theta,data)
      recLevels <- 1
      rechargeNames <- "recharge"
      colInd <- list(1:ncol(recovs))
    }
    
    nbG0covs <- ncol(g0covs)-1
    nbRecovs <- ncol(recovs)-1
    data[rechargeNames]<- list(rep(0,nrow(data)))
    
    if(is.null(par)) par <- list(g0=rep(0,nbG0covs+1),theta=rep(0,nbRecovs+1))
    for(i in 1:nbAnimals){
      for(iLevel in 1:recLevels){
        idInd <- which(data$ID==unique(data$ID)[i])
        if(nbRecovs){
          if(!all(names(par$theta)==colnames(recovs)) | !all(names(par$g0)==colnames(g0covs))) stop("column name mismatch in hierarchical recharge model -- please report to brett.mcclintock@noaa.gov")
          g0 <- par$g0 %*% t(g0covs[(i-1)*recLevels+iLevel,,drop=FALSE])
          theta <- par$theta
          data[[rechargeNames[iLevel]]][idInd] <- cumsum(c(g0,theta[colInd[[iLevel]]]%*%t(recovs[idInd[-length(idInd)],colInd[[iLevel]]])))
        }
      }
    }
    #if(is.null(hierRecharge)) newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    if(!is.null(covs)){
      covs[rechargeNames] <- lapply(data[rechargeNames],mean)
    }
  } else {
    nbG0covs <- 0
    nbRecovs <- 0
    g0covs <- NULL
    recovs <- NULL
  }
  
  tmpCovs <- model.matrix(newformula,data)
  if(is.null(covs)) {
    covs <- tmpCovs
  }
  nbCovs <- ncol(tmpCovs)-1 # substract intercept column
  
  list(newdata=data[,which(!(names(data) %in% dataNames)),drop=FALSE],newformula=newformula,formulaStates=formulaStates,recharge=recharge,hierRecharge=hierRecharge,covs=covs,nbCovs=nbCovs,nbG0covs=nbG0covs,nbRecovs=nbRecovs,g0covs=g0covs,recovs=recovs,aInd=aInd)
}

#GETSPLINEFORMULA
getSplineFormula <- function(formula,data,covs){
  newcovs<-covs
  splineInd<-FALSE
  splineCovs<-character()
  newformula<-character()
  Terms<-terms(formula,specials=splineList)
  factors<-attr(Terms,"factors")
  specials<-rownames(factors)[unlist(attr(Terms,"specials"))]
  for(k in rownames(factors)){
    if(k %in% specials){
      splineInd<-TRUE
      splineCovs<-unique(c(splineCovs,all.vars(as.formula(paste0("~",k)))))
      splineExpr<-qdapRegex::rm_between(k, "(", ",", extract=TRUE)[[1]]
      sp<-eval(substitute(eval(parse(text=k))),data,parent.frame())
      tmpcovs<-predict(sp,eval(substitute(eval(parse(text=splineExpr))),covs,parent.frame()))
      tmp<-colnames(model.matrix(as.formula(paste0("~",k)),data)[,-1])
      tmp<-gsub("[()]","",tmp)
      tmp<-gsub(" ","_",tmp)
      tmp<-gsub(",","_",tmp)
      tmp<-gsub("=","_",tmp)
      colnames(tmpcovs)<-tmp
      newcovs<-cbind(newcovs,tmpcovs)
      for(l in colnames(factors)){
        if(factors[k,l]){
          tmp<-colnames(model.matrix(as.formula(paste0("~",l)),data))[-1]
          tmp<-gsub("[()]","",tmp)
          tmp<-gsub(" ","_",tmp)
          tmp<-gsub(",","_",tmp)
          tmp<-gsub("=","_",tmp)
          newformula<-c(newformula,tmp)
        }
      }
    } else {
      if(length(specials)) {
        tmpspec<-specials
        tmpspec<-gsub("(","\\(",tmpspec,fixed=TRUE)
        tmpspec<-gsub(")","\\)",tmpspec,fixed=TRUE)
        lfact <- grep(paste(tmpspec,collapse="|"),colnames(factors), value=TRUE,invert=TRUE)
      } else lfact <- colnames(factors)
      for(l in lfact){
        if(factors[k,l]){
          newformula<-c(newformula,l)
        }
      }
    }
  }
  if(!splineInd) newformula <- formula
  else newformula <- as.formula(paste0("~",paste0(unique(newformula),collapse="+")))
  
  return(list(formula=newformula,covs=newcovs))
}

#TRMATRIX_RCPP
trMatrix_rcpp <- function(nbStates, beta, covs, betaRef) {
  .Call('_momentuHMM_trMatrix_rcpp', PACKAGE = 'momentuHMM', nbStates, beta, covs, betaRef)
}

#GET_GAMMA
get_gamma <- function(beta,covs,nbStates,i,j,betaRef,betaCons,workBounds=NULL,mixture=1){
  tmpBeta <- rep(NA,length(betaCons))
  tmpBeta[unique(c(betaCons))] <- beta
  beta <- w2wn(matrix(tmpBeta[betaCons],nrow(betaCons),ncol(betaCons)),workBounds)
  gamma <- trMatrix_rcpp(nbStates,beta[(mixture-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],covs,betaRef)[,,1]
  gamma[i,j]
}

#W2WN
w2wn <- function(wpar,workBounds,k=0){
  
  ind1<-which(is.finite(workBounds[,1]) & is.infinite(workBounds[,2]))
  ind2<-which(is.finite(workBounds[,1]) & is.finite(workBounds[,2]))
  ind3<-which(is.infinite(workBounds[,1]) & is.finite(workBounds[,2]))
  
  wpar[ind1] <- exp(wpar[ind1])+workBounds[ind1,1]
  wpar[ind2] <- (workBounds[ind2,2]-workBounds[ind2,1]) * stats::plogis(wpar[ind2])+workBounds[ind2,1]
  wpar[ind3] <- -(exp(-wpar[ind3]) - workBounds[ind3,2])
  
  if(k) wpar <- wpar[k]
  return(wpar)
}

#NEWFORMULAS
newFormulas <- function(formula,nbStates,hierarchical=FALSE)
{
  stateForms<- terms(formula, specials = c(paste0(rep(c("state","toState"),each=nbStates),1:nbStates),"recharge"))
  newformula<-formula
  formulaStates <- vector('list',nbStates*(nbStates-1))
  formulaStates[1:(nbStates*(nbStates-1))] <- list(newformula)
  formterms<-attr(terms.formula(newformula),"term.labels")
  recharge <- NULL
  
  if(nbStates>1){
    if(length(unlist(attr(stateForms,"specials")))){
      #newForm<-attr(stateForms,"term.labels")[-unlist(attr(stateForms,"specials"))]
      newForm <- attr(stateForms,"term.labels")[-which(colSums(attr(stateForms,"factors")[unlist(attr(stateForms,"specials")),,drop=FALSE])>0)]#attr(drop.terms(stateForms,which(attr(stateForms,"factors")[unlist(attr(stateForms,"specials")),]>0)),"term.labels")
      for(i in 1:nbStates){
        if(!is.null(attr(stateForms,"specials")[[paste0("state",i)]])){
          for(j in 1:(nbStates-1)){
            #newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("state",i)]]]))
            stateFact <- attr(stateForms,"factors")[attr(stateForms,"specials")[[paste0("state",i)]],,drop=FALSE]
            newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),colnames(stateFact)[which(stateFact>0)]))
          }
        }
        if(!is.null(attr(stateForms,"specials")[[paste0("toState",i)]])){
          betaInd<-matrix(0,nbStates,nbStates,byrow=TRUE)
          diag(betaInd)<-NA
          betaInd[!is.na(betaInd)] <- seq(1:(nbStates*(nbStates-1)))
          betaInd<-t(betaInd)[,i]
          betaInd<-betaInd[!is.na(betaInd)]
          for(j in betaInd){
            #newForm<-c(newForm,gsub(paste0("toState",i),paste0("betaCol",j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("toState",i)]]]))
            stateFact <- attr(stateForms,"factors")[attr(stateForms,"specials")[[paste0("toState",i)]],,drop=FALSE]
            newForm<-c(newForm,gsub(paste0("toState",i),paste0("betaCol",j),colnames(stateFact)[which(stateFact>0)]))
          }
        }
      }
      if(!is.null(attr(stateForms,"specials")$recharge)){
        reSub <- rownames(attr(stateForms,"factors"))[attr(stateForms,"specials")$recharge]
        recharge <- lapply(reSub,function(x) eval(parse(text=x)))
        #recharge <- stateFormulas(as.formula(paste("~",paste(attr(stateForms,"term.labels")[which(grepl("recharge",attr(stateForms,"term.labels")))]))),1)[[1]]
        hierInd <- FALSE
        for(k in 1:length(reSub)){
          #if(!hierarchical){
          reFact <- attr(stateForms,"factors")[attr(stateForms,"specials")$recharge[k],,drop=FALSE]
          iLevel <- gsub("\") * 1)","",gsub("I((level == \"","",names(which(attr(stateForms,"factors")[,colnames(reFact)[which(reFact>0)]]==2)),fixed=TRUE),fixed=TRUE)
          if(length(iLevel)){
            hierInd <- TRUE
            if(!hierarchical){
              newForm <- c(newForm,gsub(reSub[k],paste0("recharge",iLevel),colnames(reFact)[which(reFact>0)],fixed=TRUE))
            }
            names(recharge)[k] <- paste0("level",iLevel)
          } else {
            if(!hierarchical){
              newForm <- c(newForm,gsub(reSub[k],"recharge",colnames(reFact)[which(reFact>0)],fixed=TRUE))
            }
          }
          #}
        }
        if(length(newForm)){
          newformula<-as.formula(paste("~",ifelse(attr(terms(formula),"intercept"),"","0+"),paste(newForm,collapse="+")))
        } else {
          newformula <- ~1
        }
        if(!hierInd) recharge <- recharge[[1]]
      } else {
        newformula<-as.formula(paste("~",paste(newForm,collapse="+")))
      }
    }
    formulaStates<-stateFormulas(newformula,nbStates*(nbStates-1),spec="betaCol")
    if(length(unlist(attr(terms(newformula, specials = c(paste0("betaCol",1:(nbStates*(nbStates-1))),"cosinor")),"specials")))){
      allTerms<-unlist(lapply(formulaStates,function(x) attr(terms(x),"term.labels")))
      newformula<-as.formula(paste("~",paste(allTerms,collapse="+")))
      formterms<-attr(terms.formula(newformula),"term.labels")
    } else if(is.null(recharge)) {
      formterms<-attr(terms.formula(newformula),"term.labels")
      newformula<-formula
    } else {
      formterms<-attr(terms.formula(newformula),"term.labels")
    }
    if(hierarchical){
      recharge <- expandRechargeFormulas(recharge)
    }
  }  
  return(list(formulaStates=formulaStates,formterms=formterms,newformula=newformula,recharge=recharge))
}

#MUFFWARN
muffWarn <- function(w) {
  
  if(any(grepl("zero-length arrow is of indeterminate angle and so skipped",w)))
    
    invokeRestart("muffleWarning")
  
}

#STATEFORMULAS
stateFormulas <- function(formula,nbStates,spec="state",angleMean=FALSE,data=NULL){
  
  Terms <- terms(formula, specials = c(paste0(spec,1:nbStates),"cosinor","angleFormula","recharge"))
  if(any(grepl("angleStrength\\(",attr(Terms,"term.labels")))) stop("'angleStrength' is defunct in momentuHMM >=1.4.2. Please use 'angleFormula' instead")
  if(any(attr(Terms,"order")>1)){
    if(any(grepl("angleFormula\\(",attr(Terms,"term.labels")[attr(Terms,"order")>1]))) stop("interactions with angleFormula are not allowed")
  }
  
  stateFormula<-list()
  if(length(unlist(attr(Terms,"specials"))) | angleMean){
    varnames <- attr(Terms,"term.labels")
    mainpart <- varnames
    cosInd <- survival::untangle.specials(Terms,"cosinor",order=1:10)$terms
    if(length(cosInd) & angleMean) stop("cosinor models are not supported for angle means")
    angInd <- survival::untangle.specials(Terms,"angleFormula",order=1:10)$terms
    if(length(angInd) & !angleMean) stop("angleFormula models are only allowed for angle means")
    recInd <- survival::untangle.specials(Terms,"recharge",order=1:10)$terms
    if(length(recInd)>1) stop("only a single recharge model is permitted")
    stateInd <- numeric()
    for(j in 1:nbStates){
      tmpInd <- survival::untangle.specials(Terms,paste0(spec,j),order=1)$terms
      if(length(tmpInd)) stateInd<-c(stateInd,tmpInd)
    }
    if(length(cosInd) | length(angInd) | length(recInd) | length(stateInd)){
      mainpart <- varnames[-c(cosInd,angInd,recInd,stateInd)]
    }
    if(angleMean & length(mainpart)){
      tmpmainpart <- mainpart
      mainpart <- character()
      if(any(grepl("cos",tmpmainpart)) | any(grepl("sin",tmpmainpart))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
      for(j in 1:length(tmpmainpart)){
        mainpart <- c(mainpart,paste0(c("sin","cos"),"(",tmpmainpart[j],")"))
      }
    }
    for(j in varnames[cosInd])
      mainpart<-c(mainpart,paste0(gsub("cosinor","cosinorCos",j)),paste0(gsub("cosinor","cosinorSin",j)))
    for(j in varnames[recInd])
      mainpart<-c(mainpart,paste0(gsub(")\\s*$","",gsub("recharge\\(","",j))))
    if(length(angInd)){
      stmp <- prodlim::strip.terms(Terms[attr(Terms,"specials")$angleFormula],specials="angleFormula",arguments=list(angleFormula=list("strength"=NULL,"by"=NULL)))
      if(any(grepl("cos",attr(stmp,"term.labels"))) | any(grepl("sin",attr(stmp,"term.labels")))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
      for(jj in attr(stmp,"term.labels")){
        if(is.null(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength)){
          tmpForm <- ~ - 1
          strengthInd <- FALSE
        }
        else {
          if(!is.na(suppressWarnings(as.numeric((attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))))) stop("angleFormula has invalid strength argument")
          tmpForm <- as.formula(paste0("~-1+",attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))
          if(any(attr(terms(tmpForm),"order")>1)) stop("angleFormula strength argument for ",jj," cannot include term interactions; use the 'by' argument")
          strengthInd <- TRUE
        }
        group <- attr(stmp,"stripped.arguments")$angleFormula[[jj]]$by
        if(!is.null(data)){
          DMterms <- attr(terms(tmpForm),"term.labels")
          factorterms<-names(data)[unlist(lapply(data,is.factor))]
          factorcovs<-paste0(rep(factorterms,times=unlist(lapply(data[factorterms],nlevels))),unlist(lapply(data[factorterms],levels)))
          for(cov in DMterms){
            form<-formula(paste("~",cov))
            varform<-all.vars(form)
            if(any(varform %in% factorcovs)){
              factorvar<-factorcovs %in% varform
              tmpcov<-rep(factorterms,times=unlist(lapply(data[factorterms],nlevels)))[which(factorvar)]
              tmpcovj <- cov
              for(j in 1:length(tmpcov)){
                tmpcovj <- gsub(factorcovs[factorvar][j],tmpcov[j],tmpcovj)
              }
              form <- formula(paste("~ 0 + ",tmpcovj))
            }
            if(strengthInd){
              if(any(model.matrix(form,data)<0)) stop(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength," must be >=0 in order to be used in angleFormula")
              if(any(unlist(lapply(data[all.vars(form)],function(x) inherits(x,"factor"))))) stop("angleFormula strength argument cannot be a factor; use the 'by' argument for factors")
            }
          }
          if(!is.null(group)){
            form<-formula(paste("~",group))
            varform<-all.vars(form)
            if(any(varform %in% factorcovs)){
              factorvar<-factorcovs %in% varform
              varform<-rep(factorterms,times=unlist(lapply(data[factorterms],nlevels)))[which(factorvar)]
            }
            if(any(!(varform %in% names(data)))) stop("angleFormula 'by' argument ",varform[which(!(varform %in% names(data)))]," not found in data")
            if(any(!unlist(lapply(data[varform],function(x) inherits(x,"factor"))))) stop("angleFormula 'by' argument must be of class factor")
          }
        }
        
        if(!is.null(group)) group <- paste0(group,":")
        mainpart<-c(mainpart,paste0(group,ifelse(strengthInd,paste0(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength,":"),""),c("sin","cos"),"(",jj,")"))
      }
    }
    
    for(j in 1:nbStates){
      tmplabs<-attr(Terms,"term.labels")[attr(Terms,"specials")[[paste0(spec,j)]]]
      if(length(tmplabs)){
        tmp<- terms(as.formula(paste("~",substr(tmplabs,nchar(paste0(spec,j))+1,nchar(tmplabs)),collapse="+")),specials=c("cosinor","angleFormula","recharge"))
        
        tmpnames<-attr(tmp,"term.labels")
        if(any(grepl("angleStrength\\(",tmpnames))) stop("'angleStrength' is defunct in momentuHMM >=1.4.2. Please use 'angleFormula' instead")
        if(any(attr(tmp,"order")>1)){
          if(any(grepl("angleFormula\\(",tmpnames[attr(tmp,"order")>1]))) stop("interactions with angleFormula are not allowed")
        }
        mp<-tmpnames
        if(!is.null(unlist(attr(tmp,"specials"))) | angleMean){
          cosInd <- survival::untangle.specials(tmp,"cosinor",order=1:10)$terms
          if(length(cosInd) & angleMean) stop("cosinor models are not supported for angle means")
          angInd <- survival::untangle.specials(tmp,"angleFormula",order=1:10)$terms
          if(length(angInd) & !angleMean) stop("angleFormula models are only allowed for angle means")
          recInd <- survival::untangle.specials(tmp,"recharge",order=1:10)$terms
          if(length(recInd)) stop("recharge models cannot be state-dependent")
          if(length(cosInd) | length(angInd)){
            mp <- c(tmpnames[-c(cosInd,angInd)])
          }
          if(angleMean & length(mp)){
            tmpmp <- mp
            mp <- character()
            if(any(grepl("cos",tmpmp)) | any(grepl("sin",tmpmp))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
            for(jj in 1:length(tmpmp)){
              mp <- c(mp,paste0(c("sin","cos"),"(",tmpmp[jj],")"))
            }
          }
          for(i in tmpnames[cosInd])
            mp<-c(mp,paste0(gsub("cosinor","cosinorCos",i)),paste0(gsub("cosinor","cosinorSin",i)))
          for(i in tmpnames[recInd])
            mp<-c(mp,paste0(gsub(")\\s*$","",gsub("recharge\\(","",i))))
          if(length(angInd)){
            stmp <- prodlim::strip.terms(tmp[attr(tmp,"specials")$angleFormula],specials="angleFormula",arguments=list(angleFormula=list("strength"=NULL,"by"=NULL)))
            if(any(grepl("cos",attr(stmp,"term.labels"))) | any(grepl("sin",attr(stmp,"term.labels")))) stop("sorry, the strings 'cos' and 'sin' are reserved and cannot appear in mean angle formulas and/or covariate names")
            for(jj in attr(stmp,"term.labels")){
              if(is.null(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength)){
                tmpForm <- ~ - 1
                strengthInd <- FALSE
              }
              else {
                if(!is.na(suppressWarnings(as.numeric((attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))))) stop("angleFormula has invalid strength argument")
                tmpForm <- as.formula(paste0("~-1+",attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength))
                if(any(attr(terms(tmpForm),"order")>1)) stop("angleFormula strength argument cannot include term interactions; use the 'by' argument")
                strengthInd <- TRUE
              }
              group <- attr(stmp,"stripped.arguments")$angleFormula[[jj]]$by
              if(!is.null(data)){
                if(strengthInd){
                  if(any(model.matrix(tmpForm,data)<0)) stop(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength," must be >=0 in order to be used in angleFormula")
                  if(inherits(data[[attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength]],"factor")) stop("angleFormula strength argument cannot be a factor; use the 'by' argument for factors")
                }
                if(!is.null(group)){
                  varform <- all.vars(as.formula(paste0("~",group)))
                  if(any(!(varform %in% names(data)))) stop("angleFormula 'by' argument ",varform[which(!(varform %in% names(data)))]," not found in data")
                  if(any(!unlist(lapply(data[varform],function(x) inherits(x,"factor"))))) stop("angleFormula 'by' argument must be of class factor")
                }
              }
              if(!is.null(group)) group <- paste0(group,":")
              mp<-c(mp,paste0(group,ifelse(strengthInd,paste0(attr(stmp,"stripped.arguments")$angleFormula[[jj]]$strength,":"),""),c("sin","cos"),"(",jj,")"))
            }
          }
        }
      } else {
        tmp <- Terms
        mp <- character()
      }
      stateFormula[[j]]<-as.formula(paste("~",paste(c(attr(tmp,"intercept"),mainpart,mp),collapse = " + "),collapse=" + "))
    }
  } else {
    for(j in 1:nbStates){
      stateFormula[[j]] <- formula
    }
  }
  stateFormula
}

#GET_STAT
get_stat <- function(beta,covs,nbStates,i,betaRef,betaCons,workBounds=matrix(c(-Inf,Inf),length(betaCons),2,byrow=TRUE),mixture=1,ref=1:nbStates) {
  
  gamma <- get_gamma(beta,covs,nbStates,1:nbStates,1:nbStates,betaRef,betaCons,workBounds,mixture)
  
  solve(t(diag(length(ref))-gamma[ref,ref]+1),rep(1,length(ref)))[i]
  
}

#DELTS_BC
delta_bc <- function(m){
  
  if(is.momentuHMM(m) | is.miSum(m)){
    if(!is.null(m$conditions$fit)){
      if(!m$conditions$fit) warning("The given model hasn't been fitted.")
    } else m$conditions$fit <- TRUE
    if(is.null(m$conditions$workBounds)){
      distnames <- names(m$conditions$dist)
      
      parCount<- lapply(m$conditions$fullDM,ncol)
      for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
        parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
      }
      parindex <- c(0,cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
      names(parindex) <- distnames
      
      workBounds <- vector('list',length(distnames))
      names(workBounds) <- distnames
      if(is.miSum(m)){
        beta <- m$Par$beta$beta$est
        delta <- m$Par$beta$delta$est
      } else {
        beta <- m$CIbeta$beta$est
        delta <- m$CIbeta$delta$est
      }
      beta <- list(beta=beta,g0=m$mle$g0,theta=m$mle$theta)
      m$conditions$workBounds <- getWorkBounds(workBounds,distnames,m$mod$estimate,parindex,parCount,m$conditions$DM,beta,delta)
    }
    if(length(m$stateNames)>1 && is.null(m$conditions$betaCons)){
      if(is.miSum(m) & !is.null(m$Par$beta$beta)) m$conditions$betaCons <- matrix(1:length(m$Par$beta$beta$est),nrow(m$Par$beta$beta$est),ncol(m$Par$beta$beta$est))
      else if(is.momentuHMM(m) & !is.null(m$mle$beta)) m$conditions$betaCons <- matrix(1:length(m$mle$beta),nrow(m$mle$beta),ncol(m$mle$beta))
    }
    if(is.null(m$conditions$betaRef)) m$conditions$betaRef <- as.integer(1:length(m$stateNames))
    if(is.momentuHMM(m)){
      if(is.null(m$mod$wpar)) m$mod$wpar <- m$mod$estimate
      if(is.null(m$mod$Sigma) & !is.null(m$mod$hessian)) m$mod$Sigma <- MASS::ginv(m$mod$hessian)
    } else {
      ####### compatability hack for change to MIcombine in momentuHMM >= 1.4.3 ######
      if(is.null(m$conditions$optInd)){
        for(i in names(m$conditions$dist)){
          m$conditions$cons[[i]]<-rep(1,length(m$conditions$cons[[i]]))
          m$conditions$workcons[[i]]<-rep(0,length(m$conditions$workcons[[i]]))
          m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
        }
      }
      ################################################################################
    }
    if(is.null(m$conditions$mixtures)) m$conditions$mixtures <- 1
    if(is.null(m$covsPi)) m$covsPi <- matrix(1,length(unique(m$data$ID)),1)
    if(is.null(attr(m$data,"coords")) & !is.null(m$data$x) & !is.null(m$data$y)) attr(m$data,"coords") <- c("x","y")
  } else if(!is.miHMM(m) & any(unlist(lapply(m,is.momentuHMM)))){
    m <- HMMfits(m)
  }
  m
}


#IS.MOMENTUM
is.momentuHMM <- function(x)
  inherits(x,"momentuHMM")


#IS.MISUM
is.miSum <- function(x)
  inherits(x,"miSum")


#GETCOVS
getCovs <-function(m,covs,ID,checkHier=TRUE){
  
  if(is.null(covs)){
    
    if(inherits(m,"hierarchical")) covs <- as.data.frame(lapply(m$data,function(x) x[which.max(!is.na(x))]))
    
    else covs <- m$data[which(m$data$ID %in% ID),][1,]
    
    for(j in names(m$data)[which(unlist(lapply(m$data,function(x) any(class(x) %in% meansList))))]){
      
      if(inherits(m$data[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(m$data[[j]][which(m$data$ID %in% ID)][!is.na(m$data[[j]][which(m$data$ID %in% ID)])])
      
      else covs[[j]]<-mean(m$data[[j]][which(m$data$ID %in% ID)],na.rm=TRUE)
      
    }
    
  } else {
    
    if(!is.data.frame(covs)) stop('covs must be a data frame')
    
    if(nrow(covs)>1) stop('covs must consist of a single row')
    
    if(is.null(recharge))
      
      if(!all(names(covs) %in% names(m$data))) stop('invalid covs specified')
    
    else 
      
      if(!all(names(covs) %in% c(names(m$data),"recharge"))) stop('invalid covs specified')
    
    if(any(names(covs) %in% "ID")) covs$ID<-factor(covs$ID,levels=unique(m$data$ID))
    
    if(checkHier && inherits(m,"hierarchical") && any(names(covs) %in% "level")) stop("covs$level cannot be specified for hierarchical models")
    
    for(j in names(m$data)[which(names(m$data) %in% names(covs))]){
      
      if(inherits(m$data[[j]],"factor")) covs[[j]] <- factor(covs[[j]],levels=levels(m$data[[j]]))
      
      if(is.na(covs[[j]])) stop("check covs value for ",j)
      
    }  
    
    for(j in names(m$data)[which(!(names(m$data) %in% names(covs)))]){
      
      if(any(class(m$data[[j]]) %in% meansList)){
        
        if(inherits(m$data[[j]],"angle")) covs[[j]] <- CircStats::circ.mean(m$data[[j]][!is.na(m$data[[j]])])
        
        else covs[[j]]<-mean(m$data[[j]],na.rm=TRUE)
        
      } else {
        
        if(inherits(m,"hierarchical")) covInd <- which.max(!is.na(m$data[[j]]))
        
        else covInd <- 1
        
        covs[[j]] <- m$data[[j]][covInd]
        
      }
      
    }
    
  }
  
  covs
  
}


#GET_NCMEAN
get_ncmean <- function(distnames,fullDM,circularAngleMean,nbStates){
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    if(!isFALSE(circularAngleMean[[i]])) {
      meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  list(nc=nc,meanind=meanind)
}
