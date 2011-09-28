
#############################################
### Convert list of platforms and responses
### to DR-suitable data
#############################################

convertToDr <- function(dataIntersection, response, nPlatforms = length(data)) {
	int <- rownames(dataIntersection[[1]])
	samples <- colnames(dataIntersection[[1]])
	dataForDr <- list()
	responseForCommonPatients <- response[(response[,1]%in%samples),]
	rownames(responseForCommonPatients) <- responseForCommonPatients[,1]
	###r = responseForCommonPatients[samples, ]
	for (i in c(1:length(int))) {
		if (i%%1000 == 0) { print(i)}
		dataForDr[[i]] <- data.frame(x = t(dataIntersection[[1]][int[i],samples]) )
		if (nPlatforms >=2) {
			for (j in c(2:length(dataIntersection))) {
				dataForDr[[i]] <- cbind(dataForDr[[i]],
							  t(dataIntersection[[j]][int[i],samples]) )
			}
		}
		dataForDr[[i]] <- cbind(dataForDr[[i]],
					  responseForCommonPatients[samples, 2])
		colnames(dataForDr[[i]]) <- c(names(dataIntersection), "response")
	}
	names(dataForDr) <- int
	return(dataForDr)
}						


##########################################
### Compute dimension reduction stats
##########################################

computeScore <- function(data, columns, method) {
	dat <- data[, c(columns, ncol(data))]
	stats <- 0
	if (method == "dev") {
		dFit <- glm(response~. , dat, family=binomial('logit'))
		stats <- dFit$null.deviance - dFit$deviance
	}
	if (method == "aic") {
		fullFit <- glm(response~. , dat, family=binomial('logit'))
		nullFit <- glm(response~1 , dat, family=binomial('logit'))
		stats <- nullFit$aic - fullFit$aic
	}
	if (method == "bic") {
		fullFit <- glm(response~. , dat, family=binomial('logit'))
		nullFit <- glm(response~1 , dat, family=binomial('logit'))
		stats <- AIC(nullFit, k = log(fullFit$df.null+1))-
			AIC(fullFit, k = log(fullFit$df.null+1))
	}
	return(stats)
}

computeDrStat <- function(data, columns= c(1:(ncol(data)-1)), method = "dev", integrate = TRUE) {
	if (integrate) {
		stats <- lapply(data, computeScore, columns, method)
		listOut <- list()
		listOut$integrated <- unlist(stats)
	} else {
		listOut <- list()
		for (p in columns) {
			stats <- lapply(data, computeScore, columns = p, method)
			listOut[[p]] <- unlist(stats)
			names(listOut)[p] <- colnames(data[[1]])[p]
		}
	}
	return(listOut)
}

#######################################################
###set a function that run geneSetTest after evaluating if genes are present in data
###also turns to absolute values if needed
#######################################################

runGSE <- function(fgs, data, absolute=TRUE, gseFunc=NULL, ...) {
	
###absolute value
	#print(summary(data))
	if (absolute) {
		data <- abs(data)
		#print(summary(data))
	}

###compute Pvalues after checking if genes are available
	gsePval <- NULL
	if (sum(names(data)%in%fgs)==length(data) | sum(names(data)%in%fgs)==0 ) {
		gsePval <- NA
          } else {
		  if(is.null(gseFunc)) {
###with geneSetTest
			  selection <- names(data)%in%fgs
			  gsePval <- geneSetTest(selection, data, ...)
		  } else {
###with user defined function
			  selection <- names(data)%in%fgs
			  gsePval <- gseFunc(selection, data, ...)
		  }
	  }

###output
	return(gsePval)
}


#######################################################
###set a function that performs runGSE  over a list of FGS
#######################################################

runGSEforFGSlist <- function(fgsList, data, ...) {
	out <- sapply(fgsList,runGSE, data, ...)
	return(out)
}


#######################################################
###set a function that makes the test over a list of lists of FGS
#######################################################

runGSEforFGSlistBatch <- function(listFgsList, data, ...) {
	out <- lapply(listFgsList, runGSEforFGSlist, data, ...)
	return(out)
}


#######################################################
###set a function that run the test over a list of lists of FGS over a list of data sets
#######################################################

runBatchGSE <- function(dataList, fgsList, ...) {
	func <- function(x, y, ...){out <- runGSEforFGSlistBatch(y, x, ...)}
	out <- lapply(dataList, func, y = fgsList, ...)
	return(out)
}


#######################################################
###set a function that integrates several pvalues from GSE analysis
#######################################################

combineGSE <- function(gseOut, method){
	nms <- names(gseOut)
	fgs <- unique(unlist(lapply(gseOut,names)))
	tmp <- unlist(gseOut,recursive=F)
	tmp <- lapply(fgs,function(x,y) data.frame(y[grep(x,names(y))]), y=tmp)
	names(tmp) <- fgs
	tmp <- lapply(tmp, myMean, method)
	out <- list(combinedScore=tmp)
	return(out)
}


#######################################################
###computes mean by row in a data.frame using various  methods, NOT EXPORTED
#######################################################

myMean <- function(matM, method, ...){
###geometric mean
	geoMeanM <- function(matM, ...) {
		out <- apply(matM,1,function(x){exp(mean(log(x), ...))})
	}
###random
	randomM <- function(matM, ...) {
		out <- apply(matM,1,sample,size=1, ...)
	}
###mean
	meanM <- function(matM, ...) {
		out <- apply(matM,1,mean, ...)
	}
###median
	medianM <- function(matM, ...) {
		out <- apply(matM,1,median)
	}
###max
	maxM <- function(matM, ...) {
		out <- apply(matM,1,max)
	}
###min
	minM <- function(matM, ...) {
		out <- apply(matM,1,min)
	}
###select method to remove redundancy
	choices <- c("geometricMean","mean","median","random","min","max")
	method <- match.arg(method, choices)
	if(method%in%choices){
		out <- switch(method,
			      geometricMean = geoMeanM(matM, ...),
			      mean = meanM(matM, ...),
			      random = randomM(matM, ...),
			      median = medianM(matM, ...),
			      min = minM(matM, ...),
			      max = maxM(matM, ...)
			      )
	} else {stop(paste("Select an available method (", paste(method, collapse=", "), ")", sep=""))}
	return(out)
}


#######################################################
###set a function that integrates several pvalues from GSE analysis
#######################################################

adjustPvalGSE <- function(gseOut, proc = "BH" , alpha = 0.05, na.rm = FALSE){
	out <- lapply(gseOut, function(x, ...) lapply(x,pAdjust, proc, alpha, na.rm) )
	return(out)
}


#######################################################
###adjust on vector of p-values; NOT EXPORTED
#######################################################

pAdjust <- function(rawp, ...) {
	Adj <- mt.rawp2adjp(rawp, ...)
	Adj <- Adj$adjp[order(Adj$index),]
	rownames(Adj) <- names(rawp)
	return(Adj)
}



