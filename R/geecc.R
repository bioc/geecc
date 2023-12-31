

.loglinmodel <- function( CT, null.model=~ factor1 + factor2 + factor3 ){

	res <- MASS::loglm(formula = as.formula(null.model) , data = CT)
	obj <- list(statistic=setNames(c(res$pearson), c("Pearson"))
		, parameter=setNames(res$df, "df")
		, estimate=setNames(getOddsRatio(CT), "odds ratio")
		, p.value=(c(summary(res)$tests[2,3]))
	)
return(obj)
}



.my_separator <- function(){return("_##_")}




.getContingencyCube <- function( CT, x1__Ix_1_, x2__Ix_2_, x1__Ix_2_, x2__Ix_1_, x__1, x__2 ){
	CT[1,1,1] <- length( intersect( x1__Ix_1_, x__1 ) )
	CT[1,1,2] <- length( intersect( x1__Ix_1_, x__2 ) )
	CT[1,2,1] <- length( intersect( x1__Ix_2_, x__1 ) )
	CT[2,1,1] <- length( intersect( x2__Ix_1_, x__1 ) )

	CT[1,2,2] <- length( intersect( x1__Ix_2_, x__2 ) )
	CT[2,1,2] <- length( intersect( x2__Ix_1_, x__2 ) )
	CT[2,2,1] <- length( intersect( x2__Ix_2_, x__1 ) )
	CT[2,2,2] <- length( intersect( x2__Ix_2_, x__2 ) )
return(CT)
}
.getContingencyTable <- function( CT, x1__Ix_1_, x2__Ix_2_, x1__Ix_2_, x2__Ix_1_ ){
	CT[1:4] <- c( length( ( x1__Ix_1_ ) ), length( ( x2__Ix_1_ ) ), length( ( x1__Ix_2_ ) ), length( ( x2__Ix_2_ ) ) )

return(CT)
}



.getItemsInEachCategory <- function(obj){ return(obj@categories) }
.getNamesOfEachCategory <- function(obj){ return(lapply(obj@categories, names)) }
.getNumberOfFactorLevels <- function(obj){ return(sapply(obj@categories, length)) }
.getOptOfCategory <- function(obj){ return(obj@options) }


#test when choice between hypergeom and chisq
.performTest_approx <- function(frm, CT, CT_t, minExpectedValues, approx, nthreads, test.direction){
	if( minExpectedValues < approx ){
		return(hypergeom.test(CT_t, nthreads=nthreads, alternative=test.direction))
	}
	return(.loglinmodel(CT, null.model=frm))
}

#test when no approx possible
.performTest <- function(frm, CT){
	res <- list()
	res <- .loglinmodel(CT, null.model=frm)
	#res$estimate <- getOddsRatio(CT, za=TRUE)
return(res)
}


runConCub <- function( obj, filter, nthreads=2, subset=NULL, verbose=list(output.step=0, show.cat1=FALSE, show.cat2=FALSE, show.cat3=FALSE)){

	NCATS <- length(obj@categories)
	nms_categories <- names(obj@categories)
	
	# check if concub-filter object passed fits to concub object
	if( (length(nms_categories) != length(filter@names)) || (intersect(nms_categories, filter@names) != nms_categories)  ){
		warning("Names from concubfilter-object (", paste0(nms_categories, collapse=",") ,") and concub-object (", paste0(filter@names, collapse=","), ") do not match.")
		return(obj)
	}

	N <- length(obj@population) #size of population
	N_factor <- .getNumberOfFactorLevels(obj)

	items_factor <- .getItemsInEachCategory(obj) # items of each categorie (list of lists)
	opt_factor <- .getOptOfCategory(obj)
	rng_factor <- vector('list', NCATS)


	sub_categories <- .getNamesOfEachCategory(obj) # names of categories
	if( !is.null(subset) ){ for(nm in names(subset)){sub_categories[[ nm ]] <- intersect(subset[[nm]], sub_categories[[ nm ]]); } } # use intersect to avoid incorp of variables that have been dropped because, e.g., they are empty
	len_sub_categories <- sapply(sub_categories, length)


	ttt <- .getTypeOf_transformTable(form=obj@null.model, nms=nms_categories)
	do.strat <- FALSE
	if(NCATS==3){do.strat <- obj@options[[ 3 ]][['strat']]}
	frm <- obj@null.model
	if( do.strat == TRUE & ttt[[1]] == 'mi' ){
		frm <- update(frm, as.formula(paste("~.-", names(obj@categories)[3])))
	}
	message(paste("Testing: counts ~ ", as.character(obj@null.model)[2], " (", ttt[[1]],")", sep=""), sep="")
	if( any(sapply(opt_factor, function(x){return(x$grouping != "none")}) ) ){
		message(paste("Grouping: ", paste(names(opt_factor), ":", sapply(opt_factor, function(x){return(x$grouping)}), sep="", collapse=", "), sep=""), sep="")
	}

	####
	## grouping of ordinal variables
	## 
    RNG <- setNames( vector("list", NCATS), nms_categories ) # stores the variable names for each category that should be taken together
	if( sum( sapply( opt_factor, function(x){ x[['grouping']] != 'none' } ) ) > 1 ){ warning("Grouping-option for multiple categories. Results might be hard to interprete.") }
	for( g in nms_categories ){
		opt_grouping <- opt_factor[[g]][[ 'grouping' ]]
		opt_width <- opt_factor[[g]][[ 'width' ]]
		local_len <- len_sub_categories[g]
		local_rng <- sub_categories[[g]]
		RNG[[ g ]] <- setNames( vector("list", local_len), local_rng )

		for( g2 in 1:local_len ){
			term <- local_rng[g2]
			term2cum <- term;

 			if( opt_grouping == 'none' ){ term2cum <- term; }
			if( opt_grouping == 'cumf' ){
				if( g2 == local_len ){ next; } #don't include last
				term2cum <- local_rng[ 1:g2 ]
			}
			if( opt_grouping == 'cumr' ){
				if( g2 == 1 ){ next; } #don't include first
				term2cum <- local_rng[ g2:local_len ]
			}
			if( opt_grouping == 'sw' ){
				l <- opt_width; start <- max(c(1, g2-l)); stop <- min(c(g2+l, local_len))
				term2cum <- local_rng[start : stop]
			}

 			RNG[[ g ]][[ term ]] <- term2cum
 		}

		if( opt_grouping == 'cumf' ){ sub_categories[[g]] <- sub_categories[[g]][1:(local_len-1)] }
		if( opt_grouping == 'cumr' ){ sub_categories[[g]] <- sub_categories[[g]][2:local_len] }
	}
	##
	####

	## prepare progress output for terminal
	Len_term3 <- c();
	if(NCATS == 3 && verbose$show.cat3){
		L <- length(sub_categories[[3]])
		Len_term3 <- setNames(vector("list", L), sub_categories[[3]]);
		for( l in 1:L ){
			x <- ifelse( l==1, sub_categories[[3]][L], sub_categories[[3]][l-1] )
			Len_term3[[l]] <- paste(paste(rep("\b", times=nchar(x)), collapse=""), "\t", sub_categories[[3]][l], collapse="", sep="")
		}
	}
	##

	my_separator <- .my_separator()
	tmp0_nms1 <- c(t(outer( sub_categories[[1]], sub_categories[[2]], paste, sep=my_separator )))
	if(NCATS==3){tmp0_nms1 <- c( t(outer(tmp0_nms1, sub_categories[[3]], paste, sep=my_separator) ))}
	
	
	
	## pre-calc some sets
	tmp1 <- items_factor[[ nms_categories[2] ]]; 
	tmp2 <- RNG[[ nms_categories[2] ]]
 	PreCalc__x_1_ <- setNames( vector("list", length(sub_categories[[ nms_categories[2] ]])), names( sub_categories[[ nms_categories[2] ]] ) )
	for( t in as.integer(1:length(sub_categories[[ nms_categories[2] ]])) ){
        PreCalc__x_1_[[t]] <- unique(.special1( tmp1, tmp2, t )) # 8-fold faster and 60-fold less memory, compared to version 1.1.10
	}
	PreCalc__x__1 <- PreCalc__x__2 <- list()
	if( NCATS == 3 ){ 
		loc_nm3 <- nms_categories[3]
        tmp1 <- items_factor[[ loc_nm3 ]]; 
        tmp2 <- RNG[[ loc_nm3 ]]
        PreCalc__x__1 <- setNames( vector("list", length(sub_categories[[ loc_nm3 ]])), names( sub_categories[[ loc_nm3 ]] ) )

        for( t in as.integer(1:length(sub_categories[[ loc_nm3 ]])) ){
            PreCalc__x__1[[t]] <- unique(.special1( tmp1, tmp2, t )) # how fast?
        }
	}
	##
	
	ITER <- setNames(vector('list', length(tmp0_nms1)), tmp0_nms1)
	loc_nm1 <- nms_categories[1]
	for( g1 in 1:length(sub_categories[[ loc_nm1 ]]) ){
		term1 <- sub_categories[[ loc_nm1 ]][ g1 ]; if(verbose$show.cat1){cat(term1, sep="")}; notterm1 <- paste('not_', term1, sep="")

		x1__ <- unique( unlist(items_factor[[ loc_nm1 ]][ RNG[[ loc_nm1 ]][[ term1 ]] ], use.names=FALSE) )
		x2__ <- setdiffPresort( obj@population, x1__ ) # already sorted

		RES_CAT2 <- list()
		loc_nm2 <- nms_categories[2]
		for( g2 in 1:length(sub_categories[[ loc_nm2 ]]) ){
			term2 <- sub_categories[[ loc_nm2 ]][ g2 ]; if(verbose$show.cat2){cat("\r\t", term2, "\t", sep="")}

			x_1_ <- PreCalc__x_1_[[g2]]
			x_2_ <- setdiffPresort(obj@population, PreCalc__x_1_[[g2]]) # already sorted; could be pre-calculted, but calculated again and again to avoid memory overflow

			x1__Ix_1_ <- intersect( x1__, x_1_ );
			x1__Ix_2_ <- intersectPresort( x_2_, x1__ );
			x2__Ix_1_ <- intersectPresort( x2__, x_1_ );
			x2__Ix_2_ <- setdiffPresort( obj@population, c( x1__Ix_1_, x1__Ix_2_, x2__Ix_1_ ) ) # already sorted
            
			RES_CAT3 <- list()
			if(NCATS==2){
				res <- list()
                subpop <- as.character(( x1__Ix_1_)); len_subpop <- length(subpop)
 				CT <- array( c(len_subpop, vapply( list(x2__Ix_1_, x1__Ix_2_, x2__Ix_2_), length, FUN.VALUE=123 )), dim=c(2,2), dimnames=list( factor1=c( term1, notterm1 ), factor2=c( term2, paste('not_', term2, sep="") ) ) )
				names(dimnames(CT)) <- nms_categories

				
				ExpectedValues <- hypergea::getExpectedValues(CT)

				if( !( ( skip.zeroobs(filter) && len_subpop == 0 )
					|| ( skip.min.obs(filter) >= len_subpop )
					|| all( skip.min.group(filter) - c( length( x1__ ), length(x_1_ ) ) > 0 ) ) ){ #perform test when none of the conditions is fulfilled
					res <- .performTest_approx(frm, CT, CT, minExpectedValues=min(ExpectedValues), approx=obj@approx, nthreads, test.direction(filter))
				}

				res <- c(res, list(subpop=subpop))
				RES_CAT2[[ term2 ]] <- res
			}
			if(NCATS==3){ # ToDo: optimze for speed  etc
				RES_CAT3 <- list()

				loc_nm3 <- nms_categories[3]
				for( g3 in 1:length(sub_categories[[ loc_nm3 ]]) ){
					res <- list()
					term3 <- sub_categories[[ loc_nm3 ]][ g3 ]; if(verbose$show.cat3){message(Len_term3[[ term3 ]])}

					x__1 <- PreCalc__x__1[[g3]]
					x__2 <- setdiffPresort(obj@population, PreCalc__x__1[[g3]]) 
					
					subpop <- as.character(intersect( x1__Ix_1_, x__1 )); len_subpop <- length(subpop)
					
					CT <- array(NA, dim=c(2,2,2), dimnames=list( factor1=c( term1, notterm1 ), factor2=c( term2, paste('not_', term2, sep="") ), factor3=c( term3, paste('not_', term3, sep="") ) ))
					names(dimnames(CT)) <- nms_categories
					res <- list("estimate"=1, p.value=1, subpop=subpop)


					if( skip.zeroobs(filter) && len_subpop == 0 || skip.min.obs(filter) >= len_subpop ){ RES_CAT3[[ term3 ]] <- res; next; } # skip if not enough items in x000
					group.len <- sapply( list( x1__, x_1_, x__1 ), length )
					if( any( skip.min.group(filter) >= group.len ) ){ RES_CAT3[[ term3 ]] <- res; next; }
					
					
					CT <- .getContingencyCube( CT, x1__Ix_1_, x2__Ix_2_, x1__Ix_2_, x2__Ix_1_, x__1, x__2 )
					ExpectedValues <- getExpectedValues(CT)
					
					if( ttt[[1]] %in% c("mi") ){
						res <- .performTest_approx(frm, CT, CT, minExpectedValues=min(ExpectedValues), approx=obj@approx, nthreads=nthreads, test.direction(filter))
					}else{
						if( ttt[[1]] %in% c("sp.1", "sp.2", "sp.3") ){
							CT_t <- .transformTable(CT, x=ttt)
							ExpectedValues_t <- getExpectedValues(CT_t)
							res <- .performTest_approx(frm, CT, CT_t, minExpectedValues=min(ExpectedValues_t), approx=obj@approx, nthreads=nthreads, test.direction(filter))
						}else{
							res <- .performTest(frm, CT)
						}
					}

					RES_CAT3[[ term3 ]] <- c(res, list(subpop=subpop))
				}  # END loop g3

				tmp_nms3 <- names( RES_CAT3 )
				RES_CAT2[ paste( term2, tmp_nms3, sep=my_separator ) ] <- RES_CAT3[ tmp_nms3 ]
			}
			
			if( verbose[['output.step']] > 0 && ( (g2 %% verbose[['output.step']] == 0) || (g2 == length(rng_factor[[2]])) ) ){
				if(verbose$show.cat2){cat("\n\tpassed: second category (", nms_categories[2], "), variable ", term2, " (", g2, ")", sep="")}
				if(verbose$show.cat2){cat("\n")}
			}

		}  # END loop g2
		tmp_nms2 <- names( RES_CAT2 )
		ITER[ paste( term1, tmp_nms2, sep=my_separator ) ] <- RES_CAT2[ tmp_nms2 ]

		if( verbose$show.cat1 ){cat("\n")}
		if(g1%%10==0){gc()}
 	} # END loop g1
 	cat("\n")
 	
 	obj@test.result <- ITER

return(obj)
}



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

.filter_drop.wrongdir.layer <- function(filter, NCATS, ODDSRATIO, PVAL, nms_categories){
	if( !(filter@test.direction %in% c("both", "two.sided")) ){

		local__filter_drop.wrongdir.slice <- filter@drop.wrongdir.layer
		local__filter_test.direction <- filter@test.direction

		cond <- paste( "all( c(x) < 1, na.rm=TRUE )" )
		if(local__filter_test.direction %in% c("less", "under")){ cond <- paste( "all( c(x) > 1, na.rm=TRUE )" ) }

		for(d in 1:NCATS){
			dn <- nms_categories[d]
			if( local__filter_drop.wrongdir.slice[dn] == TRUE ){
				keep <- apply(ODDSRATIO, d, function(x){if( eval(parse(text=cond)) ){return(FALSE)}else{return(TRUE)} })
				ODDSRATIO <- .keepLayer(ODDSRATIO, keep, NCATS, d)
				PVAL <- .keepLayer(PVAL, keep, NCATS, d)
				if((.haszerodim(ODDSRATIO))){ return(NULL) }
			}
		}
	}
return(list(ODDSRATIO, PVAL))
}

.filter_drop.insignif.layer <- function( filter, NCATS, ODDSRATIO, PVAL, nms_categories ){
	local__filter_drop.insignif.slice <- drop.insignif.layer(filter)
	local__filter_p.value <- p.value(filter)
	for(d in 1:NCATS){
		dn <- nms_categories[d]
		if( local__filter_drop.insignif.slice[dn] == TRUE ){
			keep <- apply(PVAL, d, function(x){if( all( x > local__filter_p.value, na.rm=TRUE ) ){return(FALSE)}else{return(TRUE)} })
			ODDSRATIO <- .keepLayer(ODDSRATIO, keep, NCATS, d)
			PVAL <- .keepLayer(PVAL, keep, NCATS, d)
			if((.haszerodim(ODDSRATIO))){ return(NULL) }
		}
	}
return(list(ODDSRATIO, PVAL))
}
.filter_drop.lowl2or.layer <- function( filter, NCATS, ODDSRATIO, PVAL, nms_categories ){
	local__1 <- drop.lowl2or.layer(filter)
	local__val <- minimum.l2or(filter)
	for(d in 1:NCATS){
		dn <- nms_categories[d]
		if( local__1[dn] == TRUE ){
			keep <- apply(abs(log2(ODDSRATIO)), d, function(x){if( all( x < local__val, na.rm=TRUE ) ){return(FALSE)}else{return(TRUE)} })
			ODDSRATIO <- .keepLayer(ODDSRATIO, keep, NCATS, d)
			PVAL <- .keepLayer(PVAL, keep, NCATS, d)
			if(.haszerodim(ODDSRATIO)){ return(NULL) }
		}
	}
return(list(ODDSRATIO, PVAL))
}


.haszerodim <- function(x){ if(any(dim(x) == 0)){ warning(paste0("Obtained dimension of zero: ", paste0(dim(x), collapse=","), ". Weaken your filters.")); return(TRUE) }else{return(FALSE)} }
.keepLayer <- function(x, keep, NCATS, d){
	if(NCATS==2){
		if(d==1){x <- x[keep,, drop=FALSE]; }
		if(d==2){x <- x[,keep, drop=FALSE]; }
	}

	if(NCATS==3){
		if(d==1){x <- x[keep,,, drop=FALSE]; }
		if(d==2){x <- x[,keep,, drop=FALSE]; }
		if(d==3){x <- x[,,keep, drop=FALSE]; }
	}
return(x)
}

.filter_dropLayer <- function( filter, NCATS, ODDSRATIO, PVAL, nms_categories ){
	for(d in 1:NCATS){
		dn <- nms_categories[d]
		KEEP <- c()
		if( drop.lowl2or.layer(filter)[dn] ){
			keep <- apply(abs(log2(ODDSRATIO)), d, function(x){if( all( x < minimum.l2or(filter), na.rm=TRUE ) ){return(FALSE)}else{return(TRUE)} })
			KEEP <- cbind(KEEP, keep)
		}
	
		if( drop.insignif.layer(filter)[dn]  ){
			keep <- apply(PVAL, d, function(x){if( all( x > p.value(filter), na.rm=TRUE ) ){return(FALSE)}else{return(TRUE)} })
			KEEP <- cbind(KEEP, keep)
		}
		
		if( !(test.direction(filter) %in% c("both", "two.sided")) ){
			if( drop.wrongdir.layer(filter)[dn]  ){
				cond <- paste( "all( c(x) < 1, na.rm=TRUE )" ) # test.direction(filter) %in% c("greater", "over")
				if(test.direction(filter) %in% c("less", "under")){ cond <- paste( "all( c(x) > 1, na.rm=TRUE )" ) }
				keep <- apply(ODDSRATIO, d, function(x){if( eval(parse(text=cond)) ){return(FALSE)}else{return(TRUE)} })
				KEEP <- cbind(KEEP, keep)
			}
		}

		if( is.matrix(KEEP) ){
			keep <- apply(KEEP, 1, all)
			ODDSRATIO <- .keepLayer(ODDSRATIO, keep, NCATS, d)
			if(.haszerodim(ODDSRATIO)){ return(NULL) }
			PVAL <- .keepLayer(PVAL, keep, NCATS, d)
		}
	}
return( list(ODDSRATIO, PVAL) )
}

filterConCub <- function(obj, filter, p.adjust.method='none', verbose=1){
	NCATS <- length(obj@categories)
	nms_categories <- names(obj@categories)

	ResultTest <- obj@test.result

	##############################################################
	my_separator <- .my_separator()

	all_lab <- names(ResultTest)

	split_all_lab <- strsplit(all_lab, split=my_separator)
	dmnms <- matrix(unlist(split_all_lab), nrow=length(split_all_lab), byrow=TRUE)
	dmnms0 <- setNames(apply(dmnms, 2, unique), nms_categories)

	## guarantee correct order (stored in sub_categories)
	sub_categories <- .getNamesOfEachCategory(obj)
	for( i in 1:length( dmnms0 ) ){	dmnms0[[i]] <- intersect( sub_categories[[i]], dmnms0[[i]] ) }

 	M <- array( NA, dim=c(unlist(sapply(dmnms0, length))), dimnames=dmnms0)
 	PVAL <- ODDSRATIO <- M
 	if(NCATS==2){
		for( i in 1:length( split_all_lab ) ){
			ii <- split_all_lab[[i]]
			if( is.null(ResultTest[[ all_lab[i] ]][[ 'p.value' ]][1]) ){
				PVAL[ ii[1], ii[2] ] <- NA
				ODDSRATIO[ ii[1], ii[2] ] <- NA
			}else{
				PVAL[ ii[1], ii[2] ] <- ResultTest[[ all_lab[i] ]][[ 'p.value' ]][1]
				ODDSRATIO[ ii[1], ii[2] ] <- ResultTest[[ all_lab[i] ]][[ 'estimate' ]][1]
			}
		}
 	}
 	if(NCATS==3){
		for( i in 1:length( split_all_lab ) ){
			ii <- split_all_lab[[i]]
			if( is.null(ResultTest[[ all_lab[i] ]][[ 'p.value' ]][1]) ){
				PVAL[ ii[1], ii[2], ii[3] ] <- NA
				ODDSRATIO[ ii[1], ii[2], ii[3] ] <- NA
			}else{
				PVAL[ ii[1], ii[2], ii[3] ] <- ResultTest[[ all_lab[i] ]][[ 'p.value' ]][1]
				ODDSRATIO[ ii[1], ii[2], ii[3] ] <- ResultTest[[ all_lab[i] ]][[ 'estimate' ]][1]
			}
		}
 	}

 	rm(M); rm(split_all_lab)

	.checkValidOddsRatio <- function(x){
		BOOL <- rep(FALSE, times=length(x))
		for( i in 1:length(x) ){
			if( is.na(x[i]) )	{ BOOL[i] <- TRUE; next } #or==?
			if( is.nan(x[i]) )	{ BOOL[i] <- TRUE; next } #or==0/0
			if( x[i] == 0 )		{ BOOL[i] <- TRUE; next } #or==0/a
			if( x[i] == Inf )	{ BOOL[i] <- TRUE; next } #or==a/0
		}
	return(BOOL)
	}
	.printDim <- function(){ return( paste( paste(names(dimnames(ODDSRATIO)), dim(ODDSRATIO), sep="="), collapse=", ") ) }

	if(verbose>0){ message( "Dimension before filtering: ", .printDim(), sep="" ) }

	
	###
	### adjust p-values for multiple testing
	###
	if( !(p.adjust.method %in% p.adjust.methods) ){p.adjust.method <- 'none'}
	PVAL <- array(p.adjust(PVAL, method=p.adjust.method), dim=dim(PVAL), dimnames=dimnames(PVAL))


	###
	### drop slices with wrong direction
	###
	tmp <- .filter_drop.wrongdir.layer( filter, NCATS, ODDSRATIO, PVAL, nms_categories )
	if(is.null(tmp)){ return(obj) }
	ODDSRATIO <- tmp[[1]]; PVAL <- tmp[[2]]
	if(verbose>1){ message( "Dimension after filtering for direction: ", .printDim(), "\n", sep="" ) }

    
	###
	### drop slices
	###
	tmp <- .filter_dropLayer( filter, NCATS, ODDSRATIO, PVAL, nms_categories )
	if(is.null(tmp)){ return(obj) }
	ODDSRATIO <- tmp[[1]]; PVAL <- tmp[[2]]
    if(verbose>1){ message( "Dimension after filtering for P-value and low abs(log2(or)): ", .printDim(), "\n", sep="" ) }


    ###
	### drop slices with zero/NaN/Inf-valued odds ratio
	###
#     if(verbose>1){ message( "Removing uninformative layers\n", sep="" ) }
# 	for(d in 1:NCATS){
# 		keep <- apply(ODDSRATIO, d, function(x){ if( all(.checkValidOddsRatio(c(x))) ){return(FALSE)}else{return(TRUE)} })
# 		ODDSRATIO <- .keepLayer(ODDSRATIO, keep, NCATS, d)
# 		PVAL <- .keepLayer(PVAL, keep, NCATS, d)
# 	}

	if(verbose>0){ message( "Dimension after filtering: ", .printDim(), sep="" ) }
	obj@test.result.filter <- list( odds.ratio=ODDSRATIO, p.value=PVAL )
return(obj)
}

plotConCub <- function(obj, filter, fix.cat=1, show=list(), dontshow=list(), args_heatmap.2=list(), col=list(range=NULL), alt.names=list(), t=FALSE){


	NCATS <- length(obj@categories)
	if(is.numeric(fix.cat)){fix.cat <- names(obj@categories)[ fix.cat ]}

	nms_categories <- names(obj@categories)
	N_factor <- .getNumberOfFactorLevels(obj)

	ResultTest <- obj@test.result.filter
	if( length(ResultTest) == 0 ){ warning(paste("Call filterConCub before plotting", sep="")); return(NULL) }

	ODDSRATIO_0 <- ResultTest[['odds.ratio']]
	PVAL_0 <- ResultTest[['p.value']]

	if( length(show) > 0 ){
		arr.ind <- lapply(dimnames(ODDSRATIO_0), function(x){return(setNames(rep(FALSE, times=length(x)), x))} )
		for( i in 1:NCATS ){
			tmp <- show[[ nms_categories[i] ]]
			if( !is.null( tmp ) && length(tmp) > 0 ){
				arr.ind2 <- arr.ind; arr.ind2[[ nms_categories[i] ]][ which( names(arr.ind[[ nms_categories[i] ]]) %in% tmp ) ] <- TRUE
				if( NCATS == 2 ){ ODDSRATIO_0 <- ODDSRATIO_0[arr.ind2[[1]], arr.ind2[[2]], drop=FALSE]; PVAL_0 <- PVAL_0[arr.ind2[[1]], arr.ind2[[2]], drop=FALSE] }
				if( NCATS == 3 ){ ODDSRATIO_0 <- ODDSRATIO_0[arr.ind2[[1]], arr.ind2[[2]], arr.ind2[[3]], drop=FALSE]; PVAL_0 <- PVAL_0[arr.ind2[[1]], arr.ind2[[2]], arr.ind2[[3]], drop=FALSE] }
			}
		}
	}
	if( length(dontshow) > 0 ){
		POS <- setNames(vector("list", NCATS), names(dimnames(ODDSRATIO_0)))
		for( i in 1:NCATS ){
			pos <- which( dimnames(ODDSRATIO_0)[[i]] %in% dontshow[[ nms_categories[i] ]] )
			if( length(pos) > 0 ){ POS[[nms_categories[i]]] <- pos }else{ POS[[nms_categories[i]]] <- length(dimnames(ODDSRATIO_0)[[i]])+1 }
		}
		if( NCATS == 2 ){ ODDSRATIO_0 <- ODDSRATIO_0[-POS[[1]], -POS[[2]], drop=FALSE]; PVAL_0 <- PVAL_0[-POS[[1]], -POS[[2]], drop=FALSE] }
		if( NCATS == 3 ){ ODDSRATIO_0 <- ODDSRATIO_0[-POS[[1]], -POS[[2]], -POS[[3]], drop=FALSE]; PVAL_0 <- PVAL_0[-POS[[1]], -POS[[2]], -POS[[3]], drop=FALSE] }
	}

	COL_RANGE <- col$range
	tmp <- log2(ODDSRATIO_0); tmp[which( is.nan(tmp) | tmp == -Inf | tmp == Inf | is.na(tmp) )] <- 0
	if( is.null(col$range) ){ COL_RANGE <- range( c(tmp), na.rm=TRUE ); COL_RANGE <- c( -max(abs(COL_RANGE)), +max(abs(COL_RANGE)) ); }
	COL_LEVELS <- seq( COL_RANGE[1], COL_RANGE[2], length.out=1000 )
	COLORS <- c(colorpanel(floor(length(c(COL_LEVELS))/2), low='violet', mid='blue', high='white'), gplots::colorpanel(ceiling(length(COL_LEVELS)/2), low='white', mid='orange', high='red'))

	CollectPlots <- list()
	for( cnt in 1:length(obj@categories[[fix.cat]]) ){

		if(NCATS==2){
			ODDSRATIO_1 <- ODDSRATIO_0;PVAL_1 <- PVAL_0
		}
		if(NCATS==3){
			FixCounter <- names(obj@categories[[fix.cat]])[cnt]

			if( !(FixCounter %in% dimnames(ODDSRATIO_0)[[fix.cat]]) && !(FixCounter %in% names(alt.names[[fix.cat]]))){ #continue in case that current factor level was in 'dontshow'
				next;
			}
			ODDSRATIO_1 <- ODDSRATIO_0[FixCounter, , ];PVAL_1 <- PVAL_0[FixCounter, , ]
			if( fix.cat==nms_categories[2] ){ODDSRATIO_1 <- ODDSRATIO_0[, FixCounter, ]; PAVL_1 <- PVAL_0[, FixCounter, ]};
			if( fix.cat==nms_categories[3] ){ODDSRATIO_1 <- ODDSRATIO_0[, , FixCounter]; PAVL_1 <- PVAL_0[, , FixCounter]};
		}
		# don't show NA-rows and NA-cols; they are not removed by 'filter'-method when dontshow-option is used
		if( !obj@keep.empty.vars[[ names(dimnames(ODDSRATIO_1))[1] ]] ){
			keep <- apply( ODDSRATIO_1, 1, function(x){ !all(is.na(x)) } ); ODDSRATIO_1 <- ODDSRATIO_1[keep, , drop=FALSE]; PVAL_1 <- PVAL_1[keep, ,drop=FALSE]
		}
		if( !obj@keep.empty.vars[[ names(dimnames(ODDSRATIO_1))[2] ]] ){
			keep <- apply( ODDSRATIO_1, 2, function(x){ !all(is.na(x)) } ); ODDSRATIO_1 <- ODDSRATIO_1[, keep, drop=FALSE]; PVAL_1 <- PVAL_1[, keep, drop=FALSE]
		}
		tmp <- .filter_dropLayer( filter, 2, ODDSRATIO_1, PVAL_1, names(dimnames(ODDSRATIO_1)) )
		ODDSRATIO_1 <- tmp[[1]]; PVAL_1 <- tmp[[2]]
		rm(tmp, keep)

		if( length(alt.names) > 0 ){
			for( i in 1:length(alt.names) ){
				if( is.null(alt.names[[i]]) || !(names(alt.names)[i] %in% nms_categories) ){next;}
				ODDSRATIO_1 <- .translate(ODDSRATIO_1, margin=names(alt.names)[i], translation=alt.names[[i]])
				PVAL_1 <- .translate(PVAL_1, margin=names(alt.names)[i], translation=alt.names[[i]])
			}
		}


		if(t==TRUE){ ODDSRATIO_1 <- t(ODDSRATIO_1); PVAL_1 <- t(PVAL_1) }

		PVALSTAR <- pval2star(PVAL_1)
		ODDSRATIO_2 <- log2(ODDSRATIO_1); ODDSRATIO_2[which( is.nan(ODDSRATIO_2) | ODDSRATIO_2 == -Inf | ODDSRATIO_2 == Inf | is.na(ODDSRATIO_2) )] <- 0

		if( test.direction(filter) %in% c('over', "greater") ){ PVALSTAR[ ODDSRATIO_2 < 0 ] <- "" }
		if( test.direction(filter) %in% c('under', "less") ){ PVALSTAR[ ODDSRATIO_2 > 0 ] <- "" }

		#posMin <- apply( ODDSRATIO_2, 1, function(x){ if( abs(min(x)) > abs(max(x)) ){return(which.min(x)[1])}else{return(which.max(x)[1])} } )

		cexRow <- 1/log10(nrow(ODDSRATIO_2))
		my_args_heatmap.2 <- list( x=ODDSRATIO_2
			, Rowv=NA, Colv=NA, dendrogram='none', main='', revC=TRUE
			, margins=c(11,min(15, max(nchar(rownames(ODDSRATIO_2))))), symm=FALSE, cexRow=cexRow, cexCol=1, scale='none'
			, trace='none', col=COLORS
			, cellnote=PVALSTAR, notecol='black', notecex=cexRow
			, na.color="black"
			, symbreaks=TRUE
			)
		isct <- intersect(names(args_heatmap.2), names(my_args_heatmap.2))
		sdf1 <- setdiff(names(args_heatmap.2), names(my_args_heatmap.2))
		sdf2 <- setdiff(names(my_args_heatmap.2), names(args_heatmap.2))
		local_args_heatmap.2 <- c(args_heatmap.2, my_args_heatmap.2[sdf2])


		lab <- paste(fix.cat, names(obj@categories[[fix.cat]])[cnt], sep=", ")
		hm2 <- try({
		do.call("heatmap.2", local_args_heatmap.2)
		title(main=ifelse(NCATS==3, paste(fix.cat, names(obj@categories[[fix.cat]])[cnt], sep=", "), ""))
		})
		CollectPlots[[ lab ]] <- list(heat=hm2, p.value=PVAL_1, log2.odds.ratio=ODDSRATIO_2)
		
		if(NCATS==2){break;}
	}

	obj@test.result.filter.heatmap <- CollectPlots

invisible(obj)
}


