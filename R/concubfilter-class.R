

setClass( Class="concubfilter"
	, representation=representation(
		nfact="numeric"
		, names="character"
		, p.value="numeric"
 		, test.direction="character"
 		, minimum.l2or="numeric"

 		, skip.min.group="numeric"
 		, skip.min.obs="numeric"
 		, skip.zeroobs="logical"

 		, drop.insignif.layer="logical"
 		, drop.wrongdir.layer="logical"
 		, drop.lowl2or.layer="logical"
 #TODO: , drop.min.obs.layer
		)
	, prototype=prototype(
		nfact=2
		, names=c("factor1", "factor2")
		, p.value=0.1
		, minimum.l2or=0.0
		, test.direction="two.sided"

		, skip.min.group=rep(0, times=2)
 		, skip.min.obs=1
 		, skip.zeroobs=TRUE

 		, drop.insignif.layer=c(FALSE, FALSE)
 		, drop.wrongdir.layer=c(FALSE, FALSE)
 		, drop.lowl2or.layer=c(FALSE, FALSE)
	)
	, validity=function(object){
		## set names for vectors
		slot_names <- slotNames(object)
		for( sn in c( "skip.min.group", "drop.insignif.layer", "drop.wrongdir.layer", "drop.lowl2or.layer" ) ){
			v <- slot(object, sn)
			if(is.null( attr(v, "names") )){
				stop(paste0("Missing names for concubfilter-slot ", dQuote(sn)))
				#slot(.Object, sn) <- setNames(slot(.Object, sn), .Object@names)
			}
		}
	return(TRUE)
	}
)

slot_names <- slotNames("concubfilter")

for( sn in setdiff(slot_names, c("names", "nfact")) ){
	dq0 <- paste0("\"", "concubfilter", "\"")
	dq1 <- paste0("\"", sn, "<-", "\"")
	dq2 <- paste0("\"", sn, "\"")
	str1 <- paste0( "setGeneric(", dq1, ", function(x, value){standardGeneric(", dq1, ")})" )
	str2 <- paste0( "setReplaceMethod(", dq2, ", ", dq0, ", function(x, value){ x@", sn, " <- value; validObject(x); return(x) })" )
	str3 <- paste0( "setGeneric(", dq2, ", function(x){standardGeneric(", dq2, ")})" )
	str4 <- paste0( "setMethod(", dq2, ", ", dq0, ", function(x){ return(x@", sn, ") })" )

	eval( parse( text=str1 ) )
	eval( parse( text=str2 ) )
	eval( parse( text=str3 ) )
	eval( parse( text=str4 ) )
}


# ## create Rd entry
# for( sn in setdiff(slot_names, c("names", "nfact")) ){
# 	cat(paste("\\alias{", sn, c("", "<-", ",concubfilter-method", "<-,concubfilter-method"), sep=""), sep="}\n")
# }


setMethod(f="initialize"
	, signature="concubfilter"
	, definition=function(.Object, names, p.value, test.direction, minimum.l2or
			, skip.min.group, skip.min.obs, skip.zeroobs
			, drop.insignif.layer, drop.wrongdir.layer, drop.lowl2or.layer){
		if( !missing(names) ){
			if( length(names) <= 1 || length(names) > 3){ stop(paste("Wrong number of names (names=", paste(collapse=",", names), ").", sep="")) }
			.Object@names <- names
			.Object@nfact <- length(names)
		}else{stop("names required")}
		
		if( !missing(p.value) ){
			if( !( p.value >= 0 && p.value <= 1 ) ){ stop(paste("p.value must be between 0 and 1.")) }
		}
		if( !missing(minimum.l2or) ){
			if( !( minimum.l2or < 0 ) ){ stop(paste("minimum.l2or must be non-negative.")) }
		}

		if( !missing(test.direction) ){
			if( !( test.direction %in% c("two.sided", "greater", "less") ) ){ stop(paste("Wrong setting for test.direction.")) }
		}

		##
		## skip test
		##
		if( !missing(skip.min.group) ){
			.Object@skip.min.group <- skip.min.group
		}else{ .Object@skip.min.group <- rep(0, times=.Object@nfact) }
		if( !missing(skip.min.obs) ){
			if( skip.min.obs < 0 ){ skip.min.obs <- as.integer(0); .Object@skip.min.obs <- skip.min.obs }
		}else{ .Object@skip.min.obs <- 1 }
		if( !missing(skip.zeroobs) ){
			if( !is.logical(skip.zeroobs) ){ stop(paste("skip.zeroobs must be TRUE or FALSE")) }
		}

		##
		## reduce amount of test results
		##
		if( !missing(drop.insignif.layer) ){
			print(drop.insignif.layer)
			print(.Object@drop.insignif.layer)
			.Object@drop.insignif.layer <- setNames(drop.insignif.layer, .Object@names)
			print(.Object@drop.insignif.layer)
		}else{ .Object@drop.insignif.layer <- rep(FALSE, times=.Object@nfact) }
		if( !missing(drop.wrongdir.layer) ){
			.Object@drop.wrongdir.layer <- setNames(drop.wrongdir.layer, .Object@names)
		}else{ .Object@drop.wrongdir.layer <- rep(FALSE, times=.Object@nfact) }
		if( !missing(drop.lowl2or.layer) ){
			.Object@drop.lowl2or.layer <- setNames(drop.lowl2or.layer, .Object@names)
		}else{ .Object@drop.lowl2or.layer <- rep(FALSE, times=.Object@nfact) }
		
		
		for( sn in c( "skip.min.group", "drop.insignif.layer", "drop.wrongdir.layer", "drop.lowl2or.layer" ) ){
			v <- slot(.Object, sn)
			if(is.null( attr(v, "names") )){slot(.Object, sn) <- setNames(slot(.Object, sn), .Object@names)}
		}

		validObject(.Object)
	return(.Object)
	}
)


setMethod("show", "concubfilter", function(object){
	.pasteNamedVector <- function(x){ paste( attr(x, "names"), x, sep="=", collapse="," ) }
	.printsepline <- function(){ cat(paste( rep("#", times=20), collapse="" ), "\n", sep="") }

	cat("\n", sep="")
	.printsepline()
 	cat("# ", "current filter settings", "\n", sep="")
	.printsepline()
	cat( "Number of categories: ", object@nfact, "\n", sep="" )
	cat( "Maximum P-value: ", object@p.value, "\n", sep="" )
	cat( "Minimum absolute log2 odds ratio: ", object@minimum.l2or, "\n", sep="" )
	cat( "Test direction: ", object@test.direction, "\n", sep="" )

	cat( "\nSkip test in case of", "\n", sep="" )
	cat( "\tno items: ", object@skip.zeroobs, "\n", sep="" )
	cat( "\tnumber of items less than: ", object@skip.min.obs, "\n", sep="" )
	cat( "\tminimum marginal (to skip small gene sets): ", .pasteNamedVector(object@skip.min.group), "\n", sep="" )

	cat( "\nLayers to be dropped in case of", "\n", sep="" )
	cat( "\tinsignificant P-values: ", .pasteNamedVector( object@drop.insignif.layer ), "\n", sep="" )
	cat( "\twrong direction: ", .pasteNamedVector( object@drop.wrongdir.layer), "\n", sep="" )
	cat( "\tsmall abs(log2(or)): ", .pasteNamedVector( object@drop.lowl2or.layer), "\n", sep="" )
	.printsepline()
})

