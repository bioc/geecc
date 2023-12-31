

setClass( "concub",
	representation(
		categories="list", population="character", keep.empty.vars="list", options="list"
		, approx="numeric", null.model="formula"
		, test.result="list", test.result.filter="list", test.result.filter.heatmap="list"
	)
)


setMethod("initialize", signature="concub"
	, definition=function(.Object, categories, population, keep.empty.vars, options, approx, null.model){

	test_result <- list()
	MAX_NUM_FACT_SUPPORT <- 3
	num_categories <- 0
	nms_categories <- ""

	if( missing(categories) ){ stop("List of categories (parameter ", sQuote("categories"), ") required.") }
	else{
		num_categories <- length(categories)
		if( num_categories > MAX_NUM_FACT_SUPPORT ){ warning("Only two or three categories are supported"); num_categories <- MAX_NUM_FACT_SUPPORT; }
		nms_categories <- names(categories);
		if( is.null(nms_categories) ){
			nms_categories <- paste0( "category", 1:min(num_categories, MAX_NUM_FACT_SUPPORT) ); names(categories) <- nms_categories
			message(paste0("Changed names of categories: ", paste(nms_categories, collapse=",")))
		}
		.Object@categories <- categories
	}
	
	tmp <- setNames( vector("list", num_categories), nms_categories ); tmp <- lapply(tmp, function(x){FALSE})
	if(missing(keep.empty.vars)){
		keep.empty.vars <- tmp 
	}else{
		keep.empty.vars[ setdiff( nms_categories, names(keep.empty.vars) ) ] <- FALSE
	}
	.Object@keep.empty.vars <- keep.empty.vars
	
	if(missing(population)){
        population <- sortAscii(unique(unlist( sapply( categories, function(x){ unique(unlist(x, use.names=FALSE)) } ), use.names=FALSE )));
	}else{
		if(is.character(population)){
            population <- sortAscii(unique(population));
			categories <- lapply( categories, function(x1){ lapply( x1, function( x ){ isct <- intersectPresort( population, x ); if(length(isct)==0){return(NULL)};return( isct ) } ) } )
			#categories <- lapply( categories, function(x1){ lapply( x1, function( x ){ isct <- intersect( population, x ); if(length(isct)==0){return(NULL)};return( isct ) } ) } )
			for( i in 1:num_categories ){
				ii <- nms_categories[i]
				# if declared and FALSE, then remove
				if( !keep.empty.vars[[ ii ]] ){
					categories[[ ii ]] <- categories[[ ii ]][ !sapply( categories[[ ii ]], is.null) ]
				} #otherwise keep everything
			}
			.Object@categories <- categories
		}else{
			if( class(population) %in% c("eSet", "ExpressionSet", "DGEList") ){ 
				cls <- class(population)
				population <- sortAscii( unique(rownames(population) ));
				if( cls == "eSet" ){ population <- population[ grep(population, pattern="^AFFX", invert=TRUE) ] }
			}
		}
	}
	.Object@population <- population
	
	
	# change order of categories for speed-up during computations
	# categories like GO-terms typically contain thousands of factor-levels, other ones like phylostrata 10-20
	tmp <- order( sapply( categories, length ), decreasing=FALSE );
	categories <- categories[ tmp ]
	if( num_categories==MAX_NUM_FACT_SUPPORT ){ categories <- categories[ c(1,3,2) ] } # swap
    if( any( names(categories) != names(.Object@categories) ) ){
        nms_categories <- names(categories)
        message("Changed order of categories: ", paste(nms_categories, collapse=","))
        .Object@categories <- categories 
    }


	if( missing(null.model) ){ null.model <- as.formula( paste("~", paste(nms_categories, collapse="+")) )
	}else{
		bool <- .checkFormula(null.model) # stop if invalid formula
		if( !bool ){return(NULL)}
	}
 	.Object@null.model <- update(null.model, null.model)
 	if( !missing(approx) ){ .Object@approx <- max(c(approx, 0)) }else{ .Object@approx <- 0 }

	default_factor_opt <- list( grouping=c("none", 'cumf', 'cumr', 'sw')[1], width=1, strat=FALSE )
	my_opt <- setNames(vector("list", num_categories), nms_categories)
	if( !missing(options) ){ #set grouping options for sets
		for( i in seq_len(num_categories) ){
			ii <- nms_categories[i]
			my_opt[[ ii ]] <- default_factor_opt
			for( nm in intersect(names(default_factor_opt), names(options[[ ii ]])) ){ my_opt[[ ii ]][[nm]] <- options[[ ii ]][[nm]] }
		}
		.Object@options <- my_opt
	}else{ my_opt <- lapply(my_opt, function(x){ default_factor_opt }); .Object@options <- my_opt }


	# set slots filled step by step at runtime to empty list
	.Object@test.result <- list()
	.Object@test.result.filter <- list()
	.Object@test.result.filter.heatmap <- list()

	validObject(.Object)
	return(.Object)
})

setMethod("show", "concub", function(object){
	.printsepline <- function(){ cat(paste( rep("#", times=20), collapse="" ), "\n", sep="") }

	cat("\n", sep="")
	.printsepline()
 	cat("# ", "settings", "\n", sep="")
	.printsepline()

	satmod <- paste0( "count ~ ", paste0(names(object@categories), collapse="*") )
	cat("Comparing null-model '", paste0("count ~ ", as.character(object@null.model)[2]), "' against alternative model '", satmod, "' \n", sep="")
	cat("Using chi-squared approximation"); if(object@approx>0){ cat(" unless expected value greater than ", object@approx, "\n", sep="") }
	cat("\n\n")

	x <- object@categories
	for( i in seq_along(x) ){
		Lxi <- length(x[[i]])
		cat("Category ", i, " (", names(x)[i], ") with ", Lxi, " variables\n", sep="")
		print(lapply(x[[i]][1:min(5, Lxi)], head))
		print(sapply(x[[i]][1:min(5, Lxi)], length))
		if( Lxi > 5 ){cat("[... output truncated after 5 items]\n", sep="")}
		cat("\n")
	}
	Lpop <- length(object@population)
	cat("Population provided or guessed from categories (", Lpop, " items):\n", sep="")
	print(head(object@population, 20))
	if(Lpop>20){cat("[... output truncated after 20 items]\n", sep="")}
	cat("\n\n")

# 	cat("")
# 	if( length(object@test.result) > 0 ){ #no sense to show this unless good formatation
# 		print(object@test.result[!sapply(object@test.result, is.null)][1:10])
# 	}
#
# 	if(length(object@test.result.filter)>0){
# 		print(head(object@test.result.filter[[1]]))
# 		print(head(object@test.result.filter[[2]]))
# 	}
})





setGeneric(name="getTable", def=function(object, na.rm=TRUE, dontshow=list()){ standardGeneric("getTable") })
setMethod(f="getTable", signature="concub",
	definition=function(object, na.rm=TRUE, dontshow=list()){

		if(is.null(object) || length(object@test.result.filter) == 0){ warning("Empty list in concub-object."); return(NULL); }

		items_categories <- .getItemsInEachCategory(object)
		len_sub_categories <- lapply(items_categories, function(x){sapply(x, length)})

		or <- object@test.result.filter[['odds.ratio']]
		pval <- object@test.result.filter[['p.value']]
		tmp <- dimnames(or)
		my_separator <- .my_separator()

		cat.names <- names(tmp)
		n.cat <- paste("n", names(tmp), sep=".")
		cl <- c(cat.names, n.cat, "n.tags", "p.value", "log2.odds.ratio", "tags")

		labs0 <- as.matrix(expand.grid(tmp))
		labs1 <- apply(labs0, 1, paste, collapse=my_separator)
		tab2 <- matrix(NA, nrow=length(labs1), ncol=length(cl), dimnames=list(labs1, cl))
		tab2 <- as.data.frame(tab2)
		n <- expand.grid(len_sub_categories);
		rownames(n) <- apply(expand.grid(lapply(object@categories, names)), 1, paste, collapse=my_separator)
		colnames(n) <- n.cat
		tab2[labs1, cat.names] <- labs0
		tab2[labs1, n.cat] <- n[labs1, ]

		tab2[labs1, 'p.value'] <- pval[ labs0 ]
		tab2[labs1, 'log2.odds.ratio'] <- log2(or[ labs0 ])
		tab2[labs1, 'n.tags'] <- sapply(object@test.result[ labs1 ], function(x){length(x$subpop)})
		tab2[labs1, 'tags'] <- sapply(object@test.result[ labs1 ], function(x){paste(x$subpop, collapse=",", sep=",")})


		rownames(tab2) <- NULL
		if( na.rm ){ na.pos <- which(is.na(tab2[, "p.value"])); if( length(na.pos)>0 ){ tab2 <- tab2[-na.pos, ] } }

		if( !is.null(dontshow) && length(dontshow) > 0 ){
			for( nm in names(dontshow) ){
				if( !is.null(dontshow[[nm]]) && length(dontshow[[nm]]) > 0 ){
					tab2 <- tab2[ (tab2[, nm] %in% dontshow[[nm]]), , drop=FALSE]
				}
			}
		}
		for( cln in names(tmp) ){ tab2[, cln] <- as.character(tab2[, cln]) }
		for( cln in c(paste("n", names(tmp), sep="."), 'n.tags') ){ tab2[, cln] <- as.integer(as.character(tab2[, cln])) }
		for( cln in c('p.value', 'log2.odds.ratio') ){ tab2[, cln] <- as.numeric(as.character(tab2[, cln])) }

return(tab2)
})

