
.getLoadedList <- function( dbase, call ){

	if( is.list(dbase) ){ return(dbase) }	# dbase is list (keys=GO, value=genes)
	if( is.data.frame(dbase) || is.matrix(dbase) ){ 		# dbase is matrix or dataframe (col1=GO, col2=genes)
		loadedList <- split( as.character(dbase[, 2]), factor(as.character(dbase[, 1])) )
		return(loadedList)
	}
	if( length(grep( class(dbase), pattern="bimap", ignore.case=TRUE )) > 0 ){ # dbase is Bimap
			if(!require(AnnotationDbi)){ warning(paste0("Library ", sQuote("AnnotationDbi"), " cannot be loaded. Returning empty list.")); return(list()) }
			loadedList <- AnnotationDbi::as.list(dbase)
			return(loadedList)
	}
	
	warning(paste0("Data structure of passed object dbase=", call, " not supported. Returning empty list.")); return(list())
}

GO2list <- function(dbase, go.cat=NULL, rm=NULL, keep=NULL){

	loadedList <- .getLoadedList(dbase, call=deparse(substitute(dbase)))

	if(!is.null(go.cat)){ #remove GO categories
		if(!require(GO.db)){ warning(paste0("Library ", sQuote("GO.db"), " cannot be loaded. Returning empty list.")); return(list()) }
		my_keep <- unlist(sapply(names(loadedList), function( l ){ ifelse( ( Ontology( l ) %in% toupper(go.cat) ), TRUE, FALSE) }))
		loadedList <- loadedList[ my_keep ]
	}
	if(!is.null(rm)){ #remove GO terms
		loadedList <- loadedList[ setdiff( names(loadedList), toupper(rm) ) ]
	}
	if(!is.null(keep)){ #keep GO terms
		loadedList <- loadedList[ intersect( names(loadedList), toupper(keep) ) ]
	}
return(loadedList)
}

KEGG2list <- function(dbase, rm=NULL, keep=NULL){

	loadedList <- .getLoadedList(dbase, call=deparse(substitute(dbase)))

	if(!is.null(rm)){ #remove KEGG pathways
		loadedList <- loadedList[ setdiff( toupper(names(loadedList)), toupper(rm) ) ]
	}
	if(!is.null(keep)){ #keep KEGG pathways
		loadedList <- loadedList[ intersect( toupper(names(loadedList)), toupper(keep) ) ]
	}

return(loadedList)
}





.translate <- function(x, margin=2, translation=NULL){
	if(is.null(translation)){return(x)}
	dimnames(x)[[margin]] <- translation[dimnames(x)[[margin]]]
return(x)
}



######
### Formula
###

.checkFormula <- function(form){
	ch2 <- paste(as.character(form), collapse=" ")
	if( !(class(form) == "formula") ){ warning(paste( "", ch2, " has to be a valid formula" )); return(FALSE); }
return(TRUE)
}
###
### Formula
######

#####
## detect null-model

## mutual independence
## (X, Y, Z) =^= X+Y+Z			-- strsplit('+').len == 3

## single pairwise, joint independence		-- strsplit('+').len == 4
## (X, (YZ)) =^= X+Y+Z+Y*Z == X+Y+Z+Y:Z
## (Y, (XZ)) =^= X+Y+Z+X*Z == X+Y+Z+X:Z
## (Z, (XY)) =^= X+Y+Z+X*Y == X+Y+Z+X:Y

## two pairwise, conditional independence		-- strsplit('+').len == 5
## ((YX), (ZX)) =^= X+Y+Z+X*Y+X*Z == X+Y+Z+X:Y+X:Z	Y independent from Z given X
## ((XY), (ZY)) =^= X+Y+Z+X*Y+Y*Z == X+Y+Z+X:Y+Y:Z	X independent from Z given Y
## ((XZ), (YZ)) =^= X+Y+Z+X*Z+Y*Z == X+Y+Z+X:Z+Y:Z	Y independent from X given Z

## homogeneous association, no three-way interaction	-- strsplit('+').len == 6
## ((XY), (XZ), (YZ)) =^= X+Y+Z+X*Y+X*Z+Y*Z == X+Y+Z+X:Y+X:Z+Y:Z

## saturated model			-- strsplit('+').len == 7
## ((XYZ)) =^= X*Y*Z == X+Y+Z+X:Z+Y:Z+Y:X+Y:Z:X


# nms: all factor names
.getTypeOf_transformTable <- function(form= ~ factor1 + factor2 + factor3, nms=c('factor1','factor2','factor3')){
	form <- update(form, form)
	len <- length(nms)

	ch <- as.character( form )
	splt <- strsplit(ch[2], split=" \\+ ")[[1]] #split independent factors

	ttt <- ""; ordered_names <- nms
	if( ch[1] == "~" ){
		l <- length(splt)
		if( l == len & length(intersect(splt, nms)) == len ){ ttt <- "mi"; ordered_names <- nms }
		else{
			if( l == len+1 ){ #single pairwise (joint independence), single in front because of 'update'
				splt2 <- strsplit( splt[len+1], split=":" )[[1]] # split joint factors
				if( length(splt2) == len-1 ){
					sdf <- setdiff(nms, splt2)
					if( sdf == nms[1] ){ ttt <- "sp.1"; ordered_names <- c(nms[1], nms[c(2,3)]) }
					if( sdf == nms[2] ){ ttt <- "sp.2"; ordered_names <- c(nms[2], nms[c(1,3)]) }
					if( sdf == nms[3] ){ ttt <- "sp.3"; ordered_names <- c(nms[3], nms[c(1,2)]) }
				}
			}else{ #no exact hypergeom test

				if( l == len+2 ){ #two pairwise (conditional independence)
					ttt <- "ci"
				}
				if( l == len+3 ){ #homogeneous association
					ttt <- "ha"
				}
				if( l == len+4 ){ #saturated model
					ttt <- "sm"
				}
				if( ttt == "" ){ttt <- "any"}
				#message( paste0( "Currently comparisons of type ", form[2], " are not supported by exact hypergeometric test. Running approximate test." ) )
			}
		}
	}else{
		stop( paste( "Formula ", paste(ch, collapse="~"), " can not be resolved" ) )
	}
return(list(ttt, ordered_names))
}

# x is the list returned from .getTypeOf_transformTable()
.transformTable <- function(CT, x){
	CT_tmp <- matrix(NA, nrow=2, ncol=4, dimnames=list())
	type <- x[[1]]
	ordered_names <- x[[2]]
	ct <- c()
	if( type == "mi" ){ ct <- CT }
	else{
		if( type %in% c("sp.1", "sp.2", "sp.3")){
			CT_tmp <- as.matrix(ftable(CT, row.vars=ordered_names[1], col.vars=ordered_names[2:3]))
			colnames(CT_tmp) <- gsub(colnames(CT_tmp), pattern="_", replacement=":")
			names(dimnames(CT_tmp))[2] <- gsub(names(dimnames(CT_tmp))[2], pattern="_", replacement=":")
			ct <- CT_tmp
		}else{ct <- CT}
	}
return(ct)
}
##############



.pval2star <- function(x){
	y <- matrix("", nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
	Thr <- c(0.0001, 0.001, 0.01, 0.05, 0.1)
	Sym <- c('****', '***', '**', '*', '.')

	alreadyFound <- c()
	for(i in 1:5){
		pos <- which( x <= Thr[i] )
		if(length(pos)>0){
			pos0 <- setdiff(pos, alreadyFound)
			if(length(pos0)>0){y[pos0] <- Sym[i]}
			alreadyFound <- c(alreadyFound, pos0)
		}
	}
return(y)
}

