

.getLoadedList <- function( dbase, call ){

	if( is.list(dbase) ){ return(dbase) }	# dbase is list (keys=GO, value=genes)
	if( is.data.frame(dbase) || is.matrix(dbase) ){ 		# dbase is matrix or dataframe (col1=GO, col2=genes)
		loadedList <- split( as.character(dbase[, 2]), factor(as.character(dbase[, 1])) )
		return(loadedList)
	}
	if( length(grep( class(dbase), pattern="bimap", ignore.case=TRUE )) > 0 ){ # dbase is Bimap
			if(!requireNamespace("AnnotationDbi")){ warning(paste0("Library ", sQuote("AnnotationDbi"), " cannot be loaded. Returning empty list.")); return(list()) }
			loadedList <- AnnotationDbi::as.list(dbase)
			return(loadedList)
	}
	
	warning(paste0("Data structure of passed object dbase=", call, " not supported. Returning empty list.")); return(list())
}

GO2list <- function(dbase, go.cat=NULL, rm=NULL, keep=NULL){

	loadedList <- .getLoadedList(dbase, call=deparse(substitute(dbase)))

	if(!is.null(go.cat)){ #remove GO categories
		if(!requireNamespace("GO.db")){ warning(paste0("Library ", sQuote("GO.db"), " cannot be loaded. Returning empty list.")); return(list()) }
		my_keep <- vapply(names(loadedList), function( l ){ ifelse( ( AnnotationDbi::Ontology( l ) %in% go.cat ), TRUE, FALSE) }, FUN.VALUE=TRUE )
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


get_gochildren <- function(nms){
	CAT <- c("BP", "CC", "MF")
}


GO2offspring <- function(x){
	if(!is.list(x)){ warning(paste0("Parameter ", sQuote("x"), " must be a list. Returning input unmodified.")); return(x) }
	if(!requireNamespace("GO.db")){ warning(paste0("Library ", sQuote("GO.db"), " cannot be loaded. Returning input list unmodified.")); return(x) }
	if(!requireNamespace("AnnotationDbi")){ warning(paste0("Library ", sQuote("AnnotationDbi"), " cannot be loaded. Returning input list unmodified.")); return(x) }
	
	R <- 1:length(x)
	nms <- names(x)
	bp <- AnnotationDbi::as.list(GO.db::GOBPOFFSPRING)[ nms ]
 	cc <- AnnotationDbi::as.list(GO.db::GOCCOFFSPRING)[ nms ]
 	mf <- AnnotationDbi::as.list(GO.db::GOMFOFFSPRING)[ nms ]
	res <- c( bp, cc, mf )
	off <- res[ !sapply(res, is.null) ]
	
	xx <- sapply( R, function( r ){ v <- unique( unlist(c( x[[ nms[r] ]], x[ off[[ nms[r] ]] ] )) ); return( v[!is.na(v)] ) }  )
	names(xx) <- nms
return(xx)
}


GO2level <- function(x, go.level=-1, relation=c("is_a")){
	if( !is.numeric(go.level) ){ warning(paste0(dQuote('go.level'), " needs to be ", sQuote("-1"), " or positive integer. Returning input list unmodified.")); return(x); }
	if( go.level==-1 || go.level==0 ){ return(x) }
	
	if(!requireNamespace("GO.db")){ warning(paste0("Library ", sQuote("GO.db"), " cannot be loaded. Returning input list unmodified.")); return(x) }
	if(!requireNamespace("AnnotationDbi")){ warning(paste0("Library ", sQuote("AnnotationDbi"), " cannot be loaded. Returning input list unmodified.")); return(x) }

	
	u <- unlist(AnnotationDbi::Term(GO.db::GOTERM))
	u2 <- u[grep("biological_process|molecular_function|cellular_component", u)]
	goid_gene_ontology <- "GO:0003673"
	roots_id2term <- u2
	roots_term2id <- setNames( names(u2), u2 )
	#print(u2)
	CAT <- c("BP", "CC", "MF")
	
	my_get <- function(id, type="CHILDREN"){
		res <- list()
		for( ct in CAT ){
			tmp <- AnnotationDbi::as.list( get(paste0("GO", ct, type)) )[ id ];
			res[[ ct ]] <- lapply(tmp, function(v){ v[names(v) %in% relation] })
		}
		res <- unique(unlist(res))
	return( res[ !sapply(res, is.null) ] )
	}

	
	allAncestors <- my_get(names(x), type="ANCESTOR");
	goid <- my_get( names(roots_id2term), type="CHILDREN" )

	next.level <- 2
	while( next.level<=go.level ){
		goid <- my_get( c(goid), type="CHILDREN" )
		next.level <- next.level+1
	}
	
	#summarize GO terms
	goid <- intersect(goid, allAncestors);
	off <- my_get( goid, type="OFFSPRING" )
	xx <- sapply( off, function( v ){ unique(unlist(x[ v ])) } )

return(xx)
}



sortAscii <- function(x){
    tmp <- Sys.getlocale("LC_COLLATE"); tmp2 <- Sys.setlocale("LC_COLLATE", "C"); x <- sort(x); tmp2 <- Sys.setlocale("LC_COLLATE", tmp);
return(x)
}
intersectPresort <- function(pop, x){ cf_intersect6(pop, x) }
# setdiffPresort <- function(pop, x){ cf_setdiff1(pop, x) }
# .special1 <- function( x,y ,t ){ cf_special1( x,y ,t ) }



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
## (X, (YZ)) =^= X+Y+Z+Y*Z == X+Y+Z+Y:Z == X+Y*Z
## (Y, (XZ)) =^= X+Y+Z+X*Z == X+Y+Z+X:Z == Y+X*Z
## (Z, (XY)) =^= X+Y+Z+X*Y == X+Y+Z+X:Y == Z+X*Y

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



pval2star <- function(x){
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


