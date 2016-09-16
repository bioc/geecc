#include <Rcpp.h>
#include <vector>


using namespace Rcpp;

// [[Rcpp::export()]]

	   StringVector cf_special1(const Rcpp::List &items_factor, const Rcpp::List &RNG, IntegerVector t){
		
		int N=0;
			Rcpp::StringVector aux( RNG[t[0]-1] );
			// http://gallery.rcpp.org/articles/subsetting/
			Rcpp::List aux2=as<Rcpp::List>( items_factor[aux] ); // get all form RNG[t]
			for( int i=0; i<aux2.size(); ++i ){ N += (as<StringVector>(aux2[i])).size(); }
			//Rcpp::Rcout << "N=" << N << std::endl;
			
			Rcpp::StringVector output(N);
			std::size_t index = 0;
			for( int i=0; i<aux2.size(); ++i ){ // http://stackoverflow.com/questions/30175104/how-to-effectively-combine-a-list-of-numericvectors-into-one-large-numericvector
                StringVector el=aux2[i];
                std::copy(el.begin(), el.end(), output.begin() + index);

      index += el.size();
            }
	   
	   return(output);
	   }
 

// [[Rcpp::export()]]
 
StringVector cf_intersect6(const StringVector &pop, const StringVector &x){
    int Lx=x.size(), Lpop=pop.size(), j=0;
   StringVector MED(Lx);
    for( int i=0; i<Lx; ++i ){
        int left=0, right=Lpop-1, med=0;
        while( left<=right ){
            med=floor((right + left) / 2);
            int compare=strcmp( x[i], pop[med] );
            if( compare==0 ){
                MED[j]=x[i]; ++j; break;
            }
            if( compare<0 ){right=med-1;}else{ left=med+1; }
        }
    }
return(StringVector(MED.begin(), MED.end()-(Lx-j)));
}

// [[Rcpp::export()]]

StringVector cf_setdiff1(const StringVector &pop, const StringVector &x){
    int Lx=x.size(), Lpop=pop.size();
    std::vector<bool> keep=std::vector<bool>(Lpop, true);

    int ndrop=0;
    for( int i=0; i<Lx; ++i ){
        int left=0, right=Lpop-1, med=0;
        while( left<=right ){
            med=floor((right + left) / 2);
            int compare=strcmp( x[i], pop[med] );
            if( compare==0 ){
                keep[med]=false; ++ndrop; break;
            }
            if( compare<0 ){right=med-1;}else{ left=med+1; }
        }
    }
    //Rprintf("Lpop=%i, ndrop=%i, diff=%i\\n", Lpop, ndrop, Lpop-ndrop);
   StringVector MED(Lpop-ndrop);    
    int j=0;
    for( int i=0; i<Lpop; ++i ){
      //  if( Lpop-i < 20 ){     Rprintf("i=%i, j=%i\\n", i, j); Rcpp::Rcout << pop[i] << std::endl;}
        if( keep[i] ){ 
            MED[j]=pop[i];
            ++j;
        }
    }
return(MED);
}


