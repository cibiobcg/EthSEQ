#include <Rcpp.h>
using namespace Rcpp;

// 
// [[Rcpp::export]]
List splitGenoDF(DataFrame& df) {
  int idxGT;
  int nr = df.nrows(), nc= df.size() ;
  List out(nc) ;
  for( int i=0; i<nc; i++) {
    if(i<8) {
      out[i] = df[i] ;
    }
    if(i==8){
      List cv = df[i] ;
      // char format[cv.size()][2];
      std::string sub2 = cv[0] ;
      char* c = const_cast<char*>(sub2.c_str());
      char *pch;
      pch = strtok(c,":");
      int k=0;
      while (pch != NULL) {
        if(strstr(pch,"GT")!=NULL){
          idxGT = k;
          break;
        } else {
          k++;
        }
        pch = strtok(NULL,":");
      }
      StringVector ins(cv.size(),"GT");
      out[i] = ins;
    }
    if(i>8){
      StringVector cv = df[i] ;
      for(int k=0;k<cv.size();k++){
            char *pch;
            pch = strtok(cv[k],":");
            // pch = strtok(line,"\t");
            int l=0;
            char *str_tokens[10];
            while (pch != NULL)
            {
              str_tokens[l++]=pch;
              pch = strtok(NULL,"\t");
            }
        cv[k] = str_tokens[idxGT];
      }
      out[i] = cv;
    }
  }
  out.attr("class") = df.attr("class") ;
  out.attr("row.names") = df.attr("row.names") ;
  return out;
}