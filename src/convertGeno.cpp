#include <string>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

int indexGT(char* c){
  char *pch;
  pch = strtok(c,":");
  int k=0, idxGT=-1;
  while (pch != NULL) {
    if(strstr(pch,"GT")!=NULL){
      idxGT = k;
    } else {
      k++;
    }
    pch = strtok(NULL,":");
  }
  if(idxGT==-1){
    throw std::invalid_argument("GT format not found.");
  }
  return idxGT;
}

int conversionFormat(char* geno){
  if(strstr(geno,"0/1")!=NULL || strstr(geno,"1/0")!=NULL || strstr(geno,"1|0")!=NULL || strstr(geno,"0|1")!=NULL) {
    return 1;
  } else if(strstr(geno,"1/1") || strstr(geno,"1|1") ) {
    return 2;
  } else if (strstr(geno,"0/0") || strstr(geno,"0|0") ) {
    return 0;
  } else {
    return 3;
  }
}


// [[Rcpp::export]]
bool convertGeno(DataFrame& df, StringVector& ab) {
  int nr = df.nrows(), nc= df.size(), idxGT;
  int mafCnt[nr][3];
  bool maf[nr];
  String ba = "B/A"; 
  // List out(nc) ;
  for(int i=8; i<nc; i++) {
    // if(i<8) {
    //   out[i] = df[i] ;
    // }
    if(i==8){
      List cv = df[i] ;
      std::string sub2 = cv[0] ;
      char* c = const_cast<char*>(sub2.c_str());
      // char *pch;
      // pch = strtok(c,":");
      // int k=0;
      // while (pch != NULL) {
      //   if(strstr(pch,"GT")!=NULL){
      //     idxGT = k;
      //   } else {
      //     k++;
      //   }
      //   pch = strtok(NULL,":");
      // }
      // if(idxGT==-1){
      //   return 1;
      // }
      try {
        idxGT = indexGT(c);
      } catch (std::invalid_argument& e) {
        fprintf(stdout,"%s\n",e.what());
        return false;
      }
      StringVector cv8 = df[i] ;
      for(int k=0;k<cv8.size();k++){
        cv8[k] = "GT";
        mafCnt[k][0] = 0;
        mafCnt[k][1] = 0;
        mafCnt[k][2] = 0;
      }
      // StringVector ins(cv.size(),"GT");
      // out[i] = ins;
    }
    
    if(i>8){
      StringVector cv = df[i];
      for(int k=0;k<cv.size();k++){
        char *pch;
        pch = strtok(cv[k],":");
        // pch = strtok(line,"\t");
        int l=0;
        char *str_tokens[k+1];
        while (pch != NULL)
        {
          str_tokens[l++]=pch;
          pch = strtok(NULL,":");
        }
        int genoINT = conversionFormat(str_tokens[idxGT]);
        if(genoINT==0){
          mafCnt[k][0]++;
        } else if(genoINT==1){
          mafCnt[k][1]++;
        } else if (genoINT==2) {
          mafCnt[k][2]++;
        }
        cv[k] = genoINT;
      }
      // cv = as<NumericVector>(cv);
      // double cnt0 = std::count(cv.begin(),cv.end(),0);
      // double cnt1 = std::count(cv.begin(),cv.end(),1);
      // double cnt2 = std::count(cv.begin(),cv.end(),2);
      // double maf = (cnt0*2+cnt1)/((cnt0+cnt1+cnt2)*2);
      // fprintf(stdout,"%f\n",maf);
      
      // out[i] = cv;
    }
  }
  
  for(int k=0;k<nr;k++) {
    maf[k] = (double(mafCnt[k][0])*2+double(mafCnt[k][1]))/((double(mafCnt[k][0])+double(mafCnt[k][1])+double(mafCnt[k][2]))*2)>0.5;
  }
  for(int i=9; i<nc; i++) {
    StringVector cv = df[i];
    for(int k=0;k<cv.size();k++){
      if(maf[k]){
        ab[k] = ba;
        // ab.insert(k,"B/A");
        if(cv[k]=="0") {
          cv[k]="2";
        } else if (cv[k]=="2"){
          cv[k]="0";
        }
      }
    }
  }
  
  // out.attr("class") = df.attr("class") ;
  // out.attr("row.names") = df.attr("row.names") ;
  return true;
}
