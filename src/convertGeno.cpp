#include <string>
#include <algorithm>
#include <vector>
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
    throw std::invalid_argument("GT format not found. Please check your VCF file format.");
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
  String ba = "B/A"; 
  for(int i=8; i<nc; i++) {
    if(i==8){
      List cv = df[i] ;
      std::string sub2 = cv[0] ;
      char* c = const_cast<char*>(sub2.c_str());
      try {
        idxGT = indexGT(c);
      } catch (std::invalid_argument& e) {
        Rprintf("%s\n",e.what());
        return false;
      }
      StringVector cv8 = df[i];
      for(int k=0;k<nr;k++){
        cv8[k] = "GT";
        mafCnt[k][0] = 0;
        mafCnt[k][1] = 0;
        mafCnt[k][2] = 0;
      }
    }
    
    if(i>8){
      StringVector cv = df[i];
      for(int k=0;k<nr;k++){
        char *pch;
        pch = strtok(cv[k],":");
        int l=0;
        char *str_tokens[15];
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
    }
  }
  // for(int k=0;k<nr;k++) {
  //   // maf[k] = (mafCnt[k][0]*2+mafCnt[k][1])/((mafCnt[k][0]+mafCnt[k][1]+mafCnt[k][2])*2)>0.5;
  //   maf.at(k) = (mafCnt[k][0]*2+mafCnt[k][1])/((mafCnt[k][0]+mafCnt[k][1]+mafCnt[k][2])*2)>0.5
  // }
  for(int i=9; i<nc; i++) {
    StringVector cv = df[i];
    for(int k=0;k<cv.size();k++){
      // if(maf[k]){
      if((mafCnt[k][0]*2+mafCnt[k][1])/((mafCnt[k][0]+mafCnt[k][1]+mafCnt[k][2])*2)>0.5){
        ab[k] = ba;
        if(cv[k]=="0") {
          cv[k]="2";
        } else if (cv[k]=="2"){
          cv[k]="0";
        }
      }
    }
  }
  return true;
}
