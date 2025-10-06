#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector altern_Craw(NumericVector x, int k,
    NumericVector i0, NumericVector i1){

    int n = x.size();
    int n0 = i0.size();
    int n1 = i1.size();
    NumericVector cts(2*k+3);

    for(int i = 0; i < 2*k+3; i++){
      cts[i] = 0;
    }
    cts[1] = 1;
    cts[2] = 1;
    for(int i = 0; i < n; i++){
      if(x[i]==-1){
        for(int j = 0; j < n0; j++){
          cts[i0[j]+3-1] = cts[i0[j]+3-1] + cts[i0[j]+3-2-1];
        }
        } else {
        for(int j = 0; j < n1; j++){
          cts[i1[j]+3-1] = cts[i1[j]+3-1] + cts[i1[j]+3-2-1];
        }
        }
    }
    return(cts);
    }
 
// [[Rcpp::export]]   
NumericMatrix altern_Craw2(NumericVector x, int k){

    int n = x.size();
    NumericMatrix A(k+1,2);

    A(0,0) = 1;
    A(0,1) = 1;
    for(int i = 0; i < n; i++){
      if(x[i]==-1){
        for(int j = 1; j < k+1; j++){
          A(j,0) = A(j,0) + A(j-1,1);
        }
        } else {
        for(int j = 1; j < k+1; j++){
          A(j,1) = A(j,1) + A(j-1,0);
        }
        }
    }
    return(A);
    }
    
// [[Rcpp::export]]
double compKdepth(NumericVector x, int K){

    int n = x.size();
    double cts;
    int ci;

    cts = 0;
    for(int i = 0; i < n-K; i++){
      ci = 0;
      if(x[i]!=x[i+1]){
        ci = 1;
        }
      while((ci>0)&(ci<K)){
      if(x[i+ci]!=x[i+ci+1]){
        ci++;
        } else {
        ci = 0;
        }
      }
      if(ci==K) cts++;
    }
    return(cts);
    }

// [[Rcpp::export]]    
NumericVector kSD_Craw2(NumericVector ord, 
    NumericVector sgn, int k, bool echo){

    int n = ord.size();
    NumericVector curo(3), curs(3);
    NumericVector D(2);
    double temp;
    bool cont;

    cont = true;
    for(int i = 0; i < n-2; i++){
      if(cont){
      curs[0] = sgn[i];
      curs[1] = sgn[i+1];
      curs[2] = sgn[i+2];
      curo[0] = ord[i];
      curo[1] = ord[i+1];
      curo[2] = ord[i+2];
      if(curo[0]>curo[1]){
        temp = curs[0];
        curs[0] = curs[1];
        curs[1] = temp;
        temp = curo[0];
        curo[0] = curo[1];
        curo[1] = temp;
      }
      if(curo[2]<curo[1]){
        temp = curs[1];
        curs[1] = curs[2];
        curs[2] = temp;
        temp = curo[1];
        curo[1] = curo[2];
        curo[2] = temp;      
      }
      if(curo[1]<curo[0]){
        temp = curs[0];
        curs[0] = curs[1];
        curs[1] = temp;
        temp = curo[0];
        curo[0] = curo[1];
        curo[1] = temp;      
      }
      if(((curs[0]==-1)&(curs[1]==1)&(curs[2]==-1))|
      ((curs[0]==1)&(curs[1]==-1)&(curs[2]==1))){
        D[0] = D[0] + 1;
        if(echo) Rcout << i+1 << "," << i+2 << "," << i+3 << std::endl;
        if((k==2)&(i<n-3)){
          curs[0] = sgn[i+1];
          curs[1] = sgn[i+2];
          curs[2] = sgn[i+3];
          curo[0] = ord[i+1];
          curo[1] = ord[i+2];
          curo[2] = ord[i+3];
          if(curo[0]>curo[1]){
            temp = curs[0];
            curs[0] = curs[1];
            curs[1] = temp;
            temp = curo[0];
            curo[0] = curo[1];
            curo[1] = temp;
          }
          if(curo[2]<curo[1]){
            temp = curs[1];
            curs[1] = curs[2];
            curs[2] = temp;
            temp = curo[1];
            curo[1] = curo[2];
            curo[2] = temp;      
          }
          if(curo[1]<curo[0]){
            temp = curs[0];
            curs[0] = curs[1];
            curs[1] = temp;
            temp = curo[0];
            curo[0] = curo[1];
            curo[1] = temp;      
          }
          if(((curs[0]==-1)&(curs[1]==1)&(curs[2]==-1))|
          ((curs[0]==1)&(curs[1]==-1)&(curs[2]==1))){
            D[1] = D[1] + 1;          
          } else {
           cont = false;
          }
        }
      }
      } else {
      cont = true;
      }
    }
    return(D);
    }

// [[Rcpp::export]] 
bool triangle(double si, double sj, double sk, double ri, double rj, double rk){
     if((si==sj)&(sj==sk)) return(false); else {
      if(si==sj){
      /* sk must be in the middle */
        if((ri<rj)&(ri<rk)&(rk<rj)){
          return(true); 
          } else { 
        if((ri>rj)&(ri>rk)&(rk>rj)) return(true);
          }
        return(false);
      }
      if(si==sk){
      /* sj must be in the middle */
        if((ri<rk)&(ri<rj)&(rj<rk)){
          return(true); 
          } else { 
        if((ri>rk)&(ri>rj)&(rj>rk)) return(true);
          }
        return(false);
      }
      if(sj==sk){
      /* si must be in the middle */
        if((rj<rk)&(rj<ri)&(ri<rk)){
          return(true); 
          } else { 
        if((rj>rk)&(rj>ri)&(ri>rk)) return(true);
          }
        return(false);
      }
     }
     return(false);
}

// [[Rcpp::export]]  
double kSD_Cfull2(NumericVector x, NumericVector rnk){

    int n = x.size();
    double cnt;
    
    cnt = 0;
    for(int i = 0; i<(n-3); i++){
      for(int j = i+1; j<(n-2); j++){
        for(int k = j+1; k<(n-1); k++){
          if(triangle(x[i],x[j],x[k],rnk[i],rnk[j],rnk[k])){
            for(int l = k+1; l<n; l++){
              if(triangle(x[j],x[k],x[l],rnk[j],rnk[k],rnk[l])){
              /* if(echo) Rcout << i+1 << "," << j+1 << "," << k+1 << 
              "," << l+1 << std::endl;  */
              cnt++;
              }
            }
          }
        }
      }
    }
    return cnt;
    }
    
// [[Rcpp::export]] 
bool testtriangle(double ri, double rj, double rk, double rl, int tp){
     /* types (tp) coding: 
     1 - 1101
     2 - 1011
     3 - 1010
     4 - 1001
     5 - 0101
     6 - 0110
     7 - 0010
     8 - 0100
     */
     int imin, imax;
     
     /* finding minimum and mamximum index */
     NumericVector r = {ri,rj,rk,rl};
     
     imin = which_min(r)+1; 
     imax = which_max(r)+1;
     
     /* if((ri<rj)&(ri<rk)&(ri<rl)) imin = 1; else {
     if((rj<ri)&(rj<rk)&(rj<rl)) imin = 2; else {
     if((rk<ri)&(rk<rj)&(rk<rl)) imin = 3; else imin = 4;
     }
     }
     if((ri>rj)&(ri>rk)&(ri>rl)) imax = 1; else {
     if((rj>ri)&(rj>rk)&(rj>rl)) imax = 2; else {
     if((rk>ri)&(rk>rj)&(rk>rl)) imax = 3; else imax = 4;
     }
     }*/
     
     /* depending on the type, testing the triangles */
     if(tp==1){ /* 1101 - possible only 12 */
         if(((imax==1)&(imin==2))|((imax==2)&(imin==1))) return(true); else return(false);      
     }
     if(tp==2){ /* 1011 - possible only 34 */
         if(((imax==3)&(imin==4))|((imax==4)&(imin==3))) return(true); else return(false);      
     }
     if(tp==3){ /* 1010 - possible only 14 */
         if(((imax==1)&(imin==4))|((imax==4)&(imin==1))) return(true); else return(false);      
     }          
     if(tp==4){ /* 1001 - possible only 23 */
         if(((imax==2)&(imin==3))|((imax==3)&(imin==2))) return(true); else return(false);      
     }
     if(tp==5){ /* 0101 - possible only 14 */
         if(((imax==1)&(imin==4))|((imax==4)&(imin==1))) return(true); else return(false);      
     }
     if(tp==6){ /* 0110 - possible only 23 */
         if(((imax==2)&(imin==3))|((imax==3)&(imin==2))) return(true); else return(false);      
     }
     if(tp==7){ /* 0010 - possible only 12 */
         if(((imax==1)&(imin==2))|((imax==2)&(imin==1))) return(true); else return(false);      
     }
     if(tp==8){ /* 0100 - possible only 34 */
         if(((imax==3)&(imin==4))|((imax==4)&(imin==3))) return(true); else return(false);      
     }
     /* Rcout << ri+1 << "," << rj+1 << "," << rk+1 << 
              "," << rl+1 << " ,type: " << tp << std::endl; */
     return(false);
}    
    
// [[Rcpp::export]]  
double kSD_Cfull(NumericVector x, NumericVector rnk){

    /* x is now a sequence of 0 and 1 that is ordered in their appearance on the line
       rnk is the index of each of the observations w.r.t. their time of observation
    */

    int n = x.size();
    NumericVector V1(n), W1(n);
    double n2 = n*(n-1)/2;
    double n3 = n*(n-1)*(n-2)/6;
    NumericMatrix V2(n2,2), V3(n2,2), W2(n2,2), W3(n2,2);
    NumericMatrix V4(n3,3), V5(n3,3), V6(n3,3), W4(n3,3), W5(n3,3), W6(n3,3);
    double iv1 = 0, iv2 = 0, iv3 = 0, iv4 = 0, iv5 = 0, iv6 = 0;
    double iw1 = 0, iw2 = 0, iw3 = 0, iw4 = 0, iw5 = 0, iw6 = 0;
    double cnt = 0;
    bool tst;
    
    for(int i = 0; i<n; i++){
      if(x[i]==1){  /* if + comes */
        /* put V1 to V2 (before adding 1 to V1) */
        for(int ii = 0; ii<iv1; ii++){
          V2(iv2,0) = V1(ii);
          V2(iv2,1) = i;
          iv2++;
        }
        /* add 1 to V1 */
        V1(iv1) = i;
        iv1++;
        /* put V3 to V5 (before putting W1 to V3) */
        for(int ii = 0; ii<iv3; ii++){
          V5(iv5,0) = V3(ii,0);
          V5(iv5,1) = V3(ii,1);
          V5(iv5,2) = i;
          iv5++;
        }    
        /* send V4 to test (before putting W2 to V4) */ 
        for(int ii = 0; ii<iv4; ii++){
          tst = testtriangle(rnk[V4(ii,0)], rnk[V4(ii,1)], rnk[V4(ii,2)], rnk[i], 2);
          if(tst) cnt++;
        }   
        /* put W1 to V3 */
        for(int ii = 0; ii<iw1; ii++){
          V3(iv3,0) = W1(ii);
          V3(iv3,1) = i;
          iv3++;
        }
        /* put W2 to V4 */
        for(int ii = 0; ii<iw2; ii++){
          V4(iv4,0) = W2(ii,0);
          V4(iv4,1) = W2(ii,1);
          V4(iv4,2) = i;
          iv4++;
        }
        /* put W3 to V6 */
        for(int ii = 0; ii<iw3; ii++){
          V6(iv6,0) = W3(ii,0);
          V6(iv6,1) = W3(ii,1);
          V6(iv6,2) = i;
          iv6++;
        }  
        /* send W4 to test */ 
        for(int ii = 0; ii<iw4; ii++){
          tst = testtriangle(rnk[W4(ii,0)], rnk[W4(ii,1)], rnk[W4(ii,2)], rnk[i], 1);
          if(tst) cnt++;
        }   
        /* send W5 to test */ 
        for(int ii = 0; ii<iw5; ii++){
          tst = testtriangle(rnk[W5(ii,0)], rnk[W5(ii,1)], rnk[W5(ii,2)], rnk[i], 4);
          if(tst) cnt++;
        }  
        /* send W6 to test */ 
        for(int ii = 0; ii<iw6; ii++){
          tst = testtriangle(rnk[W6(ii,0)], rnk[W6(ii,1)], rnk[W6(ii,2)], rnk[i], 5);
          if(tst) cnt++;
        }  
      } else {
      /* if - comes */
        /* put W1 to W3 (before adding 1 to W1) */
        for(int ii = 0; ii<iw1; ii++){
          W3(iw3,0) = W1(ii);
          W3(iw3,1) = i;
          iw3++;
        }
        /* add 1 to W1 */
        W1(iw1) = i;
        iw1++;
        /* put W2 to W5 (before putting V1 to W2) */
        for(int ii = 0; ii<iw2; ii++){
          W5(iw5,0) = W2(ii,0);
          W5(iw5,1) = W2(ii,1);
          W5(iw5,2) = i;
          iw5++;
        }    
        /* send W6 to test (before putting V3 to W6) */ 
        for(int ii = 0; ii<iw6; ii++){
          tst = testtriangle(rnk[W6(ii,0)], rnk[W6(ii,1)], rnk[W6(ii,2)], rnk[i], 8);
          if(tst) cnt++;
        }   
        /* put V1 to W2 */
        for(int ii = 0; ii<iv1; ii++){
          W2(iw2,0) = V1(ii);
          W2(iw2,1) = i;
          iw2++;
        }
        /* put V2 to W4 */
        for(int ii = 0; ii<iv2; ii++){
          W4(iw4,0) = V2(ii,0);
          W4(iw4,1) = V2(ii,1);
          W4(iw4,2) = i;
          iw4++;
        }
        /* put V3 to W6 */
        for(int ii = 0; ii<iv3; ii++){
          W6(iw6,0) = V3(ii,0);
          W6(iw6,1) = V3(ii,1);
          W6(iw6,2) = i;
          iw6++;
        }  
        /* send V4 to test */ 
        for(int ii = 0; ii<iv4; ii++){
        tst = testtriangle(rnk[V4(ii,0)], rnk[V4(ii,1)], rnk[V4(ii,2)], rnk[i], 3);
          if(tst) cnt++;
        }   
        /* send V5 to test */ 
        for(int ii = 0; ii<iv5; ii++){
          tst = testtriangle(rnk[V5(ii,0)], rnk[V5(ii,1)], rnk[V5(ii,2)], rnk[i], 6);
          if(tst) cnt++;
        }  
        /* send V6 to test */ 
        for(int ii = 0; ii<iv6; ii++){
          tst = testtriangle(rnk[V6(ii,0)], rnk[V6(ii,1)], rnk[V6(ii,2)], rnk[i], 7);
          if(tst) cnt++;
        }       
      }   
    }
    return(cnt);
    }

// [[Rcpp::export]]  
bool C123(double a, double b, double c){
  if((a<b)&(b<c)) return(true);
  if((a>b)&(b>c)) return(true);
  return(false);
}

// [[Rcpp::export]]  
bool C132(double a, double b, double c){
  if((a<c)&(c<b)) return(true);
  if((a>c)&(c>b)) return(true);
  return(false);
}

// [[Rcpp::export]]  
bool C312(double a, double b, double c){
  if((c<a)&(a<b)) return(true);
  if((c>a)&(a>b)) return(true);
  return(false);
}

// [[Rcpp::export]]  
NumericVector kSD_Cfull4(NumericVector s, NumericVector r){

    /* x is now a sequence of -1 and +1 that is ordered as in the original sample
       r is the rank of each element in the ordered sequence of +1 and -1
    */

    int m = s.size();
    NumericVector p(m), n(m); /* p for positive, n for negative */
    NumericMatrix pp(m*(m-1)/2,2), nn(m*(m-1)/2,2), np(m*(m-1)/2,2), pn(m*(m-1)/2,2);
    NumericMatrix ppn(m*(m-1)*(m-2)/6,3), nnp(m*(m-1)*(m-2)/6,3), pnp(m*(m-1)*(m-2)/6,3);
    NumericMatrix npn(m*(m-1)*(m-2)/6,3), pnn(m*(m-1)*(m-2)/6,3), npp(m*(m-1)*(m-2)/6,3);
    double ip = 0, in = 0, ipp = 0, inn = 0, inp = 0, ipn = 0;
    double ippn = 0, innp = 0, ipnp = 0, inpn = 0, ipnn = 0, inpp = 0;
    double cnt = 0;
    
    for(int i=0; i<m; i++){
      if(s[i]==1){
          if(innp>0){
            for(int j = 0; j<innp; j++){
              if(C312(r[nnp(j,1)],r[nnp(j,2)],r[i])){
                cnt++;
              }
            }
          }
          if(ipnp>0){
            for(int j = 0; j<ipnp; j++){
              if(C312(r[pnp(j,1)],r[pnp(j,2)],r[i])){
                cnt++;
              }
            }
          } 
          if(ippn>0){
            for(int j = 0; j<ippn; j++){
              if(C123(r[ppn(j,1)],r[ppn(j,2)],r[i])){
                cnt++;
              }
            }
          } 
          if(inpn>0){
            for(int j = 0; j<inpn; j++){
              if(C123(r[npn(j,1)],r[npn(j,2)],r[i])){
                cnt++;
              }
            }
          } 
          if(ipnn>0){
            for(int j = 0; j<ipnn; j++){
              if(C132(r[pnn(j,1)],r[pnn(j,2)],r[i])){
                cnt++;
              }
            }
          }           
          if(ip>0){
            for(int j = 0; j<ip; j++){
              pp(ipp,0) = p[j];
              pp(ipp,1) = i;
              ipp++;
              } /* add vector p to pp */
            }
          if(in>0){
            for(int j = 0; j<in; j++){
              np(inp,0) = n[j];
              np(inp,1) = i;
              inp++;
              } /* add vector n to np */
          }
          if(inn>0){
            for(int j = 0; j<inn; j++){
              if(C132(r[nn(j,0)],r[nn(j,1)],r[i])){
                nnp(innp,0) = nn(j,0);
                nnp(innp,1) = nn(j,1);
                nnp(innp,2) = i;
                innp++;
                }
              } /* add vector nn to nnp */
          }
          if(ipn>0){ 
            for(int j = 0; j<ipn; j++){
              if(C123(r[pn(j,0)],r[pn(j,1)],r[i])){
                pnp(ipnp,0) = pn(j,0);
                pnp(ipnp,1) = pn(j,1);
                pnp(ipnp,2) = i;
                ipnp++;
                }
              }  /* add vector pn to pnp */
          }
          if(inp>0){
            for(int j = 0; j<inp; j++){
              if(C312(r[np(j,0)],r[np(j,1)],r[i])){
                npp(inpp,0) = np(j,0);
                npp(inpp,1) = np(j,1);
                npp(inpp,2) = i;
                inpp++;
                }
              } /* add vector np to npp */
          }
          p[ip] = i;           
          ip++;  /* add i to p*/
      } else {
      /* s[i] == -1 */
          if(innp>0){
            for(int j = 0; j<innp; j++){
              if(C123(r[nnp(j,1)],r[nnp(j,2)],r[i])){
                cnt++;
              }
            }
          }
          if(ipnp>0){
            for(int j = 0; j<ipnp; j++){
              if(C123(r[pnp(j,1)],r[pnp(j,2)],r[i])){
                cnt++;
              }
            }
          } 
          if(ippn>0){
            for(int j = 0; j<ippn; j++){
              if(C312(r[ppn(j,1)],r[ppn(j,2)],r[i])){
                cnt++;
              }
            }
          } 
          if(inpn>0){
            for(int j = 0; j<inpn; j++){
              if(C312(r[npn(j,1)],r[npn(j,2)],r[i])){
                cnt++;
              }
            }
          } 
          if(inpp>0){
            for(int j = 0; j<inpp; j++){
              if(C132(r[npp(j,1)],r[npp(j,2)],r[i])){
                cnt++;
              }
            }
          }           
          if(ip>0){
            for(int j = 0; j<ip; j++){
              pn(ipn,0) = p[j];
              pn(ipn,1) = i;
              ipn++;
              } /* add vector p to pn */
            }
          if(in>0){
            for(int j = 0; j<in; j++){
              nn(inn,0) = n[j];
              nn(inn,1) = i;
              inn++;
              } /* add vector n to nn */
          }
          if(ipp>0){
            for(int j = 0; j<ipp; j++){
              if(C132(r[pp(j,0)],r[pp(j,1)],r[i])){
                ppn(ippn,0) = pp(j,0);
                ppn(ippn,1) = pp(j,1);
                ppn(ippn,2) = i;
                ippn++;
                }
              } /* add vector pp to ppn */
          }
          if(ipn>0){ 
            for(int j = 0; j<ipn; j++){
              if(C312(r[pn(j,0)],r[pn(j,1)],r[i])){
                pnn(ipnn,0) = pn(j,0);
                pnn(ipnn,1) = pn(j,1);
                pnn(ipnn,2) = i;
                ipnn++;
                }
              }  /* add vector pn to pnn */
          }
          if(inp>0){
            for(int j = 0; j<inp; j++){
              if(C123(r[np(j,0)],r[np(j,1)],r[i])){
                npn(inpn,0) = np(j,0);
                npn(inpn,1) = np(j,1);
                npn(inpn,2) = i;
                inpn++;
                }
              } /* add vector np to npn */
          }
          n[in] = i;           
          in++;  /* add i to n*/
      }
    }
    NumericVector res(2);
    res(0) =  innp+inpn+ipnn+ippn+ipnp+inpp;
    res(1) = cnt;
    return(res);
    }