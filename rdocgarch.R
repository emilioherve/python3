
setwd("C:/Users/Emile Ndoumbe/Dropbox/emile1/docgarch")
PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
Sys.which("make")

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(tgp)
library(readr)
library(rugarch)
library(rmgarch)
library(xts)
library(scales)


###### Functions for Doc in mean and Doc in volatility

"matrix.sqrt.inv.pc" = function(A)
{
  evd = eigen(A, symmetric=TRUE)
  Asqrtinv = diag(1/sqrt(evd$values)) %*% t(evd$vectors)
  return(Asqrtinv)
}


"givens.rotation" = function(theta=0, d=2, which=c(1,2))
{
  # David S. Matteson
  # 2008.04.28
  # For a given angle theta, returns a d x d Givens rotation matrix
  #
  # Ex: for i < j , d = 2:  (c -s)
  #                         (s  c)
  c = cos(theta)
  s = sin(theta)
  M = diag(d)
  a = which[1]
  b = which[2]
  M[a,a] =  c
  M[b,b] =  c
  M[a,b] = -s
  M[b,a] =  s
  M
}

"theta2W" = function(theta)
{
  # David S. Matteson
  # 2010.02.17
  # For a vector of angles theta, returns W, a d x d Givens rotation matrix:
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2 
  ##  if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
  d = (sqrt(8*length(theta)+1)+1)/2
  if(d - floor(d) != 0){stop("theta must have length: d(d-1)/2")}
  W = diag(d)
  index = 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      Q.ij = givens.rotation(theta[index], d, c(i,j))
      W = Q.ij %*% W 
      index = index + 1
    }
  }
  W
}

src1 <- '
Rcpp::NumericVector theta_r(theta_s);  // creates Rcpp vector from SEXP
arma::colvec theta(theta_r.begin(), theta_r.size(), false);  // reuses memory and avoids extra copy

Rcpp::NumericMatrix Z_r(Z_s);                 // creates Rcpp matrix from SEXP
int n = Z_r.nrow(), d = Z_r.ncol();
arma::mat Z(Z_r.begin(), n, d, false);       // reuses memory and avoids extra copy

Rcpp::NumericVector L_r(L_s);                 // creates Rcpp matrix from SEXP
int L = as<int>(L_r);

Rcpp::NumericVector PHI_r(PHI_s);  // creates Rcpp vector from SEXP
arma::colvec PHI(PHI_r.begin(), PHI_r.size(), false);  // reuses memory and avoids extra copy

//int p = theta_r.size();

arma::mat Iden(d,d);
Iden.eye();

int k = 0; //index for theta
arma::mat W(d,d);
W.eye();

for(int j=1;j<=(d-1);j++){
	for(int i=(j+1);i<=d;i++){				
    arma::mat Q(d,d);
    Q.eye();
    double c = cos(theta[k]);
    double s = sin(theta[k]);	
    Q(i-1,i-1) = c;  
		Q(j-1,j-1) = c;
		Q(i-1,j-1) = -s;
    Q(j-1,i-1) = s;  
    W = Q * W;
		k++;		
	}//end inner for		
}//end outer for 

arma::mat tW = trans(W);
arma::mat S = Z * tW;
arma::mat ColMeans = mean(S, 0);

int z=0; //index for out
double ans= 0.0; 
int b = (d*(d-1))/2 + L*d*(d-1); //number of entries in output
arma::colvec out(b);

//Lag 0. A total of d(d-1)/2 entries in out will be filled. i and j are column numbers. k is row number 
for(int i=0;i<d;i++){	
	for(int j=i+1;j<d;j++){
		
		ans = 0.0; //next column so restart 
		ans = dot(S.col(i), S.col(j));							
		out[z] = ans/n - ColMeans(i)*ColMeans(j); //adjust by cross mean and store in out 			
		z++;							
	}//end j for
}//end i for  -- END Lag 0

//z now contains the index of the next position for outptr
/**Lag 1 and above. */
  
	for(int lag=1;lag<=L;lag++){		
		for(int i=0;i<d;i++){	
			for(int j=0;j<d;j++){			
				ans=0.0;
				
				if(i!=j){      
        // A.submat(p, r, q, s)	A.submat(first_row, first_col, last_row, last_col)       
        ans = dot(S.submat(0,i,n-lag-1,i), S.submat(lag,j,n-1,j));							
        out[z] = ans/n - ColMeans(i)*ColMeans(j); //adjust by cross mean and store in out 			
        z++;
					}else{ 
						//Do nothing				
					} //end if-else 			
			} //end j for	-- Finished column in 2nd matrix. Move to next one. 
		} //end i for -- Finished with column in 1st matrix. Move to next one. 
} //end lag for. Finished calculations for this lag. Next lag. 
	
double obj = dot((out % PHI), out);

return Rcpp::wrap(obj);
'

DOCinMean = cxxfunction(signature(theta_s="numeric", Z_s="matrix", L_s="integer", PHI_s="numeric"), 
                        src1, Rcpp=TRUE, verbose=TRUE, plugin = "RcppArmadillo",
                        libargs = c("-L/opt/local/lib -lgsl -lgslcblas -lm"), #-lgsl 
                        cppargs = c("-I/opt/local/include"), 
                        includes = c("#include <stdio.h>")
)


src <- '
Rcpp::NumericVector theta_r(theta_s);  // creates Rcpp vector from SEXP
arma::colvec theta(theta_r.begin(), theta_r.size(), false);  // reuses memory and avoids extra copy

Rcpp::NumericMatrix Z_r(Z_s);                 // creates Rcpp matrix from SEXP
int n = Z_r.nrow(), d = Z_r.ncol();
arma::mat Z(Z_r.begin(), n, d, false);       // reuses memory and avoids extra copy

Rcpp::NumericVector L_r(L_s);                 // creates Rcpp matrix from SEXP
int L = as<int>(L_r);

Rcpp::NumericVector PHI_r(PHI_s);  // creates Rcpp vector from SEXP
arma::colvec PHI(PHI_r.begin(), PHI_r.size(), false);  // reuses memory and avoids extra copy

Rcpp::NumericVector C_r(C_s);                 // creates Rcpp matrix from SEXP
double C = C_r(0);

arma::mat Iden(d,d);
Iden.eye();

int k = 0; //index for theta
arma::mat W(d,d);
W.eye();

for(int j=1;j<=(d-1);j++){
	for(int i=(j+1);i<=d;i++){				
    arma::mat Q(d,d);
    Q.eye();
    double c = cos(theta[k]);
    double s = sin(theta[k]);	
    Q(i-1,i-1) = c;  
		Q(j-1,j-1) = c;
		Q(i-1,j-1) = -s;
    Q(j-1,i-1) = s;  
    W = Q * W;
		k++;		
	}//end inner for		
}//end outer for 

arma::mat tW = trans(W);
arma::mat S = Z * tW;

arma::mat SH = abs(S);
for(int j=0; j<d; j++){
  for(int i=0; i<n; i++){
    double sTemp = SH(i,j);
    if(sTemp <= C){
      SH(i,j) = sTemp*sTemp;
    }
    if(sTemp > C){
      SH(i,j) = 2*C*sTemp - C*C;
    }
  }
}

arma::mat ColMeans = mean(SH, 0);

int z=0; //index for out
double ans= 0.0; 
int b = (d*(d-1))/2 + L*d*(d-1); //number of entries in output
arma::colvec out(b);

//Lag 0. A total of d(d-1)/2 entries in out will be filled. i and j are column numbers. k is row number 
for(int i=0;i<d;i++){	
	for(int j=i+1;j<d;j++){
		
		ans = 0.0; //next column so restart 
		ans = dot(SH.col(i), SH.col(j));							
		out[z] = ans/n - ColMeans(i)*ColMeans(j); //adjust by cross mean and store in out 			
		z++;							
	}//end j for
}//end i for  -- END Lag 0

//z now contains the index of the next position for outptr
/**Lag 1 and above. */
  
	for(int lag=1;lag<=L;lag++){		
		for(int i=0;i<d;i++){	
			for(int j=0;j<d;j++){			
				ans=0.0;
				
				if(i!=j){      
        // A.submat(p, r, q, s)	A.submat(first_row, first_col, last_row, last_col)       
        ans = dot(SH.submat(0,i,n-lag-1,i), SH.submat(lag,j,n-1,j));							
        out[z] = ans/n - ColMeans(i)*ColMeans(j); //adjust by cross mean and store in out 			
        z++;
					}else{ 
						//Do nothing				
					} //end if-else 			
			} //end j for	-- Finished column in 2nd matrix. Move to next one. 
		} //end i for -- Finished with column in 1st matrix. Move to next one. 
} //end lag for. Finished calculations for this lag. Next lag. 
	
double obj = dot((out % PHI), out);

return Rcpp::wrap(obj);
'

DOCinVar = cxxfunction(signature(theta_s="numeric", Z_s="matrix", L_s="integer", PHI_s="numeric", C_s="numeric"), 
                       src, Rcpp=TRUE, verbose=TRUE, plugin = "RcppArmadillo",
                       libargs = c("-L/opt/local/lib -lgsl -lgslcblas -lm"), #-lgsl 
                       cppargs = c("-I/opt/local/include"), 
                       includes = c("#include <stdio.h>")
)

###### Doc test function

DOC.test = function(A, m){
  N = dim(A)[1]
  k = dim(A)[2]
  temp = numeric(m)
  pval = numeric(m)
  out = as.data.frame(matrix(0,m+1,4))
  names(out) = c("m", "Q(m)", "d.f.", "p-value")
  out[,1] = 0:m
  
  Q.temp = N*sum(cor(A)[lower.tri(cor(A), diag = FALSE)]^2)
  out[1,2] = Q.temp 
  df.temp = k*(k-1)/2
  out[1,3] = df.temp
  out[1,4] = 1-pchisq(Q.temp, df.temp)
  
  for(j in 1:m){
    ccf = cor(A[-(1:j),],A[-((N-j+1):N),])
    Q.temp = Q.temp + N*(N+2)*sum(ccf[lower.tri(ccf, diag = FALSE)]^2)/(N-j)
    Q.temp = Q.temp + N*(N+2)*sum(ccf[upper.tri(ccf, diag = FALSE)]^2)/(N-j)
    out[(j+1),2] = Q.temp 
    df.temp = df.temp+ k*(k-1)
    out[(j+1),3] = df.temp
    out[(j+1),4] = 1-pchisq(Q.temp, df.temp)
  }	
  round(out,3)
}

###### Import data files

mydir="month"
myfiles=list.files(path = mydir,pattern = "*.csv",full.names=TRUE)

###### Estimation Doc in mean, Doc in volatility, garch(1,1) series and DOC tes

V1=list()
D1=list()
for (i in 1:63) {
  a=read_csv(myfiles[i])
  f=a[,c(1,2)]
  x=a[,-c(1,2)]
  x=as.matrix(x)
  x=scale(x)
  # PCA "by hand"
  U.hat = matrix.sqrt.inv.pc(cov(x))
  Z.hat = x %*% t(U.hat)
  
  k = 10
  Z.hat =  Z.hat[,1:k]
  size=dim(Z.hat)
  L=1
  P=size[1]
  d=size[2]
  # Weight matrix
  p = d*(d-1)/2 ; p2 = 2*p
  phi = 1 - (0:L)/(L+1) ; phi.total = sum(phi) ; phi = phi / phi.total
  PHI = c( rep(phi[1],p), rep(phi[-1],each = p2) ) 
  # DOC in mean estimation
  
  out = nlminb(start = rep(0.66,p), objective = DOCinMean,  # Change this to DOCinMean as necessary
               gradient = NULL, hessian = NULL, scale = 1, control = list(iter.max = 1000, eval.max = 1000), 
               lower = -pi, upper = pi, Z_s = Z.hat, L_s = 1, PHI_s = PHI)
  theta.hat = out$par  
  W.hat = (theta2W(out$par))  # Separating Matrix
  S.hat = Z.hat %*% t(W.hat)
  # Residual vector error
  e=Z.hat-S.hat
  #Uncorrelated residual errors
  U = matrix.sqrt.inv.pc(cov(e))
  Z = e %*% t(U)
  #Doc in volatility estimation
  out = nlminb(start =rep(0.66,p), objective = DOCinVar,  # Change this to DOCinMean as necessary
               gradient = NULL, hessian = NULL, scale = 1, control = list(iter.max = 1000, eval.max = 1000), 
               lower = -pi, upper = pi, Z_s = Z, L_s = 1, PHI_s = PHI, C_s= 2.25)
  theta = out$par  
  W = (theta2W(out$par))  # Separating Matrix
  S = Z %*% t(W)
  A=S^2
  m=1
  B= matrix(NA,ncol =d, nrow = P)
  for (j in 1:d) {
    h=data.frame(S)
    ug_spec= ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder=c(0, 0)))
    ugfit = ugarchfit(spec = ug_spec, data = h[,j],solver = "hybrid")
    sigma1 = rugarch::sigma(ugfit)
    B[,j]= sigma1
    
    
  }
  E=data.frame(f,B)
  V1[[i]]=E
  D1[[i]]=DOC.test(A,m)
}
D1
###### Export data

a1=data.frame(V1[1])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file1.csv")

a2=data.frame(V1[2])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file2.csv")


a3=data.frame(V1[3])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file3.csv")

a4=data.frame(V1[4])
write.csv(a4,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file4.csv")

a5=data.frame(V1[5])
write.csv(a5,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file5.csv")

a6=data.frame(V1[6])
write.csv(a6,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file6.csv")


a7=data.frame(V1[7])
write.csv(a7,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file7.csv")

a8=data.frame(V1[8])
write.csv(a8,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file8.csv")

a9=data.frame(V1[9])
write.csv(a9,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file9.csv")

a10=data.frame(V1[10])
write.csv(a10,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file10.csv")


a1=data.frame(V1[11])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file11.csv")

a2=data.frame(V1[12])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file12.csv")


a3=data.frame(V1[13])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file13.csv")

a4=data.frame(V1[14])
write.csv(a4,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file14.csv")

a5=data.frame(V1[15])
write.csv(a5,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file15.csv")

a6=data.frame(V1[16])
write.csv(a6,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file16.csv")


a7=data.frame(V1[17])
write.csv(a7,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file17.csv")

a8=data.frame(V1[18])
write.csv(a8,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file18.csv")

a9=data.frame(V1[19])
write.csv(a9,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file19.csv")

a10=data.frame(V1[20])
write.csv(a10,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file20.csv")


a1=data.frame(V1[21])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file21.csv")

a2=data.frame(V1[22])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file22.csv")


a3=data.frame(V1[23])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file23.csv")

a4=data.frame(V1[24])
write.csv(a4,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file24.csv")

a5=data.frame(V1[25])
write.csv(a5,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file25.csv")

a6=data.frame(V1[26])
write.csv(a6,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file26.csv")


a7=data.frame(V1[27])
write.csv(a7,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file27.csv")

a8=data.frame(V1[28])
write.csv(a8,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file28.csv")

a9=data.frame(V1[29])
write.csv(a9,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file29.csv")

a10=data.frame(V1[30])
write.csv(a10,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file30.csv")


a1=data.frame(V1[31])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file31.csv")

a2=data.frame(V1[32])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file32.csv")


a3=data.frame(V1[33])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file33.csv")

a4=data.frame(V1[34])
write.csv(a4,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file34.csv")

a5=data.frame(V1[35])
write.csv(a5,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file35.csv")

a6=data.frame(V1[36])
write.csv(a6,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file36.csv")


a7=data.frame(V1[37])
write.csv(a7,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file37.csv")

a8=data.frame(V1[38])
write.csv(a8,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file38.csv")

a9=data.frame(V1[39])
write.csv(a9,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file39.csv")

a10=data.frame(V1[40])
write.csv(a10,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file40.csv")



a1=data.frame(V1[41])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file41.csv")

a2=data.frame(V1[42])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file42.csv")


a3=data.frame(V1[43])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file43.csv")

a4=data.frame(V1[44])
write.csv(a4,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file44.csv")

a5=data.frame(V1[45])
write.csv(a5,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file45.csv")

a6=data.frame(V1[46])
write.csv(a6,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file46.csv")


a7=data.frame(V1[47])
write.csv(a7,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file47.csv")

a8=data.frame(V1[48])
write.csv(a8,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file48.csv")

a9=data.frame(V1[49])
write.csv(a9,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file49.csv")

a10=data.frame(V1[50])
write.csv(a10,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file50.csv")


a1=data.frame(V1[51])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file51.csv")

a2=data.frame(V1[52])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file52.csv")


a3=data.frame(V1[53])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file53.csv")

a4=data.frame(V1[54])
write.csv(a4,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file54.csv")

a5=data.frame(V1[55])
write.csv(a5,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file55.csv")

a6=data.frame(V1[56])
write.csv(a6,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file56.csv")


a7=data.frame(V1[57])
write.csv(a7,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file57.csv")

a8=data.frame(V1[58])
write.csv(a8,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file58.csv")

a9=data.frame(V1[59])
write.csv(a9,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file59.csv")

a10=data.frame(V1[60])
write.csv(a10,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file60.csv")


a1=data.frame(V1[61])
write.csv(a1,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file61.csv")

a2=data.frame(V1[62])
write.csv(a2,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file62.csv")


a3=data.frame(V1[63])
write.csv(a3,"C:/Users/Emile Ndoumbe/Dropbox/emile1/garch_month/file63.csv")



