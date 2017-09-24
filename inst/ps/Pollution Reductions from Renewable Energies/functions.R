# This function computes MEA(W,D)

mea = function(beta,W,D) {
  mea = 0
  
  for (i in 1:3) {
    for (j in 0:3) {
      
      mea = mea - i*beta[i,j+1]*W^(i-1)*D^j
      
    }
  }
  return(mea)
}


# Plotting MEA(W,D)

plot_mea = function(pol,data) {
  frml = formula(paste("txok_",pol,"_tr ~ 
                          wind_tr + wind2_tr + wind3_tr +
                       wind_load_tr + wind2_load_tr + wind3_load_tr + 
                       wind_load2_tr + wind2_load2_tr + wind3_load2_tr +
                       wind_load3_tr + wind2_load3_tr + wind3_load3_tr +
                       load_tr + load2_tr + load3_tr +
                       spp_load_tr + spp_load2_tr + spp_load3_tr", sep=""))
  mgem = lm(frml ,data)
  
  beta = matrix(coef(mgem)[2:13], nrow=3, ncol=4)
  
  if (pol=="co2") {
    title = "Marginal CO2 Emissions Avoided"
  } else if (pol=="nox") {
    title = "Marginal NOx Emissions Avoided"
  } else {
    title = "Marginal SO2 Emissions Avoided"
  }
  
  library(ggplot2)
  ggplot(data.frame(x = c(23500, 55000)), aes(x = x)) +
    stat_function(fun = mea, args = list(beta = beta, W = 125), aes(colour=" 125")) +
    stat_function(fun = mea, args = list(beta = beta, W = 820), aes(colour=" 820")) +
    stat_function(fun = mea, args = list(beta = beta, W = 1935), aes(colour="1935")) +
    stat_function(fun = mea, args = list(beta = beta, W = 3285), aes(colour="3285")) +
    xlab("Demand [MWh]") +
    ylab("MEA [tons / MWh]") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual("Wind [MWh]", values = c("red", "blue", "green", "orange"))
}


# This function computes AEA(K)
aea = function(k, xdata, Ddata, beta, xbin=40, Dbin=40) {
  
  library(GenKern)
  dist_est = KernSur(xdata, Ddata, xbin, Dbin)
  
  x = dist_est$xords
  D = dist_est$yords
  f = dist_est$zden

  df = data.frame(K=k, AEA=rep(NA, times=length(k)))
  
  aea_num = 0
  aea_den = 0
  
  for (i in 1:xbin) {
    for (j in 1:Dbin) {
      
      mea = 0
      
      for (m in 1:3) {
        for (n in 0:3) {
          mea = mea - m*beta[m,n+1]*(x[i]*k)^(m-1)*D[j]^n
          
        }
      }
      
      aea_num = aea_num + x[i]*mea*f[i,j]
      aea_den = aea_den + x[i]*f[i,j]
      df[which(df[,1]==k),2] = aea_num/aea_den
      
    }
  }
  return(df)
}


# This function computes the beta_i,j, which are required in MEA(x*K,D).

beta_ij = function(data) {
  
  mgem_co2 = lm(txok_co2_tr ~
                  wind_tr + wind2_tr + wind3_tr +
                  wind_load_tr + wind2_load_tr + wind3_load_tr + 
                  wind_load2_tr + wind2_load2_tr + wind3_load2_tr +
                  wind_load3_tr + wind2_load3_tr + wind3_load3_tr +
                  load_tr + load2_tr + load3_tr +
                  spp_load_tr + spp_load2_tr + spp_load3_tr, data)
  beta_co2 = matrix(coef(mgem_co2)[2:13], nrow=3, ncol=4)
  
  mgem_nox = lm(txok_nox_tr ~
                  wind_tr + wind2_tr + wind3_tr +
                  wind_load_tr + wind2_load_tr + wind3_load_tr + 
                  wind_load2_tr + wind2_load2_tr + wind3_load2_tr +
                  wind_load3_tr + wind2_load3_tr + wind3_load3_tr +
                  load_tr + load2_tr + load3_tr +
                  spp_load_tr + spp_load2_tr + spp_load3_tr, data)
  beta_nox = matrix(coef(mgem_nox)[2:13], nrow=3, ncol=4)
  
  mgem_so2 = lm(txok_so2_tr ~
                  wind_tr + wind2_tr + wind3_tr +
                  wind_load_tr + wind2_load_tr + wind3_load_tr + 
                  wind_load2_tr + wind2_load2_tr + wind3_load2_tr +
                  wind_load3_tr + wind2_load3_tr + wind3_load3_tr +
                  load_tr + load2_tr + load3_tr +
                  spp_load_tr + spp_load2_tr + spp_load3_tr, data)
  beta_so2 = matrix(coef(mgem_so2)[2:13], nrow=3, ncol=4)
  
  beta = list("co2"=beta_co2, "nox"=beta_nox, "so2"=beta_so2)
  return(beta)
}


# standard errors of AEA(k) for wind

se_aea = function(K, xdata, Ddata, regdata, xbin=40, ybin=40) {
  
  library(GenKern)
  dist_est = KernSur(xdata, Ddata, xbin, ybin)
  
  x = matrix(rep(dist_est$xords, each = 40))
  D = matrix(rep(dist_est$yords, times = 40))
  f = matrix(t(dist_est$zden), ncol=1)
  
  for (k in K) {
    
    weights = matrix(nrow=1600, ncol=12)
    
    for (i in 1:1600) {
      for (m in 0:3) {
        for (n in 1:3) {
          weights[i,m*3+n]=n*(k*x[i])^(n-1)*D[i]^m
        }
      }
    }
    
    name = paste("weightsmea_", k, sep="")
    assign(name, weights)
  }
  
  
  
  for (i in c("co2", "nox", "so2")) {
    
    frml = formula(paste("txok_",i,"_tr ~
                  wind_tr + wind2_tr + wind3_tr +
                  wind_load_tr + wind2_load_tr + wind3_load_tr + 
                  wind_load2_tr + wind2_load2_tr + wind3_load2_tr +
                  wind_load3_tr + wind2_load3_tr + wind3_load3_tr +
                  load_tr + load2_tr + load3_tr +
                  spp_load_tr + spp_load2_tr + spp_load3_tr", sep=""))
    mgem = lm(frml, regdata)
    
    library(sandwich)
    vcov_mgem = NeweyWest(mgem, lag=24, prewhite=FALSE)[2:13,2:13]*(43794-18-1)/(43794-18-1825-9660)

    name = paste("vcov_mgem_", i, sep="")
    assign(name, vcov_mgem)
  }
  
  
  df = data.frame(K=K, co2=rep(NA, times=length(K)), nox=rep(NA, times=length(K)), so2=rep(NA, times=length(K)))  
  
  for (i in c("co2", "nox", "so2")) {
    for (k in K) {
      
      weights = get(paste("weightsmea_", k, sep=""))
      vcov_mgem = get(paste("vcov_mgem_", i, sep=""))
      
      se = sqrt(t(x*f) %*% weights %*% vcov_mgem %*% t(weights) %*% (x*f) / (t(x) %*% f)^2)
      df[which(df[,1] == k),which(colnames(df)==i)] = se[1,1]
      
    }
  }
  
  return(df)
}




