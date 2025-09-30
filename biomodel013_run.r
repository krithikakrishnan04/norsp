# sample code for BioModel "013" 

library(igraph)
library(MASS)
require(RColorBrewer)
library(leaps)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(RCy3)

total = 0 
speciesOut <- c()
speciesIn <- c()
species.df.collect <- c()


  
  model_name = "013"
  jmat = read.csv(file = "jac013.csv")
  fmat = read.csv(file="fTable013.csv")
  
  # get perturbed species details 
  
  path <- paste("maxSpeciesPerturbed",model_name,".txt",sep="")
  
  speciesMax = read.table(path,sep=",",header = T)
  perturbed_name = speciesMax$Name
  
  # process jacobian 
  jmat = jmat[,-1] 
  jnames = colnames(jmat)
  rownames(jmat) = jnames
  
  J0 = jmat  # remove nodes with no out-degree
  mus <- apply(abs(J0),2,sum)
  inds_rm <- which(mus==0)
  if (length(inds_rm)>0){ 
    J <- J0[-inds_rm,-inds_rm]
  } else {
    J <- J0
  }
  
  A <- abs(sign(J))
  G = graph_from_adjacency_matrix(abs(t(sign(J))),weighted=T)
  
  igraph::V(G)$comp <- igraph::components(G)$membership
  table(V(G)$comp)
  main <- induced_subgraph(G,V(G)$comp==1) # removing unconnected nodes
  # remove self loops
  G = simplify(main)
  # process f table for inconsistent names
  fmat[,1] <-  colnames(jmat)
  
  # plot graph
  setEPS()
  postscript(paste('infuenceMat_',model_name,".eps"),
             width = 10.0, height = 10.0, fonts=c("serif", "Palatino"))
  V(G)$color <- ifelse(V(G)$name %in% speciesMax$Name , "red", "lightblue")
  plot(G, vertex.size=15, vertex.label.family="serif",
       edge.label.family="Palatino",vertex.label.cex=0.9,edge.arrow.size= 0.6)
  legend('topright',legend=c("Perturbed","Unperturbed"),fill=c("red","lightblue"))
  dev.off()
  
  # adding reduced number of vertices to a new species list 
  speciesNames <- V(G)$name
  
  if((speciesMax$Name %in% speciesNames) == T){
    speciesMax$Index = which(speciesNames %in% speciesMax$Name)
    perturbed = speciesMax$Index
    speciesIn = c(speciesIn,speciesMax$Name)
    total = total + 1 # count total models with perturbed nodes present in subset selection 
  }else{
    perturbed = 0
    speciesOut = c(speciesOut,speciesMax$Name)
  }
  
  inds_conv <- match(speciesNames, V(G)$name)
  J <- J[speciesNames,colnames(J) %in% speciesNames]
  Wp <- t(J)
  I <- diag(length(V(G)))
  
  # propagation models
  alpha = .9 # fixed 
  

  
  ##############################################################################
  # Undirected Propagation model ###############################################
  ##############################################################################
  
  W <- t(abs(t(sign(J))) + abs(sign(J)))
  W[W>0] <- 1
  
  D1 <- matrix(0,length(V(G)),length(V(G)))
  D2 <- matrix(0,length(V(G)),length(V(G)))
  # out degree prob
  diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
  # in degree prob
  diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
  # Wp[i,j] is Wp divided by out degree of i and in degree of j
  Wp <- D1 %*% W %*% D2
  
  I <- diag(length(V(G)))
  Mprince_undir <- solve(I-alpha * Wp) * (1-alpha)
  
  colnames(Mprince_undir) = speciesNames
  rownames(Mprince_undir) = speciesNames
  ##############################################################################
  # Directed Propagation model 
  ##############################################################################
  
  W <- abs(t(sign(J)))
  D1 <- matrix(0,length(V(G)),length(V(G)))
  D2 <- matrix(0,length(V(G)),length(V(G)))
  # out degree prob
  diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
  # in degree prob
  diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
  # Wp[i,j] is Wp divided by out degree of i and in degree of j
  Wp <- D1 %*% W %*% D2
  
  I <- diag(length(V(G)))
  Mprince_dir <- solve(I-alpha * Wp) * (1-alpha)
  
  colnames(Mprince_dir) = speciesNames
  rownames(Mprince_dir) = speciesNames
  

  ##############################################################################
  ## Directed signed Propagation model 
  ##############################################################################
  
  W <- t(sign(J))
  D1 <- matrix(0,length(V(G)),length(V(G)))
  D2 <- matrix(0,length(V(G)),length(V(G)))
  # out degree prob
  diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
  # in degree prob
  diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
  # Wp[i,j] is Wp divided by out degree of i and in degree of j
  Wp <- D1 %*% W %*% D2
  
  I <- diag(length(V(G)))
  Mprince_dir_sign <- solve(I-alpha * Wp) * (1-alpha)
  
  colnames(Mprince_dir_sign) = speciesNames
  rownames(Mprince_dir_sign) = speciesNames
  
  ##############################################################################
  ## subset selection ##########################################################
  ##############################################################################
  
  # retrieve F data
  colnames(fmat)[1] = "Species"
  fmat.sub = fmat[fmat$Species %in% speciesNames,]
  
  x_control <- fmat.sub[,2]
  x_treat <- fmat.sub[,3]
  ##############################################################################
  # prediction function for subset selection
  predict.regsubsets = function(object, newdata, id){
    form = as.formula(object$call[[2]])
    mat = model.matrix(form, newdata)
    coefi = coef(object,id=id)
    xvars=names(coefi)
    if (id > 1){
      pred = mat[,xvars]%*%coefi
    }
    if (id == 1){
      pred = t(mat[,xvars] * coefi)
    }
    return(pred)
  }
  ##############################################################################
  n = length(speciesNames)
  
  s_true_list <- list(Mprince_undir,Mprince_dir,Mprince_dir_sign) # exact sensitivity matrix
  names(s_true_list) <- c("undir","dir","dir_signed")
  method <- names(s_true_list)
  total.models <- c()
  tally.tab.methodwise <- c()
  F_hat_sel_compile <- c()
  delta_hat_compile <- c()
  ##############################################################################
  F_0 = x_control
  F_t = x_treat
  y = F_t - F_0
  ##############################################################################
  # initial error
  F_init_err <- sqrt(sum(F_t-F_0)^2)/sqrt(sum(F_t)^2) # relative initial error
  ##############################################################################
  thresh = .9 # threshold for model selection on adj R2 
  
  F_err_store <- c()
  F_err_raw_store <- c()
  tot_sel_store <-c()
  hits_store <-c()
  miss_rate_store <-c()
  ## removed loops which run through all 3 propagation methods
  smat = s_true_list["undir"] # change s_true_list["dir"] or s_true_list["dir_signed"] 
    
  
    s_true <- as.data.frame(smat)
    models.fr <- regsubsets(y~., data=s_true, method="forward", intercept = FALSE, nvmax = n)
    res.sum.fr <- summary(models.fr)
    
    res.sum.fr$adjr2 <- unlist(lapply(res.sum.fr$adjr2,function(x) replace(x, is.na(x), 0)))
  
    best_adj_r2_val = max(res.sum.fr$adjr2)
    
    if (best_adj_r2_val < thresh){ # general case
      best_adj_r2 = length(res.sum.fr$adjr2) # we take the largest model if we do not reach the threshold
    }
    
    if (best_adj_r2_val > thresh){ # really good -- we take the first one that exceeds the threshold
      temp <- which(res.sum.fr$adjr2 > thresh) # find all above the threshold
      best_adj_r2 = temp[1] # take the first one
    }
    
  
    for (kk in 1:n){
      pred_val <- predict(models.fr, newdata = s_true, id = kk) # predicted values
      if (kk == 1){ # if the model only has one variable, we need to take the transpose
        pred_val <- t(pred_val)
      }
      
      F_hat_sel = F_0 + pred_val
      F_err <- sqrt(sum(F_t-F_hat_sel)^2)/sqrt(sum(F_t)^2) # relative error
      F_err_raw <- sqrt(sum(F_t-F_hat_sel)^2) # raw error 
  
  
      F_err_store <- c(F_err_store, F_err)
      F_err_raw_store <- c(F_err_raw_store, F_err_raw) 
      
    }
    
    # a binary vector to track information
    pert <- rep(0,n)
    pert[perturbed] <- 1
    
    picks <- which(res.sum.fr$outmat[best_adj_r2,] == "*") # which ones were selected
    selected <- rep(0,n)
    selected[picks] <- 1
    # Does delta_hat = delta_true via subset selection ?
    delta_hat <- rep(0,n)
    delta_hat[picks] <- coef(models.fr, best_adj_r2) # the selected cofficients
    
    # track total number selected
    tot <- length(picks)
    hits <- length(intersect(perturbed, picks))
    
    ## represent the error of those selected as a misclassification rate
    tot_sel_store <- c(tot_sel_store, tot)
    hits_store <- c(hits_store, hits)
    miss_rate <- (1/n)*sum(abs(pert-selected))
   
  
  
  plot(models.fr,main="Forward Selection")

  plot(res.sum.fr$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "b")
  points(best_adj_r2, res.sum.fr$adjr2[best_adj_r2],col = "red")
  
  data.mat.without.F.error.raw <- data.frame(x_label = 1:n,
                                             "F error" = F_err_store)
  
  data_long <- melt(data.mat.without.F.error.raw, id = "x_label")
  
  F.error <- ggplot(data_long,
                    aes(x = x_label,
                        y = value,
                        color = variable)) +  geom_line(color="blue") +
    theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
    ylab('F Error')+
    xlab('Number of nodes (n)') + theme_grey(base_size = 22)
  
  F.error <- F.error + theme(legend.position="none")
  F.error + geom_hline(yintercept=F_init_err, linetype="dashed", color = "red")
  
