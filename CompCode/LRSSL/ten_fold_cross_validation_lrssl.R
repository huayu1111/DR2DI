get.knn <- function(v, k){
  ind <- order(v, decreasing = T)
  return(ind[1:k])
}

get.knn.graph <- function(S, k){
  n <- nrow(S)
  S.knn <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    ind <- get.knn(S[i,], k)
    S.knn[i,ind] <- 1
    S.knn[ind,i] <- 1
    #S.kkn[i,ind] <- S[i,ind]
    #S.knn[ind,i] <- S[i,ind]
  }
  return(S.knn)
}

lrssl <- function(Xs, Ls, Y, mx, ml, mu, lam, gam, max.iter, eps){
  n <- nrow(Y)
  c <- ncol(Y)
  Gs <- vector("list", mx)
  e1ds <- vector("list", mx)
  ds <- rep(0, mx)
  for(i in 1:mx){
    ds[i] <- nrow(Xs[[i]])
    Gs[[i]] <- matrix(runif(ds[i]*c), nrow = ds[i], ncol = c)
    e1ds[[i]] <- rep(1, ds[i])
  }
  alpha <- rep(1, ml)/ml
  
  check.step <- 20;
  Gs.old <- vector("list", mx)
  As <- vector("list", mx)
  As.pos <- vector("list", mx)
  As.neg <- vector("list", mx)
  Bs <- vector("list", mx)
  Bs.pos <- vector("list", mx)
  Bs.neg <- vector("list", mx)
  L <- matrix(0, nrow = n, ncol = n)
  for(i in 1:ml){
    L <- L + alpha[i]^gam*Ls[[i]]
  }
  
  t <- 0
  while(t < max.iter){
    t <- t + 1
    Q <- Y
    for(i in 1:mx){
      Gs.old[[i]] <- Gs[[i]]
      Q <- Q + mu*t(Xs[[i]])%*%Gs[[i]]
    }
    P <- solve(L + (1 + mx*mu)*diag(1, n, n))
    F.mat <- P%*%Q
    
    for(i in 1:mx){
      As[[i]] <- Xs[[i]]%*%(mu*diag(1,n,n)-mu^2*t(P))%*%t(Xs[[i]]) + lam*(e1ds[[i]]%*%t(e1ds[[i]]))
      As.pos[[i]] <- (As[[i]] + abs(As[[i]]))/2
      As.neg[[i]] <- (abs(As[[i]]) - As[[i]])/2
      Bs[[i]] <- mu*Xs[[i]]%*%P%*%Y
      for(j in 1:mx){
        if(i == j){
          next
        }else{
          Bs[[i]] <- Bs[[i]] + mu^2*Xs[[i]]%*%t(P)%*%t(Xs[[j]])%*%Gs[[j]]
        }
      }
      Bs.pos[[i]] <- (Bs[[i]] + abs(Bs[[i]]))/2
      Bs.neg[[i]] <- (abs(Bs[[i]]) - Bs[[i]])/2
    }
    for(i in 1:mx){
      Gs[[i]] <- Gs[[i]]*sqrt((Bs.pos[[i]] + As.neg[[i]]%*%Gs[[i]])/(Bs.neg[[i]] + As.pos[[i]]%*%Gs[[i]]))
    }
    
    for(i in 1:ml){
      alpha[i] <- (1/sum(diag(t(F.mat)%*%Ls[[i]]%*%F.mat)))^(1/(gam - 1))
    }
    alpha <- alpha/sum(alpha)
    
    L <- matrix(0, nrow = n, ncol = n)
    for(i in 1:ml){
      L <- L + alpha[i]^gam*Ls[[i]]
    }
    
    diff.G <- rep(0, mx)
    for(i in 1:mx){
      diff.G[i] <- norm(Gs[[i]] - Gs.old[[i]], "f")/norm(Gs.old[[i]], "f")
    }
    
    if(t%%check.step == 0){
      mesg <- sprintf("t = %i, diffG mean = %.4e", t, mean(diff.G))
      print(mesg)
    }
    if(mean(diff.G) < eps)
      break
  }
  return(list(Gs = Gs, F.mat = F.mat, alpha = alpha, diff.G = diff.G, t = t))
}
##############################read training data###################################################
chem.mat.all <-read.table("./drug_pubchem_info.txt", sep = " ", header = T, row.names = 1)
dom.mat.all <- read.table("./drug_target_domain_info.txt", sep = " ", header = T, row.names = 1)
go.mat.all <- read.table("./drug_target_go_info.txt", sep = " ", header = T, row.names = 1)
phar.mat <- read.table("./drug_dis_mat.txt", sep = "\t", header = T, row.names = 1)
dis.sim <- read.table("./disease_sim_matrix.txt", sep = "\t", header = T, row.names = 1)
matchIndexes <- match(colnames(phar.mat),rownames(dis.sim))
dis.sim <- dis.sim[matchIndexes,matchIndexes]
drug.num <- floor(nrow(phar.mat)/10)
chem.mat <- 0
dom.mat <- 0
go.mat <- 0
chem.mat.new <- 0
dom.mat.new <- 0
go.mat.new <- 0
Y <- 0
res <- matrix(0,nrow=1,ncol=2)

for(iterer in seq(1,10)){
	if(iterer == 10){
		chem.mat <- chem.mat.all[seq(1,(9*drug.num)),]
		dom.mat <- dom.mat.all[seq(1,(9*drug.num)),]
		go.mat <- go.mat.all[seq(1,(9*drug.num)),]
		chem.mat.new <- chem.mat.all[seq(((iterer-1)*drug.num+1),nrow(chem.mat.all)),]
		dom.mat.new <- dom.mat.all[seq(((iterer-1)*drug.num+1),nrow(chem.mat.all)),]
		go.mat.new <- go.mat.all[seq(((iterer-1)*drug.num+1),nrow(chem.mat.all)),]
		Y <- as.matrix(phar.mat[seq(1,(9*drug.num)),])
	}else if(iterer==1){
		chem.mat <- chem.mat.all[seq((iterer*drug.num+1),nrow(chem.mat.all)),]
		dom.mat <- dom.mat.all[seq((iterer*drug.num+1),nrow(dom.mat.all)),]
		go.mat <- go.mat.all[seq((iterer*drug.num + 1),nrow(go.mat.all)),]
		chem.mat.new <- chem.mat.all[seq(((iterer-1)*drug.num+1),(iterer*drug.num)),]
		dom.mat.new <- dom.mat.all[seq(((iterer-1)*drug.num+1),(iterer*drug.num)),]
		go.mat.new <- go.mat.all[seq(((iterer-1)*drug.num+1),(iterer*drug.num)),]
		Y <- as.matrix(phar.mat[seq((iterer*drug.num+1),nrow(chem.mat.all)),])
	}else{
		chem.mat <- chem.mat.all[c(seq(((iterer-2)*drug.num+1),((iterer-1)*drug.num)),seq((iterer*drug.num + 1),nrow(chem.mat.all))),]
		dom.mat <- dom.mat.all[c(seq(((iterer-2)*drug.num+1),((iterer-1)*drug.num)),seq((iterer*drug.num + 1),nrow(dom.mat.all))),]
		go.mat <- go.mat.all[c(seq(((iterer-2)*drug.num+1),((iterer-1)*drug.num)),seq((iterer*drug.num + 1),nrow(go.mat.all))),]
		chem.mat.new <- chem.mat.all[seq(((iterer-1)*drug.num+1),(iterer*drug.num)),]
		dom.mat.new <- dom.mat.all[seq(((iterer-1)*drug.num+1),(iterer*drug.num)),]
		go.mat.new <- go.mat.all[seq(((iterer-1)*drug.num+1),(iterer*drug.num)),]
		Y <- as.matrix(phar.mat[c(seq(((iterer-2)*drug.num+1),((iterer-1)*drug.num)),seq((iterer*drug.num + 1),nrow(chem.mat.all))),])
	}
	
	X1 <- as.matrix(t(chem.mat))
	X2 <- as.matrix(t(dom.mat))
	X3 <- as.matrix(t(go.mat))

	n <- nrow(Y)
	c <- ncol(Y)

	Xs <- list(X1, X2, X3)
	##############################caculate knn graph###################################################
	k <- 10

	S1 <- t(X1)%*%X1/sqrt(colSums(X1)%*%t(colSums(X1)))
	S2 <- t(X2)%*%X2/sqrt(colSums(X1)%*%t(colSums(X1)))
	S3 <- t(X3)%*%X3/sqrt(colSums(X3)%*%t(colSums(X3)))

	S1 <- S1 - diag(diag(S1))
	S2 <- S2 - diag(diag(S2))
	S3 <- S3 - diag(diag(S3))

	S4 <- matrix(0, nrow = n, ncol = n)
	for(i in 1:(n-1)){
	  for(j in (i+1):n){
		s <- dis.sim[Y[i,]==1,Y[j,]==1]
		S4[i,j] <- max(s)
		#S4[i,j] <- mean(s)
	  }
	}
	S4 <- S4 + t(S4)

	S1.knn <- get.knn.graph(S1, k)
	S2.knn <- get.knn.graph(S2, k)
	S3.knn <- get.knn.graph(S3, k)
	S4.knn <- get.knn.graph(S4, k)

	D1 <- diag(colSums(S1.knn))
	D2 <- diag(colSums(S2.knn))
	D3 <- diag(colSums(S3.knn))
	D4 <- diag(colSums(S4.knn))

	L1 <- D1 - S1.knn
	L2 <- D2 - S2.knn
	L3 <- D3 - S3.knn
	L4 <- D4 - S4.knn

	Ls <- list(L1, L2, L3, L4)
	##############################Training##############################################################
	mx <- 3
	ml <- 4
	mu <- 0.01
	lam <- 0.01
	gam <- 2

	train.res <- lrssl(Xs, Ls, Y, mx, ml, mu, lam, gam, 100, 1e-6)
	###############################Indication prediction for new drugs#################################


	X1.new <- as.matrix(chem.mat.new)
	X2.new <- as.matrix(dom.mat.new)
	X3.new <- as.matrix(go.mat.new)

	G1 <- train.res$Gs[[1]]
	G2 <- train.res$Gs[[2]]
	G3 <- train.res$Gs[[3]]

	Y.pred1 <- X1.new%*%G1
	Y.pred2 <- X2.new%*%G2
	Y.pred3 <- X3.new%*%G3

	Y.pred1 <- Y.pred1/rowSums(Y.pred1)
	Y.pred2 <- Y.pred2/rowSums(Y.pred2)
	Y.pred3 <- Y.pred3/rowSums(Y.pred3)

	alpha <- train.res$alpha
	Y.pred <- alpha[1]*Y.pred1 + alpha[2]*Y.pred2 + alpha[3]*Y.pred3
	if(iterer==10){
		Y = phar.mat[((iterer-1)*drug.num+1):nrow(chem.mat.all),];
		pos.indexes <- which(Y==1,arr.ind=TRUE);
		neg.indexes <- which(Y==0,arr.ind=TRUE);
		random.indexes <- sample(nrow(neg.indexes),nrow(pos.indexes))
		neg.indexes <- neg.indexes[random.indexes,]
		cross_res <- cbind(c(Y.pred[pos.indexes],Y.pred[neg.indexes]),c(rep(1,nrow(pos.indexes)),rep(0,nrow(neg.indexes))))
		if(nrow(res) == 1){
			res = cross_res;
		}else{
		  cat(dim(res),"\n")
		  cat(dim(cross_res),"\n")
			res = rbind(res,cross_res);
		}
	}else if(iterer==1){
		Y = phar.mat[((iterer-1)*drug.num+1):(iterer*drug.num),];
		pos.indexes <- which(Y==1,arr.ind=TRUE);
		neg.indexes <- which(Y==0,arr.ind=TRUE);
		random.indexes <- sample(nrow(neg.indexes),nrow(pos.indexes))
		neg.indexes <- neg.indexes[random.indexes,]
		cross_res = cbind(c(Y.pred[pos.indexes],Y.pred[neg.indexes]),c(rep(1,nrow(pos.indexes)),rep(0,nrow(neg.indexes))))
		if(nrow(res)==1){
			res = cross_res;
		}else{
		  cat(dim(res),"\n")
		  cat(dim(cross_res),"\n")
			res = rbind(res,cross_res);
		}
	}else{
		Y = phar.mat[((iterer-1)*drug.num+1):(iterer*drug.num),];
		pos.indexes <- which(Y==1,arr.ind=TRUE);
		neg.indexes <- which(Y==0,arr.ind=TRUE);
		random.indexes <- sample(nrow(neg.indexes),nrow(pos.indexes))
		neg.indexes <- neg.indexes[random.indexes,]
		cross_res = cbind(c(Y.pred[pos.indexes],Y.pred[neg.indexes]),c(rep(1,nrow(pos.indexes)),rep(0,nrow(neg.indexes))))
		if(nrow(res)==1){
			res = cross_res;
		}else{
		  cat(dim(res),"\n")
		  cat(dim(cross_res),"\n")
			res = rbind(res,cross_res);
		}
	}
}

