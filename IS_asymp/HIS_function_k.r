HIS.est <- function(y, beta, k) {
	N <- ncol(y)
	T <- nrow(y)
	b.orth <- Null(beta)
	# VECM, known beta
	z <- y %*% beta
	if (ncol(beta) == 1) lz <- c(NA, z[-nrow(z)]) else {
		lz <- rbind(NA, z[-nrow(z),]) }
	ly <- rbind(NA, y[-nrow(y),])
	dy <- y - ly
	fm.data <- as.data.frame(cbind(dy, lz))
	fm <- lm(dy ~ -1 + lz)
	if (k > 1) {
		for (j in 1:(k-1)) {
			if (j == 1)	ldy1 <- rbind(NA, dy[-nrow(dy),]) else {
				eval(parse(text=sprintf("ldy%s <- rbind(NA, ldy%s[-nrow(ldy%s),])", j, j-1, j-1)))
				}
			eval(parse(text=sprintf("fm.data <- cbind(fm.data,ldy%s)", j)))
			}
		colnames(fm.data) <- c(paste("dy", 1:N, sep=""), paste("lz", 1:(N-1), sep=""),
			paste(rep(paste("ldy", 1:(k-1), sep=""), each=N), 1:N, sep=""))
		dy <- as.matrix(fm.data[,1:N])
		fm.data <- fm.data[,-(1:N)]	
		fm <- lm(dy ~ -1 + ., data=fm.data)
		}
	alpha.hat <- t(coef(fm)[1:(N-1),])
	if (ncol(beta) == 1) alpha.hat <- matrix(coef(fm)[1:(N-1),], nrow=N)
	alpha.bar <- alpha.hat %*% solve(t(alpha.hat) %*% alpha.hat)
	if (k == 1) Gamma.hat <- matrix(0, nrow=N, ncol=N)
	if (k == 2)	Gamma.hat <- coef(fm)[-(1:(N-1)),]
	if (k > 2) {
		Gamma.hat <- coef(fm)[-(1:(N-1)),]
		temp.mat <- matrix(NA, nrow=N, ncol=N)
		for (j in 1:N) {
			temp.mat[j,] <- colSums(Gamma.hat[seq(j,by=N,length=k-1),])
			}
		Gamma.hat <- temp.mat
		}
	Psi.hat <- diag(N) - Gamma.hat
	u.hat <- resid(fm)
	Omega.hat <- cov(u.hat)

	a.hat.orth <- Null(alpha.hat)
	a.hat.orth <- a.hat.orth/sum(a.hat.orth)
	C.hat <- b.orth %*% solve(t(a.hat.orth) %*% (diag(N) - Gamma.hat) %*% b.orth) %*% t(a.hat.orth)

	Q1 <- (t(C.hat) %*% t(Psi.hat) - diag(N)) %*% alpha.bar
	Q.star <- cbind(Q1,t(C.hat))	

	ldy1 <- rbind(NA, dy[-nrow(dy),])
	s.hat <- cbind(lz,ldy1)
	M.ss <- 1/T* (t(na.omit(s.hat)) %*% na.omit(s.hat))
	kappa <- Q.star %*% solve(M.ss) %*% t(Q.star)
	if (k > 2) {
		d2y <- diff(dy)
		w2 <- rbind(NA,d2y)
		if (k > 3) {
			for (j in 2:(k-1)) {
				d2y.lag <- rbind(matrix(NA,nrow=j, ncol=N),head(d2y,-(j-1)))
				w2 <- cbind(w2, d2y.lag)
				}
			}
		# adjust to same length
		temp <- na.omit(cbind(s.hat,w2))
		s.hat <- temp[,1:(2*N-1)]
		w2 <- temp[,-c(1:(2*N-1))]

		M.ss <- 1/T* (t(s.hat) %*% s.hat)
		M.sw2 <- 1/T* (t(s.hat) %*% w2)
		M.w2w2 <- 1/T* (t(w2) %*% w2)
		M.w2s <- 1/T* (t(w2) %*% s.hat)
		M.ss.w2 <- M.ss - M.sw2 %*% solve(M.w2w2) %*% M.w2s

		kappa <- Q.star %*% solve(M.ss.w2) %*% t(Q.star)
	}

	# information share
	xi.hat <- solve(t(a.hat.orth) %*% (diag(N) - Gamma.hat) %*% b.orth) %*% t(a.hat.orth)
	F.hat <- t(chol(Omega.hat))
	IS.num <- xi.hat %*% F.hat
	IS.denom <- xi.hat %*% Omega.hat %*% t(xi.hat)
	IS <- IS.num^2 / as.vector(IS.denom)
	names(IS) <- paste("y", 1:N, sep="")

	n_vech <- N*(N+1)/2
	Vmat <- t( apply(u.hat, 1, function(e){
		M <- tcrossprod(e) - Omega.hat
		M[ lower.tri(M, diag=TRUE) ]
		}) )           # T x n_vech
	S.Omega <- crossprod(Vmat) / T   # n_vech x n_vech

	# define f_chol: vech(Omega) -> vec( L )  where L = t(chol(Omega))
	f_chol <- function(v){
		# recover N from length(v), 
		# solve \[ \ell = \text{length}(v) \;=\;\tfrac12\,N(N+1)\,. \]
		N0 <- (-1 + sqrt(1 + 8*length(v))) / 2
		Om <- matrix(0, N0, N0)
		Om[ lower.tri(Om, diag=TRUE) ] <- v
		Om <- Om + t(Om) - diag(diag(Om))
		L  <- t(chol(Om))     
		as.numeric(L)          # stacks columns into a length-N0^2 vector
	}

	v0 <- Omega.hat[ lower.tri(Omega.hat, diag=TRUE) ]
	H.hat <- jacobian(f_chol, v0)
	S.F <- H.hat %*% S.Omega %*% t(H.hat)

	A <- t(F.hat)
	B <- diag(N) %x% xi.hat
	S.xi <- ((xi.hat %*% Omega.hat %*% t(xi.hat)) %x% kappa)
	# symmetrically-distributed errors, cross-terms zero
	Z <- (A %*% S.xi %*% t(A) + B %*% S.F %*% t(B)) / T

	theta <- xi.hat %*% F.hat
	denom <- sum(theta^2)

	# Jacobian
	Jg <- matrix(NA, nrow=N, ncol=N)
	for (j in 1:N) {
		for (i in 1:N) {
			if (i == j) {
				theta.sum <- denom - theta[i]^2
				Jg[i,j] <- 2*theta[j]*theta.sum/denom^2
			} else {
				Jg[i,j] <- -2*theta[j]^2*theta[i]/denom^2 }
			}
		}

	# covariance matrix
	S.IS <- t(Jg) %*% Z %*% Jg

	results <- list(IS=IS, S.IS=S.IS, xi.hat=xi.hat, S.xi=S.xi/T, 
			Omega.hat=Omega.hat, S.Omega=S.Omega/T, F.hat=F.hat, S.F=S.F/T)
	return(results)
	}





