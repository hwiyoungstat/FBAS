#' FBAS
#'
#' @import CVXR
#' @param x Exposure
#' @param y Outcome
#' @param m mediator
#' @param lambda regularizer
#' @param rho lagrangian penalty
#'
#' @return
#' @export


# Mediation Lasso
FBAS <- function(x,y,m,lambda,rho){
  Eps.rel <- 1e-3
  Eps.abs <- 1e-1
  Max.ite <- 100

  N <- dim(m)[1]
  P <- dim(m)[2]
  y.ast <- scale(resid(lm(y~x)))
  x <- scale(x)
  y <- scale(y)
  m <- apply(m,2,scale)
  sgn.tau <- sign(coef(lm(y~x))[2])


  Proj_compM <- (diag(N)-x%*%t(x)/(N-1))%*%m

  w <- rep(1,P)
  u <- rep(0.01,P)
  t <- rep(0,P)

  I <- diag(P)

  MXY <- (t(m)%*%x%*%t(y.ast)%*%Proj_compM)/((N-1)^2)
#  Quad <- 1/4*(2*rho*I - MXY - t(MXY))
  if(sgn.tau>0){Quad <- 1/4*(2*rho*I - MXY - t(MXY))}else{Quad <- 1/4*(2*rho*I + MXY + t(MXY))}

  for(i in 1:Max.ite){
    # Use CVXR for SDR
    Cost <- rbind(cbind(Quad,rho/2*(t-u)),c(rho/2*(t-u),0))
    W <- Variable(P+1,P+1,PSD=TRUE)
    obj <- matrix_trace(Cost%*%W)
    constr <- list(matrix_trace((diag(c(rep(0,P),1)))%*%W)==1,
                   matrix_trace((rbind(cbind(t(m)%*%m,0),0))%*%W)==1,
    #              matrix_trace((rbind(cbind(t(Proj_compM)%*%Proj_compM,0),0))%*%W)==1,
                   matrix_trace(rbind(cbind(diag(rep(0,P)),t(m)%*%rep(1,N)/(2*N)),c(t(rep(1,N))%*%m/(2*N),0))%*%W)==0,
                   matrix_trace((diag(c(rep(1,P),0)))%*%W)<=1) # L2 constraint for stability
    prob <- Problem(Minimize(obj),constr)
    result <- solve(prob)
    Optimal.W <- result[[1]]
    sign <- 2*(Optimal.W[P+1,P+1]>0)-1 # if t=1 -> w: optimal, if t=-1 -> -w; optimal

    Eigen <- eigen(Optimal.W)
    w.new <- sign*sqrt(Eigen$values[1])*Eigen$vectors[c(1:P),1]

    ## Update U
    u.new <- soft(w.new+t,lambda/rho)

    # Define Residual
    n.resid.prim <- sqrt(t(c(w.new-u.new))%*%c(w.new-u.new))
    n.resid.dual <- sqrt(t(rho*c(u-u.new))%*%(rho*c(u-u.new)))


    # Update rho
    Mu <- 10
    tau.inc <- 2
    tau.dec <- 2

    rho.new <- ifelse(n.resid.prim>Mu*n.resid.dual,tau.inc*rho,
                      ifelse(n.resid.dual>Mu*n.resid.prim,rho/tau.dec,rho))

    scale <- rho.new/rho

    t.new <- t + w.new - u.new

    # Stopping criterion
    Eps.prime <- sqrt(P)*Eps.abs + Eps.rel*max(sqrt(t(w.new)%*%w.new),sqrt(t(u.new)%*%u.new))
    Eps.dual <- sqrt(P)*Eps.abs + Eps.rel*norm(sqrt(t(t.new)%*%t.new))


    if(n.resid.prim<Eps.prime & n.resid.dual<Eps.dual) break



    # Update for the next step
    w <- w.new
    u <- u.new
    t <- t.new
  }

  result <- list("U"=u.new)
  return(result)
}

