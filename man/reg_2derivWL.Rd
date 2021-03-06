% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npregderiv.R
\name{reg_2derivWL}
\alias{reg_2derivWL}
\title{Estimating the second derivative of a regression function by the method of Wang and Lin (2015)}
\usage{
reg_2derivWL(xdat, ydat, k, u)
}
\arguments{
\item{xdat}{numerical vector of the increasing and evenly spaced design points}

\item{ydat}{numerical vector of the corresponding data points}

\item{k}{integer value of the smoothing parameter}

\item{u}{an estimation point (scalar) that must be one of the \eqn{xdat} points from the interior region}
}
\value{
The computed second derivative estimate (scalar) at the point \eqn{u}.
}
\description{
The design data \eqn{xdat} must be increasing and evenly spaced. The second derivative estimate is computed at \eqn{u} that must be one of the design points \eqn{xdat} from the interior region. The interior region excludes the first \eqn{k} points with the smallest \eqn{xdat} values (left boundary) and the last \eqn{k} points with the largest \eqn{xdat} values (right boundary).
}
\details{
The method's smoothing parameter \eqn{k} shows how many data points are used from each side of an estimation point \eqn{u} to compute the second derivative estimate at that point. Choose \eqn{k<n/4}, where \eqn{n} is the sample size. The value of \eqn{k=0.1n} can be taken as a starting point.
}
\examples{
# EXAMPLE 1 (simulated data).
# The regression function is taken from p. 2631 of the article of Wang and Lin (2015).
m=function(x)               # Regression function from the article of Wang and Lin (2015)
  32*exp(-8*(1-2*x)^2)*(1-2*x)
m2=function(x)              # second derivative of m
  -2048*exp(-8*(1-2*x)^2)*(128*x^3-192*x^2+90*x-13)
N=100                       # sample size
xi=1:N/N                    # equidistant design points
sigma=0.2                   # noise level used in the article
yi=m(xi)+rnorm(N,sd=sigma)  # generating data
K=0.1*N                     # selected value of the smoothing parameter k
x_inter=xi[(K+1):(N-K)]     # interior estimation region
N_inter=length(x_inter)     # number of points where the second derivative estimate is computed
M2=array(0,N_inter)         # Vector of estimated second derivatives at the points x_inter
for (i in 1:N_inter)
  M2[i]=reg_2derivWL(xi,yi,K,x_inter[i])
op=par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(xi,yi,pch=20,cex=1.7,main="Regression function and generated data",xlab="",ylab="",
     cex.axis=1.8,cex.main=1.5)
lines(xi,m(xi),'l',col="red",lwd=2)
title(xlab="argument",ylab="regression function",line=2.5,cex.lab=1.8)
legend(0.35,5,lty="solid",col="red",lwd=2,legend="true regression function",bty="n",
       cex=1.45,seg.len=0.5)
plot(x_inter,M2,'l',lwd=2,main="Actual and estimated second derivative",xlab="",ylab="",
     cex.axis=1.8,cex.main=1.5)
lines(x_inter,m2(x_inter),'l',lwd=2,col="red")
title(xlab="argument",ylab="second derivative",line=2.5,cex.lab=1.8)
legend(0.01,800,lty=c("solid","solid"),col=c("red","black"),lwd=c(2,2),
      legend=c("true second derivative", "estimate"),bty="n",cex=1.45,seg.len=0.5)
par(op)

\donttest{
# EXAMPLE 2: Estimating the second derivative of the deflection function for the data on
# deflection after film deposition scanned at 90 degrees.
# See Savchuk O., and Volinsky, A. (2020) for the experiment's description.
xdesign=wafer40$x_dat
ydata=wafer40$y_after_90
n_data=length(xdesign)      # sample size (original)
n_new=1034                  # reducing the number of points where the second 
                            # derivative is estimated.
K=ceiling(0.1*n_data)       # value of the smoothing parameter
index=(K+1)+0:(n_new-1)*6   # cutting about 10\% of points from each side and reducing the 
                            # number of points where the second derivative is estimated
x_inter=xdesign[index]      # the values of argument where the second derivative is to be
                            # estimated
y_inter=ydata[index]
Der2=array(0,n_new)
for (i in 1:n_new)
  Der2[i]=reg_2derivWL(xdesign,ydata,K,x_inter[i])
op=par(no.readonly=TRUE)
par(mfrow=c(1,2))
plot(xdesign*1000,ydata*10^6,pch=20,main="Deflection AFTER film deposition. Angle: 90 degrees.",
     xlab="",ylab="",cex.axis=1.8,cex.main=1.5)
title(ylab=expression(paste("deflection, ", mu, "m")), xlab="position, mm", line=2.25, 
      cex.lab=1.8)
plot(x_inter*1000,Der2,'l',lwd=2,
     main="2-nd derivative AFTER film deposition. Angle: 90 degrees.",xlab="",ylab="",
     cex.axis=1.8, cex.main=1.5)
title(ylab="second derivative of deflection, 1/m", xlab="position, mm", line=2.5,
      cex.lab=1.8)
par(op)
}
}
\references{
Wang, W.W., Lin, L. (2015) \url{http://www.jmlr.org/papers/v16/wang15b.html} ``Derivative Estimation Based on Difference Sequence via Locally
Weighted Least Squares Regression'', \emph{Journal of Machine Learning Research}, 16, 2617-2641.

Savchuk, O., Volinsky, A.A. (2020). ``Nonparametric estimation of SiC film residual stress
from the wafer profile''.
}
\seealso{
\code{\link{reg_1derivWL}}, \code{\link{wafer40}}
}
