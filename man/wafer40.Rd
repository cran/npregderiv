% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npregderiv.R
\docType{data}
\name{wafer40}
\alias{wafer40}
\title{The data set on the \eqn{n}=7755 substrate's deflection measurements before and after the thin film deposition in two radial directions.}
\format{
The data frame has five variables:
\describe{
\item{x_dat}{design points (separated by 0.000005 m);}
\item{y_before_0}{deflections measurements (in m) before film deposition at 0 degrees;}
\item{y_after_0}{deflections measurements (in m) after film deposition at 0 degrees;}
\item{y_before_90}{deflections measurements (in m) before film deposition at 90 degrees;}
\item{y_after_90}{deflections measurements (in m) after film deposition at 90 degrees.}
}
}
\source{
Savchuk, O., Volinsky, A.A. (2020). ``Nonparametric estimation of SiC film residual stress from the wafer profile''.
}
\usage{
wafer40
}
\description{
The data frame on the deflection measurements is described and analyzed in the article of Savchuk and Volinsky (2020).
}
\examples{
# EXAMPLE: Plotting scatter plot of deflection in the case of 0 degrees before the film
# deposition. Note: position is displayed in mm, whereas the deflection is displayed in
# micrometers. No smoothing is used. The original n=7755 data points are plotted.
dev.new()
plot(wafer40$x_dat*1000,wafer40$y_before_0*1000000,pch=20,
     main="Deflection BEFORE film deposition. Angle: 0 degrees.",xlab="",ylab="",
     cex.axis=2,cex.main=2)
title(ylab=expression(paste("deflection, ", mu, "m")), xlab="position, mm",
      line=2.25, cex.lab=2)
}
\seealso{
\code{\link{reg_2derivWL}}.
}
\keyword{datasets}
