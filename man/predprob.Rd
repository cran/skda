\name{predprob}
\alias{predprob}
\title{
Conditional class probability prediction for SKDA
}
\description{
This function predicts the conditional class probability for the sparse (multicategory) kernel discriminant analysis (\code{\link{skda}}).
 }
\usage{
predprob(x, y, lam, xnew, method = "Bayes")
}
\arguments{
  \item{x}{
     a matrix ( n X p ) that contains predictors.
  }
  \item{y}{
     a vector that contains the categorical response coded as 1, 2, ..., K.
  }
  \item{lam}{
    a length-K vector of inverse smoothing bandwidths.
  }
  \item{xnew}{
    a matrix (with p columns) that contains predictors of new observations.
  }
  \item{method}{
     method (mle or Bayes) to be used in the KDA classifier.
  }
}
\value{
    \item{phat}{
    a matrix (with K columns) that contains predicted conditional class probabilities.
    }
}
\references{
L. A. Stefanski, Y. Wu, and K. White (2013) Variable selection in nonparametric classification via measurement error model selection likelihoods
\emph{Journal of the American Statistical Association}, \bold{??},
???-???.



}
\author{L. A. Stefanski, Y. Wu, and K. White}
\seealso{
\code{\link{skda}}, 
\code{\link{cvskda}}
}
\examples{


        # a binary example
          rm(list=ls())
          n=200
          p=20
          r=0.8
          x=matrix(rnorm(n*p,mean=0,sd=1),n,p)
          y=rbinom(n,1,0.5)
          x[,3]=x[,3]+(2*y-1)*r
          x[,11]=x[,11]+(2*y-1)*r

          y=y+1

          ind1=which(y>1.5)
          ind0=which(y<1.5)
          plot(-4:4, -4:4, type = "n")
          points(x[ind1, 3], x[ind1,11],col="blue")
          points(x[ind0, 3], x[ind0,11],col="red")

          lam=skda(x,y,3)$lam
         

          xtest=matrix(rnorm(10*n*p,mean=0,sd=1),10*n,p)
          ytest=rbinom(10*n,1,0.5)
          xtest[,3]=xtest[,3]+(2*ytest-1)*r
          xtest[,11]=xtest[,11]+(2*ytest-1)*r
          ytest=ytest+1

          ptesthat=predprob(x,y,lam, xtest)
          ytesthat=apply(ptesthat, 1,which.max)
          error=sum(abs(ytesthat-ytest)>0.5)
          error/length(ytest)





        # a three-class example
          rm(list=ls())
          n=200
          p=20
          r=2
          x=matrix(rnorm(n*p,mean=0,sd=1),n,p)
          y=ceiling(runif(n,0,3))
          thetas=c(0, 2*pi/3, 4*pi/3)
          x[,3]=x[,3]+r*cos(thetas[y])
          x[,11]=x[,11]+r*sin(thetas[y])

          ind1=which(y==1)
          ind2=which(y==2)
          ind3=which(y==3)

          plot(-6:6, -6:6, type = "n")
          points(x[ind1, 3], x[ind1,11],col="blue")
          points(x[ind2, 3], x[ind2,11],col="red")
          points(x[ind3, 3], x[ind3,11],col="black")



          lam=skda(x,y,3)$lam
         

          xtest=matrix(rnorm(10*n*p,mean=0,sd=1),10*n,p)
          ytest=ceiling(runif(10*n,0,3))
          xtest[,3]=xtest[,3]+r*cos(thetas[ytest])
          xtest[,11]=xtest[,11]+r*sin(thetas[ytest])


          ptesthat=predprob(x,y,lam, xtest)
          ytesthat=apply(ptesthat, 1,which.max)
          error=sum(abs(ytesthat-ytest)>0.5)
          error/length(ytest)



}
\keyword{file}

