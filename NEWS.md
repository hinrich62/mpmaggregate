# mpmaggregate 0.1.1
* Initial release. 
* Adds mpm_aggregate() and leslie_aggregate() for standard and 
* elasticity-consistent aggregation in the lambda and R0 frameworks.

\section{mpmaggregate 0.1.1}{

  \subsection{New features}{
    \itemize{
      \item Added \code{generation_time()}, a function to compute generation time
      from matrix population models using either a \eqn{\lambda}-based or
      \eqn{R_0}-based framework.

      \item Added \code{elasticity()}, a function to compute elasticities of the
      dominant population growth rate (\eqn{\lambda}) or net reproductive rate
      (\eqn{R_0}) with respect to entries of matrix population models, using
      framework-specific input validation.
    }
  }

}
