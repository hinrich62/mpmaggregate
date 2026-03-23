pkgname <- "mpmaggregate"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "mpmaggregate-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('mpmaggregate')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("dominant_eigen")
### * dominant_eigen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dominant_eigen
### Title: Dominant eigen-elements of a population projection matrix
### Aliases: dominant_eigen

### ** Examples

A <- matrix(c(
  0,   1,   2,
  0.5, 0,   0,
  0,   0.8, 0.9
), nrow = 3, byrow = TRUE)

dom <- dominant_eigen(A)
dom$lambda
sum(dom$w)
sum(dom$v * dom$w)

#reducible example
A<-rbind(c(1,1,0),
         c(1,0,0),
         c(0,1,0))

dom <- dominant_eigen(A)
dom$lambda
sum(dom$w)
sum(dom$v * dom$w)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dominant_eigen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generation_time")
### * generation_time

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generation_time
### Title: Generation time from a matrix population model
### Aliases: generation_time

### ** Examples

matU <- matrix(c(0.2, 0.0,
                0.3, 0.4), nrow = 2, byrow = TRUE)
matF <- matrix(c(0.0, 1.2,
                0.0, 0.0), nrow = 2, byrow = TRUE)
generation_time(matF = matF, matU = matU, framework = "lambda")
generation_time(matF = matF, matU = matU, framework = "R0")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generation_time", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_leslie")
### * is_leslie

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_leslie
### Title: Test whether a matrix is a Leslie matrix
### Aliases: is_leslie

### ** Examples

L <- matrix(c(
  0, 2, 1,
  0.6, 0, 0,
  0, 0.7, 0
), nrow = 3, byrow = TRUE)
is_leslie(L)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_leslie", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("leslie_aggregate")
### * leslie_aggregate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: leslie_aggregate
### Title: Aggregate Leslie-to-Leslie matrix population model
### Aliases: leslie_aggregate

### ** Examples

# A simple 3x3 Leslie matrix (fertility rates in first row; survival probabilities on subdiagonal)
A <- matrix(c(
  0.0, 1.2, 1.8,
  0.5, 0.0, 0.0,
  0.0, 0.7, 0.0
), nrow = 3, byrow = TRUE)

# Aggregate to 2 age groups
res_std <- leslie_aggregate(
  matA = A,
  ngroups = 2,
  framework = "lambda",
  criterion = "standard"
)
res_std$matAk_agg
res_std$effectiveness

res_el <- leslie_aggregate(
  matA = A,
  ngroups = 2,
  framework = "lambda",
  criterion = "elasticity"
)
res_el$matAk_agg
res_el$effectiveness




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("leslie_aggregate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("leslie_disaggregate")
### * leslie_disaggregate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: leslie_disaggregate
### Title: Disaggregate a Leslie matrix population model to a compatible
###   dimensionality
### Aliases: leslie_disaggregate

### ** Examples

L <- matrix(c(
  0, 2, 1,
  0.6, 0, 0,
  0, 0.7, 0
), nrow = 3, byrow = TRUE)
leslie_disaggregate(L, m = 2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("leslie_disaggregate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("leslie_dominant_eigen")
### * leslie_dominant_eigen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: leslie_dominant_eigen
### Title: Dominant eigen-elements of a Leslie matrix
### Aliases: leslie_dominant_eigen

### ** Examples

leslie_dominant_eigen(1)

L <- matrix(c(
  0, 2, 1,
  0.6, 0, 0,
  0, 0.7, 0
), nrow = 3, byrow = TRUE)
dom <- leslie_dominant_eigen(L)
dom$lambda
sum(dom$w)
sum(dom$v * dom$w)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("leslie_dominant_eigen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("leslie_reproductive_values")
### * leslie_reproductive_values

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: leslie_reproductive_values
### Title: Reproductive values for a Leslie matrix
### Aliases: leslie_reproductive_values

### ** Examples

leslie_reproductive_values(1)

L <- matrix(c(
  0, 2, 1,
  0.6, 0, 0,
  0, 0.7, 0
), nrow = 3, byrow = TRUE)
v <- leslie_reproductive_values(L)
w <- leslie_stable_age(L)
sum(v * w)  # should be 1




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("leslie_reproductive_values", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("leslie_stable_age")
### * leslie_stable_age

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: leslie_stable_age
### Title: Stable age distribution for a Leslie matrix
### Aliases: leslie_stable_age

### ** Examples

leslie_stable_age(1)

L <- matrix(c(
  0, 2, 1,
  0.6, 0, 0,
  0, 0.7, 0
), nrow = 3, byrow = TRUE)
leslie_stable_age(L)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("leslie_stable_age", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mpm_aggregate")
### * mpm_aggregate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mpm_aggregate
### Title: Aggregate a general-to-general matrix population model
### Aliases: mpm_aggregate

### ** Examples

# Example aggregation of a 3x3 projection matrix to 2x2 using groups:
# group 1 = stage 1; group 2 = stages 2 and 3.
A <- matrix(c(
  0.0, 1.0, 2.0,
  0.5, 0.0, 0.0,
  0.0, 0.8, 0.9
), nrow = 3, byrow = TRUE)

res_std <- mpm_aggregate(
  matA = A,
  groups = list(c(1), c(2, 3)),
  framework = "lambda",
  criterion = "standard"
)
res_std$matA_agg
res_std$effectiveness

res_el <- mpm_aggregate(
  matA = A,
  groups = list(c(1), c(2, 3)),
  framework = "lambda",
  criterion = "elasticity"
)
res_el$matA_agg
res_el$effectiveness




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mpm_aggregate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mpm_elasticity")
### * mpm_elasticity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mpm_elasticity
### Title: Elasticity of lambda or R_0 with respect to entries of 'matA'
### Aliases: mpm_elasticity

### ** Examples

## Lambda framework: matA provided directly
matA <- matrix(
  c(0.2, 1.2,
    0.3, 0.4),
  nrow = 2, byrow = TRUE
)
out_lambda <- mpm_elasticity(matA = matA, framework = "lambda")
str(out_lambda)

## R0 framework: matA constructed from matF and matU
matU <- matrix(
  c(0.2, 0.0,
    0.3, 0.4),
  nrow = 2, byrow = TRUE
)
matF <- matrix(
  c(0.0, 1.2,
    0.0, 0.0),
  nrow = 2, byrow = TRUE
)
out_R0 <- mpm_elasticity(matF = matF, matU = matU, framework = "R0")
str(out_R0)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mpm_elasticity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mpm_partition")
### * mpm_partition

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mpm_partition
### Title: Create a partitioning matrix for MPM aggregation
### Aliases: mpm_partition

### ** Examples

g <- list(c(1), c(2:3))
P <- mpm_partition(g, n = 3)
P




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mpm_partition", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("reproductive_values")
### * reproductive_values

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: reproductive_values
### Title: Compute reproductive values
### Aliases: reproductive_values

### ** Examples

A <- matrix(c(
  0, 1,
  0.5, 0
), nrow = 2, byrow = TRUE)

w <- stable_stage(A)
v <- reproductive_values(A)
v
sum(v * w)  # should be 1

A<-rbind(c(1,1,0),
         c(1,0,0),
         c(0,1,0))

w <- stable_stage(A)
w
sum(w)

w <- stable_stage(A)
v <- reproductive_values(A)
v
sum(v * w)  # should be 1




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("reproductive_values", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spectral_radius")
### * spectral_radius

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spectral_radius
### Title: Compute the spectral radius of a matrix
### Aliases: spectral_radius

### ** Examples

A <- matrix(c(
  0, 1,
  0.5, 0
), nrow = 2, byrow = TRUE)

spectral_radius(A)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spectral_radius", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stable_stage")
### * stable_stage

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stable_stage
### Title: Compute the stable stage distribution
### Aliases: stable_stage

### ** Examples

#irreducible example
A <- matrix(c(
  0, 1,
  0.5, 0
), nrow = 2, byrow = TRUE)

w <- stable_stage(A)
w
sum(w)
#reducible example
A<-rbind(c(1,1,0),
         c(1,0,0),
         c(0,1,0))

w <- stable_stage(A)
w
sum(w)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stable_stage", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
