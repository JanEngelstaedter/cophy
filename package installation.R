devtools::document()
devtools::check()
devtools::build()

devtools::install("~/Dropbox/Projects/Host-shift dynamics/Package\ development/cophy.dev/cophy", build_vignettes = TRUE)
library(cophy)
browseVignettes(package = "cophy")

citation("ape")
packageVersion("ape")
packageVersion("cophy")

library(ape)


htree <- rphylo(10,1,0.5, fossils = TRUE)
htree <- rcoal(10)
plot(htree)
cop <- rcophylo(HTree = htree)
plot(cop)


set.seed(2)
htree <- rphylo_H(tmax = 5, lambda = 1, mu = 0.4)
coph <- rcophylo(HTree = htree, PStartT = 1)
plot(coph)

coph[[2]]
