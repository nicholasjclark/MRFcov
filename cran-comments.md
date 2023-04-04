## Patch updates for release v1.0.39
This update contains improvements code (changing `ncol` to `NCOL` where appropriate) and uses Bernoulli draws for more realistic binary predictions

## Test environments
* Windows install: R 4.2.1
* win-builder: R-devel

## R CMD check results
There were no ERRORs or WARNINGs. 

An error was found when multiple processes were spawned during `donttest` runs on `CRAN`: `Error in .check_ncores(length(names)) : 3 simultaneous processes spawned`. This has been fixed

Maintainer: 'Nicholas J Clark <nicholas.j.clark1214@gmail.com>'

