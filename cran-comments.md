## Patch updates for release v1.0.38
This update contains improvments to the nonparanormal mapping of predictions for discrete variables

## Test environments
* OS X install: R 3.6.0
* win-builder: R-devel

## R CMD check results
There were no ERRORs or WARNINGs. 

An ERROR appears in the win-builder check for r-devel-windows-ix86+x86_64, likely owing to failures on r-devel-windows-ix86+x86_64 due to missing dependencies that need compilation. A test using `rhub::check( platform="windows-x86_64-devel", env_vars=c(R_COMPILE_AND_INSTALL_PACKAGES = "always"))` passes the win-builder check

Previously the version of this update was not sufficient for CRAN requirements: Insufficient package version (submitted: 1.0.4, existing: 1.0.37). This has been changed (new version is 1.0.38)

A URL in the one of the vignettes previously did not have https. The CRAN note was
* Found the following (possibly) invalid URLs: URL: http://www.jmlr.org/papers/v10/liu09a.html (moved to https://www.jmlr.org/papers/v10/liu09a.html) From: inst/doc/Gaussian_Poisson_CRFs.html Status: 200

This has been changed to the suggested https format

Previously a DOI was not formatted according to CRAN requirements: Found the following URLs which should use \doi (with the DOI name only): File 'Bird.parasites.Rd': http://dx.doi.org/10.5061/dryad.pp6k4

This has been fixed


Maintainer: 'Nicholas J Clark <nicholas.j.clark1214@gmail.com>'

