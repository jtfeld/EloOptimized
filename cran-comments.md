## Test environments
* local Windows 10 install, R 4.1.0
* ubuntu 16.04 (on travis-ci), R-release.
* macOS-latest (on github actions), R-release, R-4.0, R-devel.
* r-hub
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Joseph Feldblum <feldblum@umich.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Elo (2:18, 10:5)
  Foerster (10:32)
  al (10:51)
  et (10:48)
  
  These words are spelled correctly
  
## Resubmission
This is a resubmission. In this version I have:

* Changed the Description field so that it does not start with 
  'This package'.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* I have run R CMD check on the 7 downstream dependencies.
  (Summary at https://github.com/jtfeld/EloOptimized/blob/master/revdep/README.md). 
  
* FAILURE SUMMARY

* No failures.
