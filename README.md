# Seg3IsBack
Github Repository for The R Project GSoc 2019 (Segmentor3IsBack package)

## CRAN Checks Fix
Currently, the [CRAN Package Check Results for Package Segmentor3IsBack](https://cloud.r-project.org/web/checks/check_results_Segmentor3IsBack.html) returns **WARN** for the following environments:
  * r-devel-linux-x86_64-debian-clang
  * r-devel-linux-x86_64-fedora-clang
  * r-release-osx-x86_64
  * r-oldrel-osx-x86_64

The **WARN** message being:
_Found the following significant warnings:
  BinNegative.cpp:108:5: warning: assigning field to itself [-Wself-assign-field]_

The Fix changes the statment 'T = T ;' (line: 108, src/BinNegative.cpp) to a comment so as to not change the internal working of the programme.
