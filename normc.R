# normc.R: Residualize c with d' to get normalized c values

# Created: 02/16/2020 By: Evan Layher (normc.R)
# Revised: 07/18/2020 By: Evan Layher (updated for online upload)
#
## --- LICENSE INFORMATION --- ##
## Modified BSD-2 License - for Non-Commercial Use Only

## Copyright (c) 2020, The Regents of the University of California
## All rights reserved.

## Redistribution and use in source and binary forms, with or without modification, are 
## permitted for non-commercial use only provided that the following conditions are met:

## 1. Redistributions of source code must retain the above copyright notice, this list 
##    of conditions and the following disclaimer.

## 2. Redistributions in binary form must reproduce the above copyright notice, this list 
##    of conditions and the following disclaimer in the documentation and/or other 
##    materials provided with the distribution.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
## EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
## OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
## SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
## INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
## TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
## OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
## WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## For permission to use for commercial purposes, please contact UCSB’s Office of 
## Technology & Industry Alliances at 805-893-5180 or info@tia.ucsb.edu.
## --------------------------- ##
#
# REQUIRED INPUTS: 2 vectors of equal length
# [1] cVec: vector of c values
# [2] dVec: vector of d' values
#
# Outputs:
# [1] cn: normalized c values
# [2] resid: residuals
# [3] slope: regression slope
# [4] int: regression intercept
# [5] r: Pearson's r
# [6] r2: Person's r squared

normc <- function(cVec = NULL, dVec = NULL) {
  
  # Ensure valid inputs (1) vectors must contain values; (2) vectors must be equal length
  if (is.null(cVec) | is.null(dVec)) {stop('Must specify input vectors "cVec" and "dVec"')}
  if (length(cVec) != length(dVec)) {stop('Vector lengths must be equal: "cVec" (', length(cVec), ') "dVec" (', length(dVec), ')')}
  
  mC <- mean(cVec) # mean of c vector
  mD <- mean(dVec) # mean of d' vector
  
  cSq <- sum((cVec - mC) ^ 2) # sum of squares: c (y)
  dSq <- sum((dVec - mD) ^ 2) # sum of squares: d' (x)
  cdSq <- sum((dVec - mD) * (cVec - mC)) # sum of squares: c * d' (y)
  
  slope <- cdSq / dSq # regression slope
  int <- mC - (slope * mD) # regression intercept
  r <- cdSq / sqrt((cSq * dSq)) # Pearson's r
  r2 <- r ^ 2 # explained variance (%)
  ln <- slope * dVec + int # y = mx + b
  resid <- cVec - ln # residuals
  cn <- resid + mC # normalized c (residuals + mean c)
  
  return(cbind(cn, resid, slope, int, r, r2))
} # normc
