#!/bin/bash
rm -f src/RcppExports.cpp
R CMD INSTALL --no-multiarch --with-keep.source ../bayesEpi
