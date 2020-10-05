# CSPP Sounder-QL

## Introduction

The Sounder Quicklooks (Sounder-QL) package is a python package for plotting
navigated satellite based sounder datasets derived from one of several software
packages that create level-2 sounder retrievals. It was developed within the
Community Satellite Processing Package, a project funded by the NOAA/NASA Joint
Polar Satellite System ([JPSS](https://www.jpss.noaa.gov)) project.

The CSPP Sounder-QL software is available from the CSPP website:

http://cimss.ssec.wisc.edu/cspp

Please use the ‘Contact Us’ form on the website to submit any questions or
comments about CSPP.

## Overview

This software package contains the scripts for generating navigated temperature,
water vapor, dewpoint and relative humidity plots, and skew-T diagrams. These
plots are generated from level-2 product files produced by the following
software packages...

* International ATOVS Processing Package (CSPP-IAPP)
* Microwave Integrated Retrieval System (MIRS)
* CSPP Hyperspectral Retrieval (HSRTV) Package
* NOAA Unique CrIS/ATMS Product System (NUCAPS)
* Hyper-Spectral Enterprise Algorithm Package (HEAP)

The package obtained from the CSPP website includes a self-contained
python distribution, `ShellB3`.

## System requirements

System requirements for the CSPP Sounder-QL v2.0 software are as follows:

- Intel or AMD CPU with 64-bit instruction support
- 1 GB RAM
- CentOS-6 (or newer) 64-bit Linux (or other compatible 64-bit Linux
distribution)
- 3 GB of disk space (plus space for your own DB data and CSPP-Sounder-QL
products)

## Acknowlegements

This package uses:
* [MetPy](https://unidata.github.io/MetPy/latest/index.html#) >= 0.12.2
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/#) >= 0.18.0

## Disclaimer

Original scripts and automation included as part of this package are distributed
under the GNU GENERAL PUBLIC LICENSE agreement version 3. Binary executable
files included as part of this software package are copyrighted and licensed by
their respective organizations, and distributed consistent with their licensing
terms. The University of Wisconsin-Madison Space Science and Engineering
Center (SSEC) makes no warranty of any kind with regard to the CSPP software
or any accompanying documentation, including but not limited to the implied
warranties of merchantability and fitness for a particular purpose. SSEC does
not indemnify any infringement of copyright, patent, or trademark through
the use or modification of this software. There is no expressed or implied
warranty made to anyone as to the suitability of this software for any purpose.
All risk of use is assumed by the user. Users agree not to hold SSEC, the
University of Wisconsin-Madison, or any of its employees or assigns liable for
any consequences resulting from the use of the CSPP software.

