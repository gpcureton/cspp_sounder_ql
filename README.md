# CSPP-Sounder-QL

## Introduction

### Overview

This document contains instructions for installation and operation of the Community Satellite
Processing Package (CSPP) release of the Sounder Quicklooks Package (Sounder-QL)
software for plotting navigated sounder datasets derived from one of several software packages
that create level-2 sounder retrievals.

The CSPP-Sounder-QL software is available from the CSPP website:

http://cimss.ssec.wisc.edu/cspp

Software, test data, and documentation may be downloaded from this web site. Please use the
‘Contact Us’ form on the website to submit any questions or comments about CSPP.

### System requirements

System requirements for the CSPP-Sounder-QL v1.0 software are as follows:

- Intel or AMD CPU with 64-bit instruction support,
- 1 GB RAM
- CentOS-6 64-bit Linux (or other compatible 64-bit Linux distribution),
- 3 GB of disk space (plus space for your own DB data and CSPP-Sounder-QL products).
Linux terminal commands included in these instructions are for the bash shell.

Linux terminal commands included in these instructions are for the bash shell.

### Acknowledgements

This package relies heavily on various projects, which include:

- MetPy v0.7.0 (https://doi.org/10.5065/D6WW7G29)
- Matplotlib v2.1.2 (https://doi.org/10.5281/zenodo.1154287)
- xarray v0.10.8 (https://doi.org/10.5281/zenodo.598201)
- netCDF4-python v1.3.1 (http://unidata.github.io/netcdf4-python/)
- h5py v2.7.1 (https://github.com/h5py)

### Disclaimer

Original scripts and automation included as part of this package are distributed under the GNU
GENERAL PUBLIC LICENSE agreement version 3. Binary executable files included as part of
this software package are copyrighted and licensed by their respective organizations, and
distributed consistent with their licensing terms.
The University of Wisconsin-Madison Space Science and Engineering Center (SSEC) makes no
warranty of any kind with regard to the CSPP software or any accompanying documentation,
including but not limited to the implied warranties of merchantability and fitness for a particular
purpose. SSEC does not indemnify any infringement of copyright, patent, or trademark through
the use or modification of this software.
There is no expressed or implied warranty made to anyone as to the suitability of this software
for any purpose. All risk of use is assumed by the user. Users agree not to hold SSEC, the
University of Wisconsin-Madison, or any of its employees or assigns liable for any
consequences resulting from the use of the CSPP software.

## Installation and Configuration

### Overview

This software package contains the CSPP-Sounder-QL plotting scripts, which generate
navigated temperature, dewpoint and relative humidity plots, and skew-T diagrams. These plots
are generated from level-2 product files produced by the following software packages...

* CSPP-International ATOVS Processing Package (CSPP-IAPP)
* Microwave Integrated Retrieval System (MIRS)
* CSPP Hyperspectral Retrieval (HSRTV) Package
* NOAA Unique CrIS/ATMS Product System (NUCAPS)

### Installation of CSPP-Sounder-QL Software

Create the installation directory CSPP_Sounder_QL_v1.2, and the working and tarball directories,
and change into the tarball directory...

```bash
mkdir -p CSPP_Sounder_QL_v1.2/tarballs CSPP_Sounder_QL_v1.2/Work
cd CSPP_Sounder_QL_v1.2/tarballs
```

Download the CSPP-Sounder-QL software and test data tarballs (and associated sha1
checksum files) from http://cimss.ssec.wisc.edu/cspp/...

```bash
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.2.tar.gz
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.2.tar.gz.sha1
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.2_SAMPLE_DATA.tar.gz
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.2_SAMPLE_DATA.tar.gz.sha1
```

Using the sha1 files, confirm the integrity of the downloaded tarballs...

```bash
sha1sum -c CSPP_Sounder_QL_v1.2.tar.gz.sha1
sha1sum -c CSPP_Sounder_QL_v1.2_SAMPLE_DATA.tar.gz.sha1
```

Unpack the tarballs, and change back to the top-level directory CSPP_Sounder_QL_v1.0

```bash
tar xzvf CSPP_Sounder_QL_v1.2.tar.gz -C ../
tar xzvf CSPP_Sounder_QL_v1.2_SAMPLE_DATA.tar.gz -C ../Work
cd ..
```

Set up the run environment...

```bash
export CSPP_SOUNDER_QL_HOME=$(readlink -f Sounder_QL_1_2)
. $CSPP_SOUNDER_QL_HOME/cspp_sounder_ql_env.sh
```

The CSPP-Sounder-QL package is now ready to use, and test data files are included in the
`Sounder_QL_1_2/input` and `Sounder_QL_1_0/truth_output` directories.


In order to test the installation, execute the main driver scripts now with no arguments:

```bash
ql_level2_image.sh
ql_level2_skewt.sh
```

If the installation has been successful to this point, you will be presented with command line
switches and options for the CSPP-Sounder-QL pacakge scripts. This does not mean that your
installation is *processing* data correctly, but it does mean that the
environment has been successfully set to start using CSPP-Sounder-QL.


### Implementation Notes

When installed as shown in section 2.2, the directory/file structure of CSPP-Sounder-QL
(pruned for the purposes of this document) is:

```
Sounder_QL_1_2
    ├── bin
    │   ├── ql_level2_image.sh
    │   └── ql_level2_skewt.sh
    ├── cspp_sounder_ql_env.sh
    ├── cspp_sounder_ql_runtime.sh
    ├── scripts
    │   ├── sounder_image.py
    │   ├── sounder_skewt.py
    │   └── thermo.py
    └── vendor
        └── env
```

## Using CSPP-Sounder-QL

### ql_level2_image.sh

The bash script `ql_level2_image.sh` checks for environment variable `CSPP_SOUNDER_QL_HOME`
and then invokes the Python script `$CSPP_SOUNDER_QL_HOME/scripts/sounder_image.py` ,
which contains the logic to ingest the various level-2 product type, and plot either the navigated
temperature, dewpoint or relative humidity. The script `ql_level2_image.sh` requires, at a
minimum, the name of the input level-2 file, and the type of the input sounder file, in that order...

The bash script `ql_level2_image.sh` checks for environment variable `CSPP_SOUNDER_QL_HOME`
and then invokes the Python script `$CSPP_SOUNDER_QL_HOME/scripts/sounder_image.py`, which
contains the logic to ingest the various level-2 product type, and plot either the navigated
temperature, dewpoint or relative humidity. The script `ql_level2_image.sh` requires, at a
minimum, the type of the input sounder file, and the name of the input level-2 file(s), in that order...

```text
ql_level2_image.sh FILE_TYPE [INPUT_FILE [INPUT_FILE ...]]
```

Various command line options are available for ql_level2_image.sh as shown below:

```text
> ql_level2_image.sh
usage: sounder_image.py [--dset {temp,wvap,dwpt,relh,ctp,ctt}]
                        [--pressure PRESSURE] [-o OUTPUT_FILE] [-V] [-v] [-q]
                        [-h] [-x]
                        [{iapp,mirs,hsrtv,nucaps}]
                        [input_file [input_file ...]]

Create a plot of temperature, dewpoint or something else at a particular pressure level. Supports IAPP, MIRS, HSRTV and NUCAPS files.

positional arguments:
  {iapp,mirs,hsrtv,nucaps}
                        The type of the input sounder data file. Possible values are:
                                'iapp'    (International TOVS Processing Package),
                                'mirs'    (Microwave Integrated Retrieval System),
                                'hsrtv'   (CrIS, AIRS and IASI Hyperspectral Retrieval Software),
                                'nucaps'  (NOAA Unique Combined Atmospheric Processing System),
  input_file            One or more input files.

optional arguments:
  --dset {temp,wvap,dwpt,relh,ctp,ctt}
                        The sounder dataset to plot. Possible values are... 'temp', 'wvap', 'dwpt', 'relh', 'ctp', 'ctt'. [default: temp]
  --pressure PRESSURE   The pressure level (in mbar) to plot. [default: 850.0]
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        The filename of the output Skew-T png file. [default: None]
  -V, --version         Print the CSPP Sounder QL package version
  -v, --verbosity       Each occurrence increases verbosity one level from INFO. -v=DEBUG
  -q, --quiet           Silence all console output
  -h, --help            Show this help message and exit
  -x, --expert          Display all help options, including the expert ones.
```

By default, if no optional parameters are given, the script will plot the temperature dataset at the
pressure level closest to 850mb.

### ql_level2_skewt.sh

The bash script `ql_level2_skewt.sh` checks for environment variable `CSPP_SOUNDER_QL_HOME`
and then invokes the Python script `$CSPP_SOUNDER_QL_HOME/scripts/sounder_skewt.py`,
which contains the logic to ingest the various level-2 product types, and plot the Skew-T plot.
The script `ql_level2_skewt.sh` requires, at a
minimum, the type of the input sounder file, and the name of the input level-2 file(s), in that order...

```text
ql_level2_skewt.sh FILE FILE_TYPE
```

Various command line options are available for `ql_level2_image.sh` as shown below:

```text
> ql_level2_skewt.sh
usage: sounder_skewt.py [-o OUTPUT_FILE] [-V] [-v] [-q] [-h] [-x]
                        [{iapp,mirs,hsrtv,nucaps}] input_file

Create a Skew-T plot from input sounder data. Supports IAPP, MIRS, HSRTV and NUCAPS files.

positional arguments:
  {iapp,mirs,hsrtv,nucaps}
                        The type of the input sounder data file. Possible values are:
                                'iapp'    (International TOVS Processing Package),
                                'mirs'    (Microwave Integrated Retrieval System),
                                'hsrtv'   (CrIS, AIRS and IASI Hyperspectral Retrieval Software),
                                'nucaps'  (NOAA Unique Combined Atmospheric Processing System),
  input_file            A single input file.

optional arguments:
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        The filename of the output Skew-T png file. [default: None]
  -V, --version         Print the CSPP Sounder QL package version
  -v, --verbosity       Each occurrence increases verbosity one level from INFO. -v=DEBUG
  -q, --quiet           Silence all console output
  -h, --help            Show this help message and exit
  -x, --expert          Display all help options, including the expert ones.
```

### Example Plots

Below we have example plots of from each of the packages HSRTV, IAPP, MIRS and NUCAPS. For each package there are the dewpoint temperature,
relative humidity, skew-T, temperature and water vapor mixing ratio. The Skew-T is calculated at the center of the swaths (latitude and longitude is indicated).

```bash
export CSPP_SOUNDER_QL_HOME=$(readlink -f Sounder_QL_1_2)
. $CSPP_SOUNDER_QL_HOME/cspp_sounder_ql_env.sh
cd Work
mkdir hsrtv iapp mirs nucaps
```

#### HSRTV
```bash
for files in $(ls sample_data/input/hsrtv/*.h5); \
do \
    out_prefix=hsrtv/$(basename $files)
    ql_level2_image.sh hsrtv $files --dset temp --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh hsrtv $files --dset dwpt --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh hsrtv $files --dset wvap --plot-min=0.   --plot-max=14.  -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh hsrtv $files --dset relh --plot-min=0.   --plot-max=100. -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh hsrtv $files --dset ctp -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh hsrtv $files --dset ctt -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_skewt.sh hsrtv $files -O $out_prefix ; \
done
```

#### IAPP
```bash
for files in $(ls sample_data/input/iapp/*_iapp.nc); \
do \
    out_prefix=iapp/$(basename $files)
    ql_level2_image.sh iapp $files --dset temp --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 60 -O $out_prefix ; \
    ql_level2_image.sh iapp $files --dset dwpt --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 60 -O $out_prefix ; \
    ql_level2_image.sh iapp $files --dset wvap --plot-min=0.   --plot-max=14.  -m 'l' --scatter-plot -P 60 -O $out_prefix ; \
    ql_level2_image.sh iapp $files --dset relh --plot-min=0.   --plot-max=100. -m 'l' --scatter-plot -P 60 -O $out_prefix ; \
    ql_level2_image.sh iapp $files --dset ctp -m 'l' --scatter-plot -P 60 -O $out_prefix ; \
    ql_level2_image.sh iapp $files --dset ctt -m 'l' --scatter-plot -P 60 -O $out_prefix ; \
    ql_level2_skewt.sh iapp $files -O $out_prefix ; \
done
```

#### MIRS
```bash
for files in $(ls sample_data/input/mirs/NPR-MIRS-SND*.nc); \
do \
    out_prefix=mirs/$(basename $files)
    ql_level2_image.sh mirs $files --dset temp --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh mirs $files --dset dwpt --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh mirs $files --dset wvap --plot-min=0.   --plot-max=14.  -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_image.sh mirs $files --dset relh --plot-min=0.   --plot-max=100. -m 'l' --scatter-plot -P 5 -O $out_prefix ; \
    ql_level2_skewt.sh mirs $files -O $out_prefix ; \
done

for files in $(ls sample_data/input/mirs/SND*.nc); \
do \
    out_prefix=mirs/$(basename $files)
    ql_level2_image.sh mirs $files --dset temp --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_image.sh mirs $files --dset dwpt --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_image.sh mirs $files --dset wvap --plot-min=0.   --plot-max=14.  -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_image.sh mirs $files --dset relh --plot-min=0.   --plot-max=100. -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_skewt.sh mirs $files -O $out_prefix ; \
done
```

#### NUCAPS
```bash
for files in $(ls sample_data/input/nucaps/IASI*.nc); \
do \
    out_prefix=nucaps/$(basename $files)
    ql_level2_image.sh nucaps $files --dset temp --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_image.sh nucaps $files --dset dwpt --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_image.sh nucaps $files --dset wvap --plot-min=0.   --plot-max=14.  -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_image.sh nucaps $files --dset relh --plot-min=0.   --plot-max=100. -m 'l' --scatter-plot -P 25 -O $out_prefix ; \
    ql_level2_skewt.sh nucaps $files -O $out_prefix ; \
done

files=sample_data/input/nucaps/NUCAPS-EDR_v1r0_npp*.nc
out_prefix=nucaps/$(basename $(echo $files | gawk '{print $1}'))
ql_level2_image.sh nucaps $files --dset temp --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 25 -O $out_prefix
ql_level2_image.sh nucaps $files --dset dwpt --plot-min=250. --plot-max=300. -m 'l' --scatter-plot -P 25 -O $out_prefix
ql_level2_image.sh nucaps $files --dset wvap --plot-min=0.   --plot-max=14.  -m 'l' --scatter-plot -P 25 -O $out_prefix
ql_level2_image.sh nucaps $files --dset relh --plot-min=0.   --plot-max=100. -m 'l' --scatter-plot -P 25 -O $out_prefix
ql_level2_skewt.sh nucaps $(echo $files | gawk '{print $1}') -O $out_prefix
```
