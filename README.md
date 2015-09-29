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

Included with the CSPP-Sounder-QL package is a self-contained Python distribution, `ShellB3`.

### Installation of CSPP-Sounder-QL Software

Create the installation directory CSPP_Sounder_QL_v1.0, and the working and tarball directories,
and change into the tarball directory...

```bash
mkdir -p CSPP_Sounder_QL_v1.0/tarballs CSPP_Sounder_QL_v1.0/Work
cd CSPP_Sounder_QL_v1.0/tarballs
```

Download the CSPP-Sounder-QL software and test data tarballs (and associated sha1
checksum files) from http://cimss.ssec.wisc.edu/cspp/...

```bash
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.0.tar.gz
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.0.tar.gz.sha1
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.0_SAMPLE_DATA.tar.gz
wget -c http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/downloads/CSPP_Sounder_QL_v1.0_SAMPLE_DATA.tar.gz.sha1
```

Using the sha1 files, confirm the integrity of the downloaded tarballs...

```bash
sha1sum -c CSPP_Sounder_QL_v1.0.tar.gz.sha1
sha1sum -c CSPP_Sounder_QL_v1.0_SAMPLE_DATA.tar.gz.sha1
```

Unpack the tarballs, and change back to the top-level directory CSPP_Sounder_QL_v1.0

```bash
tar xzvf CSPP_Sounder_QL_v1.0.tar.gz -C ../
tar xzvf CSPP_Sounder_QL_v1.0_SAMPLE_DATA.tar.gz -C ../Work
cd ..
```

Set up the run environment...

```bash
export CSPP_SOUNDER_QL_HOME=$(readlink -f Sounder_QL_1_0)
. $CSPP_SOUNDER_QL_HOME/cspp_sounder_ql_env.sh
```

The CSPP-Sounder-QL package is now ready to use, and test data files are included in the
`Sounder_QL_1_0/input` and `Sounder_QL_1_0/truth_output` directories.


In order to test the installation, execute the main driver scripts now with no arguments:

```bash
$CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh
$CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh
```

If the installation has been successful to this point, you will be presented with command line
switches and options for the CSPP-Sounder-QL pacakge scripts. This does not mean that your
installation is *processing* data correctly, but it does mean that the
environment has been successfully set to start using CSPP-Sounder-QL.


### Implementation Notes

When installed as shown in section 2.2, the directory/file structure of CSPP-Sounder-QL
(pruned for the purposes of this document) is:

```
Sounder_QL_1_0
            │
            ├── cspp_sounder_ql_env.sh
            ├── cspp_sounder_ql_runtime.sh
            ├── cspp_sounder_ql_unset.sh
            │
            ├── common
            │   ├── cspp_cfg
            │   └── ShellB3
            │
            └── scripts
                ├── ql_level2_image.sh
                ├── ql_level2_skewt.sh
                ├── sounder_image.py
                ├── sounder_skewt.py
                ├── skewt.py
                └── thermo.py
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
minimum, the name of the input level-2 file, and the type of the input sounder file, in that order...

```bash
bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh FILE FILE_TYPE`
```

Various command line options are available for ql_level2_image.sh as shown below:

```
bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh

usage: sounder_image.py [-h] [--dset {temp,wvap,dwpt,relh}] [-S STRIDE]
                        [--pressure PRESSURE] [--plotMin PLOTMIN]
                        [--plotMax PLOTMAX] [-d DPI] [--lat_0 LAT_0]
                        [--lon_0 LON_0] [--latMin LATMIN] [--latMax LATMAX]
                        [--lonMin LONMIN] [--lonMax LONMAX] [--scatter_plot]
                        [-P POINTSIZE] [-s SCALE] [-m {c,l,i}]
                        [-o OUTPUT_FILE] [-O OUTPUTFILEPREFIX] [-v] [-q]
                        input_file {IAPP,MIRS,HSRTV,NUCAPS}

Create a contour plot of temperature, dewpoint or something else at a
particular pressure level. Supports IAPP, MIRS, HSRTV and NUCAPS files.

positional arguments:
  input_file            The fully qualified path to the input file(s). May be
                        a file glob (which should be in quotes).
  {IAPP,MIRS,HSRTV,NUCAPS}
                        The type of the input sounder data file. Possible
                        values are... 'IAPP', 'MIRS', 'HSRTV', 'NUCAPS'.

optional arguments:
  -h, --help            show this help message and exit
  --dset {temp,wvap,dwpt,relh}
                        The sounder dataset to plot. Possible values are...
                        'temp', 'wvap', 'dwpt', 'relh', 'ctp','ctt'. [default: temp]
  -S STRIDE, --stride STRIDE
                        Sample every STRIDE pixels in the band data. [default:
                        1]
  --pressure PRESSURE   The pressure level (in mbar) to plot.. [default:
                        850.0]
  --plotMin PLOTMIN     Minimum value to plot.
  --plotMax PLOTMAX     Maximum value to plot.
  -d DPI, --dpi DPI     The resolution in dots per inch of the output png
                        file. [default: 200]
  --lat_0 LAT_0         Center latitude of plot.
  --lon_0 LON_0         Center longitude of plot.
  --latMin LATMIN       Minimum latitude to plot.
  --latMax LATMAX       Maximum latitude to plot.
  --lonMin LONMIN       Minimum longitude to plot.
  --lonMax LONMAX       Maximum longitude to plot.
  --scatter_plot        Generate the plot using a scatterplot approach.
  -P POINTSIZE, --pointSize POINTSIZE
                        Size of the plot point used to represent each pixel.
                        [default: 1]
  -s SCALE, --scale SCALE
                        The scaling factor for the default viewport size of (w
                        x h) = (4200000.0 x 4500000.0) meters. [default: 1]
  -m {c,l,i}, --map_res {c,l,i}
                        The map coastline resolution. Possible values are 'c'
                        (coarse),'l' (low) and 'i' (intermediate). [default:
                        'c']
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        The filename of the output png file.
  -O OUTPUTFILEPREFIX, --output_file_prefix OUTPUTFILEPREFIX
                        String to prefix to the automatically generated png
                        name. [default: None]
  -v, --verbose         each occurrence increases verbosity 1 level from INFO.
                        -v=DEBUG
  -q, --quiet           Silence all console output
```

By default, if no optional parameters are given, the script will plot the temperature dataset at the
pressure level closest to 850mb.

### ql_level2_skewt.sh

The bash script `ql_level2_skewt.sh` checks for environment variable `CSPP_SOUNDER_QL_HOME`
and then invokes the Python script `$CSPP_SOUNDER_QL_HOME/scripts/sounder_skewt.py`,
which contains the logic to ingest the various level-2 product type, and plot either the navigated
temperature, dewpoint or relative humidity. The script `ql_level2_skewt.sh` requires, at a
minimum, the name of the input level-2 file, and the type of the input sounder file, in that order...

```bash
bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh FILE FILE_TYPE
```

Various command line options are available for `ql_level2_image.sh` as shown below:

```
bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh

usage: sounder_skewt.py [-h] [--temp_min TEMP_MIN] [--temp_max TEMP_MAX]
                        [-d DPI] [--lat_0 LAT_0] [--lon_0 LON_0]
                        [-o OUTPUT_FILE] [-O OUTPUTFILEPREFIX] [-v] [-q]
                        input_file {IAPP,MIRS,HSRTV,NUCAPS}

Create a Skew-T plot from input sounder data.

positional arguments:
  input_file            The fully qualified path to a single level-1d input
                        file.
  {IAPP,MIRS,HSRTV,NUCAPS}
                        The type of the input sounder data file. Possible
                        values are... 'IAPP', 'MIRS', 'HSRTV', 'NUCAPS'.

optional arguments:
  -h, --help            show this help message and exit
  --temp_min TEMP_MIN   Minimum temperature value to plot. [default: -40 deg
                        C]
  --temp_max TEMP_MAX   Maximum temperature value to plot. [default: 40 deg C]
  -d DPI, --dpi DPI     The resolution in dots per inch of the output png
                        file. [default: 200]
  --lat_0 LAT_0         Plot the sounder footprint closest to lat_0.
  --lon_0 LON_0         Plot the sounder footprint closest to lon_0.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        The filename of the output Skew-T png file. [default:
                        None]
  -O OUTPUTFILEPREFIX, --output_file_prefix OUTPUTFILEPREFIX
                        String to prefix to the automatically generated png
                        name. [default: None]
  -v, --verbose         each occurrence increases verbosity 1 level from INFO.
                        -v=DEBUG
  -q, --quiet           Silence all console output
```

### Example Plots

Below we have example plots of from each of the packages HSRTV, IAPP, MIRS and NUCAPS, broken down by satellite. Each row has the dewpoint temperature,
relative humidity, skew-T, temperature and water vapor mixing ratio. The Skew-T is calculated at the center of the swaths (latitude and longitude is indicated).

#### HSRTV

- Suomi-NPP

    ```bash
    for files in $(ls ../../input/hsrtv/CrIS*.h5); \
    do \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset temp --plotMin=250. --plotMax=300. -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset dwpt --plotMin=250. --plotMax=300. -m 'l' --scatter_plot -P 8; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset wvap --plotMin=0.   --plotMax=14.  -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset relh --plotMin=0.   --plotMax=100. -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset ctp -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset ctt -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh $files HSRTV ; \
    done 
    ```

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/CrIS_d20150403_t190649.atm_prof_rtv.h5.HSRTV_dwpt_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/CrIS_d20150403_t190649.atm_prof_rtv.h5.HSRTV_relh_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/CrIS_d20150403_t190649.atm_prof_rtv.h5.HSRTV_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/CrIS_d20150403_t190649.atm_prof_rtv.h5.HSRTV_temp_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/CrIS_d20150403_t190649.atm_prof_rtv.h5.HSRTV_wvap_852mb.png "" "width:200px")


- Metop-B

    ```bash
    for files in $(ls ../../input/hsrtv/IASI*.h5); \
    do \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset temp --plotMin=250. --plotMax=300. -m 'l' ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset dwpt --plotMin=250. --plotMax=300. -m 'l' ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset wvap --plotMin=0.   --plotMax=14.  -m 'l' ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset relh --plotMin=0.   --plotMax=100. -m 'l' ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset ctp -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files HSRTV --dset ctt -m 'l' --scatter_plot -P 8 ; \
        bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh $files HSRTV ; \
    done 
    ```

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/IASI_d20150331_t153700_M01.atm_prof_rtv.h5.HSRTV_dwpt_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/IASI_d20150331_t153700_M01.atm_prof_rtv.h5.HSRTV_relh_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/IASI_d20150331_t153700_M01.atm_prof_rtv.h5.HSRTV_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/IASI_d20150331_t153700_M01.atm_prof_rtv.h5.HSRTV_temp_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/hsrtv/IASI_d20150331_t153700_M01.atm_prof_rtv.h5.HSRTV_wvap_852mb.png "" "width:200px")

#### IAPP

```bash
for files in $(ls ../../input/iapp/*_iapp.nc); \
do \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files IAPP --dset temp --plotMin=250. --plotMax=300. -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files IAPP --dset dwpt --plotMin=250. --plotMax=300. -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files IAPP --dset wvap --plotMin=0.   --plotMax=14.  -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files IAPP --dset relh --plotMin=0.   --plotMax=100. -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files IAPP --dset ctp -m 'l' --scatter_plot -P 8 ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files IAPP --dset ctt -m 'l' --scatter_plot -P 8 ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh $files IAPP ; \
done 
```

- Metop-A

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopa_L2_d20150402_t1549110_e1554374_c20150402160354864396_iapp.nc.IAPP_dwpt_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopa_L2_d20150402_t1549110_e1554374_c20150402160354864396_iapp.nc.IAPP_relh_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopa_L2_d20150402_t1549110_e1554374_c20150402160354864396_iapp.nc.IAPP_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopa_L2_d20150402_t1549110_e1554374_c20150402160354864396_iapp.nc.IAPP_temp_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopa_L2_d20150402_t1549110_e1554374_c20150402160354864396_iapp.nc.IAPP_wvap_850mb.png "" "width:200px")

- Metop-B

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopb_L2_d20150330_t0341557_e0353397_c20150330040308043882_iapp.nc.IAPP_dwpt_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopb_L2_d20150330_t0341557_e0353397_c20150330040308043882_iapp.nc.IAPP_relh_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopb_L2_d20150330_t0341557_e0353397_c20150330040308043882_iapp.nc.IAPP_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopb_L2_d20150330_t0341557_e0353397_c20150330040308043882_iapp.nc.IAPP_temp_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/metopb_L2_d20150330_t0341557_e0353397_c20150330040308043882_iapp.nc.IAPP_wvap_850mb.png "" "width:200px")

- NOAA-15

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa15_L2_d20060727_t2100044_e2159548_c20150402223923492393_iapp.nc.IAPP_dwpt_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa15_L2_d20060727_t2100044_e2159548_c20150402223923492393_iapp.nc.IAPP_relh_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa15_L2_d20060727_t2100044_e2159548_c20150402223923492393_iapp.nc.IAPP_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa15_L2_d20060727_t2100044_e2159548_c20150402223923492393_iapp.nc.IAPP_temp_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa15_L2_d20060727_t2100044_e2159548_c20150402223923492393_iapp.nc.IAPP_wvap_850mb.png "" "width:200px")


- NOAA-16

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa16_L2_d20140105_t0519406_e0713540_c20150402224320161020_iapp.nc.IAPP_dwpt_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa16_L2_d20140105_t0519406_e0713540_c20150402224320161020_iapp.nc.IAPP_relh_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa16_L2_d20140105_t0519406_e0713540_c20150402224320161020_iapp.nc.IAPP_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa16_L2_d20140105_t0519406_e0713540_c20150402224320161020_iapp.nc.IAPP_temp_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa16_L2_d20140105_t0519406_e0713540_c20150402224320161020_iapp.nc.IAPP_wvap_850mb.png "" "width:200px")


- NOAA-18

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa18_L2_d20150331_t1023110_e1034432_c20150331104437295133_iapp.nc.IAPP_dwpt_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa18_L2_d20150331_t1023110_e1034432_c20150331104437295133_iapp.nc.IAPP_relh_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa18_L2_d20150331_t1023110_e1034432_c20150331104437295133_iapp.nc.IAPP_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa18_L2_d20150331_t1023110_e1034432_c20150331104437295133_iapp.nc.IAPP_temp_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa18_L2_d20150331_t1023110_e1034432_c20150331104437295133_iapp.nc.IAPP_wvap_850mb.png "" "width:200px")


- NOAA-19

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa19_L2_d20150331_t0905533_e0918221_c20150331092749438101_iapp.nc.IAPP_dwpt_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa19_L2_d20150331_t0905533_e0918221_c20150331092749438101_iapp.nc.IAPP_relh_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa19_L2_d20150331_t0905533_e0918221_c20150331092749438101_iapp.nc.IAPP_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa19_L2_d20150331_t0905533_e0918221_c20150331092749438101_iapp.nc.IAPP_temp_850mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/iapp/noaa19_L2_d20150331_t0905533_e0918221_c20150331092749438101_iapp.nc.IAPP_wvap_850mb.png "" "width:200px")


#### MIRS

```bash
for files in $(ls ../../input/mirs/SND*.nc); \
do \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files MIRS --dset temp --plotMin=250. --plotMax=300. -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files MIRS --dset dwpt --plotMin=250. --plotMax=300. -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files MIRS --dset wvap --plotMin=0.   --plotMax=14.  -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh $files MIRS --dset relh --plotMin=0.   --plotMax=100. -m 'l' ; \
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh $files MIRS ; \
done 
```

- Metop-A

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M2.D15091.S1610.E1615.B0000001.WE.LR.ORB.nc.MIRS_dwpt_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M2.D15091.S1610.E1615.B0000001.WE.LR.ORB.nc.MIRS_relh_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M2.D15091.S1610.E1615.B0000001.WE.LR.ORB.nc.MIRS_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M2.D15091.S1610.E1615.B0000001.WE.LR.ORB.nc.MIRS_temp_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M2.D15091.S1610.E1615.B0000001.WE.LR.ORB.nc.MIRS_wvap_840mb.png "" "width:200px")

- Metop-B

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M1.D15091.S1516.E1528.B0000001.WE.LR.ORB.nc.MIRS_dwpt_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M1.D15091.S1516.E1528.B0000001.WE.LR.ORB.nc.MIRS_relh_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M1.D15091.S1516.E1528.B0000001.WE.LR.ORB.nc.MIRS_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M1.D15091.S1516.E1528.B0000001.WE.LR.ORB.nc.MIRS_temp_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.M1.D15091.S1516.E1528.B0000001.WE.LR.ORB.nc.MIRS_wvap_840mb.png "" "width:200px")

- NOAA-18

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NN.D15091.S2139.E2151.B5083131.WE.LR.ORB.nc.MIRS_dwpt_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NN.D15091.S2139.E2151.B5083131.WE.LR.ORB.nc.MIRS_relh_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NN.D15091.S2139.E2151.B5083131.WE.LR.ORB.nc.MIRS_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NN.D15091.S2139.E2151.B5083131.WE.LR.ORB.nc.MIRS_temp_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NN.D15091.S2139.E2151.B5083131.WE.LR.ORB.nc.MIRS_wvap_840mb.png "" "width:200px")

- NOAA-19

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NP.D15089.S1904.E1916.B3164343.WE.LR.ORB.nc.MIRS_dwpt_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NP.D15089.S1904.E1916.B3164343.WE.LR.ORB.nc.MIRS_relh_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NP.D15089.S1904.E1916.B3164343.WE.LR.ORB.nc.MIRS_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NP.D15089.S1904.E1916.B3164343.WE.LR.ORB.nc.MIRS_temp_840mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/mirs/SND_SX.NP.D15089.S1904.E1916.B3164343.WE.LR.ORB.nc.MIRS_wvap_840mb.png "" "width:200px")


#### NUCAPS

- Suomi-NPP

    ```bash
    files="../../input/nucaps/NUCAPS-EDR_v1r0_npp_s2015040319*.nc"
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh "$files" NUCAPS --dset temp --plotMin=250. --plotMax=300. -m 'l'
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh "$files" NUCAPS --dset dwpt --plotMin=250. --plotMax=300. -m 'l'
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh "$files" NUCAPS --dset wvap --plotMin=0.   --plotMax=14.  -m 'l'
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_image.sh "$files" NUCAPS --dset relh --plotMin=0.   --plotMax=100. -m 'l'
    bash $CSPP_SOUNDER_QL_HOME/scripts/ql_level2_skewt.sh "$files" NUCAPS
    ```

![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/nucaps/NUCAPS-EDR_v1r0_npp_s201504031906499_e201504031907197_c201504032037040.nc.NUCAPS_dwpt_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/nucaps/NUCAPS-EDR_v1r0_npp_s201504031906499_e201504031907197_c201504032037040.nc.NUCAPS_relh_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/nucaps/NUCAPS-EDR_v1r0_npp_s201504031906499_e201504031907197_c201504032037040.nc.NUCAPS_SkewT.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/nucaps/NUCAPS-EDR_v1r0_npp_s201504031906499_e201504031907197_c201504032037040.nc.NUCAPS_temp_852mb.png "" "width:200px")![Alt Text](http://www.ssec.wisc.edu/~geoffc/CSPP_Sounder_QL/quicklooks/nucaps/NUCAPS-EDR_v1r0_npp_s201504031906499_e201504031907197_c201504032037040.nc.NUCAPS_wvap_852mb.png "" "width:200px")

