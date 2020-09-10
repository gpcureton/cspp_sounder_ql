#!/usr/bin/env python
# encoding: utf-8
"""
args.py

 * DESCRIPTION: This file uses the argparse module to gather up the command line inputs, and
 feeding them to the main loop.

Created by Geoff Cureton on 2018-08-10.
Copyright (c) 2018 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
import sys
import argparse
from datetime import datetime
from collections import OrderedDict
import logging

import log_common
from utils import check_and_convert_path, create_dir

LOG = logging.getLogger(__name__)

def argument_parser_common():
    '''
    Common options for command line applications.
    '''

    # Do we need the expert help messages...
    is_expert = False
    if '--expert' in sys.argv:
        expert_index = sys.argv.index('--expert')
        sys.argv[expert_index] = '--help'
        is_expert = True
    elif '-x' in sys.argv:
        expert_index = sys.argv.index('-x')
        sys.argv[expert_index] = '--help'
        is_expert = True
    else:
        pass
    flags = {}
    flags['is_expert'] = is_expert

    cspp_sounder_ql_version = 'cspp-sounder-ql-2.0'
    
    # Set the common help strings.
    help_strings = {}
    help_strings['work_dir'] = """Work directory which all activity will occur in, defaults to""" \
        """ current directory."""
    help_strings['debug'] = """Always retain intermediate files. [default: %(default)s]"""
    help_strings['verbosity'] = """Each occurrence increases verbosity one level from the""" \
        """ default, where the\nlevels are [ERROR, WARNING, INFO, DEBUG].  [default: INFO]"""
    help_strings['version'] = """Print the CSPP Sounder QL package version"""
    help_strings['help'] = """Show this help message and exit"""
    help_strings['expert'] = """Display all help options, including the expert ones."""

    # Construct the dictionaries containing the parser options.

    parser_lists = OrderedDict()
    parser_dicts = OrderedDict()

    parser_lists['work_dir'] = ['-W', '--work-dir']
    parser_dicts['work_dir'] = {
        'action': "store",
        'dest': 'work_dir',
        'metavar': 'WORK_DIR',
        'default': '.',
        'help': help_strings['work_dir']
    }

    parser_lists['debug'] = ['-d', '--debug']
    parser_dicts['debug'] = {
        'action': 'store_true',
        'dest': 'debug',
        'default': False,
        'help': help_strings['debug']
    }

    parser_lists['verbosity'] = ["-v", "--verbosity"]
    parser_dicts['verbosity'] = {
        'action': "count",
        'dest': 'verbosity',
        'default': 2,
        'help': help_strings['verbosity']
    }

    parser_lists['version'] = ['-V', '--version']
    parser_dicts['version'] = {
        'action': 'version',
        'version': cspp_sounder_ql_version,
        'help': help_strings['version']
    }

    parser_lists['help'] = ["-h", "--help"]
    parser_dicts['help'] = {
        'action': "help",
        'help': help_strings['help']
    }

    parser_lists['expert'] = ['-x', '--expert']
    parser_dicts['expert'] = {
        'action': "store_true",
        'dest': 'is_expert',
        'default': False,
        'help': help_strings['expert']
    }

    return parser_lists, parser_dicts, flags

def argument_parser_image():
    '''
    geocat specific options.
    '''

    # Get the common options...
    common_parser_lists, common_parser_dicts, flags = argument_parser_common()

    dataChoices = OrderedDict([
            ('iapp', 'International TOVS Processing Package'),
            ('mirs', 'Microwave Integrated Retrieval System'),
            ('hsrtv', 'CrIS, AIRS and IASI Hyperspectral Retrieval Software'),
            ('heap', 'Hyper-Spectral Enterprise Algorithm Package'),
            ('nucaps', 'NOAA Unique Combined Atmospheric Processing System')
            ])
    prodChoices=['temp', 'wvap', 'dwpt','relh']
    # prodChoices=['temp','temp_gdas','wvap','wvap_gdas','dwpt','relh','relh_gdas']
                 # 'ctp','ctt']
                 # 'ctp','ctt','2temp','cold_air_aloft']
    map_res_choice = ['c','l','i']
    plot_type_choice = ['image']
    # plot_type_choice = ['image','slice']
    map_proj_choice = {
            'eckiv':   'Eckert IV Equal Area',
            'lcc':     'Lambert Conformal',
            'ortho':   'Orthographic',
            'geos':    'Geostationary',
            'npstere': 'North-Polar Stereographic',
            'spstere': 'South-Polar Stereographic',
            'plat':    'PlateCarree',
            'sinu':    'Sinusoidal',
            'moll':    'Mollweide'
            }

    defaults = {
            'datatype':None,
            'dataset':'temp',
            'stride':1,
            'pressure':800.,
            # 'elevation':None,
            'plotMin'  : None,
            'plotMax'  : None,
            'dpi':200,
            'lat_0':None,
            'lon_0':None,
            # 'yMin'  : None,
            # 'yMax'  : None,
            # 'bounding_lat':0.,
            'image_size' : [7.5, 7.5],
            'map_axis' : [0.05, 0.15, 0.90, 0.75 ],
            'cbar_axis' : [0.10 , 0.05, 0.8 , 0.05],
            'scatter_plot':False,
            'global':False,
            'pointSize':1,
            'font_scale':1.,
            'cmap':None,
            'plot_type'  : 'image',
            # 'scale':1,
            'map_res':'c',
            'proj':'lcc',
            'outputFilePrefix' : None,
            'llcrnrx':None,
            'llcrnry':None,
            'urcrnrx':None,
            'urcrnry':None,
            'output_file':None,
            'unnavigated':False
            }

    help_strings = {}
    help_strings['inputs'] = """One or more input files or directories."""
    help_strings['work_dir'] = """The work directory."""
    help_strings['datatype'] = """The type of the input sounder data file. Possible values """ \
            """are:\n {}""".format(
            "".join(["\t'{0:8s} ({1}),\n".format(tups[0]+"'",tups[1]) for
                tups in zip(dataChoices.keys(),dataChoices.values())]))
    help_strings['dataset'] = """The sounder dataset to plot. Possible values are...\n""" \
            """{},\n{}. [default: {}]""".format(prodChoices[:6].__str__()[1:-1],
                    prodChoices[6:].__str__()[1:-1], defaults["dataset"])
    help_strings['stride'] = """Sample every STRIDE pixels in the band data. """ \
            """[default: {}]""".format(defaults["stride"])
    help_strings['pressure'] = """The pressure level (in mbar) to plot [default: {}]. """ \
            """""".format(defaults["pressure"])
            # """Mutually exclusive with --elevation.""".format(defaults["pressure"])
    # help_strings['elevation'] = """The elevation level (in feet) to plot [default: {}]. """ \
            # """Mutually exclusive with --pressure.""".format(defaults["elevation"])
    help_strings['plotMin'] = "Minimum value to plot.".format(defaults["plotMin"])
    help_strings['plotMax'] = "Maximum value to plot.".format(defaults["plotMax"])
    help_strings['dpi'] = """The resolution in dots per inch of the output png file. """ \
            """[default: {}]""".format(defaults["dpi"])
    help_strings['lat_0'] = "Center latitude of plot."
    help_strings['lon_0'] = "Center longitude of plot."
    # help_strings['latMin'] = "Minimum latitude to plot."
    # help_strings['latMax'] = "Maximum latitude to plot."
    # help_strings['lonMin'] = "Minimum longitude to plot."
    # help_strings['lonMax'] = "Maximum longitude to plot."
    # help_strings['yMin'] = "Minimum y-axis limit on slice plots."
    # help_strings['yMax'] = "Maximum y-axis limit on slice plots."
    # help_strings['footprint'] = """The cross-track footprint to use when plotting a profile slice."""
    # help_strings['bounding_lat'] = """The minimum/maximum latitude to plot for the polar """ \
            # """projections. [default: {}]""".format(defaults["bounding_lat"])
    # help_strings['viewport'] = """Lower-left and upper-right coordinates """ \
            # """[*llcrnrx*, *llcrnry*, *urcrnrx*, *urcrnry*]\nof the projection viewport, in the """ \
            # """range [-0.5,+0.5] (for navigated plots only)"""
    help_strings['image_size'] = """The size of the output image [*width*, *height*] in inches. """ \
            """[default: '{}']""".format(defaults["image_size"])
    help_strings['map_axis'] = """Set the map axes at position """ \
            """[*left*, *bottom*, *width*, *height*] where all\nquantities are in fractions of """ \
            """figure width and height. [default: '{}']""".format(defaults["map_axis"])
    help_strings['cbar_axis'] = """Set the colorbar axes at position """ \
            """[*left*, *bottom*, *width*, *height*] where all\nquantities are in fractions of """ \
            """figure width and height. [default: '{}']""".format(defaults["cbar_axis"])
    help_strings['scatter_plot'] = "Generate the plot using a scatterplot approach."
    help_strings['global'] = "Plot the maximum extent of the desired projection, rather than just """ \
            """the extent\nof the data. [default: '{}']""".format(defaults["global"])
    help_strings['logscale'] = """Plot the dataset using a logarithmic scale."""
    help_strings['no_logscale'] = """Plot the dataset using a linear scale."""
    help_strings['pointSize'] = """Size of the plot point used to represent each pixel. """ \
            """[default: {}]""".format(defaults["pointSize"])
    help_strings['font_scale'] = """The scale factor to apply to the default font size for the""" \
            """plot labels. [default: {}]""".format(defaults["font_scale"])
    help_strings['cmap'] = """The matplotlib colormap to use. [default: '{}']""".format(defaults["cmap"])
    help_strings['plot_type'] = """The type of plot. Possible values are {}.""".format(
            plot_type_choice.__str__()[1:-1])
    help_strings['plot_title'] = """The plot title. Must be placed in double quotes."""
    help_strings['cbar_title'] = """The colourbar title. Must be placed in double quotes."""
    # help_strings['scale'] = """The scaling factor for the default viewport size of (w x h) = """ \
            # """(4200000.0 x 4500000.0)\nmeters. [default: {}]""".format(defaults["scale"])
    help_strings['map_res'] = """The map coastline resolution. Possible values are 'c' (coarse), """ \
            """'l' (low) and\n'i' (intermediate). [default: '{}']""".format(defaults["map_res"])
    help_strings['proj'] = """The map projection. Possible values are:\n {}[default: '{}']""".format(
            "".join(["\t'{0:8s} ({1}),\n".format(tups[0]+"'",tups[1]) for
                tups in zip(map_proj_choice.keys(),map_proj_choice.values())]),
            defaults["proj"])

    help_strings['output_file'] = """The filename of the output png file.
                      """
    help_strings['outputFilePrefix'] = """String to prefix to the automatically generated png """ \
            """name. [default: {}]""".format(defaults["outputFilePrefix"])

    # Construct the dictionaries containing the parser options.

    parser_lists = OrderedDict()
    parser_dicts = OrderedDict()

    is_expert = flags['is_expert']

    # Mandatory/positional arguments

    parser_lists['datatype'] = []
    parser_dicts['datatype'] = {
        'dest': "datatype",
        'action': "store",
        'choices': dataChoices.keys(),
        'nargs': '?',
        'help': help_strings['datatype']
    }

    parser_lists['inputs'] = []
    parser_dicts['inputs'] = {
        'action': 'store',
        'dest': 'inputs',
        'type': str,
        'nargs': '*',
        'help': help_strings['inputs']
    }

    # Optional arguments

    parser_lists['dataset'] = ['--dset']
    parser_dicts['dataset'] = {
        'dest': "dataset",
        'action': "store",
        'default': defaults['dataset'],
        'choices': prodChoices,
        'help': help_strings['dataset']
    }

    parser_lists['plot_type'] = ['--plot-type']
    parser_dicts['plot_type'] = {
        'action': "store",
        'dest': "plot_type",
        'default': defaults["plot_type"],
        'type': str,
        'choices': plot_type_choice,
        'help': help_strings['plot_type']
    }

    parser_lists['stride'] = ['-S','--stride']
    parser_dicts['stride'] = {
        'dest': "stride",
        'action': "store",
        'default': defaults["stride"],
        'type': int,
        'help': help_strings['stride'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['pressure'] = ['--pressure']
    parser_dicts['pressure'] = {
        'dest': "pressure",
        'action': "store",
        'default': defaults["pressure"],
        'type': float,
        'help': help_strings['pressure'] if is_expert else argparse.SUPPRESS
    }

    # parser_lists['elevation'] = ['--elevation']
    # parser_dicts['elevation'] = {
        # 'action': "store",
        # 'dest': "elevation",
        # 'default': defaults["elevation"],
        # 'type': float,
        # 'help': help_strings['elevation'] if is_expert else argparse.SUPPRESS
    # }

    parser_lists['plotMin'] = ['--plotMin']
    parser_dicts['plotMin'] = {
        'action': "store",
        'dest': "plotMin",
        'default': defaults["plotMin"],
        'type': float,
        'help': help_strings['plotMin'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['plotMax'] = ['--plotMax']
    parser_dicts['plotMax'] = {
        'action': "store",
        'dest': "plotMax",
        'default': defaults["plotMax"],
        'type': float,
        'help': help_strings['plotMax'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['dpi'] = ['--dpi']
    parser_dicts['dpi'] = {
        'action': "store",
        'dest': "dpi",
        'default': defaults["dpi"],
        'type': float,
        'help': help_strings['dpi'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['lat_0'] = ['--lat-0']
    parser_dicts['lat_0'] = {
        'action': "store",
        'dest': "lat_0",
        'type': float,
        'help': help_strings['lat_0'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['lon_0'] = ['--lon-0']
    parser_dicts['lon_0'] = {
        'action': "store",
        'dest': "lon_0",
        'type': float,
        'help': help_strings['lon_0'] if is_expert else argparse.SUPPRESS
    }

    # parser_lists['latMin'] = ['--latMin']
    # parser_dicts['latMin'] = {
        # 'action': "store",
        # 'dest': "latMin",
        # 'type': float,
        # 'help': help_strings['latMin'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['latMax'] = ['--latMax']
    # parser_dicts['latMax'] = {
        # 'action': "store",
        # 'dest': "latMax",
        # 'type': float,
        # 'help': help_strings['latMax'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['lonMin'] = ['--lonMin']
    # parser_dicts['lonMin'] = {
        # 'action': "store",
        # 'dest': "lonMin",
        # 'type': float,
        # 'help': help_strings['lonMin'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['lonMax'] = ['--lonMax']
    # parser_dicts['lonMax'] = {
        # 'action': "store",
        # 'dest': "lonMax",
        # 'type': float,
        # 'help': help_strings['lonMax'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['yMin'] = ['--yMin']
    # parser_dicts['yMin'] = {
        # 'action': "store",
        # 'dest': "yMin",
        # 'type': float,
        # 'help': help_strings['yMin'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['yMax'] = ['--yMax']
    # parser_dicts['yMax'] = {
        # 'action': "store",
        # 'dest': "yMax",
        # 'type': float,
        # 'help': help_strings['yMax'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['footprint'] = ['--footprint']
    # parser_dicts['footprint'] = {
        # 'action': "store",
        # 'dest': "footprint",
        # 'type': int,
        # 'help': help_strings['footprint'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['bounding_lat'] = ['--bounding-lat']
    # parser_dicts['bounding_lat'] = {
        # 'action': "store",
        # 'dest': "bounding_lat",
        # 'default': defaults["bounding_lat"],
        # 'type': float,
        # 'help': help_strings['bounding_lat'] if is_expert else argparse.SUPPRESS
    # }

    # parser_lists['viewport'] = ['--viewport']
    # parser_dicts['viewport'] = {
        # 'action': "store",
        # 'dest': "viewport",
        # 'type': float,
        # 'nargs': 4,
        # 'metavar': ('LLCRNRX', 'LLCRNRY', 'URCRNRX', 'URCRNRY'),
        # 'help': help_strings['viewport'] if is_expert else argparse.SUPPRESS
    # }

    parser_lists['image_size'] = ['--image-size']
    parser_dicts['image_size'] = {
        'action': "store",
        'dest': "image_size",
        'default': defaults["image_size"],
        'type': float,
        'nargs': 2,
        'metavar': ('WIDTH', 'HEIGHT'),
        'help': help_strings['image_size'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['map_axis'] = ['--map-axis']
    parser_dicts['map_axis'] = {
        'action': "store",
        'dest': "map_axis",
        'default': defaults["map_axis"],
        'type': float,
        'nargs': 4,
        'metavar': ('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
        'help': help_strings['map_axis'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['cbar_axis'] = ['--cbar-axis']
    parser_dicts['cbar_axis'] = {
        'action': "store",
        'dest': "cbar_axis",
        'default': defaults["cbar_axis"],
        'type': float,
        'nargs': 4,
        'metavar': ('LEFT', 'BOTTOM', 'WIDTH', 'HEIGHT'),
        'help': help_strings['cbar_axis'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['scatter_plot'] = ['--scatter-plot']
    parser_dicts['scatter_plot'] = {
        'action': "store_true",
        'dest': "doScatterPlot",
        'default': defaults["scatter_plot"],
        'help': help_strings['scatter_plot']
    }

    parser_lists['global'] = ['--global']
    parser_dicts['global'] = {
        'action': "store_true",
        'dest': "plot_global",
        'default': defaults["global"],
        'help': help_strings['global']
    }

    parser_lists['logscale'] = ['--logscale']
    parser_dicts['logscale'] = {
        'action': "store_true",
        'dest': "logscale",
        'help': help_strings['logscale']
    }

    parser_lists['no_logscale'] = ['--no-logscale']
    parser_dicts['no_logscale'] = {
        'action': "store_true",
        'dest': "no_logscale",
        'help': help_strings['no_logscale'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['pointSize'] = ['-P','--pointSize']
    parser_dicts['pointSize'] = {
        'action': "store",
        'dest': "pointSize",
        'default': defaults["pointSize"],
        'type': float,
        'help': help_strings['pointSize'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['font_scale'] = ['--font-scale']
    parser_dicts['font_scale'] = {
        'action': "store",
        'dest': "font_scale",
        'default': defaults["font_scale"],
        'type': float,
        'help': help_strings['font_scale'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['cmap'] = ['--cmap']
    parser_dicts['cmap'] = {
        'action': "store",
        'dest': "cmap",
        'default': defaults["cmap"],
        'type': str,
        'help': help_strings['cmap'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['plot_title'] = ['--plot-title']
    parser_dicts['plot_title'] = {
        'action': "store",
        'dest': "plot_title",
        'type': str,
        'help': help_strings['plot_title'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['cbar_title'] = ['--cbar-title']
    parser_dicts['cbar_title'] = {
        'action': "store",
        'dest': "cbar_title",
        'type': str,
        'help': help_strings['cbar_title'] if is_expert else argparse.SUPPRESS
    }

    # parser_lists['scale'] = ['-s','--scale']
    # parser_dicts['scale'] = {
        # 'action': "store",
        # 'dest': "scale",
        # 'default': defaults["scale"],
        # 'type': float,
        # 'help': help_strings['scale'] if is_expert else argparse.SUPPRESS
    # }

    parser_lists['map_res'] = ['-m','--map-res']
    parser_dicts['map_res'] = {
        'action': "store",
        'dest': "map_res",
        'default': defaults["map_res"],
        'type': str,
        'choices': map_res_choice,
        'help': help_strings['map_res'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['proj'] = ['--proj']
    parser_dicts['proj'] = {
        'action': "store",
        'dest': "proj",
        'default': defaults["proj"],
        'type': str,
        'choices': map_proj_choice,
        'help': help_strings['proj'] if is_expert else argparse.SUPPRESS
    }

    parser_lists['output_file'] = ['-o','--output-file']
    parser_dicts['output_file'] = {
        'action': "store",
        'dest': "output_file",
        'default': defaults["output_file"],
        'type': str,
        'help': help_strings['output_file']
    }

    parser_lists['outputFilePrefix'] = ['-O','--output-file-prefix']
    parser_dicts['outputFilePrefix'] = {
        'action': "store",
        'dest': "outputFilePrefix",
        'default': defaults["outputFilePrefix"],
        'type': str,
        'help': help_strings['outputFilePrefix']
    }

    # Initialise the parser.
    desc = '''Create a plot of temperature, dewpoint or something else at a particular ''' \
            '''pressure level.\nSupports IAPP, MIRS, HSRTV, HEAP and NUCAPS files.'''
            # '''pressure or elevation level.\nSupports IAPP, MIRS, HSRTV, HEAP and NUCAPS files.'''
    epilog = ''
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False,
                                     epilog=epilog
                                     )

    # Add the quicklook image options...
    image_options = parser_lists.keys()
    for option in image_options:
        parser.add_argument(*parser_lists[option], **parser_dicts[option])

    # Add the common options...
    common_options = common_parser_lists.keys()
    for option in common_options:
        parser.add_argument(*common_parser_lists[option], **common_parser_dicts[option])

    # Parse the options lists
    args = parser.parse_args()

    # Set up the logging
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[args.verbosity if args.verbosity < 4 else 3]

    # Create the work directory if it doesn't exist
    work_dir = os.path.abspath(os.path.expanduser(args.work_dir))
    work_dir = create_dir(work_dir)
    work_dir = check_and_convert_path("WORK_DIR", work_dir)

    dt = datetime.utcnow()
    timestamp = dt.strftime('%Y%m%dT%H%M%S')
    logname = "cspp_sounder_ql." + timestamp + ".log"
    logfile = os.path.join(work_dir, logname)
    log_common.configure_logging(level, FILE=logfile)
    # log_common.configure_logging(level, FILE=None)

    LOG.debug('work directory : {}'.format(work_dir))

    #docleanup = True
    #if args.debug is True:
        #docleanup = False


    return args, work_dir


def argument_parser_skewt():
    '''
    geocat specific options.
    '''

    # Get the common options...
    common_parser_lists, common_parser_dicts, flags = argument_parser_common()

def argument_parser_slice():
    '''
    geocat specific options.
    '''

    # Get the common options...
    common_parser_lists, common_parser_dicts, flags = argument_parser_common()
