#!/usr/bin/env python
# encoding: utf-8
'''
utils.py

 * DESCRIPTION: This is the collection of utilities for parsing text, running external processes and
 binaries, checking environments and whatever other mundane tasks aren't specific to this project.

Created Oct 2011 by R.K.Garcia <rayg@ssec.wisc.edu>
Copyright (c) 2011 University of Wisconsin Regents.
Licensed under GNU GPLv3.
'''

import os
import sys
import re
import string
import logging
import traceback
import time
import fileinput
import shutil
from copy import copy
import uuid
from subprocess import Popen, call, PIPE
from datetime import datetime, timedelta, timezone
from threading import Thread
from queue import Queue, Empty
from six import string_types
from pathlib import Path

LOG = logging.getLogger(__name__)

PROFILING_ENABLED = os.environ.get('CSPP_PROFILE', None) is not None
STRACE_ENABLED = os.environ.get('CSPP_STRACE', None) is not None


def gunzip(local, preserve_original=True):
    """ Uncompress a gzip'd file named "local" (probably ending with ".gz")

    If local ends with ".gz", the output file is the same name without
    the ".gz".  Otherwise the output file is the input file.

    Input file is always removed, even if decompression fails. (But not
    if removing the file fails.)
    """
    import gzip
    from os import unlink
    BUFFER_SIZE = 1024 * 1024 * 16
    if local[-3:] == ".gz":
        output = local[:-3]
    else:
        raise RuntimeError("FIXME: gunzip without a .gz extension isn't currently supported")

    with gzip.open(local, 'rb') as inf:

        # Rely on Unix behavior of open file handles
        # remaining valid when a file is removed
        # Do this before uncompression to handle case
        # where output == local.
        if not preserve_original:
            unlink(local)

        with open(output, 'wb') as out:
            while True:
                buffer = inf.read(BUFFER_SIZE)
                if len(buffer) == 0:
                    break
                out.write(buffer)

    return output


def split_search_path(s):
    '''
    Break a colon-separated list of directories into a list of directories, else empty-list
    '''
    if not s:
        return []

    back_list = []
    for path in s.split(':'):
        back_list.append(Path(path).expanduser().absolute())

    return back_list


def print_ld_lib_path(env_str, path_env):
    for idx, pth in enumerate(path_env.split(':')):
        LOG.debug('{} ({}) = {}'.format(env_str, idx, pth))


def _replaceAll(intputfile, searchExp, replaceExp):
    '''
    Replace all instances of 'searchExp' with 'replaceExp' in 'intputfile'
    '''
    for line in fileinput.input(intputfile, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp, replaceExp)
        sys.stdout.write(line)
    fileinput.close()


def cleanup(objs_to_remove):
    """
    cleanup directories / files
    """
    for file_obj in objs_to_remove:
        try:
            if Path(file_obj).is_dir():
                LOG.debug('\tRemoving directory: {}'.format(file_obj))
                shutil.rmtree(file_obj)
            elif Path(file_obj).is_file():
                LOG.debug('\tRemoving file: {}'.format(file_obj))
                os.unlink(file_obj)
        except Exception:
            LOG.warn("\tUnable to remove {}".format(file_obj))
            LOG.debug(traceback.format_exc())


class AscLineParser(object):
    def time_range(self, ascLine):
        '''
        :param ascLine:
        :return:
        '''
        day, time = self.extract_time_range_tokens(ascLine)
        return self.time_from_tokens(day, time)

    def extract_time_range_tokens(self, ascLine):
        return ascLine.split('"')[3:4][0].split(' ')

    def time_from_tokens(self, day, time):
        dt = datetime.strptime(day + time, '%Y-%m-%d%H:%M:%S.%f')
        return dt


def getURID(URID_timeObj=None):
    '''
    Create a new URID to be used in making the asc filenames
    '''

    URID_dict = {}

    if URID_timeObj is None:
        URID_timeObj = datetime.now(tz=timezone.utc)

    creationDateStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.%f")
    creationDate_nousecStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.000000")

    tv_sec = int(URID_timeObj.strftime("%s"))
    tv_usec = int(URID_timeObj.strftime("%f"))
    hostId_ = uuid.getnode()
    thisAddress = id(URID_timeObj)

    all = tv_sec + tv_usec + hostId_ + thisAddress

    URID = '-'.join(('{0:08x}'.format(tv_sec)[:8],
                     '{0:05x}'.format(tv_usec)[:5],
                     '{0:08x}'.format(hostId_)[:8],
                     '{0:08x}'.format(all)[:8]))

    URID_dict['creationDateStr'] = creationDateStr
    URID_dict['creationDate_nousecStr'] = creationDate_nousecStr
    URID_dict['tv_sec'] = tv_sec
    URID_dict['tv_usec'] = tv_usec
    URID_dict['hostId_'] = hostId_
    URID_dict['thisAddress'] = thisAddress
    URID_dict['URID'] = URID

    return URID_dict


def link_files(dest_path, files):
    '''
    Link ancillary files into a destination directory.
    '''
    files_linked = 0
    for src_file in files:
        src = Path(src_file).name
        dest_file = Path(*(dest_path, src))
        if not Path(dest_file).is_file():
            LOG.debug("Link {0} -> {1}".format(src_file, dest_file))
            Path(dest_file).symlink_to(src_file)
            files_linked += 1
        else:
            LOG.warn('link already exists: {}'.format(dest_file))
            files_linked += 1
    return files_linked


class CsppEnvironment(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def check_and_convert_path(key, a_path, check_write=False):
    '''
    Make sure the path or paths specified exist
    Return the path or list of absolute paths that exist
    '''
    if isinstance(a_path, Path):
        a_path = str(a_path)

    abs_locations = []
    if ":" in a_path:
        paths = a_path.split(":")
        # isinstance(s, string_types)
    elif isinstance(a_path, string_types):
        paths = [a_path]
    else:
        paths = a_path

    for path in paths:
        if not Path(path).exists():
            if key:
                msg = "Environment variable {} refers to a path that does not exist.  {}={}".format(
                    key, key, path)
            else:
                msg = "Required path {} does not exist.".format(path)

            raise CsppEnvironment(msg)
            # sys.exit(2)
            return None
        else:
            LOG.debug("Found: {} at {} {}".format(key, path, Path(path).absolute()))
            abs_locations.append(Path(path).absolute())

        if check_write:
            if not os.access(path, os.W_OK):
                msg = "Path exists but is not writable {}={}".format(key, path)
                raise CsppEnvironment(msg)

    # return a string if only one and an array if more
    if len(abs_locations) == 1:
        return abs_locations[0]
    else:
        # return abs_locations
        # return a :-joined string for use in an env variable
        return ':'.join(abs_locations)


def check_existing_env_var(varname, default_value=None, flag_warn=False):
    '''
    If variable exists then use value, otherwise use default
    '''
    value = None
    if varname in os.environ:
        value = os.environ.get(varname)
    else:
        if default_value is not None:
            value = default_value
        else:
            if flag_warn:
                LOG.warn("{} is not set, please update environment and re-try".format(varname))
                LOG.warn("Environment variable missing. {}".format(varname))
            else:
                LOG.debug("{} is not set, please update environment and re-try".format(varname))
                LOG.debug("Environment variable missing. {}".format(varname))

    return value


def check_and_convert_env_var(varname, check_write=False, default_value=None, flag_warn=False):
    value = check_existing_env_var(varname, default_value=default_value, flag_warn=flag_warn)
    path = check_and_convert_path(varname, value, check_write=check_write)
    return value, path


def what_package_am_i():
    path = Path(__file__).absolute().parent
    cspp_x = path.split("/common")
    cspp_x_home = cspp_x[0]

    return cspp_x_home


def _ldd_verify(exe):
    '''
    check that a program is ready to run
    '''
    rc = call(['ldd', exe], stdout=os.tmpfile(), stderr=os.tmpfile())
    return rc == 0


def env(**kv):
    '''
    augment environment with new values
    '''
    zult = dict(os.environ)
    zult.update(kv)

    return zult


def _convert_datetime(s):
    '''
    converter which takes strings from ASC and converts to computable datetime objects
    '''
    pt = s.rfind('.')
    micro_s = s[pt + 1:]
    micro_s += '0' * (6 - len(micro_s))
    # when = dt.datetime.strptime(s[:pt], '%Y-%m-%d %H:%M:%S').replace(microsecond = int(micro_s))
    when = datetime.strptime(s[:pt], '%Y-%m-%d %H:%M:%S').replace(microsecond=int(micro_s))
    return when


def _convert_isodatetime(s):
    '''
    converter which takes strings from ASC and converts to computable datetime objects
    '''
    pt = s.rfind('.')
    micro_s = s[pt + 1:]
    micro_s += '0' * (6 - len(micro_s))
    # when = dt.datetime.strptime(s[:pt], '%Y-%m-%d %H:%M:%S').replace(microsecond = int(micro_s))
    when = datetime.strptime(s[:pt], '%Y-%m-%dT%H:%M:%S').replace(microsecond=int(micro_s))
    return when


def make_time_stamp_d(timeObj):
    '''
    Returns a timestamp ending in deciseconds
    '''
    dateStamp = timeObj.strftime("%Y-%m-%d")
    # seconds = repr(int(round(timeObj.second + float(timeObj.microsecond) / 1000000.)))
    deciSeconds = int(round(float(timeObj.microsecond) / 100000.))
    deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
    timeStamp = "{}.{}".format(timeObj.strftime("%H:%M:%S"), deciSeconds)
    return "{} {}".format(dateStamp, timeStamp)


def make_time_stamp_m(timeObj):
    '''
    Returns a timestamp ending in milliseconds
    '''
    dateStamp = timeObj.strftime("%Y-%m-%d")
    # seconds = repr(int(round(timeObj.second + float(timeObj.microsecond) / 1000000.)))
    milliseconds = int(round(float(timeObj.microsecond) / 1000.))
    milliseconds = repr(000 if milliseconds > 999 else milliseconds)
    timeStamp = "{}.{}".format(timeObj.strftime("%H:%M:%S"), str(milliseconds).zfill(3))
    return "{} {}".format(dateStamp, timeStamp)


def execution_time(startTime, endTime):
    '''
    Converts a time duration in seconds to days, hours, minutes etc...
    '''

    time_dict = {}

    delta = endTime - startTime
    days, remainder = divmod(delta, 86400.)
    hours, remainder = divmod(remainder, 3600.)
    minutes, seconds = divmod(remainder, 60.)

    time_dict['delta'] = delta
    time_dict['days'] = int(days)
    time_dict['hours'] = int(hours)
    time_dict['minutes'] = int(minutes)
    time_dict['seconds'] = seconds

    return time_dict


class NonBlockingStreamReader:
    '''
    Implements a reader for a data stream (associated with a subprocess) which
    does not block the process. This is done by writing the stream to a queue
    (decoupling the stream from the reading), and then slurping data off of the
    queue and passing it to wherever it's needed.
    '''

    def __init__(self, stream):
        '''
        stream: the stream to read from.
                Usually a process' stdout or stderr.
        '''

        self.stream = stream
        self.queue = Queue()

        def _populateQueue(stream, queue):
            '''
            Collect lines from 'stream' and put them in 'queue'.
            '''

            try:
                while True:
                    line = stream.readline()
                    if line:
                        queue.put(line)
                    else:
                        raise UnexpectedEndOfStream
                        pass
            except UnexpectedEndOfStream:
                LOG.debug("The process output stream has ended.")
            except ValueError:
                LOG.debug("ValueError: The process output stream has ended.")

        self.thread = Thread(target=_populateQueue, args=(self.stream, self.queue))
        self.thread.daemon = True
        self.thread.start()  # start collecting lines from the stream

    def readline(self, timeout=None):
        try:
            return self.queue.get(block=timeout is not None,
                                  timeout=timeout)
        except Empty:
            # print "Need to close the thread"
            return None


class UnexpectedEndOfStream(Exception):
    pass


def execute_binary_captured_inject_io(work_dir, cmd, err_dict, log_execution=True, log_stdout=True,
                                      log_stderr=True, **kv):
    '''
    Execute an external script, capturing stdout and stderr without blocking the
    called script.
    '''

    LOG.debug('executing {} with kv={}'.format(cmd, kv))
    pop = Popen(cmd,
                cwd=work_dir,
                env=env(**kv),
                shell=True,
                stdin=PIPE,
                stdout=PIPE,
                stderr=PIPE,
                close_fds=True)

    # wrap pop.std* streams with NonBlockingStreamReader objects:
    nbsr_stdout = NonBlockingStreamReader(pop.stdout)
    nbsr_stderr = NonBlockingStreamReader(pop.stderr)

    error_keys = err_dict['error_keys']
    del(err_dict['error_keys'])

    # get the output
    out_str = ""
    while pop.poll() is None and nbsr_stdout.thread.is_alive() and nbsr_stderr.thread.is_alive():

        '''
        Trawl through the stdout stream
        '''
        output_stdout = nbsr_stdout.readline(0.01)  # 0.01 secs to let the shell output the result

        if output_stdout is not None:

            # Gather the stdout stream for output to a log file.
            time_obj = datetime.now(tz=timezone.utc)
            time_stamp = make_time_stamp_m(time_obj)
            out_str += "{} (INFO)  : {}".format(time_stamp, output_stdout)

            # Search stdout for exe error strings and pass them to the logger.
            for error_key in error_keys:
                error_pattern = err_dict[error_key]['pattern']
                if error_pattern in output_stdout.decode():
                    output_stdout = string.replace(output_stdout, "\n", "")
                    err_dict[error_key]['count'] += 1

                    if err_dict[error_key]['count_only']:
                        if err_dict[error_key]['count'] < err_dict[error_key]['max_count']:
                            LOG.warn(string.replace(output_stdout, "\n", ""))
                        if err_dict[error_key]['count'] == err_dict[error_key]['max_count']:
                            LOG.warn(string.replace(output_stdout, "\n", ""))
                            LOG.warn(
                                'Maximum number of "{}" messages reached,' +
                                'further instances will be counted only'.format(error_key))
                    else:
                        LOG.warn(string.replace(output_stdout, "\n", ""))
                    break

        '''
        Trawl through the stderr stream
        '''
        output_stderr = nbsr_stderr.readline()  # 0.1 secs to let the shell output the result

        if output_stderr is not None:

            # Gather the stderr stream for output to a log file.
            time_obj = datetime.now(tz=timezone.utc)
            time_stamp = make_time_stamp_m(time_obj)
            out_str += "{} (WARNING) : {}".format(time_stamp, output_stderr)

        '''
        Check to see if the stdout and stderr streams are ended
        '''
        if not nbsr_stdout.thread.is_alive():
            LOG.debug("stdout thread has ended for {}".format(cmd.split(" ")[-1]))
        if not nbsr_stderr.thread.is_alive():
            LOG.debug("stderr thread has ended for {}".format(cmd.split(" ")[-1]))

    # Flush the remaining content in the stdout and stderr streams
    while True:
        try:
            # 0.01 secs to let the shell output the result
            output_stdout = nbsr_stdout.readline(0.01)
            # 0.1 secs to let the shell output the result
            output_stderr = nbsr_stderr.readline()

            if output_stdout is not None or output_stderr is not None:

                if output_stdout is not None:
                    # Gather the stdout stream for output to a log file.
                    time_obj = datetime.now(tz=timezone.utc)
                    time_stamp = make_time_stamp_m(time_obj)
                    out_str += "{} (INFO)  : {}".format(time_stamp, output_stdout)

                if output_stderr is not None:
                    # Gather the stderr stream for output to a log file.
                    time_obj = datetime.now(tz=timezone.utc)
                    time_stamp = make_time_stamp_m(time_obj)
                    out_str += "{} (WARNING)  : {}".format(time_stamp, output_stderr)
            else:
                break

        except IOError:
            pass

    # Poll for the return code. A "None" value indicates that the process hasn’t terminated yet.
    # A negative value -N indicates that the child was terminated by signal N
    max_rc_poll_attempts = 20
    rc_poll_attempts = 0
    continue_polling = True
    while continue_polling:
        if rc_poll_attempts == max_rc_poll_attempts:
            LOG.warn(
                'Maximum number of attempts ({}) of obtaining return code for {} reached,' +
                'setting to zero.'.format(rc_poll_attempts, cmd.split(" ")[-1],))
            rc = 0
            break

        rc = pop.returncode
        LOG.debug("{} : pop.returncode = {}".format(cmd.split(" ")[-1], rc))
        if rc is not None:
            continue_polling = False

        rc_poll_attempts += 1
        time.sleep(0.5)

    LOG.debug("{}: rc = {}".format(cmd, rc))

    return rc, out_str


def get_return_code(num_unpacking_problems, num_xml_files_to_process,
                    num_no_output_runs, noncritical_problem, environment_error):
    '''
    based on problems encountered, print final disposition message, return
    return code to be passed back to caller. Non-zero return code indicates a
    critical problem was encountered.
    '''
    # considered a noncritical problem if there were any runs that crashed,
    # produced no output, where Geo failed, where ADL logs indicated a problem,
    # or where output SDRs failed the imaginary quality check

    # critical problems: set non-zero return code and log error messages
    rc = 0
    if num_unpacking_problems > 0:
        rc |= 2
        LOG.error('Failed to unpack input data.')
    # skipping this check if no XML files to process
    if num_xml_files_to_process and (num_xml_files_to_process <= num_no_output_runs):
        rc |= 1
        LOG.error('Failed to generate any SDR granules.')
    if environment_error:
        rc |= 8
        LOG.error("Environment error.")

    # if critical error was encountered, print failure message and return error code
    if rc != 0:
        LOG.error('Failure. Refer to previous error messages')
        LOG.info('Failure. Refer to previous error messages')
        return rc

    # otherwise no errors or only non-critical errors: print success message and return 0
    if noncritical_problem:
        LOG.info('Normal Completion. Encountered some problems (refer to previous error messages).')
    else:
        LOG.info('Normal Completion.')
    return rc


def create_dir(dir):
    '''
    Create a directory
    '''
    returned_dir = Path(copy(dir))
    LOG.debug("We want to create the dir {} ...".format(dir))

    try:
        if returned_dir is not None:
            returned_dir_path = returned_dir.parent
            returned_dir_base = returned_dir.name
            LOG.debug("returned_dir_path = {}".format(returned_dir_path))
            LOG.debug("returned_dir_base = {}".format(returned_dir_base))
            # Check if a directory and has write permissions...
            if not returned_dir.is_dir() and os.access(returned_dir_path, os.W_OK):
                LOG.debug("Creating directory {} ...".format(returned_dir))
                os.makedirs(returned_dir)
                # Check if the created dir has write permissions
                if not os.access(returned_dir, os.W_OK):
                    msg = "Created dir {} is not writable.".format(returned_dir)
                    raise CsppEnvironment(msg)
            elif returned_dir.is_dir():
                LOG.debug("Directory {} exists...".format(returned_dir))
                if not os.access(returned_dir, os.W_OK):
                    msg = "Existing dir {} is not writable.".format(returned_dir)
                    raise CsppEnvironment(msg)
            else:
                raise CsppEnvironment("Cannot create {}".format(returned_dir))
    except CsppEnvironment:
        LOG.debug("Unable to create {}".format(returned_dir))
        LOG.debug(traceback.format_exc())
        returned_dir = None
    except OSError:
        LOG.debug("Unable to create new dir '{}' in {}".format(
            returned_dir_base, returned_dir_path))
        LOG.debug(traceback.format_exc())
        returned_dir = None
    except Exception:
        LOG.warning("General error for {}".format(returned_dir))
        LOG.debug(traceback.format_exc())
        returned_dir = None

    LOG.debug('Final returned_dir = {}'.format(returned_dir))
    return returned_dir
