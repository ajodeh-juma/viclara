#!/usr/bin/env python

"""
utility functions
"""

import os
import re
import sys
import errno
import shutil
import subprocess
from os import path
from itertools import groupby


def find_executable(names, default=None):
    """
    find an executable PATH from the given list of names.
    Raises an error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.

    :param names: list of given executable names
    :param default:
    :return: <str> path to the first executable found in PATH from the given list of names.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of {} in PATH={}".format(names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception
    return exe


def run_shell_command(cmd, raise_errors=False, extra_env=None):
    """
    run the given command string via Bash with error checking

    :param cmd: command given to the bash shell for executing
    :param raise_errors: bool to raise error if running command fails/succeeds
    :param extra_env: mapping that provides keys and values which are overlayed onto the default subprocess environment.
    :return:
    """

    env = os.environ.copy()

    if extra_env:
        env.update(extra_env)

    try:
        subprocess.check_call("set -euo pipefail; " + cmd,
                                  shell=True,
                                  universal_newlines=True,
                                  executable="/bin/bash"
                                  )
    except (subprocess.CalledProcessError, OSError) as error:
        rc = error.returncode
        if rc == 127:
            extra = "Are you sure this program is installed?"
        else:
            extra = " "
        print("Error occurred: shell exited with return code: {}\ncommand running: {} {}".format(
            error.returncode, cmd, extra)
        )
        if raise_errors:
            raise
        else:
            return False
    except OSError as error:
        print("Unable to run shell command using {}! tool requires {} to be installed.".format(
            error.filename, error.filename)
        )
        if raise_errors:
            raise
        else:
            return False
    else:
        return True


def mkdir(directory):
    """
    recursivley create a directory if it does not exist

    :param directory: path to the directory to be created
    :return: directory
    """
    directory = os.path.abspath(directory)
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    return directory


def fasta_iterator(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    Author: Brent Pedersen
    https://www.biostars.org/p/710/
    """
    fa_fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fa_fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq
    fa_fh.close()


def copy_file(src, dest):
    """
    copy file from one directory to another

    :param src: path of the source file
    :param dest: path of the destination
    :return:
    """
    # copy to directory
    try:
        shutil.copy(src, dest)

    # source and destination are same
    except shutil.SameFileError:
        print("source\t{} and destination\t{} represents the same file.".format(src, dest))

    # permission error
    except PermissionError:
        print("Permission denied!")

    # any other errors
    except Exception as error:
        print(error)
    return dest