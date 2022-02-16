#!/usr/bin/env python

"""
build kraken database for classifying reads

example python3 build_kraken_db.py \\
                --task download-taxonomy \\ 
                --db virus_db \\
                --library viral \\ 
                -t 4 -v
"""

# --- standard python imports ---#
import os
import logging
import argparse
import textwrap
# from pathlib import Path

# --- third-party imports ---#

# --- project specific imports ---#
from utils import mkdir
from utils import find_executable
from utils import run_shell_command

log_level = logging.DEBUG
logging.basicConfig(level=log_level,
                    format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p'
                    )


def parse_args():
    """
    command line arguments

    :return:
    """
    parser = argparse.ArgumentParser(prog="build_kraken_db.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=textwrap.dedent('''\
                                             build Kraken2 database
                                             ---------------------------------------------------------------------------

                                             task: operation to be performed by kraken2
                                             library: library to download sequences from
                                             db: database name/path
                                             '''),
                                     argument_default=argparse.SUPPRESS)
    required_group = parser.add_argument_group('required arguments')

    required_group.add_argument('--task', default='', action='store', metavar='<str>', required=True,
                                choices=["download-taxonomy", "download-library", "special",
                                         "add-to-library", "build", "clean", "standard"],
                                help="operation to be performed by kraken2 "
                                     "['download-taxonomy', 'download-library', 'special', 'add-to-library', 'build', "
                                     "'clean', 'standard'] "
                                )

    required_group.add_argument('--db', required=True, type=str, dest="db",
                                help="database name (path)"
                                )
    required_group.add_argument('--library', metavar='<str>',
                                choices=["archaea", "bacteria", "plasmid",
                                         "viral", "human", "fungi", "plant", "protozoa",
                                         "nr", "nt", "env_nr", "env_nt", "UniVec",
                                         "UniVec_Core"],
                                help="library to download sequences from - ['archaea', 'bacteria', 'plasmid','viral', "
                                     "'human', 'fungi', 'plant', 'protozoa','nr', 'nt', 'env_nr', 'env_nt', 'UniVec',"
                                     "'UniVec_Core']")

    parser.add_argument('-t', '--threads', metavar='<int>', type=int,
                        dest='n_threads', default=4,
                        help="number of threads/cpus to use; specifying the value "
                             "auto' will cause the number of available CPU cores "
                             "on your system, if determinable, to be used"
                        )
    parser.add_argument('--extra-args', default="", type=str, metavar='<str>',
                        help='extra arguments to be passed directly to the kraken2 --download-taxonomy executable '
                             '\n(e.g., --extra-args="--use-ftp")'
                        )
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        dest="verbose",
                        help="verbose output (repeat for increased verbosity)"
                        )
    parser.add_argument('-q', '--quiet',
                        action='count',
                        default=0,
                        dest='quiet',
                        help="quiet output (show errors only)"
                        )
    return parser


def download_taxonomy(db, task="download-taxonomy", n_threads=4, extra_args=""):
    """
    download the accession number to taxon maps, as well as the taxonomic name and tree information from NCBI.

    :param task: <str> operation to be performed
    :param db: <str> path/location to the database
    :param n_threads: <int> number of cpus/threads
    :param extra_args: <str> extra arguments passed to the taxonomy executable
    :return:
    """

    # locate the executable
    kraken2 = find_executable(["kraken2-build"])

    # create db directory
    db = os.path.join(os.path.realpath(db),'')
    mkdir(db)
    taxonomy_dir = os.path.join(db, 'taxonomy')
    # mkdir(taxonomy_dir)
    # prelim = os.path.join(taxonomy_dir, 'prelim_map.txt')
    if os.path.exists(taxonomy_dir) and len(os.listdir(taxonomy_dir)) > 0:
        logging.critical('taxonomy already downloaded \n\t{}'.format(taxonomy_dir))

    else:
        call = ["{} --{} --use-ftp --threads {} --db {} {}".format(kraken2, task, n_threads, db, extra_args)]
        cmd = " ".join(call)
        # run the shell command
        logging.info("downloading taxonomy")
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return db


def download_library(db, library, task="download-library", n_threads=4, extra_args=""):
    """
    download sets of standard genomes/proteins while providing a reference.

    :param db: <str> path/location to the database
    :param library: <str> name of the reference library whose sequences will be downloaded
    :param task: <str> operation to be performed
    :param n_threads: <int> number of cpus/threads
    :param extra_args: <str> extra arguments passed to the taxonomy executable
    :return:
    """

    # locate the executable
    kraken2 = find_executable(["kraken2-build"])

    # create db directory
    db = os.path.join(os.path.realpath(db),'')
    mkdir(db)

    lib_dir = os.path.join(db, 'library', library)
    fnames = ['assembly_summary.txt', 'manifest.txt', 'prelim_map.txt', 'library.fna.masked', 'library.fna']
    ind = list(os.path.join(lib_dir, f) for f in fnames)
    i_bool = [os.path.isfile(filename) for filename in ind]
    result = all(i_bool)

    if result is True:
        logging.critical("library files exist \n\t{}".format(
            '\n\t'.join([filename for filename in ind if os.path.isfile(filename)])))
    else:
        call = ["{} --{} {} --threads {} --use-ftp --db {} {}".format(kraken2, task, library, n_threads, db, extra_args)
                ]
        cmd = " ".join(call)

        # run the shell command
        logging.info("downloading library")
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
        # if not p and os.path.exists(os.path.join(lib_dir, 'assembly_summary.txt')):
        #     assembly_summary = os.path.join(lib_dir, 'assembly_summary.txt')
        #     fix_ftp_paths(assembly_summary)
        #     p = run_shell_command(cmd=cmd, logfile=f_out, raise_errors=False, extra_env=None)
    return db


def build(db, task="build", n_threads=4, extra_args=""):
    """
    build the database using 'kraken2-build --build' once library has been installed

    :param db: <str> path/location to the database
    :param task: <str> operation to be performed
    :param n_threads: <int> number of cpus/threads
    :param extra_args: <str> extra arguments passed to the taxonomy executable
    :return:
    """

    # locate the executable
    kraken2 = find_executable(["kraken2-build"])

    # create db directory
    db = os.path.join(os.path.realpath(db),'')
    mkdir(db)

    # check if indices exist
    indices = ['opts.k2d', 'hash.k2d', 'taxo.k2d', 'seqid2taxid.map']
    ind = list(os.path.join(db, ext) for ext in indices)
    i_bool = [os.path.isfile(filename) for filename in ind]
    result = all(i_bool)

    if result is True:
        logging.critical("indices exist \n\t{}".format(
            '\n\t'.join([filename for filename in ind if os.path.isfile(filename)])))
    else:
        # run process
        call = ["{} --{} --threads {} --db {} {}".format(kraken2, task, n_threads, db, extra_args)
                ]
        cmd = " ".join(call)

        # run the shell command
        logging.info("building database")
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return db


def build_standard(db, task='standard', n_threads=4, extra_args=""):
    """
    build standard kraken2 database
    :param db: <str> path/location to the database
    :param task: <str> operation to be performed
    :param n_threads: <int> number of cpus/threads
    :param extra_args: <str> extra arguments passed to the taxonomy executable
    :return:
    """

    # locate the executable
    kraken2_build = find_executable(["kraken2-build"])

    # create db directory
    db = os.path.join(os.path.realpath(db),'')
    mkdir(db)

    lib_dir = os.path.join(db, 'library')
    if os.path.exists(lib_dir) and len(os.listdir(lib_dir)) > 0:
        logging.critical('library already downloaded -> \n\t{}'.format(lib_dir))

    else:
        # run process
        call = ["{} --{} --threads {} --use-ftp --db {} {}".format(kraken2_build, task, n_threads, db, extra_args)]
        cmd = " ".join(call)

        # run the shell command
        logging.info("building standard database")
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)

    return db


def main():
    parser = parse_args()
    args = parser.parse_args()

    if args.task == 'download-taxonomy':
        download_taxonomy(db=args.db,
                          task=args.task,
                          n_threads=args.n_threads,
                          extra_args=args.extra_args
                          )
    if args.task == 'download-library':
        download_library(db=args.db,
                         library=args.library,
                         task=args.task,
                         n_threads=args.n_threads,
                         extra_args=args.extra_args
                         )
    if args.task == 'build':
        build(db=args.db,
              task=args.task,
              n_threads=args.n_threads,
              extra_args=args.extra_args
              )

    if args.task == 'standard':
        build_standard(db=args.db,
                       task=args.task,
                       n_threads=args.n_threads,
                       extra_args=args.extra_args
                       )


if __name__ == '__main__':
    main()
