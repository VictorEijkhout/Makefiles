#!/usr/bin/env python3
# -*- python -*-

import os
import subprocess
import sys

##
## Argument handling
##
args = sys.argv

import argparse
parser = argparse.ArgumentParser\
    ( prog="install_ladder",
      description="Install according to ladder file",
      add_help=True )
parser.add_argument( '-f','--force',action='store_true',default=False )
parser.add_argument( 'ladder_file', nargs='*', help="ladder file" )
arguments = parser.parse_args()

force_install = arguments.force 
if force_install:
    print( "Forcing all install" )
ladder_file   = arguments.ladder_file[0]
#ladder_file   = ladder_file.decode("utf-8")
print( f"Using ladder file: {ladder_file}" )

from utils import *
from MrPackMod.process import process_initiate,process_terminate, process_execute, nonnull

non_packages = [ "mkl","nvpl","blaslapack", "mpi", ]
install_process = process_initiate()
with open(ladder_file,"r") as ladder:
    os.chdir( ".." )
    cwd = os.getcwd()
    for package in ladder.readlines():
        package = package.strip()
        if package in non_packages:
            print( f" .. skip install of non-package: {package}" )
            continue
        print( f"Package: {package}" )
        if force_install:
            do_install = True
        else:
            package_var = f"TACC_{package.upper()}_DIR"
            installation_dir = os.getenv( package_var,"" )
            if installation_dir != "":
                print( f" .. loaded at {installation_dir}: {os.path.isdir( installation_dir )}" )
            else:
                avails = process_execute( f"module -t avail {package}",terminal=None )
                if nonnull( avails):
                    print( f" .. availability: {avails}; loading" )
                    process_execute( f"module -t load {package}",
                                     process=install_process )
                else:
                    print( " .. proceeding with installation" )
                    packageloc     = package_dir(package)
                    package_config = configuration_file( package )
                    process_execute\
                        ( f"cd {cwd}/{packageloc}/ && mpm.py -c {package_config} install",
                          process=install_process,immediate=True )
                    process_execute( f"module -t load {package}",process=install_process )

process_terminate(install_process)
