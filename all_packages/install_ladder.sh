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
ladder_file   = arguments.ladder_file[0]
#ladder_file   = ladder_file.decode("utf-8")
print(ladder_file)

with open(ladder_file,"r") as ladder:
    os.chdir( ".." )
    cwd = os.getcwd()
    for package in ladder.readlines():
        package = package.strip()
        package_var = f"TACC_{package.upper()}_DIR"
        package_dir = os.getenv( package_var,"" )
        do_install  = force_install \
            or package_dir == "" \
            or not os.path.isdir( package_dir )
        if do_install:
            print( f"Installing package: {package}" )
            list_reqs = subprocess.Popen\
                ( f"cd {cwd}/{packageloc}/ && mpm.py -c {package_config} dependencies",
                  shell=True,env=os.environ,
                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
            reqs = list_reqs.communicate()[0].strip().decode("utf-8").split()
