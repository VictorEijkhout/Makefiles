#!/usr/bin/env python3
# -*- python -*-

import re
import os
import subprocess
import sys

##
## Argument handling
##
args = sys.argv
import argparse
parser = argparse.ArgumentParser\
    ( prog="build_ladder",
      description="Build ladder file",
      add_help=True )
parser.add_argument( '-a','--all',action='store_true',default=False )
parser.add_argument( 'packages', nargs='*', help="root level packages" )
arguments = parser.parse_args()

all_packages = []
after = {}
before = {}
top = []
ignored_packages = \
    [ ".git", ".gitignore", "a1example", "all_packages", 
      "benchpro", "demangle", "demangler", "testing",
      # system stuff:
      "clang", "gcc", "intel", "intel-mpi-binding-kit", "cuda", "mkl", "mpich",
      "gtest", "julia", "tau", "utils",
      # stuff I'm not building:
      "adaptivecpp", "foam", "mapl", "facebook_nle",
      "nethack", "netcdfc", "netcdfx", "netcdff",
      "gklib-karypis", "metis-karypis", # karypis stuff is abandonware
      "neuralfortran", "octopus-auto", "opensycl", "parmetis-git", "plascom",
      "python", "pylauncher",
      "rangev3", "roms", "vtkhdf", "wannier", "wannier90",
      # stuff I should build
      "alps", "athenapk", "cantera",
      "blaspp", "lapackpp", "mfemcuda", # cuda stuff
      "corrfunc", # needs python
      "cesm", "claymore", "ecbuild", "fargparse", "geos",
      "gftl", "gftlshared", # these go together
      "libceed", "ratel", # go together
      "moose", "hpx", "ncl", "nclncarg",
      "qt5", "gnuplot", # gnuplot dpeneds on qt5
      "nanobind", "openblas", "osubenchmark", "plastcom", "petscchaco", "ratel", "rmp",
      "yafyaml", "pflogger", # go together
      "zoltan", # whole build seems sick
    ]

##
## Determine packages to build
##
os.chdir( f"{os.getenv('HOME')}/Software" )

build_all = arguments.all
packages  = arguments.packages
if build_all:
    if len(package)>0:
        print( f"Warning: flagg -a/--all override explicit package specification." )
    packages = [ d for d in os.listdir(".") if os.path.isdir(d) and d not in ignored_packages ]
    print( f"Packages from directory listing: {packages}" )
else:
    print( f"Packages from args: {packages}" )

##
## modules that are built from a different directory
##
variants = ["parallelnetcdf", "parpack", "phdf5",
            "ptscotch", ]
packagedirs = { "parallelnetcdf":"netcdf", "parpack":"arpack", "phdf5":"hdf5",
                "ptscotch":"scotch", }
packagetgts = { "parallelnetcdf":"par", "parpack":"par", "phdf5":"par",
                "ptscotch":"par32", }
def package_dir( package ):
    directories = { 'parallelnetcdf':'netcdf', 'phdf5':'hdf5',
                   }
    if package in directories.keys():
        return directories[package]
    else: return package

def configuration_file( package ):
    config_files = { 'boost':'Configuration.seq',
                     'dealii':'Configuration.real',
                     'hdf5':'Configuration.seq', 'phdf5':'Configuration.par',
                     'netcdf':'Configuration.seq', 'parallenetcdf':'Configuration.par',
                     'sundials':'Configuration.par',
                     }
    if package in config_files.keys():
        config_file = config_files[package]
    else: config_file = "Configuration"
    # test on existence is done in the calling environment
    return config_file

##
## Loop over package directories and variants
##
def closure( packages,before ):
    # assume we are home/Software
    cwd = os.getcwd()
    changes = False
    for p in packages:
        if p in before.keys(): continue
        # only look at uninvestigated packages
        changes = True
        before[p] = []
        packageloc    = package_dir(p)
        package_config = configuration_file( packageloc )
        config_full_path = f"{cwd}/{packageloc}/{package_config}"
        if not os.path.exists( config_full_path ):
            raise Exception( f"No such config file: {config_full_path}" )
        list_reqs = subprocess.Popen\
            ( f"cd {cwd}/{packageloc}/ && mpm.py -c {package_config} dependencies",
              shell=True,env=os.environ,
              stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
        reqs = list_reqs.communicate()[0].strip().decode("utf-8").split()
        print( f"{p} <- {reqs}" )
        for q in reqs:
            if q not in packages:
                packages.append(q)
            before[p].append(q)
    if changes:
        return closure(packages,before)
    else:
        return packages,before

packages,before = closure( packages,{} )
print( f"{packages} with relations: {before}" )
sys.exit(0)

print( f"Before: {before}" )
#print( f"Top: {top}" )
#print( f"After: {after}" )

from collections import defaultdict, deque

def topological_sort(strings, before):
    # Build graph and in-degree count
    graph = defaultdict(list)
    in_degree = defaultdict(int)

    # Initialize in-degree for all nodes
    for s in strings:
        in_degree[s] = 0

    # Add edges based on the 'before' relationships
    for s1, predecessors in before.items():
        #print( f"{s1} <= {predecessors}" )
        for s2 in predecessors:
            graph[s2].append(s1)   # s2 -> s1 (s2 must come before s1)
            in_degree[s1] += 1     # s1 has one more dependency

    # Start with all nodes that have no incoming edges
    queue = deque([s for s in strings if in_degree[s] == 0])
    #print( f"start queue: {queue}" )
    sorted_strings = []

    while queue:
        current = queue.popleft()
        sorted_strings.append(current)
        for neighbor in graph[current]:
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    if len(sorted_strings) != len(strings):
        print( f"Missing: { set(strings)-set(sorted_strings) }" )
        raise ValueError("Cycle detected or incomplete input: cannot topologically sort.")

    return sorted_strings

print( "sorting" )
sorted = topological_sort( all_packages,before )
print(sorted)

ladder = "all_packages/all_ladder.sh"
with open( ladder,"w" ) as ladderfile:
    for s in sorted:
        ladderfile.write( s+'\n' )
print( f"Written ladder to: {ladder}" )

import sys
sys.exit(0)

for package in packages :
    if package in ignored_packages : continue
    all_packages.append(package)
    if os.path.isdir(package):
        # there is a directory for this package
        packagedir = package
        packagetgt = ""
    else:
        # this package is built from a different directory
        packagedir = packagedirs[package]
        packagetgt = f"TARGET={packagetgts[package]}"
    # find prerequisites by running "make listmodules"
    list_prereqs = subprocess.Popen\
        ( f"cd {packagedir}/ && make listmodules {packagetgt}",
          shell=True,env=os.environ,
          stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    prereqs = list_prereqs.communicate()[0].strip().decode("utf-8").split()
    # filter ignored packages
    prereqs = set(prereqs) - set(ignored_packages)
    # remove explicit version numbers
    prereqs = [ p.split("/")[0] for p in prereqs ]
    print( f"Package: {package}, prereqs: {prereqs}" )
    before[package] = prereqs
    if len(prereqs)==0:
        top.append(package)
    else:
        for p in prereqs:
            if p not in after:
                after[p] = []
            after[p].append(package)
