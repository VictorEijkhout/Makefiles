#!/usr/bin/env python3
# -*- python -*-

import re
import os
import subprocess

print( "inventory" )
packages = [ d for d in os.listdir(".") if os.path.isdir(d) ]
all_packages = []
after = {}
before = {}
top = []
ignored_packages = \
    [ ".git", ".gitignore", "a1example", "all_packages", 
      "benchpro", "demangle", "demangler", "testing",
      # system stuff:
      "clang", "gcc", "intel", "intel-mpi-binding-kit", "cuda", "mkl", "mpich",
      "gtest", "julia", "tau",
      # stuff I'm not building:
      "adaptivecpp", "foam", "mapl", "nethack", "netcdfx", "facebook_nle", 
      "gklib-karypis", "metis-karypis", # karypis stuff is abandonware
      "octopus-auto", "opensycl", "parmetis-git", "python", "pylauncher", "wannier",
      # stuff I should build
      "alps", "athenapk",
      "blaspp", "lapackpp", "mfemcuda", # cuda stuff
      "cesm", "claymore", "ecbuild", "fargparse", "geos",
      "gftl", "gftlshared", # these go together
      "libceed", "ratel", # go together
      "gmp", "moose", "hpx", "ncl", "nclncarg",
      "qt5", "gnuplot", # gnuplot dpeneds on qt5
      "nanobind", "openblas", "osubenchmark", "petscchaco", "ratel", "rmp",
      "yafyaml", "pflogger", # go together
      "zoltan", # whole build seems sick
    ]

##
## modules that are built from a different directory
##
variants = ["parallelnetcdf", "parpack", "phdf5", ]
packagedirs = { "parallelnetcdf":"netcdf", "parpack":"arpack", "phdf5":"hdf5", }
packagetgts = { "parallelnetcdf":"par", "parpack":"par", "phdf5":"par", }
packages = packages+variants

##
## Loop over package directorie and variants
##
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
    # print( f"Package: {package}, prereqs: {prereqs}" )
    before[package] = prereqs
    if len(prereqs)==0:
        top.append(package)
    else:
        for p in prereqs:
            if p not in after:
                after[p] = []
            after[p].append(package)

print( f"Top: {top}" )
print( f"Before: {before}" )
print( f"After: {after}" )

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
        for s2 in predecessors:
            graph[s2].append(s1)   # s2 -> s1 (s2 must come before s1)
            in_degree[s1] += 1     # s1 has one more dependency

    # Start with all nodes that have no incoming edges
    queue = deque([s for s in strings if in_degree[s] == 0])
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

with open( "all_packages/all_ladder.sh","w" ) as ladderfile:
    for s in sorted:
        ladderfile.write( s+'\n' )

import sys
sys.exit(0)

##
## old code
## broken
##

def list_packages_and_prereqs( package_list,writefun,indent="" ):
    global done,before,after
    for package in package_list:
        if package in done: continue
        print( f"{package} visiting" )
        # first do all prereqs
        if package in before.keys(): # and before[package] is not []:
            # package has prerequisites
            print( f"{package} before: {before[package]}" )
            for prereq in before[package]:
                list_packages_and_prereqs( before[package],writefun,indent=indent )
            before.pop(package)
        # then do package
        if package in done: continue
        print( f"{package} listed" )
        writefun( f"{indent}{package}" )
        done.append(package)
        # then do dependents
        if package in after.keys():
            print( f"{package} after: {after[package]}" )
            list_packages_and_prereqs( after[package],writefun,indent=indent+"  " )
        print( f"{package} fully finished" )

#done = []
#list_packages_and_prereqs( top,lambda x:print(x) )

done = []
with open( "all_packages/all_ladder.sh","w" ) as ladderfile:
    list_packages_and_prereqs( top,lambda x:ladderfile.write(x+'\n') )
