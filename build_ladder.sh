#!/usr/bin/env python3
# -*- python -*-

import os
import subprocess

packages = os.listdir(".")
after = {}
top = []
for package in packages:
    if package in [ ".git", "a1example", "all_packages", "demangle", "demangler",
                    "clang", "gcc", "intel", "intel-mpi-binding-kit",
                    "python", "pylauncher",
                   ]: continue
    if not os.path.isdir(package): continue
    list_prereqs = subprocess.Popen\
        ( f"cd {package} && make listmodules",
          shell=True,env=os.environ,
          stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    prereqs = list_prereqs.communicate()[0].strip().decode("utf-8").split()
    # print( f"Package: {package}, prereqs: {prereqs}" )
    if len(prereqs)==0:
        top.append(package)
    else:
        for p in prereqs:
            if p not in after:
                after[p] = []
            after[p].append(package)

print( f"Top: {top}" )
print( f"After: {after}" )

def list_packages_and_prereqs( package_list,writefun,indent="" ):
    global done
    for package in package_list:
        if package in done: continue
        done.append(package)
        writefun( f"{indent}{package}" )
        if package in after.keys():
            list_packages_and_prereqs( after[package],writefun,indent=indent+"    " )

done = []
list_packages_and_prereqs( top,lambda x:print(x) )
done = []
with open( "all_packages/all_ladder.sh","w" ) as ladderfile:
    list_packages_and_prereqs( top,lambda x:ladderfile.write(x+'\n') )
