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
                     'netcdf':'Configuration.seq', 'parallelnetcdf':'Configuration.par',
                     'sundials':'Configuration.par',
                     }
    if package in config_files.keys():
        config_file = config_files[package]
    else: config_file = "Configuration"
    # test on existence is done in the calling environment
    return config_file
