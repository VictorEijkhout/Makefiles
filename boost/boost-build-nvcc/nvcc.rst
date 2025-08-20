Boost.Build ``nvcc`` Toolset
============================

.. contents::

Overview
--------

This project provides a Boost.Build toolset for the NVIDIA ``nvcc``
compiler.

Supported Development Platforms
```````````````````````````````

Every effort has been made to ensure that the toolsets would support
proper operation on all platforms supported by Boost.Build, including
being run from Cygwin on Windows.

The toolset searches for the compiler in the standard location on each
platform.

Usage
-----

Once the toolset is configured it should work like any other
Boost.Build toolset within the constraints of the NVIDIA C/C++
Compiler.

A ``using`` directive without parameters searches for the code
generation tools in the normal places.  The first to match wins.

::

   using nvcc ;

Specifying the version performs the same search as above but stops
with the first toolset found that provides that version number.

::

   using nvcc : 9.2 ;

Specifying the path will use the path specified.  If the version does
not match the desired version, it is an error.

::

   using nvcc : 9.2 : /Developer/NVIDIA/CUDA-9.2/bin/nvcc ;

GPU Support
-----------

The ``nvcc`` toolset adds support for compiling and assembling the
various CUDA-specific file types.

* ``.cu`` - mixed host and device code in C++ and C
* ``.ptx`` - PTX assembling code

The toolset provides two new features in support of GPU architecture
and code.

* ``gpu-architecture``
* ``gpu-code``

These features can be used as any other feature in Boost.Build.

Run the following commands to set GPU architecture, GPU code, or both.

.. code::

   b2 gpu-architecture=compute_70

.. code::

   b2 gpu-architecture=compute_70 gpu-code=sm_70

.. code::

   b2 gpu-code=sm_70

Note that ``nvcc`` also allows multiple GPU architectures or codes to
be set using a comma-separated format.  Currently, the ``nvcc``
toolset does not support this format.
