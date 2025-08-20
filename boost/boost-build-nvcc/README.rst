Boost.Build Toolset for NVIDIA nvcc Compiler
============================================

.. contents::

Overview
--------

This project contains a Boost.Build toolset for the NVIDIA C/C++
Compuler (``nvcc``) compiler.

Toolset
-------

See the `documentation <nvcc.rst>`__ for a description of how to use
NVIDIA compilers.

See the implementation files for the toolset for implementation
details for toolset.

* `NVIDIA C/C++ Compiler <nvcc.jam>`__

Testing
-------

This project can be tested using the `Boost.Build Toolset Tester
<https://github.com/tee3/boost-build-toolset-tester>`__.  There is a
branch ``devel-nvcc``, which configures the toolset.

.. code:: sh

   $ b2 --test-config=user-config.jam \
         toolset=nvcc \
         link=static
