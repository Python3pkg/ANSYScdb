ANSYScdb v0.14
========

Python/Cython module to convert ANSYS CDB ASCII block files to numpy and
vtk unstructured grid objects.

Installation
------------

From PyPI https://pypi.python.org/pypi/ANSYScdb

``pip install ANSYScdb``

From source

``python setup.py install``

License
-------

ANSYScdb is licensed under the MIT license. The full statement is
provided in the file named ``LICENSE``.

Dependencies
------------

Required: ``numpy``, ``cython``. Optional: ``vtk``

Minimum requirements are numpy and cython to parse an ANSYS cdb. To
convert the raw data to a VTK unstructured grid, vtk 5.0 or greater must
be installed with Python bindings.

Tests
-----

Test installation with the following from Python

.. code:: python

    from ANSYScdb import Tests

    # Load a hexahedral beam from a cdb
    Tests.Load.HexBeam() # returns a cdb object

    # Load and display the same beam
    Tests.Display.HexBeam() # returns a cdb object

Example Code
------------

Assumes you have downloaded the example CDB files. Otherwise, replace
‘Beam.cdb’ with your own blocked \*.cdb file.

.. code:: python

    #Load module
    from ANSYScdb import CDB_Reader

    #Load ANSYS cdb file
    cdb = CDB_Reader.Read('HexBeam.cdb')

    # Parse the raw data into a VTK unstructured grid
    cdb.ParseVTK() # Should the cython reader fail, use: cdb.ParseVTK(use_cython=False)

    # Plot the resulting unstructured grid
    cdb.Plot()
