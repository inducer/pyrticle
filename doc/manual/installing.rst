.. highlight:: sh

Installing Pyrticle
===================

This tutorial will walk you through the process of installing
:mod:`pyrticle`. To follow, you need a few basic standard software 
packages:

* A UNIX-like machine with web access.
* A C++ compiler, preferably a Version 4.x gcc.
* A working `Python <http://www.python.org>`_ installation, 
  Version 2.4 or newer. If you are using Python 2.4, 
  you should also install `PySqlite <http://pysqlite.org>`_.
  (Sqlite is included in Python 2.5 and up.)
* The `Basic Linear Algebra Subroutines (BLAS) <http://netlib.org/blas>`_
  or a tuned implementation thereof, and
* `LAPACK <http://netlib.org/lapack>`_.

:mod:`pyrticle` has a number of prerequisites that need to be
installed for it to be usable. Once these prerequisites are installed,
installing this package is fairly straightforward, because it will
simply reuse the build configuration that you create for its 
prerequisites.

.. note::

    Some of the :mod:`pyrticle`'s prerequisites have overlapping 
    dependencies. For example, nearly every package depends on
    Boost and/or numpy. If multiple packages call for the same 
    dependency, installing that dependency only once is sufficient.

Step 1: Install :mod:`hedge`
----------------------------

Navigate to the `installation documentation
<http://documen.tician.de/hedge/installing.html>`_ for the
`hedge <http://mathema.tician.de/software/hedge>`_ Discontinuous
Galerkin solver and follow the instructions there. Since pyrticle
requires the BLAS, you may want to make sure that you also compile
hedge with support for the BLAS, as documented.

Step 2: Install :mod:`pylo`
---------------------------

To make its simulation state accessible to visualization software,
pyrticle relies on the `Silo <https://wci.llnl.gov/codes/silo/>`_
data exchange format. In order to write Silo files, :mod:`pyrticle`
requires `Pylo <http://mathema.tician.de/software/pylo>`_, a Python
interface to the standard Silo data access library.

To install Pylo, please navigate to its `installation tutorial
<http://documen.tician.de/pylo/installing.html>`_.

Step 3: Install Boost.Unordered
-------------------------------

Download this `snapshot <http://tiker.net/tmp/boost-unordered.tar.gz>`_ of the
Boost Unordered Container Library, and unpack it::

    $ tar xfz boost-unordered.tar.gz

Then add the :file:`boost-unordered` directory that was unpacked as the root of
the archive's contents to the `BOOST_INC_DIR` variable in your
:file:`.aksetup-defaults.py` file, located in your home directory. That is, if
this variable used to read::

    BOOST_INC_DIR = ['/home/andreas/pool/include/boost-1_35']

change it to read::

    BOOST_INC_DIR = [
      '/home/andreas/pool/include/boost-1_35',
      '/home/andreas/boost-unordered',
      ]

Step 4: Tell :mod:`pyrticle` about LAPACK
-----------------------------------------

Add the following lines to the file :file:`.aksetup-defaults.py` in your home
directory::

    HAVE_LAPACK = True
    LAPACK_LIB_DIR = ['/where/your/lapack/is/installed']
    LAPACK_LIBNAME = ['your_lapack_libname'] # without leading lib and trailing .a/.so

You will need to customize these lines for your LAPACK install locations.

Step 5: Obtain and Unpack :mod:`pyrticle`
-----------------------------------------

We will assume here that you have obtained a source code snapshot of 
:mod:`pyrticle`. Place this file in the directory from where you wish
to install :mod:`pyrticle` and type::

    $ tar xfz pyrticle.git-VERSION.tar.gz

or::

    $ tar xf pyrticle.git-VERSION.tar.gz

(Whenever you see the "`$`" dollar sign in this tutorial, this means
you should enter this at your shell prompt. You don't have to be
`root`. A few spots are marked with "sudo" to show that these *do*
require root privileges *if* you are using a Python interpreter that
is install globally.)

.. warning::

    While the file will be named `.tar.gz`, your browser may
    already have removed the gzip compression, and left an archive in
    uncompressed tar format. One of the above commands will work, but
    which depends how you downloaded the archive.

Step 6: Build and Install :mod:`pyrticle`
-----------------------------------------

Actually compiling and installing :mod:`pyrticle` should now be fairly
simple::

    $ cd pyrticle.git 
    $ sudo python setup.py install

Get some coffee while :mod:`pyrticle` is installed. Note that depending on the
software versions you have installed, this may print some (actually quite a
bunch of) warnings. Most of these are fixed in newer versions of Python and
Boost. This is no cause for concern--the resulting software will work just
fine.  If you get no errors, congratulations! You have successfully built
:mod:`pyrticle`.

Success! So what now?
---------------------

Next, we suggest that you continue on to the next section and try your hand at
running some simple simulations.
