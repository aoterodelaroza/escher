
Binary instalation *escherpy* dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the easiest and fastest way.

Download ``Anaconda`` python binaries, install it and start using *escherpy*.

You can continue with 

* :ref:`Instalation of the *escherpy* module in Fedora 19 <sec:install_escherpyf19>`



Instalation of *escherpy* dependencies in Fedora 19 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation instructions to run the *escherpy* module.
More or less the same with other distributions.

- Install gcc make automake

::

    sudo yum install make automake gcc gcc-c++

- Install cmake-2.8.12.1 or >=2.8

::

    sudo yum install cmake
    #./bootstrap
    #  make

- Install git

::

    sudo yum install git

- Download VTK

::

    git clone git://vtk.org/VTK.git
    cd VTK
    I'm sure release candidate 3 works. I cannot assure others
    will work. It must be v6 at least.
    git checkout v6.0.0.rc3

- Install dependencies of VTK

::

    sudo yum install freeglut freeglut-devel
    sudo yum install libXt-devel
    sudo yum install python-devel

- Install VTK

::

    ccmake .
    Type c to configure
    Turn ON BUILD_SHARED_LIBS
    It is advisable that you change 
    CMAKE_INSTALL_PREFIX to a folder in your $HOME
    Turn ON VTK_WRAP_PYTHON
    Type c to configure
    Again to be sure:
    Type c to configure
    Type g to generate configuration file
    make
    make install
    echo 'export PYTHONPATH=$PYTHONPATH:/path/to/VTK/bin/' >> $HOME/.bashrc
    echo 'export PYTHONPATH=$PYTHONPATH:/path/to/VTK/lib/' >> $HOME/.bashrc
    echo 'export PYTHONPATH=$PYTHONPATH:/path/to/VTK/Wrapping/Python/' >> $HOME/.bashrc
    echo 'export LD_LIBARY_PATH=$LD_LIBARY_PATH:/path/to/prefix/lib/' >> $HOME/.bashrc

check if it works:

::

    python -c 'import vtk'

- Install numpy

::

    sudo yum install numpy

Forget about the *parse* module by now.
The *parse* python module is planned to be used in the future:

::

    sudo pip install parse

or 

::

    sudo easy_install parse


.. _`sec:install_escherpyf19`:

Instalation of the *escherpy* module in Fedora 19
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Download escher

::

    git clone azufre.quimica.uniovi.es:/home/daniel/src/escher.git
    git clone azufre.quimica.uniovi.es:/home/daniel/src/escher_data.git

Define the following environment variables:

::

    echo 'export ESCHER_DATA=/path/to/escher_data' >> $HOME/.bashrc
    echo 'export PYTHONPATH=$PYTHONPATH:/path/to/escher/escherpy' >> $HOME/.bashrc

Optionally, you can let the log write detailed info about the run,
otherwise set it to ``OFF``:

::

    echo 'export ESCHER_DEBUG="ON"' >> $HOME/.bashrc

Run 

::

    source ~/.bashrc

to refresh your environment variables.

Go to ``$ESCHER_HOME/escherpy/test/`` and run ``make test`` or
``make testcoverage`` to see the module statements covered by the tests.

Go to ``$ESCHER_HOME/escherpy/examples/`` and run ``./testall.py``, or
individually ``./testme.py``:

::

    ./testall.py

If the following error occurs
::

    Traceback (most recent call last):

    File "./testme.py", line 5, in <module>
        import escherpy as esc
    ImportError: No module named escherpy 

the escher path has not been added properly to your ``PYTHONPATH``.

After that, several 3D files will appear and a Blender script
to visualize the VRML file. It can be run with:

::

    blender -P vrml.bpy


Use of escherpy
~~~~~~~~~~~~~~~~~~


Complete set of instructions.

::

    import escherpy as esc

    mol = esc.Molecule()

    mol. structfile = 'file path'
    mol.readstruct()
    mol.rot(60., 30., 20.)
    mol.stickball()

    mol. cpfile = 'file path'
    mol.readcps()
    mol.cpball()

    mol.basinfile = 'file path'
    mol.readsurf('Na')

    mol.isovalue = 0.1
    mol.isosurface('file path')

    mol. densfile = 'file path'
    mol. gradfile = 'file path'
    mol.nciplot()

    mol.show()

