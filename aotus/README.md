Advanced Options and Tables in Universal Scripting
==================================================

The AOTUS library provides a Fortran wrapper around the C-API of the
[Lua](http://www.lua.org) scripting language, allowing a convenient usage of Lua
scripts as configuration files in Fortran applications.

Please visit the
[FORD generated documentation](https://geb.inf.tu-dresden.de/doxy/aotus)
for more information on its usage and details.

*This library is released under a simplified MIT licence, please have a look
 into the COPYRIGHT file for details.*


How To Build
------------

[Waf](http://waf.io/) is used as build system, but alternatively
you may also use a Makefile, see below.

Run:

~~~~~~~~~~~{.sh}
./waf configure build
~~~~~~~~~~~

to build the Aotus library.
If you want to select a specific Fortran compiler, set the environment variable
*FC*.
And for a specific C compiler, set the environment variable *CC*.
The Fortran compiler flags are set with the help of fcopts, which provide
a set of compiler flag combinations for various compilers.
They are found in the `fortran_compiler.py` file in the root directory of the
project.

By running:

~~~~~~~~~~~{.sh}
./waf --help
~~~~~~~~~~~

you get a list of available options to the waf script.

### Generating a Makefile

If you do not want to use waf for the actual compilation of the project,
you can generate a `Makefile` with waf by running the `dump` command.
To resemble more closely the usual auto-tool steps, there is a
`configure` script to run the configuration phase of waf and the
generation of the Makefile.
This will create a Makefile in the build subdirectory, so you can
compile aotus with:

~~~~~~~~~~~{.sh}
./configure
cd build
make
~~~~~~~~~~~

You can get a list of available options by running

~~~~~~~~~~~{.sh}
./configure --help
~~~~~~~~~~~


### Build with the smeka system

The aotus library can also built by using the provided `Makefile.smeka`,
which utilizes the [smeka](https://github.com/zerothi/smeka) build system.

This build system requires the compilation in a subdirectory.
A minimal compilation (defaulting to the GNU compiler suite,
`gfortran`/`gcc`) is achieved by running:

~~~~~{.sh}
mkdir build
cd build
{
echo 'TOP_DIR = ..'
echo 'include $(TOP_DIR)/Makefile.smeka'
} > Makefile
make
~~~~~

and the resulting *libaotus.a* will be put in that directory.  
To control the compiler flags you may create a file `setup.make`,
which can define the usual `FC`, `FFLAGS`, `CC`, `CFLAGS`, `INCLUDES`
and `LIBS` variables which are used for compilation, and linking.
For instance, to compile with just the `-g` flag, you can use:
~~~~~{.sh}
{
echo FFLAGS = -g
echo CFLAGS = -g
} > setup.make
~~~~~

Note that the default compiler is still the GNU suite.

Aotus may be built with two variations of extended precisions. They may
individually be turned on by setting these variables:
~~~~~{makefile}
EXTDOUBLE = 1 # defaults to 0 == do not build with extended double
QUADRUPLE = 1 # defaults to 0 == do not build with quadruple double
~~~~~

To install the library together with the modules in a given directory,
you can run:
~~~~~{.sh}
make install PREFIX=$HOME/aotus
~~~~~
This will create `bin/`, `include/` and `lib/` directories within the
`PREFIX` path with the typical libraries and modules.


_NOTE_: Currently `liblua.a` is not added to the archive, this may easily be
achieved by performing this shell command (after having run `make`)
~~~~~{.sh}
rm libaotus.a
ar -ru libaotus.a *.o
ranlib libaotus.a
~~~~~

What is Built
-------------

For your convenience the Lua library is included in version 5.4.4 (released
2022-01-26).
Its objects are completely gathered into the final *libaotus* library, so you
only need to link against this single static library to gain the
configuration features of Aotus in your Fortran application.
Due to the compiler specific module information required by any application
using libaotus, the suggested approach to incorporate libaotus is to include
its building in the build process of the final application. This is straight
forward if waf is used for the complete project. But also in other build
environments it should not be too hard to make use of the generated *build*
directory.
Yet, if you would rather install the *libaotus.a* and the module files into a
*$PREFIX* directory, you can make use of:

~~~~~~~~~~~{.sh}
./waf install
~~~~~~~~~~~

The default build process will also create some unit test executables and
execute them to ensure functionality of the various parts in the library.

The documentation can be built with
[FORD](https://github.com/Fortran-FOSS-Programmers/ford) by running:

~~~~~~~~~~~{.sh}
ford aot_mainpage.md
~~~~~~~~~~~

This will build a docu directory with the resulting documentation.
Note that this requires
[FORD to be installed](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Installation)
beforehand.
This documentation is also available online at
[our server](https://geb.inf.tu-dresden.de/doxy/aotus).

### Example

There is an example program built, called aotus_sample, which you will find in
the *build* directory.
It can be used with the provided *config.lua* in the *sample* directory, where
also the source of this small program is found.

Getting Started
---------------
The central module in this library is the
[aotus_module](https://geb.inf.tu-dresden.de/doxy/aotus/module/aotus_module.html).
Its documentation and the
[Aotus overview](https://geb.inf.tu-dresden.de/doxy/aotus/page/index.html) would
be good starting points.

Related Projects
----------------

Some projects with similar goals or related information:

* [f2k3-lua](https://github.com/MaikBeckmann/f2k3-lua/tree/simple)
* [FortLua](https://github.com/adolgert/FortLua)

License
=======

Aotus is licensed under the terms of the MIT license reproduced below.
This means that Aotus is free software and can be used for both academic and
commercial purposes at absolutely no cost. You are free to do with the code
whatever you want.
The only requirement is that some credit to the authors is given by putting this
copyright notice somewhere in your project.
The MIT license is chosen for full compatibility with Lua.

For the license of the underlying Lua library have a look at
http://www.lua.org/license.html.

---
Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
                        Aachen and others.
              2014-2018 University of Siegen.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

---
