
librsb README file
================================================================================


librsb - Recursive Sparse Blocks  Matrix computations library
--------------------------------------------------------------------------------

 A library for sparse matrix computations featuring the Recursive Sparse Blocks
 (RSB) matrix format, geared towards cache efficient and multi-threaded
 (that is, shared memory parallel) operations on large sparse matrices.
 It provides the most common operations necessary to iterative solvers, like
 matrix-vector and matrix-matrix multiplication, triangular solution, scaling of
 rows/columns, diagonal extraction / setting, blocks extraction, computation of
 norm, formats conversion.
 The RSB format is especially well suited for symmetric and transposed 
 multiplication variants.
 Much of the source code is machine-generated; the supported numerical types can
 be chosen by the user at build time.
 This library is dual-interfaced: it can be used via the native (`RSB`) 
 interface (with function identifiers prefixed by `rsb_` or `RSB_`), and a 
 Sparse BLAS one (function identifiers prefixed by `BLAS_`).
 The `RSB` interface can be used from C/C++ (`rsb.h` header) or via modern 
 Fortran ISO-C-BINDING (`rsb` module in `rsb.F90`).
 New in version 1.3 is the C++-only API in `rsb.hpp` (`Rsb`-prefixed classes).
 The Sparse BLAS interface is usable from C/C++ via the `blas_sparse.h` header, 
 and from Fortran via the `blas_sparse` module.


--------------------------------------------------------------------------------

 This (README) is the first document you should read about librsb.
 It contains basic instructions to configure, build, install, and use librsb.
 The reference documentation for programming with librsb is contained in the
 `./doc/` source package subdirectory and when installed, placed in the
 appropriate system directories as both Unix man pages (`./doc/man/`) and HTML
 (`./doc/html/`).
 If you have used a previous release version of librsb, see the `NEWS` file 
 for a succinct list of user-relevant changes.


## <a id="introduction"></a>INTRODUCTION ##
--------------------------------------------------------------------------------

 librsb is a library for sparse matrix algebra computations.
 It is stand-alone: it does not require any other library to build or work.
 It is shared memory parallel, using constructs from the OpenMP standard.
 It is focused on high performance SpMV and SpMM operations.
 A part of the library code is automatically generated from templates and
 macros, on the basis of the numerical types a user wishes to have supported.
 The configure script provides many build time options, especially with respect
 to debug and additional verbosity (see `configure --help`).
 Defaults should be fine for most users, though.

 - [INTRODUCTION](#introduction)
 - [MAIN ASPECTS,FEATURES](#mainaspectsfeatures)
 - [QUICK INSTALL AND TESTING EXAMPLE](#quickinstallandtestingexample)
 - [LIBRARY CONFIGURATION, GENERATION, BUILD ](#libraryconfigurationgenerationbuild)
 - [INSTALLATION, USAGE](#installationusage)
 - [EXECUTION AND ENVIRONMENT VARIABLES](#executionandenvironmentvariables)
 - [DOCUMENTATION, EXAMPLES AND PROGRAMMING GUIDELINES](#documentationexamplesandprogrammingguidelines)
 - [CONFIGURE, BUILD AND BENCHMARK EXAMPLE](#configurebuildandbenchmarkexample)
 - [COMPATIBILITY](#compatibility)
 - [FAQ](#faq)
 - [POSSIBLE FUTURE ENHANCEMENTS](#possiblefutureenhancements)
 - [ABOUT THE INTERNALS](#abouttheinternals)
 - [BUGS](#bugs)
 - [CONTACTS](#contacts)
 - [CREDITS](#credits)
 - [LICENSE](#license)


## <a id="mainaspectsfeatures"></a>MAIN ASPECTS,FEATURES ##
--------------------------------------------------------------------------------

 * conceived for efficient multithreaded SpMV/SpMM
 * thread level (shared memory) parallelism by using OpenMP
 * threads/structure auto-tuning feature for additional performance
 * support for multiple numerical data types which can be turned
   on/off individually (e.g.:`double`, `float`, `int`, `char`, `float complex`, 
   `double complex`) at configure time
 * a Sparse BLAS interface for matrix assembly, computation, destruction
 * code generators for the inner CSR, COO computational kernels
 * based on a recursive memory layout of submatrices
 * enough functionality to implement the most common iterative methods 
 * basic input sanitizing (index types overflow checks, etc)
 * parallel matrix assembly and conversion routines
 * auxiliary functions for matrix I/O (using the "Matrix Market" format:
   real, integer, complex and pattern are supported)
 * implemented as a building block for solvers like e.g. PSBLAS, MaPHyS++
 * dual implementation of inner kernels: with "full-" and "half-word" indices
 * basic (unoptimized) sparse matrices multiplication and summation
 * interactive usage possible by using the `sparsersb` plugin for GNU Octave 
   or the PyRSB Python package
 * complete with examples and a test suite
 * see the `NEWS` text file for a detailed list of changes in each release


## <a id="quickinstallandtestingexample"></a>QUICK INSTALL AND TESTING EXAMPLE ##
--------------------------------------------------------------------------------

	# unpack the archives or get them from the repositories
	./autogen.sh	# only necessary if  configure  file does not exist
	./configure --prefix=$HOME/local/librsb/ # with installation destination
        # see also ./configure --help for many other options
        # see      ./configure --help=recursive for even more options
	# librsb has been configured
	make help	# provide information
	make		# build the library and test programs
	# librsb has been built
        make  qtests	# perform brief sanity tests
        make qqtests	# the same, but with less output
        make   tests	# perform extended sanity tests
	ls examples/*.c   # editable C       examples; build them with 'make'
	ls examples/*.cpp # editable C++     examples; build them with 'make'
	ls examples/*.F90 # editable Fortran examples; build them with 'make'
	make install	# install to $HOME/local/librsb/
	make itests     # perform post-installation tests on examples
	# librsb has been installed and can be used

	# for instance, try using one of the librsb examples as a base:
	mkdir -p ~/rsb-test/ && cp examples/hello.c ~/rsb-test/myrsb.c
	# adapt hello.c to your needs and recompile:
	cd ~/rsb-test/
	export PATH=$PATH:$HOME/local/librsb/bin/
	gcc `librsb-config --I_opts`.  -c myrsb.c 
 	gcc -o myrsb myrsb.o `librsb-config --static --ldflags --extra_libs`
 	./myrsb         # run your program


LIBRARY CONFIGURATION, GENERATION, BUILD
--------------------------------------------------------------------------------

 A good portion of this library is C code generated via M4 preprocessor macros.
 Certain generation options (notably, supported numerical types) can be set via
 the `./configure`  script.
 The location of the M4 executable can be specified via `./configure M4=...`
 After invoking `./configure` and before running `make` it is possible to invoke
 `make cleanall` to make sure that auto-generated code is deleted first.
 
 The `configure` script attempts at detecting the system cache memory hierarchy
 parameters and prints them out.
 If detection fails, you can pass them via the `--with-memhinfo=...`  option.
 For instance, declaring having 3 MB of L3, 256 kB of L2, and 32 kB of L1:
    `./configure --with-memhinfo="L3:12/64/3072K,L2:8/64/256K,L1:8/64/32K" `
 These values need not be exact: also approximate can help achieving better
 performance.
 You can override the configure-time default at run time: see API documentation.
 Read a few sections further for a description of the memory hierarchy info
 string format.

 If you want to enable Fortran examples, be sure of running `./configure` with
 the `--enable-fortran-examples` option.  You can specify the desired Fortran
 compiler and compilation flags via the `FC` and `FCFLAGS` variables.

 Set the `CPPFLAGS` variable at configure time to provide additional compilation
 flags; e.g. `configure` to detect necessary headers in non-standard location.
 Similarly, the `LDFLAGS` variable can be set to contain link time options; so
 you can use it to specify libraries to be linked to librsb examples.
 Invoke `./configure --help` for details of other relevant environment
 variables.
 
 The `./configure` script will emit information about the current build options.
 If all went fine, you can invoke `make` to build the library and the examples.

 To check for library consistency, run:

	make qtests # takes a short time, spots most problems
or

	make tests  # takes longer, qtests is usually enough
 
 If these tests terminate with an error code, it may be that it has been caused
 by a bug in librsb, so please report it (see BUGS).
 These tests also check error reporting capabilities, so don't be scared by 
 error messages appearing in the running output.


## <a id="installationusage"></a>INSTALLATION, USAGE ##
--------------------------------------------------------------------------------
 
 Once built, the library can be installed with e.g.:

	sudo make install	# 'sudo' may be needed for system-wide locations

 This installs header files, binary library files, documentation, examples,
 and the `librsb-config` program.
 Then, application C programs should include the `rsb.h` header file with
 `#include <rsb.h>`
 C++ programs can use that as well, or `#include <rsb.hpp>`.
 Path to the librsb headers and extra options can be obtained via
 `librsb-config --I_opts`.

 To link to the librsb library and its dependencies one can use the output of
 `librsb-config --static --ldflags --extra_libs`.
 
 Users of `pkg-config` can manually copy the `librsb.pc` file to the appropriate
 directory to use `pkg-config` in a way similar to `librsb-config`.


## <a id="executionandenvironmentvariables"></a>EXECUTION AND ENVIRONMENT VARIABLES ##
--------------------------------------------------------------------------------
 
 By default, librsb reads the environment variable
 `RSB_USER_SET_MEM_HIERARCHY_INFO`, which allows to override settings of memory
 hierarchy set at configure-time or detected at runtime.

 Its value is specified as n concatenated strings of the form:

	L<l>:<a_l>/<b_l>/<c_l>

 These strings are separated by a comma (","), and each of them is made
 up from substrings where:

	<n> is the cache memories hierarchy height, from 1 upwards.
	<l> is the cache level, from 1 upwards.
	<a_l> is the cache associativity
	<b_l> is the cache block size (cache line length)
	<c_l> is the cache capacity (size)

 The `<a_l>`, `<b_l>`, `<c_l>` substrings consist of an integer number with an
 optional multiplier character among {K,M,G} (to specify respectively 2^10,
 2^20 or 2^30).
 Any value is permitted, a long as it is positive. Higher level cache
 capacities are required to be larger than lower level ones.
 Please note that currently, only the cache capacity value is being used.
 Example strings and usage in the BASH shell:

	RSB_USER_SET_MEM_HIERARCHY_INFO="L2:4/64/512K,L1:8/64/32K"  <your program>
	RSB_USER_SET_MEM_HIERARCHY_INFO="L1:8/128/2M"  <your program>

 Experimenting with this environment variable can help tuning performance.

 Setting this environment variable may be also needed if automatic detection
 fails (e.g. very recent systems).

 A default value for this memory hierarchy info string can be set at configure
 time by using the  `--with-memhinfo`  configure option.

 If you don't know values for these parameters, you can run the
 `./scripts/linux-sys-cache.sh`
 script to try to get a guess on a Linux system.
 On other systems, please consult the available documentation.
 E.g.: On Mac OS 10.6 it was possible to get this information by invoking
  `sysctl -a | grep cache`.
  
 You can control the active threads count with OpenMP variable 
 `OMP_NUM_THREADS`, but also override it via `RSB_NUM_THREADS`.


## <a id="documentationexamplesandprogrammingguidelines"></a>DOCUMENTATION, EXAMPLES AND PROGRAMMING GUIDELINES ##
--------------------------------------------------------------------------------

 The `<rsb.h>` header file specifies the C API of librsb.
 Header `<rsb.hpp>` provides two C++ classes for librsb and most of the
 functionality of `<rsb.h>`.
 So these two files are good starting points to learn to use librsb.
 
 To generate librsb API documentation, you need to run configure with the
 `--enable-doc-build` option. This requires the `doxygen` tool.
 Documentation will be placed in `doc/man/` and `doc/html`, and then
 installed with `make install`.
 Once installed and indexed, librsb man pages can be listed with
 `apropos rsb` or `man -k rsb`.

 The latest API documentation can also be found on http://librsb.sourceforge.net

 The `examples` directory has a number of working example programs.

 The library declares symbols prefixed by `rsb_`.
 To avoid name clashes, you should avoid declaring identifiers prefixed that way
 in programs using librsb.

 If `configure` has been invoked with the `--enable-sparse-blas-interface`, then
 the corresponding `BLAS_`- and `blas_`- prefixed symbols will also be built.


## <a id="configurebuildandbenchmarkexample"></a>CONFIGURE, BUILD AND BENCHMARK EXAMPLE ##
--------------------------------------------------------------------------------

 The performance-critical parts of librsb are written in C and C++.
 Compilation flags for C and C++ are specified by the CFLAGS and CXXFLAGS
 variables
 Fortran flags can be set by FCFLAGS, but are of minor importance.
 A sensible setup can be:

	./configure CC=gcc CFLAGS='-Ofast' CXX=g++ CXXFLAGS='-Ofast -std=c++17' FC=gfortran --prefix=/opt/librsb

 You should be able to specify other combination of compilers, like icx+icpx+ifc
 or icc+icpc+ifort. If you have luck, even mix compilers from different suites.

 If you wish to build the librsb benchmark program with the MKL library (say,
 for performance comparison purposes), you can use one of the following three
 recipes as a base.

 With gcc, 64 bit (and a few more options you may want to adapt):

	export MKLROOT=/opt/intel/mkl
	./configure \
	  CC=gcc CFLAGS='-Ofast' \
	  CXX=g++ CXXFLAGS='-Ofast -std=c++17'         \
          FC=gfortran \
	  --with-mkl="-static -L${MKLROOT}/lib/intel64 \
	  -Wl,--start-group,-lmkl_intel_lp64,-lmkl_gnu_thread,-lmkl_core,--end-group \
	  -fopenmp -lpthread"                        \
	  --with-memhinfo=L2:4/64/512K,L1:8/64/24K   \
	  --with-mkl-include=/opt/intel/mkl/include/

 With icc, 64 bit)

	export MKLROOT=/opt/intel/mkl
	./configure \
	  CC=icc CFLAGS='-Ofast' \
	  CXX=icpc CXXFLAGS='-Ofast -std=c++17'    \
          FC=ifort \
	--with-mkl="-static -L${MKLROOT}/lib/intel64 -openmp -lpthread \
	-Wl,--start-group,-lmkl_intel_lp64,-lmkl_intel_thread,-lmkl_core,--end-group" \
	--with-memhinfo=L2:4/64/512K,L1:8/64/24K   \
	--with-mkl-include=/opt/intel/mkl/include/

  or 32 bit:

	./configure \
	  CC=gcc  CFLAGS='-Ofast'
	  CXX=g++ CXXFLAGS='-Ofast -std=c++17'        \
          FC=gfortran \
	 --with-memhinfo=L2:4/64/512K,L1:8/64/24K     \
	 --with-mkl="-static -L/opt/intel/mkl/lib/ia32/ -lmkl_solver \
	 -Wl,--start-group,-lmkl_intel,-lmkl_gnu_thread,-lmkl_core,--end-group \
	 -fopenmp -lpthread" \
	 --with-mkl-include=/opt/intel/mkl/include/

Once you chose your configure options, you want to build.

	make cleanall # (optional) delete generated sources
	make          # build library, documentation, tests, examples
	make qtests   # optional
	make install  # optional
	make itests   # optional

If you modified configure options impacting code generation, may want to run
`make cleanall`, which deletes stale sources and ensures sources re-generation
before the 'make'.

To speed up 'make', consider using the parallelism option, e.g.:

	make -j 3     # use 3 build threads

 After the build, say you want to benchmark the library for SpMV/SpMM.
 You have a Matrix Market file representing a matrix, `A.mtx`,
 you want to use 1 and 4 cores, and type Z (double complex).
 With SpMV, and SpMM with 2 right hand sides laid in C order.
 Then running:

	./rsbench -oa -Ob --bench -f A.mtx -qH -R -n1,4 -T z --verbose --nrhs 1,2 --by-rows

 will output performance and timing results in a tabular form.

 If not specifying a type (argument to the `-T` option), all will be used.
 If configured in at build time, choices may be `-T D` (where `D` is the BLAS
 prefix for `double`), `-T Z` (`Z` stands for `double complex`) and so on.
 You can also pass `-T :` to specify all of the configured types.

 For more options and configuration information, invoke:

	./rsbench -oa -Ob --help
	./rsbench --help
	./rsbench --version
	./rsbench -I
	./rsbench -C

 An example Matrix Market matrix file contents:

	%%MatrixMarket matrix coordinate pattern general
	% This is a comment.
	% See other examples in the distributed *.mtx files.
	2 2 3
	1 1
	2 1
	2 2


## <a id="compatibility"></a>COMPATIBILITY ##
--------------------------------------------------------------------------------

 This library was developed mostly on Debian Linux and using only free software.
 
 This library has been built and tested on Unix machines.
 Microsoft Windows users can try building librsb under the Cygwin environment.

 (Ancient comment kept for the nostalgic reader)
 Some tricks may have to be used on IBM AIX. For instance, adding the
 `--without-xdr` or the `--without-zlib` switch to `./configure`.
 Your mileage may vary.
 AIX's `make` program may give problems; use the GNU version `gmake` instead;
 the same shall be done with the M4 interpreter.


## <a id="faq"></a>FAQ ##
--------------------------------------------------------------------------------

 Q: **Can you provide me good configure defaults for an optimized build ?**

 A: Default `./configure` options are appropriate for an optimized build.
    A good starting point for `gcc` is `./configure CC=gcc CFLAGS='-O3'`. 
    However, if you need complex arithmetic and are using GCC, I'd advise using
    -Ofast. On many versions of GCC I observed sub-optimal complex arithmetic
    performance with `-O3`, regardless use of e.g. `-mtune=native`.
    For more, consult your compiler documentation (e.g. `man gcc`, `man icx`, 
    `man icc`),
    and learn about the best flags for your specific platform.
    Striping your executable (`make install-strip` for librsb's `rsbench`) may
    help.


 Q: **I am a beginner and I wish librsb to be very verbose when I invoke
    library interface functions incorrectly.
    Can you provide me good configure defaults for such a "debug" build ?**

 A: Yes: via the `./scripts/configure_for_debug.sh` script.


 Q: **I have problems linking, seems like some Fortran library is missing. 
    What should I do?**

 A: Did you try `./configure --enable-fortran-linker` ?


 Q: **I have machine X, compiler Y, compiling flags Z; is SpMV performance P 
    with matrix M good ?**

 A: In general, hard to tell. You can `make hinfo.log` and send to me (see
    CONTACTS) the `hinfo.log` file and your matrix in Matrix Market format
    (well, please don't send matrices by email but rather upload them
    somewhere on the web and send an URL to them).
    The `hinfo.log` file will contain useful compile and machine information.
    Then I *may* get an idea about the performance you should get with that
    matrix on that computer.


 Q: **What is the Sparse BLAS ?**

 A: It's a programming interface specification:

 * [sparseblas_2001]:
   BLAS Technical Forum Standard, Chapter 3, Sparse BLAS.
   http://www.netlib.org/blas/blast-forum/chapter3.pdf

 * [dhp_2002]:
   An Overview of the Sparse Basic Linear Algebra Subprograms:
    The New Standard from the BLAS Technical Forum.
   IAIN S. DUFF, CERFACS and Rutherford Appleton Laboratory.
   MICHAEL A. HEROUX, Sandia National Laboratories.
   ROLDAN POZO, National Institute of Standards and Technology.

 * [dv_2002]:
   Algorithm 818:
    A Reference Model Implementation of the Sparse BLAS in Fortran 95.
   IAIN S. DUFF, CERFACS, France and Atlas Centre, RAL, England.
   CHRISTOF VÖMEL, CERFACS, France.


 Q: **Is there an easy way to profile librsb usage in my application ?**

 A: Yes: build with `--enable-librsb-stats` and extract time elapsed in librsb
    via e.g.:`rsb_lib_get_opt(RSB_IO_WANT_LIBRSB_ETIME,&etime);`.


 Q: **Why another sparse matrix library ?**

 A: This library is the fruit of the author's PhD work, focused on researching
    improved multi threaded and cache friendly matrix storage schemes for e.g.
    PSBLAS.


 Q: **What are the key features of this library when compared to other ones ?**

 A: Recursive storage, a code generator, parallel BLAS operations
    (including matrix assembly, matrix-matrix multiplication, transposed
     matrix-vector multiply), a rich test suite, a Sparse BLAS
     interface and a free software licensing.
 

 Q: **How do I detect librsb from my package's `configure` script ?**

 A: Please check out `examples/configure.ac` and `examples/makefile.am`:


 Q: **How is correctness checked in the librsb test suite ?**

 A: In different ways, often configuration dependent.
    Check out the `qtests` and `tests` targets in the Makefile.

 Q: **Why did you originally write the library in C and not in C++ ?**

 A: C is pretty easy to interface with C++ and Fortran.
    Using a debugger in C++ can be a headache.
    Also C's `restrict` keyword.
    

 Q: **Why did you use C and not Fortran ?**

 A: This library is slightly system-oriented, and system calls interfacing is
    much easier in C. Also C's pointer arithmetic was desirable.


 Q: **Is there a quick and easy way to perform an artificial performance
    test with huge matrices without having to program ?**

 A: Express your matrix in a Matrix Market format as the `A.mtx` file and
    use it as e.g.:**

	./rsbench -oa -Ob --bench -f A.mtx --verbose --nrhs 1,4 --by-rows


 Q: **I've found a bug! What should I do ?**

 A: First please make sure it is really a bug:
    read the documentation, check, double check.
    Then you can write a description of the problem, with a minimal program
    source code and data to replicate it.
    Then you can jump to the CONTACTS details section.


 Q: **Is it possible to build matrices of, say, `long double` or
    `long double complex` or `int` or `short int` ?**

 A: Yes. Invoke the configure script accordingly, e.g.:
      `--enable-matrix-types="long double"` or
      `--enable-matrix-types="double,long double,double complex"`
    If this breaks code compilation, feel free to contact the author
    (see the CONTACTS section).
    Note that some combinations may break `make qtests` or other test recipes.


 Q: **Is there a way to compare the performance of this library to some other
    high performance libraries ?**

 A: If you build `rsbench` with support for the Intel MKL library, then you
    can do performance comparisons with e.g.:
    `# ./rsbench -oa -Ob -qH -R --gen-diag 100 --compare-competitors --verbose`
    or use the following script:
    `# bench/dense.sh ' '`
    Or even better, check out the `--write-performance-record` feature.
    For details see the output of:
    `# rsbench -oa -Ob --help`


 Q: **Is there a non-threaded (serial) version of librsb ?**

 A: Yes: you can configure the library to work serially (with no OpenMP).
    See `./configure --help`. 


 Q: **Is this library thread-safe ?**

 A: Probably yes: no static buffers are being used, and reentrant C standard
    library functions are invoked.


 Q: **Does the librsb library run on GPUs?**

 A: Not yet.


 Q: **I built and compiled the code without enabling any BLAS type (S,D,C,Z), 
     and both `make qtests` and `make tests` ran successfully outside the
     `./examples` directory, but `make tests` breaks within `./examples` 
     directory.**

 A: Well, the tests passed because the examples testing was simply skipped.
    The example programs need at least one of these types to work.


 Q: **At build time I get many "unused variable" warnings. Why ?**

 A: librsb accommodates many code generation and build time configuration
    options. Some combinations may turn off compilation of certain parts of the
    code, leading some variables to be unused.


 Q: **Are there papers to read about the RSB format and algorithms ?**

 A: Yes, the following:

 *  Michele Martone, Simone Bacchio.
    Portable performance on multi-threaded Sparse BLAS operations with PyRSB
    Proceedings of SciPy 2021 (Scientific Python Conference).
    https://doi.org/10.25080/majora-1b6fd038-00e

 *  Michele Martone.
    Efficient Multithreaded Untransposed, Transposed or Symmetric Sparse
    Matrix-Vector Multiplication with the Recursive Sparse Blocks Format.
    Parallel Computing 40(7): 251-270 (2014).
    http://dx.doi.org/10.1016/j.parco.2014.03.008

 *  Michele Martone.
    Cache and Energy Efficiency of Sparse Matrix-Vector Multiplication for
    Different BLAS Numerical Types with the RSB Format.
    Proceedings of the ParCo 2013 conference, September 2013, Munich, Germany.
    PARCO 2013: 193-202.
    http://dx.doi.org/10.3233/978-1-61499-381-0-193

 *  Michele Martone, Marcin Paprzycki, Salvatore Filippone.
    An Improved Sparse Matrix-Vector Multiply Based on Recursive Sparse Blocks  
    Layout.
    LSSC 2011: 606-613.
    http://dx.doi.org/10.1007/978-3-642-29843-1_69

 *  Michele Martone, Salvatore Filippone, Salvatore Tucci, Marcin Paprzycki,
    Maria Ganzha.
    Utilizing Recursive Storage in Sparse Matrix-Vector.
    Multiplication - Preliminary Considerations. CATA 2010: 300-305
    
 *  Michele Martone, Salvatore Filippone, Marcin Paprzycki, Salvatore Tucci.
    Assembling Recursively Stored Sparse Matrices. IMCSIT 2010: 317-325.
    http://www.proceedings2010.imcsit.org/pliks/205.pdf

 *  Michele Martone, Salvatore Filippone, Pawel Gepner, Marcin Paprzycki,
    Salvatore Tucci.
    Use of Hybrid Recursive CSR/COO Data Structures in Sparse Matrices-Vector
    Multiplication.
    IMCSIT 2010: 327-335.
    http://dx.doi.org/10.1109/SYNASC.2010.72

 *  Michele Martone, Salvatore Filippone, Marcin Paprzycki, Salvatore Tucci.
    On BLAS Operations with Recursively Stored Sparse Matrices.
    SYNASC 2010: 49-56.
    http://dx.doi.org/10.1109/SYNASC.2010.72

 *  Michele Martone, Salvatore Filippone, Marcin Paprzycki, Salvatore Tucci.
    On the Usage of 16 Bit Indices in Recursively Stored Sparse Matrices.
    SYNASC 2010: 57-64.
    http://dx.doi.org/10.1109/SYNASC.2010.77


 Q: **I have make-related problems on my FreeBSD/OpenBSD/NetBSD.
    What should I do ?**

 A: Did you try using GNU Make implementation via
    e.g.: `./configure MAKE=gmake ...`?


 Q: **Running `make` keep crashing when recreating rsb.F90!
      And I do not need Fortran..**

 A: Ok, disable it: `./configure FC=' ' ...`.


 Q: **I keep getting C++ compile or link problems.
     How do I do?**

 A: Use: `./configure CXX=' ' ...`.
    This turns off rsb_spmm/usmm faster kernels, and several test programs.


 Q: **I have M4-related problems on IBM SP5/SP6/... . What should I do ?**

 A: A fix is to use a GNU M4 implementation 
    e.g.: `./configure M4=/opt/freeware/bin/m4 ...` (or `M4=gm4`)
    BTW, if it's an old POWER machine consider donating it to a computer museum,
    like https://museo.freaknet.org/
   

 Q: **Can I skip the OpenMP flags guessing mechanism and specify flags my own?**

 A: Yes, e.g.: `--enable-openmp 'ac_cv_prog_c_openmp=-DMY_OMP_FLAGS'`
    If you use e.g. `OPENMP_CFLAGS=-DMY_OMP_FLAGS` it won't skip the tests, and
    may introduce unwanted flags.


## <a id="possiblefutureenhancements"></a>POSSIBLE FUTURE ENHANCEMENTS ##
--------------------------------------------------------------------------------

 * auxiliary functions for numerical vectors
 * CSC,BCSR,BCSC and other block-level formats
 * performance prediction/estimation facilities (experimental)
 * types of the blocks, nonzeroes, and coordinates indices can be user specified
 * automatic matrix blocking selection (for BCSR/BCSC) 
 * an arbitrary subset of block size kernels can be specified to be generated
 * recursive storage variants of blocked formats (non uniform blocking)
 * more auto-tuning and prediction control


## <a id="abouttheinternals"></a>ABOUT THE INTERNALS ##
--------------------------------------------------------------------------------

 The following good practices are being followed during development of librsb:

 - only symbols beginning with `rsb_` or `blas_` or `BLAS_` are being exported
 - internal functions are usually prefixed by `rsb__`
 - no library internal function is expected to call any API function


## <a id="bugs"></a>BUGS ##
--------------------------------------------------------------------------------

 If you encounter any bug (e.g.: mismatch of library/program behaviour and
 documentation, please let me know about it by sending me (see CONTACTS) all
 relevant information (code snippet, originating data/matrix, `config.log` ), in
 such a way that I can replicate the bug behaviour on my machines.
 If the bug occurred when using rsb interfaced to some proprietary library,
 please make sure the bug is in librsb.

 It may be of great help to you to build the library with the debug compile
 options on (e.g.: `CFLAGS='-O0 -ggdb'`), and with appropriate library verbosity
 levels, e.g. (`--enable-internals-error-verbosity=1` and
  `--enable-interface-error-verbosity=1` options to configure) to better 
 understand the program behaviour before sending a report.

 Make sure you have the latest version of the library when reporting a bug. 


## <a id="contacts"></a>CONTACTS ##
--------------------------------------------------------------------------------

 You are welcome to contact the librsb author:

  Michele Martone `<michelemartone AT users DOT sourceforge DOT net>`
 
 Please specify "librsb" in the "Subject:" line of your emails.

 More information and downloads on  http://sourceforge.net/projects/librsb

 Mailing list: https://lists.sourceforge.net/lists/listinfo/librsb-users
 

## <a id="credits"></a>CREDITS	(in alphabetical order) ##
--------------------------------------------------------------------------------

For librsb-1.3:

 * Juan Durillo Barrionuevo patiently tested and helped smoothing the build
   system on the Mac OS.
 * Salvatore Cielo provided testing on the Mac OS.
 * Rafael Laboissiere for continued help in aspects of build, documentation,
   testing and distribution.
 * Markus Muetzel suggested use of putenv() as Windows alternative to setenv().
 * Nisarg Patel provided testing on the Mac OS.
 * Matthieu A. Simonin for testing, discussion and suggestions.

For librsb-1.2:

 * Marco Atzeri provided testing, patches to build librsb under cygwin over
   nearly each release, and spotted a few bugs.
 * Fabio Cassini spotted an unintended conversion via sparsersb and +.
 * John Donoghue spotted a rendering corner case bug.
 * Sebastian Koenig spotted a computational bug in -rc6.
 * Rafael Laboissiere helped a lot improving the documentation and the build 
   system.
 * Mu-Chu Lee provided a patch to fix sorting code crashing with > 10^9 nnz.
 * Constanza Manassero spotted an inconsistency in the usmm/ussm interface.
 * Markus Muetzel helped debugging rsb_mtx_rndr().
 * Dmitri Sergatskov spotted a double free in rsb_mtx_rndr() and convinced about
   the necessity of sanitizing memory usage.

For librsb-1.1:

 * Gilles Gouaillardet provided a patch for OpenMP-encapsulated I/O.
 * Marco Restelli provided with testing and detailed comments and suggestions.

For librsb-1.0:

 * Francis Casson helped with testing and documentation reviewing during 
   the first release.
 * Nitya Hariharan helped revising early versions of the documentation.


ACKNOWLEDGEMENTS
--------------------------------------------------------------------------------

  In 2019-2022, librsb has been developed under the PRACE-6IP-WP8 (GRANT
  AGREEMENT NUMBER 823767 -- PRACE-6IP) EU project LyNcs (partners:
  Computation-based Science and Technology Research Centre (CaSToRC) of The
  Cyprus Institute and Inria).


## <a id="license"></a>LICENSE ##
--------------------------------------------------------------------------------

 This software is distributed under the terms of the Lesser GNU Public License
 version 3 (LGPLv3) or later.
 See the COPYING file for a copy of the LGPLv3.

 librsb is free software.
 To support it, consider writing to the author and acknowledging use of librsb 
 in your publications.
 All that would be very appreciated.


--------------------------------------------------------------------------------
