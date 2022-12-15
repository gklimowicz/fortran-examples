title: Overview

*Avanced Options in Tables and Universal Scripting*

![Aotus](logo.svg)

This is a Fortran wrapper for the [Lua](http://www.lua.org) scripting language,
dedicated to enable flexible configuration of Fortran applications with this
full fledged scripting language.
Aotus is also the [Night monkey](http://en.wikipedia.org/wiki/Night_monkey)
(living in south america).
Thus, it can be understood to be interacting with the moon (Lua, provided by
Pontifical Catholic University of Rio de Janeiro in Brazil) and providing the
means to retrieve data from its scripts.

![Powered by Lua](powered-by-lua.gif)

The most prominent data structure in Lua are tables, which provide the
possibility to store complex data structures.
Thus, the configuration is usually done by global variables in the Lua script or
tables.

Aotus provides several layers, encapsulating the bare Lua C-API:

 * [[lua_fif]]: this just provides the ISO_C_BINDING interface declarations to
   the Lua API
 * [[flu_binding]]: this is the actual Fortran binding wrapped around lua_fif,
   to provide a more Fortran like interface. Especially the **flu_state** type
   is declared, which maintains the handle for the Lua context.
 * [[aot_table_module]]: provides some convenience functions to work on Lua
   tables in Fortran
 * [[aot_fun_module]]: provides some convenience functions to work with Lua
   functions in Fortran
 * [[aotus_module]]: provides the high end level to easily retrieve data from a
   Lua script

On top of those there is an additional [[aot_vector_module]] that allows the
direct reading of values into arrays of rank 1.

Finally, there is an additional module which allows [output](Output.html) of
Fortran values into nested Lua tables.

The library can be compiled by various modern Fortran compilers as outlined in
[Compiler Support](CompilerSupport.html).

**An example, showing the usage of this library in a Fortran application, is
given in aotus_sample.f90. You can find it in the *build* directory. The
corresponding Lua script is given in config.lua, which you can find in the
*sample* directory**

Note on usage in parallel environments: Aotus itself is not providing parallel
facilities. But it can be nicely used in parallel aswell. However, for
massively parallel systems, it is advisable to minimize the access to config
files. To avoid excessive filesystem meta accesses it is recommended to
load required files only on one process.
An implementation of this for MPI can be found in TreElMs
[distconf](https://geb.inf.tu-dresden.de/doxy/treelm/proc/tem_open_distconf_array.html).


Reading a Lua Script
--------------------

You need a handle for the Lua context of type [[flu_state]]. You can get that by
opening and processing a Lua script with [[open_config_file]]:

    :::fortran
    call open_config_file(L, filename, errCode, errString)

Instead of reading a script from a file, you can also execute a string directly
by using [[open_config_chunk]] with the same interface, but replacing the
filename with the chunk of Lua code to be executed.
The arguments errCode and errString are optional and return errors, which might
occur while loading or executing the Lua code.

It is also possible to load already processed scripts in byte code by using the
[[open_config_buffer]] routine, which expects an array of characters with the
script in byte code.

In the end, after getting all configuration values, close it again with
[[close_config]]:

    :::fortran
    call close_config(L)

Retrieving Variables from the Script
------------------------------------

From configuration files you usually want to obtain some parameters to steer
your application. The most important interface for this functionality is
[[aot_get_val]], it is a generic interface to different functions, which allow
you to obtain global values, values within tables and complete vectors as
tables. Sometimes, especially in the context of evaluating functions, you might
also need [[aot_top_get_val]], which always tries to obtain the topmost value
on the Lua stack.

The [[aot_get_val]] interface is provided by the [[aotus_module]] and
generically looks like this:

    :::fortran
    call aot_get_val(val, errCode, L, key, default)

Where:

* **val**: is the variable the value should be stored in, all intrinsic
  data types should be supported, reals with kind single and double precision
* **errCode**: Returns an error code with various bits set for different errors,
  which might happen while retrieving the variable.
  They can be checked by btest, and the different error codes are encoded in
  parameters:
    * **aoterr_fatal**: Something went irrecoverably wrong
    * **aoterr_nonExistent**: The requested variable is not set in the Lua script
    * **aoterr_wrongType**: The requested variable in the Lua script does not meet
      the requested data type
    * For example you can check for a fatal error by using
      *btest(errCode, aoterr_fatal)*
* **L**: is the Lua context of type flu_state
* **key**: is a string identifying the variable you want to retrieve
* _optional_ **default**: A default value to put into the variable, if the
  variable is not provided in the Lua script

In general we get the following shape for the interface:

    :::fortran
    call aot_{top}_get_val(<outputs>, <id>, default)

Where *outputs* is: **val** and **errCode** and *id* is at least the Lua
context (**L**) for the *aot_top* variant. For global variables there has
to be a **key** in the *id* and for tables there has to be a **thandle**.
In tables the **key** might be replaced by a **pos** argument.

### Tables

The interface to work with tables is trying to resemble IO, thus you could think
of a table as a file, which you can open and read values out of by
referencing its unit (handle).
Opening and closing tables is provided by the [[aot_table_module]].

To work with a table, you first need to get a handle to identify the table.
For globally defined tables this can be done by using

    :::fortran
    call aot_table_open(L, thandle, key)

Where:

* **L**: is the Lua context of type flu_state
* **thandle**: a handle to reference this table
* **key**: is the name of the globally defined table to retrieve

For a table within an already opened table use:

    :::fortran
    call aot_table_open(L, parent, thandle, key, pos)

Where the additional arguments are:

* **parent**: the handle of the table, which this table should be looked up in
* _optional_ **pos**: Referring to the table to retrieve by position instead of
  name, it is optional as well as the **key**, and one of them has to be present

The handle will be 0, if the variable does not exist, or is not a table.
After you have the handle to the table, you can access its components with

    :::fortran
    call aot_table_get_val(val, errCode, L, thandle, key, pos, default)

Which is essentially the same interface as for global variables, except for
the optional argument **pos**, by which the unnamed entries in the table are
accessible by their position and the handle to the table, where the component
is to be looked up.
Both **pos** and **key** are optional, providing the ability to access the
variables either purely by their order, or their name.
If both are provided, the key takes precedence.
The handling of positional or named addressing is a little bit similar to
the Fortran convention, that is, as soon as there is a named component in the
table, all following components should also be named. Positional references are
only valid up to this position.

After all values are read from the table, the table should be closed again by
calling

    :::fortran
    call aot_table_close(L, thandle)

### Type dependent actions

Sometimes there might be different types possible for a given setting, and
different actions need to be taken for each possible type.
In this case you might first want to check the type of a given variable
before proceeding and reading the actual value.
This can be achieve by [[aot_type_of]], which is a function that will
put the requested variable onto the top of the stack and return the
Lua data type of it:

    :::fortran
    luatype = aot_type_of(L, thandle, key, pos)

With:

* **L**: is the Lua context
* **thandle**: handle of the parent table
* **key**: name of the Lua variable to get to the top of the stack
* **pos**: positional addressing of the variable inside thandle

Afterwards, the actual value can then be read by the corrsponding
[[aot_top_get_val]].

### Functions

Again functions try to resemble the usage of files, however in this case its
slightly more complicated, as you first need to "write" the input parameters
into the function, then execute it and finally retrieve the results.

To use a function, that is globally defined, open it with:

    :::fortran
    call aot_fun_open(L, fun, key)

Where:

* **L**: is the Lua context
* **fun**: is the handle to the opened function
* **key**: is the name of the function you want to access

To access a function, which is within a table, use:

    :::fortran
    call aot_fun_open(L, parent, fun, key, pos)

Where the additional arguments are:

* **parent**: the table handle of the table, the function should be looked up in
* _optional_ **pos**: Refer to the function by position instead of name (**key**
  is also optional)

After the function is opened, its arguments need to be filled with:

    :::fortran
    call aot_fun_put(L, fun, arg)

Where:

* **L**: Lua context
* **fun**: handle to the function to put the arguments in
* **arg**: argument to provide to the function

When all arguments are written, the function needs to be executed with:

    :::fortran
    call aot_fun_do(L, fun, nresults)

Where:

* **L**: Lua context
* **fun**: opened function to execute
* **nresults**: number of results you want to retrieve from that function

After the function is executed, the results can be read, using:

    :::fortran
    call aot_top_get_val(val, errCode, L, default)

Where:

* **val**: value to return
* **ErrCode**: an error code, if something went wrong
* **L**: Lua context
* _optional_ **default**: a default value to use, if no value can be retrieved

You will get the results in reversed order if there are multiple results.
That is, because the first call to [[aot_top_get_val]] will return the last result returned
by the function, the next the second last, and so on.

You may then go on and put new arguments into the function, execute
it and retrieve the corresponding results.

After you are done with the evaluation of the function it has to be closed
with:

    :::fortran
    call aot_fun_close(L, fun)

This should cover the most typical tasks for configurations.


### Vectors

In addition to the scalar retrieval routines, there is an aot_vector_module
provided, which provides interfaces with arrays of rank 1 to access vectorial
data in form of tables.
They follow the same interfaces, as the scalar routines, however the values,
error codes and defaults have to be one-dimensional arrays.

    :::fortran
    call aot_get_val(val, errCode, L, key, default)

Where:

* **L**: is the Lua context of type flu_state
* **key**: is a string identifying the variable you want to retrieve
* **val(:)**: is the array the vector should be stored in
* **errCode(:)**: Returns an error code with various bits set for different errors,
  which might happen, while retrieving each component of the vector.
  They can be checked by btest, and the different error codes are encoded in
  parameters:
    * **aoterr_fatal**: Something went irrecoverably wrong
    * **aoterr_nonExistent**: The requested component is not set in the Lua script
    * **aoterr_wrongType**: The requested component in the Lua script does not
      meet the requested data type
    * For example you can check for a fatal error by using
      "btest(errCode(1), aoterr_fatal)"
* _optional_ **default(:)**: A default vector to put into the variable, if the
  variable is not provided in the Lua script, it also fills vectors, which are
  only partially defined in the configuration script.

The interface for vectors within other tables is defined accordingly.
In the interface described above, the conf_val vector has a given size and
just the values need to be filled.
However, it might be necessary to retrieve arrays, of which the size is not
known beforehand, and should depend on the table definition in the configuration.
For these cases there are additional routines defined, which take allocatable
arrays as input.
These are then allocated and filled according to the configuration or the
provided default vector.
In these interfaces you have to provide an additional parameter **maxlength**
that limits the size of the vector to allocate.
If the Lua variable in question is not a table but rather a scalar of the
correct type, an array of length 1 will be created for this single value.
For undefined vectors zero sized arrays are returned for the values and the
error codes.

Using Fortran Code in Lua Scripts
---------------------------------

Aotus also provides the interface to register functions in Lua to allow the
usage of your Fortran implementations in Lua.
A short description for this feature is provided in
[Using Fortran in Lua](FortranInLua.html).
