title: Output

The [[aot_out_module]] provides some convenience functions to
output data into nested Lua tables.
It purely relies on Fortran formatted IO, and therefore could also be used
without the Lua library.
Despite its high independence from the rest of the library it can be considered
as a counter part to the reading functions.
That is Aotus allows you to use Lua as a closed system, where you store values
from one Fortran application in a script and read it back in another Fortran
application using Lua.
This yields a great flexibility and allows you to design easily extendable and
portable file exchanges.

Writing a Lua Script
--------------------

The output mechanism relies on a file-handle described by the **aot_out_type**.
It is initialized with:

    :::fortran
    call aot_out_open(put_conf, filename, outUnit, indentation)

Where:

 * **put_conf**: is the handle to the opened output file
 * _optional_ **filename**: name of the file to write to
 * _optional_ **outUnit**: Pre-Connected file unit to write to
 * _optional_ **indentation**: Number of spaces to use for indentation
   (defaults to 4).

You either have to state the **filename** or the **outUnit**, if both are given,
the filename will be used and outUnit ignored.
When you specify a filename to open, the according file will be replaced, if it
exists and opened for writing.

After putting all values into the script, the output has to be finished with

    :::fortran
    call aot_out_close(put_conf)

This will close the given handle and, if the file was opened it is closed as
well.

The concept of writing into the Lua script is very similar to the reading.
Values are written by calling [[aot_out_val]], which is a generic interface for
the various intrinsic types in Fortran, except complex numbers.
The interface for scalar values looks like this:

    :::fortran
    call aot_out_val(put_conf, val, vname, advance_previous)

Where:

 * **put_conf**: is the handle to the output file
 * **val**: value to write
 * _optional_ **vname**: A name to assign to this value.
 * _optional_ **advance_previous**: Flag if the value should be put on a new
   line.

Arrays will be put into tables internally, and the interface looks slightly
different:

    :::fortran
    call aot_out_val(put_conf, val, vname, advance_previous, max_per_line)

Where all parameters have the same meaning as in the scalar interface, and the
additional optional parameter **max_per_line** indicates how many entries of
the array should be put on a common line. If none is provided some default
depending on the data type will be used.
The opening bracket of the table will always be on the same line as the first
entry, and the closing bracket on the same line as the last entry.


Writing to Tables
-----------------

To put values into arbitrarily nested tables, there are opening and closing
calls for the tables, all values in between these two calls will be put into
this table.
Again this handling of tables is similar to the interface used in the input
routines.
Opening a table is done with:

    :::fortran
    call aot_out_open_table(put_conf, tname, advance_previous)

Where:

 * **put_conf**: is the handle to the script to write into
 * _optional_ **tname**: Optionally a name might be assigned to the table
 * _optional_ **advance_previous**: A flag to indicate if this table should be
   started on a new line (default = .true.)

To close the table again the following call has to be used:

    :::fortran
    call aot_out_close_table(put_conf, advance_previous)

Where:

 * **put_conf**: is the handle to the script to write into
 * _optional_ **advance_previous**: A flag to indicate the closing bracket
   should be put on a new line (default = .true.)

Example
-------

Here is a short example how this could be used in Fortran:

    :::fortran
    program aot_out_test
      use aot_out_module
    
      implicit none
    
      type(aot_out_type) :: dummyOut
    
    
      call aot_out_open(put_conf = dummyOut, filename = 'dummy.lua')
    
      call aot_out_open_table(dummyOut, 'screen')
      ! Screen table
      call aot_out_val(dummyOut, 123, 'width')
      call aot_out_val(dummyOut, 456, 'height')
    
      call aot_out_val(dummyOut, [100.0, 0.0], vname='origin')
      ! End of screen table
      call aot_out_close_table(dummyOut)
    
      call aot_out_val(dummyOut, [0, 1, 2, 3], vname='testarray')
    
      call aot_out_close(dummyOut)
    end program aot_out_test

And here is the resulting Lua script:

    :::lua
    screen = {
        width = 123,
        height = 456,
        origin = { 100.00000000, 0.00000000 } 
    }
    testarray = { 0, 1, 2, 3 }
