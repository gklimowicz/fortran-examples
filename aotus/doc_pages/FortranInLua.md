title: Using Fortran in Lua

The [[flu_binding]] also includes the interfaces to functions, required to
incorporate your Fortran code in Lua scripts.
A related project to this approach is
[Numeric Lua](http://numlua.luaforge.net/), that exposes standard numeric
libraries to Lua.
If you intend to use your Fortran implementations in such a way, embedding them
in a Lua script you can do this by using Aotus and compiling your Fortran code
together with it into a library.

To register Fortran functions for usage in Lua, you need to define a
*luaopen_libname* function and call [[flu_register]], for each function to
expose to Lua.

Here is an example for this setup:

    :::fortran
    module libtest
    
      use, intrinsic :: iso_c_binding
      use flu_binding
    
      implicit none
    
    contains
    
      function luaopen_libtest(lua_state) result(val) bind(c)
    
          type(c_ptr), value :: lua_state
          integer(c_int) :: val
          type(flu_state) :: fL
    
          fL = flu_copyptr(lua_state)
    
          call flu_register(fL, "hello_f90", hello_f90)
    
          val = 0
    
      end function luaopen_libtest
    
      function hello_f90(lua_state) result(val) bind(c)
    
          use aotus_module
    
          type(c_ptr), value :: lua_state
          integer(c_int) :: val
          integer :: x, err
          type(flu_state) :: fL
    
          fL = flu_copyptr(lua_state)
          write (6,*) "HELLO FROM FORTRAN! :-)"
          call flu_pushinteger(fL, 102)
          call flu_setglobal(fL, "af90")
          val = 0
    
      end function hello_f90
    end module libtest

After compilation of this code it could be used in Lua like this:

    :::sh
    $ lua -e "require('libtest'); hello_f90(); print(af90)"

Resulting in:

    HELLO FROM FORTRAN! :-)
    102
