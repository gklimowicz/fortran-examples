(* circe2/thoMatrix.mli --  *)
(* Copyright (C) 2001-2014 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
   Circe2 is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by 
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
   Circe2 is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)  

val copy : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  'a array array -> 'a array array
val map : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  ('a -> 'b) -> 'a array array -> 'b array array
val iter : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  ('a -> unit) -> 'a array array -> unit
val fold_left : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  ('a -> 'b -> 'a) -> 'a -> 'b array array -> 'a
val sum_float :  ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  float array array -> float
val size : 'a array array -> int

val transpose : 'a array array -> 'a array array

val suite : OUnit.test

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
