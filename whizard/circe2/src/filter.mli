(* circe2/filter.mli --  *)
(* Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
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

type t
val unit : t
val gaussian : float -> t
val apply : ?inf:int -> ?sup:int -> t -> float array -> float array
val apply1 : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  t -> float array array -> float array array
val apply2 : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  t -> float array array -> float array array
val apply12 : ?inf1:int -> ?sup1:int -> ?inf2:int -> ?sup2:int ->
  t -> t -> float array array -> float array array

exception Out_of_bounds of int * int 

val suite : OUnit.test

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)


