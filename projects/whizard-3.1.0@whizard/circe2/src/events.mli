(* circe2/events.mli --  *)
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

(* We're dealing with Fortran style \texttt{DOUBLE PRECISION}
   arrays exclusively. *)

type t =
    (float, Bigarray.float64_elt, Bigarray.fortran_layout) Bigarray.Array2.t

(* Read an ASCII representation of a big array from a channel or a file.
   The array is read in pieces of [chunk] columns each; the default value
   for [chunk] is 100000.  The number of rows is given by the integer
   argument, while the number of columns is determined by the number of
   lines in the file.  If the [file] argument is present the resulting
   bigarray is mapped to a file.  *) 
val of_ascii_channel : ?file:string -> ?chunk:int -> int -> in_channel -> t
val of_ascii_file : ?file:string -> ?chunk:int -> int -> string -> t

(* Map a file containing a binary representation of a big array.  The
   number of rows is again given by the argument and the number of
   columns is determined by the size of the file.  The first version
   does a read-only (or rather copy-on-write) map, while the second
   version allows modifications. *)
val of_binary_file : int -> string -> t
val shared_map_binary_file : int -> string -> t

(* Selfexplaining, hopefully \ldots *)
val to_ascii_channel : out_channel -> t -> unit
val to_ascii_file : string -> t -> unit
val to_binary_file : string -> t -> unit

(* Rescale the entries. *)
val rescale : float -> float -> t -> unit

(* Utilities for reading ASCII representations. *)
val next_float : Lexing.lexbuf -> float

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
