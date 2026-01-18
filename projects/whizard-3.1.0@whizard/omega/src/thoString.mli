(* thoString.mli --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>

   WHIZARD is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   WHIZARD is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)

(* This is a very simple library if stroing manipulation functions missing
   in O'Caml's standard library. *)

(* [strip_prefix prefix string] returns [string] with 0 or 1
   occurences of a leading [prefix] removed. *)
val strip_prefix : string -> string -> string

(* [strip_prefix_star prefix string] returns [string] with any number
   of leading occurences of [prefix] removed. *)
val strip_prefix_star : char -> string -> string

(* [strip_prefix prefix string] returns [string] with a leading
   [prefix] removed, raises [Invalid_argument] if there's no match. *)
val strip_required_prefix : string -> string -> string

(* [strip_from_first c s] returns [s] with everything starting from
   the first [c] removed.  [strip_from_last c s] returns [s] with
   everything starting from the last [c] removed. *)
val strip_from_first : char -> string -> string
val strip_from_last : char -> string -> string

(* [index_string pattern string] returns the index of the first
   occurence of [pattern] in [string], if any.  Raises [Not_found], if
   [pattern] is not in [string]. *)
val index_string : string -> string -> int

(* This silently fails if the argument contains both single and double quotes! *)
val quote : string -> string

(* The corresponding functions from [String] have become obsolescent
   with O'Caml~4.0.3.  Quanrantine them here. *)
val uppercase : string -> string
val lowercase : string -> string

(* Ignore the case in comparisons. *)
val compare_caseless :  string -> string -> int

(* Match the regular expression
   \texttt{\lbrack A-Za-z\rbrack\lbrack A-Za-z0-9\_\rbrack*} *)
val valid_fortran_id : string -> bool

(* Replace any invalid character by ['_'] and prepend ["N_"] iff
   the string doesn't start with a letter. *)
val sanitize_fortran_id : string -> string

module Test : sig val suite : OUnit.test end

