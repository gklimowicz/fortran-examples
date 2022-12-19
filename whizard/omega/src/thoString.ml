(* thoString.ml --

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

let strip_prefix p s =
  let lp = String.length p
  and ls = String.length s in
  if lp > ls then
    s
  else
    let rec strip_prefix' i =
      if i >= lp then
	String.sub s i (ls - i)
      else if p.[i] <> s.[i] then
	s
      else
	strip_prefix' (succ i)
    in
    strip_prefix' 0

let strip_prefix_star p s =
  let ls = String.length s in
  if ls < 1 then
    s
  else
    let rec strip_prefix_star' i =
      if i < ls then begin
	if p <> s.[i] then
	  String.sub s i (ls - i)
	else
	  strip_prefix_star' (succ i)
      end else
	""
    in
    strip_prefix_star' 0

let strip_required_prefix p s =
  let lp = String.length p
  and ls = String.length s in
  if lp > ls then
    invalid_arg ("strip_required_prefix: expected `" ^ p ^ "' got `" ^ s ^ "'")
  else
    let rec strip_prefix' i =
      if i >= lp then
	String.sub s i (ls - i)
      else if p.[i] <> s.[i] then
	invalid_arg ("strip_required_prefix: expected `" ^ p ^ "' got `" ^ s ^ "'")
      else
	strip_prefix' (succ i)
    in
    strip_prefix' 0

let strip_from_first c s =
  try
    String.sub s 0 (String.index s c)
  with
  | Not_found -> s

let strip_from_last c s =
  try
    String.sub s 0 (String.rindex s c)
  with
  | Not_found -> s

let index_string pat s =
  let lpat = String.length pat
  and ls = String.length s in
  if lpat = 0 then
    0
  else
    let rec index_string' n =
      let i = String.index_from s n pat.[0] in
      if i + lpat > ls then
        raise Not_found
      else
        if String.compare pat (String.sub s i lpat) = 0 then
          i
        else
          index_string' (succ i)
    in
    index_string' 0

let quote s =
  if String.contains s ' ' || String.contains s '\n' then begin
    if String.contains s '"' then
      "'" ^ s ^ "'"
    else
      "\"" ^ s ^ "\""
  end else
    s

let uppercase = String.uppercase_ascii
let lowercase = String.lowercase_ascii

let compare_caseless s1 s2 =
  String.compare (lowercase s1) (lowercase s2)

let is_alpha c =
  ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z')

let is_numeric c =
  '0' <= c && c <= '9'

let is_alphanum c =
  is_alpha c || is_numeric c || c = '_'

let valid_fortran_id s =
  let rec valid_fortran_id' n =
    if n < 0 then
      false
    else if n = 0 then
      is_alpha s.[0]
    else if is_alphanum s.[n] then
      valid_fortran_id' (pred n)
    else
      false in
  valid_fortran_id' (pred (String.length s))

let sanitize_fortran_id s =
  let sanitize s =
    String.map (fun c -> if is_alphanum c then c else '_') s in
  if String.length s <= 0 then
    invalid_arg "ThoString.sanitize_fortran_id: empty"
  else if is_alpha s.[0] then
    sanitize s
  else
    "N_" ^ sanitize s

module Test =
  struct

    open OUnit

    let fortran_empty =
      "empty" >::
	(fun () -> assert_equal false (valid_fortran_id ""))

    let fortran_digit =
      "0" >::
	(fun () -> assert_equal false (valid_fortran_id "0"))

    let fortran_digit_alpha =
      "0abc" >::
	(fun () -> assert_equal false (valid_fortran_id "0abc"))

    let fortran_underscore =
      "_" >::
	(fun () -> assert_equal false (valid_fortran_id "_"))

    let fortran_underscore_alpha =
      "_ABC" >::
	(fun () -> assert_equal false (valid_fortran_id "_ABC"))

    let fortran_questionmark =
      "A?C" >::
	(fun () -> assert_equal false (valid_fortran_id "A?C"))

    let fortran_valid =
      "A_xyz_0_" >::
	(fun () -> assert_equal true (valid_fortran_id "A_xyz_0_"))

    let sanitize_digit =
      "0" >::
	(fun () -> assert_equal "N_0" (sanitize_fortran_id "0"))

    let sanitize_digit_alpha =
      "0abc" >::
	(fun () -> assert_equal "N_0abc" (sanitize_fortran_id "0abc"))

    let sanitize_underscore =
      "_" >::
	(fun () -> assert_equal "N__" (sanitize_fortran_id "_"))

    let sanitize_underscore_alpha =
      "_ABC" >::
	(fun () -> assert_equal "N__ABC" (sanitize_fortran_id "_ABC"))

    let sanitize_questionmark =
      "A?C" >::
	(fun () -> assert_equal "A_C" (sanitize_fortran_id "A?C"))

    let sanitize_valid =
      "A_xyz_0_" >::
	(fun () -> assert_equal "A_xyz_0_" (sanitize_fortran_id "A_xyz_0_"))

    let suite_fortran =
      "valid_fortran_id" >:::
        [fortran_empty;
         fortran_digit;
         fortran_digit_alpha;
         fortran_underscore;
         fortran_underscore_alpha;
         fortran_questionmark;
         fortran_valid]

    let suite_sanitize =
      "sanitize_fortran_id" >:::
        [sanitize_digit;
         sanitize_digit_alpha;
         sanitize_underscore;
         sanitize_underscore_alpha;
         sanitize_questionmark;
         sanitize_valid]

    let suite =
      "ThoString" >:::
	[suite_fortran;
         suite_sanitize]

  end

