(* circe2/float.ml --  *)
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

open Printf

module type T =
  sig
    type t
    (* Difference between~$1.0$ and the minimum float greater than~$1.0$ *)
    val epsilon : t
    val to_string : t -> string
    val input_binary_float : in_channel -> float
    val input_binary_floats : in_channel -> float array -> unit
  end

module Double =
  struct

    type t = float

  (* Difference between~$1.0$ and the minimum float greater than~$1.0$
     \begin{dubious}
       This is the hard coded value for double precision on
       Linux/Intel.  We should determine this \emph{machine dependent}
       value during configuration.
     \end{dubious} *)
    let epsilon = 2.2204460492503131e-16
    let little_endian = true

    let to_string x =
      let s = Bytes.of_string (sprintf "%.17E" x) in
      for i = 0 to Bytes.length s - 1 do
        let c = Bytes.get s i in
        if c = 'e' || c = 'E' then
          Bytes.set s i 'D'
      done;
      Bytes.to_string s

    (* Identity floatingpoint numbers that are indistinguishable from
       integers for more concise printing. *)

    type int_or_float =
      | Int of int
      | Float of float

    let float_min_int = float min_int
    let float_max_int = float max_int

    let soft_truncate x =
      let eps = 2.0 *. abs_float x *. epsilon in
      if x >= 0.0 then begin
        if x > float_max_int then
          Float x
        else if x -. floor x <= eps then
          Int (int_of_float x)
        else if ceil x -. x <= eps then
          Int (int_of_float x + 1)
        else
          Float x
      end else begin
        if x < float_min_int then
          Float x
        else if ceil x -. x <= eps then
          Int (int_of_float x)
        else if x -. floor x <= eps then
          Int (int_of_float x - 1)
        else
          Float x
      end

    let to_short_string x =
      match soft_truncate x with
      | Int i -> string_of_int i ^ "D0"
      | Float x -> to_string x

    (* Suggested by Xavier Leroy: *)

    let output_float_big_endian oc f =
      let n = ref (Int64.bits_of_float f) in
      for i = 0 to 7 do
        output_byte oc (Int64.to_int (Int64.shift_right_logical !n 56));
        n := Int64.shift_left !n 8
      done

    let output_float_little_endian oc f =
      let n = ref (Int64.bits_of_float f) in
      for i = 0 to 7 do
        output_byte oc (Int64.to_int !n);
        n := Int64.shift_right_logical !n 8
      done
        
    let input_float_big_endian ic =
      let n = ref Int64.zero in
      for i = 0 to 7 do
        let b = input_byte ic in
        n := Int64.logor (Int64.shift_left !n 8) (Int64.of_int b)
      done;
      Int64.float_of_bits !n

    let input_float_little_endian ic =
      let n = ref Int64.zero in
      for i = 0 to 7 do
        let b = input_byte ic in
        n := Int64.logor !n (Int64.shift_left (Int64.of_int b) (i*8))
      done;
      Int64.float_of_bits !n

    let input_binary_float = input_float_little_endian

    let input_binary_floats ic array =
      for i = 0 to Array.length array - 1 do
        array.(i) <- input_binary_float ic
      done

  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
