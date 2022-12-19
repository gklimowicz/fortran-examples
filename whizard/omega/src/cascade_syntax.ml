(* cascade_syntax.ml --

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

(* Concerning the Gaussian propagators, we admit the following: In
   principle, they would allow for flavor sums like the off-shell 
   lines, but for all practical purposes they are used only for
   determining the significance of a specified intermediate state. 
   So we select them in the same manner as on-shell states. *)

(* [False] is probably redundant.  *)

type ('flavor, 'p, 'constant) t =
  | True
  | False
  | On_shell of 'flavor list * 'p
  | On_shell_not of 'flavor list * 'p
  | Off_shell of 'flavor list * 'p
  | Off_shell_not of 'flavor list * 'p
  | Gauss of 'flavor list * 'p
  | Gauss_not of 'flavor list * 'p
  | Any_flavor of 'p
  | And of ('flavor, 'p, 'constant) t list
  | X_Flavor of 'flavor list
  | X_Vertex of 'constant list * 'flavor list list

let mk_true () = True
let mk_false () = False
let mk_on_shell f p = On_shell (f, p)
let mk_on_shell_not f p = On_shell_not (f, p)
let mk_off_shell f p = Off_shell (f, p)
let mk_off_shell_not f p = Off_shell_not (f, p)
let mk_gauss f p = Gauss (f, p)
let mk_gauss_not f p = Gauss_not (f, p)
let mk_any_flavor p = Any_flavor p

let mk_and c1 c2 =
  match c1, c2 with
  | c, True | True, c -> c
  | c, False | False, c -> False
  | And cs, And cs' -> And (cs @ cs')
  | And cs, c | c, And cs -> And (c::cs)
  | c, c' -> And [c; c']

let mk_x_flavor f = X_Flavor f
let mk_x_vertex c fs = X_Vertex (c, fs)
    
let to_string flavor_to_string momentum_to_string coupling_to_string cascades =
  let flavors_to_string fs =
    String.concat ":" (List.map flavor_to_string fs)
  and couplings_to_string cs =
    String.concat ":" (List.map coupling_to_string cs) in
  let rec to_string' = function
    | True -> "true"
    | False -> "false"
    | On_shell (fs, p) ->
        momentum_to_string p ^ " = " ^ flavors_to_string fs
    | On_shell_not (fs, p) ->
        momentum_to_string p ^ " = !" ^ flavors_to_string fs
    | Off_shell (fs, p) ->
        momentum_to_string p  ^ " ~ " ^ flavors_to_string fs
    | Off_shell_not (fs, p) ->
        momentum_to_string p  ^ " ~ !" ^ flavors_to_string fs
    | Gauss (fs, p) ->
        momentum_to_string p ^ " # " ^ flavors_to_string fs
    | Gauss_not (fs, p) ->
        momentum_to_string p ^ " # !" ^ flavors_to_string fs
    | Any_flavor p ->
        momentum_to_string p ^ " ~ ?"
    | And cs ->
        String.concat " && " (List.map (fun c -> "(" ^ to_string' c ^ ")") cs)
    | X_Flavor fs ->
        "!" ^ String.concat ":" (List.map flavor_to_string fs)
    | X_Vertex (cs, fss) ->
        "^" ^ couplings_to_string cs ^
        "[" ^ (String.concat "," (List.map flavors_to_string fss)) ^ "]"
  in
  to_string' cascades

let int_list_to_string p =
  String.concat "+" (List.map string_of_int (List.sort compare p))

exception Syntax_Error of string * int * int

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)

