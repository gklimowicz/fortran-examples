(* cascade_syntax.mli --

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

val mk_true : unit -> ('flavor, 'p, 'constant) t
val mk_false : unit -> ('flavor, 'p, 'constant) t
val mk_on_shell : 'flavor list -> 'p -> ('flavor, 'p, 'constant) t
val mk_on_shell_not : 'flavor list -> 'p -> ('flavor, 'p, 'constant) t
val mk_off_shell : 'flavor list -> 'p -> ('flavor, 'p, 'constant) t
val mk_off_shell_not : 'flavor list -> 'p -> ('flavor, 'p, 'constant) t
val mk_gauss : 'flavor list -> 'p -> ('flavor, 'p, 'constant) t
val mk_gauss_not : 'flavor list -> 'p -> ('flavor, 'p, 'constant) t
val mk_any_flavor : 'p -> ('flavor, 'p, 'constant) t
val mk_and : ('flavor, 'p, 'constant) t ->
  ('flavor, 'p, 'constant) t -> ('flavor, 'p, 'constant) t
val mk_x_flavor : 'flavor list -> ('flavor, 'p, 'constant) t
val mk_x_vertex : 'constant list -> 'flavor list list ->
  ('flavor, 'p, 'constant) t

val to_string : ('flavor -> string) -> ('p -> string) ->
  ('constant -> string) -> ('flavor, 'p, 'constant) t -> string

exception Syntax_Error of string * int * int

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)

