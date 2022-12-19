(* tree2.ml --

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

(* Dependency trees for wavefunctions. *)

type ('n, 'e) t = 
  | Node of ('e * 'n * ('n, 'e) t list) list
  | Leaf of 'n

let leaf node = Leaf node

let sort_children (edge, node, children) =
  (edge, node, List.sort compare children)

let cons fusions = Node (List.sort compare (List.map sort_children fusions))

let is_singleton = function
  | Leaf _ -> true
  | _ -> false

let rec to_string n2s e2s = function
  | Leaf n -> n2s n
  | Node [children] ->
     children_to_string n2s e2s children
  | Node children2 ->
     "{ " ^
       String.concat " | " (List.map (children_to_string n2s e2s) children2) ^
         " }"

and children_to_string n2s e2s (e, n, children) =
  "(" ^ (match e2s e with "" -> "" | s -> s ^ ">") ^ n2s n ^ ":" ^
    (String.concat "," (List.map (to_string n2s e2s) children)) ^ ")"

  
let rec to_channel ch n2s e2s = function
  | Leaf n -> Printf.fprintf ch "%s" (n2s n)
  | Node [] -> Printf.fprintf ch "{ }";
  | Node [children] -> children_to_channel ch n2s e2s children
  | Node (children::children2) ->
     Printf.fprintf ch "{ ";
     children_to_channel ch n2s e2s children;
     List.iter
       (fun children ->
         Printf.fprintf ch " \\\n   | ";
         children_to_channel ch n2s e2s children)
       children2;
     Printf.fprintf ch " }"

and children_to_channel ch n2s e2s (e, n, children) =
  Printf.fprintf ch "(";
  begin match e2s e with
  | "" -> ()
  | s -> Printf.fprintf ch "%s>" s
  end;
  Printf.fprintf ch "%s:" (n2s n);
  begin match children with
  | [] -> ()
  | [child] -> to_channel ch n2s e2s child
  | child::children ->
     to_channel ch n2s e2s child;
     List.iter
       (fun child ->
         Printf.fprintf ch ",";
         to_channel ch n2s e2s child)
       children
  end;
  Printf.fprintf ch ")"

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
