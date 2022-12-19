(* partition.ml --

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

(* All unordered pairs of integers with the same sum~$n$ in a given
   range~$\{n_1,\ldots,n_2\}$:
   \begin{equation}
     \text{\ocwlowerid{pairs}}: (n, n_1, n_2) \to
        \bigl\{ (i,j) \,\vert\, i+j=n
                \land n_1\le i \le j \le n_2 \bigr\}
   \end{equation} *)

let rec pairs' acc n1 n2 =
  if n1 > n2 then
    List.rev acc
  else
    pairs' ((n1, n2) :: acc) (succ n1) (pred n2)

let pairs sum min_n1 max_n2 =
  let n1 = max min_n1 (sum - max_n2) in
  let n2 = sum - n1 in
  if n2 <= max_n2 then
    pairs' [] n1 n2
  else
    []

let rec tuples d sum n_min n_max =
  if d <= 0 then
    invalid_arg "tuples"
  else if d > 1 then
    tuples' d sum n_min n_max n_min
  else if sum >= n_min && sum <= n_max then
    [[sum]]
  else
    []

and tuples' d sum n_min n_max n =
  if n > n_max then
    []
  else
    List.fold_right (fun l ll -> (n :: l) :: ll)
      (tuples (pred d) (sum - n) (max n_min n) n_max)
      (tuples' d sum n_min n_max (succ n))

(* \begin{dubious}
     When I find a little spare time, I can provide a dedicated implementation,
     but we \emph{know} that [Impossible] is \emph{never} raised and the present
     approach is just as good (except for a possible tiny inefficiency).
   \end{dubious} *)
exception Impossible of string
let impossible name = raise (Impossible name)

let triples sum n_min n_max =
  List.map (function [n1; n2; n3] -> (n1, n2, n3) | _ -> impossible "triples")
    (tuples 3 sum n_min n_max)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)



