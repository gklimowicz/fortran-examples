(* partition.mli --

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

(* [pairs n n1 n2] returns all (unordered) pairs of integers with the
    sum~$n$ in the range from~[n1] to~[n2]. *)
val pairs : int -> int -> int -> (int * int) list
val triples : int -> int -> int -> (int * int * int) list

(* [tuples d n n_min n_max] returns
   all~$\lbrack n_1; n_2; \ldots; n_d\rbrack$
   with~$n_{\min}\le n_1\le n_2\le\ldots\le n_d\le n_{\max}$ and
   \begin{equation}
     \sum_{i=1}^d n_i = n
   \end{equation} *)
val tuples : int -> int -> int -> int -> int list list

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)

