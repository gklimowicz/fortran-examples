(* young.mli --

   Copyright (C) 2022- by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>

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

(* Caveat: the following are not optimized for large Young diagrams and
   tableaux.  They are straightforward implementations of the
   definitions, since we are unlikely to meet large diagrams.

   To make matters worse, native integer arithmetic will overflow
   already for diagrams with more than 20 cells.
   Since the [Num] library has been removed from the O'Caml
   distribution with version 4.06, we can not use it as
   a shortcut.  Requiring Whizard/O'Mega users to install
   [Num] or its successor [Zarith] is probably not worth
   the effort. *)

(* \ytableausetup{centertableaux,smalltableaux} *)

(* \thocwmodulesection{Young Diagrams} *)
   
(* Young diagrams can be represented by a non-increasing list
   of positive integers, corresponding to the number of boxes
   in each row:
   \begin{equation}
     \ydiagram{5,4,4,2} \Longleftrightarrow \lbrack 5;4;4;2 \rbrack
   \end{equation} *)
type diagram = int list

(* Check that the diagram is valid, i.\,e.~the number of boxes
   is non-increasing from top to bottom. *)
val valid_diagram : diagram -> bool

(* Count the number of cells. *)
val num_cells_diagram : diagram -> int

(* Conjugate a diagram:
   \begin{equation}
     \ydiagram{5,4,4,2} \mapsto \ydiagram{4,4,3,3,1}
   \end{equation} *)
val conjugate_diagram : diagram -> diagram

(* The product of all the ``hook lengths'' in the diagram, e.\,g.
   \begin{equation}
     \ydiagram{5,4,4,2}
     \mapsto \ytableaushort{87541,6532,5421,21}
     \mapsto 8 \cdot 7 \cdot 6 \cdot 5^3 \cdot 4^2 \cdot 3 \cdot 2^3
     = 16128000
   \end{equation}
   where the intermediate step is only for illustration and does not
   represent a Young tableau! *)
val hook_lengths_product : diagram -> int

(* Number of standard tableaux corresponding to the diagram.
   Also, the dimension of the representation of~$S_n$ described
   by this diagram
   \begin{equation}
     d = \frac{n!}{\prod_{i=1}^n h_i}
   \end{equation}
   with~$n$ the number of cells and~$h_i$ the hook length of
   the $i$th cell. *)
val num_standard_tableaux : diagram -> int

(* Normalization of the projector on the representation of $\mathrm{GL(N)}$
   described by the diagram
   \begin{equation}
     \alpha = \frac{\prod_{R} |R|!\prod_{C} |C|!}{\prod_{i=1}^n h_i}
   \end{equation}
   with~$|R|$ and~$|C|$ the lengths of the row~$R$ and column~$C$,
   respectively.  Returned as a pair of numerator and denominator,
   because it is not guaranteed to be integer. *)
val normalization : diagram -> int * int

(* \thocwmodulesection{Young Tableaux} *)
(* There is an obvious representation as a list of lists:
   \begin{equation}
     \ytableaushort{023,14}
     \Longleftrightarrow
     \lbrack \lbrack 0; 2; 3 \rbrack;
             \lbrack 1; 4 \rbrack \rbrack
   \end{equation} *)
type 'a tableau = 'a list list

(* Ignoring the contents of the cells of a Young tableau produces
   a unique corresponding Young diagram.
   \begin{equation}
     \ytableaushort{023,14}
     \mapsto \ydiagram{3,2}
   \end{equation} *)
val diagram_of_tableau : 'a tableau -> diagram

(* The number of columns must be non-increasing.  Obviously,
   [valid_tableau] is the composition of [diagram_of_tableau]
   and [valid_diagram].*)
val valid_tableau : 'a tableau -> bool

(* A tableau is called \textit{semistandard}, iff the entries
   don't increase along rows and strictly increase along columns.
   Therefore, the conjugate of a semistandard tableau is \emph{not}
   necessarily semistandard. *)
val semistandard_tableau : 'a tableau -> bool

(* A tableau is called \textit{standard}, iff it is semistandard
   and the entries are an uninterrupted sequence of natural numbers.
   If the optional [offset] is specified, it must match the smallest
   of these numbers.  Some authors expect [offset=1], but we want
   to be able to start from 0 as well.
   The conjugate of a standard tableau is again a standard tableau. *)
val standard_tableau : ?offset:int -> int tableau -> bool

(* The contents of the cells and their number. *)
val cells_tableau : 'a tableau -> 'a list
val num_cells_tableau : 'a tableau -> int

(* Conjugate a Young tableau
   \begin{equation}
     \ytableaushort{023,14}
     \mapsto \ytableaushort{01,24,3}
   \end{equation} *)
val conjugate_tableau : 'a tableau -> 'a tableau

(* \thocwmodulesection{Unit Tests} *)
module type Test =
  sig
    val suite : OUnit.test
    val suite_long : OUnit.test
  end

module Test : Test

