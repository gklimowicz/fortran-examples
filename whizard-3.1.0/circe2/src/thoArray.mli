(* circe2/thoArray.mli --  *)
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

exception Out_of_bounds of int * int 

(* Interpret optional array boundaries.
   Assuming that [Array.length a] $\mapsto n$, we have
   \begin{itemize}
     \item\relax [decode_inf a] $\mapsto 0$
     \item\relax [decode_sup a] $\mapsto n-1$
     \item\relax [decode_inf ~inf:i a] $\mapsto i$ for $0\le i \le n-1$
     \item\relax [decode_sup ~sup:i a] $\mapsto i$ for $0\le i \le n-1$
     \item\relax [decode_inf ~inf:(-i) a] $\mapsto n-i$ for $1\le i \le n$
     \item\relax [decode_sup ~sup:(-i) a] $\mapsto n-i$ for $1\le i \le n$
     \item\relax [decode_inf ~inf:i a] raises [Out_of_bounds]
        for $i\ge n \lor i<-n$
     \item\relax [decode_sup ~sup:i a] raises [Out_of_bounds]
        for $i\ge n \lor i<-n$
   \end{itemize}
   In particular
   \begin{itemize}
     \item\relax [decode_inf ~inf:(-2) a] $\mapsto n-2$,
       i.\,e.~the idex of the next-to-last element.
     \item\relax [decode_sup ~sup:(-1) a] $\mapsto n-1$,
       i.\,e.~the idex of the last element.  
   \end{itemize} *)
val decode_inf : ?inf:int -> 'a array -> int
val decode_sup : ?sup:int -> 'a array -> int

(* Just like the functions from [Array] of the same name,
   but acting only on the subarray specified by the optional
   [~inf] and [~sup], interpreted as above.
   E.\,g.~[copy ~inf:1 ~sup:(-2) a] chops off the first
   and last elements. *)
val map : ?inf:int -> ?sup:int -> ('a -> 'b) -> 'a array -> 'b array
val copy : ?inf:int -> ?sup:int -> 'a array -> 'a array
val iter : ?inf:int -> ?sup:int -> ('a -> unit) -> 'a array -> unit
val fold_left : ?inf:int -> ?sup:int ->
  ('a -> 'b -> 'a) -> 'a -> 'b array -> 'a

(* A convenience function. *)
val sum_float : ?inf:int -> ?sup:int -> float array -> float

val suite : OUnit.test

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
