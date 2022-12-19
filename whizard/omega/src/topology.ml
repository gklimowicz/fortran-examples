(* topology.ml --

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

module type T =
  sig
    type partition
    val partitions : int -> partition list
    type 'a children
    val keystones : 'a list -> ('a list * 'a list children list) list
    val max_subtree : int -> int
    val inspect_partition : partition -> int list
  end

(* \thocwmodulesection{Factorizing Diagrams for $\phi^3$} *)

module Binary =
  struct
    type partition = int * int * int
    let inspect_partition (n1, n2, n3) = [n1; n2; n3]

(* One way~\cite{ALPHA:1997} to lift the degeneracy is to select the
   vertex that is closest to the center
   (see table~\ref{tab:partition}):
   \begin{equation}
   \label{eq:partition}
     \text{\ocwlowerid{partitions}}: n \to
        \bigl\{ (n_1,n_2,n_3) \,\vert\, n_1 + n_2 + n_3 = n
                \land n_1 \le n_2 \le n_3 \le \lfloor n/2 \rfloor \bigr\}
   \end{equation}
   Other, less symmetric, approaches are possible.  The simplest
   of these is: choose the vertex adjacent to a fixed
   external line~\cite{HELAC:2000}.  They will be made available
   for comparison in the future.
   \begin{table}
     \begin{center}
       \begin{tabular}{ r | l }
         [n]& [partitions n] \\\hline
          4 & (1,1,2) \\
          5 & (1,2,2) \\
          6 & (1,2,3), (2,2,2) \\
          7 & (1,3,3), (2,2,3) \\
          8 & (1,3,4), (2,2,4), (2,3,3) \\
          9 & (1,4,4), (2,3,4), (3,3,3) \\
         10 & (1,4,5), (2,3,5), (2,4,4), (3,3,4) \\
         11 & (1,5,5), (2,4,5), (3,3,5), (3,4,4) \\
         12 & (1,5,6), (2,4,6), (2,5,5), (3,3,6), (3,4,5), (4,4,4) \\
         13 & (1,6,6), (2,5,6), (3,4,6), (3,5,5), (4,4,5) \\
         14 & (1,6,7), (2,5,7), (2,6,6), (3,4,7), (3,5,6), (4,4,6), (4,5,5) \\
         15 & (1,7,7), (2,6,7), (3,5,7), (3,6,6), (4,4,7), (4,5,6), (5,5,5) \\
         16 & (1,7,8), (2,6,8), (2,7,7), (3,5,8), (3,6,7), (4,4,8), (4,5,7), (4,6,6), (5,5,6) 
       \end{tabular}
     \end{center}
     \caption{\label{tab:partition} [partitions n] for moderate values
       of [n].}
   \end{table} *)

(* An obvious consequence of~$n_1 + n_2 + n_3 = n$
   and~$n_1 \le n_2 \le n_3$ is $n_1\le\lfloor n/3 \rfloor$: *)
    let rec partitions' n n1 =
      if n1 > n / 3 then
        []
      else
        List.map (fun (n2, n3) -> (n1, n2, n3))
          (Partition.pairs (n - n1) n1 (n / 2)) @ partitions' n (succ n1)

    let partitions n = partitions' n 1

(* \begin{figure}
     \begin{center}
        \hfil\\
        \begin{fmfgraph*}(25,20)
          \fmfstraight
          \fmfbottomn{b}{2}
          \fmftopn{t}{1}
          \fmf{plain}{t1,v}
          \fmf{plain}{b1,v}
          \fmf{plain}{b2,v}
          \fmfv{d.sh=circle,d.f=empty,d.si=18pt,l=$n$,l.d=0}{b1}  
          \fmfv{d.sh=circle,d.f=empty,d.si=18pt,l=$n$,l.d=0}{b2}  
          \fmfv{d.sh=circle,d.f=empty,d.si=18pt,l=$n$,l.d=0}{t1}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{v}  
        \end{fmfgraph*}
        \qquad\qquad\qquad\qquad
        \begin{fmfgraph*}(25,20)
          \fmfstraight
          \fmfbottomn{b}{3}
          \fmftopn{t}{1}
          \fmf{plain}{b1,t1}
          \fmf{plain}{b2,t1}
          \fmf{plain}{b3,t1}
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n$,l.d=0}{b1}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n$,l.d=0}{b2}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n$,l.d=0}{b3}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{t1}  
        \end{fmfgraph*}
     \end{center} 
     \caption{\label{fig:nnn} Topologies with a blatant three-fold
       permutation symmetry, if the number of external lines is a
       multiple of three}
   \end{figure}
   \begin{figure}
     \begin{center}
        \begin{fmfgraph*}(15,20)
          \fmfstraight
          \fmfbottomn{b}{2}
          \fmftopn{t}{1}
          \fmf{plain}{b1,v}
          \fmf{plain}{b2,v}
          \fmf{plain,tension=2}{t1,v}
          \fmfv{d.sh=circle,d.f=empty,d.si=18pt,l=$n$,l.d=0}{t1}  
          \fmfv{d.sh=circle,d.f=empty,d.si=18pt,l=$n'$,l.d=0}{b1}  
          \fmfv{d.sh=circle,d.f=empty,d.si=18pt,l=$n'$,l.d=0}{b2}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{v}
        \end{fmfgraph*}
        \qquad\qquad\qquad\qquad
        \begin{fmfgraph*}(25,20)
          \fmfstraight
          \fmfbottomn{b}{3}
          \fmftopn{t}{1}
          \fmf{plain}{b1,t1}
          \fmf{plain}{b2,t1}
          \fmf{plain}{b3,t1}
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n$,l.d=0}{b1}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=30pt,l=$n'$,l.d=0}{b2}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=30pt,l=$n'$,l.d=0}{b3}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{t1}
          \fmfshift{(0,.2h)}{b1}
        \end{fmfgraph*}
        \qquad\qquad
        \begin{fmfgraph*}(25,20)
          \fmfstraight
          \fmfbottomn{b}{3}
          \fmftopn{t}{1}
          \fmf{plain}{b1,t1}
          \fmf{plain}{b2,t1}
          \fmf{plain}{b3,t1}
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n'$,l.d=0}{b1}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n'$,l.d=0}{b2}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=30pt,l=$n$,l.d=0}{b3}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{t1}
          \fmfshift{(0,.2h)}{b1,b2}
        \end{fmfgraph*}
     \end{center} 
     \caption{\label{fig:n1n2n2} Topologies with a blatant two-fold symmetry.}
   \end{figure}
   \begin{figure}
     \begin{center}
        \hfil\\
        \begin{fmfgraph*}(25,20)
          \fmfstraight
          \fmfbottomn{b}{3}
          \fmftopn{t}{1}
          \fmf{plain}{b1,t1}
          \fmf{plain}{b2,t1}
          \fmf{plain}{b3,t1}
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n_1$,l.d=0}{b1}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=30pt,l=$n_2$,l.d=0}{b2}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=35pt,l=$n_3$,l.d=0}{b3}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{t1}  
          \fmfshift{(0,.30h)}{b1}
          \fmfshift{(0,.15h)}{b2}
        \end{fmfgraph*}
        \qquad\qquad
        \begin{fmfgraph*}(25,20)
          \fmfstraight
          \fmfbottomn{b}{3}
          \fmftopn{t}{1}
          \fmf{plain}{b1,t1}
          \fmf{plain}{b2,t1}
          \fmf{plain}{b3,t1}
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n$,l.d=0}{b1}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=25pt,l=$n$,l.d=0}{b2}  
          \fmfv{d.sh=triangle,d.f=empty,d.si=35pt,l=$2n$,l.d=0}{b3}  
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{t1}  
          \fmfshift{(0,.20h)}{b1}
          \fmfshift{(0,.20h)}{b2}
        \end{fmfgraph*}
     \end{center} 
     \caption{\label{fig:n1n2n3} If~$n_3=n_1+n_2$, the apparently
       asymmetric topologies on the left hand side have a non obvious
       two-fold symmetry, that exchanges the two halves.  Therefore,
       the topologies on the right hand side have a four fold symmetry.}
   \end{figure} *)

    type 'a children = 'a Tuple.Binary.t

(* There remains one peculiar case, when the number of external lines is
   even and~$n_3=n_1+n_2$ (cf.~figure~\ref{fig:n1n2n3}).
   Unfortunately, this reflection symmetry is not respected by the equivalence
   classes. E.\,g.
   \begin{equation}
     \{1\}\{2,3\}\{4,5,6\}\mapsto\bigl\{
       \{4\}\{5,6\}\{1,2,3\}; \{5\}\{4,6\}\{1,2,3\}; \{6\}\{4,5\}\{1,2,3\} \bigr\}
   \end{equation}
   However, these reflections will always exchange the two halves
   and a representative can be chosen by requiring that one fixed
   momentum remains in one half.  We choose to filter out the half
   of the partitions where the element~[p] appears in the second
   half, i.\,e.~the list of length~[n3].

   Finally, a closed expression for the number of Feynman diagrams
   in the equivalence class $(n_1,n_2,n_3)$ is
   \begin{equation}
     N(n_1,n_2,n_3) =
       \frac{(n_1+n_2+n_3)!}{S(n_1,n_2,n_3)}
       \prod_{i=1}^{3} \frac{(2n_i-3)!!}{n_i!}
   \end{equation}
   where the symmetry factor from the above arguments is
   \begin{equation}
   \label{eq:S(1,2,3)}
     S(n_1,n_2,n_3) =
       \begin{cases}
          3!      & \text{for $n_1 = n_2 = n_3$} \\
          2\cdot2 & \text{for $n_3 = 2n_1 = 2n_2$} \\
          2       & \text{for $n_1 = n_2 \lor n_2 = n_3$} \\
          2       & \text{for $n_1 + n_2 = n_3$} 
       \end{cases}
   \end{equation}
   Indeed, the sum of all Feynman diagrams
   \begin{equation}
   \label{eq:keystone-check}
     \sum_{\substack{n_1 + n_2 + n_3 = n\\
                     1 \le n_1 \le n_2 \le n_3 \le \lfloor n/2 \rfloor}}
        N(n_1,n_2,n_3) = (2n-5)!!
   \end{equation}
   can be checked numerically for large values of $n=n_1+n_2+n_3$,
   verifying the symmetry factor (see table~\ref{tab:keystone-check}).
   \begin{dubious}
     P.\,M.~claims to have seen similar formulae in the context of
     Young tableaux.  That's a good occasion to read the new edition
     of Howard's book \ldots
   \end{dubious}
   \begin{table}
     \begin{center}
       \begin{tabular}{ r | r | l }
         $n$ & $(2n-5)!!$ & $\sum N(n_1,n_2,n_3)$ \\\hline
          4  &         3 & $3\cdot(1,1,2)$ \\
          5  &        15 & $15\cdot(1,2,2)$ \\
          6  &       105 & $90\cdot(1,2,3) + 15\cdot(2,2,2)$ \\
          7  &       945 & $630\cdot(1,3,3) + 315\cdot(2,2,3)$ \\
          8  &     10395 & $6300\cdot(1,3,4) + 1575\cdot(2,2,4) + 2520\cdot(2,3,3)$ \\
          9  &    135135 & $70875\cdot(1,4,4) + 56700\cdot(2,3,4) + 7560\cdot(3,3,3)$ \\
         10  &   2027025 & $992250\cdot(1,4,5) + 396900\cdot(2,3,5)$ \\
             &           & \quad$\mbox{}+ 354375\cdot(2,4,4) + 283500\cdot(3,3,4)$ \\
         11  &  34459425 & $15280650\cdot(1,5,5) + 10914750\cdot(2,4,5)$ \\
             &           & \quad$\mbox{}+ 4365900\cdot(3,3,5) + 3898125\cdot(3,4,4)$ \\
         12  & 654729075 & $275051700\cdot(1,5,6) + 98232750\cdot(2,4,6)$ \\
             &           & \quad$\mbox{}+ 91683900\cdot(2,5,5)+ 39293100\cdot(3,3,6)$ \\
             &           & \quad$\mbox{}+ 130977000\cdot(3,4,5) + 19490625\cdot(4,4,4)$
       \end{tabular}
     \end{center}
     \caption{\label{tab:keystone-check} Equation~(\ref{eq:keystone-check}) for
       small values of $n$.}
   \end{table} *)

(* Return a list of all inequivalent partitions of the list~[l] in three
   lists of length [n1], [n2] and [n3], respectively. Common first lists
   are factored. This is nothing more than a typedafe wrapper around
   [Combinatorics.factorized_keystones].  *)

    exception Impossible of string
    let tuple_of_list2 = function
      | [x1; x2] -> Tuple.Binary.of2 x1 x2
      | _ -> raise (Impossible "Topology.tuple_of_list")

    let keystone (n1, n2, n3) l =
      List.map (fun (p1, p23) -> (p1, List.rev_map tuple_of_list2 p23))
        (Combinatorics.factorized_keystones [n1; n2; n3] l)

    let keystones l =
      ThoList.flatmap (fun n123 -> keystone n123 l) (partitions (List.length l))

    let max_subtree n = n / 2

  end
    
(* \thocwmodulesection{Factorizing Diagrams for $\sum_n\lambda_n\phi^n$} *)

(* \begin{figure}
     \begin{center}
        \begin{fmfgraph}(25,20)
          \fmfleftn{l}{3}
          \fmfrightn{r}{3}
          \fmf{plain}{l1,v4}
          \fmf{plain}{l2,v4}
          \fmf{plain}{l3,v4}
          \fmf{plain}{r1,v1}
          \fmf{plain}{r2,v1}
          \fmf{plain}{v1,v2}
          \fmf{plain}{r3,v2}
          \fmf{plain}{v2,v4}
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{v4}  
          \fmfdot{v1,v2}
        \end{fmfgraph}
        \qquad\qquad
        \begin{fmfgraph}(25,20)
          \fmfleftn{l}{3}
          \fmfrightn{r}{3}
          \fmf{plain}{l1,v4}
          \fmf{plain}{l2,v4}
          \fmf{plain}{l3,v4}
          \fmf{plain}{r1,v1}
          \fmf{plain}{r2,v1}
          \fmf{plain}{v1,v2}
          \fmf{plain}{r3,v2}
          \fmf{plain}{v2,v4}
          \fmfv{d.sh=circle,d.f=empty,d.si=5thin}{v2}  
          \fmfdot{v1,v4}
        \end{fmfgraph}
     \end{center} 
     \caption{\label{fig:n1n2n3n4} Degenerate $(1,1,1,3)$ and $(1,2,3)$.}
   \end{figure} *)

(* Mixed $\phi^n$ adds new degeneracies, as in figure~\ref{fig:n1n2n3n4}.
   They appear if and only if one part takes exactly half of the external
   lines and can relate central vertices of different arity. *)

module Nary (B : Tuple.Bound) =
  struct
    type partition = int list
    let inspect_partition p = p

    let partition d sum =
      Partition.tuples d sum 1 (sum / 2)

    let rec partitions' d sum =
      if d < 3 then
        []
      else
        partition d sum @ partitions' (pred d) sum

    let partitions sum = partitions' (succ (B.max_arity ())) sum

(* \begin{table}
     \begin{center}
       \begin{tabular}{ r | r | l }
         $n$ & $\sum$    & $\sum$ \\\hline
          4  &         4 & $1\cdot(1,1,1,1) + 3\cdot(1,1,2)$ \\
          5  &        25 & $10\cdot(1,1,1,2) + 15\cdot(1,2,2)$ \\
          6  &       220 & $40\cdot(1,1,1,3) + 45\cdot(1,1,2,2)
                            + 120\cdot(1,2,3) + 15\cdot(2,2,2)$ \\
          7  &      2485 & $840\cdot(1,1,2,3) + 105\cdot(1,2,2,2)
                            + 1120\cdot(1,3,3) + 420\cdot(2,2,3)$ \\
          8  &     34300 & $5250\cdot(1,1,2,4) + 4480\cdot(1,1,3,3) + 3360\cdot(1,2,2,3)$\\
             &           & \quad$\mbox{}+ 105\cdot(2,2,2,2) + 14000\cdot(1,3,4)$\\
             &           & \quad$\mbox{}+ 2625\cdot(2,2,4) + 4480\cdot(2,3,3)$ \\
          9  &    559405 & $126000\cdot(1,1,3,4) + 47250\cdot(1,2,2,4) + 40320\cdot(1,2,3,3)$\\
             &           & \quad$\mbox{}+ 5040\cdot(2,2,2,3) + 196875\cdot(1,4,4)$\\
             &           & \quad$\mbox{}+ 126000\cdot(2,3,4) + 17920\cdot(3,3,3)$ \\
         10  &  10525900 & $1108800\cdot(1,1,3,5) + 984375\cdot(1,1,4,4) + 415800\cdot(1,2,2,5)$\\
             &           & \quad$\mbox{}+ 1260000\cdot(1,2,3,4) + 179200\cdot(1,3,3,3)
                                        + 78750\cdot(2,2,2,4)$\\
             &           & \quad$\mbox{}+ 100800\cdot(2,2,3,3) + 3465000\cdot(1,4,5)
                                        + 1108800\cdot(2,3,5)$\\
             &           & \quad$\mbox{}+ 984375\cdot(2,4,4) + 840000\cdot(3,3,4)$
       \end{tabular}
     \end{center}
     \caption{\label{tab:keystone-check4}%
       $\mathcal{L}=\lambda_3\phi^3+\lambda_4\phi^4$}
   \end{table}
   \begin{table}
     \begin{center}
       \begin{tabular}{ r | r | l }
         $n$ & $\sum$    & $\sum$ \\\hline
          4  &         4 & $1\cdot(1,1,1,1) + 3\cdot(1,1,2)$ \\
          5  &        26 & $1\cdot(1,1,1,1,1) + 10\cdot(1,1,1,2) + 15\cdot(1,2,2)$ \\
          6  &       236 & $1\cdot(1,1,1,1,1,1) + 15\cdot(1,1,1,1,2) + 40\cdot(1,1,1,3)$\\
             &           & \quad$\mbox{}+ 45\cdot(1,1,2,2) + 120\cdot(1,2,3) + 15\cdot(2,2,2)$ \\
          7  &      2751 & $21\cdot(1,1,1,1,1,2) + 140\cdot(1,1,1,1,3) + 105\cdot(1,1,1,2,2)$\\
             &           & \quad$\mbox{}+ 840\cdot(1,1,2,3) + 105\cdot(1,2,2,2)
                                        + 1120\cdot(1,3,3) + 420\cdot(2,2,3)$ \\
          8  &     39179 & $224\cdot(1,1,1,1,1,3) + 210\cdot(1,1,1,1,2,2) + 910\cdot(1,1,1,1,4)$\\
             &           & \quad$\mbox{}+ 2240\cdot(1,1,1,2,3) + 420\cdot(1,1,2,2,2)
                                        + 5460\cdot(1,1,2,4)$\\
             &           & \quad$\mbox{}+ 4480\cdot(1,1,3,3) + 3360\cdot(1,2,2,3)
                                        + 105\cdot(2,2,2,2)$\\
             &           & \quad$\mbox{}+ 14560\cdot(1,3,4) + 2730\cdot(2,2,4) + 4480\cdot(2,3,3)$
       \end{tabular}
     \end{center}
     \caption{\label{tab:keystone-check6}%
       $\mathcal{L}=\lambda_3\phi^3+\lambda_4\phi^4+\lambda_5\phi^5+\lambda_6\phi^6$}
   \end{table} *)

    module Tuple = Tuple.Nary(B)
    type 'a children = 'a Tuple.t

    let keystones' l =
      let n = List.length l in
      ThoList.flatmap (fun p -> Combinatorics.factorized_keystones p l)
        (partitions n)
     
    let keystones l =
      List.map (fun (bra, kets) -> (bra, List.map Tuple.of_list kets))
        (keystones' l)

    let max_subtree n = n / 2

  end
    
module Nary4 = Nary (struct let max_arity () = 3 end)

(* \thocwmodulesection{Factorizing Diagrams for $\phi^4$} *)

module Ternary =
  struct
    type partition = int * int * int * int
    let inspect_partition (n1, n2, n3, n4) = [n1; n2; n3; n4]
    type 'a children = 'a Tuple.Ternary.t
    let collect4 acc = function
      | [x; y; z; u] -> (x, y, z, u) :: acc
      | _ -> acc
    let partitions n =
      List.fold_left collect4 [] (Nary4.partitions n)
    let collect3 acc = function
      | [x; y; z] -> Tuple.Ternary.of3 x y z :: acc
      | _ -> acc
    let keystones l =
      List.map (fun (bra, kets) -> (bra, List.fold_left collect3 [] kets))
        (Nary4.keystones' l)
    let max_subtree = Nary4.max_subtree
  end
    
(* \thocwmodulesection{Factorizing Diagrams for $\phi^3+\phi^4$} *)

module Mixed23 =
  struct
    type partition =
      | P3 of int * int * int
      | P4 of int * int * int * int
    let inspect_partition = function
      | P3 (n1, n2, n3) -> [n1; n2; n3]
      | P4 (n1, n2, n3, n4) -> [n1; n2; n3; n4]
    type 'a children = 'a Tuple.Mixed23.t
    let collect34 acc = function
      | [x; y; z] -> P3 (x, y, z) :: acc
      | [x; y; z; u] -> P4 (x, y, z, u) :: acc
      | _ -> acc
    let partitions n =
      List.fold_left collect34 [] (Nary4.partitions n)
    let collect23 acc = function
      | [x; y] -> Tuple.Mixed23.of2 x y :: acc
      | [x; y; z] -> Tuple.Mixed23.of3 x y z :: acc
      | _ -> acc
    let keystones l =
      List.map (fun (bra, kets) -> (bra, List.fold_left collect23 [] kets))
        (Nary4.keystones' l)
    let max_subtree = Nary4.max_subtree
  end
    
(* \thocwmodulesection{%
     Diagnostics: Counting Diagrams and Factorizations for $\sum_n\lambda_n\phi^n$} *)

module type Integer =
  sig
    type t
    val zero : t
    val one : t
    val ( + ) : t -> t -> t
    val ( - ) : t -> t -> t
    val ( * ) : t -> t -> t
    val ( / ) : t -> t -> t
    val pred : t -> t
    val succ : t -> t
    val ( = ) : t -> t -> bool
    val ( <> ) : t -> t -> bool
    val ( < ) : t -> t -> bool
    val ( <= ) : t -> t -> bool
    val ( > ) : t -> t -> bool
    val ( >= ) : t -> t -> bool
    val of_int : int -> t
    val to_int : t -> int
    val to_string : t -> string
    val compare : t -> t -> int
    val factorial : t -> t
  end

(* O'Caml's native integers suffice for all applications, but in
   appendix~\ref{sec:count}, we want to use big integers for numeric
   checks in high orders: *)

module Int : Integer =
  struct
    type t = int
    let zero = 0
    let one = 1
    let ( + ) = ( + )
    let ( - ) = ( - )
    let ( * ) = ( * )
    let ( / ) = ( / )
    let pred = pred
    let succ = succ
    let ( = ) = ( = )
    let ( <> ) = ( <> )
    let ( < ) = ( < )
    let ( <= ) = ( <= )
    let ( > ) = ( > )
    let ( >= ) = ( >= )
    let of_int n = n
    let to_int n = n
    let to_string = string_of_int
    let compare = compare
    let factorial = Combinatorics.factorial
  end

module type Count =
  sig
    type integer
    val diagrams : ?f:(integer -> bool) -> integer -> integer -> integer
    val diagrams_via_keystones : integer -> integer -> integer
    val keystones : integer list -> integer
    val diagrams_per_keystone : integer -> integer list -> integer
  end

module Count (I : Integer) =
  struct
    let description = ["(still inoperational) phi^n topology"]

    type integer = I.t
    open I
    let two = of_int 2
    let three = of_int 3

(* If [I.t] is an abstract datatype, the polymorphic [Pervasives.min]
   can fail.  Provide our own version using the specific comparison
   ``[(<=)]''. *)

    let min x y =
      if x <= y then
        x
      else
        y

(* \thocwmodulesubsection{Counting Diagrams for $\sum_n\lambda_n\phi^n$} *)

(* Classes of diagrams are defined by the number of vertices and their
   degrees.  We could use fixed size arrays, but we will use a map
   instead.  For efficiency, we also maintain the number of external
   lines and the total number of propagators. *)

    module IMap = Map.Make (struct type t = integer let compare = compare end)

    type diagram_class = { ext : integer; prop : integer; v : integer IMap.t }

(*i
    let to_string cl =
      IMap.fold
        (fun d n s ->
          s ^ Printf.sprintf ", #%s=%s" (to_string d) (to_string n)) cl.v
        (Printf.sprintf "#ext=%s, #prop=%s"
           (to_string cl.ext) (to_string cl.prop))
i*)

(* The numbers of external lines, propagators and vertices are determined
   by the degrees and multiplicities of vertices:
   \begin{subequations}
   \begin{align}
     E(\{n_3,n_4,\ldots\}) &= 2 + \sum_{d=3}^{\infty} (d-2)n_d \\
     P(\{n_3,n_4,\ldots\}) &= \sum_{d=3}^{\infty} n_d  - 1
                            = V(\{n_3,n_4,\ldots\}) - 1 \\
     V(\{n_3,n_4,\ldots\}) &= \sum_{d=3}^{\infty} n_d
   \end{align}
   \end{subequations} *)

    let num_ext v =
      List.fold_left (fun sum (d, n)  -> sum + (d - two) * n) two v

    let num_prop v =
      List.fold_left (fun sum (_, n)  -> sum + n) (zero - one) v

(* The sum of all vertex degrees must be equal to the number of propagator end
   points.  This can be verified easily:
   \begin{equation}
     2 P(\{n_3,n_4,\ldots\}) + E(\{n_3,n_4,\ldots\}) = \sum_{d=3}^{\infty} dn_d
   \end{equation} *)

    let add_degree map (d, n) =
      if d < three then
	invalid_arg "add_degree: d < 3"
      else if n < zero then
	invalid_arg "add_degree: n <= 0"
      else if n = zero then
	map
      else
	IMap.add d n map

    let create_class v =
      { ext = num_ext v;
        prop = num_prop v;
        v = List.fold_left add_degree IMap.empty v }

    let multiplicity cl d =
      if d >= three then
        try
          IMap.find d cl.v
        with
        | Not_found -> zero
      else
        invalid_arg "multiplicity: d < 3"

(* Remove one vertex of degree [d], maintaining the invariants.  Raises
   [Zero] if all vertices of degree [d] are exhausted.  *)

    exception Zero

    let remove cl d =
      let n = pred (multiplicity cl d) in
      if n < zero then
        raise Zero
      else
        { ext = cl.ext - (d - two);
          prop = pred cl.prop;
          v = if n = zero then
            IMap.remove d cl.v
          else
            IMap.add d n cl.v }

(* Add one vertex of degree [d], maintaining the invariants.  *)

    let add cl d =
      { ext = cl.ext + (d - two);
        prop = succ cl.prop;
        v = IMap.add d (succ (multiplicity cl d)) cl.v }

(* Count the number of diagrams. Any diagram can be obtained recursively either
   from a diagram with one ternary vertex less by insertion if a ternary vertex
   in an internal or external propagator or from a diagram with a higher order
   vertex that has its degree reduced by one:
   \begin{multline}
     D(\{n_3,n_4,\ldots\}) = \\
      \left(P(\{n_3-1,n_4,\ldots\})+E(\{n_3-1,n_4,\ldots\})\right)
      D(\{n_3-1,n_4,\ldots\}) \\
      {} + \sum_{d=4}^{\infty} (n_{d-1} + 1) D(\{n_3,n_4,\ldots,n_{d-1}+1,n_d-1,\ldots\})
   \end{multline} *)

    let rec class_size cl =
      if cl.ext = two || cl.prop = zero then
        one
      else
        IMap.fold (fun d _ s -> class_size_n cl d + s) cl.v (class_size_3 cl)

(* Purely ternary vertices recurse among themselves: *)

    and class_size_3 cl =
      try
        let d' = remove cl three in
        (d'.ext + d'.prop) * class_size d'
      with
      | Zero -> zero
            
(* Vertices of higher degree recurse one step towards lower degrees: *)

    and class_size_n cl d =
      if d > three then begin
        try
          let d' = pred d in
          let cl' = add (remove cl d) d' in
          multiplicity cl' d' * class_size cl'
        with
        | Zero -> zero
      end else
        zero

(* Find all $\{n_3,n_4,\ldots,n_d\}$ with
   \begin{equation}
     E(\{n_3,n_4,\ldots,n_d\}) - 2 = \sum_{i=3}^cl (i-2)n_i  = \ocwlowerid{sum}
   \end{equation}
   The implementation is a variant of [tuples] above. *)

    let rec distribute_degrees' d sum =
      if d < three then
        invalid_arg "distribute_degrees"
      else if d = three then
        [[(d, sum)]]
      else
        distribute_degrees'' d sum (sum / (d - two))

    and distribute_degrees'' d sum n =
      if n < zero then
        []
      else
        List.fold_left (fun ll l -> ((d, n) :: l) :: ll)
          (distribute_degrees'' d sum (pred n))
          (distribute_degrees' (pred d) (sum - (d - two) * n))
          
(* Actually, we need to find all $\{n_3,n_4,\ldots,n_d\}$ with
   \begin{equation}
     E(\{n_3,n_4,\ldots,n_d\}) = \ocwlowerid{sum}
   \end{equation} *)

    let distribute_degrees d sum = distribute_degrees' d (sum - two)

(* Finally we can count all diagrams by adding all possible ways of
   splitting the degrees of vertices. We can also count diagrams where
   \emph{all} degrees satisfy a predicate [f]: *)

    let diagrams ?(f = fun _ -> true) deg n =
      List.fold_left (fun s d ->
        if List.for_all (fun (d', n') -> f d' || n' = zero) d then
          s + class_size (create_class d)
        else
          s)
        zero (distribute_degrees deg n)

(* The next two are duplicated from [ThoList] and [Combinatorics],
   in order to use the specific comparison functions.  *)

    let classify l =
      let rec add_to_class a = function
        | [] -> [of_int 1, a]
        | (n, a') :: rest ->
            if a = a' then
              (succ n, a) :: rest
            else
              (n, a') :: add_to_class a rest
      in
      let rec classify' cl = function
        | [] -> cl
        | a :: rest -> classify' (add_to_class a cl) rest
      in
      classify' [] l

    let permutation_symmetry l =
      List.fold_left (fun s (n, _) -> factorial n * s) one (classify l)

    let symmetry l =
      let sum = List.fold_left (+) zero l in
      if List.exists (fun x -> two * x = sum) l then
	two * permutation_symmetry l
      else
	permutation_symmetry l

(* The number of Feynman diagrams built of vertices with maximum
   degree~$d_{\max}$ in a partition $N_{d,n}=\{n_1,n_2,\ldots,n_d\}$
   with $n = n_1 + n_2 + \cdots + n_d$ and
   \begin{equation}
     \tilde F(d_{\max},N_{d,n}) =
       \frac{n!}{ |\mathcal{S}(N_{d,n})| \sigma(n_d,n)}
       \prod_{i=1}^{d} \frac{F(d_{\max},n_i+1)}{n_i!}
   \end{equation}
   with~$|\mathcal{S}(N)|$ the size of the symmetric group of~$N$,
   $\sigma(n,2n) = 2$ and $\sigma(n,m) = 1$ otherwise. *)

    let keystones p =
      let sum = List.fold_left (+) zero p in
      List.fold_left (fun acc n -> acc / (factorial n)) (factorial sum) p
        / symmetry p
        
    let diagrams_per_keystone deg p =
      List.fold_left (fun acc n -> acc * diagrams deg (succ n)) one p
        
(* We must find
   \begin{equation}
     F(d_{\max},n) =
       \sum_{d=3}^{d_{\max}}
       \sum_{\substack{N = \{n_1,n_2,\ldots,n_d\}\\
                       n_1 + n_2 + \cdots + n_d = n\\
                       1 \le n_1 \le n_2 \le \cdots \le n_d \le \lfloor n/2 \rfloor}}
        \tilde F(d_{\max},N)
   \end{equation} *)

    let diagrams_via_keystones deg n =
      let module N = Nary (struct let max_arity () = to_int (pred deg) end) in
      List.fold_left
        (fun acc p -> acc + diagrams_per_keystone deg p * keystones p)
        zero (List.map (List.map of_int) (N.partitions (to_int n)))

  end

(* \thocwmodulesection{Emulating HELAC} *)

(* In~\cite{HELAC:2000}, one leg is singled out:  *)

module Helac (B : Tuple.Bound) =
  struct
    module Tuple = Tuple.Nary(B)

    type partition = int list
    let inspect_partition p = p

    let partition d sum =
      Partition.tuples d sum 1 (sum - d + 1)

    let rec partitions' d sum =
      let d' = pred d in
      if d' < 2 then
        []
      else
        List.map (fun p -> 1::p) (partition d' (pred sum)) @ partitions' d' sum

    let partitions sum = partitions' (succ (B.max_arity ())) sum

    type 'a children = 'a Tuple.t

    let keystones' l =
      match l with
      | [] -> []
      | head :: tail ->
          [([head],
            ThoList.flatmap (fun p -> Combinatorics.partitions (List.tl p) tail)
              (partitions (List.length l)))]
     
    let keystones l =
      List.map (fun (bra, kets) -> (bra, List.map Tuple.of_list kets))
        (keystones' l)

    let max_subtree n = pred n
  end
    
(* \begin{dubious}
     The following is not tested, but it is no rocket science either \ldots
   \end{dubious} *)

module Helac_Binary =
  struct
    type partition = int * int * int
    let inspect_partition (n1, n2, n3) = [n1; n2; n3]

    let partitions sum =
      List.map (fun (n2, n3) -> (1, n2, n3))
        (Partition.pairs (sum - 1) 1 (sum - 2))

    type 'a children = 'a Tuple.Binary.t

    let keystones' l =
      match l with
      | [] -> []
      | head :: tail ->
          [([head],
            ThoList.flatmap (fun (_, p2, _) -> Combinatorics.split p2 tail)
              (partitions (List.length l)))]
     
    let keystones l =
      List.map (fun (bra, kets) ->
        (bra, List.map (fun (x, y) -> Tuple.Binary.of2 x y) kets))
        (keystones' l)

    let max_subtree n = pred n

  end
    
(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)



