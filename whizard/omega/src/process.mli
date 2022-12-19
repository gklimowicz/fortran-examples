(* process.mli --

   Copyright (C) 1999-2022 by

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

module type T =
  sig

    type flavor

(* \begin{dubious}
     Eventually this should become an abstract type:
   \end{dubious} *)
    type t = flavor list * flavor list

    val incoming : t -> flavor list
    val outgoing : t -> flavor list

(* [parse_decay s] decodes a decay description ["a -> b c ..."], where
    each word is split into a bag of flavors separated by [':']s. *)
    type decay
    val parse_decay : string -> decay
    val expand_decays : decay list -> t list

(* [parse_scattering s] decodes a scattering description ["a b -> c d ..."],
    where each word is split into a bag of flavors separated by [':']s. *)
    type scattering
    val parse_scattering : string -> scattering
    val expand_scatterings : scattering list -> t list

(* [parse_process s] decodes process descriptions
   \begin{subequations}
   \begin{align}
     \text{\texttt{"a b c d"}} &\Rightarrow \text{[Any [a; b; c; d]]} \\
     \text{\texttt{"a -> b c d"}} &\Rightarrow \text{[Decay (a, [b; c; d])]} \\
     \text{\texttt{"a b -> c d"}} &\Rightarrow \text{[Scattering (a, b, [c; d])]}
   \end{align}
   \end{subequations}
   where each word is split into a bag of flavors separated by `\texttt{:}'s. *)
    type any
    type process = Any of any | Decay of decay | Scattering of scattering
    val parse_process : string -> process

(* [remove_duplicate_final_states partition processes] removes duplicates from
   [processes], which differ only by a permutation of final state particles.
   The permutation must respect the partitioning given by the offset 1 integers
   in [partition]. *)
    val remove_duplicate_final_states : int list list -> t list -> t list

(* [diff set1 set2] returns the processes in [set1] with the processes in [set2]
   removed.  [set2] does not need to be a subset of [set1]. *)
    val diff : t list -> t list -> t list

(* \begin{dubious}
     Not functional yet.  Interface subject to change.  Should be moved to
     [Fusion.Multi], because we will want to cross \emph{colored} matrix
     elements.
   \end{dubious} *)

(* Factor amplitudes that are related by crossing symmetry. *)

    val crossing : t list -> (flavor list * int list * t) list

  end

module Make (M : Model.T) : T with type flavor = M.flavor

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
