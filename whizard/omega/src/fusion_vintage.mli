(* fusion.mli --

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

    val options : Options.t

(* JRR's implementation of Majoranas needs a special case. *)
    val vintage : bool

(* Wavefunctions are an abstract data type, containing a momentum~[p]
   and additional quantum numbers, collected in~[flavor]. *)
    type wf
    val conjugate : wf -> wf

(* Obviously, [flavor] is not restricted to the physical notion of
   flavor, but can carry spin, color, etc. *)
    type flavor
    val flavor : wf -> flavor
    type flavor_sans_color
    val flavor_sans_color : wf -> flavor_sans_color

(* Momenta are represented by an abstract datatype (defined
   in~[Momentum]) that is optimized for performance.  They can be
   accessed either abstractly or as lists of indices of the external
   momenta.  These indices are assigned sequentially by [amplitude] below. *)
    type p
    val momentum : wf -> p
    val momentum_list : wf -> int list

(* At tree level, the wave functions are uniquely specified by [flavor]
   and momentum.  If loops are included, we need to distinguish among
   orders.  Also, if we build a result from an incomplete sum of diagrams,
   we need to add a distinguishing mark.  At the moment, we assume that a
   [string] that can be attached to the symbol suffices.  *)
    val wf_tag : wf -> string option

(* Coupling constants *)
    type constant

(* and right hand sides of assignments.  The latter are formed from a sign from
   Fermi statistics, a coupling (constand and Lorentz structure) and wave
   functions. *)
    type coupling
    type rhs
    type 'a children
    val sign : rhs -> int
    val coupling : rhs -> constant Coupling.t

    val coupling_tag : rhs -> string option

    type exclusions
    val no_exclusions : exclusions

(* In renormalized perturbation theory, couplings come in different orders
   of the loop expansion.  Be prepared: [val order : rhs -> int] *)

(* \begin{dubious}
     This is here only for the benefit of [Target] and shall become
     [val children : rhs -> wf children] later \ldots
   \end{dubious} *)
    val children : rhs -> wf list

(* Fusions come in two types: fusions of wave functions to off-shell wave
   functions:
   \begin{equation*}
     \phi(p+q) = \phi(p)\phi(q)
   \end{equation*} *)
    type fusion
    val lhs : fusion -> wf
    val rhs : fusion -> rhs list

(* and products at the keystones:
   \begin{equation*}
     \phi(-p-q)\cdot\phi(p)\phi(q)
   \end{equation*} *)
    type braket
    val bra : braket -> wf
    val ket : braket -> rhs list

(* [amplitude goldstones incoming outgoing] calculates the
   amplitude for scattering of [incoming] to [outgoing].  If
   [goldstones] is true, also non-propagating off-shell Goldstone
   amplitudes are included to allow the checking of Slavnov-Taylor
   identities. *)
    type amplitude
    type amplitude_sans_color
    type selectors
    val amplitudes : bool -> exclusions -> selectors ->
      flavor_sans_color list -> flavor_sans_color list -> amplitude list
    val amplitude_sans_color : bool -> exclusions -> selectors ->
      flavor_sans_color list -> flavor_sans_color list -> amplitude_sans_color

    val dependencies : amplitude -> wf -> (wf, coupling) Tree2.t

(* We should be precise regarding the semantics of the following functions, since
   modules implementating [Target] must not make any mistakes interpreting the
   return values.  Instead of calculating the amplitude
   \begin{subequations}
   \begin{equation}
   \label{eq:physical-amplitude}
     \Braket{f_3,p_3,f_4,p_4,\ldots|T|f_1,p_1,f_2,p_2}
   \end{equation}
   directly, O'Mega calculates the---equivalent, but more symmetrical---crossed
   amplitude 
   \begin{equation}
     \Braket{\bar f_1,-p_1,\bar f_2,-p_2,f_3,p_3,f_4,p_4,\ldots|T|0}
   \end{equation}
   Internally, all flavors are represented by their charge conjugates
   \begin{equation}
   \label{eq:internal-amplitude}
     A(f_1,-p_1,f_2,-p_2,\bar f_3,p_3,\bar f_4,p_4,\ldots)
   \end{equation}
   \end{subequations}
   The correspondence of vertex and term in the lagrangian
   \begin{equation}
     \parbox{26\unitlength}{%
       \fmfframe(5,3)(5,3){%
         \begin{fmfgraph*}(15,20)
           \fmfleft{v}
           \fmfright{p,A,e}
           \fmflabel{$\mathrm{e}^-$}{e}
           \fmflabel{$\mathrm{e}^+$}{p}
           \fmflabel{$\mathrm{A}$}{A}
           \fmf{fermion}{p,v,e}
           \fmf{photon}{A,v}
           \fmfdot{v}
         \end{fmfgraph*}}}: \bar\psi\fmslash{A}\psi
   \end{equation}
   suggests to denote the \emph{outgoing} particle by the flavor of the
   \emph{anti}particle and the \emph{outgoing} \emph{anti}particle by the
   flavor of the particle, since this choice allows to represent the vertex
   by a triple
   \begin{equation}
     \bar\psi\fmslash{A}\psi: (\mathrm{e}^+,A,\mathrm{e}^-)
   \end{equation}
   which is more intuitive than the alternative $(\mathrm{e}^-,A,\mathrm{e}^+)$.
   Also, when thinking in terms of building wavefunctions from the outside in,
   the outgoing \emph{antiparticle} is represented by a \emph{particle}
   propagator and vice versa\footnote{Even if this choice will appear slightly
   counter-intuitive on the [Target] side, one must keep in mind that much more
   people are expected to prepare [Model]s.}.
   [incoming] and [outgoing] are the physical flavors as
   in~(\ref{eq:physical-amplitude}) *)
    val incoming : amplitude -> flavor list
    val outgoing : amplitude -> flavor list

(* [externals] are flavors and momenta as in~(\ref{eq:internal-amplitude}) *)
    val externals : amplitude -> wf list

    val variables : amplitude -> wf list
    val fusions : amplitude -> fusion list
    val brakets : amplitude -> braket list
    val on_shell : amplitude -> (wf -> bool)
    val is_gauss : amplitude -> (wf -> bool)
    val constraints : amplitude -> string option
    val symmetry : amplitude -> int

    val allowed : amplitude -> bool

(*i
(* \thocwmodulesubsection{Performance Hacks} *)

    val initialize_cache : string -> unit
    val set_cache_name : string -> unit
i*)

(* \thocwmodulesubsection{Diagnostics} *)

    val check_charges : unit -> flavor_sans_color list list
    val count_fusions : amplitude -> int
    val count_propagators : amplitude -> int
    val count_diagrams : amplitude -> int

    val forest : wf -> amplitude -> ((wf * coupling option, wf) Tree.t) list
    val poles : amplitude -> wf list list
    val s_channel : amplitude -> wf list

    val tower_to_dot : out_channel -> amplitude -> unit
    val amplitude_to_dot : out_channel -> amplitude -> unit

(* \thocwmodulesubsection{WHIZARD} *)

    val phase_space_channels : out_channel -> amplitude_sans_color -> unit
    val phase_space_channels_flipped : out_channel -> amplitude_sans_color -> unit

  end

(* There is more than one way to make fusions.  *)

module type Maker =
    functor (P : Momentum.T) -> functor (M : Model.T) ->
      T with type p = P.t
      and type flavor = Colorize.It(M).flavor
      and type flavor_sans_color = M.flavor
      and type constant = M.constant
      and type selectors = Cascade.Make(M)(P).selectors

(*i If we want or need to expose [Make], here's how to do it:

module type Stat =
  sig
    type flavor
    type stat
    exception Impossible
    val stat : flavor -> int -> stat
    val stat_fuse : stat -> stat -> flavor -> stat
    val stat_sign : stat -> int
  end

module type Stat_Maker = functor (M : Model.T) ->
  Stat with type flavor = M.flavor

module Make : functor (PT : Tuple.Poly) (Stat : Stat_Maker)
                      (T : Topology.T with type 'a children = 'a PT.t) -> Maker

i*)

(* Straightforward Dirac fermions vs. slightly more complicated
   Majorana fermions: *)

module Binary : Maker
module Binary_Majorana : Maker

module Mixed23 : Maker
module Mixed23_Majorana : Maker

module Nary : functor (B : Tuple.Bound) -> Maker
module Nary_Majorana : functor (B : Tuple.Bound) -> Maker

(* We can also proceed \'a la~\cite{HELAC:2000}.  Empirically,
   this will use slightly~($O(10\%)$) fewer fusions than the
   symmetric factorization.  Our implementation uses
   significantly~($O(50\%)$) fewer fusions than reported
   by~\cite{HELAC:2000}.  Our pruning of the DAG might
   be responsible for this.  *)

module Helac : functor (B : Tuple.Bound) -> Maker
module Helac_Majorana : functor (B : Tuple.Bound) -> Maker

(* \thocwmodulesection{Multiple Amplitudes} *)

module type Multi =
  sig
    exception Mismatch
    val options : Options.t

    type flavor
    type process = flavor list * flavor list
    type amplitude
    type fusion
    type wf
    type exclusions
    val no_exclusions : exclusions
    type selectors
    type amplitudes

    (* Construct all possible color flow amplitudes for a given process. *)
    val amplitudes : bool -> int option ->
      exclusions -> selectors -> process list -> amplitudes
    val empty : amplitudes

(*i
    (* Precompute the vertex table cache. *)
    val initialize_cache : string -> unit
    val set_cache_name : string -> unit
i*)

    (* The list of all combinations of incoming and outgoing particles
       with a nonvanishing scattering amplitude. *)
    val flavors : amplitudes -> process list

    (* The list of all combinations of incoming and outgoing particles that
       don't lead to any color flow with non vanishing scattering amplitude. *)
    val vanishing_flavors : amplitudes -> process list

    (* The list of all color flows with a nonvanishing scattering amplitude. *)
    val color_flows : amplitudes -> Color.Flow.t list

    (* The list of all valid helicity combinations. *)
    val helicities : amplitudes -> (int list * int list) list

    (* The list of all amplitudes. *)
    val processes : amplitudes -> amplitude list

    (* [(process_table a).(f).(c)] returns the amplitude for the [f]th
       allowed flavor combination and the [c]th allowed color flow as
       an [amplitude option]. *)
    val process_table : amplitudes -> amplitude option array array

    (* The list of all non redundant fusions together with the amplitudes
       they came from. *)
    val fusions : amplitudes -> (fusion * amplitude) list

    (* If there's more than external flavor state, the wavefunctions are
       \emph{not} uniquely specified by [flavor] and [Momentum.t].  This
       function can be used to determine how many variables must be allocated. *)
    val multiplicity : amplitudes -> wf -> int

    (* This function can be used to disambiguate wavefunctions with the same
       combination of [flavor] and [Momentum.t]. *)
    val dictionary : amplitudes -> amplitude -> wf -> int

    (* [(color_factors a).(c1).(c2)] power of~$N_C$ for the given product
       of color flows. *)
    val color_factors : amplitudes -> Color.Flow.factor array array

    (* A description of optional diagram selectors. *)
    val constraints : amplitudes -> string option

  end

module type Multi_Maker = functor (Fusion_Maker : Maker) ->
  functor (P : Momentum.T) ->
    functor (M : Model.T) ->
      Multi with type flavor = M.flavor
      and type amplitude = Fusion_Maker(P)(M).amplitude
      and type fusion = Fusion_Maker(P)(M).fusion
      and type wf = Fusion_Maker(P)(M).wf
      and type selectors = Fusion_Maker(P)(M).selectors

module Multi : Multi_Maker

(* \thocwmodulesection{Tags} *)

(* It appears that there are useful applications for tagging couplings
   and wave functions, e.\,g.~skeleton expansion and diagram selections.
   We can abstract this in a [Tags] signature: *)

module type Tags =
  sig
    type wf
    type coupling
    type 'a children
    val null_wf : wf
    val null_coupling : coupling
    val fuse : coupling -> wf children -> wf
    val wf_to_string : wf -> string option
    val coupling_to_string : coupling -> string option
  end

module type Tagger =
    functor (PT : Tuple.Poly) -> Tags with type 'a children = 'a PT.t

module type Tagged_Maker =
    functor (Tagger : Tagger) ->
      functor (P : Momentum.T) -> functor (M : Model.T) ->
        T with type p = P.t
        and type flavor = Colorize.It(M).flavor
        and type flavor_sans_color = M.flavor
        and type constant = M.constant

module Tagged_Binary : Tagged_Maker

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
