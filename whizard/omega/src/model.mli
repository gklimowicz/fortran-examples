(* model.mli --

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

(* \thocwmodulesection{General Quantum Field Theories} *)

module type T =
  sig

(* [flavor] abstractly encodes all quantum numbers. *) 
    type flavor

(* [Color.t] encodes the ($\textrm{SU}(N)$) color representation. *) 
    val color : flavor -> Color.t
    val nc : unit -> int

(* The set of conserved charges. *)
    module Ch : Charges.T
    val charges : flavor -> Ch.t

(* The PDG particle code for interfacing with Monte Carlos. *)
    val pdg : flavor -> int

(* The Lorentz representation of the particle. *)
    val lorentz : flavor -> Coupling.lorentz

(* The propagator for the particle, which \emph{can} depend
   on a gauge parameter. *)
    type gauge
    val propagator : flavor -> gauge Coupling.propagator

(* \emph{Not} the symbol for the numerical value, but the
   scheme or strategy.  *)
    val width : flavor -> Coupling.width

(* Charge conjugation, with and without color.  *)
    val conjugate : flavor -> flavor

(* Returns $1$ for fermions, $-1$ for anti-fermions, $2$ for Majoranas
   and $0$ otherwise.  *)
    val fermion : flavor -> int

(* The Feynman rules.  [vertices] and [(fuse2, fuse3, fusen)] are
   redundant, of course.  However, [vertices] is required for building
   functors for models and [vertices] can be recovered from
   [(fuse2, fuse3, fusen)] only at great cost. *)

(* \begin{dubious}
     Nevertheless: [vertices] is a candidate for removal, b/c we can
     build a smarter [Colorize] functor acting on [(fuse2, fuse3, fusen)].
     It can support an arbitrary numer of color lines.  But we have to test
     whether it is efficient enough.  And we have to make sure that this
     wouldn't break the UFO interface.
   \end{dubious} *)
    type constant 

    (* Later: [type orders] to count orders of couplings *)

    val max_degree : unit -> int
    val vertices : unit ->
      ((((flavor * flavor * flavor) * constant Coupling.vertex3 * constant) list)
         * (((flavor * flavor * flavor * flavor) * constant Coupling.vertex4 * constant) list)
         * (((flavor list) * constant Coupling.vertexn * constant) list))
    val fuse2 : flavor -> flavor -> (flavor * constant Coupling.t) list
    val fuse3 : flavor -> flavor -> flavor -> (flavor * constant Coupling.t) list
    val fuse : flavor list -> (flavor * constant Coupling.t) list

    (* Later: [val orders : constant -> orders] counting orders of couplings *)

(* The list of all known flavors. *)
    val flavors : unit -> flavor list

(* The flavors that can appear in incoming or outgoing states, grouped
   in a way that is useful for user interfaces. *)
    val external_flavors : unit -> (string * flavor list) list

(* The Goldstone bosons corresponding to a gauge field, if any. *)
    val goldstone : flavor -> (flavor * constant Coupling.expr) option

(* The dependent parameters. *)
    val parameters : unit -> constant Coupling.parameters
        
(* Translate from and to convenient textual representations of flavors.  *)
    val flavor_of_string : string -> flavor
    val flavor_to_string : flavor -> string

(* \TeX{} and \LaTeX{} *) 
    val flavor_to_TeX : flavor -> string

(* The following must return unique symbols that are acceptable as
   symbols in all programming languages under consideration as targets.
   Strings of alphanumeric characters (starting with a letter) should
   be safe.  Underscores are also usable, but would violate strict
   Fortran77. *)
    val flavor_symbol : flavor -> string
    val gauge_symbol : gauge -> string
    val mass_symbol : flavor -> string
    val width_symbol : flavor -> string
    val constant_symbol : constant -> string

(* Model specific options. *)
    val options : Options.t

(* \textit{Not ready for prime time} or other warnings to
   be written to the source files for the amplitudes. *)

    val caveats : unit -> string list

  end

(* In addition to hardcoded models, we can have models that are
   initialized at run time. *)

(* \thocwmodulesection{Mutable Quantum Field Theories} *)

module type Mutable =
  sig
    include T

    val init : unit -> unit

(* Export only one big initialization function to discourage
   partial initializations.  Labels make this usable. *)

    val setup :
        color:(flavor -> Color.t) ->
        nc:(unit -> int) ->
        pdg:(flavor -> int) ->
        lorentz:(flavor -> Coupling.lorentz) ->
        propagator:(flavor -> gauge Coupling.propagator) ->
        width:(flavor -> Coupling.width) ->
        goldstone:(flavor -> (flavor * constant Coupling.expr) option) ->
        conjugate:(flavor -> flavor) ->
        fermion:(flavor -> int) ->
        vertices:
          (unit ->
           ((((flavor * flavor * flavor) * constant Coupling.vertex3 * constant) list)
            * (((flavor * flavor * flavor * flavor) * constant Coupling.vertex4 * constant) list)
            * (((flavor list) * constant Coupling.vertexn * constant) list))) ->
        flavors:((string * flavor list) list) ->
        parameters:(unit -> constant Coupling.parameters) ->
        flavor_of_string:(string -> flavor) ->
        flavor_to_string:(flavor -> string) ->
        flavor_to_TeX:(flavor -> string) ->
        flavor_symbol:(flavor -> string) ->
        gauge_symbol:(gauge -> string) ->
        mass_symbol:(flavor -> string) ->
        width_symbol:(flavor -> string) ->
        constant_symbol:(constant -> string) ->
        unit
  end

(* \thocwmodulesection{Gauge Field Theories} *)

(* The following signatures are used only for model building.  The diagrammatics
   and numerics is supposed to be completely ignorant about the detail of the
   models and expected to rely on the interface [T] exclusively.
   \begin{dubious}
     In the end, we might have functors [(M : T) -> Gauge], but we will
     need to add the quantum numbers to [T].
   \end{dubious} *)

module type Gauge =
  sig
    include T

(* Matter field carry conserved quantum numbers and can be replicated
   in generations without changing the gauge sector.  *)
    type matter_field

(* Gauge bosons proper.  *)
    type gauge_boson

(* Higgses, Goldstones and all the rest:  *)
    type other

(* We can query the kind of field *)
    type field =
      | Matter of matter_field
      | Gauge of gauge_boson
      | Other of other
    val field : flavor -> field

(* and we can build new fields of a given kind: *)
    val matter_field : matter_field -> flavor
    val gauge_boson : gauge_boson -> flavor
    val other : other -> flavor
  end

(* \thocwmodulesection{Gauge Field Theories with Broken Gauge Symmetries} *)

(* Both are carefully crafted as subtypes of [Gauge] so that
   they can be used in place of [Gauge] and [T] everywhere: *)

module type Broken_Gauge =
  sig
    include Gauge

    type massless
    type massive
    type goldstone

    type kind =
      | Massless of massless
      | Massive of massive
      | Goldstone of goldstone
    val kind : gauge_boson -> kind

    val massless : massive -> gauge_boson
    val massive : massive -> gauge_boson
    val goldstone : goldstone -> gauge_boson

  end

module type Unitarity_Gauge =
  sig
    include Gauge

    type massless
    type massive

    type kind =
      | Massless of massless
      | Massive of massive
    val kind : gauge_boson -> kind

    val massless : massive -> gauge_boson
    val massive : massive -> gauge_boson

  end

module type Colorized =
  sig

    include T

    type flavor_sans_color
    val flavor_sans_color : flavor -> flavor_sans_color
    val conjugate_sans_color : flavor_sans_color -> flavor_sans_color

    val amplitude : flavor_sans_color list -> flavor_sans_color list ->
      (flavor list * flavor list) list
    val flow : flavor list -> flavor list -> Color.Flow.t

  end

module type Colorized_Gauge =
  sig

    include Gauge

    type flavor_sans_color
    val flavor_sans_color : flavor -> flavor_sans_color
    val conjugate_sans_color : flavor_sans_color -> flavor_sans_color

    val amplitude : flavor_sans_color list -> flavor_sans_color list ->
      (flavor list * flavor list) list
    val flow : flavor list -> flavor list -> Color.Flow.t

  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
