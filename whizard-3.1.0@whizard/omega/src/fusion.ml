(* fusion.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>
       Marco Sekulla <marco.sekulla@kit.edu>

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

(* Avoid refering to [Pervasives.compare], because [Pervasives] will
   become [Stdlib.Pervasives] in O'Caml 4.07 and [Stdlib] in O'Caml 4.08. *)
let pcompare = compare

module type T =
  sig
    val options : Options.t
    val vintage : bool
    type wf
    val conjugate : wf -> wf
    type flavor
    type flavor_sans_color
    val flavor : wf -> flavor
    val flavor_sans_color : wf -> flavor_sans_color
    type p
    val momentum : wf -> p
    val momentum_list : wf -> int list
    val wf_tag : wf -> string option
    type constant
    type coupling
    type rhs
    type 'a children
    val sign : rhs -> int
    val coupling : rhs -> constant Coupling.t
    val coupling_tag : rhs -> string option
    type exclusions
    val no_exclusions : exclusions
    val children : rhs -> wf list
    type fusion
    val lhs : fusion -> wf
    val rhs : fusion -> rhs list
    type braket
    val bra : braket -> wf
    val ket : braket -> rhs list
    type amplitude
    type amplitude_sans_color
    type selectors
    val amplitudes : bool -> exclusions -> selectors ->
      flavor_sans_color list -> flavor_sans_color list -> amplitude list
    val amplitude_sans_color : bool -> exclusions -> selectors ->
      flavor_sans_color list -> flavor_sans_color list -> amplitude_sans_color
    val dependencies : amplitude -> wf -> (wf, coupling) Tree2.t
    val incoming : amplitude -> flavor list
    val outgoing : amplitude -> flavor list
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
    val initialize_cache : string -> unit
    val set_cache_name : string -> unit
i*)
    val check_charges : unit -> flavor_sans_color list list
    val count_fusions : amplitude -> int
    val count_propagators : amplitude -> int
    val count_diagrams : amplitude -> int
    val forest : wf -> amplitude -> ((wf * coupling option, wf) Tree.t) list
    val poles : amplitude -> wf list list
    val s_channel : amplitude -> wf list
    val tower_to_dot : out_channel -> amplitude -> unit
    val amplitude_to_dot : out_channel -> amplitude -> unit
    val phase_space_channels : out_channel -> amplitude_sans_color -> unit
    val phase_space_channels_flipped : out_channel -> amplitude_sans_color -> unit
  end

module type Maker =
    functor (P : Momentum.T) -> functor (M : Model.T) ->
      T with type p = P.t
      and type flavor = Colorize.It(M).flavor
      and type flavor_sans_color = M.flavor
      and type constant = M.constant
      and type selectors = Cascade.Make(M)(P).selectors

(* \thocwmodulesection{Fermi Statistics} *)

module type Stat =
  sig

    (* This will be [Model.T.flavor]. *)
    type flavor

    (* A record of the fermion lines in the 1POW. *)
    type stat

    (* Vertices with an odd number of fermion fields. *)
    exception Impossible

    (* External lines. *)
    val stat : flavor -> int -> stat

    (* [stat_fuse (Some flines) slist f] combines the fermion lines
       in the elements of [slist] according to the connections listed
       in [flines].
       On the other hand, [stat_fuse None slist f] corresponds to
       the legacy mode with \emph{at most} two fermions.
       The resulting flavor [f] of the 1POW can be ignored for models
       with only Dirac fermions, except for debugging, since
       the direction of the arrows is unambiguous.
       However, in the case of Majorana fermions and/or fermion number
       violating interactions, the flavor [f] must be used. *)
    val stat_fuse :
      Coupling.fermion_lines option -> stat list -> flavor -> stat

    (* Analogous to [stat_fuse], but for the finalizing keystone
       instead of the 1POW.  *) 
    val stat_keystone :
      Coupling.fermion_lines option -> stat list -> flavor -> stat

    (* Compute the sign corresponding to the fermion lines in
       a 1POW or keystone. *)
    val stat_sign : stat -> int

    (* Debugging and consistency checks \ldots *)
    val stat_to_string : stat -> string
    val equal : stat -> stat -> bool
    val saturated : stat -> bool

end

module type Stat_Maker = functor (M : Model.T) ->
  Stat with type flavor = M.flavor

(* \thocwmodulesection{Dirac Fermions} *)

let dirac_log silent logging = logging
let dirac_log silent logging = silent

exception Majorana

module Stat_Dirac (M : Model.T) : (Stat with type flavor = M.flavor) =
  struct 
    type flavor = M.flavor

(* \begin{equation}
     \gamma_\mu\psi(1)\,G^{\mu\nu}\,\bar\psi(2)\gamma_\nu\psi(3)
         - \gamma_\mu\psi(3)\,G^{\mu\nu}\,\bar\psi(2)\gamma_\nu\psi(1)
   \end{equation} *)

    type stat =
      | Fermion of int * (int option * int option) list
      | AntiFermion of int * (int option * int option) list
      | Boson of (int option * int option) list

    let lines_to_string lines =
      ThoList.to_string
        (function
         | Some i, Some j -> Printf.sprintf "%d>%d" i j
         | Some i, None -> Printf.sprintf "%d>*" i
         | None, Some j -> Printf.sprintf "*>%d" j
         | None, None -> "*>*")
        lines

    let stat_to_string = function
      | Boson lines -> Printf.sprintf "Boson %s" (lines_to_string lines)
      | Fermion (p, lines) ->
         Printf.sprintf "Fermion (%d, %s)" p (lines_to_string lines)
      | AntiFermion (p, lines) ->
         Printf.sprintf "AntiFermion (%d, %s)" p (lines_to_string lines)

    let equal s1 s2 =
      match s1, s2 with
      | Boson l1, Boson l2 ->
         List.sort compare l1 = List.sort compare l2
      | Fermion (p1, l1), Fermion (p2, l2)
      | AntiFermion (p1, l1), AntiFermion (p2, l2) ->
         p1 = p2 && List.sort compare l1 = List.sort compare l2
      | _ -> false

    let saturated = function
      | Boson _ -> true
      | _ -> false

    let stat f p =
      match M.fermion f with
      | 0 -> Boson []
      | 1 -> Fermion (p, [])
      | -1 -> AntiFermion (p, [])
      | 2 -> raise Majorana
      | _ -> invalid_arg "Fusion.Stat_Dirac: invalid fermion number"

    exception Impossible

    let stat_fuse_pair_legacy f s1 s2 =
      match s1, s2 with
      | Boson l1, Boson l2 -> Boson (l1 @ l2)
      | Boson l1, Fermion (p, l2) -> Fermion (p, l1 @ l2)
      | Boson l1, AntiFermion (p, l2) -> AntiFermion (p, l1 @ l2)
      | Fermion (p, l1), Boson l2 -> Fermion (p, l1 @ l2)
      | AntiFermion (p, l1), Boson l2 -> AntiFermion (p, l1 @ l2)
      | AntiFermion (pbar, l1), Fermion (p, l2) ->
          Boson ((Some pbar, Some p) :: l1 @ l2)
      | Fermion (p, l1), AntiFermion (pbar, l2) ->
          Boson ((Some pbar, Some p) :: l1 @ l2)
      | Fermion _, Fermion _ | AntiFermion _, AntiFermion _ ->
          raise Impossible

    let stat_fuse_legacy s1 s23__n f =
      List.fold_right (stat_fuse_pair_legacy f) s23__n s1

    let stat_fuse_legacy_logging s1 s23__n f =
      let s = stat_fuse_legacy s1 s23__n f in
      Printf.eprintf
        "stat_fuse_legacy: %s <- %s -> %s\n"
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string (s1 :: s23__n))
        (stat_to_string s);
      s

    let stat_fuse_legacy =
      dirac_log stat_fuse_legacy stat_fuse_legacy_logging

    module IMap = Map.Make (struct type t = int let compare = compare end)

    type partial =
      { stat : stat (* the [stat] accumulated so far *);
        fermions : int IMap.t (* a map from the indices in the vertex to open fermion lines *);
        antifermions : int IMap.t (* a map from the indices in the vertex to open antifermion lines *);
        n : int (* the number of incoming propagators *) }

    let partial_to_string p =
      Printf.sprintf
        "{ fermions=%s, antifermions=%s, state=%s, #=%d }"
        (ThoList.to_string
           (fun (i, f) -> Printf.sprintf "%d@%d" f i)
           (IMap.bindings p.fermions))
        (ThoList.to_string
           (fun (i, f) -> Printf.sprintf "%d@%d" f i)
           (IMap.bindings p.antifermions))
        (stat_to_string p.stat)
        p.n

    let add_lines l = function
      | Boson l' -> Boson (List.rev_append l l')
      | Fermion (n, l') -> Fermion (n, List.rev_append l l')
      | AntiFermion (n, l') -> AntiFermion (n, List.rev_append l l')

    let partial_of_slist slist =
      List.fold_left
        (fun acc s ->
          let n = succ acc.n in
          match s with
          | Boson l ->
             { acc with
               stat = add_lines l acc.stat;
               n }
          | Fermion (p, l) ->
             { acc with
               fermions = IMap.add n p acc.fermions;
               stat = add_lines l acc.stat;
               n }
          | AntiFermion (p, l) ->
             { acc with
               antifermions = IMap.add n p acc.antifermions;
               stat = add_lines l acc.stat;
               n } )
        { stat = Boson [];
          fermions = IMap.empty;
          antifermions = IMap.empty;
          n = 0 }
        slist

    let find_opt p map =
      try Some (IMap.find p map) with Not_found -> None

    let match_fermion_line p (i, j) =
      if i <= p.n && j <= p.n then
        match find_opt i p.fermions, find_opt j p.antifermions with
        | (Some _ as f), (Some _ as fbar) ->
           { p with
             stat = add_lines [fbar, f] p.stat;
             fermions = IMap.remove i p.fermions;
             antifermions = IMap.remove j p.antifermions }
        | _ ->
           invalid_arg "match_fermion_line: mismatched boson"
      else if i <= p.n then
        match find_opt i p.fermions, p.stat with
        | Some f, Boson l ->
           { p with
             stat = Fermion (f, l);
             fermions = IMap.remove i p.fermions }
        | _ ->
           invalid_arg "match_fermion_line: mismatched fermion"
      else if j <= p.n then
        match find_opt j p.antifermions, p.stat with
        | Some fbar, Boson l ->
           { p with
             stat = AntiFermion (fbar, l);
             antifermions = IMap.remove j p.antifermions }
        | _ ->
           invalid_arg "match_fermion_line: mismatched antifermion"
      else
        failwith "match_fermion_line: impossible"

    let match_fermion_line_logging p (i, j) =
      Printf.eprintf
        "match_fermion_line %s (%d, %d)"
        (partial_to_string p) i j;
      let p' = match_fermion_line p (i, j) in
      Printf.eprintf " >> %s\n" (partial_to_string p');
      p'

    let match_fermion_line =
      dirac_log match_fermion_line match_fermion_line_logging

    let match_fermion_lines flines s1 s23__n =
      let p = partial_of_slist (s1 :: s23__n) in
      List.fold_left match_fermion_line p flines

    let stat_fuse_new flines s1 s23__n f =
      (match_fermion_lines flines s1 s23__n).stat

    let stat_fuse_new_checking flines s1 s23__n f =
      let stat = stat_fuse_new flines s1 s23__n f in
      if List.length flines < 2 then
        begin
          let legacy = stat_fuse_legacy s1 s23__n f in
          if not (equal stat legacy) then
            failwith
              (Printf.sprintf
                 "Fusion.Stat_Dirac.stat_fuse_new: %s <> %s!"
                 (stat_to_string stat)
                 (stat_to_string legacy))
        end;
      stat

    let stat_fuse_new_logging flines s1 s23__n f =
      Printf.eprintf
        "stat_fuse_new: connecting fermion lines %s in %s <- %s\n"
        (UFO_Lorentz.fermion_lines_to_string flines)
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string (s1 :: s23__n));
      stat_fuse_new_checking flines s1 s23__n f

    let stat_fuse_new =
      dirac_log stat_fuse_new stat_fuse_new_logging

    let stat_fuse flines_opt slist f =
      match slist with
      | [] -> invalid_arg "Fusion.Stat_Dirac.stat_fuse: empty"
      | s1 :: s23__n ->
         begin match flines_opt with
         | Some flines -> stat_fuse_new flines s1 s23__n f
         | None -> stat_fuse_legacy s1 s23__n f
         end

    let stat_fuse_logging flines_opt slist f =
      Printf.eprintf
        "stat_fuse: %s <- %s\n"
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string slist);
      stat_fuse flines_opt slist f

    let stat_fuse =
      dirac_log stat_fuse stat_fuse_logging

    let stat_keystone_legacy s1 s23__n f =
      let s2 = List.hd s23__n
      and s34__n = List.tl s23__n in
      stat_fuse_legacy s1 [stat_fuse_legacy s2 s34__n (M.conjugate f)] f

    let stat_keystone_legacy_logging s1 s23__n f =
      let s = stat_keystone_legacy s1 s23__n f in
      Printf.eprintf
        "stat_keystone_legacy: %s (%s) %s -> %s\n"
        (stat_to_string s1)
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string s23__n)
        (stat_to_string s);
      s

    let stat_keystone_legacy =
      dirac_log stat_keystone_legacy stat_keystone_legacy_logging

    let stat_keystone flines_opt slist f =
      match slist with
      | [] -> invalid_arg "Fusion.Stat_Dirac.stat_keystone: empty"
      | [s] -> invalid_arg "Fusion.Stat_Dirac.stat_keystone: singleton"
      | s1 :: (s2 :: s34__n as s23__n) ->
         begin match flines_opt with
         | None -> stat_keystone_legacy s1 s23__n f
         | Some flines ->
            (* The fermion line indices in [flines] must match
               the lines on one side of the keystone. *)
            let stat =
              stat_fuse_legacy s1 [stat_fuse_new flines s2 s34__n f] f in
            if saturated stat then
              stat
            else
              failwith
                (Printf.sprintf
                   "Fusion.Stat_Dirac.stat_keystone: incomplete %s!"
                   (stat_to_string stat))
         end

    let stat_keystone_logging flines_opt slist f =
      let s = stat_keystone flines_opt slist f in
      Printf.eprintf
        "stat_keystone:        %s (%s) %s -> %s\n"
        (stat_to_string (List.hd slist))
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string (List.tl slist))
        (stat_to_string s);
      s

    let stat_keystone =
      dirac_log stat_keystone stat_keystone_logging

(* \begin{figure}
     \begin{displaymath}
       \parbox{26\unitlength}{%
         \begin{fmfgraph*}(25,15)
           \fmfstraight
           \fmfleft{f}
           \fmfright{f1,f2,f3}
           \fmflabel{$\psi(1)$}{f1}
           \fmflabel{$\bar\psi(2)$}{f2}
           \fmflabel{$\psi(3)$}{f3}
           \fmflabel{$0$}{f}
           \fmf{fermion}{f1,v1,f}
           \fmffreeze
           \fmf{fermion,tension=0.5}{f3,v2,f2}
           \fmf{photon}{v1,v2}
           \fmfdot{v1,v2}
         \end{fmfgraph*}}
       \qquad\qquad-\qquad
       \parbox{26\unitlength}{%
         \begin{fmfgraph*}(25,15)
           \fmfstraight
           \fmfleft{f}
           \fmfright{f1,f2,f3}
           \fmflabel{$\psi(1)$}{f1}
           \fmflabel{$\bar\psi(2)$}{f2}
           \fmflabel{$\psi(3)$}{f3}
           \fmflabel{$0$}{f}
           \fmf{fermion}{f3,v1,f}
           \fmffreeze
           \fmf{fermion,tension=0.5}{f1,v2,f2}
           \fmf{photon}{v1,v2}
           \fmfdot{v1,v2}
         \end{fmfgraph*}}
     \end{displaymath} 
     \caption{\label{fig:stat_fuse} Relative sign from Fermi statistics.}
   \end{figure} *)

(* \begin{equation}
     \epsilon \left(\left\{ (0,1), (2,3) \right\}\right)
       = - \epsilon \left(\left\{ (0,3), (2,1) \right\}\right)
   \end{equation} *)

    let permutation lines =
      let fout, fin = List.split lines in
      let eps_in, _ = Combinatorics.sort_signed fin
      and eps_out, _ = Combinatorics.sort_signed fout in
      (eps_in * eps_out)

(* \begin{dubious}
     This comparing of permutations of fermion lines is a bit tedious
     and takes a macroscopic fraction of time.  However, it's less than
     20\,\%, so we don't focus on improving on it yet.
   \end{dubious} *)

    let stat_sign = function
      | Boson lines -> permutation lines
      | Fermion (p, lines) -> permutation ((None, Some p) :: lines)
      | AntiFermion (pbar, lines) -> permutation ((Some pbar, None) :: lines)

  end

(* \thocwmodulesection{Tags} *)

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

(* No tags is one option for good tags \ldots *)

module No_Tags (PT : Tuple.Poly) =
  struct
    type wf = unit
    type coupling = unit
    type 'a children = 'a PT.t
    let null_wf = ()
    let null_coupling = ()
    let fuse () _ = ()
    let wf_to_string () = None
    let coupling_to_string () = None
  end

(* \begin{dubious}
     Here's a simple additive tag that can grow into something useful
     for loop calculations.
   \end{dubious} *)

module Loop_Tags (PT : Tuple.Poly) =
  struct
    type wf = int
    type coupling = int
    type 'a children = 'a PT.t
    let null_wf = 0
    let null_coupling = 0
    let fuse c wfs = PT.fold_left (+) c wfs
    let wf_to_string n = Some (string_of_int n)
    let coupling_to_string n = Some (string_of_int n)
  end

module Order_Tags (PT : Tuple.Poly) =
  struct
    type wf = int
    type coupling = int
    type 'a children = 'a PT.t
    let null_wf = 0
    let null_coupling = 0
    let fuse c wfs = PT.fold_left (+) c wfs
    let wf_to_string n = Some (string_of_int n)
    let coupling_to_string n = Some (string_of_int n)
  end
    
(* \thocwmodulesection{[Tagged], the [Fusion.Make] Functor} *)

module Tagged (Tagger : Tagger) (PT : Tuple.Poly)
    (Stat : Stat_Maker) (T : Topology.T with type 'a children = 'a PT.t)
    (P : Momentum.T) (M : Model.T) =
  struct 

    let vintage = false

    type cache_mode = Cache_Use | Cache_Ignore | Cache_Overwrite
    let cache_option = ref Cache_Ignore
    type qcd_order = 
      | QCD_order of int
    type ew_order = 
      | EW_order of int
    let qcd_order = ref (QCD_order 99)
    let ew_order = ref (EW_order 99)

    let options = Options.create
        [
(*i
          "ignore-cache", Arg.Unit (fun () -> cache_option := Cache_Ignore),
          " ignore cached model tables (default)";
          "use-cache", Arg.Unit (fun () -> cache_option := Cache_Use),
          " use cached model tables";
          "overwrite-cache", Arg.Unit (fun () -> cache_option := Cache_Overwrite),
          " overwrite cached model tables";
i*)
	  "qcd", Arg.Int (fun n -> qcd_order := QCD_order n), 
	  " set QCD order n [>= 0, default = 99] (ignored)";
	  "ew", Arg.Int (fun n -> ew_order := EW_order n), 
	  " set QCD order n [>=0, default = 99] (ignored)"]

    exception Negative_QCD_order
    exception Negative_EW_order
    exception Vanishing_couplings      
    exception Negative_QCD_EW_orders

    let int_orders = 
      match !qcd_order, !ew_order with
	| QCD_order n, EW_order n' when n < 0 &&  n' >= 0 -> 
	    raise Negative_QCD_order
	| QCD_order n, EW_order n' when n >= 0 &&  n' < 0 -> 
	    raise Negative_EW_order
	| QCD_order n, EW_order n' when n < 0 && n' < 0 -> 
	    raise Negative_QCD_EW_orders
	| QCD_order n, EW_order n' -> (n, n')

    open Coupling

    module S = Stat(M)

    type stat = S.stat
    let stat = S.stat
    let stat_sign = S.stat_sign

(* \begin{dubious}
     This will do \emph{something} for 4-, 6-, \ldots fermion vertices,
     but not necessarily the right thing \ldots
   \end{dubious} *)

    (* \begin{dubious}
         This is copied from [Colorize] and should be factored!
       \end{dubious} *)

    (* \begin{dubious}
         In the long run, it will probably be beneficial to apply
         the permutations in [Modeltools.add_vertexn]!
       \end{dubious} *)

    module PosMap =
      Partial.Make (struct type t = int let compare = compare end)

    let partial_map_undoing_permutation l l' =
      let module P = Permutation.Default in
      let p = P.of_list (List.map pred l') in
      PosMap.of_lists l (P.list p l)

    let partial_map_undoing_fuse fuse =
      partial_map_undoing_permutation
        (ThoList.range 1 (List.length fuse))
        fuse

    let undo_permutation_of_fuse fuse =
      PosMap.apply_with_fallback
        (fun _ -> invalid_arg "permutation_of_fuse")
        (partial_map_undoing_fuse fuse)

    let fermion_lines = function
      | Coupling.V3 _ | Coupling.V4 _ -> None
      | Coupling.Vn (Coupling.UFO (_, _, _, fl, _), fuse, _) ->
         Some (UFO_Lorentz.map_fermion_lines (undo_permutation_of_fuse fuse) fl)

    type constant = M.constant

(* \thocwmodulesubsection{Wave Functions} *)

(* \begin{dubious}
     The code below is not yet functional.  Too often, we assign to
     [Tags.null_wf] instead of calling [Tags.fuse].
   \end{dubious} *)

(* We will need two types of amplitudes: with color and without color.  Since
   we can build them using the same types with only [flavor] replaced, it pays
   to use a functor to set up the scaffolding. *)

    module Tags = Tagger(PT)

(* In the future, we might want to have [Coupling] among the functor
   arguments.  However, for the moment, [Coupling] is assumed to be
   comprehensive. *)

    module type Tagged_Coupling =
      sig
        type sign = int
        type t =
            { sign : sign;
              coupling : constant Coupling.t;
              coupling_tag : Tags.coupling }
        val sign : t -> sign
        val coupling : t -> constant Coupling.t
        val coupling_tag : t -> string option
      end

    module Tagged_Coupling : Tagged_Coupling =
      struct
        type sign = int
        type t =
            { sign : sign;
              coupling : constant Coupling.t;
              coupling_tag : Tags.coupling }
        let sign c = c.sign
        let coupling c = c.coupling
        let coupling_tag_raw c = c.coupling_tag
        let coupling_tag rhs = Tags.coupling_to_string (coupling_tag_raw rhs)
      end

(* \thocwmodulesubsection{Amplitudes: Monochrome and Colored} *)

    module type Amplitude =
      sig

        module Tags : Tags

        type flavor
        type p

        type wf =
            { flavor : flavor;
              momentum : p;
              wf_tag : Tags.wf }

        val flavor : wf -> flavor
        val conjugate : wf -> wf
        val momentum : wf -> p
        val momentum_list : wf -> int list
        val wf_tag : wf -> string option
	val wf_tag_raw : wf -> Tags.wf
        val order_wf : wf -> wf -> int
        val external_wfs : int -> (flavor * int) list -> wf list

        type 'a children
        type coupling = Tagged_Coupling.t
        type rhs = coupling * wf children
        val sign : rhs -> int
        val coupling : rhs -> constant Coupling.t
        val coupling_tag : rhs -> string option
	type exclusions
	val no_exclusions : exclusions
	    
        val children : rhs -> wf list

        type fusion = wf * rhs list
        val lhs : fusion -> wf
        val rhs : fusion -> rhs list

        type braket = wf * rhs list
        val bra : braket -> wf
        val ket : braket -> rhs list

        module D :
            DAG.T with type node = wf and type edge = coupling and type children = wf children

        val wavefunctions : braket list -> wf list

        type amplitude =
            { fusions : fusion list;
              brakets : braket list;
              on_shell : (wf -> bool);
              is_gauss : (wf -> bool);
              constraints : string option;
              incoming : flavor list;
              outgoing : flavor list;
              externals : wf list;
              symmetry : int;
              dependencies : (wf -> (wf, coupling) Tree2.t);
              fusion_tower : D.t;
              fusion_dag : D.t }

        val incoming : amplitude -> flavor list
        val outgoing : amplitude -> flavor list
        val externals : amplitude -> wf list
        val variables : amplitude -> wf list
        val fusions : amplitude -> fusion list
        val brakets : amplitude -> braket list
        val on_shell : amplitude -> (wf -> bool)
        val is_gauss : amplitude -> (wf -> bool)
        val constraints : amplitude -> string option
        val symmetry : amplitude -> int
        val dependencies : amplitude -> wf -> (wf, coupling) Tree2.t
        val fusion_dag : amplitude -> D.t

      end

    module Amplitude (PT : Tuple.Poly) (P : Momentum.T) (M : Model.T) :
        Amplitude
        with type p = P.t
        and type flavor = M.flavor
        and type 'a children = 'a PT.t
        and module Tags = Tags =
      struct

        type flavor = M.flavor
        type p = P.t

        module Tags = Tags

        type wf =
            { flavor : flavor;
              momentum : p;
              wf_tag : Tags.wf }

        let flavor wf = wf.flavor
        let conjugate wf = { wf with flavor = M.conjugate wf.flavor }
        let momentum wf = wf.momentum
        let momentum_list wf = P.to_ints wf.momentum
        let wf_tag wf = Tags.wf_to_string wf.wf_tag
        let wf_tag_raw wf = wf.wf_tag

        let external_wfs rank particles =
          List.map
            (fun (f, p) ->
              { flavor = f;
                momentum = P.singleton rank p;
                wf_tag = Tags.null_wf })
            particles

(* Order wavefunctions so that the external come first, then the pairs, etc.
   Also put possible Goldstone bosons \emph{before} their gauge bosons. *)

        let lorentz_ordering f =
          match M.lorentz f with
          | Coupling.Scalar -> 0
          | Coupling.Spinor -> 1
          | Coupling.ConjSpinor -> 2
          | Coupling.Majorana -> 3
          | Coupling.Vector -> 4
          | Coupling.Massive_Vector -> 5
          | Coupling.Tensor_2 -> 6
          | Coupling.Tensor_1 -> 7
          | Coupling.Vectorspinor -> 8
          | Coupling.BRS Coupling.Scalar -> 9
          | Coupling.BRS Coupling.Spinor -> 10
          | Coupling.BRS Coupling.ConjSpinor -> 11
          | Coupling.BRS Coupling.Majorana -> 12
          | Coupling.BRS Coupling.Vector -> 13
          | Coupling.BRS Coupling.Massive_Vector -> 14
          | Coupling.BRS Coupling.Tensor_2 -> 15
          | Coupling.BRS Coupling.Tensor_1 -> 16
          | Coupling.BRS Coupling.Vectorspinor -> 17
          | Coupling.BRS _ -> invalid_arg "Fusion.lorentz_ordering: not needed"
          | Coupling.Maj_Ghost -> 18
    (*i   | Coupling.Ward_Vector -> 19  i*)

        let order_flavor f1 f2 =
          let c = compare (lorentz_ordering f1) (lorentz_ordering f2) in
          if c <> 0 then
            c
          else
            compare f1 f2

(* Note that [Momentum().compare] guarantees that wavefunctions will be
   ordered according to \emph{increasing} [Momentum().rank] of their
   momenta. *)

        let order_wf wf1 wf2 =
          let c = P.compare wf1.momentum wf2.momentum in
          if c <> 0 then
            c
          else
            let c = order_flavor wf1.flavor wf2.flavor in
            if c <> 0 then
              c
            else
              compare wf1.wf_tag wf2.wf_tag

(* This \emph{must} be a pair matching the [edge * node children] pairs of
   [DAG.Forest]! *)

        type coupling = Tagged_Coupling.t
        type 'a children = 'a PT.t
        type rhs = coupling * wf children
        let sign (c, _) = Tagged_Coupling.sign c
        let coupling (c, _) = Tagged_Coupling.coupling c
        let coupling_tag (c, _) = Tagged_Coupling.coupling_tag c
	type exclusions =
	  { x_flavors : flavor list;
	    x_couplings : coupling list }
	let no_exclusions = { x_flavors = []; x_couplings = [] }
        let children (_, wfs) = PT.to_list wfs

        type fusion = wf * rhs list
        let lhs (l, _) = l
        let rhs (_, r) = r

        type braket = wf * rhs list
        let bra (b, _) = b
        let ket (_, k) = k

        module D = DAG.Make
            (DAG.Forest(PT)
               (struct type t = wf let compare = order_wf end)
               (struct type t = coupling let compare = compare end))

        module WFSet =
          Set.Make (struct type t = wf let compare = order_wf end)

        let wavefunctions brakets =
          WFSet.elements (List.fold_left (fun set (wf1, wf23) ->
            WFSet.add wf1 (List.fold_left (fun set' (_, wfs) ->
              PT.fold_right WFSet.add wfs set') set wf23)) WFSet.empty brakets)
          
        type amplitude =
            { fusions : fusion list;
              brakets : braket list;
              on_shell : (wf -> bool);
              is_gauss : (wf -> bool);
              constraints : string option;
              incoming : flavor list;
              outgoing : flavor list;
              externals : wf list;
              symmetry : int;
              dependencies : (wf -> (wf, coupling) Tree2.t);
              fusion_tower : D.t;
              fusion_dag : D.t }

        let incoming a = a.incoming
        let outgoing a = a.outgoing
        let externals a = a.externals
        let fusions a = a.fusions
        let brakets a = a.brakets
        let symmetry a = a.symmetry
        let on_shell a = a.on_shell
        let is_gauss a = a.is_gauss
        let constraints a = a.constraints
        let variables a = List.map lhs a.fusions
        let dependencies a = a.dependencies
        let fusion_dag a = a.fusion_dag

      end

    module A = Amplitude(PT)(P)(M)

(* Operator insertions can be fused only if they are external. *)
    let is_source wf =
      match M.propagator wf.A.flavor with
      | Only_Insertion -> P.rank wf.A.momentum = 1
      | _ -> true

(* [is_goldstone_of g v] is [true] if and only if [g] is the Goldstone boson
   corresponding to the gauge particle [v]. *)
    let is_goldstone_of g v =
      match M.goldstone v with
      | None -> false
      | Some (g', _) -> g = g'

(* \begin{dubious}
     In the end, [PT.to_list] should become redudant!
   \end{dubious} *)
    let fuse_rhs rhs = M.fuse (PT.to_list rhs)

(* \thocwmodulesubsection{Vertices} *)

(* Compute the set of all vertices in the model from the allowed
   fusions and the set of all flavors:
   \begin{dubious}
     One could think of using [M.vertices] instead of [M.fuse2],
     [M.fuse3] and [M.fuse] \ldots
   \end{dubious} *)

    module VSet = Map.Make(struct type t = A.flavor let compare = compare end)

    let add_vertices f rhs m =
      VSet.add f (try rhs :: VSet.find f m with Not_found -> [rhs]) m

    let collect_vertices rhs =
      List.fold_right (fun (f1, c) -> add_vertices (M.conjugate f1) (c, rhs))
        (fuse_rhs rhs)

(* The set of all vertices with common left fields factored. *)

(*   I used to think that constant initializers are a good idea to allow
     compile time optimizations.  The down side turned out to be that the
     constant initializers will be evaluated \emph{every time} the functor
     is applied.   \emph{Relying on the fact that the functor will be
     called only once is not a good idea!} *)

    type vertices = (A.flavor * (constant Coupling.t * A.flavor PT.t) list) list

(* \begin{dubious}
     This is \emph{very} inefficient for [max_degree > 6].  Find a better
     approach that avoids precomputing the huge lookup table!
   \end{dubious}
   \begin{dubious}
     I should revive the above Idea to use [M.vertices] instead directly,
     instead of rebuilding it from [M.fuse2],
     [M.fuse3] and [M.fuse]!
   \end{dubious} *)

    let vertices_nocache max_degree flavors : vertices =
      VSet.fold (fun f rhs v -> (f, rhs) :: v)
        (PT.power_fold
           ~truncate:(pred max_degree)
           collect_vertices flavors VSet.empty) []

(* Performance hack: *)

    type vertex_table =
            ((A.flavor * A.flavor * A.flavor) * constant Coupling.vertex3 * constant) list
          * ((A.flavor * A.flavor * A.flavor * A.flavor)
               * constant Coupling.vertex4 * constant) list
          * (A.flavor list * constant Coupling.vertexn * constant) list

(*i
    module VCache =
      Cache.Make (struct type t = vertex_table end) (struct type t = vertices end)

    let vertices_cache = ref None
    let hash () = VCache.hash (M.vertices ())

(* \begin{dubious}
     Can we do better than the executable name provided by [Config.cache_prefix]???
     We need a better way to avoid collisions among the caches for different models
     in the same program.
   \end{dubious} *)

    let cache_name =
      ref (Config.cache_prefix ^ "." ^ Config.cache_suffix)

    let set_cache_name name = 
      cache_name := name

    let initialize_cache dir =
      Printf.eprintf
        " >>> Initializing vertex table %s.  This may take some time ... "
        !cache_name;
      flush stderr;
      VCache.write_dir (hash ()) dir !cache_name
        (vertices_nocache  (M.max_degree ()) (M.flavors()));
      Printf.eprintf "done. <<< \n"

    let vertices max_degree flavors : vertices =
      match !vertices_cache with 
      | None -> 
          begin match !cache_option with
          | Cache_Use ->
              begin match VCache.maybe_read (hash ()) !cache_name with
              | VCache.Hit result -> result
              | VCache.Miss ->
                  Printf.eprintf
                    " >>> Initializing vertex table %s.  This may take some time ... "
                    !cache_name;
                  flush stderr;
                  let result = vertices_nocache max_degree flavors in
                  VCache.write (hash ()) !cache_name (result);
                  vertices_cache := Some result;
                  Printf.eprintf "done. <<< \n";
                  flush stderr;
                  result
              | VCache.Stale file ->
                  Printf.eprintf
                    " >>> Re-initializing stale vertex table %s in file %s.  "
                    !cache_name file;
                  Printf.eprintf "This may take some time ... ";
                  flush stderr;
                  let result = vertices_nocache max_degree flavors in
                  VCache.write (hash ()) !cache_name (result);
                  vertices_cache := Some result;
                  Printf.eprintf "done. <<< \n";
                  flush stderr;
                  result
              end
          | Cache_Overwrite ->
              Printf.eprintf
                " >>> Overwriting vertex table %s.  This may take some time ... "
                !cache_name;
              flush stderr;
              let result = vertices_nocache max_degree flavors in
              VCache.write (hash ()) !cache_name (result);
              vertices_cache := Some result;
              Printf.eprintf "done. <<< \n";
              flush stderr;
              result
          | Cache_Ignore ->
              let result = vertices_nocache max_degree flavors in
              vertices_cache := Some result;
              result
          end
      | Some result -> result
i*)
    let vertices = vertices_nocache

    let vertices' max_degree flavors =
      Printf.eprintf ">>> vertices %d ..." max_degree;
      flush stderr;
      let v = vertices max_degree flavors in
      Printf.eprintf " done.\n";
      flush stderr;
      v

(* Note that we must perform any filtering of the vertices \emph{after}
   caching, because the restrictions \emph{must not} influence the
   cache (unless we tag the cache with model and restrictions).  *)

(*i
    let unpack_constant = function
      | Coupling.V3 (_, _, cs) -> cs
      | Coupling.V4 (_, _, cs) -> cs
      | Coupling.Vn (_, _, cs) -> cs

    let coupling_and_flavors_to_string (c, fs) =
      M.constant_symbol (unpack_constant c) ^ "[" ^
	String.concat ", " (List.map M.flavor_to_string (PT.to_list fs)) ^ "]"

    let fusions_to_string (f, cfs) =
      M.flavor_to_string f ^ " <- { " ^
	String.concat " | " (List.map coupling_and_flavors_to_string cfs) ^
	" }"

    let vertices_to_string vertices =
      String.concat "; " (List.map fusions_to_string vertices)
  i*)

    let filter_vertices select_vtx vertices =
      List.fold_left
	(fun acc (f, cfs) ->
	  let f' = M.conjugate f in
	  let cfs =
	    List.filter
	      (fun (c, fs) -> select_vtx c f' (PT.to_list fs))
	      cfs
	  in
	  match cfs with
	  | [] -> acc
	  | cfs -> (f, cfs) :: acc)
	[] vertices

(* \thocwmodulesubsection{Partitions} *)

(* Vertices that are not crossing invariant need special treatment so
   that they're only generated for the correct combinations of momenta.

   NB: the [crossing] checks here are a bit redundant, because  [CM.fuse] below
   will bring the killed vertices back to life and will have to filter once more.
   Nevertheless, we keep them here, for the unlikely case that anybody ever wants
   to use uncolored amplitudes directly.

   NB: the analogous problem does not occur for [select_wf], because this applies
   to momenta instead of vertices. *)

(* \begin{dubious}
     This approach worked before the colorize, but has become \emph{futile},
     because [CM.fuse] will bring the killed vertices back to life.  We need
     to implement the same checks there again!!!
   \end{dubious}  *)

(* \begin{dubious}
     Using [PT.Mismatched_arity] is not really good style \ldots

   Tho's approach doesn't work since he does not catch charge conjugated processes or
   crossed processes. Another very strange thing is that O'Mega seems always to run in the
   q2 q3 timelike case, but not in the other two. (Property of how the DAG is built?).    
   For the $ZZZZ$ vertex I add the same vertex again, but interchange 1 and 3 in the 
   [crossing] vertex

   \end{dubious} *)

    let kmatrix_cuts c momenta =
      match c with
      | V4 (Vector4_K_Matrix_tho (disc, _), fusion, _) 
      | V4 (Vector4_K_Matrix_jr (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end
      | V4 (Vector4_K_Matrix_cf_t0 (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end          
      | V4 (Vector4_K_Matrix_cf_t1 (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end          
      | V4 (Vector4_K_Matrix_cf_t2 (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end
      | V4 (Vector4_K_Matrix_cf_t_rsi (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end
      | V4 (Vector4_K_Matrix_cf_m0 (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end
      | V4 (Vector4_K_Matrix_cf_m1 (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end
      | V4 (Vector4_K_Matrix_cf_m7 (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F431|F342|F432) 
          | 1, false, true, false, (F134|F143|F234|F243)
          | 1, false, false, true, (F314|F413|F324|F423) ->
              true
          | 2, true, false, false, (F123|F213|F124|F214)
          | 2, false, true, false, (F312|F321|F412|F421)
          | 2, false, false, true, (F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true 
          | _ -> false 
          end    
      | V4 (DScalar2_Vector2_K_Matrix_ms (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F432|F123|F214) 
          | 1, false, true, false, (F134|F243|F312|F421)
          | 1, false, false, true, (F314|F423|F132|F241) ->
              true
          | 2, true, false, false, (F431|F342|F213|F124)
          | 2, false, true, false, (F143|F234|F321|F412)
          | 2, false, false, true, (F413|F324|F231|F142) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true
          | 4, true, false, false, (F142|F413|F231|F324)
          | 4, false, true, false, (F214|F341|F123|F432)
          | 4, false, false, true, (F124|F431|F213|F342) ->
              true
          | 5, true, false, false, (F143|F412|F321|F234)
          | 5, false, true, false, (F314|F241|F132|F423)
          | 5, false, false, true, (F134|F421|F312|F243) ->
              true
          | 6, true, false, false, (F134|F132|F314|F312|F241|F243|F421|F423)
          | 6, false, true, false, (F213|F413|F231|F431|F124|F324|F142|F342)
          | 6, false, false, true, (F143|F123|F341|F321|F412|F214|F432|F234) ->
              true
          | 7, true, false, false, (F134|F312|F421|F243)
          | 7, false, true, false, (F413|F231|F142|F324)
          | 7, false, false, true, (F143|F321|F412|F432) ->
              true
          | 8, true, false, false, (F132|F314|F241|F423)
          | 8, false, true, false, (F213|F431|F124|F342)
          | 8, false, false, true, (F123|F341|F214|F234) ->
              true
          | _ -> false
          end
      | V4 (DScalar2_Vector2_m_0_K_Matrix_cf (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F432|F123|F214) 
          | 1, false, true, false, (F134|F243|F312|F421)
          | 1, false, false, true, (F314|F423|F132|F241) ->
              true
          | 2, true, false, false, (F431|F342|F213|F124)
          | 2, false, true, false, (F143|F234|F321|F412)
          | 2, false, false, true, (F413|F324|F231|F142) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true
          | 4, true, false, false, (F142|F413|F231|F324)
          | 4, false, true, false, (F214|F341|F123|F432)
          | 4, false, false, true, (F124|F431|F213|F342) ->
              true
          | 5, true, false, false, (F143|F412|F321|F234)
          | 5, false, true, false, (F314|F241|F132|F423)
          | 5, false, false, true, (F134|F421|F312|F243) ->
              true
          | 6, true, false, false, (F134|F132|F314|F312|F241|F243|F421|F423)
          | 6, false, true, false, (F213|F413|F231|F431|F124|F324|F142|F342)
          | 6, false, false, true, (F143|F123|F341|F321|F412|F214|F432|F234) ->
              true
          | 7, true, false, false, (F134|F312|F421|F243)
          | 7, false, true, false, (F413|F231|F142|F324)
          | 7, false, false, true, (F143|F321|F412|F432) ->
              true
          | 8, true, false, false, (F132|F314|F241|F423)
          | 8, false, true, false, (F213|F431|F124|F342)
          | 8, false, false, true, (F123|F341|F214|F234) ->
              true
          | _ -> false
          end 
      | V4 (DScalar2_Vector2_m_1_K_Matrix_cf (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F432|F123|F214) 
          | 1, false, true, false, (F134|F243|F312|F421)
          | 1, false, false, true, (F314|F423|F132|F241) ->
              true
          | 2, true, false, false, (F431|F342|F213|F124)
          | 2, false, true, false, (F143|F234|F321|F412)
          | 2, false, false, true, (F413|F324|F231|F142) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true
          | 4, true, false, false, (F142|F413|F231|F324)
          | 4, false, true, false, (F214|F341|F123|F432)
          | 4, false, false, true, (F124|F431|F213|F342) ->
              true
          | 5, true, false, false, (F143|F412|F321|F234)
          | 5, false, true, false, (F314|F241|F132|F423)
          | 5, false, false, true, (F134|F421|F312|F243) ->
              true
          | 6, true, false, false, (F134|F132|F314|F312|F241|F243|F421|F423)
          | 6, false, true, false, (F213|F413|F231|F431|F124|F324|F142|F342)
          | 6, false, false, true, (F143|F123|F341|F321|F412|F214|F432|F234) ->
              true
          | 7, true, false, false, (F134|F312|F421|F243)
          | 7, false, true, false, (F413|F231|F142|F324)
          | 7, false, false, true, (F143|F321|F412|F432) ->
              true
          | 8, true, false, false, (F132|F314|F241|F423)
          | 8, false, true, false, (F213|F431|F124|F342)
          | 8, false, false, true, (F123|F341|F214|F234) ->
              true
          | _ -> false
          end  
      | V4 (DScalar2_Vector2_m_7_K_Matrix_cf (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 1, true, false, false, (F341|F432|F123|F214) 
          | 1, false, true, false, (F134|F243|F312|F421)
          | 1, false, false, true, (F314|F423|F132|F241) ->
              true
          | 2, true, false, false, (F431|F342|F213|F124)
          | 2, false, true, false, (F143|F234|F321|F412)
          | 2, false, false, true, (F413|F324|F231|F142) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true
          | 4, true, false, false, (F142|F413|F231|F324)
          | 4, false, true, false, (F214|F341|F123|F432)
          | 4, false, false, true, (F124|F431|F213|F342) ->
              true
          | 5, true, false, false, (F143|F412|F321|F234)
          | 5, false, true, false, (F314|F241|F132|F423)
          | 5, false, false, true, (F134|F421|F312|F243) ->
              true
          | 6, true, false, false, (F134|F132|F314|F312|F241|F243|F421|F423)
          | 6, false, true, false, (F213|F413|F231|F431|F124|F324|F142|F342)
          | 6, false, false, true, (F143|F123|F341|F321|F412|F214|F432|F234) ->
              true
          | 7, true, false, false, (F134|F312|F421|F243)
          | 7, false, true, false, (F413|F231|F142|F324)
          | 7, false, false, true, (F143|F321|F412|F432) ->
              true
          | 8, true, false, false, (F132|F314|F241|F423)
          | 8, false, true, false, (F213|F431|F124|F342)
          | 8, false, false, true, (F123|F341|F214|F234) ->
              true
          | _ -> false
          end    
      | V4 (DScalar4_K_Matrix_ms (disc, _), fusion, _) ->
          let s12, s23, s13 =
            begin match PT.to_list momenta with
            | [q1; q2; q3] -> (P.Scattering.timelike (P.add q1 q2),
                               P.Scattering.timelike (P.add q2 q3),
                               P.Scattering.timelike (P.add q1 q3))
            | _ -> raise PT.Mismatched_arity
            end in
          begin match disc, s12, s23, s13, fusion with
          | 0, true, false, false, (F341|F431|F342|F432|F123|F213|F124|F214)
          | 0, false, true, false, (F134|F143|F234|F243|F312|F321|F412|F421)
          | 0, false, false, true, (F314|F413|F324|F423|F132|F231|F142|F241) ->
              true
          | 3, true, false, false, (F143|F413|F142|F412|F321|F231|F324|F234)
          | 3, false, true, false, (F314|F341|F214|F241|F132|F123|F432|F423)
          | 3, false, false, true, (F134|F431|F124|F421|F312|F213|F342|F243) ->
              true
          | 4, true, false, false, (F142|F413|F231|F324)
          | 4, false, true, false, (F214|F341|F123|F432)
          | 4, false, false, true, (F124|F431|F213|F342) ->
              true
          | 5, true, false, false, (F143|F412|F321|F234)
          | 5, false, true, false, (F314|F241|F132|F423)
          | 5, false, false, true, (F134|F421|F312|F243) ->
              true
          | 6, true, false, false, (F134|F132|F314|F312|F241|F243|F421|F423)
          | 6, false, true, false, (F213|F413|F231|F431|F124|F324|F142|F342)
          | 6, false, false, true, (F143|F123|F341|F321|F412|F214|F432|F234) ->
              true
          | 7, true, false, false, (F134|F312|F421|F243)
          | 7, false, true, false, (F413|F231|F142|F324)
          | 7, false, false, true, (F143|F321|F412|F432) ->
              true
          | 8, true, false, false, (F132|F314|F241|F423)
          | 8, false, true, false, (F213|F431|F124|F342)
          | 8, false, false, true, (F123|F341|F214|F234) ->
              true
          | _ -> false
          end
      | _ -> true


(* Counting QCD and EW orders. *)

    let qcd_ew_check orders = 
      if fst (orders) <= fst (int_orders) &&
	 snd (orders) <= snd (int_orders) then
	true
      else
	false


(* Match a set of flavors to a set of momenta.  Form the direct product for
   the lists of momenta two and three with the list of couplings and flavors
   two and three.  *)

    let flavor_keystone select_p dim (f1, f23) (p1, p23) =
      ({ A.flavor = f1;
         A.momentum = P.of_ints dim p1;
         A.wf_tag = A.Tags.null_wf },
       Product.fold2 (fun (c, f) p acc ->
         try
           let p' = PT.map (P.of_ints dim) p in
           if select_p (P.of_ints dim p1) (PT.to_list p') && kmatrix_cuts c p' then
             (c, PT.map2 (fun f'' p'' -> { A.flavor = f'';
                                           A.momentum = p'';
                                           A.wf_tag = A.Tags.null_wf }) f p') :: acc
           else
             acc
         with
         | PT.Mismatched_arity -> acc) f23 p23 [])

(*i
    let cnt = ref 0

    let gc_stat () =
      let minor, promoted, major = Gc.counters () in
      Printf.sprintf "(%12.0f, %12.0f, %12.0f)" minor promoted major

    let flavor_keystone select_p n (f1, f23) (p1, p23) =
      incr cnt;
      Gc.set { (Gc.get()) with Gc.space_overhead = 20 };
      Printf.eprintf "%6d@%8.1f: %s\n" !cnt (Sys.time ()) (gc_stat ());
      flush stderr;
      flavor_keystone select_p n (f1, f23) (p1, p23)
i*)

(* Produce all possible combinations of vertices (flavor keystones)
   and momenta by forming the direct product.  The semantically equivalent
   [Product.list2 (flavor_keystone select_wf n) vertices keystones] with
   \emph{subsequent} filtering would be a \emph{very bad} idea, because
   a potentially huge intermediate list is built for large models.
   E.\,g.~for the MSSM this would lead to non-termination by thrashing
   for $2\to4$ processes on most PCs. *)

    let flavor_keystones filter select_p dim vertices keystones =
      Product.fold2 (fun v k acc ->
        filter (flavor_keystone select_p dim v k) acc) vertices keystones []

(* Flatten the nested lists of vertices into a list of attached lines. *)

    let flatten_keystones t =
      ThoList.flatmap (fun (p1, p23) ->
        p1 :: (ThoList.flatmap (fun (_, rhs) -> PT.to_list rhs) p23)) t

(* \thocwmodulesubsection{Subtrees} *)

(* Fuse a tuple of wavefunctions, keeping track of Fermi statistics.
   Record only the the sign \emph{relative} to the children.
   (The type annotation is only for documentation.) *)

    let fuse select_wf select_vtx wfss : (A.wf * stat * A.rhs) list =
      if PT.for_all (fun (wf, _) -> is_source wf) wfss then
        try
          let wfs, ss = PT.split wfss in
          let flavors = PT.map A.flavor wfs
          and momenta = PT.map A.momentum wfs
(*i       and wf_tags = PT.map A.wf_tag_raw wfs i*) in
          let p = PT.fold_left_internal P.add momenta in
(*i	  let wft = PT.fold_left Tags.fuse wf_tags in i*)
          List.fold_left
            (fun acc (f, c) ->
              if select_wf f p (PT.to_list momenta)
		&& select_vtx c f (PT.to_list flavors)
		&& kmatrix_cuts c momenta then
                (* [let _ = 
                  Printf.eprintf
                    "Fusion.fuse: %s <- %s\n"
                    (M.flavor_to_string f)
                    (ThoList.to_string M.flavor_to_string (PT.to_list flavors)) in] *)
                let s = S.stat_fuse (fermion_lines c) (PT.to_list ss) f in
                let flip =
                  PT.fold_left (fun acc s' -> acc * stat_sign s') (stat_sign s) ss in
                ({ A.flavor = f;
                   A.momentum = p;
                   A.wf_tag = A.Tags.null_wf }, s,
                 ({ Tagged_Coupling.sign = flip;
                    Tagged_Coupling.coupling = c;
                    Tagged_Coupling.coupling_tag = A.Tags.null_coupling }, wfs)) :: acc
              else
                acc)
            [] (fuse_rhs flavors)
        with
        | P.Duplicate _ | S.Impossible -> []
      else
        []

(* \begin{dubious}
     Eventually, the pairs of [tower] and [dag] in [fusion_tower']
     below could and should be replaced by a graded [DAG].  This will
     look like, but currently [tower] containts statistics information
     that is missing from [dag]:
     \begin{quote}
       \verb+Type node = flavor * p is not compatible with type wf * stat+
     \end{quote}
     This should be easy to fix.  However, replacing [type t = wf]
     with [type t = wf * stat] is \emph{not} a good idea because the variable
     [stat] makes it impossible to test for the existance of a particular
     [wf] in a [DAG].
   \end{dubious}
   \begin{dubious}
     In summary, it seems that [(wf * stat) list array * A.D.t] should be
     replaced by [(wf -> stat) * A.D.t].
   \end{dubious} *)
    module GF =
      struct
        module Nodes =
          struct
            type t = A.wf
            module G = struct type t = int let compare = compare end
            let compare = A.order_wf
            let rank wf = P.rank wf.A.momentum
          end
        module Edges = struct type t = A.coupling let compare = compare end
        module F = DAG.Forest(PT)(Nodes)(Edges)
        type node = Nodes.t
        type edge = F.edge
        type children = F.children
        type t = F.t
        let compare = F.compare
        let for_all = F.for_all
        let fold = F.fold
      end

    module D' = DAG.Graded(GF)

    let tower_of_dag dag =
      let _, max_rank = D'.min_max_rank dag in
      Array.init max_rank (fun n -> D'.ranked n dag)

(* The function [fusion_tower']
   recursively builds the tower of all fusions from bottom up to a chosen
   level. The argument [tower] is an array of lists, where the $i$-th sublist
   (counting from 0) represents all off shell wave functions depending on
   $i+1$~momenta and their Fermistatistics.
   \begin{equation}
     \begin{aligned}
         \Bigl\lbrack
                 & \{ \phi_1(p_1), \phi_2(p_2), \phi_3(p_3), \ldots \}, \\
                 & \{ \phi_{12}(p_1+p_2), \phi'_{12}(p_1+p_2), \ldots,
                      \phi_{13}(p_1+p_3), \ldots, \phi_{23}(p_2+p_3), \ldots \}, \\
                 & \ldots \\
                 & \{ \phi_{1\cdots n}(p_1+\cdots+p_n),
                      \phi'_{1\cdots n}(p_1+\cdots+p_n), \ldots \} \Bigr\rbrack
     \end{aligned}
   \end{equation}
   The argument [dag] is a DAG representing all the fusions calculated so far.
   NB: The outer array in [tower] is always very short, so we could also
   have accessed a list with [List.nth].   Appending of new members at the
   end brings no loss of performance.  NB: the array is supposed to be
   immutable.  *)

(* The towers must be sorted so that the combinatorical functions can
   make consistent selections.
   \begin{dubious}
     Intuitively, this seems to be correct.  However, one could have
     expected that no element appears twice and that this ordering is
     not necessary \ldots
   \end{dubious} *)
    let grow select_wf select_vtx tower =
      let rank = succ (Array.length tower) in
      List.sort pcompare
        (PT.graded_sym_power_fold rank
           (fun wfs acc -> fuse select_wf select_vtx wfs @ acc) tower [])

    let add_offspring dag (wf, _, rhs) =
      A.D.add_offspring wf rhs dag

    let filter_offspring fusions =
      List.map (fun (wf, s, _) -> (wf, s)) fusions

    let rec fusion_tower' n_max select_wf select_vtx tower dag : (A.wf * stat) list array * A.D.t =
      if Array.length tower >= n_max then
        (tower, dag)
      else
        let tower' = grow select_wf select_vtx tower in
        fusion_tower' n_max select_wf select_vtx
          (Array.append tower [|filter_offspring tower'|])
          (List.fold_left add_offspring dag tower')

(* Discard the tower and return a map from wave functions to Fermistatistics
   together with the DAG. *)

    let make_external_dag wfs =
      List.fold_left (fun m (wf, _) -> A.D.add_node wf m) A.D.empty wfs

    let mixed_fold_left f acc lists =
      Array.fold_left (List.fold_left f) acc lists

    module Stat_Map =
      Map.Make (struct type t = A.wf let compare = A.order_wf end)

    let fusion_tower height select_wf select_vtx wfs : (A.wf -> stat) * A.D.t =
      let tower, dag =
        fusion_tower' height select_wf select_vtx [|wfs|] (make_external_dag wfs) in
      let stats = mixed_fold_left
          (fun m (wf, s) -> Stat_Map.add wf s m) Stat_Map.empty tower in
      ((fun wf -> Stat_Map.find wf stats), dag)

(* Calculate the minimal tower of fusions that suffices for calculating
   the amplitude.  *)

    let minimal_fusion_tower n select_wf select_vtx wfs : (A.wf -> stat) * A.D.t =
      fusion_tower (T.max_subtree n) select_wf select_vtx wfs

(* Calculate the complete tower of fusions.  It is much larger than required,
   but it allows a complete set of gauge checks.  *)
    let complete_fusion_tower select_wf select_vtx wfs : (A.wf -> stat) * A.D.t =
      fusion_tower (List.length wfs - 1) select_wf select_vtx wfs

(* \begin{dubious}
     There is a natural product of two DAGs using [fuse].  Can this be
     used in a replacement for [fusion_tower]?  The hard part is to avoid
     double counting, of course.  A straight forward solution
     could do a diagonal sum (in order to reject flipped offspring representing
     the same fusion) and rely on the uniqueness in [DAG] otherwise.
     However, this will (probably) slow down the procedure significanty,
     because most fusions (including Fermi signs!) will be calculated before
     being rejected by [DAG().add_offspring].
   \end{dubious} *)

(* Add to [dag] all Goldstone bosons defined in [tower] that correspond
   to gauge bosons in [dag].  This is only required for checking
   Slavnov-Taylor identities in unitarity gauge.  Currently, it is not used,
   because we use the complete tower for gauge checking. *)
    let harvest_goldstones tower dag =
      A.D.fold_nodes (fun wf dag' ->
        match M.goldstone wf.A.flavor with
        | Some (g, _) ->
            let wf' = { wf with A.flavor = g } in
            if A.D.is_node wf' tower then begin
              A.D.harvest tower wf' dag'
            end else begin
              dag'
            end
        | None -> dag') dag dag

(* Calculate the sign from Fermi statistics that is not already included
   in the children. *)

    let strip_fermion_lines = function
      | (Coupling.V3 _ | Coupling.V4 _ as v) -> v
      | Coupling.Vn (Coupling.UFO (c, l, s, fl, col), f, x) ->
         Coupling.Vn (Coupling.UFO (c, l, s, [], col), f, x)

    let num_fermion_lines_v3 = function
      | FBF _ | PBP _ | BBB _ | GBG _ -> 1
      | _ -> 0

    let num_fermion_lines = function
      | Coupling.Vn (Coupling.UFO (c, l, s, fl, col), f, x) -> List.length fl
      | Coupling.V3 (v3, _, _) -> num_fermion_lines_v3 v3
      | Coupling.V4 _ -> 0

    let stat_keystone v stats wf1 wfs =
      let wf1' = stats wf1
      and wfs' = PT.map stats wfs in
      let f = A.flavor wf1 in
      let slist = wf1' :: PT.to_list wfs' in
      let stat = S.stat_keystone (fermion_lines v) slist f in
      (* We can compare with the legacy implementation only if there
         are no fermion line ambiguities possible, i.\,e.~for
         at most one line. *)
      if num_fermion_lines v < 2 then
        begin
          let legacy = S.stat_keystone None slist f in
          if not (S.equal stat legacy) then
            failwith
              (Printf.sprintf
                 "Fusion.stat_keystone: %s <> %s!"
                 (S.stat_to_string legacy)
                 (S.stat_to_string stat));
          if not (S.saturated legacy) then
            failwith
              (Printf.sprintf
                 "Fusion.stat_keystone: legacy incomplete: %s!"
                 (S.stat_to_string legacy))
        end;
      if not (S.saturated stat) then
        failwith
          (Printf.sprintf
             "Fusion.stat_keystone: incomplete: %s!"
             (S.stat_to_string stat));
      stat_sign stat
        * PT.fold_left (fun acc wf -> acc * stat_sign wf) (stat_sign wf1') wfs'

    let stat_keystone_logging v stats wf1 wfs =
      let sign = stat_keystone v stats wf1 wfs in
      Printf.eprintf
        "Fusion.stat_keystone: %s * %s -> %d\n"
        (M.flavor_to_string (A.flavor wf1))
        (ThoList.to_string
           (fun wf -> M.flavor_to_string (A.flavor wf))
           (PT.to_list wfs))
        sign;
      sign

(* Test all members of a list of wave functions are defined by the DAG
   simultaneously: *)
    let test_rhs dag (_, wfs) =
      PT.for_all (fun wf -> is_source wf && A.D.is_node wf dag) wfs

(* Add the keystone [(wf1,pairs)] to [acc] only if it is present in [dag]
   and calculate the statistical factor depending on [stats]
   \emph{en passant}: *)
    let filter_keystone stats dag (wf1, pairs) acc =
      if is_source wf1 && A.D.is_node wf1 dag then
        match List.filter (test_rhs dag) pairs with
        | [] -> acc
        | pairs' -> (wf1, List.map (fun (c, wfs) ->
            ({ Tagged_Coupling.sign = stat_keystone c stats wf1 wfs;
               Tagged_Coupling.coupling = c;
               Tagged_Coupling.coupling_tag = A.Tags.null_coupling },
             wfs)) pairs') :: acc
      else
        acc

(* \begin{figure}
     \begin{center}
       \thocwincludegraphics{width=\textwidth}{bhabha0}\\
       \hfil\\
       \thocwincludegraphics{width=\textwidth}{bhabha}
     \end{center}
     \caption{\label{fig:bhabha}
       The DAGs for Bhabha scattering before and after weeding out unused
       nodes. The blatant asymmetry of these DAGs is caused by our
       prescription for removing doubling counting for an even number
       of external lines.}
   \end{figure}
   \begin{figure}
     \begin{center}
       \thocwincludegraphics{width=\textwidth}{epemudbarmunumubar0}\\
       \hfil\\
       \thocwincludegraphics{width=\textwidth}{epemudbarmunumubar}
     \end{center}
     \caption{\label{fig:epemudbarmunumubar}
       The DAGs for $e^+e^-\to u\bar d \mu^-\bar\nu_\mu$ before and after
       weeding out unused nodes.}
   \end{figure}
   \begin{figure}
     \begin{center}
       \thocwincludegraphics{width=\textwidth}{epemudbardubar0}\\
       \hfil\\
       \thocwincludegraphics{width=\textwidth}{epemudbardubar}
     \end{center}
     \caption{\label{fig:epemudbardubar}
       The DAGs for $e^+e^-\to u\bar d d\bar u$ before and after weeding
       out unused nodes.}
   \end{figure} *)

(* \thocwmodulesubsection{Amplitudes} *)

    module C = Cascade.Make(M)(P)
    type selectors = C.selectors

    let external_wfs n particles =
      List.map (fun (f, p) ->
        ({ A.flavor = f;
           A.momentum = P.singleton n p;
           A.wf_tag = A.Tags.null_wf },
         stat f p)) particles

(* \thocwmodulesubsection{Main Function} *)

    module WFMap = Map.Make (struct type t = A.wf let compare = compare end)

(* [map_amplitude_wfs f a] applies the function [f : wf -> wf] to all
   wavefunctions appearing in the amplitude [a]. *)
    let map_amplitude_wfs f a =
      let map_rhs (c, wfs) = (c, PT.map f wfs) in
      let map_braket (wf, rhs) = (f wf, List.map map_rhs rhs)
      and map_fusion (lhs, rhs) = (f lhs, List.map map_rhs rhs) in
      let map_dag = A.D.map f (fun node rhs -> map_rhs rhs) in
      let tower = map_dag a.A.fusion_tower
      and dag = map_dag a.A.fusion_dag in
      let dependencies_map =
        A.D.fold (fun wf _ -> WFMap.add wf (A.D.dependencies dag wf)) dag WFMap.empty in
      { A.fusions = List.map map_fusion a.A.fusions;
        A.brakets = List.map map_braket a.A.brakets;
        A.on_shell = a.A.on_shell;
        A.is_gauss = a.A.is_gauss;
        A.constraints = a.A.constraints;
        A.incoming = a.A.incoming;
        A.outgoing = a.A.outgoing;
        A.externals = List.map f a.A.externals;
        A.symmetry = a.A.symmetry;
        A.dependencies = (fun wf -> WFMap.find wf dependencies_map);
        A.fusion_tower = tower;
        A.fusion_dag = dag }

(*i
(* \begin{dubious}
     Just a silly little test:
   \end{dubious} *)

    let hack_amplitude =
      map_amplitude_wfs (fun wf -> { wf with momentum = P.split 2 16 wf.momentum })
i*)

(* This is the main function that constructs the amplitude for sets
   of incoming and outgoing particles and returns the results in
   conveniently packaged pieces.  *)

    let amplitude goldstones selectors fin fout =

      (* Set up external lines and match flavors with numbered momenta. *)
      let f = fin @ List.map M.conjugate fout in
      let nin, nout = List.length fin, List.length fout in
      let n = nin + nout in
      let externals = List.combine f (ThoList.range 1 n) in
      let wfs = external_wfs n externals in
      let select_p = C.select_p selectors in
      let select_wf =
        match fin with
        | [_] -> C.select_wf selectors P.Decay.timelike
        | _ -> C.select_wf selectors P.Scattering.timelike in
      let select_vtx = C.select_vtx selectors in

      (* Build the full fusion tower (including nodes that are never
         needed in the amplitude). *)
      let stats, tower =

        if goldstones then
          complete_fusion_tower select_wf select_vtx wfs
        else
          minimal_fusion_tower n select_wf select_vtx wfs in

      (* Find all vertices for which \emph{all} off shell wavefunctions
         are defined by the tower. *)

      let brakets =
        flavor_keystones (filter_keystone stats tower) select_p n
          (filter_vertices select_vtx
	     (vertices (min n (M.max_degree ())) (M.flavors ())))
          (T.keystones (ThoList.range 1 n)) in

      (* Remove the part of the DAG that is never needed in the amplitude. *)
      let dag =
        if goldstones then
          tower
        else
          A.D.harvest_list tower (A.wavefunctions brakets) in

      (* Remove the leaf nodes of the DAG, corresponding to external lines. *)
      let fusions =
        List.filter (function (_, []) -> false | _ -> true) (A.D.lists dag) in

      (* Calculate the symmetry factor for identical particles in the
         final state. *)
      let symmetry =
        Combinatorics.symmetry fout in

      let dependencies_map =
        A.D.fold (fun wf _ -> WFMap.add wf (A.D.dependencies dag wf)) dag WFMap.empty in
      
      (* Finally: package the results: *)
      { A.fusions = fusions;
        A.brakets = brakets;
        A.on_shell = (fun wf -> C.on_shell selectors (A.flavor wf) wf.A.momentum);
        A.is_gauss = (fun wf -> C.is_gauss selectors (A.flavor wf) wf.A.momentum);
        A.constraints = C.description selectors;
        A.incoming = fin;
        A.outgoing = fout;
        A.externals = List.map fst wfs;        
        A.symmetry = symmetry;
        A.dependencies = (fun wf -> WFMap.find wf dependencies_map);
        A.fusion_tower = tower;
        A.fusion_dag = dag }

(* \thocwmodulesubsection{Color} *)

    module CM = Colorize.It(M)
    module CA = Amplitude(PT)(P)(CM)

    let colorize_wf flavor wf =
      { CA.flavor = flavor;
        CA.momentum = wf.A.momentum;
        CA.wf_tag = wf.A.wf_tag }

    let uncolorize_wf wf =
      { A.flavor = CM.flavor_sans_color wf.CA.flavor;
        A.momentum = wf.CA.momentum;
        A.wf_tag = wf.CA.wf_tag }

(* \begin{dubious}
     At the end of the day, I shall want to have some sort of
     \textit{fibered DAG} as abstract data type, with a projection
     of colored nodes to their uncolored counterparts.
   \end{dubious} *)

    module CWFBundle = Bundle.Make
        (struct
          type elt = CA.wf
          let compare_elt = compare
          type base = A.wf
          let compare_base = compare
          let pi wf =
            { A.flavor = CM.flavor_sans_color wf.CA.flavor;
              A.momentum = wf.CA.momentum;
              A.wf_tag = wf.CA.wf_tag }
        end)

(* \begin{dubious}
     For now, we can live with simple aggregation:
   \end{dubious} *)

    type fibered_dag = { dag : CA.D.t; bundle : CWFBundle.t }

(* Not yet(?) needed: [module CS = Stat (CM)] *)

    let colorize_sterile_nodes dag f wf fibered_dag = 
      if A.D.is_sterile wf dag then
        let wf', wf_bundle' = f wf fibered_dag in
        { dag = CA.D.add_node wf' fibered_dag.dag;
          bundle = wf_bundle' }
      else
        fibered_dag

    let colorize_nodes f wf rhs fibered_dag =
      let wf_rhs_list', wf_bundle' = f wf rhs fibered_dag in
      let dag' =
        List.fold_right
          (fun (wf', rhs') -> CA.D.add_offspring wf' rhs')
          wf_rhs_list' fibered_dag.dag in
      { dag = dag';
        bundle = wf_bundle' }

(* O'Caml (correctly) infers the type
   [val colorize_dag : (D.node -> D.edge * D.children -> fibered_dag ->
                        (CA.D.node * (CA.D.edge * CA.D.children)) list * CWFBundle.t) ->
                       (D.node -> fibered_dag -> CA.D.node * CWFBundle.t) ->
                       D.t -> CWFBundle.t -> fibered_dag]. *)

    let colorize_dag f_node f_ext dag wf_bundle =
      A.D.fold (colorize_nodes f_node) dag
        (A.D.fold_nodes (colorize_sterile_nodes dag f_ext) dag
           { dag = CA.D.empty; bundle = wf_bundle })

    let colorize_external wf fibered_dag = 
      match CWFBundle.inv_pi wf fibered_dag.bundle with
      | [c_wf] -> (c_wf, fibered_dag.bundle)
      | [] -> failwith "colorize_external: not found"
      | _ -> failwith "colorize_external: not unique"

    let fuse_c_wf rhs =
      let momenta = PT.map (fun wf -> wf.CA.momentum) rhs in
      List.filter
        (fun (_, c) -> kmatrix_cuts c momenta)
        (CM.fuse (List.map (fun wf -> wf.CA.flavor) (PT.to_list rhs)))

    let colorize_coupling c coupling =
        { coupling with Tagged_Coupling.coupling = c }

    let colorize_fusion wf (coupling, children) fibered_dag =
      let match_flavor (f, _) = (CM.flavor_sans_color f = A.flavor wf)
      and find_colored wf' = CWFBundle.inv_pi wf' fibered_dag.bundle in
      let fusions =
        ThoList.flatmap
          (fun c_children ->
            List.map 
              (fun (f, c) ->
                (colorize_wf f wf, (colorize_coupling c coupling, c_children)))
              (List.filter match_flavor (fuse_c_wf c_children)))
          (PT.product (PT.map find_colored children)) in
      let bundle =
        List.fold_right
          (fun (c_wf, _) -> CWFBundle.add c_wf)
          fusions fibered_dag.bundle in
      (fusions, bundle)

    let colorize_braket1 (wf, (coupling, children)) fibered_dag =
      let find_colored wf' = CWFBundle.inv_pi wf' fibered_dag.bundle in
      Product.fold2
        (fun bra ket acc ->
          List.fold_left
            (fun brakets (f, c) ->
              if CM.conjugate bra.CA.flavor = f then
                (bra, (colorize_coupling c coupling, ket)) :: brakets
              else
                brakets)
            acc (fuse_c_wf ket))
        (find_colored wf) (PT.product (PT.map find_colored children)) []

    module CWFMap =
      Map.Make (struct type t = CA.wf let compare = CA.order_wf end)

    module CKetSet =
      Set.Make (struct type t = CA.rhs let compare = compare end)

    (* Find a set of kets in [map] that belong to [bra].
       Return the empty set, if nothing is found. *)

    let lookup_ketset bra map =
      try CWFMap.find bra map with Not_found -> CKetSet.empty

    (* Return the set of kets belonging to [bra] in [map],
       augmented by [ket]. *)

    let addto_ketset bra ket map =
      CKetSet.add ket (lookup_ketset bra map)

    (* Augment or update [map] with a new [(bra, ket)] relation. *)

    let addto_ketset_map map (bra, ket) =
      CWFMap.add bra (addto_ketset bra ket map) map

    (* Take a list of [(bra, ket)] pairs and group the [ket]s
       according to [bra].  This is very similar to
       [ThoList.factorize] on page~\pageref{ThoList.factorize},
       but the latter keeps duplicate copies, while we keep
       only one, with equality determined by [CA.order_wf]. *)

    (* \begin{dubious}
         Isn't [Bundle]~\ref{Bundle} the correct framework for this?
       \end{dubious} *)

    let factorize_brakets brakets =
      CWFMap.fold
        (fun bra ket acc -> (bra, CKetSet.elements ket) :: acc)
        (List.fold_left addto_ketset_map CWFMap.empty brakets)
        []

    let colorize_braket (wf, rhs_list) fibered_dag =
      factorize_brakets
        (ThoList.flatmap
           (fun rhs -> (colorize_braket1 (wf, rhs) fibered_dag))
           rhs_list)

    let colorize_amplitude a fin fout =
      let f = fin @ List.map CM.conjugate fout in
      let nin, nout = List.length fin, List.length fout in
      let n = nin + nout in
      let externals = List.combine f (ThoList.range 1 n) in
      let external_wfs = CA.external_wfs n externals in
      let wf_bundle = CWFBundle.of_list external_wfs  in

      let fibered_dag =
        colorize_dag
          colorize_fusion colorize_external a.A.fusion_dag wf_bundle in

      let brakets =
        ThoList.flatmap
          (fun braket -> colorize_braket braket fibered_dag)
          a.A.brakets in

      let dag = CA.D.harvest_list fibered_dag.dag (CA.wavefunctions brakets) in

      let fusions =
        List.filter (function (_, []) -> false | _ -> true) (CA.D.lists dag) in

      let dependencies_map =
        CA.D.fold
          (fun wf _ -> CWFMap.add wf (CA.D.dependencies dag wf))
          dag CWFMap.empty in

      { CA.fusions = fusions;
        CA.brakets = brakets;
        CA.constraints = a.A.constraints;
        CA.incoming = fin;
        CA.outgoing = fout;
        CA.externals = external_wfs;
        CA.fusion_dag = dag;
        CA.fusion_tower = dag; 
        CA.symmetry = a.A.symmetry;
        CA.on_shell = (fun wf -> a.A.on_shell (uncolorize_wf wf));
        CA.is_gauss = (fun wf -> a.A.is_gauss (uncolorize_wf wf));
        CA.dependencies = (fun wf -> CWFMap.find wf dependencies_map) }

    let allowed amplitude =
      match amplitude.CA.brakets with
      | [] -> false
      | _ -> true

    let colorize_amplitudes a =
      List.fold_left
        (fun amps (fin, fout) ->
          let amp = colorize_amplitude a fin fout in
          if allowed amp then
            amp :: amps
          else
            amps)
        [] (CM.amplitude a.A.incoming a.A.outgoing)

    let amplitudes goldstones exclusions selectors fin fout =
      colorize_amplitudes (amplitude goldstones selectors fin fout)

    let amplitude_sans_color goldstones exclusions selectors fin fout =
      amplitude goldstones selectors fin fout

    type flavor = CA.flavor
    type flavor_sans_color = A.flavor
    type p = A.p
    type wf = CA.wf
    let conjugate = CA.conjugate
    let flavor = CA.flavor
    let flavor_sans_color wf = CM.flavor_sans_color (CA.flavor wf)
    let momentum = CA.momentum
    let momentum_list = CA.momentum_list
    let wf_tag = CA.wf_tag

    type coupling = CA.coupling

    let sign = CA.sign
    let coupling = CA.coupling
    let coupling_tag = CA.coupling_tag
    type exclusions = CA.exclusions
    let no_exclusions = CA.no_exclusions

    type 'a children = 'a CA.children
    type rhs = CA.rhs
    let children = CA.children

    type fusion = CA.fusion
    let lhs = CA.lhs
    let rhs = CA.rhs

    type braket = CA.braket
    let bra = CA.bra
    let ket = CA.ket   

    type amplitude = CA.amplitude
    type amplitude_sans_color = A.amplitude
    let incoming = CA.incoming
    let outgoing = CA.outgoing
    let externals = CA.externals
    let fusions = CA.fusions
    let brakets = CA.brakets
    let symmetry = CA.symmetry
    let on_shell = CA.on_shell
    let is_gauss = CA.is_gauss
    let constraints = CA.constraints
    let variables a = List.map lhs (fusions a)
    let dependencies = CA.dependencies

(* \thocwmodulesubsection{Checking Conservation Laws} *)

    let check_charges () =
      let vlist3, vlist4, vlistn = M.vertices () in
      List.filter
        (fun flist -> not (M.Ch.is_null (M.Ch.sum (List.map M.charges flist))))
        (List.map (fun ((f1, f2, f3), _, _) -> [f1; f2; f3]) vlist3
         @ List.map (fun ((f1, f2, f3, f4), _, _) -> [f1; f2; f3; f4]) vlist4
         @ List.map (fun (flist, _, _) -> flist) vlistn)

(* \thocwmodulesubsection{Diagnostics} *)

    let count_propagators a =
      List.length a.CA.fusions

    let count_fusions a =
      List.fold_left (fun n (_, a) -> n + List.length a) 0 a.CA.fusions
        + List.fold_left (fun n (_, t) -> n + List.length t) 0 a.CA.brakets
        + List.length a.CA.brakets

(* \begin{dubious}
     This brute force approach blows up for more than ten particles.
     Find a smarter algorithm.
   \end{dubious} *)

    let count_diagrams a =
      List.fold_left (fun n (wf1, wf23) ->
        n + CA.D.count_trees wf1 a.CA.fusion_dag *
          (List.fold_left (fun n' (_, wfs) ->
            n' + PT.fold_left (fun n'' wf ->
              n'' * CA.D.count_trees wf a.CA.fusion_dag) 1 wfs) 0 wf23))
        0 a.CA.brakets

    exception Impossible

    let forest' a =
      let below wf = CA.D.forest_memoized wf a.CA.fusion_dag in
      ThoList.flatmap
        (fun (bra, ket) ->
          (Product.list2 (fun bra' ket' -> bra' :: ket')
             (below bra)
             (ThoList.flatmap
                (fun (_, wfs) ->
                  Product.list (fun w -> w) (PT.to_list (PT.map below wfs)))
                ket)))
        a.CA.brakets

    let cross wf =
      { CA.flavor = CM.conjugate wf.CA.flavor;
        CA.momentum = P.neg wf.CA.momentum;
        CA.wf_tag = wf.CA.wf_tag }

    let fuse_trees wf ts =
      Tree.fuse (fun (wf', e) -> (cross wf', e))
        wf (fun t -> List.mem wf (Tree.leafs t)) ts
      
    let forest wf a =
      List.map (fuse_trees wf) (forest' a)

(*i
(* \begin{dubious}
     The following duplication should be replaced by polymorphism
     or a functor.
   \end{dubious} *)

    let forest_uncolored' a =
      let below wf = A.D.forest_memoized wf a.A.fusion_dag in
      ThoList.flatmap
        (fun (bra, ket) ->
          (Product.list2 (fun bra' ket' -> bra' :: ket')
             (below bra)
             (ThoList.flatmap
                (fun (_, wfs) ->
                  Product.list (fun w -> w) (PT.to_list (PT.map below wfs)))
                ket)))
        a.A.brakets

    let cross_uncolored wf =
      { A.flavor = M.conjugate wf.A.flavor;
        A.momentum = P.neg wf.A.momentum;
        A.wf_tag = wf.A.wf_tag }

    let fuse_trees_uncolored wf ts =
      Tree.fuse (fun (wf', e) -> (cross_uncolored wf', e))
        wf (fun t -> List.mem wf (Tree.leafs t)) ts
      
    let forest_sans_color wf a =
      List.map (fuse_trees_uncolored wf) (forest_uncolored' a)
i*)

    let poles_beneath wf dag =
      CA.D.eval_memoized (fun wf' -> [[]])
        (fun wf' _ p -> List.map (fun p' -> wf' :: p') p)
        (fun wf1 wf2 ->
          Product.fold2 (fun wf' wfs' wfs'' -> (wf' @ wfs') :: wfs'') wf1 wf2 [])
        (@) [[]] [[]] wf dag

    let poles a =
      ThoList.flatmap (fun (wf1, wf23) ->
        let poles_wf1 = poles_beneath wf1 a.CA.fusion_dag in
        (ThoList.flatmap (fun (_, wfs) ->
          Product.list List.flatten
            (PT.to_list (PT.map (fun wf ->
              poles_wf1 @ poles_beneath wf a.CA.fusion_dag) wfs)))
           wf23))
        a.CA.brakets

    module WFSet =
      Set.Make (struct type t = CA.wf let compare = CA.order_wf end)

    let s_channel a =
      WFSet.elements
        (ThoList.fold_right2
           (fun wf wfs ->
             if P.Scattering.timelike wf.CA.momentum then
               WFSet.add wf wfs
             else
               wfs) (poles a) WFSet.empty)
      
(* \begin{dubious}
     This should be much faster!  Is it correct?  Is it faster indeed?
   \end{dubious} *)

    let poles' a =
      List.map CA.lhs a.CA.fusions

    let s_channel a =
      WFSet.elements
        (List.fold_right
           (fun wf wfs ->
             if P.Scattering.timelike wf.CA.momentum then
               WFSet.add wf wfs
             else
               wfs) (poles' a) WFSet.empty)
      
(* \thocwmodulesubsection{Pictures} *)

(* Export the DAG in the \texttt{dot(1)} file format so that we can
   draw pretty pictures to impress audiences \ldots *)

    let p2s p =
      if p >= 0 && p <= 9 then
        string_of_int p
      else if p <= 36 then
        String.make 1 (Char.chr (Char.code 'A' + p - 10))
      else
        "_"

    let variable wf =
      CM.flavor_symbol wf.CA.flavor ^
      String.concat "" (List.map p2s (P.to_ints wf.CA.momentum))

    module Int = Map.Make (struct type t = int let compare = compare end)

    let add_to_list i n m =
      Int.add i (n :: try Int.find i m with Not_found -> []) m

    let classify_nodes dag =
      Int.fold (fun i n acc -> (i, n) :: acc)
        (CA.D.fold_nodes (fun wf -> add_to_list (P.rank wf.CA.momentum) wf)
           dag Int.empty) []

    let dag_to_dot ch brakets dag =
      Printf.fprintf ch "digraph OMEGA {\n";
      CA.D.iter_nodes (fun wf ->
        Printf.fprintf ch "  \"%s\" [ label = \"%s\" ];\n"
          (variable wf) (variable wf)) dag;
      List.iter (fun (_, wfs) ->
        Printf.fprintf ch "  { rank = same;";
        List.iter (fun n ->
          Printf.fprintf ch " \"%s\";" (variable n)) wfs;
        Printf.fprintf ch " };\n") (classify_nodes dag);
      List.iter (fun n ->
        Printf.fprintf ch " \"*\" -> \"%s\";\n" (variable n))
        (flatten_keystones brakets);
      CA.D.iter (fun n (_, ns) ->
        let p = variable n in
        PT.iter (fun n' ->
          Printf.fprintf ch "  \"%s\" -> \"%s\";\n" p (variable n')) ns) dag;
      Printf.fprintf ch "}\n"

    let tower_to_dot ch a =
      dag_to_dot ch a.CA.brakets a.CA.fusion_tower

    let amplitude_to_dot ch a =
      dag_to_dot ch a.CA.brakets a.CA.fusion_dag

(* \thocwmodulesubsection{Phasespace} *)


    let variable wf =
      M.flavor_to_string wf.A.flavor ^
        "[" ^ String.concat "/" (List.map p2s (P.to_ints wf.A.momentum)) ^ "]"

    let below_to_channel transform ch dag wf =
      let n2s wf = variable (transform wf)
      and e2s c = "" in
      Tree2.to_channel ch n2s e2s (A.D.dependencies dag wf)

    let bra_to_channel transform ch dag wf =
      let tree = A.D.dependencies dag wf in
      if Tree2.is_singleton tree then
        let n2s wf = variable (transform wf)
        and e2s c = "" in
        Tree2.to_channel ch n2s e2s tree
      else
        failwith "Fusion.phase_space_channels: wrong topology!"

    let ket_to_channel transform ch dag ket =
      Printf.fprintf ch "(";
      begin match A.children ket with
      | [] -> ()
      | [child] -> below_to_channel transform ch dag child
      | child :: children ->
         below_to_channel transform ch dag child;
         List.iter
           (fun child ->
             Printf.fprintf ch ",";
             below_to_channel transform ch dag child)
           children
      end;
      Printf.fprintf ch ")"

    let phase_space_braket transform ch (bra, ket) dag =
      bra_to_channel transform ch dag bra;
      Printf.fprintf ch ": {";
      begin match ket with
      | [] -> ()
      | [ket1] ->
         Printf.fprintf ch " ";
         ket_to_channel transform ch dag ket1
      | ket1 :: kets ->
         Printf.fprintf ch " ";
         ket_to_channel transform ch dag ket1;
         List.iter
           (fun k ->
             Printf.fprintf ch " \\\n   | ";
             ket_to_channel transform ch dag k)
           kets
      end;
      Printf.fprintf ch " }\n"

(*i Food for thought:

    let braket_to_tree2 dag (bra, ket) =
      let bra' = A.D.dependencies dag bra in
      if Tree2.is_singleton bra' then
        Tree2.cons
          [(fst ket, bra, List.map (A.D.dependencies dag) (A.children ket))]
      else
        failwith "Fusion.phase_space_channels: wrong topology!"

    let phase_space_braket transform ch (bra, ket) dag =
      let n2s wf = variable (transform wf)
      and e2s c = "" in
      Printf.fprintf
        ch "%s\n" (Tree2.to_string n2s e2s (braket_to_tree2 dag (bra, ket)))
i*)

    let phase_space_channels_transformed transform ch a =
      List.iter
        (fun braket -> phase_space_braket transform ch braket a.A.fusion_dag)
        a.A.brakets

    let phase_space_channels ch a =
      phase_space_channels_transformed (fun wf -> wf) ch a

    let exchange_momenta_list p1 p2 p =
      List.map
        (fun pi ->
          if pi = p1 then
            p2
          else if pi = p2 then
            p1
          else
            pi)
        p

    let exchange_momenta p1 p2 p =
      P.of_ints (P.dim p) (exchange_momenta_list p1 p2 (P.to_ints p))

    let flip_momenta wf =
      { wf with A.momentum = exchange_momenta 1 2 wf.A.momentum }

    let phase_space_channels_flipped ch a =
      phase_space_channels_transformed flip_momenta ch a

  end

module Make = Tagged(No_Tags)

module Binary = Make(Tuple.Binary)(Stat_Dirac)(Topology.Binary)
module Tagged_Binary (T : Tagger) =
  Tagged(T)(Tuple.Binary)(Stat_Dirac)(Topology.Binary)

(* \thocwmodulesection{Fusions with Majorana Fermions} *)

let majorana_log silent logging = logging
let majorana_log silent logging = silent
let force_legacy = true
let force_legacy = false

module Stat_Majorana (M : Model.T) : (Stat with type flavor = M.flavor) =
  struct 

    exception Impossible

    type flavor = M.flavor

    (* \thocwmodulesubsection{Keeping Track of Fermion Lines} *)

    (* JRR's algorithm doesn't use lists of pairs representing
       directed arrows as in [Stat_Dirac().stat] above, but a list
       of integers denoting the external leg a fermion line connects
       to: *)
    type stat =
      | Fermion of int * int list
      | AntiFermion of int * int list
      | Boson of int list
      | Majorana of int * int list        

    let sign_of_permutation lines = fst (Combinatorics.sort_signed lines)   

    let lines_equivalent l1 l2 =
      sign_of_permutation l1 = sign_of_permutation l2
      
    let stat_to_string s =
      let open Printf in
      let l2s = ThoList.to_string string_of_int in
      match s with
      | Boson lines -> sprintf "B%s" (l2s lines)
      | Fermion (p, lines) -> sprintf "F(%d, %s)" p (l2s lines)
      | AntiFermion (p, lines) -> sprintf "A(%d, %s)" p (l2s lines)
      | Majorana (p, lines) -> sprintf "M(%d, %s)" p (l2s lines)

    (* Writing all cases explicitely is tedious, but allows exhaustiveness
       checking.  *)
    let equal s1 s2 =
      match s1, s2 with
      | Boson l1, Boson l2 ->
         lines_equivalent l1 l2
      | Majorana (p1, l1), Majorana (p2, l2)
      | Fermion (p1, l1), Fermion (p2, l2)
      | AntiFermion (p1, l1), AntiFermion (p2, l2) ->
         p1 = p2 && lines_equivalent l1 l2
      | Boson _, (Fermion _ | AntiFermion _ | Majorana _ )
      | (Fermion _ | AntiFermion _ | Majorana _ ), Boson _
      | Majorana _, (Fermion _ | AntiFermion _)
      | (Fermion _ | AntiFermion _), Majorana _
      | Fermion _ , AntiFermion _
      | AntiFermion _ , Fermion _ -> false

    (* The final amplitude must not be fermionic! *)
    let saturated = function
      | Boson _ -> true
      | Fermion _ | AntiFermion _ | Majorana _ -> false

    (* [stat f p] interprets the numeric fermion numbers of flavor [f]
       at external leg [p] at creates a leaf: *)
    let stat f p =
      match M.fermion f with
      | 0 -> Boson []
      | 1 -> Fermion (p, [])
      | -1 -> AntiFermion (p, [])
      | 2 -> Majorana (p, [])
      | _ -> invalid_arg "Fusion.Stat_Majorana: invalid fermion number"

(* The formalism of~\cite{Denner:Majorana} does not distinguish
   spinors from conjugate spinors, it is only important to know in which direction
   a fermion line is calculated. So the sign is made by the calculation together
   with an aditional one due to the permuation of the pairs of endpoints of
   fermion lines in the direction they are calculated. We propose a
   ``canonical'' direction from the right to the left child at a fusion point
   so we only have to keep in mind which external particle hangs at each side.
   Therefore we need not to have a list of pairs of conjugate spinors and
   spinors but just a list in which the pairs are right-left-right-left
   and so on. Unfortunately it is unavoidable to have couplings with clashing 
   arrows in supersymmetric theories so we need transmutations from fermions 
   in antifermions and vice versa as well. *)   

    (* \thocwmodulesubsection{Merge Fermion Lines for Legacy Models with Implied Fermion Connections} *)

    (* In the legacy case with at most one fermion line, it was straight
       forward to determine the kind of outgoing line from the 
       corresponding flavor.  In the general case, it is not
       possible to maintain this constraint, when constructing
       the $n$-ary fusion from binary ones. *)

    (* We can break up the process into two steps however:
       first perform unconstrained fusions pairwise \ldots *)

    let stat_fuse_pair_unconstrained s1 s2 =
      match s1, s2 with
      | Boson l1, Boson l2 -> Boson (l1 @ l2)
      | (Majorana (p1, l1) | Fermion (p1, l1) | AntiFermion (p1, l1)),
        (Majorana (p2, l2) | Fermion (p2, l2) | AntiFermion (p2, l2)) ->
          Boson ([p2; p1] @ l1 @ l2)
      | Boson l1, Majorana (p, l2) -> Majorana (p, l1 @ l2)
      | Boson l1, Fermion (p, l2)  -> Fermion (p, l1 @ l2)
      | Boson l1, AntiFermion (p, l2) -> AntiFermion (p, l1 @ l2)
      | Majorana (p, l1), Boson l2 -> Majorana (p, l1 @ l2)
      | Fermion (p, l1), Boson l2 -> Fermion (p, l1 @ l2)
      | AntiFermion (p, l1), Boson l2 -> AntiFermion (p, l1 @ l2)

    (* \ldots{} and only apply the constraint to the outgoing leg. *)

    let constrain_stat_fusion s f =
      match s, M.lorentz f with
      | (Majorana (p, l) | Fermion (p, l) | AntiFermion (p, l)),
        (Coupling.Majorana | Coupling.Vectorspinor | Coupling.Maj_Ghost) ->
         Majorana (p, l)
      | (Majorana (p, l) | Fermion (p, l) | AntiFermion (p, l)),
        Coupling.Spinor -> Fermion (p, l)
      | (Majorana (p, l) | Fermion (p, l) | AntiFermion (p, l)),
        Coupling.ConjSpinor -> AntiFermion (p, l)
      | (Majorana _ | Fermion _ | AntiFermion _ as s),
        (Coupling.Scalar | Coupling.Vector | Coupling.Massive_Vector
         | Coupling.Tensor_1 | Coupling.Tensor_2 | Coupling.BRS _) ->
         invalid_arg
           (Printf.sprintf
              "Fusion.stat_fuse_pair_constrained: expected boson, got %s"
              (stat_to_string s))
      | Boson l as s,
        (Coupling.Majorana | Coupling.Vectorspinor | Coupling.Maj_Ghost
         | Coupling.Spinor | Coupling.ConjSpinor) ->
         invalid_arg
           (Printf.sprintf
              "Fusion.stat_fuse_pair_constrained: expected fermion, got %s"
              (stat_to_string s))
      | Boson l,
        (Coupling.Scalar | Coupling.Vector | Coupling.Massive_Vector
         | Coupling.Tensor_1 | Coupling.Tensor_2 | Coupling.BRS _) ->
         Boson l

    let stat_fuse_pair_legacy f s1 s2 =
      stat_fuse_pair_unconstrained s1 s2

    let stat_fuse_pair_legacy_logging f s1 s2 =
      let stat = stat_fuse_pair_legacy f s1 s2 in
      Printf.eprintf
        "stat_fuse_pair_legacy: (%s, %s) -> %s = %s\n"
        (stat_to_string s1) (stat_to_string s2) (stat_to_string stat)
        (M.flavor_to_string f);
      stat

    let stat_fuse_pair_legacy =
      majorana_log stat_fuse_pair_legacy stat_fuse_pair_legacy_logging

    (* Note that we are using [List.fold_left], therefore
       we perform the fusions as
       $f(f(\ldots(f(s_1,s_2),s_3),\ldots),s_n)$.  Had we used
       [List.fold_right] instead, we would compute
       $f(s_1,f(s_2,\ldots f(s_{n-1},s_n))).$   For our Dirac
       algorithm, this makes no difference, but JRR's Majorana
       algorithm depends on the order! *)

    (* Also not that we \emph{must not} apply [constrain_stat_fusion]
       here, because [stat_fuse_legacy] will be used in
       [stat_keystone_legacy] again, where we always expect
       [Boson _]. *)
    let stat_fuse_legacy s1 s23__n f =
      List.fold_left (stat_fuse_pair_legacy f) s1 s23__n

    (*i
    let stat_fuse_legacy' s1 s23__n f =
      match List.rev (s1 :: s23__n) with
      | s1 :: s23__n -> List.fold_left (stat_fuse_pair_legacy f) s1 s23__n
      | [] -> failwith "stat_fuse_legacy: impossible"

    let stat_fuse_legacy' s1 s23__n f =
      List.fold_right (stat_fuse_pair_legacy f) s23__n s1
i*)

    let stat_fuse_legacy_logging s1 s23__n f =
      let stat = stat_fuse_legacy s1 s23__n f in
      Printf.eprintf
        "stat_fuse_legacy:      %s -> %s = %s\n"
        (ThoList.to_string stat_to_string (s1 :: s23__n))
        (stat_to_string stat)
        (M.flavor_to_string f);
      stat

    let stat_fuse_legacy =
      majorana_log stat_fuse_legacy stat_fuse_legacy_logging

    (* \thocwmodulesubsection{Merge Fermion Lines using Explicit Fermion Connections} *)

    (* We need to match the fermion lines in the incoming propagators
       using the connection information in the vertex.  This used to
       be trivial in the old omega, because there was at most one
       fermion line in a vertex. *)
    module IMap = Map.Make (struct type t = int let compare = compare end)

    (* From version 4.05 on, this is just [IMap.find_opt]. *)
    let imap_find_opt p map =
      try Some (IMap.find p map) with Not_found -> None

    (* Partially combined [stat]s of the incoming propagators and keeping
       track of the fermion lines, while we're scanning them. *)
    type partial =
      { stat : stat (* the [stat] accumulated so far *);
        fermions : int IMap.t (* a map from the indices in the vertex to open (anti)fermion lines *);
        n : int (* the number of incoming propagators *) }

    (* We will
       perform two passes:
       \begin{enumerate}
         \item collect the saturated fermion lines in a [Boson], while
           building a map from the indices in the vertex to the open
           fermion lines
         \item connect the open fermion lines using the [int -> int] map
           [fermions].
       \end{enumerate} *)

    let empty_partial =
      { stat = Boson [];
        fermions = IMap.empty;
        n = 0 }

    (* Only for debugging: *)
    let partial_to_string p =
      Printf.sprintf
        "{ fermions=%s, stat=%s, #=%d }"
        (ThoList.to_string
           (fun (i, particle) -> Printf.sprintf "%d@%d" particle i)
           (IMap.bindings p.fermions))
        (stat_to_string p.stat)
        p.n

    (* Add a list of saturated fermion lines at the top of the list
       of lines in a [stat]. *)
    let add_lines l = function
      | Boson l' -> Boson (l @ l')
      | Fermion (n, l') -> Fermion (n, l @ l')
      | AntiFermion (n, l') -> AntiFermion (n, l @ l')
      | Majorana (n, l') -> Majorana (n, l @ l')

    (* Process one line in the first pass: add the saturated fermion lines
       to the partial stat [p.stat]
       and add a pointer to an open fermion line in case of a fermion. *)
    let add_lines_to_partial p stat =
      let n = succ p.n in
      match stat with
      | Boson l ->
         { fermions = p.fermions;
           stat = add_lines l p.stat;
           n }
      | Majorana (f, l) ->
         { fermions = IMap.add n f p.fermions;
           stat = add_lines l p.stat;
           n }
      | Fermion (p, l) ->
         invalid_arg
           "add_lines_to_partial: unexpected Fermion"
      | AntiFermion (p, l) ->
         invalid_arg
           "add_lines_to_partial: unexpected AntiFermion"

    (* Do it for all lines: *)
    let partial_of_slist stat_list =
      List.fold_left add_lines_to_partial empty_partial stat_list

    let partial_of_rev_slist stat_list =
      List.fold_left add_lines_to_partial empty_partial (List.rev stat_list)

    (* The building blocks for a single step of the second pass:
       saturate a fermion line or pass it through. *)

    (* The indices [i] and [j] refer to incoming lines: add a saturated
       line to [p.stat] and remove the corresponding open lines from
       the map. *)
    let saturate_fermion_line p i j =
      match imap_find_opt i p.fermions, imap_find_opt j p.fermions with
      | Some f, Some f' ->
         { stat = add_lines [f'; f] p.stat;
           fermions = IMap.remove i (IMap.remove j p.fermions);
           n = p.n }
      | Some _, None ->
         invalid_arg "saturate_fermion_line: no open outgoing fermion line"
      | None, Some _ ->
         invalid_arg "saturate_fermion_line: no open incoming fermion line"
      | None, None ->
         invalid_arg "saturate_fermion_line: no open fermion lines"

    (* The index [i] refers to an incoming line: add the open line
       to [p.stat] and remove it from the map. *)
    let pass_through_fermion_line p i =
      match imap_find_opt i p.fermions, p.stat with
      | Some f, Boson l ->
         { stat = Majorana (f, l);
           fermions = IMap.remove i p.fermions;
           n = p.n }
      | Some _ , (Majorana _ | Fermion _ | AntiFermion _) ->
         invalid_arg "pass_through_fermion_line: more than one open line"
      | None, _ ->
         invalid_arg "pass_through_fermion_line: expected fermion not found"

    (* Ignoring the direction of the fermion line reproduces JRR's algorithm. *)
    let sort_pair (i, j) =
      if i < j then
        (i, j)
      else
        (j, i)

    (* The index [p.n + 1] corresponds to the outgoing line: *)
    let is_incoming p i =
      i <= p.n

    let match_fermion_line p (i, j) =
      let i, j = sort_pair (i, j) in
      if is_incoming p i && is_incoming p j then
        saturate_fermion_line p i j
      else if is_incoming p i then
        pass_through_fermion_line p i
      else if is_incoming p j then
        pass_through_fermion_line p j
      else
        failwith "match_fermion_line: both lines outgoing"

    let match_fermion_line_logging p (i, j) =
      Printf.eprintf
        "match_fermion_line     %s [%d->%d]"
        (partial_to_string p) i j;
      let p' = match_fermion_line p (i, j) in
      Printf.eprintf " >> %s\n" (partial_to_string p');
      p'

    let match_fermion_line =
      majorana_log match_fermion_line match_fermion_line_logging

    (* Combine the passes \ldots *)
    let match_fermion_lines flines s1 s23__n =
      List.fold_left match_fermion_line (partial_of_slist (s1 :: s23__n)) flines

    (* \ldots{} and keep only the [stat]. *)
    let stat_fuse_new flines s1 s23__n _ =
      (match_fermion_lines flines s1 s23__n).stat

    (* If there is at most a single fermion line, we can compare [stat]
       against the result of [stat_fuse_legacy] for checking
       [stat_fuse_new] (admittedly, this case is rather trivial) \ldots *)
    let stat_fuse_new_check stat flines s1 s23__n f =
      if List.length flines < 2 then
        begin
          let legacy = stat_fuse_legacy s1 s23__n f in
          if not (equal stat legacy) then
            failwith
              (Printf.sprintf
                 "stat_fuse_new: %s <> %s!"
                 (stat_to_string stat)
                 (stat_to_string legacy))
        end

    (* \ldots{} do it, but only when we are writing debugging output. *)
    let stat_fuse_new_logging flines s1 s23__n f =
      let stat = stat_fuse_new flines s1 s23__n f in
      Printf.eprintf
        "stat_fuse_new: %s: %s -> %s = %s\n"
        (UFO_Lorentz.fermion_lines_to_string flines)
        (ThoList.to_string stat_to_string (s1 :: s23__n))
        (stat_to_string stat)
        (M.flavor_to_string f);
      stat_fuse_new_check stat flines s1 s23__n f;
      stat

    let stat_fuse_new =
      majorana_log stat_fuse_new stat_fuse_new_logging

    (* Use [stat_fuse_new], whenever fermion connections are
       available.  NB: [Some []] is \emph{not} the same as [None]! *)
    let stat_fuse flines_opt slist f =
      match slist with
      | [] -> invalid_arg "stat_fuse: empty"
      | s1 :: s23__n ->
         constrain_stat_fusion
           (match flines_opt with
            | Some flines -> stat_fuse_new flines s1 s23__n f
            | None -> stat_fuse_legacy s1 s23__n f)
           f

    let stat_fuse_logging flines_opt slist f =
      let stat = stat_fuse flines_opt slist f in
      Printf.eprintf
        "stat_fuse:             %s -> %s = %s\n"
        (ThoList.to_string stat_to_string slist)
        (stat_to_string stat)
        (M.flavor_to_string f);
      stat

    let stat_fuse =
      majorana_log stat_fuse stat_fuse_logging

    (* \thocwmodulesubsection{Final Step using Implied Fermion Connections} *)

    let stat_keystone_legacy s1 s23__n f =
      stat_fuse_legacy s1 s23__n f

    let stat_keystone_legacy_logging s1 s23__n f =
      let s = stat_keystone_legacy s1 s23__n f in
      Printf.eprintf
        "stat_keystone_legacy: %s (%s) %s -> %s\n"
        (stat_to_string s1)
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string s23__n)
        (stat_to_string s);
      s

    let stat_keystone_legacy =
      majorana_log stat_keystone_legacy stat_keystone_legacy_logging

    (* \thocwmodulesubsection{Final Step using Explicit Fermion Connections} *)

    let stat_keystone_new flines slist f =
      match slist with
      | [] -> invalid_arg "stat_keystone: empty"
      | [s] -> invalid_arg "stat_keystone: singleton"
      | s1 :: s2 :: s34__n ->
         let stat =
           stat_fuse_pair_unconstrained s1 (stat_fuse_new flines s2 s34__n f) in
         if saturated stat then
           stat
         else
           failwith
             (Printf.sprintf
                "stat_keystone: incomplete %s!"
                (stat_to_string stat))

    let stat_keystone_new_check stat slist f =
      match slist with
      | [] -> invalid_arg "stat_keystone_check: empty"
      | s1 :: s23__n ->
         let legacy = stat_keystone_legacy s1 s23__n f in
         if not (equal stat legacy) then
           failwith
             (Printf.sprintf
                "stat_keystone_check: %s <> %s!"
                (stat_to_string stat)
                (stat_to_string legacy))

    let stat_keystone flines_opt slist f =
      match flines_opt with
      | Some flines -> stat_keystone_new flines slist f
      | None ->
         begin match slist with
         | [] -> invalid_arg "stat_keystone: empty"
         | s1 :: s23__n -> stat_keystone_legacy s1 s23__n f
         end

    let stat_keystone_logging flines_opt slist f =
      let stat = stat_keystone flines_opt slist f in
      Printf.eprintf
        "stat_keystone:        %s (%s) %s -> %s\n"
        (stat_to_string (List.hd slist))
        (M.flavor_to_string f)
        (ThoList.to_string stat_to_string (List.tl slist))
        (stat_to_string stat);
      stat_keystone_new_check stat slist f;
      stat

    let stat_keystone =
      majorana_log stat_keystone stat_keystone_logging

    (* Force the legacy version w/o checking against the
       new implementation for comparing generated code
       against the hard coded models: *)

    let stat_fuse flines_opt slist f =
      if force_legacy then
        stat_fuse_legacy (List.hd slist) (List.tl slist) f
      else
        stat_fuse flines_opt slist f

    let stat_keystone flines_opt slist f =
      if force_legacy then
        stat_keystone_legacy (List.hd slist) (List.tl slist) f
      else
        stat_keystone flines_opt slist f

    (* \thocwmodulesubsection{Evaluate Signs from Fermion Permuations} *)

    let stat_sign = function
      | Boson lines -> sign_of_permutation lines
      | Fermion (p, lines) -> sign_of_permutation (p :: lines)
      | AntiFermion (pbar, lines) -> sign_of_permutation (pbar :: lines)
      | Majorana (pm, lines) -> sign_of_permutation (pm :: lines)  

    let stat_sign_logging stat =
      let sign = stat_sign stat in
      Printf.eprintf
        "stat_sign: %s -> %d\n"
        (stat_to_string stat) sign;
      sign

    let stat_sign =
      majorana_log stat_sign stat_sign_logging

  end

module Binary_Majorana =
  Make(Tuple.Binary)(Stat_Majorana)(Topology.Binary)

module Nary (B: Tuple.Bound) =
  Make(Tuple.Nary(B))(Stat_Dirac)(Topology.Nary(B))
module Nary_Majorana (B: Tuple.Bound) =
  Make(Tuple.Nary(B))(Stat_Majorana)(Topology.Nary(B))

module Mixed23 =
  Make(Tuple.Mixed23)(Stat_Dirac)(Topology.Mixed23)
module Mixed23_Majorana =
  Make(Tuple.Mixed23)(Stat_Majorana)(Topology.Mixed23)

module Helac (B: Tuple.Bound) =
  Make(Tuple.Nary(B))(Stat_Dirac)(Topology.Helac(B))
module Helac_Majorana (B: Tuple.Bound) =
  Make(Tuple.Nary(B))(Stat_Majorana)(Topology.Helac(B))

module B2 = struct let max_arity () = 2 end
module B3 = struct let max_arity () = 3 end
module Helac_Binary = Helac(B2)
module Helac_Binary_Majorana = Helac(B2)
module Helac_Mixed23 = Helac(B3)
module Helac_Mixed23_Majorana = Helac(B3)

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
    val amplitudes : bool -> int option ->
      exclusions -> selectors -> process list -> amplitudes
    val empty : amplitudes
(*i
    val initialize_cache : string -> unit
    val set_cache_name : string -> unit
i*)
    val flavors : amplitudes -> process list
    val vanishing_flavors : amplitudes -> process list
    val color_flows : amplitudes -> Color.Flow.t list
    val helicities : amplitudes -> (int list * int list) list
    val processes : amplitudes -> amplitude list
    val process_table : amplitudes -> amplitude option array array
    val fusions : amplitudes -> (fusion * amplitude) list
    val multiplicity : amplitudes -> wf -> int
    val dictionary : amplitudes -> amplitude -> wf -> int
    val color_factors : amplitudes -> Color.Flow.factor array array
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

module Multi (Fusion_Maker : Maker) (P : Momentum.T) (M : Model.T) =
  struct

    exception Mismatch

    type progress_mode =
      | Quiet
      | Channel of out_channel
      | File of string

    let progress_option = ref Quiet

    module CM = Colorize.It(M)
    module F = Fusion_Maker(P)(M)
    module C = Cascade.Make(M)(P)

(* \begin{dubious}
     A kludge, at best \ldots
   \end{dubious} *)

    let options = Options.extend F.options
        [ "progress", Arg.Unit (fun () -> progress_option := Channel stderr),
          "report progress to the standard error stream";
          "progress_file", Arg.String (fun s -> progress_option := File s),
          "report progress to a file" ]

    type flavor = M.flavor
    type p = F.p
    type process = flavor list * flavor list
    type amplitude = F.amplitude
    type fusion = F.fusion
    type wf = F.wf
    type exclusions = F.exclusions
    let no_exclusions = F.no_exclusions
    type selectors = F.selectors

    type flavors = flavor list array
    type helicities = int list array
    type colors = Color.Flow.t array

    type amplitudes' = amplitude array array array

    type amplitudes =
        { flavors : process list;
          vanishing_flavors : process list;
          color_flows : Color.Flow.t list;
          helicities : (int list * int list) list; 
          processes : amplitude list;
          process_table : amplitude option array array;
          fusions : (fusion * amplitude) list;
          multiplicity : (wf -> int);
          dictionary : (amplitude -> wf -> int);
          color_factors : Color.Flow.factor array array;
          constraints : string option }

    let flavors a = a.flavors
    let vanishing_flavors a = a.vanishing_flavors
    let color_flows a = a.color_flows
    let helicities a = a.helicities
    let processes a = a.processes
    let process_table a = a.process_table
    let fusions a = a.fusions
    let multiplicity a = a.multiplicity
    let dictionary a = a.dictionary
    let color_factors a = a.color_factors
    let constraints a = a.constraints

    let sans_colors f =
      List.map CM.flavor_sans_color f

    let colors (fin, fout) =
      List.map M.color (fin @ fout)

    let process_sans_color a =
      (sans_colors (F.incoming a), sans_colors (F.outgoing a))

    let color_flow a =
      CM.flow (F.incoming a) (F.outgoing a)

    let process_to_string fin fout =
      String.concat " " (List.map M.flavor_to_string fin)
      ^ " -> " ^ String.concat " " (List.map M.flavor_to_string fout)

    let count_processes colored_processes =
      List.length colored_processes

    module FMap =
      Map.Make (struct type t = process let compare = compare end)

    module CMap =
      Map.Make (struct type t = Color.Flow.t let compare = compare end)

(* Recently [Product.list] began to guarantee lexicographic order for sorted
   arguments.  Anyway, we still force a lexicographic order. *)

    let rec order_spin_table1 s1 s2 =
      match s1, s2 with
      | h1 :: t1, h2 :: t2 ->
          let c = compare h1 h2 in
          if c <> 0 then
            c
          else
            order_spin_table1 t1 t2
      | [], [] -> 0
      | _ -> invalid_arg "order_spin_table: inconsistent lengths"
      
    let order_spin_table (s1_in, s1_out) (s2_in, s2_out) =
      let c = compare s1_in s2_in in
      if c <> 0 then
        c
      else
        order_spin_table1 s1_out s2_out
          
    let sort_spin_table table =
      List.sort order_spin_table table

    let id x = x

    let pair x y = (x, y)

(* \begin{dubious}
     Improve support for on shell Ward identities: [Coupling.Vector -> [4]] for one
     and only one external vector.
   \end{dubious} *)

    let rec hs_of_lorentz = function
      | Coupling.Scalar -> [0]
      | Coupling.Spinor | Coupling.ConjSpinor
      | Coupling.Majorana | Coupling.Maj_Ghost -> [-1; 1]
      | Coupling.Vector -> [-1; 1]
      | Coupling.Massive_Vector -> [-1; 0; 1]
      | Coupling.Tensor_1 -> [-1; 0; 1]
      | Coupling.Vectorspinor -> [-2; -1; 1; 2]
      | Coupling.Tensor_2 -> [-2; -1; 0; 1; 2]
      | Coupling.BRS f -> hs_of_lorentz f

    let hs_of_flavor f =
      hs_of_lorentz (M.lorentz f)

    let hs_of_flavors (fin, fout) =
      (List.map hs_of_flavor fin, List.map hs_of_flavor fout)

    let rec unphysical_of_lorentz = function
      | Coupling.Vector -> [4]
      | Coupling.Massive_Vector -> [4]
      | _ -> invalid_arg "unphysical_of_lorentz: not a vector particle"

    let unphysical_of_flavor f =
      unphysical_of_lorentz (M.lorentz f)

    let unphysical_of_flavors1 n f_list =
      ThoList.mapi
        (fun i f -> if i = n then unphysical_of_flavor f else hs_of_flavor f)
        1 f_list
      
    let unphysical_of_flavors n (fin, fout) =
      (unphysical_of_flavors1 n fin, unphysical_of_flavors1 (n - List.length fin) fout)

    let helicity_table unphysical flavors =
      let hs =
        begin match unphysical with
        | None -> List.map hs_of_flavors flavors
        | Some n ->  List.map (unphysical_of_flavors n) flavors
        end in
      if not (ThoList.homogeneous hs) then
        invalid_arg "Fusion.helicity_table: not all flavors have the same helicity states!"
      else
        match hs with
        | [] -> []
        | (hs_in, hs_out) :: _ ->
            sort_spin_table (Product.list2 pair (Product.list id hs_in) (Product.list id hs_out))

    module Proc = Process.Make(M)

    module WFMap = Map.Make (struct type t = F.wf let compare = compare end)
    module WFSet2 =
      Set.Make (struct type t = F.wf * (F.wf, F.coupling) Tree2.t let compare = compare end)
    module WFMap2 =
      Map.Make (struct type t = F.wf * (F.wf, F.coupling) Tree2.t let compare = compare end)
    module WFTSet =
      Set.Make (struct type t = (F.wf, F.coupling) Tree2.t let compare = compare end)

(* All wavefunctions are unique per amplitude.  So we can use per-amplitude
   dependency trees without additional \emph{internal} tags to identify identical
   wave functions. *)

(* \textbf{NB:} we miss potential optimizations, because we assume all coupling to
   be different, while in fact we have horizontal/family symmetries and non abelian
   gauge couplings are universal anyway. *)

    let disambiguate_fusions amplitudes =
      let fusions =
        ThoList.flatmap (fun amplitude ->
          List.map
            (fun fusion -> (fusion, F.dependencies amplitude (F.lhs fusion)))
            (F.fusions amplitude))
          amplitudes in
      let duplicates =
        List.fold_left
          (fun map (fusion, dependencies) ->
            let wf = F.lhs fusion in
            let set = try WFMap.find wf map with Not_found -> WFTSet.empty in
            WFMap.add wf (WFTSet.add dependencies set) map)
          WFMap.empty fusions in
      let multiplicity_map =
        WFMap.fold (fun wf dependencies acc ->
          let cardinal = WFTSet.cardinal dependencies in
          if cardinal <= 1 then
            acc
          else
            WFMap.add wf cardinal acc)
          duplicates WFMap.empty
      and dictionary_map =  
        WFMap.fold (fun wf dependencies acc ->
          let cardinal = WFTSet.cardinal dependencies in
          if cardinal <= 1 then
            acc
          else
            snd (WFTSet.fold
                   (fun dependency (i', acc') ->
                     (succ i', WFMap2.add (wf, dependency) i' acc'))
                   dependencies (1, acc)))
          duplicates WFMap2.empty in
      let multiplicity wf = 
        WFMap.find wf multiplicity_map
      and dictionary amplitude wf =
        WFMap2.find (wf, F.dependencies amplitude wf) dictionary_map in
      (multiplicity, dictionary)

    let eliminate_common_fusions1 seen_wfs amplitude =
      List.fold_left
        (fun (seen, acc) f ->
          let wf = F.lhs f in
          let dependencies = F.dependencies amplitude wf in
          if WFSet2.mem (wf, dependencies) seen then
            (seen, acc)
          else
            (WFSet2.add (wf, dependencies) seen, (f, amplitude) :: acc))
        seen_wfs (F.fusions amplitude)

    let eliminate_common_fusions processes =
      let _, rev_fusions =
        List.fold_left
          eliminate_common_fusions1
          (WFSet2.empty, []) processes in
      List.rev rev_fusions

(*i
    let eliminate_common_fusions processes =
      ThoList.flatmap
        (fun amplitude ->
          (List.map (fun f -> (f, amplitude)) (F.fusions amplitude)))
        processes
i*)

(* \thocwmodulesubsection{Calculate All The Amplitudes} *)

    let amplitudes goldstones unphysical exclusions select_wf processes =

(* \begin{dubious}
     Eventually, we might want to support inhomogeneous helicities.  However,
     this makes little physics sense for external particles on the mass shell,
     unless we have a model with degenerate massive fermions and bosons.
   \end{dubious} *)

      if not (ThoList.homogeneous (List.map hs_of_flavors processes)) then
        invalid_arg "Fusion.Multi.amplitudes: incompatible helicities";

      let unique_uncolored_processes =
        Proc.remove_duplicate_final_states (C.partition select_wf) processes in

      let progress =
        match !progress_option with
        | Quiet -> Progress.dummy
        | Channel oc -> Progress.channel oc (count_processes unique_uncolored_processes)
        | File name -> Progress.file name (count_processes unique_uncolored_processes) in

      let allowed =
        ThoList.flatmap
          (fun (fi, fo) ->
            Progress.begin_step progress (process_to_string fi fo);
            let amps = F.amplitudes goldstones exclusions select_wf fi fo in
            begin match amps with
            | [] -> Progress.end_step progress "forbidden"
            | _ -> Progress.end_step progress "allowed"
            end;
            amps) unique_uncolored_processes in
 
      Progress.summary progress "all processes done";
          
      let color_flows =
        ThoList.uniq (List.sort compare (List.map color_flow allowed))
      and flavors =
        ThoList.uniq (List.sort compare (List.map process_sans_color allowed)) in

      let vanishing_flavors =
        Proc.diff processes flavors in

      let helicities =
        helicity_table unphysical flavors in

      let f_index = 
        fst (List.fold_left
               (fun (m, i) f -> (FMap.add f i m, succ i))
               (FMap.empty, 0) flavors)
      and c_index = 
        fst (List.fold_left
               (fun (m, i) c -> (CMap.add c i m, succ i))
               (CMap.empty, 0) color_flows) in

      let table =
        Array.make_matrix (List.length flavors) (List.length color_flows) None in
      List.iter
        (fun a ->
          let f = FMap.find (process_sans_color a) f_index
          and c = CMap.find (color_flow a) c_index in
          table.(f).(c) <- Some (a))
        allowed;

      let cf_array = Array.of_list color_flows in
      let ncf = Array.length cf_array in
      let color_factor_table = Array.make_matrix ncf ncf Color.Flow.zero in

      for i = 0 to pred ncf do
        for j = 0 to i do
          color_factor_table.(i).(j) <-
            Color.Flow.factor cf_array.(i) cf_array.(j);
          color_factor_table.(j).(i) <-
            color_factor_table.(i).(j)
        done
      done;

      let fusions = eliminate_common_fusions allowed
      and multiplicity, dictionary = disambiguate_fusions allowed in
      
      { flavors = flavors;
        vanishing_flavors = vanishing_flavors;
        color_flows = color_flows;
        helicities = helicities;
        processes = allowed;
        process_table = table;
        fusions = fusions;
        multiplicity = multiplicity;
        dictionary = dictionary;
        color_factors = color_factor_table;
        constraints = C.description select_wf }

(*i
    let initialize_cache = F.initialize_cache
    let set_cache_name = F.set_cache_name
i*)

    let empty =
      { flavors = [];
        vanishing_flavors = [];
        color_flows = [];
        helicities = [];
        processes = [];
        process_table = Array.make_matrix 0 0 None;
        fusions = [];
        multiplicity = (fun _ -> 1);
        dictionary = (fun _ _ -> 1);
        color_factors = Array.make_matrix 0 0 Color.Flow.zero;
        constraints = None }

  end
