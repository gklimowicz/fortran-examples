(* cascade.ml --

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

    type constant
    type flavor
    type p

    type t
    val of_string_list : int -> string list -> t
    val to_string : t -> string

    type selectors
    val to_selectors : t -> selectors
    val no_cascades : selectors

    val select_wf : selectors -> (p -> bool) -> flavor -> p -> p list -> bool
    val select_p : selectors -> p -> p list -> bool
    val on_shell : selectors -> flavor -> p -> bool
    val is_gauss : selectors -> flavor -> p -> bool

    val select_vtx : selectors -> constant Coupling.t ->
      flavor -> flavor list -> bool

    val partition : selectors -> int list list
    val description : selectors -> string option

  end

module Make (M : Model.T) (P : Momentum.T) :
    (T with type flavor = M.flavor and type constant = M.constant and type p = P.t) =
  struct

    module CS = Cascade_syntax

    type constant = M.constant
    type flavor = M.flavor
    type p = P.t

(* Since we have
   \begin{equation}
      p \le q \Longleftrightarrow (-q) \le (-p)
   \end{equation}
   also for $\le$ as set inclusion [lesseq], only four of the eight
   combinations are independent
   \begin{equation}
     \begin{aligned}
        p &\le q    &&\Longleftrightarrow & (-q) &\le (-p) \\
        q &\le p    &&\Longleftrightarrow & (-p) &\le (-q) \\
        p &\le (-q) &&\Longleftrightarrow & q    &\le (-p) \\
     (-q) &\le p    &&\Longleftrightarrow & (-p) &\le q 
     \end{aligned}
   \end{equation}  *)

    let one_compatible p q =
      let neg_q = P.neg q in
      P.lesseq p q ||
      P.lesseq q p ||
      P.lesseq p neg_q ||
      P.lesseq neg_q p

(* 'tis wasteful \ldots (at least by a factor of two, because every momentum
    combination is generated, including the negative ones. *)

    let all_compatible p p_list q =
      let l = List.length p_list in
      if l <= 2 then
        one_compatible p q
      else
        let tuple_lengths = ThoList.range 2 (succ l / 2) in
        let tuples = ThoList.flatmap (fun n -> Combinatorics.choose n p_list) tuple_lengths in
        let momenta = List.map (List.fold_left P.add (P.zero (P.dim q))) tuples in
        List.for_all (one_compatible q) momenta

(* The following assumes that the [flavor list] is always very short.  Otherwise
   one should use an efficient set implementation. *)

    type wf =
      | True
      | False
      | On_shell of flavor list * P.t
      | On_shell_not of flavor list * P.t
      | Off_shell of flavor list * P.t
      | Off_shell_not of flavor list * P.t
      | Gauss of flavor list * P.t
      | Gauss_not of flavor list * P.t
      | Any_flavor of P.t
      | And of wf list

    module Constant = Modeltools.Constant (M)

    type vtx =
        { couplings : M.constant list;
          fields : flavor list }

    type t =
        { wf : wf;
          (* TODO: The following lists should be sets for efficiency. *)
          flavors : flavor list;
          vertices : vtx list }

    let default =
      { wf = True;
        flavors = [];
        vertices = [] }

    let of_string s = 
      Cascade_parser.main Cascade_lexer.token (Lexing.from_string s)

(* \begin{dubious}
     If we knew that we're dealing with a scattering, we could apply
     [P.flip_s_channel_in] to all momenta, so that $1+2$ accepts the particle
     and not the antiparticle.  Right now, we don't have this information.
   \end{dubious} *)

    let only_wf wf = { default with wf = wf }

    let cons_and_wf c wfs =
      match c.wf, wfs with
      | True, wfs -> wfs
      | False, _ -> [False]
      | wf, [] -> [wf]
      | wf, wfs -> wf :: wfs

    let and_cascades_wf c =
      match List.fold_right cons_and_wf c [] with
      | [] -> True
      | [wf] -> wf
      | wfs -> And wfs

    let uniq l =
      ThoList.uniq (List.sort compare l)

    let import dim cascades =
      let rec import' = function
        | CS.True ->
            only_wf True
        | CS.False ->
            only_wf False
        | CS.On_shell (f, p) ->
            only_wf
              (On_shell (List.map M.flavor_of_string f, P.of_ints dim p))
        | CS.On_shell_not (f, p) ->
            only_wf
              (On_shell_not (List.map M.flavor_of_string f, P.of_ints dim p))
        | CS.Off_shell (fs, p) ->
            only_wf
              (Off_shell (List.map M.flavor_of_string fs, P.of_ints dim p))
        | CS.Off_shell_not (fs, p) ->
            only_wf
              (Off_shell_not (List.map M.flavor_of_string fs, P.of_ints dim p))
        | CS.Gauss (f, p) ->
            only_wf
              (Gauss (List.map M.flavor_of_string f, P.of_ints dim p))
        | CS.Gauss_not (f, p) ->
            only_wf
              (Gauss (List.map M.flavor_of_string f, P.of_ints dim p))
        | CS.Any_flavor p ->
            only_wf (Any_flavor (P.of_ints dim p))
        | CS.And cs ->
            let cs = List.map import' cs in
            { wf = and_cascades_wf cs;
              flavors = uniq (List.concat
                                (List.map (fun c -> c.flavors) cs));
              vertices = uniq (List.concat
                                 (List.map (fun c -> c.vertices) cs)) }
        | CS.X_Flavor fs ->
            let fs = List.map M.flavor_of_string fs in
            { default with flavors = uniq (fs @ List.map M.conjugate fs) }
        | CS.X_Vertex (cs, fss) ->
            let cs = List.map Constant.of_string cs
            and fss = List.map (List.map M.flavor_of_string) fss in
            let expanded =
              List.map
                (fun fs -> { couplings = cs; fields = fs })
                (match fss with
                | [] -> [[]] (* Subtle: \emph{not} an empty list! *)
                | fss -> Product.list (fun fs -> fs) fss) in
            { default with vertices = expanded }
      in
      import' cascades

    let of_string_list dim strings =
      match List.map of_string strings with
      | [] -> default
      | first :: next ->
          import dim (List.fold_right CS.mk_and next first)

    let flavors_to_string fs =
      (String.concat ":" (List.map M.flavor_to_string fs))

    let momentum_to_string p =
      String.concat "+" (List.map string_of_int (P.to_ints p))

    let rec wf_to_string = function
      | True ->
          "true"
      | False ->
          "false"
      | On_shell (fs, p) ->
          momentum_to_string p ^ " = " ^ flavors_to_string fs
      | On_shell_not (fs, p) ->
          momentum_to_string p ^ " = !" ^ flavors_to_string fs
      | Off_shell (fs, p) ->
          momentum_to_string p  ^ " ~ " ^ flavors_to_string fs
      | Off_shell_not (fs, p) ->
          momentum_to_string p  ^ " ~ !" ^ flavors_to_string fs
      | Gauss (fs, p) ->
          momentum_to_string p ^ " # " ^ flavors_to_string fs
      | Gauss_not (fs, p) ->
          momentum_to_string p ^ " # !" ^ flavors_to_string fs
      | Any_flavor p ->
          momentum_to_string p ^ " ~ ?"
      | And cs ->
          String.concat " && " (List.map (fun c -> "(" ^ wf_to_string c ^ ")") cs)

    let vertex_to_string v =
      "^" ^ String.concat ":" (List.map M.constant_symbol v.couplings) ^
      "[" ^ String.concat "," (List.map M.flavor_to_string v.fields) ^ "]"

    let vertices_to_string vs =
      (String.concat " && " (List.map vertex_to_string vs))

    let to_string = function
      | { wf = True; flavors = []; vertices = [] } ->
          ""
      | { wf = True; flavors = fs; vertices = [] } ->
          "!" ^ flavors_to_string fs
      | { wf = True; flavors = []; vertices = vs } ->
          vertices_to_string vs
      | { wf = True; flavors = fs; vertices = vs } ->
          "!" ^ flavors_to_string fs ^ " && " ^ vertices_to_string vs
      | { wf = wf; flavors = []; vertices = [] } ->
          wf_to_string wf
      | { wf = wf; flavors = []; vertices = vs } ->
          vertices_to_string vs ^ " && " ^ wf_to_string wf
      | { wf = wf; flavors = fs; vertices = [] } ->
          "!" ^ flavors_to_string fs ^ " && " ^ wf_to_string wf
      | { wf = wf; flavors = fs; vertices = vs } ->
          "!" ^ flavors_to_string fs ^
          " && " ^ vertices_to_string vs ^
          " && " ^ wf_to_string wf

    type selectors =
        { select_p : p -> p list -> bool;
          select_wf : (p -> bool) -> flavor -> p -> p list -> bool;
          on_shell : flavor -> p -> bool;
          is_gauss : flavor -> p -> bool;
          select_vtx : constant Coupling.t -> flavor -> flavor list -> bool;
          partition : int list list;
          description : string option }

    let no_cascades =
      { select_p = (fun _ _ -> true);
        select_wf = (fun _ _ _ _ -> true);
        on_shell = (fun _ _ -> false);
        is_gauss = (fun _ _ -> false);
        select_vtx = (fun _ _ _ -> true);
        partition = [];
        description = None }

    let select_p s = s.select_p
    let select_wf s = s.select_wf
    let on_shell s = s.on_shell
    let is_gauss s = s.is_gauss
    let select_vtx s = s.select_vtx
    let partition s = s.partition
    let description s = s.description

    let to_select_p cascades p p_in =
      let rec to_select_p' = function
        | True -> true
        | False -> false
        | On_shell (_, momentum) | On_shell_not (_, momentum)
        | Off_shell (_, momentum) | Off_shell_not (_, momentum)
        | Gauss (_, momentum) | Gauss_not (_, momentum)
        | Any_flavor momentum -> all_compatible p p_in momentum
        | And [] -> false
        | And cs -> List.for_all to_select_p' cs in
      to_select_p' cascades

    let to_select_wf cascades is_timelike f p p_in =
      let f' = M.conjugate f in
      let rec to_select_wf' = function
        | True -> true
        | False -> false
        | Off_shell (flavors, momentum) ->
            if p = momentum then
              List.mem f' flavors || (if is_timelike p then false else List.mem f flavors)
            else if p = P.neg momentum then
              List.mem f flavors || (if is_timelike p then false else List.mem f' flavors)
            else
              one_compatible p momentum && all_compatible p p_in momentum
        | On_shell (flavors, momentum) | Gauss (flavors, momentum) ->
             if is_timelike p then begin
               if p = momentum then
                 List.mem f' flavors
               else if p = P.neg momentum then
                 List.mem f flavors
               else
                 one_compatible p momentum && all_compatible p p_in momentum
             end else
               false
        | Off_shell_not (flavors, momentum) ->
            if p = momentum then
              not (List.mem f' flavors || (if is_timelike p then false else List.mem f flavors))
            else if p = P.neg momentum then
              not (List.mem f flavors || (if is_timelike p then false else List.mem f' flavors))
            else
              one_compatible p momentum && all_compatible p p_in momentum
        | On_shell_not (flavors, momentum) | Gauss_not (flavors, momentum) ->
            if is_timelike p then begin
              if p = momentum then
                not (List.mem f' flavors)
              else if p = P.neg momentum then
                not (List.mem f flavors)
              else
                one_compatible p momentum && all_compatible p p_in momentum
            end else
              false
        | Any_flavor momentum ->
            one_compatible p momentum && all_compatible p p_in momentum
        | And [] -> false
        | And cs -> List.for_all to_select_wf' cs in
      not (List.mem f cascades.flavors) && to_select_wf' cascades.wf


(* In case you're wondering: [to_on_shell f p] and [is_gauss f p] only search
   for on shell conditions and are to be used in a target, not in [Fusion]! *)

    let to_on_shell cascades f p =
      let f' = M.conjugate f in
      let rec to_on_shell' = function
        | True | False | Any_flavor _
        | Off_shell (_, _) | Off_shell_not (_, _)
        | Gauss (_, _) | Gauss_not (_, _) -> false
        | On_shell (flavors, momentum) ->
            (p = momentum || p = P.neg momentum) && (List.mem f flavors || List.mem f' flavors)
        | On_shell_not (flavors, momentum) ->
            (p = momentum || p = P.neg momentum) && not (List.mem f flavors || List.mem f' flavors)
        | And [] -> false
        | And cs -> List.for_all to_on_shell' cs in
      to_on_shell' cascades


    let to_gauss cascades f p =
      let f' = M.conjugate f in
      let rec to_gauss' = function
        | True | False | Any_flavor _
        | Off_shell (_, _) | Off_shell_not (_, _)
        | On_shell (_, _) | On_shell_not (_, _) -> false
        | Gauss (flavors, momentum) ->
            (p = momentum || p = P.neg momentum) &&
            (List.mem f flavors || List.mem f' flavors)
        | Gauss_not (flavors, momentum) ->
            (p = momentum || p = P.neg momentum) &&
            not (List.mem f flavors || List.mem f' flavors)
        | And [] -> false
        | And cs -> List.for_all to_gauss' cs in
      to_gauss' cascades

    module Fields =
      struct
        type f = M.flavor
        type c = M.constant list
        let compare = compare
        let conjugate = M.conjugate
      end

    module Fusions = Modeltools.Fusions (Fields)

    let dummy3 = Coupling.Scalar_Scalar_Scalar 1
    let dummy4 = Coupling.Scalar4 1
    let dummyn = Coupling.UFO (Algebra.QC.unit, "dummy", [], [], Color.Vertex.one)

(* Translate the vertices in a pair of lists: the first is the list
   of always rejected couplings and the second the remaining
   vertices suitable as input to [Fusions.of_vertices]. *)

    let translate_vertices vertices =
      List.fold_left
        (fun (cs, (v3, v4, vn) as acc) v ->
          match v.fields with
          | [] -> (v.couplings @ cs, (v3, v4, vn))
          | [_] | [_;_] -> acc
          | [f1; f2; f3] ->
              (cs, (((f1, f2, f3), dummy3, v.couplings)::v3, v4, vn))
          | [f1; f2; f3; f4] ->
              (cs, (v3, ((f1, f2, f3, f4), dummy4, v.couplings)::v4, vn))
          | fs -> (cs, (v3, v4, (fs, dummyn, v.couplings)::vn)))
        ([], ([], [], [])) vertices

(*i
    let fusion_to_string c f fs =
      M.flavor_to_string f ^ " <- " ^ M.constant_symbol c ^ "[" ^
      String.concat " , " (List.map M.flavor_to_string fs) ^ "]"
i*)

    let unpack_constant = function
      | Coupling.V3 (_, _, cs) -> cs
      | Coupling.V4 (_, _, cs) -> cs
      | Coupling.Vn (_, _, cs) -> cs

(* Sometimes, the empty list is a wildcard and matches any coupling: *)

    let match_coupling c cs =
      List.mem c cs

    let match_coupling_wildcard c = function
      | [] -> true
      | cs -> match_coupling c cs

    let to_select_vtx cascades =
      match cascades.vertices with
      | [] ->
          (* No vertex constraints means that we always accept. *)
          (fun c f fs -> true)
      | vertices ->
          match translate_vertices vertices with
          | [], ([],[],[]) ->
              (* If [cascades.vertices] is not empty, we mustn't
                 get here \ldots *)
              failwith "Cascade.to_select_vtx: unexpected"
          | couplings, ([],[],[]) ->
              (* No constraints on the fields.  Just make sure that the
		 coupling [c] doesn't appear in the vetoed [couplings]. *)
              (fun c f fs ->
                let c = unpack_constant c in
                not (match_coupling c couplings))
          | couplings, vertices ->
              (* Make sure that [Fusions.of_vertices] is only evaluated
		 once for efficiency. *)
              let fusions = Fusions.of_vertices vertices in
              (fun c f fs ->
                let c = unpack_constant c in
		(* Make sure that none of the vetoed [couplings] matches.
		   Here an empty [couplings] list is \emph{not} a
		   wildcard. *)
                if match_coupling c couplings then
		  false
		else
		  (* Also make sure that none of the vetoed [vertices]
		     matches.  Here an empty [couplings] list \emph{is}
		     a wildcard. *)
                  not (List.exists
                         (fun (f', cs') ->
                           let cs' = unpack_constant cs' in
                           f = f' && match_coupling_wildcard c cs')
                         (Fusions.fuse fusions fs)))
        
(* \begin{dubious}
     Not a working implementation yet, but it isn't used either \ldots 
   \end{dubious} *)

    module IPowSet =
      PowSet.Make (struct type t = int let compare = compare let to_string = string_of_int end)
      
    let rec coarsest_partition' = function
        | True | False -> IPowSet.empty
        | On_shell (_, momentum) | On_shell_not (_, momentum)
        | Off_shell (_, momentum) | Off_shell_not (_, momentum)
        | Gauss (_, momentum) | Gauss_not (_, momentum)
        | Any_flavor momentum -> IPowSet.of_lists [P.to_ints momentum]
        | And [] -> IPowSet.empty
        | And cs -> IPowSet.basis (IPowSet.union (List.map coarsest_partition' cs))

    let coarsest_partition cascades =
      let p = coarsest_partition' cascades in
      if IPowSet.is_empty p then
        []
      else
        IPowSet.to_lists p

    let part_to_string part =
      "{" ^ String.concat "," (List.map string_of_int part) ^ "}"

    let partition_to_string = function
      | [] -> ""
      | parts ->
          "  grouping {" ^ String.concat "," (List.map part_to_string parts) ^ "}"

    let to_selectors = function
      | { wf = True; flavors = []; vertices = [] } -> no_cascades
      | c ->
          let partition = coarsest_partition c.wf in
          { select_p = to_select_p c.wf;
            select_wf = to_select_wf c;
            on_shell = to_on_shell c.wf;
            is_gauss = to_gauss c.wf;
            select_vtx = to_select_vtx c;
            partition = partition;
            description = Some (to_string c ^ partition_to_string partition) }

(*i
    let to_selectors cascades =
      prerr_endline (">>> " ^ to_string cascades);
      to_selectors cascades
i*)
  end
