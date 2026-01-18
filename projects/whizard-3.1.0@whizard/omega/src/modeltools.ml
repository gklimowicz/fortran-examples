(* modeltools.ml --

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

(* \thocwmodulesection{Compilation} *)

(* Flavors and coupling constants:  flavors can be tested for equality
   and charge conjugation is defined.  *)

module type Flavor =
  sig
    type f
    type c
    val compare : f -> f -> int
    val conjugate : f -> f
  end

(* Compiling fusions from a list of vertices:  *)

module type Fusions =
  sig
    type t
    type f
    type c
    val fuse2 : t -> f -> f -> (f * c Coupling.t) list
    val fuse3 : t -> f -> f -> f -> (f * c Coupling.t) list
    val fuse : t -> f list -> (f * c Coupling.t) list
    val of_vertices :
        (((f * f * f) * c Coupling.vertex3 * c) list
           * ((f * f * f * f) * c Coupling.vertex4 * c) list
           * (f list * c Coupling.vertexn * c) list) -> t
  end

module Fusions (F : Flavor) : Fusions with type f = F.f and type c = F.c =
  struct

    type f = F.f
    type c = F.c

    module F2 =
      struct
        type t = f * f
        let hash = Hashtbl.hash
        let compare (f1, f2) (f1', f2') =
          let c1 = F.compare f1 f1' in
          if c1 <> 0 then
            c1
          else
            F.compare f2 f2'
        let equal f f' = compare f f' = 0
      end

    module F3 =
      struct
        type t = f * f * f
        let hash = Hashtbl.hash
        let compare (f1, f2, f3) (f1', f2', f3') =
          let c1 = F.compare f1 f1' in
          if c1 <> 0 then
            c1
          else
            let c2 = F.compare f2 f2' in
            if c2 <> 0 then
              c2
            else
              F.compare f3 f3'
        let equal f f' = compare f f' = 0
      end

    module Fn =
      struct
        type t = f list
        let hash = Hashtbl.hash
        let compare f f' = ThoList.compare ~cmp:F.compare f f'
        let equal f f' = compare f f' = 0
      end

    module H2 = Hashtbl.Make (F2)
    module H3 = Hashtbl.Make (F3)
    module Hn = Hashtbl.Make (Fn)

    type t =
        { v3 : (f * c Coupling.t) list H2.t;
          v4 : (f * c Coupling.t) list H3.t;
          vn : (f * c Coupling.t) list Hn.t }

    let lookup_fuse2 table f1 f2 =
      try H2.find table.v3 (f1, f2) with Not_found -> []

    let lookup_fuse3 table f1 f2 f3 =
      try H3.find table.v4 (f1, f2, f3) with Not_found -> []

    let lookup_fusen table f =
      try Hn.find table.vn f with Not_found -> []

    let fuse2 table f1 f2 =
      List.rev_append
        (lookup_fusen table [f1; f2])
        (lookup_fuse2 table f1 f2)

    let fuse3 table f1 f2 f3 =
      List.rev_append
        (lookup_fusen table [f1; f2; f3])
        (lookup_fuse3 table f1 f2 f3)

    let fusen table f =
      lookup_fusen table f

    let fuse table = function
      | [] | [_] -> invalid_arg "Fusions().fuse"
      | [f1; f2] -> fuse2 table f1 f2
      | [f1; f2; f3] -> fuse3 table f1 f2 f3
      | f -> fusen table f

(* Note that a pair or a triplet can appear more than once
   (e.\,g.~$e^+e^-\to \gamma$ and~$e^+e^-\to Z$).  Therefore don't
   replace the entry, but augment it instead.  *)

    let add_fusion2 table f1 f2 fusions =
      H2.add table.v3 (f1, f2) (fusions :: lookup_fuse2 table f1 f2)

    let add_fusion3 table f1 f2 f3 fusions =
      H3.add table.v4 (f1, f2, f3) (fusions :: lookup_fuse3 table f1 f2 f3)

    let add_fusionn table f fusions =
      Hn.add table.vn f (fusions :: lookup_fusen table f)

(* \begin{dubious}
     Do we need to take into account the charge conjugation
     of the coupling constants here?
   \end{dubious} *)

(* If some flavors are identical, we must not introduce the
   same vertex more than once: *)

    open Coupling

    let permute3 (f1, f2, f3) =
      [ (f1, f2), F.conjugate f3, F12;
        (f2, f1), F.conjugate f3, F21;
        (f2, f3), F.conjugate f1, F23;
        (f3, f2), F.conjugate f1, F32;
        (f3, f1), F.conjugate f2, F31;
        (f1, f3), F.conjugate f2, F13 ]

(* Here we add identical permutations of pairs only once: *)

    module F2' = Set.Make (F2)

    let add_permute3 table v c set ((f1, f2 as f12), f, p) =
      if F2'.mem f12 set then
        set
      else begin
        add_fusion2 table f1 f2 (f, V3 (v, p, c));
        F2'.add f12 set
      end

    let add_vertex3 table (f123, v, c) =
      ignore (List.fold_left (fun set f -> add_permute3 table v c set f)
                F2'.empty (permute3 f123))

(* \begin{dubious}
     Handling all the cases explicitely is OK for cubic vertices, but starts
     to become questionable already for quartic couplings.  The advantage
     remains that we can check completeness in [Targets].
   \end{dubious} *)

    let permute4 (f1, f2, f3, f4) =
      [ (f1, f2, f3), F.conjugate f4, F123;
        (f2, f3, f1), F.conjugate f4, F231;
        (f3, f1, f2), F.conjugate f4, F312;
        (f2, f1, f3), F.conjugate f4, F213;
        (f3, f2, f1), F.conjugate f4, F321;
        (f1, f3, f2), F.conjugate f4, F132;
        (f1, f2, f4), F.conjugate f3, F124;
        (f2, f4, f1), F.conjugate f3, F241;
        (f4, f1, f2), F.conjugate f3, F412;
        (f2, f1, f4), F.conjugate f3, F214;
        (f4, f2, f1), F.conjugate f3, F421;
        (f1, f4, f2), F.conjugate f3, F142;
        (f1, f3, f4), F.conjugate f2, F134;
        (f3, f4, f1), F.conjugate f2, F341;
        (f4, f1, f3), F.conjugate f2, F413;
        (f3, f1, f4), F.conjugate f2, F314;
        (f4, f3, f1), F.conjugate f2, F431;
        (f1, f4, f3), F.conjugate f2, F143;
        (f2, f3, f4), F.conjugate f1, F234;
        (f3, f4, f2), F.conjugate f1, F342;
        (f4, f2, f3), F.conjugate f1, F423;
        (f3, f2, f4), F.conjugate f1, F324;
        (f4, f3, f2), F.conjugate f1, F432;
        (f2, f4, f3), F.conjugate f1, F243 ]

(* Add identical permutations of triplets only once: *)

    module F3' = Set.Make (F3)

    let add_permute4 table v c set ((f1, f2, f3 as f123), f, p) =
      if F3'.mem f123 set then
        set
      else begin
        add_fusion3 table f1 f2 f3 (f, V4 (v, p, c));
        F3'.add f123 set
      end

    let add_vertex4 table (f1234, v, c) =
      ignore (List.fold_left (fun set f -> add_permute4 table v c set f)
                F3'.empty (permute4 f1234))

    module Fn' = Set.Make (Fn)

    let permuten = function
      | [] -> invalid_arg "Modeltools.permuten"
      | f ->
         List.map
           (fun f' ->
             match List.split f' with
             | i :: i_list, f :: f_list ->
                (f_list, F.conjugate f, i_list @ [i])
             | _ -> failwith "Modeltools.permuten: impossible")
           (Combinatorics.permute (ThoList.enumerate 1 f))

    (* This is for debugging: it provides the same permutations
       than the legacy version. *)
    let permutations = function
      | [f1; f2; f3] ->
         [ [f1; f2; f3];
           [f2; f1; f3];
           [f2; f3; f1];
           [f3; f2; f1];
           [f3; f1; f2];
           [f1; f3; f2] ]
      | [f1; f2; f3; f4] ->
         [ [f1; f2; f3; f4];
           [f1; f2; f4; f3];
           [f1; f3; f2; f4];
           [f1; f3; f4; f2];
           [f1; f4; f2; f3];
           [f1; f4; f3; f2];
           [f2; f1; f3; f4];
           [f2; f1; f4; f3];
           [f2; f3; f1; f4];
           [f2; f3; f4; f1];
           [f2; f4; f1; f3];
           [f2; f4; f3; f1];
           [f3; f1; f2; f4];
           [f3; f1; f4; f2];
           [f3; f2; f1; f4];
           [f3; f2; f4; f1];
           [f3; f4; f1; f2];
           [f3; f4; f2; f1];
           [f4; f1; f2; f3];
           [f4; f1; f3; f2];
           [f4; f2; f1; f3];
           [f4; f2; f3; f1];
           [f4; f3; f1; f2];
           [f4; f3; f2; f1] ]
      | flist -> Combinatorics.permute flist

    let permutations = Combinatorics.permute

    let permuten = function
      | [] -> invalid_arg "Modeltools.permuten"
      | f ->
         List.map
           (fun f' ->
             match List.split (List.rev f') with
             | i_list, f :: f_list ->
             (* [Printf.eprintf
                  "permuten: %s\n"
                  (ThoList.to_string string_of_int (List.rev i_list));] *)
                (List.rev f_list, F.conjugate f, List.rev i_list)
             | _ -> failwith "Modeltools.permuten: impossible")
           (permutations (ThoList.enumerate 1 f))

    let add_permuten table v c set (f12__n, f, p) =
      if Fn'.mem f12__n set then
        set
      else begin
        add_fusionn table f12__n (f, Vn (v, p, c));
        Fn'.add f12__n set
      end

    (* \begin{dubious}
         We could apply any necessary permutations
         to objects that are hidden inside of the vertex [v] here
         instead of in [Fusion.stat_fuse] and [Colorize.fuse].
       \end{dubious} *)
    let add_vertexn table (f12__n, v, c) =
      ignore
        (List.fold_left
           (fun set f -> add_permuten table v c set f)
           Fn'.empty (permuten f12__n))

    let of_vertices (vlist3, vlist4, vlistn) =
      let table =
        { v3 = H2.create 37; v4 = H3.create 37; vn = Hn.create 37 } in
      List.iter (add_vertex3 table) vlist3;
      List.iter (add_vertex4 table) vlist4;
      List.iter (add_vertexn table) vlistn;
      table

  end

module type Constant =
  sig
    type t
    val of_string : string -> t
  end

module Constant (M : Model.T) : Constant with type t = M.constant =
  struct

    type t = M.constant

    module String_Key =
      struct
        type t = string
        let hash = Hashtbl.hash
        let equal = (=)
      end
    module String_Hash = Hashtbl.Make (String_Key)

    let table = String_Hash.create 37

    let fill_table table vs =
      List.iter
        (fun (_, _, c) ->
          String_Hash.add table (M.constant_symbol c) c)
        vs

    (* Delay loading of the tables until the first use, so that
       [M.vertices] can be initialized from a file.  *)

    let tables_filled = ref false

    let fill_tables () =
      if not !tables_filled then begin
	let (v3, v4, vn) = M.vertices () in
	fill_table table v3;
	fill_table table v4;
	fill_table table vn;
	tables_filled := true
      end

    let of_string name =
      try
	fill_tables ();
        String_Hash.find table name
      with
      | Not_found ->
          invalid_arg
            ("Constant(Model).of_string: unknown coupling constant: " ^ name)

  end

(* \thocwmodulesection{Mutable Models} *)

module Mutable (FGC : sig type f and g and c end) : Model.Mutable
       with type flavor = FGC.f and type gauge = FGC.g and type constant = FGC.c =
  struct
    type flavor = FGC.f
    type gauge = FGC.g
    type constant = FGC.c

    let init () = ()

    let options = Options.empty
    let caveats () = []

    module Ch = Charges.Null
    let charges _ = ()

    exception Uninitialized of string
    let uninitialized name =
      raise (Uninitialized name)
      
(* Note that [lookup] works, by the magic of currying, for any arity.  But
   we need to supply one argument to delay evaluation. *)

(* Also note that the references are \emph{not} shared among results
   of functor applications.  Simple module renaming causes sharing.  *)
    let declare template =
      let reference = ref template in
      let update fct = reference := fct
      and lookup arg = !reference arg in
      (update, lookup)

    let set_color, color =
      declare (fun f -> uninitialized "color")

    let set_nc, nc =
      declare (fun f -> uninitialized "nc")

    let set_pdg, pdg =
      declare (fun f -> uninitialized "pdg")

    let set_lorentz, lorentz =
      declare (fun f -> uninitialized "lorentz")

    let set_propagator, propagator =
      declare (fun f -> uninitialized "propagator")

    let set_width, width =
      declare (fun f -> uninitialized "width")

    let set_goldstone, goldstone =
      declare (fun f -> uninitialized "goldstone")

    let set_conjugate, conjugate =
      declare (fun f -> uninitialized "conjugate")

    let set_fermion, fermion =
      declare (fun f -> uninitialized "fermion")

    let set_max_degree, max_degree =
      declare (fun () -> uninitialized "max_degree")

    let set_vertices, vertices =
      declare (fun () -> uninitialized "vertices")

    let set_fuse2, fuse2 =
      declare (fun f1 f2 -> uninitialized "fuse2")

    let set_fuse3, fuse3 =
      declare (fun f1 f2 f3 -> uninitialized "fuse3")

    let set_fuse, fuse =
      declare (fun f -> uninitialized "fuse")

    let set_flavors, flavors =
      declare (fun () -> [])

    let set_external_flavors, external_flavors =
      declare (fun () -> [("uninitialized", [])])

    let set_parameters, parameters =
      declare (fun () -> uninitialized "parameters")

    let set_flavor_of_string, flavor_of_string =
      declare (fun f -> uninitialized "flavor_of_string")

    let set_flavor_to_string, flavor_to_string =
      declare (fun f -> uninitialized "flavor_to_string")

    let set_flavor_to_TeX, flavor_to_TeX =
      declare (fun f -> uninitialized "flavor_to_TeX")

    let set_flavor_symbol, flavor_symbol =
      declare (fun f -> uninitialized "flavor_symbol")

    let set_gauge_symbol, gauge_symbol =
      declare (fun g -> uninitialized "gauge_symbol")

    let set_mass_symbol, mass_symbol =
      declare (fun f -> uninitialized "mass_symbol")

    let set_width_symbol, width_symbol =
      declare (fun f -> uninitialized "width_symbol")

    let set_constant_symbol, constant_symbol =
      declare (fun c -> uninitialized "constant_symbol")

    module F = Fusions (struct
      type f = flavor
      type c = constant
      let compare = compare
      let conjugate = conjugate
    end)

    let max_degree_of_vertices (v3, v4, vn) =
      List.fold_left
        (fun acc (p, _, _) -> max acc (List.length p))
        (max (match v3 with [] -> 0 | _ -> 3) (match v4 with [] -> 0 | _ -> 4))
        vn

    let setup ~color ~nc ~pdg ~lorentz ~propagator ~width ~goldstone
        ~conjugate ~fermion ~vertices
        ~flavors ~parameters ~flavor_of_string ~flavor_to_string
        ~flavor_to_TeX ~flavor_symbol
        ~gauge_symbol ~mass_symbol ~width_symbol ~constant_symbol =
      set_color color;
      set_nc nc;
      set_pdg pdg;
      set_lorentz lorentz;
      set_propagator propagator;
      set_width width;
      set_goldstone goldstone;
      set_conjugate conjugate;
      set_fermion fermion;
      let v = vertices () in
      let max_degree = max_degree_of_vertices v in
      set_max_degree (fun () -> max_degree);
      set_vertices (fun () -> v);
      let table = F.of_vertices v in
      set_fuse2 (F.fuse2 table);
      set_fuse3 (F.fuse3 table);
      set_fuse (F.fuse table);
      set_external_flavors (fun () -> flavors);
      let flavors = ThoList.flatmap snd flavors in
      set_flavors (fun () -> flavors);
      set_parameters parameters;
      set_flavor_of_string flavor_of_string;
      set_flavor_to_string flavor_to_string;
      set_flavor_to_TeX flavor_to_TeX;
      set_flavor_symbol flavor_symbol;
      set_gauge_symbol gauge_symbol;
      set_mass_symbol mass_symbol;
      set_width_symbol width_symbol;
      set_constant_symbol constant_symbol

  end

module Static (M : Model.T) =
  struct
    type flavor = M.flavor
    type gauge = M.gauge
    type constant = M.constant
    module Ch = M.Ch
    let color = M.color
    let nc = M.nc
    let charges = M.charges
    let pdg = M.pdg
    let lorentz = M.lorentz
    let propagator = M.propagator
    let width = M.width
    let conjugate = M.conjugate
    let fermion = M.fermion
    let max_degree = M.max_degree
    let vertices = M.vertices
    let fuse2 = M.fuse2
    let fuse3 = M.fuse3
    let fuse = M.fuse
    let flavors = M.flavors
    let external_flavors = M.external_flavors
    let goldstone = M.goldstone
    let parameters = M.parameters
    let flavor_of_string = M.flavor_of_string
    let flavor_to_string = M.flavor_to_string
    let flavor_to_TeX = M.flavor_to_TeX
    let flavor_symbol = M.flavor_symbol
    let gauge_symbol = M.gauge_symbol
    let mass_symbol = M.mass_symbol
    let width_symbol = M.width_symbol
    let constant_symbol = M.constant_symbol
    let options = M.options
    let caveats = M.caveats
    let init () = ()
    let setup ~color ~nc ~pdg ~lorentz ~propagator ~width ~goldstone
        ~conjugate ~fermion ~vertices
        ~flavors ~parameters ~flavor_of_string ~flavor_to_string
        ~flavor_to_TeX ~flavor_symbol
        ~gauge_symbol ~mass_symbol ~width_symbol ~constant_symbol =
      ()
  end

(* \thocwmodulesection{Topology Only} *)

(* UFO models can have more than one Lorentz structure for a
   given flavor combination.  This messes up the phase space
   generation.  There we need to be able to ignore the redundant
   flavor combinations. *)

(* Filter vertices with more than one Lorentz structure
   for a combination of flavors.  Only the first Lorentz
   structure is kept. *)
let filter_couplings flavor_coupling_list =
  List.map
    (fun (f, c_list) -> (f, List.hd c_list))
    (ThoList.factorize flavor_coupling_list)

let triple_to_nested (a, b, c) = (a, (b, c))

let nested_to_triple (a, (b, c)) = (a, b, c)

let filter_couplings_triples fc =
  List.map
    nested_to_triple
    (filter_couplings (List.map triple_to_nested fc))

(* \begin{dubious}
     It would be clearer to replace [constant Coupling.t] by
     [unit] in the resultig model, but that would require
     much more code duplication.
   \end{dubious} *)

module Topology (M : Model.T) =
  struct
    type flavor = M.flavor
    type gauge = M.gauge
    type constant = M.constant
    module Ch = M.Ch
    let color = M.color
    let nc = M.nc
    let charges = M.charges
    let pdg = M.pdg
    let lorentz = M.lorentz
    let propagator = M.propagator
    let width = M.width
    let conjugate = M.conjugate
    let fermion = M.fermion
    let max_degree = M.max_degree
    let vertices () =
      let (v3, v4, vn) = M.vertices () in
      (filter_couplings_triples v3,
       filter_couplings_triples v4,
       filter_couplings_triples vn)
    let fuse2 f1 f2 = filter_couplings (M.fuse2 f1 f2)
    let fuse3 f1 f2 f3 = filter_couplings (M.fuse3 f1 f2 f3)
    let fuse f_list = filter_couplings (M.fuse f_list)
    let flavors = M.flavors
    let external_flavors = M.external_flavors
    let goldstone = M.goldstone
    let parameters = M.parameters
    let flavor_of_string = M.flavor_of_string
    let flavor_to_string = M.flavor_to_string
    let flavor_to_TeX = M.flavor_to_TeX
    let flavor_symbol = M.flavor_symbol
    let gauge_symbol = M.gauge_symbol
    let mass_symbol = M.mass_symbol
    let width_symbol = M.width_symbol
    let constant_symbol = M.constant_symbol
    let options = M.options
    let caveats = M.caveats
  end

module Topology3 (M : Model.T) =
  struct
    type flavor = M.flavor
    type gauge = M.gauge
    type constant = M.constant
    module Ch = M.Ch
    let color = M.color
    let nc = M.nc
    let charges = M.charges
    let pdg = M.pdg
    let lorentz = M.lorentz
    let propagator = M.propagator
    let width = M.width
    let conjugate = M.conjugate
    let fermion = M.fermion
    let max_degree = M.max_degree
    let vertices () =
      let (v3, _, vn) = M.vertices () in
      (filter_couplings_triples v3,
       [],
       filter_couplings_triples
         (List.filter (fun (f, _, _) -> List.length f < 3) vn))
    let fuse2 f1 f2 = filter_couplings (M.fuse2 f1 f2)
    let fuse3 f1 f2 f3 = []
    let fuse = function
      | [_; _] as f_list -> filter_couplings (M.fuse f_list)
      | _ -> []
    let flavors = M.flavors
    let external_flavors = M.external_flavors
    let goldstone = M.goldstone
    let parameters = M.parameters
    let flavor_of_string = M.flavor_of_string
    let flavor_to_string = M.flavor_to_string
    let flavor_to_TeX = M.flavor_to_TeX
    let flavor_symbol = M.flavor_symbol
    let gauge_symbol = M.gauge_symbol
    let mass_symbol = M.mass_symbol
    let width_symbol = M.width_symbol
    let constant_symbol = M.constant_symbol
    let options = M.options
    let caveats = M.caveats
  end
