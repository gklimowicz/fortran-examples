(* UFO.ml --

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

(* Unfortunately, \texttt{ocamlweb} will not typeset all multi character
   operators nicely. E.\,g.~\verb+f @< g+ comes out as [f @< g]. *)

let (<*>) f g x =
 f (g x)

let (<**>) f g x y =
  f (g x y)

module SMap = Map.Make (struct type t = string let compare = compare end)
module SSet = Sets.String

module CMap =
  Map.Make
    (struct
      type t = string
      let compare = ThoString.compare_caseless
    end)
module CSet = Sets.String_Caseless

let error_in_string text start_pos end_pos =
  let i = start_pos.Lexing.pos_cnum
  and j = end_pos.Lexing.pos_cnum in
  String.sub text i (j - i)

let error_in_file name start_pos end_pos =
  Printf.sprintf
    "%s:%d.%d-%d.%d"
    name
    start_pos.Lexing.pos_lnum
    (start_pos.Lexing.pos_cnum - start_pos.Lexing.pos_bol)
    end_pos.Lexing.pos_lnum
    (end_pos.Lexing.pos_cnum - end_pos.Lexing.pos_bol)

let parse_string text =
  try
    UFO_parser.file
      UFO_lexer.token
      (UFO_lexer.init_position "" (Lexing.from_string text))
  with
  | UFO_tools.Lexical_Error (msg, start_pos, end_pos) ->
     invalid_arg (Printf.sprintf "lexical error (%s) at: `%s'"
                    msg  (error_in_string text start_pos end_pos))
  | UFO_syntax.Syntax_Error (msg, start_pos, end_pos) ->
     invalid_arg (Printf.sprintf "syntax error (%s) at: `%s'"
                    msg  (error_in_string text start_pos end_pos))
  | Parsing.Parse_error ->
     invalid_arg ("parse error: " ^ text)

exception File_missing of string

let parse_file name =
  let ic =
    try open_in name with
    | Sys_error msg as exc ->
       if msg = name ^ ": No such file or directory" then
         raise (File_missing name)
       else
         raise exc in
  let result =
    begin
      try
	UFO_parser.file
	  UFO_lexer.token
	  (UFO_lexer.init_position name (Lexing.from_channel ic))
      with
      | UFO_tools.Lexical_Error (msg, start_pos, end_pos) ->
	 begin
	   close_in ic;
	   invalid_arg (Printf.sprintf
			  "%s: lexical error (%s)"
			  (error_in_file name start_pos end_pos) msg)
	 end
      | UFO_syntax.Syntax_Error (msg, start_pos, end_pos) ->
	 begin
	   close_in ic;
	   invalid_arg (Printf.sprintf
			  "%s: syntax error (%s)"
			  (error_in_file name start_pos end_pos) msg)
	 end
      | Parsing.Parse_error ->
	 begin
	   close_in ic;
	   invalid_arg ("parse error: " ^ name)
	 end
    end in
  close_in ic;
  result

(* These are the contents of the Python files after lexical
   analysis as context-free variable declarations, before
   any semantic interpretation. *)

module type Files =
  sig
    
    type t = private
      { particles : UFO_syntax.t;
	couplings : UFO_syntax.t;
	coupling_orders : UFO_syntax.t;
	vertices : UFO_syntax.t;
	lorentz : UFO_syntax.t;
	parameters : UFO_syntax.t;
	propagators : UFO_syntax.t;
	decays : UFO_syntax.t }

    val parse_directory : string -> t

  end

module Files : Files =
  struct
    
    type t =
      { particles : UFO_syntax.t;
	couplings : UFO_syntax.t;
	coupling_orders : UFO_syntax.t;
	vertices : UFO_syntax.t;
	lorentz : UFO_syntax.t;
	parameters : UFO_syntax.t;
	propagators : UFO_syntax.t;
	decays : UFO_syntax.t }

    let parse_directory dir =
      let filename stem = Filename.concat dir (stem ^ ".py") in
      let parse stem = parse_file (filename stem) in
      let parse_optional stem =
        try parse stem with File_missing _ -> [] in
      { particles = parse "particles";
	couplings = parse "couplings";
	coupling_orders = parse_optional "coupling_orders";
	vertices = parse "vertices";
	lorentz = parse "lorentz";
	parameters = parse "parameters";
	propagators = parse_optional "propagators";
	decays = parse_optional "decays" }

  end

let dump_file pfx f =
  List.iter
    (fun s -> print_endline (pfx ^ ": " ^ s))
    (UFO_syntax.to_strings f)

type charge =
  | Q_Integer of int
  | Q_Fraction of int * int

let charge_to_string = function
  | Q_Integer i -> Printf.sprintf "%d" i
  | Q_Fraction (n, d) -> Printf.sprintf "%d/%d" n d

module S = UFO_syntax

let find_attrib name attribs =
  try
    (List.find (fun a -> name = a.S.a_name) attribs).S.a_value
  with
  | Not_found -> failwith ("UFO.find_attrib: \"" ^ name ^ "\" not found")

let find_attrib name attribs =
  (List.find (fun a -> name = a.S.a_name) attribs).S.a_value

let name_to_string ?strip name =
  let stripped =
    begin match strip, List.rev name with
    | Some pfx, head :: tail ->
       if pfx = head then
	 tail
       else
	 failwith ("UFO.name_to_string: expected prefix '" ^ pfx ^
		      "', got '" ^ head ^ "'")
    | _, name -> name
    end in
  String.concat "." stripped

let name_attrib ?strip name attribs =
  match find_attrib name attribs with
  | S.Name n -> name_to_string ?strip n
  | _ -> invalid_arg ("UFO.name_attrib: " ^ name)

let integer_attrib name attribs =
  match find_attrib name attribs with
  | S.Integer i -> i
  | _ -> invalid_arg ("UFO.integer_attrib: " ^ name)

let charge_attrib name attribs =
  match find_attrib name attribs with
  | S.Integer i -> Q_Integer i
  | S.Fraction (n, d) -> Q_Fraction (n, d)
  | _ -> invalid_arg ("UFO.charge_attrib: " ^ name)

let string_attrib name attribs =
  match find_attrib name attribs with
  | S.String s -> s
  | _ -> invalid_arg ("UFO.string_attrib: " ^ name)

let string_expr_attrib name attribs =
  match find_attrib name attribs with
  | S.Name n -> [S.Macro n]
  | S.String s -> [S.Literal s]
  | S.String_Expr e -> e
  | _ -> invalid_arg ("UFO.string_expr_attrib: " ^ name)

let boolean_attrib name attribs =
  try
    match ThoString.lowercase (name_attrib name attribs) with
    | "true" -> true
    | "false" -> false
    | _ -> invalid_arg ("UFO.boolean_attrib: " ^ name)
  with
  | Not_found -> false

type value =
  | Integer of int
  | Fraction of int * int
  | Float of float
  | Expr of UFOx.Expr.t
  | Name of string list

let map_expr f default = function
  | Integer _ | Fraction (_, _) | Float _ | Name _ -> default
  | Expr e -> f e

let variables = map_expr UFOx.Expr.variables CSet.empty
let functions = map_expr UFOx.Expr.functions CSet.empty

let add_to_set_in_map key element map =
  let set = try CMap.find key map with Not_found -> CSet.empty in
  CMap.add key (CSet.add element set) map

(* Add all variables in [value] to the [map] from variables
   to the names in which they appear, indicating
   that [name] depends on these variables. *)
let dependency name value map =
  CSet.fold
    (fun variable acc -> add_to_set_in_map variable name acc)
    (variables value)
    map

let dependencies name_value_list =
  List.fold_left
    (fun acc (name, value) -> dependency name value acc)
    CMap.empty
    name_value_list

let dependency_to_string (variable, appearences) =
  Printf.sprintf
    "%s -> {%s}"
    variable (String.concat ", " (CSet.elements appearences))

let dependencies_to_strings map =
  List.map dependency_to_string (CMap.bindings map)

let expr_to_string =
  UFOx.Value.to_string <*> UFOx.Value.of_expr

let value_to_string = function
  | Integer i -> Printf.sprintf "%d" i
  | Fraction (n, d) -> Printf.sprintf "%d/%d" n d
  | Float x -> string_of_float x
  | Expr e -> "'" ^ expr_to_string e ^ "'"
  | Name n -> name_to_string n

let value_to_expr substitutions = function
  | Integer i -> Printf.sprintf "%d" i
  | Fraction (n, d) -> Printf.sprintf "%d/%d" n d
  | Float x -> string_of_float x
  | Expr e -> expr_to_string (substitutions e)
  | Name n -> name_to_string n

let value_to_coupling substitutions atom = function
  | Integer i -> Coupling.Integer i
  | Fraction (n, d) -> Coupling.Quot (Coupling.Integer n, Coupling.Integer d)
  | Float x -> Coupling.Float x
  | Expr e ->
     UFOx.Value.to_coupling atom (UFOx.Value.of_expr (substitutions e))
  | Name n -> failwith "UFO.value_to_coupling: Name not supported yet!"

let value_to_numeric = function
  | Integer i -> Printf.sprintf "%d" i
  | Fraction (n, d) -> Printf.sprintf "%g" (float n /. float d)
  | Float x -> Printf.sprintf "%g" x
  | Expr e -> invalid_arg ("UFO.value_to_numeric: expr = " ^ (expr_to_string e))
  | Name n -> invalid_arg ("UFO.value_to_numeric: name = " ^ name_to_string n)

let value_to_float = function
  | Integer i -> float i
  | Fraction (n, d) -> float n /. float d
  | Float x -> x
  | Expr e -> invalid_arg ("UFO.value_to_float: string = " ^ (expr_to_string e))
  | Name n -> invalid_arg ("UFO.value_to_float: name = " ^ name_to_string n)

let value_attrib name attribs =
  match find_attrib name attribs with
  | S.Integer i -> Integer i
  | S.Fraction (n, d) -> Fraction (n, d)
  | S.Float x -> Float x
  | S.String s -> Expr (UFOx.Expr.of_string s)
  | S.Name n -> Name n
  | _ -> invalid_arg ("UFO.value_attrib: " ^ name)

let string_list_attrib name attribs =
  match find_attrib name attribs with
  | S.String_List l -> l
  | _ -> invalid_arg ("UFO.string_list_attrib: " ^ name)

let name_list_attrib ~strip name attribs =
  match find_attrib name attribs with
  | S.Name_List l -> List.map (name_to_string ~strip) l
  | _ -> invalid_arg ("UFO.name_list_attrib: " ^ name)

let integer_list_attrib name attribs =
  match find_attrib name attribs with
  | S.Integer_List l -> l
  | _ -> invalid_arg ("UFO.integer_list_attrib: " ^ name)

let order_dictionary_attrib name attribs =
  match find_attrib name attribs with
  | S.Order_Dictionary d -> d
  | _ -> invalid_arg ("UFO.order_dictionary_attrib: " ^ name)

let coupling_dictionary_attrib ~strip name attribs =
  match find_attrib name attribs with
  | S.Coupling_Dictionary d ->
     List.map (fun (i, j, c) -> (i, j, name_to_string ~strip c)) d
  | _ -> invalid_arg ("UFO.coupling_dictionary_attrib: " ^ name)

let decay_dictionary_attrib name attribs =
  match find_attrib name attribs with
  | S.Decay_Dictionary d ->
     List.map (fun (p, w) -> (List.map List.hd p, w)) d
  | _ -> invalid_arg ("UFO.decay_dictionary_attrib: " ^ name)

(*i The following doesn't typecheck in applications, even with
    type annotations ...
let attrib_handlers : type attribs value.
      string -> string -> attribs ->
      ((string -> attribs -> value) -> string -> value) *
        ((string -> attribs -> value) -> string -> value -> value) =
  fun kind symbol attribs ->
  let required query name =
    try
      query name attribs
    with
    | Not_found ->
       invalid_arg
         (Printf.sprintf
            "fatal UFO error: mandatory attribute `%s' missing for %s `%s'!"
            name kind symbol)
  and optional query name default =
    try
      query name attribs
    with
    | Not_found -> default in
  (required, optional) i*)

let required_handler kind symbol attribs query name =
  try
    query name attribs
  with
  | Not_found ->
     invalid_arg
       (Printf.sprintf
          "fatal UFO error: mandatory attribute `%s' missing for %s `%s'!"
          name kind symbol)

let optional_handler attribs query name default =
  try
    query name attribs
  with
  | Not_found -> default

(* The UFO paper~\cite{Degrande:2011ua} is not clear on the question
   whether the \texttt{name} attribute of an instance
   must match its Python name.
   While the examples appear to imply this, there are examples of
   UFO files in the wild that violate this constraint. *)

let warn_symbol_name file symbol name =
  if name <> symbol then
    Printf.eprintf
      "UFO: warning: symbol '%s' <> name '%s' in %s.py: \
       while legal in UFO, it is unusual and can cause problems!\n"
      symbol name file

let valid_fortran_id kind name =
  if not (ThoString.valid_fortran_id name) then
    invalid_arg
      (Printf.sprintf
         "fatal UFO error: the %s `%s' is not a valid fortran id!"
         kind name)

let map_to_alist map =
  SMap.fold (fun key value acc -> (key, value) :: acc) map []

let keys map =
  SMap.fold (fun key _ acc -> key :: acc) map []

let keys_caseless map =
  CMap.fold (fun key _ acc -> key :: acc) map []

let values map =
  SMap.fold (fun _ value acc -> value :: acc) map []

module SKey =
  struct
    type t = string
    let hash = Hashtbl.hash
    let equal = (=)
  end
module SHash = Hashtbl.Make (SKey)

module type Particle =
  sig

    type t = private
      { pdg_code : int;
	name : string;
	antiname : string;
	spin : UFOx.Lorentz.r;
	color : UFOx.Color.r;
	mass : string;
	width : string;
        propagator : string option;
	texname : string;
	antitexname : string;
	charge : charge;
	ghost_number : int;
	lepton_number : int;
	y : charge;
	goldstone : bool;
	propagating : bool;   (* NOT HANDLED YET! *)
	line : string option; (* NOT HANDLED YET! *)
        is_anti : bool }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string
    val conjugate : t -> t
    val force_spinor : t -> t
    val force_conjspinor : t -> t
    val force_majorana : t -> t
    val is_majorana : t -> bool
    val is_ghost : t -> bool
    val is_goldstone : t -> bool
    val is_physical : t -> bool
    val filter : (t -> bool) -> t SMap.t -> t SMap.t

  end

module Particle : Particle =
  struct
    
    type t =
      { pdg_code : int;
	name : string;
	antiname : string;
	spin : UFOx.Lorentz.r;
	color : UFOx.Color.r;
	mass : string;
	width : string;
        propagator : string option;
	texname : string;
	antitexname : string;
	charge : charge;
	ghost_number : int;
	lepton_number : int;
	y : charge;
	goldstone : bool;
	propagating : bool;  (* NOT HANDLED YET! *)
	line : string option; (* NOT HANDLED YET! *)
        is_anti : bool }

    let to_string symbol p =
      Printf.sprintf
	"particle: %s => [pdg = %d, name = '%s'/'%s', \
                          spin = %s, color = %s, \
                          mass = %s, width = %s,%s \
                          Q = %s, G = %d, L = %d, Y = %s, \
                          TeX = '%s'/'%s'%s]"
	symbol p.pdg_code p.name p.antiname
	(UFOx.Lorentz.rep_to_string p.spin)
	(UFOx.Color.rep_to_string p.color)
	p.mass p.width
	(match p.propagator with
         | None -> ""
         | Some p -> " propagator = " ^ p ^ ",")
	(charge_to_string p.charge)
	p.ghost_number p.lepton_number
        (charge_to_string p.y)
	p.texname p.antitexname
	(if p.goldstone then ", GB" else "")

    let conjugate_charge = function
      | Q_Integer i -> Q_Integer (-i)
      | Q_Fraction (n, d) -> Q_Fraction (-n, d)

    let is_neutral p =
      (p.name = p.antiname)

    (* We \emph{must not} mess with [pdg_code] and [color] if
       the particle is neutral! *)
    let conjugate p =
      if is_neutral p then
	p
      else
	{ pdg_code = - p.pdg_code;
	  name = p.antiname;
	  antiname = p.name;
	  spin = UFOx.Lorentz.rep_conjugate p.spin;
	  color = UFOx.Color.rep_conjugate p.color;
	  mass = p.mass;
	  width = p.width;
          propagator = p.propagator;
	  texname = p.antitexname;
	  antitexname = p.texname;
	  charge = conjugate_charge p.charge;
	  ghost_number = - p.ghost_number;
	  lepton_number = - p.lepton_number;
	  y = conjugate_charge p.y;
	  goldstone = p.goldstone;
	  propagating = p.propagating;
	  line = p.line;
          is_anti = not p.is_anti }

    let of_file1 map d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Particle" ], attribs ->
         let required query name =
           required_handler "particle" symbol attribs query name
         and optional query name default =
           optional_handler attribs query name default in
         let name = required string_attrib "name"
	 and antiname = required string_attrib "antiname" in
         let neutral = (name = antiname) in
         let pdg_code = required integer_attrib "pdg_code" in
	 SMap.add symbol
	   { (* The required attributes per UFO docs. *)
             pdg_code;
	     name; antiname;
	     spin =
               UFOx.Lorentz.rep_of_int neutral (required integer_attrib "spin");
	     color =
               UFOx.Color.rep_of_int neutral (required integer_attrib "color");
	     mass = required (name_attrib ~strip:"Param") "mass";
	     width = required (name_attrib ~strip:"Param") "width";
	     texname = required string_attrib "texname";
	     antitexname = required string_attrib "antitexname";
	     charge = required charge_attrib "charge";
	     (* The optional attributes per UFO docs. *)
             ghost_number = optional integer_attrib "GhostNumber" 0;
	     lepton_number = optional integer_attrib "LeptonNumber" 0;
	     y = optional charge_attrib "Y" (Q_Integer 0);
	     goldstone = optional boolean_attrib "goldstone" false;
	     propagating = optional boolean_attrib "propagating" true;
	     line =
               (try Some (name_attrib "line" attribs) with _ -> None);
	     (* Undocumented extensions. *)
             propagator =
               (try Some (name_attrib ~strip:"Prop" "propagator" attribs) with _ -> None);
             (* O'Mega extensions. *)
             (* Instead of ``first come is particle'' rely on
                a negative PDG code to identify antiparticles. *)
             is_anti = pdg_code < 0 } map
      | [ "anti"; p ], [] ->
	 begin
	   try
	     SMap.add symbol (conjugate (SMap.find p map)) map
	   with
	   | Not_found ->
	      invalid_arg
		("Particle.of_file: " ^ p ^ ".anti() not yet defined!")
	 end
      | _ -> invalid_arg ("Particle.of_file: " ^ name_to_string d.S.kind)

    let of_file particles =
      List.fold_left of_file1 SMap.empty particles

    let is_spinor p =
      match UFOx.Lorentz.omega p.spin with
      | Coupling.Spinor | Coupling.ConjSpinor | Coupling.Majorana -> true
      | _ -> false

    (* \begin{dubious}
         TODO: this is a bit of a hack: try to expose the type
         [UFOx.Lorentz_Atom'.r] instead.
       \end{dubious} *)
    let force_spinor p =
      if is_spinor p then
        { p with spin = UFOx.Lorentz.rep_of_int false 2 }
      else
        p

    let force_conjspinor p =
      if is_spinor p then
        { p with spin = UFOx.Lorentz.rep_of_int false (-2) }
      else
        p

    let force_majorana p =
      if is_spinor p then
        { p with spin = UFOx.Lorentz.rep_of_int true 2 }
      else
        p

    let is_majorana p =
      match UFOx.Lorentz.omega p.spin with
      | Coupling.Majorana | Coupling.Vectorspinor | Coupling.Maj_Ghost -> true
      | _ -> false

    let is_ghost p =
      p.ghost_number <> 0

    let is_goldstone p =
      p.goldstone

    let is_physical p =
      not (is_ghost p || is_goldstone p)

    let filter predicate map =
      SMap.filter (fun symbol p -> predicate p) map

  end

module type UFO_Coupling =
  sig

    type t = private
      { name : string;
	value : UFOx.Expr.t;
	order : (string * int) list }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string

  end

module UFO_Coupling : UFO_Coupling =
  struct
    
    type t =
      { name : string;
	value : UFOx.Expr.t;
	order : (string * int) list }

    let order_to_string orders =
      String.concat ", "
	(List.map (fun (s, i) -> Printf.sprintf "'%s':%d" s i) orders)

    let to_string symbol c =
      Printf.sprintf
	"coupling: %s => [name = '%s', value = '%s', order = [%s]]"
	symbol c.name (expr_to_string c.value) (order_to_string c.order)

    let of_file1 map d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Coupling" ], attribs ->
         let required query name =
           required_handler "coupling" symbol attribs query name in
         let name = required string_attrib "name" in
         warn_symbol_name "couplings" symbol name;
         valid_fortran_id "coupling" name;
	 SMap.add symbol
           { name;
	     value = UFOx.Expr.of_string (required string_attrib "value");
	     order = required order_dictionary_attrib "order" } map
      | _ -> invalid_arg ("UFO_Coupling.of_file: " ^ name_to_string d.S.kind)

    let of_file couplings =
      List.fold_left of_file1 SMap.empty couplings

  end

module type Coupling_Order =
  sig

    type t = private
      { name : string;
	expansion_order : int;
	hierarchy : int }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string

  end

module Coupling_Order : Coupling_Order =
  struct

    type t =
      { name : string;
	expansion_order : int;
	hierarchy : int }

    let to_string symbol c =
      Printf.sprintf
	"coupling_order: %s => [name = '%s', \
                                expansion_order = '%d', \
                                hierarchy = %d]"
	symbol c.name c.expansion_order c.hierarchy

    let of_file1 map d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "CouplingOrder" ], attribs ->
         let required query name =
           required_handler "coupling order" symbol attribs query name in
         let name = required string_attrib "name" in
         warn_symbol_name "coupling_orders" symbol name;
	 SMap.add symbol
	   { name;
	     expansion_order = required integer_attrib "expansion_order";
	     hierarchy = required integer_attrib "hierarchy" } map
      | _ -> invalid_arg ("Coupling_order.of_file: " ^ name_to_string d.S.kind)

    let of_file coupling_orders =
      List.fold_left of_file1 SMap.empty coupling_orders
  end

module type Lorentz_UFO =
  sig

    (* If the \texttt{name} attribute of a \texttt{Lorentz} object
       does \emph{not} match the the name of the object, we need the
       latter for weeding out unused Lorentz structures (see
       [Vertex.contains] below).  Therefore, we keep it around. *)

    type t = private
      { name : string;
        symbol : string;
	spins : int list;
	structure : UFOx.Lorentz.t }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string

  end

module Lorentz_UFO : Lorentz_UFO =
  struct

    type t =
      { name : string;
        symbol : string;
	spins : int list;
	structure : UFOx.Lorentz.t }

    let to_string symbol l =
      Printf.sprintf
	"lorentz: %s => [name = '%s', spins = [%s], \
                         structure = %s]"
	symbol l.name
	(String.concat ", " (List.map string_of_int l.spins))
	(UFOx.Lorentz.to_string l.structure)

    let of_file1 map d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Lorentz" ], attribs ->
         let required query name =
           required_handler "lorentz" symbol attribs query name in
         let name = required string_attrib "name" in
         warn_symbol_name "lorentz" symbol name;
         valid_fortran_id "lorentz" symbol;
	 SMap.add symbol
	   { name;
	     symbol;
	     spins = required integer_list_attrib "spins";
	     structure =
	       UFOx.Lorentz.of_string (required string_attrib "structure") } map
      | _ -> invalid_arg ("Lorentz.of_file: " ^ name_to_string d.S.kind)

    let of_file lorentz =
      List.fold_left of_file1 SMap.empty lorentz

  end

module type Vertex =
  sig

    type lcc = private (* Lorentz-color-coupling *)
      { lorentz : string;
	color : UFOx.Color.t;
	coupling : string }

    type t = private
      { name : string;
	particles : string array;
	lcc : lcc list }

    val of_file : Particle.t SMap.t -> S.t -> t SMap.t
    val to_string : string -> t -> string
    val to_string_expanded :
      Lorentz_UFO.t SMap.t -> UFO_Coupling.t SMap.t -> t -> string
    val contains : Particle.t SMap.t -> (Particle.t -> bool) -> t -> bool
    val filter : (t -> bool) -> t SMap.t -> t SMap.t

  end

module Vertex : Vertex =
  struct
    
    type lcc =
      { lorentz : string;
	color : UFOx.Color.t;
	coupling : string }

    type t =
      { name : string;
	particles : string array;
	lcc : lcc list }

    let to_string symbol c =
      Printf.sprintf
	"vertex: %s => [name = '%s', particles = [%s], \
                        lorentz-color-couplings = [%s]"
	symbol c.name
	(String.concat
           ", " (Array.to_list c.particles))
	(String.concat
           ", "
           (List.map
              (fun lcc ->
                Printf.sprintf
                  "%s * %s * %s"
                  lcc.coupling lcc.lorentz
                  (UFOx.Color.to_string lcc.color))
              c.lcc))
        
    let to_string_expanded lorentz couplings c =
      let expand_lorentz s =
        try
          UFOx.Lorentz.to_string (SMap.find s lorentz).Lorentz_UFO.structure
        with
        | Not_found -> "?" in
      Printf.sprintf
	"expanded: [%s] -> { lorentz-color-couplings = [%s] }"
	(String.concat ", " (Array.to_list c.particles))
        (String.concat
           ", "
           (List.map
              (fun lcc ->
                Printf.sprintf
                  "%s * %s * %s"
                  lcc.coupling (expand_lorentz lcc.lorentz)
                  (UFOx.Color.to_string lcc.color))
              c.lcc))

    let contains particles predicate v =
      let p = v.particles in
      let rec contains' i =
	if i < 0 then
	  false
	else if predicate (SMap.find p.(i) particles) then
	  true
	else
	  contains' (pred i) in
      contains' (Array.length p - 1)
      
    let force_adj_identity1 adj_indices = function
      | UFOx.Color_Atom.Identity (a, b) as atom ->
         begin match List.mem a adj_indices, List.mem b adj_indices with
         | true, true -> UFOx.Color_Atom.Identity8 (a, b)
         | false, false -> atom
         | true, false | false, true ->
            invalid_arg "force_adj_identity: mixed representations!"
         end
      | atom -> atom

    let force_adj_identity adj_indices tensor =
      UFOx.Color.map_atoms (force_adj_identity1 adj_indices) tensor

    let find_adj_indices map particles =
      let adj_indices = ref [] in
      Array.iteri
        (fun i p ->
          (* We must pattern match against the O'Mega representation,
             because [UFOx.Color.r] is abstract. *)
          match UFOx.Color.omega (SMap.find p map).Particle.color with
          | Color.AdjSUN _ -> adj_indices := succ i :: !adj_indices
          | _ -> ())
        particles;
      !adj_indices

    let classify_color_indices map particles =
      let fund_indices = ref []
      and conj_indices = ref []
      and adj_indices = ref [] in
      Array.iteri
        (fun i p ->
          (* We must pattern match against the O'Mega representation,
             because [UFOx.Color.r] is abstract. *)
          match UFOx.Color.omega (SMap.find p map).Particle.color with
          | Color.SUN n ->
             if n > 0 then
               fund_indices := succ i :: !fund_indices
             else if n < 0 then
               conj_indices := succ i :: !conj_indices
             else
               failwith "classify_color_indices: SU(0)"
          | Color.AdjSUN n ->
             if n <> 0 then
               adj_indices := succ i :: !adj_indices
             else
               failwith "classify_color_indices: SU(0)"
          | _ -> ())
        particles;
      (!fund_indices, !conj_indices, !adj_indices)

    (* FIXME: would have expected the opposite order \ldots *)
    let force_identity1 (fund_indices, conj_indices, adj_indices) = function
      | UFOx.Color_Atom.Identity (a, b) as atom ->
         if List.mem a fund_indices then
           begin
             if List.mem b conj_indices then
               UFOx.Color_Atom.Identity (b, a)
             else
               invalid_arg "force_adj_identity: mixed representations!"
           end
         else if List.mem a conj_indices then
           begin
             if List.mem b fund_indices then
               UFOx.Color_Atom.Identity (a, b)
             else
               invalid_arg "force_adj_identity: mixed representations!"
           end else if List.mem a adj_indices then begin
             if List.mem b adj_indices then
               UFOx.Color_Atom.Identity8 (a, b)
             else
               invalid_arg "force_adj_identity: mixed representations!"
           end
         else
           atom
      | atom -> atom

    let force_identity indices tensor =
      UFOx.Color.map_atoms (force_identity1 indices) tensor

    (* Here we don't have the Lorentz structures available yet.
       Thus we set [fermion_lines = []] for now and correct this
       later. *)
    let of_file1 particle_map map d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Vertex" ], attribs ->
         let required query name =
           required_handler "vertex" symbol attribs query name in
         let name = required string_attrib "name" in
         warn_symbol_name "vertices" symbol name;
         let particles =
	   Array.of_list (required (name_list_attrib ~strip:"P") "particles") in
	 let color =
           let indices = classify_color_indices particle_map particles in
	   Array.of_list
	     (List.map
                (force_identity indices <*> UFOx.Color.of_string)
                (required string_list_attrib "color"))
	 and lorentz =
	   Array.of_list (required (name_list_attrib ~strip:"L") "lorentz")
	 and couplings_alist =
	   required (coupling_dictionary_attrib ~strip:"C") "couplings" in
	 let lcc =
	   List.map
	     (fun (i, j, c) ->
               { lorentz = lorentz.(j);
                 color = color.(i);
                 coupling = c })
	     couplings_alist in
	 SMap.add symbol { name; particles; lcc } map
      | _ -> invalid_arg ("Vertex.of_file: " ^ name_to_string d.S.kind)

    let of_file particles vertices =
      List.fold_left (of_file1 particles) SMap.empty vertices

    let filter predicate map =
      SMap.filter (fun symbol p -> predicate p) map

  end

module type Parameter =
  sig

    type nature = private Internal | External
    type ptype = private Real | Complex

    type t = private
      { name : string;
	nature : nature;
	ptype : ptype;
	value : value;
	texname : string;
	lhablock : string option;
	lhacode : int list option;
        sequence : int }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string

    val missing : string -> t

  end

module Parameter : Parameter =
  struct

    type nature = Internal | External
	
    let nature_to_string = function
      | Internal -> "internal"
      | External -> "external"

    let nature_of_string = function
      | "internal" -> Internal
      | "external" -> External
      | s -> invalid_arg ("Parameter.nature_of_string: " ^ s)
	 
    type ptype = Real | Complex

    let ptype_to_string = function
      | Real -> "real"
      | Complex -> "complex"

    let ptype_of_string = function
      | "real" -> Real
      | "complex" -> Complex
      | s -> invalid_arg ("Parameter.ptype_of_string: " ^ s)

    type t =
      { name : string;
	nature : nature;
	ptype : ptype;
	value : value;
	texname : string;
	lhablock : string option;
	lhacode : int list option;
        sequence : int }

    let to_string symbol p =
      Printf.sprintf
	"parameter: %s => [#%d, name = '%s', nature = %s, type = %s, \
                           value = %s, texname = '%s', \
                           lhablock = %s, lhacode = [%s]]"
	symbol p.sequence p.name
	(nature_to_string p.nature)
	(ptype_to_string p.ptype)
	(value_to_string p.value) p.texname
	(match p.lhablock with None -> "???" | Some s -> s)
	(match p.lhacode with
	| None -> ""
	| Some c -> String.concat ", " (List.map string_of_int c))
      
    let of_file1 (map, n) d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Parameter" ], attribs ->
         let required query name =
           required_handler "particle" symbol attribs query name in
         let name = required string_attrib "name" in
         warn_symbol_name "parameters" symbol name;
         valid_fortran_id "parameter" name;
	 (SMap.add symbol
	    { name;
	      nature = nature_of_string (required string_attrib "nature");
	      ptype = ptype_of_string (required string_attrib "type");
	      value = required value_attrib "value";
	      texname = required string_attrib "texname";
	      lhablock =
	        (try Some (string_attrib "lhablock" attribs) with
		   Not_found -> None);
	      lhacode =
	        (try Some (integer_list_attrib "lhacode" attribs) with
		   Not_found -> None);
              sequence = n } map, succ n)
      | _ -> invalid_arg ("Parameter.of_file: " ^ name_to_string d.S.kind)
    
    let of_file parameters =
      let map, _ = List.fold_left of_file1 (SMap.empty, 0) parameters in
      map

    let missing name =
      { name;
	nature = External;
	ptype = Real;
	value = Integer 0;
	texname = Printf.sprintf "\\texttt{%s}" name;
	lhablock = None;
	lhacode = None;
        sequence = 0 }

  end

(* Macros are encoded as a special [S.declaration] with
   [S.kind = "$"].  This is slightly hackish, but general enough
   and the overhead of a special union type is probably not worth
   the effort.  *)

module type Macro =
  sig
    type t
    val empty : t

    (* The domains and codomains are still a bit too much ad hoc,
       but it does the job. *)
    val define : t -> string -> S.value -> t
    val expand_string : t -> string -> S.value
    val expand_expr : t -> S.string_atom list -> string

    (* Only for documentation: *)
    val expand_atom : t -> S.string_atom -> string
  end

module Macro : Macro =
  struct

    type t = S.value SMap.t

    let empty = SMap.empty

    let define macros name expansion =
      SMap.add name expansion macros

    let expand_string macros name =
      SMap.find name macros

    let rec expand_atom macros = function
      | S.Literal s -> s
      | S.Macro [name] ->
         begin
           try
             begin match SMap.find name macros with
             | S.String s -> s
             | S.String_Expr expr -> expand_expr macros expr
             | _ -> invalid_arg ("expand_atom: not a string: " ^ name)
             end
           with
           | Not_found -> invalid_arg ("expand_atom: not found: " ^ name)
         end
      | S.Macro [] -> invalid_arg "expand_atom: empty"
      | S.Macro name ->
         invalid_arg ("expand_atom: compound name: " ^ String.concat "." name)

    and expand_expr macros expr =
      String.concat "" (List.map (expand_atom macros) expr)

  end

module type Propagator_UFO =
  sig

    type t = (* private *)
      { name : string;
	numerator : UFOx.Lorentz.t;
	denominator : UFOx.Lorentz.t }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string

  end

module Propagator_UFO : Propagator_UFO =
  struct

    type t =
      { name : string;
	numerator : UFOx.Lorentz.t;
	denominator : UFOx.Lorentz.t }

    let to_string symbol p =
      Printf.sprintf
	"propagator: %s => [name = '%s', numerator = '%s', \
                            denominator = '%s']"
	symbol p.name
        (UFOx.Lorentz.to_string p.numerator)
        (UFOx.Lorentz.to_string p.denominator)

    (* The \texttt{denominator} attribute is optional and
       there is a default (cf.~\texttt{arXiv:1308.1668}) *)
    let default_denominator =
      "P('mu', id) * P('mu', id) \
       - Mass(id) * Mass(id) \
       + complex(0,1) * Mass(id) * Width(id)"

    let of_string_with_error_correction symbol num_or_den s =
      try
        UFOx.Lorentz.of_string s
      with
      | Invalid_argument msg ->
         begin
           let fixed = s ^ ")" in
           try
             let tensor = UFOx.Lorentz.of_string fixed in
             Printf.eprintf
               "UFO.Propagator.of_string: added missing closing parenthesis \
                in %s of %s: \"%s\"\n"
               num_or_den symbol s;
             tensor
           with
           | Invalid_argument _ ->
              invalid_arg
                (Printf.sprintf
                   "UFO.Propagator.of_string: %s of %s: %s in \"%s\"\n"
                   num_or_den symbol msg fixed)
         end

    let of_file1 (macros, map) d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Propagator" ], attribs ->
         let required query name =
           required_handler "particle" symbol attribs query name
         and optional query name default =
           optional_handler attribs query name default in
        let name = required string_attrib "name" in
         warn_symbol_name "propagators" symbol name;
         let num_string_expr = required string_expr_attrib "numerator"
         and den_string =
	   begin match optional find_attrib "denominator"
                                (S.String default_denominator) with
	   | S.String s -> s
	   | S.Name [n] ->
              begin match Macro.expand_string macros n with
              | S.String s -> s
              | _ -> invalid_arg "Propagator.denominator"
              end
	   | _ -> invalid_arg "Propagator.denominator: "
	   end in
         let num_string = Macro.expand_expr macros num_string_expr in
         let numerator =
           of_string_with_error_correction symbol "numerator" num_string
         and denominator =
           of_string_with_error_correction symbol "denominator" den_string in
	 (macros, SMap.add symbol { name; numerator; denominator } map)
      | [ "$" ], [ macro ] ->
         begin match macro.S.a_value with
         | S.String _ as s ->
            (Macro.define macros symbol s, map);
         | S.String_Expr expr ->
            let expanded = S.String (Macro.expand_expr macros expr) in
            (Macro.define macros symbol expanded, map)
         | _ -> invalid_arg ("Propagator:of_file: not a string " ^ symbol)
         end
      | [ "$" ], [] ->
         invalid_arg ("Propagator:of_file: empty declaration " ^ symbol)
      | [ "$" ], _ ->
         invalid_arg ("Propagator:of_file: multiple declaration " ^ symbol)
      | _ -> invalid_arg ("Propagator:of_file: " ^ name_to_string d.S.kind)
       
    let of_file propagators =
      let _, propagators' =
	List.fold_left of_file1 (Macro.empty, SMap.empty) propagators in
      propagators'

  end

module type Decay =
  sig

    type t = private
      { name : string;
	particle : string;
	widths : (string list * string) list }

    val of_file : S.t -> t SMap.t
    val to_string : string -> t -> string

  end

module Decay : Decay =
  struct

    type t =
      { name : string;
	particle : string;
	widths : (string list * string) list }

    let width_to_string ws =
      String.concat ", "
	(List.map
	   (fun (ps, w) ->
	     "(" ^ String.concat ", " ps ^ ") -> '" ^ w ^ "'")
	   ws)

    let to_string symbol d =
      Printf.sprintf
	"decay: %s => [name = '%s', particle = '%s', widths = [%s]]"
	symbol d.name d.particle (width_to_string d.widths)

    let of_file1 map d =
      let symbol = d.S.name in
      match d.S.kind, d.S.attribs with
      | [ "Decay" ], attribs ->
         let required query name =
           required_handler "particle" symbol attribs query name in
         let name = required string_attrib "name" in
         warn_symbol_name "decays" symbol name;
	 SMap.add symbol
	   { name;
	     particle = required (name_attrib ~strip:"P") "particle";
	     widths = required decay_dictionary_attrib "partial_widths" } map
      | _ -> invalid_arg ("Decay.of_file: " ^ name_to_string d.S.kind)

    let of_file decays =
      List.fold_left of_file1 SMap.empty decays

  end

(* We can read the spinor representations off the
   vertices to check for consistency. *)
(* \begin{dubious}
     Note that we have to conjugate the representations!
   \end{dubious} *)

let collect_spinor_reps_of_vertex particles lorentz v sets =
  List.fold_left
    (fun sets' lcc ->
      let l = (SMap.find lcc.Vertex.lorentz lorentz).Lorentz_UFO.structure in
      List.fold_left
        (fun (spinors, conj_spinors as sets'') (i, rep) ->
          let p = v.Vertex.particles.(pred i) in
          match UFOx.Lorentz.omega rep with
          | Coupling.ConjSpinor -> (SSet.add p spinors, conj_spinors)
          | Coupling.Spinor -> (spinors, SSet.add p conj_spinors)
          | _ -> sets'')
        sets' (UFOx.Lorentz.classify_indices l))
    sets v.Vertex.lcc

let collect_spinor_reps_of_vertices particles lorentz vertices =
  SMap.fold
    (fun _ v -> collect_spinor_reps_of_vertex particles lorentz v)
    vertices (SSet.empty, SSet.empty)

let lorentz_reps_of_vertex particles v =
  ThoList.alist_of_list ~predicate:(not <*> UFOx.Lorentz.rep_trivial) ~offset:1
    (List.map
       (fun p ->
	 (* Why do we need to conjugate??? *)
         UFOx.Lorentz.rep_conjugate
           (SMap.find p particles).Particle.spin)
       (Array.to_list v.Vertex.particles))

let rep_compatible rep_vertex rep_particle =
  let open UFOx.Lorentz in
  let open Coupling in
  match omega rep_vertex, omega rep_particle with
  | (Spinor | ConjSpinor), Majorana -> true
  | r1, r2 -> r1 = r2

let reps_compatible reps_vertex reps_particles =
  List.for_all2
    (fun (iv, rv) (ip, rp) -> iv = ip && rep_compatible rv rp)
    reps_vertex reps_particles

let check_lorentz_reps_of_vertex particles lorentz v =
  let reps_particles =
    List.sort compare (lorentz_reps_of_vertex particles v) in
  List.iter
    (fun lcc ->
      let l = (SMap.find lcc.Vertex.lorentz lorentz).Lorentz_UFO.structure in
      let reps_vertex = List.sort compare (UFOx.Lorentz.classify_indices l) in
      if not (reps_compatible reps_vertex reps_particles) then begin
	Printf.eprintf "%s <> %s [%s]\n"
	  (UFOx.Index.classes_to_string
	     UFOx.Lorentz.rep_to_string reps_particles)
	  (UFOx.Index.classes_to_string
	     UFOx.Lorentz.rep_to_string reps_vertex)
          v.Vertex.name (* [(Vertex.to_string v.Vertex.name v)] *);
	(* [invalid_arg "check_lorentz_reps_of_vertex"] *) ()
      end)
    v.Vertex.lcc

let color_reps_of_vertex particles v =
  ThoList.alist_of_list ~predicate:(not <*> UFOx.Color.rep_trivial) ~offset:1
    (List.map
       (fun p -> (SMap.find p particles).Particle.color)
       (Array.to_list v.Vertex.particles))

let check_color_reps_of_vertex particles v =
  let reps_particles =
    List.sort compare (color_reps_of_vertex particles v) in
  List.iter
    (fun lcc ->
      let reps_vertex =
        List.sort compare (UFOx.Color.classify_indices lcc.Vertex.color) in
      if reps_vertex <> reps_particles then begin
	Printf.printf "%s <> %s\n"
	  (UFOx.Index.classes_to_string UFOx.Color.rep_to_string reps_particles)
	  (UFOx.Index.classes_to_string UFOx.Color.rep_to_string reps_vertex);
	invalid_arg "check_color_reps_of_vertex"
     end)
    v.Vertex.lcc

module P = Permutation.Default

module type Lorentz =
  sig

    type spins = private
      | Unused
      | Unique of Coupling.lorentz array
      | Ambiguous of Coupling.lorentz array SMap.t

    type t = private
      { name : string;
        n : int;
	spins : spins;
	structure : UFO_Lorentz.t;
        fermion_lines : Coupling.fermion_lines;
        variables : string list }

    val required_charge_conjugates : t -> t list
    val permute : P.t -> t -> t

    val of_lorentz_UFO :
      Particle.t SMap.t -> Vertex.t SMap.t ->
      Lorentz_UFO.t SMap.t -> t SMap.t

    val lorentz_to_string : Coupling.lorentz -> string
    val to_string : string -> t -> string

  end

module Lorentz : Lorentz =
  struct

    let rec lorentz_to_string = function
      | Coupling.Scalar -> "Scalar"
      | Coupling.Spinor -> "Spinor"
      | Coupling.ConjSpinor -> "ConjSpinor"
      | Coupling.Majorana -> "Majorana"
      | Coupling.Maj_Ghost -> "Maj_Ghost"
      | Coupling.Vector -> "Vector"
      | Coupling.Massive_Vector -> "Massive_Vector"
      | Coupling.Vectorspinor -> "Vectorspinor"
      | Coupling.Tensor_1 -> "Tensor_1"
      | Coupling.Tensor_2 -> "Tensor_2"
      | Coupling.BRS l -> "BRS(" ^ lorentz_to_string l ^ ")"

    (* Unlike UFO, O'Mega distinguishes bewteen spinors
       and conjugate spinors.  However, we can inspect
       the particles in the vertices in which a Lorentz
       structure is used to determine the correct
       quantum numbers.

       Most model files in the real world contain unused Lorentz
       structures.  This is not a problem, we can just ignore them. *)

    type spins =
      | Unused
      | Unique of Coupling.lorentz array
      | Ambiguous of Coupling.lorentz array SMap.t

    (* \begin{dubious}
         Use [UFO_targets.Fortran.fusion_name] below in order
         to avoid communication problems.  Or even move away
         from strings alltogether.
       \end{dubious} *)
    type t =
      { name : string;
        n : int;
	spins : spins;
	structure : UFO_Lorentz.t;
        fermion_lines : Coupling.fermion_lines;
        variables : string list }

    (* Add one charge conjugated fermion lines. *)
    let charge_conjugate1 l (ket, bra as fermion_line) =
      { name = l.name ^ Printf.sprintf "_c%x%x" ket bra;
        n = l.n;
        spins = l.spins;
        structure = UFO_Lorentz.charge_conjugate fermion_line l.structure;
        fermion_lines = l.fermion_lines;
        variables = l.variables }

    (* Add several charge conjugated fermion lines. *)
    let charge_conjugate l fermion_lines =
      List.fold_left charge_conjugate1 l fermion_lines

(*i
    let all_charge_conjugates l =
      List.map (charge_conjugate l) (ThoList.power l.fermion_lines)
i*)

    (* Add all combinations of charge conjugated fermion lines
       that don't leave the fusion. *)
    let required_charge_conjugates l =
      let saturated_fermion_lines =
        List.filter
          (fun (ket, bra) -> ket != 1 && bra != 1)
          l.fermion_lines in
      List.map (charge_conjugate l) (ThoList.power saturated_fermion_lines)

    let permute_spins p = function
      | Unused -> Unused
      | Unique s -> Unique (P.array p s)
      | Ambiguous map -> Ambiguous (SMap.map (P.array p) map)

    (* Note that we apply the \emph{inverse} permutation to
       the indices in order to match the permutation of the
       particles/spins. *)

    let permute_structure n p (l, f) =
      let permuted = P.array (P.inverse p) (Array.init n succ) in
      let permute_index i =
        if i > 0 then
          UFOx.Index.map_position (fun pos -> permuted.(pred pos)) i
        else
          i in
      (UFO_Lorentz.map_indices permute_index l,
       UFO_Lorentz.map_fermion_lines permute_index f)

    let permute p l =
      let structure, fermion_lines =
        permute_structure l.n p (l.structure, l.fermion_lines) in
      { name = l.name ^ "_p" ^ P.to_string (P.inverse p);
        n = l.n;
        spins = permute_spins p l.spins;
        structure;
        fermion_lines;
        variables = l.variables }

    let omega_lorentz_reps n alist =
      let reps = Array.make n Coupling.Scalar in
      List.iter
        (fun (i, rep) -> reps.(pred i) <- UFOx.Lorentz.omega rep)
        alist;
      reps

    let contained lorentz vertex =
      List.exists
        (fun lcc1 -> lcc1.Vertex.lorentz = lorentz.Lorentz_UFO.symbol)
        vertex.Vertex.lcc

    (* Find all vertices in with the Lorentz structure [lorentz] is
       used and build a map from those vertices to the O'Mega
       Lorentz representations inferred from UFO's Lorentz
       structure and the [particles] involved.
       Then scan the bindings and check that we have inferred
       the same Lorentz representation from all vertices. *)
    let lorentz_reps_of_structure particles vertices lorentz =
      let uses =
        SMap.fold
          (fun name v acc ->
            if contained lorentz v then
              SMap.add
                name
                (omega_lorentz_reps
                   (Array.length v.Vertex.particles)
                   (lorentz_reps_of_vertex particles v)) acc
            else
              acc) vertices SMap.empty in
      let variants =
        ThoList.uniq (List.sort compare (List.map snd (SMap.bindings uses))) in
      match variants with
      | [] -> Unused
      | [s] -> Unique s
      | _ ->
         Printf.eprintf "UFO.Lorentz.lorentz_reps_of_structure: AMBIGUOUS!\n";
         List.iter
           (fun variant ->
             Printf.eprintf
               "UFO.Lorentz.lorentz_reps_of_structure: %s\n"
               (ThoList.to_string lorentz_to_string (Array.to_list variant)))
           variants;
         Ambiguous uses

    let of_lorentz_tensor spins lorentz =
      match spins with
      | Unique s ->
         begin
           try
             Some (UFO_Lorentz.parse (Array.to_list s) lorentz)
           with
           | Failure msg ->
              begin
                prerr_endline msg;
                Some (UFO_Lorentz.dummy)
              end
         end
      | Unused ->
         Printf.eprintf
           "UFO.Lorentz: stripping unused structure %s\n"
           (UFOx.Lorentz.to_string lorentz);
         None
      | Ambiguous _ -> invalid_arg "UFO.Lorentz.of_lorentz_tensor: Ambiguous"

    (* NB: if the \texttt{name} attribute of a \texttt{Lorentz} object
       does \emph{not} match the the name of the object, the former has
       a better chance to correspond to a valid Fortran name.  Therefore
       we use it. *)

    let of_lorentz_UFO particles vertices lorentz_UFO =
      SMap.fold
        (fun name l acc ->
          let spins = lorentz_reps_of_structure particles vertices l in
          match of_lorentz_tensor spins l.Lorentz_UFO.structure with
          | None -> acc
          | Some structure ->
             SMap.add
               name
               { name = l.Lorentz_UFO.symbol;
                 n = List.length l.Lorentz_UFO.spins;
	         spins;
	         structure;
                 fermion_lines = UFO_Lorentz.fermion_lines structure;
                 variables = UFOx.Lorentz.variables l.Lorentz_UFO.structure }
               acc)
        lorentz_UFO SMap.empty

    let to_string symbol l =
      Printf.sprintf
	"lorentz: %s => [name = '%s', spins = %s, \
                         structure = %s, fermion_lines = %s]"
	symbol l.name
	(match l.spins with
         | Unique s ->
            "[" ^ String.concat
                    ", " (List.map lorentz_to_string (Array.to_list s)) ^ "]"
         | Ambiguous _ -> "AMBIGUOUS!"
         | Unused -> "UNUSED!")
	(UFO_Lorentz.to_string l.structure)
	(UFO_Lorentz.fermion_lines_to_string l.fermion_lines)

  end

(* According to arxiv:1308:1668, there should not be a factor
   of~$i$ in the numerators of propagators, but the (unused)
   \texttt{propagators.py} in most models violate this rule! *)
let divide_propagators_by_i = ref false

module type Propagator =
  sig

    type t = (* private *)
      { name : string;
        spins : Coupling.lorentz * Coupling.lorentz;
	numerator : UFO_Lorentz.t;
	denominator : UFO_Lorentz.t;
        variables : string list }

    val of_propagator_UFO : ?majorana:bool -> Propagator_UFO.t -> t
    val of_propagators_UFO : ?majorana:bool -> Propagator_UFO.t SMap.t -> t SMap.t

    val transpose : t -> t

    val to_string : string -> t -> string

  end

module Propagator : Propagator =
  struct

    type t = (* private *)
      { name : string;
        spins : Coupling.lorentz * Coupling.lorentz;
	numerator : UFO_Lorentz.t;
	denominator : UFO_Lorentz.t;
	variables : string list }

    let lorentz_rep_at rep_classes i =
      try
        UFOx.Lorentz.omega (List.assoc i rep_classes)
      with
      | Not_found -> Coupling.Scalar

    let imaginary = Algebra.QC.make Algebra.Q.null Algebra.Q.unit
    let scalars = [Coupling.Scalar; Coupling.Scalar]

    (* If~$51$ and~$52$ show up as indices, we must
       map $(1,51)\to(1001,2001)$ and $(2,52)\to(1002,2002)$,
       as per the UFO conventions for Lorentz structures. *)

    (* \begin{dubious}
         This does not work yet, because [UFOx.Lorentz.map_indices]
         affects also the position argument of [P], [Mass] and [Width].
       \end{dubious} *)

    let contains_51_52 tensor =
      List.exists
        (fun (i, _) -> i = 51 || i = 52)
        (UFOx.Lorentz.classify_indices tensor)

    let remap_51_52 = function
      | 1 -> 1001 | 51 -> 2001
      | 2 -> 1002 | 52 -> 2002
      | i -> i

    let canonicalize_51_52 tensor =
      if contains_51_52 tensor then
        UFOx.Lorentz.rename_indices remap_51_52 tensor
      else
        tensor

    let force_majorana = function
      | Coupling.Spinor | Coupling.ConjSpinor -> Coupling.Majorana
      | s -> s

    let string_list_union l1 l2 =
      Sets.String.elements
        (Sets.String.union
           (Sets.String.of_list l1)
           (Sets.String.of_list l2))

    (* In the current conventions, the factor of~$i$ is not included: *)
    let of_propagator_UFO ?(majorana=false) p =
      let numerator = canonicalize_51_52 p.Propagator_UFO.numerator in
      let lorentz_reps = UFOx.Lorentz.classify_indices numerator in
      let spin1 = lorentz_rep_at lorentz_reps 1
      and spin2 = lorentz_rep_at lorentz_reps 2 in
      let numerator_sans_i =
        if !divide_propagators_by_i then
          UFOx.Lorentz.map_coeff (fun q -> Algebra.QC.div q imaginary) numerator
        else
          numerator in
      { name = p.Propagator_UFO.name;
        spins =
          if majorana then
            (force_majorana spin1, force_majorana spin2)
          else
            (spin1, spin2);
        numerator =
          UFO_Lorentz.parse ~allow_denominator:true [spin1; spin2] numerator_sans_i;
        denominator = UFO_Lorentz.parse scalars p.Propagator_UFO.denominator;
        variables =
          string_list_union
            (UFOx.Lorentz.variables p.Propagator_UFO.denominator)
            (UFOx.Lorentz.variables numerator_sans_i) }

    let of_propagators_UFO ?majorana propagators_UFO =
      SMap.fold
        (fun name p acc -> SMap.add name (of_propagator_UFO ?majorana p) acc)
        propagators_UFO SMap.empty

    let permute12 = function
      | 1 -> 2
      | 2 -> 1
      | n -> n

    let transpose_positions t =
      UFOx.Index.map_position permute12 t

    let transpose p =
      { name = p.name;
        spins = (snd p.spins, fst p.spins);
        numerator = UFO_Lorentz.map_indices transpose_positions p.numerator;
        denominator = p.denominator;
        variables = p.variables }

    let to_string symbol p =
      Printf.sprintf
	"propagator: %s => [name = '%s', spin = '(%s, %s)', numerator/I = '%s', \
                            denominator = '%s']"
	symbol p.name
        (Lorentz.lorentz_to_string (fst p.spins))
        (Lorentz.lorentz_to_string (snd p.spins))
        (UFO_Lorentz.to_string p.numerator)
        (UFO_Lorentz.to_string p.denominator)

  end

type t =
  { particles : Particle.t SMap.t;
    particle_array : Particle.t array; (* for diagnostics *)
    couplings : UFO_Coupling.t SMap.t;
    coupling_orders : Coupling_Order.t SMap.t;
    vertices : Vertex.t SMap.t;
    lorentz_UFO : Lorentz_UFO.t SMap.t;
    lorentz : Lorentz.t SMap.t;
    parameters : Parameter.t SMap.t;
    propagators_UFO : Propagator_UFO.t SMap.t;
    propagators : Propagator.t SMap.t;
    decays : Decay.t SMap.t;
    nc : int }

let use_majorana_spinors = ref false

let fallback_to_majorana_if_necessary particles vertices lorentz_UFO =
  let majoranas =
    SMap.fold
      (fun p particle acc ->
        if Particle.is_majorana particle then
          SSet.add p acc
        else
          acc)
      particles SSet.empty in
  let spinors, conj_spinors =
    collect_spinor_reps_of_vertices particles lorentz_UFO vertices in
  let ambiguous =
    SSet.diff (SSet.inter spinors conj_spinors) majoranas in
  let no_majoranas = SSet.is_empty majoranas
  and no_ambiguities = SSet.is_empty ambiguous in
  if no_majoranas && no_ambiguities && not !use_majorana_spinors then
    (SMap.mapi
       (fun p particle ->
         if SSet.mem p spinors then
           Particle.force_spinor particle
         else if SSet.mem p conj_spinors then
           Particle.force_conjspinor particle
         else
           particle)
       particles,
     false)
  else
    begin
      if !use_majorana_spinors then
        Printf.eprintf "O'Mega: Majorana fermions requested.\n";
      if not no_majoranas then
        Printf.eprintf "O'Mega: found Majorana fermions!\n";
      if not no_ambiguities then
        Printf.eprintf
          "O'Mega: found ambiguous spinor representations for %s!\n"
          (String.concat ", " (SSet.elements ambiguous));
      Printf.eprintf
        "O'Mega: falling back to the Majorana representation for all fermions.\n";
      (SMap.map Particle.force_majorana particles,
       true)
    end

let nc_of_particles particles =
  let nc_set =
    List.fold_left
      (fun nc_set (_, p) ->
        match UFOx.Color.omega p.Particle.color with
        | Color.Singlet -> nc_set
        | Color.SUN nc -> Sets.Int.add (abs nc) nc_set
        | Color.AdjSUN nc -> Sets.Int.add (abs nc) nc_set)
      Sets.Int.empty (SMap.bindings particles) in
  match Sets.Int.elements nc_set with
  | [] -> 0
  | [n] -> n
  | nc_list ->
     invalid_arg
       ("UFO.Model: more than one value of N_C: " ^
          String.concat ", " (List.map string_of_int nc_list))

let of_file u =
  let particles = Particle.of_file u.Files.particles in
  let vertices = Vertex.of_file particles u.Files.vertices
  and lorentz_UFO = Lorentz_UFO.of_file u.Files.lorentz
  and propagators_UFO = Propagator_UFO.of_file u.Files.propagators in
  let particles, majorana =
    fallback_to_majorana_if_necessary particles vertices lorentz_UFO in
  let particle_array = Array.of_list (values particles)
  and lorentz = Lorentz.of_lorentz_UFO particles vertices lorentz_UFO
  and propagators = Propagator.of_propagators_UFO ~majorana propagators_UFO in
  let model =
    { particles;
      particle_array;
      couplings = UFO_Coupling.of_file u.Files.couplings;
      coupling_orders = Coupling_Order.of_file u.Files.coupling_orders;
      vertices;
      lorentz_UFO;
      lorentz;
      parameters = Parameter.of_file u.Files.parameters;
      propagators_UFO;
      propagators;
      decays = Decay.of_file u.Files.decays;
      nc = nc_of_particles particles } in
  SMap.iter
    (fun _ v ->
      check_color_reps_of_vertex model.particles v;
      check_lorentz_reps_of_vertex model.particles model.lorentz_UFO v)
    model.vertices;
  model

let parse_directory dir =
  of_file (Files.parse_directory dir)

let dump model =
  Printf.printf "NC = %d\n" model.nc;
  SMap.iter (print_endline <**> Particle.to_string) model.particles;
  SMap.iter (print_endline <**> UFO_Coupling.to_string) model.couplings;
  SMap.iter (print_endline <**> Coupling_Order.to_string) model.coupling_orders;
  (* [SMap.iter (print_endline <**> Vertex.to_string) model.vertices;] *)
  SMap.iter
    (fun symbol v ->
      (print_endline <**> Vertex.to_string) symbol v;
      print_endline
        (Vertex.to_string_expanded model.lorentz_UFO model.couplings v))
    model.vertices;
  SMap.iter (print_endline <**> Lorentz_UFO.to_string) model.lorentz_UFO;
  SMap.iter (print_endline <**> Lorentz.to_string) model.lorentz;
  SMap.iter (print_endline <**> Parameter.to_string) model.parameters;
  SMap.iter (print_endline <**> Propagator_UFO.to_string) model.propagators_UFO;
  SMap.iter (print_endline <**> Propagator.to_string) model.propagators;
  SMap.iter (print_endline <**> Decay.to_string) model.decays;
  SMap.iter
    (fun symbol d ->
      List.iter (fun (_, w) -> ignore (UFOx.Expr.of_string w)) d.Decay.widths)
    model.decays

exception Unhandled of string
let unhandled s = raise (Unhandled s)

module Model =
  struct

    (* NB: we could use [type flavor = Particle.t], but that would
       be very inefficient, because we will use [flavor] as a key
       for maps below. *)
    type flavor = int
    type constant = string
    type gauge = unit

    module M = Modeltools.Mutable
        (struct type f = flavor type g = gauge type c = constant end)

    let flavors = M.flavors
    let external_flavors = M.external_flavors
    let external_flavors = M.external_flavors
    let lorentz = M.lorentz
    let color = M.color
    let nc = M.nc
    let propagator = M.propagator
    let width = M.width
    let goldstone = M.goldstone
    let conjugate = M.conjugate
    let fermion = M.fermion
    let vertices = M.vertices
    let fuse2 = M.fuse2
    let fuse3 = M.fuse3
    let fuse = M.fuse
    let max_degree = M.max_degree
    let parameters = M.parameters
    let flavor_of_string = M.flavor_of_string
    let flavor_to_string = M.flavor_to_string
    let flavor_to_TeX = M.flavor_to_TeX
    let flavor_symbol = M.flavor_symbol
    let gauge_symbol = M.gauge_symbol
    let pdg = M.pdg
    let mass_symbol = M.mass_symbol
    let width_symbol = M.width_symbol
    let constant_symbol = M.constant_symbol
    module Ch = M.Ch
    let charges = M.charges

    let rec fermion_of_lorentz = function
      | Coupling.Spinor -> 1
      | Coupling.ConjSpinor -> -1
      | Coupling.Majorana -> 2
      | Coupling.Maj_Ghost -> 2
      | Coupling.Vectorspinor -> 1
      | Coupling.Vector | Coupling.Massive_Vector -> 0
      | Coupling.Scalar | Coupling.Tensor_1 | Coupling.Tensor_2 -> 0 
      | Coupling.BRS f -> fermion_of_lorentz f

    module Q = Algebra.Q
    module QC = Algebra.QC

    let dummy_tensor3 = Coupling.Scalar_Scalar_Scalar 1
    let dummy_tensor4 = Coupling.Scalar4 1

    let triplet p = (p.(0), p.(1), p.(2))
    let quartet p = (p.(0), p.(1), p.(2), p.(3))

    let half_times q1 q2 =
      Q.mul (Q.make 1 2) (Q.mul q1 q2)

    let name g =
      g.UFO_Coupling.name

    let fractional_coupling g r =
      let g = name g in
      match Q.to_ratio r with
      |  0, _ -> "0.0_default"
      |  1, 1 -> g
      | -1, 1 -> Printf.sprintf "(-%s)" g
      |  n, 1 -> Printf.sprintf "(%d*%s)" n g
      |  1, d -> Printf.sprintf "(%s/%d)" g d
      | -1, d -> Printf.sprintf "(-%s/%d)" g d
      |  n, d -> Printf.sprintf "(%d*%s/%d)" n g d

    let lorentz_of_symbol model symbol =
      try
	SMap.find symbol model.lorentz
      with
      | Not_found -> invalid_arg ("lorentz_of_symbol: " ^ symbol)

    let lorentz_UFO_of_symbol model symbol =
      try
	SMap.find symbol model.lorentz_UFO
      with
      | Not_found -> invalid_arg ("lorentz_UFO_of_symbol: " ^ symbol)

    let coupling_of_symbol model symbol =
      try
	SMap.find symbol model.couplings
      with
      | Not_found -> invalid_arg ("coupling_of_symbol: " ^ symbol)

    let spin_triplet model name =
      match (lorentz_of_symbol model name).Lorentz.spins with
      | Lorentz.Unique [|s0; s1; s2|] -> (s0, s1, s2)
      | Lorentz.Unique _ -> invalid_arg "spin_triplet: wrong number of spins"
      | Lorentz.Unused -> invalid_arg "spin_triplet: Unused"
      | Lorentz.Ambiguous _ -> invalid_arg "spin_triplet: Ambiguous"
        
    let spin_quartet model name =
      match (lorentz_of_symbol model name).Lorentz.spins with
      | Lorentz.Unique [|s0; s1; s2; s3|] -> (s0, s1, s2, s3)
      | Lorentz.Unique _ -> invalid_arg "spin_quartet: wrong number of spins"
      | Lorentz.Unused -> invalid_arg "spin_quartet: Unused"
      | Lorentz.Ambiguous _ -> invalid_arg "spin_quartet: Ambiguous"
        
    let spin_multiplet model name =
      match (lorentz_of_symbol model name).Lorentz.spins with
      | Lorentz.Unique sarray -> sarray
      | Lorentz.Unused -> invalid_arg "spin_multiplet: Unused"
      | Lorentz.Ambiguous _ -> invalid_arg "spin_multiplet: Ambiguous"

    (* If we have reason to belive that a $\delta_{ab}$-vertex is
       an effective $\tr(T_aT_b)$-vertex generated at loop
       level, like~$gg\to H\ldots$ in the SM, we should interpret
       it as such and use the expression~(6.2) from~\cite{Kilian:2012pz}. *)

    (* AFAIK, there is no way to distinguish these cases directly
       in a UFO file.  Instead we rely in a heuristic, in which
       each massless color octet vector particle or ghost is a gluon
       and colorless scalars are potential Higgses. *)

    let is_massless p =
      match ThoString.uppercase p.Particle.mass with
      | "ZERO" -> true
      | _ -> false

    let is_gluon model f =
      let p = model.particle_array.(f) in
      match UFOx.Color.omega p.Particle.color,
            UFOx.Lorentz.omega p.Particle.spin with
      | Color.AdjSUN _, Coupling.Vector -> is_massless p
      | Color.AdjSUN _, Coupling.Scalar ->
         if p.Particle.ghost_number <> 0 then
           is_massless p
         else
           false
      | _ -> false

    let is_color_singlet model f =
      let p = model.particle_array.(f) in
      match UFOx.Color.omega p.Particle.color with
      | Color.Singlet -> true
      | _ -> false

    let is_higgs_gluon_vertex model p adjoints =
      if Array.length p > List.length adjoints then
        List.for_all
          (fun (i, p) ->
            if List.mem i adjoints then
              is_gluon model p
            else
              is_color_singlet model p)
          (ThoList.enumerate 1 (Array.to_list p))
      else
        false

    let delta8_heuristics model p a b =
      if is_higgs_gluon_vertex model p [a; b] then
        Color.Vertex.delta8_loop a b
      else
        Color.Vertex.delta8 a b

    let verbatim_higgs_glue = ref false

    let translate_color_atom model p = function
      | UFOx.Color_Atom.Identity (i, j) -> Color.Vertex.delta3 j i
      | UFOx.Color_Atom.Identity8 (a, b) ->
         if !verbatim_higgs_glue then
           Color.Vertex.delta8 a b
         else
           delta8_heuristics model p a b
      | UFOx.Color_Atom.T (a, i, j) -> Color.Vertex.t a i j
      | UFOx.Color_Atom.F (a, b, c) -> Color.Vertex.f a b c
      | UFOx.Color_Atom.D (a, b, c) -> Color.Vertex.d a b c
      | UFOx.Color_Atom.Epsilon (i, j, k) -> Color.Vertex.epsilon [i; j; k]
      | UFOx.Color_Atom.EpsilonBar (i, j, k) -> Color.Vertex.epsilon_bar [i; j; k]
      | UFOx.Color_Atom.T6 (a, i, j) -> Color.Vertex.t6 a i j
      | UFOx.Color_Atom.K6 (i, j, k) -> Color.Vertex.k6 i j k
      | UFOx.Color_Atom.K6Bar (i, j, k) -> Color.Vertex.k6bar i j k

    let translate_color_term model p = function
      | [], q ->
         Color.Vertex.scale q Color.Vertex.one
      | [atom], q ->
         Color.Vertex.scale q (translate_color_atom model p atom)
      | atoms, q ->
         let atoms = List.map (translate_color_atom model p) atoms in
         Color.Vertex.scale q (Color.Vertex.multiply atoms)

    let translate_color model p terms =
      match terms with
      | [] -> invalid_arg "translate_color: empty"
      | [ term ] -> translate_color_term model p term
      | terms ->
         Color.Vertex.sum (List.map (translate_color_term model p) terms)

    let translate_coupling_1 model p lcc =
      let l = lcc.Vertex.lorentz in
      let s = Array.to_list (spin_multiplet model l)
      and fl = (SMap.find l model.lorentz).Lorentz.fermion_lines
      and c = name (coupling_of_symbol model lcc.Vertex.coupling) in
      match lcc.Vertex.color with
      | UFOx.Color.Linear color ->
         let col = translate_color model p color in
         (Array.to_list p, Coupling.UFO (QC.unit, l, s, fl, col), c)
      | UFOx.Color.Ratios _ as color ->
         invalid_arg
           ("UFO.Model.translate_coupling: invalid color structure" ^
              UFOx.Color.to_string color)
        

    let translate_coupling model p lcc =
      List.map (translate_coupling_1 model p) lcc

    let long_flavors = ref false

    module type Lookup =
      sig
        type f = private
          { flavors : flavor list;
            flavor_of_string : string -> flavor;
            flavor_of_symbol : string -> flavor;
            particle : flavor -> Particle.t;
            flavor_symbol : flavor -> string;
            conjugate : flavor -> flavor }
        type flavor_format =
          | Long
          | Decimal
          | Hexadecimal
        val flavor_format : flavor_format ref
        val of_model : t -> f
      end

    module Lookup : Lookup =
      struct

        type f =
          { flavors : flavor list;
            flavor_of_string : string -> flavor;
            flavor_of_symbol : string -> flavor;
            particle : flavor -> Particle.t;
            flavor_symbol : flavor -> string;
            conjugate : flavor -> flavor }
            
        type flavor_format =
          | Long
          | Decimal
          | Hexadecimal

        let flavor_format = ref Hexadecimal

(*i
        let match_pdf_code p1 p2 =
	  p1.Particle.pdg_code = p2.Particle.pdg_code
i*)

        let conjugate_of_particle_array particles =
          Array.init
	    (Array.length particles)
	    (fun i ->
	      let f' = Particle.conjugate particles.(i) in
	      match ThoArray.match_all f' particles with
	      | [i'] -> i'
	      | [] ->
	         invalid_arg ("no charge conjugate: " ^ f'.Particle.name)
	      | _ ->
	         invalid_arg ("multiple charge conjugates: " ^ f'.Particle.name))

        let invert_flavor_array a =
          let table = SHash.create 37 in
          Array.iteri (fun i s -> SHash.add table s i) a;
          (fun name ->
	    try
	      SHash.find table name
	    with
	    | Not_found -> invalid_arg ("not found: " ^ name))

        let digits base n =
          let rec digits' acc n =
            if n < 1 then
              acc
            else
              digits' (succ acc) (n / base) in
          if n < 0 then
            digits' 1 (-n)
          else if n = 0 then
            1
          else
            digits' 0 n

        let of_model model =
          let particle_array = Array.of_list (values model.particles) in
          let conjugate_array = conjugate_of_particle_array particle_array
          and name_array = Array.map (fun f -> f.Particle.name) particle_array
          and symbol_array = Array.of_list (keys model.particles) in
          let flavor_symbol f =
            begin match !flavor_format with
            | Long -> symbol_array.(f)
            | Decimal -> 
               let w = digits 10 (Array.length particle_array - 1) in
               Printf.sprintf "%0*d" w f
            | Hexadecimal ->
               let w = digits 16 (Array.length particle_array - 1) in
               Printf.sprintf "%0*X" w f
            end in
          { flavors = ThoList.range 0 (Array.length particle_array - 1);
            flavor_of_string = invert_flavor_array name_array;
            flavor_of_symbol = invert_flavor_array symbol_array;
            particle = Array.get particle_array;
            flavor_symbol = flavor_symbol;
            conjugate = Array.get conjugate_array }

      end

    (* \begin{dubious}
         We appear to need to conjugate all flavors.  Why???
       \end{dubious} *)
    let translate_vertices model tables =
      let vn =
        List.fold_left
          (fun acc v ->
            let p = Array.map tables.Lookup.flavor_of_symbol v.Vertex.particles
            and lcc = v.Vertex.lcc in
            let p = Array.map conjugate p in (* FIXME: why? *)
            translate_coupling model p lcc @ acc)
          [] (values model.vertices) in
      ([], [], vn)

    let propagator_of_lorentz = function
      | Coupling.Scalar -> Coupling.Prop_Scalar
      | Coupling.Spinor -> Coupling.Prop_Spinor
      | Coupling.ConjSpinor -> Coupling.Prop_ConjSpinor
      | Coupling.Majorana -> Coupling.Prop_Majorana
      | Coupling.Maj_Ghost -> invalid_arg
         "UFO.Model.propagator_of_lorentz: SUSY ghosts do not propagate"
      | Coupling.Vector -> Coupling.Prop_Feynman
      | Coupling.Massive_Vector -> Coupling.Prop_Unitarity
      | Coupling.Tensor_2 -> Coupling.Prop_Tensor_2
      | Coupling.Vectorspinor -> invalid_arg
         "UFO.Model.propagator_of_lorentz: Vectorspinor"
      | Coupling.Tensor_1 -> invalid_arg
	 "UFO.Model.propagator_of_lorentz: Tensor_1"
      | Coupling.BRS _ -> invalid_arg
         "UFO.Model.propagator_of_lorentz: no BRST"

    let filter_unphysical model =
      let physical_particles =
	Particle.filter Particle.is_physical model.particles in
      let physical_particle_array =
        Array.of_list (values physical_particles) in
      let physical_vertices =
	Vertex.filter
	  (not <*> (Vertex.contains model.particles (not <*> Particle.is_physical)))
	  model.vertices in
      { model with
        particles = physical_particles;
        particle_array = physical_particle_array;
        vertices = physical_vertices }

    let whizard_constants =
      SSet.of_list
        [ "ZERO" ]

    let filter_constants parameters =
      List.filter
        (fun p ->
          not (SSet.mem (ThoString.uppercase p.Parameter.name) whizard_constants))
        parameters

    let add_name set parameter =
      CSet.add parameter.Parameter.name set

    let hardcoded_parameters =
      CSet.of_list
        ["cmath.pi"]

    let missing_parameters input derived couplings =
      let input_parameters =
        List.fold_left add_name hardcoded_parameters input in
      let all_parameters =
        List.fold_left add_name input_parameters derived in
      let derived_dependencies =
        dependencies
          (List.map
             (fun p -> (p.Parameter.name, p.Parameter.value))
             derived) in
      let coupling_dependencies =
        dependencies
          (List.map
             (fun p -> (p.UFO_Coupling.name, Expr p.UFO_Coupling.value))
             (values couplings)) in
      let missing_input =
        CMap.filter
          (fun parameter derived_parameters ->
            not (CSet.mem parameter all_parameters))
          derived_dependencies
      and missing =
        CMap.filter
          (fun parameter couplings ->
            not (CSet.mem parameter all_parameters))
          coupling_dependencies in
      CMap.iter
        (fun parameter derived_parameters ->
          Printf.eprintf
            "UFO warning: undefined input parameter %s appears in derived \
             parameters {%s}: will be added to the list of input parameters!\n"
            parameter (String.concat "; " (CSet.elements derived_parameters)))
        missing_input;
      CMap.iter
        (fun parameter couplings ->
          Printf.eprintf
            "UFO warning: undefined parameter %s appears in couplings {%s}: \
             will be added to the list of input parameters!\n"
            parameter (String.concat "; " (CSet.elements couplings)))
        missing;
      keys_caseless missing_input @ keys_caseless missing

    let classify_parameters model =
      let compare_parameters p1 p2 =
        compare p1.Parameter.sequence p2.Parameter.sequence in
      let input, derived =
        List.fold_left
          (fun (input, derived) p ->
            match p.Parameter.nature with
            | Parameter.Internal -> (input, p :: derived)
            | Parameter.External ->
               begin match p.Parameter.ptype with
               | Parameter.Real -> ()
               | Parameter.Complex ->
                  Printf.eprintf
                    "UFO warning: invalid complex declaration of input \
                     parameter `%s' ignored!\n"
                    p.Parameter.name
               end;
               (p :: input, derived))
          ([], []) (filter_constants (values model.parameters)) in
      let additional = missing_parameters input derived model.couplings in
      (List.sort compare_parameters input @ List.map Parameter.missing additional,
       List.sort compare_parameters derived)

(*i
      List.iter
        (fun line -> Printf.eprintf "par: %s\n" line)
        (dependencies_to_strings derived_dependencies);
      List.iter
        (fun line -> Printf.eprintf "coupling: %s\n" line)
        (dependencies_to_strings coupling_dependencies);
i*)

    let translate_name map name =
      try SMap.find name map with Not_found -> name

    let translate_input map p =
      (translate_name map p.Parameter.name, value_to_float p.Parameter.value)

    let alpha_s_half e =
      UFOx.Expr.substitute "aS" (UFOx.Expr.half "aS") e

    let alpha_s_half_etc map e =
      UFOx.Expr.rename (map_to_alist map) (alpha_s_half e)

    let translate_derived map p =
      let make_atom s = s in
      let c = make_atom (translate_name map p.Parameter.name)
      and v =
        value_to_coupling (alpha_s_half_etc map) make_atom p.Parameter.value in
      match p.Parameter.ptype with
      | Parameter.Real -> (Coupling.Real c, v)
      | Parameter.Complex -> (Coupling.Complex c, v)

    let translate_coupling_constant map c =
      let make_atom s = s in
      (Coupling.Complex c.UFO_Coupling.name,
       Coupling.Quot
         (value_to_coupling
            (alpha_s_half_etc map) make_atom
            (Expr c.UFO_Coupling.value),
          Coupling.I))

    module Lowercase_Parameters =
      struct
        type elt = string
        type base = string
        let compare_elt = compare
        let compare_base = compare
        let pi = ThoString.lowercase
      end

    module Lowercase_Bundle = Bundle.Make (Lowercase_Parameters)

    let coupling_names model =
      SMap.fold
        (fun _ c acc -> c.UFO_Coupling.name :: acc)
        model.couplings []

    let parameter_names model =
      SMap.fold
        (fun _ c acc -> c.Parameter.name :: acc)
        model.parameters []

    let ambiguous_parameters model =
      let all_names =
        List.rev_append (coupling_names model) (parameter_names model) in
      let lc_bundle = Lowercase_Bundle.of_list all_names in
      let lc_set =
        List.fold_left
          (fun acc s -> SSet.add s acc)
          SSet.empty (Lowercase_Bundle.base lc_bundle)
      and ambiguities =
        List.filter
          (fun (_, names) -> List.length names > 1)
          (Lowercase_Bundle.fibers lc_bundle) in
      (lc_set, ambiguities)

    let disambiguate1 lc_set name =
      let rec disambiguate1' i =
        let name' = Printf.sprintf "%s_%d" name i in
        let lc_name' = ThoString.lowercase name' in
        if SSet.mem lc_name' lc_set then
          disambiguate1' (succ i)
        else
          (SSet.add lc_name' lc_set, name') in
      disambiguate1' 1

    let disambiguate lc_set names =
      let _, replacements =
        List.fold_left
          (fun (lc_set', acc) name ->
            let lc_set'', name' = disambiguate1 lc_set' name in
            (lc_set'', SMap.add name name' acc))
          (lc_set, SMap.empty) names in
      replacements

    let omegalib_names =
      ["u"; "ubar"; "v"; "vbar"; "eps"]

    let translate_parameters model =
      let lc_set, ambiguities = ambiguous_parameters model in
      let replacements =
        disambiguate lc_set (ThoList.flatmap snd ambiguities) in
      SMap.iter
        (Printf.eprintf
           "warning: case sensitive parameter names: renaming '%s' -> '%s'\n")
        replacements;
      let replacements =
        List.fold_left
          (fun acc name -> SMap.add name ("UFO_" ^ name) acc)
          replacements omegalib_names in
      let input_parameters, derived_parameters = classify_parameters model
      and couplings = values model.couplings in
      { Coupling.input =
          List.map (translate_input replacements) input_parameters;
        Coupling.derived =
          List.map (translate_derived replacements) derived_parameters @
            List.map (translate_coupling_constant replacements) couplings;
        Coupling.derived_arrays = [] }

    (* UFO requires us to look up the mass parameter to
       distinguish between massless and massive vectors.

       TODO: this is a candidate for another lookup table. *)

    let lorentz_of_particle p =
      match UFOx.Lorentz.omega p.Particle.spin with
      | Coupling.Vector ->
         begin match ThoString.uppercase p.Particle.mass with
         | "ZERO" -> Coupling.Vector
         | _ -> Coupling.Massive_Vector
         end
      | s -> s

    type state =
      { directory : string;
        model : t }

    let initialized = ref None

    let is_initialized_from dir =
      match !initialized with
      | None -> false
      | Some state -> dir = state.directory

    let dump_raw = ref false

    let init dir =
      let model = filter_unphysical (parse_directory dir) in
      if !dump_raw then
	dump model;
      let tables = Lookup.of_model model in
      let vertices () = translate_vertices model tables in
      let particle f = tables.Lookup.particle f in
      let lorentz f = lorentz_of_particle (particle f) in
      let propagator f =
        let p = particle f in
        match p.Particle.propagator with
        | None -> propagator_of_lorentz (lorentz_of_particle p)
        | Some s -> Coupling.Prop_UFO s in
      let gauge_symbol () = "?GAUGE?" in
      let constant_symbol s = s in
      let parameters = translate_parameters model in
      M.setup
        ~color:(fun f -> UFOx.Color.omega (particle f).Particle.color)
        ~nc:(fun () -> model.nc)
        ~pdg:(fun f -> (particle f).Particle.pdg_code)
        ~lorentz
        ~propagator
        ~width:(fun f -> Coupling.Constant)
        ~goldstone:(fun f -> None)
        ~conjugate:tables.Lookup.conjugate
        ~fermion:(fun f -> fermion_of_lorentz (lorentz f))
        ~vertices
        ~flavors:[("All Flavors", tables.Lookup.flavors)]
        ~parameters:(fun () -> parameters)
        ~flavor_of_string:tables.Lookup.flavor_of_string
        ~flavor_to_string:(fun f -> (particle f).Particle.name)
        ~flavor_to_TeX:(fun f -> (particle f).Particle.texname)
        ~flavor_symbol:tables.Lookup.flavor_symbol
        ~gauge_symbol
        ~mass_symbol:(fun f -> (particle f).Particle.mass)
        ~width_symbol:(fun f -> (particle f).Particle.width)
        ~constant_symbol;
      initialized := Some { directory = dir; model = model }

    let ufo_directory = ref Config.default_UFO_dir

    let load () =
      if is_initialized_from !ufo_directory then
	()
      else
	init !ufo_directory

    let include_all_fusions = ref false

    (*   In case of Majorana spinors, also generate
         all combinations of charge conjugated fermion lines.
         The naming convention is to append
         \texttt{\_c}$nm$ if the $\gamma$-matrices
         of the fermion line $n\to m$ has been charge conjugated
         (this could become impractical for too many fermions at
         a vertex, but shouldn't matter in real life). *)

    (* Here we alway generate \emph{all} charge conjugations, because
       we treat \emph{all} fermions as Majorana fermion, if there
       is at least one Majorana fermion in the model! *)

    let is_majorana = function
      | Coupling.Majorana | Coupling.Vectorspinor | Coupling.Maj_Ghost -> true
      | _ -> false

    let name_spins_structure spins l =
      (l.Lorentz.name, spins, l.Lorentz.structure)

    let fusions_of_model ?only model =
      let include_fusion =
        match !include_all_fusions, only with
        | true, _
        | false, None -> (fun name -> true)
        | false, Some names -> (fun name -> SSet.mem name names)
      in
      SMap.fold
        (fun name l acc ->
          if include_fusion name then
            List.fold_left
              (fun acc p ->
                let l' = Lorentz.permute p l in
                match l'.Lorentz.spins with
                | Lorentz.Unused -> acc
                | Lorentz.Unique spins ->
                   if Array.exists is_majorana spins then
                     List.map
                       (name_spins_structure spins)
                       (Lorentz.required_charge_conjugates l')
                     @ acc
                   else
                     name_spins_structure spins l' :: acc
                | Lorentz.Ambiguous _ -> failwith "fusions: Lorentz.Ambiguous")
              [] (Permutation.Default.cyclic l.Lorentz.n) @ acc
          else
            acc)
        model.lorentz []

    let fusions ?only () =
      match !initialized with
      | None -> []
      | Some { model = model } -> fusions_of_model ?only model

    let propagators_of_model ?only model =
      let include_propagator =
        match !include_all_fusions, only with
        | true, _
        | false, None -> (fun name -> true)
        | false, Some names -> (fun name -> SSet.mem name names)
      in
      SMap.fold
        (fun name p acc ->
          if include_propagator name then
            (name, p) :: acc
          else
            acc)
        model.propagators []

    let propagators ?only () =
      match !initialized with
      | None -> []
      | Some { model = model } -> propagators_of_model ?only model

    let include_hadrons = ref true

    let ufo_majorana_warnings =
      [ "***************************************************";
        "*                                                 *";
        "* CAVEAT:                                         *";
        "*                                                 *";
        "*   These amplitudes have been computed for a     *";
        "*   UFO model containing Majorana fermions.       *";
        "*   This version of O'Mega contains some known    *";
        "*   bugs for this case.  It was released early at *";
        "*   the request of the Linear Collider community. *";
        "*                                                 *";
        "*   These amplitudes MUST NOT be used for         *";
        "*   publications without prior consulation        *";
        "*   with the WHIZARD authors !!!                  *";
        "*                                                 *";
        "***************************************************" ]

    let caveats () =
      if !use_majorana_spinors then
        ufo_majorana_warnings
      else
        []

    module Whizard : sig val write : unit -> unit end =
      struct
        
        let write_header dir =
          Printf.printf "# WHIZARD Model file derived from UFO directory\n";
          Printf.printf "#   '%s'\n\n" dir;
          List.iter (fun s -> Printf.printf "# %s\n" s) (M.caveats ());
          Printf.printf "model \"%s\"\n\n" (Filename.basename dir)

        let write_input_parameters parameters =
          let open Parameter in
          Printf.printf "# Independent (input) Parameters\n";
          List.iter
            (fun p ->
              Printf.printf
                "parameter %s = %s"
                p.name (value_to_numeric p.value);
              begin match p.lhablock, p.lhacode with
              | None, None -> ()
              | Some name, Some (index :: indices) ->
                 Printf.printf " slha_entry %s %d" name index;
                 List.iter (fun i -> Printf.printf " %d" i) indices
              | Some name, None ->
                 Printf.eprintf
                   "UFO: parameter %s: slhablock %s without slhacode\n"
                   p.name name
              | Some name, Some [] ->
                 Printf.eprintf
                   "UFO: parameter %s: slhablock %s with empty slhacode\n"
                   p.name name
              | None, Some _ ->
                 Printf.eprintf
                   "UFO: parameter %s: slhacode without slhablock\n"
                   p.name
              end;
              Printf.printf "\n")
            parameters;
          Printf.printf "\n"

        let write_derived_parameters parameters =
          let open Parameter in
          Printf.printf "# Dependent (derived) Parameters\n";
          List.iter
            (fun p ->
              Printf.printf
                "derived %s = %s\n"
                p.name (value_to_expr alpha_s_half p.value))
            parameters

        let write_particles particles =
          let open Particle in
          Printf.printf "# Particles\n";
          Printf.printf "# NB: hypercharge assignments appear to be unreliable\n";
          Printf.printf "#     therefore we can't infer the isospin\n";
          Printf.printf "# NB: parton-, gauge- & handedness are unavailable\n";
          List.iter
            (fun p ->
              if not p.is_anti then begin
                  Printf.printf
                    "particle \"%s\" %d ### parton? gauge? left?\n"
                    p.name p.pdg_code;
                  Printf.printf
                    "  spin %s charge %s color %s ### isospin?\n"
                    (UFOx.Lorentz.rep_to_string_whizard p.spin)
                    (charge_to_string p.charge)
                    (UFOx.Color.rep_to_string_whizard p.color);
                  Printf.printf "  name \"%s\"\n" p.name;
                  if p.antiname <> p.name then
                    Printf.printf "  anti \"%s\"\n" p.antiname;
                  Printf.printf "  tex_name \"%s\"\n" p.texname;
                  if p.antiname <> p.name then
                    Printf.printf "  tex_anti \"%s\"\n" p.antitexname;
                  Printf.printf "  mass %s width %s\n\n" p.mass p.width
                end)
            (values particles);
          Printf.printf "\n"

        let write_hadrons () =
          Printf.printf "# Hadrons (protons and beam remnants)\n";
          Printf.printf "# NB: these are NOT part of the UFO model\n";
          Printf.printf "#     but added for WHIZARD's convenience!\n";
          Printf.printf "particle PROTON 2212\n";
          Printf.printf "  spin 1/2  charge 1\n";
          Printf.printf "  name p \"p+\"\n";
          Printf.printf "  anti pbar \"p-\"\n";
          Printf.printf "particle HADRON_REMNANT 90\n";
          Printf.printf "  name hr\n";
          Printf.printf "  tex_name \"had_r\"\n";
          Printf.printf "particle HADRON_REMNANT_SINGLET 91\n";
          Printf.printf "  name hr1\n";
          Printf.printf "  tex_name \"had_r^{(1)}\"\n";
          Printf.printf "particle HADRON_REMNANT_TRIPLET 92\n";
          Printf.printf "  color 3\n";
          Printf.printf "  name hr3\n";
          Printf.printf "  tex_name \"had_r^{(3)}\"\n";
          Printf.printf "  anti hr3bar\n";
          Printf.printf "  tex_anti \"had_r^{(\\bar 3)}\"\n";
          Printf.printf "particle HADRON_REMNANT_OCTET 93\n";
          Printf.printf "  color 8\n";
          Printf.printf "  name hr8\n";
          Printf.printf "  tex_name \"had_r^{(8)}\"\n";
          Printf.printf "\n"

        let vertex_to_string model v =
          String.concat
            " "
            (List.map
               (fun s ->
                 "\"" ^ (SMap.find s model.particles).Particle.name ^ "\"")
               (Array.to_list v.Vertex.particles))

        let write_vertices3 model vertices  =
          Printf.printf "# Vertices (for phasespace generation only)\n";
          Printf.printf "# NB: particles should be sorted increasing in mass.\n";
          Printf.printf "#     This is NOT implemented yet!\n";
          List.iter
            (fun v ->
              if Array.length v.Vertex.particles = 3 then
                Printf.printf "vertex %s\n" (vertex_to_string model v))
            (values vertices);
          Printf.printf "\n"

        let write_vertices_higher model vertices  =
          Printf.printf
            "# Higher Order Vertices (ignored by phasespace generation)\n";
          List.iter
            (fun v ->
              if Array.length v.Vertex.particles <> 3 then
                Printf.printf "# vertex %s\n" (vertex_to_string model v))
            (values vertices);
          Printf.printf "\n"

        let write_vertices model vertices  =
          write_vertices3 model vertices;
          write_vertices_higher model vertices

        let write () =
          match !initialized with
          | None -> failwith "UFO.Whizard.write: UFO model not initialized"
          | Some { directory = dir; model = model } ->
             let input_parameters, derived_parameters =
               classify_parameters model in
             write_header dir;
             write_input_parameters input_parameters;
             write_derived_parameters derived_parameters;
             write_particles model.particles;
             if !include_hadrons then
               write_hadrons ();
             write_vertices model model.vertices;
             exit 0

      end

    let options =
      Options.create
        [ ("UFO_dir", Arg.String (fun name -> ufo_directory := name),
           "UFO model directory (default: " ^ !ufo_directory ^ ")");
          ("Majorana", Arg.Set use_majorana_spinors,
           "use Majorana spinors (must come _before_ exec!)");
          ("divide_propagators_by_i", Arg.Set divide_propagators_by_i,
           "divide propagators by I (pre 2013 FeynRules convention)");
          ("verbatim_Hg", Arg.Set verbatim_higgs_glue,
           "don't correct the color flows for effective Higgs Gluon couplings");
          ("write_WHIZARD", Arg.Unit Whizard.write,
           "write the WHIZARD model file (required once per model)");
          ("long_flavors",
           Arg.Unit (fun () -> Lookup.flavor_format := Lookup.Long),
           "write use the UFO flavor names instead of integers");
          ("dump", Arg.Set dump_raw,
           "dump UFO model for debugging the parser (must come _before_ exec!)");
          ("all_fusions", Arg.Set include_all_fusions,
           "include all fusions in the fortran module");
          ("no_hadrons", Arg.Clear include_hadrons,
           "don't add any particle not in the UFO file");
          ("add_hadrons", Arg.Set include_hadrons,
           "add protons and beam remants for WHIZARD");
          ("exec", Arg.Unit load,
           "load the UFO model files (required _before_ using particles names)");
          ("help", Arg.Unit (fun () -> prerr_endline "..."),
           "print information on the model")]

  end

module type Fortran_Target =
  sig

    val fuse :
      Algebra.QC.t -> string ->
      Coupling.lorentzn -> Coupling.fermion_lines ->
      string -> string list -> string list -> Coupling.fusen -> unit

    val lorentz_module :
      ?only:SSet.t -> ?name:string ->
      ?fortran_module:string -> ?parameter_module:string ->
      Format_Fortran.formatter -> unit -> unit

  end

module Targets =
  struct

    module Fortran : Fortran_Target =
      struct

        open Format_Fortran

        let fuse = UFO_targets.Fortran.fuse

        let lorentz_functions ff fusions () =
          List.iter
            (fun (name, s, l) ->
              UFO_targets.Fortran.lorentz ff name s l)
            fusions

        let propagator_functions ff parameter_module propagators () =
          List.iter
            (fun (name, p) ->
              UFO_targets.Fortran.propagator
                ff name
                parameter_module p.Propagator.variables
                p.Propagator.spins
                p.Propagator.numerator p.Propagator.denominator)
            propagators

        let lorentz_module
              ?only ?(name="omega_amplitude_ufo")
              ?(fortran_module="omega95")
              ?(parameter_module="parameter_module") ff () =
          let printf fmt = fprintf ff fmt
          and nl = pp_newline ff in
          printf "module %s" name; nl ();
          printf "  use kinds"; nl ();
          printf "  use %s" fortran_module; nl ();
          printf "  implicit none"; nl ();
          printf "  private"; nl ();
          let fusions = Model.fusions ?only ()
          and propagators = Model.propagators () in
          List.iter
            (fun (name, _, _) -> printf "  public :: %s" name; nl ())
            fusions;
          List.iter
            (fun (name, _) -> printf "  public :: pr_U_%s" name; nl ())
            propagators;
          UFO_targets.Fortran.eps4_g4_g44_decl ff ();
          UFO_targets.Fortran.eps4_g4_g44_init ff ();
          printf "contains"; nl ();
          UFO_targets.Fortran.inner_product_functions ff ();
          lorentz_functions ff fusions ();
          propagator_functions ff parameter_module propagators ();
          printf "end module %s" name; nl ();
          pp_flush ff ()

      end

  end

module type Test =
  sig
    val suite : OUnit.test
  end

module Test : Test =
  struct

    open OUnit

    let lexer s =
      UFO_lexer.token (UFO_lexer.init_position "" (Lexing.from_string s))

    let suite_lexer_escapes =
      "escapes" >:::

        [ "single-quote" >::
            (fun () ->
              assert_equal (UFO_parser.STRING "a'b'c") (lexer "'a\\'b\\'c'"));

          "unterminated" >::
            (fun () ->
              assert_raises End_of_file (fun () -> lexer "'a\\'b\\'c")) ]

    let suite_lexer =
      "lexer" >:::
        [suite_lexer_escapes]

    let suite =
      "UFO" >:::
        [suite_lexer]

  end
