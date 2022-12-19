(* vertex.ml --

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

(* Avoid refering to [Pervasives.compare], because [Pervasives] will
   become [Stdlib.Pervasives] in O'Caml 4.07 and [Stdlib] in O'Caml 4.08. *)
let pcompare = compare

module type Test =
  sig
    val example : unit -> unit
    val suite : OUnit.test
  end

(* \thocwmodulesection{New Implementation: Next Version} *)

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
  Vertex_syntax.File.expand_includes
    (fun file -> invalid_arg ("parse_string: found include `" ^ file ^ "'"))
    (try
       Vertex_parser.file
	 Vertex_lexer.token
	 (Vertex_lexer.init_position "" (Lexing.from_string text))
     with
     | Vertex_syntax.Syntax_Error (msg, start_pos, end_pos) ->
       invalid_arg (Printf.sprintf "syntax error (%s) at: `%s'"
                      msg  (error_in_string text start_pos end_pos))
     | Parsing.Parse_error ->
       invalid_arg ("parse error: " ^ text))

let parse_file name =
  let parse_file_tree name =
    let ic = open_in name in
    let file_tree =
      begin try
        Vertex_parser.file
	  Vertex_lexer.token
	  (Vertex_lexer.init_position name (Lexing.from_channel ic))
      with
      | Vertex_syntax.Syntax_Error (msg, start_pos, end_pos) ->
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
    file_tree in
  Vertex_syntax.File.expand_includes parse_file_tree (parse_file_tree name)

let dump_file pfx f =
  List.iter
    (fun s -> print_endline (pfx ^ ": " ^ s))
    (Vertex_syntax.File.to_strings f)

module Parser_Test : Test =
  struct

    let example () =
      ()

    open OUnit

    let compare s_out s_in () =
      assert_equal ~printer:(String.concat " ")
        [s_out] (Vertex_syntax.File.to_strings (parse_string s_in))

    let parse_error error s () =
      assert_raises (Invalid_argument error) (fun () -> parse_string s)

    let syntax_error (msg, error) s () =
      parse_error ("syntax error (" ^ msg ^ ") at: `" ^ error ^ "'") s ()

    let (=>) s_in s_out =
      " " ^ s_in >:: compare s_out s_in

    let (?>) s =
      s => s

    let (=>!!!) s error =
      " " ^ s >:: parse_error error s

    let (=>!) s error =
      " " ^ s >:: syntax_error error s

    let empty =
      "empty" >::
        (fun () -> assert_equal [] (parse_string ""))

    let expr =
      "expr" >:::
        [ "\\vertex[2 * (17 + 4)]{}" => "\\vertex[42]{{}}";
          "\\vertex[2 * 17 + 4]{}"   => "\\vertex[38]{{}}";
	  "\\vertex[2" =>! ("missing `]'", "[2");
	  "\\vertex]{}" =>! ("expected `[' or `{'", "\\vertex]");
	  "\\vertex2]{}" =>! ("expected `[' or `{'", "\\vertex2");
	  "\\vertex}{}" =>! ("expected `[' or `{'", "\\vertex}");
	  "\\vertex2}{}" =>! ("expected `[' or `{'", "\\vertex2");
	  "\\vertex[(2}{}" =>! ("expected `)', found `}'", "(2}");
	  "\\vertex[(2]{}" =>! ("expected `)', found `]'", "(2]");
	  "\\vertex{2]{}" =>! ("syntax error", "2");
	  "\\vertex[2}{}" =>! ("expected `]', found `}'", "[2}");
	  "\\vertex[2{}" =>! ("syntax error", "2");
	  "\\vertex[2*]{}" =>! ("syntax error", "2") ]

    let index =
      "index" >:::
        [ "\\vertex{{a}_{1}^{2}}" => "\\vertex{a^2_1}";
          "\\vertex{a_{11}^2}"    => "\\vertex{a^2_{11}}";
          "\\vertex{a_{1_1}^2}"   => "\\vertex{a^2_{1_1}}" ]

    let electron1 =
      "electron1" >:::
        [ ?> "\\charged{e^-}{e^+}";
          "\\charged{{e^-}}{{e^+}}" => "\\charged{e^-}{e^+}" ]

    let electron2 =
      "electron2" >:::
        [ "\\charged{e^-}{e^+}\\fortran{ele}" =>
          "\\charged{e^-}{e^+}\\fortran{{ele}}";
          "\\charged{e^-}{e^+}\\fortran{electron}\\fortran{ele}" =>
          "\\charged{e^-}{e^+}\\fortran{{ele}}\\fortran{{electron}}";
          "\\charged{e^-}{e^+}\\alias{e2}\\alias{e1}" =>
          "\\charged{e^-}{e^+}\\alias{{e1}}\\alias{{e2}}";
          "\\charged{e^-}{e^+}\\fortran{ele}\\anti\\fortran{pos}" =>
          "\\charged{e^-}{e^+}\\fortran{{ele}}\\anti\\fortran{{pos}}" ]

    let particles =
      "particles" >:::
        [electron1;
         electron2]

    let parameters =
      "parameters" >:::
        [ ?> "\\parameter{\\alpha}{1/137}";
          ?> "\\derived{\\alpha_s}{1/\\ln{\\frac{\\mu}{\\Lambda}}}";
          "\\parameter{\\alpha}{1/137}\\anti\\fortran{alpha}" =>!
          ("invalid parameter attribute", "\\anti") ]

    let indices =
      "indices" >:::
        [ ?> "\\index{a}\\color{8}";
          "\\index{a}\\color[SU(2)]{3}" => "\\index{a}\\color[{SU(2)}]{3}"  ]

    let tensors =
      "tensors" >:::
        [ "\\tensor{T}\\color{3}" => "\\tensor{T}\\color{3}"]

    let vertices =
      "vertex" >:::
        [ "\\vertex{\\bar\\psi\\gamma_\\mu\\psi A_\\mu}" =>
          "\\vertex{{{\\bar\\psi\\gamma_\\mu\\psi A_\\mu}}}" ]

    module T = Vertex_syntax.Token

    let parse_token s =
      match parse_string ("\\vertex{" ^ s ^ "}") with
      | [Vertex_syntax.File.Vertex (_, v)] -> v
      | _ -> invalid_arg "only_vertex"

    let print_token pfx t =
      print_endline (pfx ^ ": " ^ T.to_string t)

    let test_stem s_out s_in () =
      assert_equal ~printer:T.to_string
        (parse_token s_out)
        (T.stem (parse_token s_in))

    let (=>>) s_in s_out =
      "stem " ^ s_in >:: test_stem s_out s_in

    let tokens =
      "tokens" >:::
        [ "\\vertex{a'}" => "\\vertex{a^\\prime}";
          "\\vertex{a''}" => "\\vertex{a^{\\prime\\prime}}";
          "\\bar\\psi''_{i,\\alpha}" =>> "\\psi";
          "\\phi^\\dagger_{i'}" =>> "\\phi";
          "\\bar{\\phi\\psi}''_{i,\\alpha}" =>> "\\psi";
          "\\vertex{\\phi}" => "\\vertex{\\phi}";
          "\\vertex{\\phi_1}" => "\\vertex{\\phi_1}";
          "\\vertex{{{\\phi}'}}" => "\\vertex{\\phi^\\prime}";
          "\\vertex{\\hat{\\bar\\psi}_1}" => "\\vertex{\\hat\\bar\\psi_1}";
          "\\vertex{{a_b}_{cd}}" => "\\vertex{a_{bcd}}";
          "\\vertex{{\\phi_1}_2}" => "\\vertex{\\phi_{12}}";
          "\\vertex{{\\phi_{12}}_{34}}" => "\\vertex{\\phi_{1234}}";
          "\\vertex{{\\phi_{12}}^{34}}" => "\\vertex{\\phi^{34}_{12}}";
          "\\vertex{\\bar{\\psi_{\\mathrm{e}}}_\\alpha\\gamma_{\\alpha\\beta}^\\mu{\\psi_{\\mathrm{e}}}_\\beta}" =>
          "\\vertex{{{\\bar\\psi_{\\mathrm e\\alpha}\\gamma^\\mu_{\\alpha\\beta}\\psi_{\\mathrm e\\beta}}}}"]

    let suite =
      "Vertex_Parser" >:::
        [empty;
         index;
         expr;
         particles;
         parameters;
         indices;
         tensors;
         vertices;
         tokens ]

  end

(* \thocwmodulesubsection{Symbol Tables} *)

module type Symbol =
  sig

    type file = Vertex_syntax.File.t
    type t = Vertex_syntax.Token.t

    (* Tensors and their indices are representations of
       color, flavor or Lorentz groups.  In the end it might
       turn out to be unnecessary to distinguish [Color] from
       [Flavor].  *)
 
    type space =
    | Color of Vertex_syntax.Lie.t
    | Flavor of t list * t list
    | Lorentz of t list

    (* A symbol (i.\,e.~a [Symbol.t = Vertex_syntax.Token.t])
       can refer either to particles, to parameters (derived and input)
       or to tensors and indices.  *)
    type kind =
    | Neutral
    | Charged
    | Anti
    | Parameter
    | Derived
    | Index of space
    | Tensor of space

    type table
    val load : file -> table
    val dump : out_channel -> table -> unit

    (* Look up the [kind] of a symbol. *)
    val kind_of_symbol : table -> t -> kind option

    (* Look up the [kind] of a symbol's stem. *)
    val kind_of_stem : table -> t -> kind option

    (* Look up the [kind] of a symbol and fall back to the
       [kind] of the symbol's stem, if necessary. *)
    val kind_of_symbol_or_stem : table -> t -> kind option

    (* A table to look up all symbols with the same [stem]. *)
    val common_stem : table -> t -> t list

    exception Missing_Space of t
    exception Conflicting_Space of t

  end

module Symbol : Symbol =
  struct

    module T = Vertex_syntax.Token
    module F = Vertex_syntax.File
    module P = Vertex_syntax.Particle
    module I = Vertex_syntax.Index
    module L = Vertex_syntax.Lie
    module Q = Vertex_syntax.Parameter
    module X = Vertex_syntax.Tensor

    type file = F.t
    type t = T.t

    type space =
    | Color of L.t
    | Flavor of t list * t list
    | Lorentz of t list
        
    let space_to_string = function
      | Color (g, r) ->
	 "color:" ^ L.group_to_string g ^ ":" ^ L.rep_to_string r
      | Flavor (_, _) -> "flavor"
      | Lorentz _ -> "Lorentz"

    type kind =
    | Neutral
    | Charged
    | Anti
    | Parameter
    | Derived
    | Index of space
    | Tensor of space

    let kind_to_string = function
      | Neutral -> "neutral particle"
      | Charged -> "charged particle"
      | Anti -> "charged anti particle"
      | Parameter -> "input parameter"
      | Derived -> "derived parameter"
      | Index space -> space_to_string space ^ " index"
      | Tensor space -> space_to_string space ^ " tensor"

    module ST = Map.Make (T)
    module SS = Set.Make (T)

    type table = 
	{ symbol_kinds : kind ST.t;
	  stem_kinds : kind ST.t;
	  common_stems : SS.t ST.t }

    let empty = 
	{ symbol_kinds = ST.empty;
	  stem_kinds = ST.empty;
	  common_stems = ST.empty }

    let kind_of_symbol table token =
      try Some (ST.find token table.symbol_kinds) with Not_found -> None

    let kind_of_stem table token =
      try
	Some (ST.find (T.stem token) table.stem_kinds)
      with
      | Not_found -> None

    let kind_of_symbol_or_stem symbol_table token =
      match kind_of_symbol symbol_table token with
      | Some _ as kind -> kind
      | None -> kind_of_stem symbol_table token

    let common_stem table token =
      try
	SS.elements (ST.find (T.stem token) table.common_stems)
      with
      | Not_found -> []

    let add_symbol_kind table token kind =
      try
	let old_kind = ST.find token table in
	if kind = old_kind then
	  table
	else
	  invalid_arg ("conflicting symbol kind: " ^
			 T.to_string token ^ " -> " ^
			   kind_to_string kind ^ " vs " ^
			     kind_to_string old_kind)
      with
      | Not_found -> ST.add token kind table

    let add_stem_kind table token kind =
      let stem = T.stem token in
      try
	let old_kind = ST.find stem table in
	if kind = old_kind then
	  table
	else begin
	    match kind, old_kind with
	    | Charged, Anti -> ST.add stem Charged table
	    | Anti, Charged -> table
	    | _, _ ->
	       invalid_arg ("conflicting stem kind: " ^
			      T.to_string token ^ " -> " ^
				T.to_string stem ^ " -> " ^
				  kind_to_string kind ^ " vs " ^
				    kind_to_string old_kind)
	  end
      with
      | Not_found -> ST.add stem kind table

    let add_kind table token kind =
      { table with
	symbol_kinds = add_symbol_kind table.symbol_kinds token kind;
	stem_kinds = add_stem_kind table.stem_kinds token kind }

    let add_stem table token =
      let stem = T.stem token in
      let set =
	try
	  ST.find stem table.common_stems
	with
	| Not_found -> SS.empty in
      { table with
	common_stems = ST.add stem (SS.add token set) table.common_stems }

    (* Go through the list of attributes, make sure that
       the [space] is declared and unique.  Return the space. *)

    exception Missing_Space of t
    exception Conflicting_Space of t

    let group_rep_of_tokens group rep =
      let group =
	match group with
	| [] -> L.default_group
	| group -> L.group_of_string (T.list_to_string group) in
	Color (group, L.rep_of_string group (T.list_to_string rep))

    let index_space index =
      let spaces =
        List.fold_left
          (fun acc -> function
          | I.Color (group, rep) -> group_rep_of_tokens group rep :: acc
          | I.Flavor (group, rep) -> Flavor (rep, group) :: acc
          | I.Lorentz t -> Lorentz t :: acc)
          [] index.I.attr in
      match ThoList.uniq (List.sort compare spaces) with
      | [space] -> space
      | [] -> raise (Missing_Space index.I.name)
      | _ -> raise (Conflicting_Space index.I.name)

    let tensor_space tensor =
      let spaces =
        List.fold_left
          (fun acc -> function
          | X.Color (group, rep) -> group_rep_of_tokens rep group :: acc
          | X.Flavor (group, rep) -> Flavor (rep, group) :: acc
          | X.Lorentz t -> Lorentz t :: acc)
          [] tensor.X.attr in
      match ThoList.uniq (List.sort compare spaces) with
      | [space] -> space
      | [] -> raise (Missing_Space tensor.X.name)
      | _ -> raise (Conflicting_Space tensor.X.name)

    (* NB: if [P.Charged (name, name)] below, only
       the [Charged] will survive, [Anti] will be shadowed. *)
    let insert_kind table = function
      | F.Particle p ->
        begin match p.P.name with
        | P.Neutral name -> add_kind table name Neutral
        | P.Charged (name, anti) ->
          add_kind (add_kind table anti Anti) name Charged
        end
      | F.Index i -> add_kind table i.I.name (Index (index_space i))
      | F.Tensor t -> add_kind table t.X.name (Tensor (tensor_space t))
      | F.Parameter p ->
        begin match p with
        | Q.Parameter name -> add_kind table name.Q.name Parameter
        | Q.Derived name -> add_kind table name.Q.name Derived
        end
      | F.Vertex _ -> table

    let insert_stem table = function
      | F.Particle p ->
        begin match p.P.name with
        | P.Neutral name -> add_stem table name
        | P.Charged (name, anti) -> add_stem (add_stem table name) anti
        end
      | F.Index i -> add_stem table i.I.name
      | F.Tensor t -> add_stem table t.X.name
      | F.Parameter p ->
        begin match p with
        | Q.Parameter name
	| Q.Derived name -> add_stem table name.Q.name
        end
      | F.Vertex _ -> table

    let insert table token =
      insert_stem (insert_kind table token) token

    let load decls =
      List.fold_left insert empty decls

    let dump oc table =
      Printf.fprintf oc "<<< Symbol Table: >>>\n";
      ST.iter
	(fun s k ->
	 Printf.fprintf oc "%s -> %s\n" (T.to_string s) (kind_to_string k))
	table.symbol_kinds;
      Printf.fprintf oc "<<< Stem Table: >>>\n";
      ST.iter
	(fun s k ->
	 Printf.fprintf oc "%s -> %s\n" (T.to_string s) (kind_to_string k))
	table.stem_kinds;
      Printf.fprintf oc "<<< Common Stems: >>>\n";
      ST.iter
	(fun stem symbols ->
	 Printf.fprintf
	   oc "%s -> %s\n"
	   (T.to_string stem)
	   (String.concat
	      ", " (List.map T.to_string (SS.elements symbols))))
	table.common_stems

  end

(* \thocwmodulesubsection{Declarations} *)

module type Declaration =
  sig

    type t

    val of_string : string -> t list
    val to_string : t list -> string

    (* For testing and debugging *)
    val of_string_and_back : string -> string

    val count_indices : t -> (int * Symbol.t) list
    val indices_ok : t -> unit

  end

module Declaration : Declaration =
  struct

    module S = Symbol
    module T = Vertex_syntax.Token

    type factor =
      { stem : T.t;
        prefix : T.prefix list;
        particle : T.t list;
        color : T.t list;
        flavor : T.t list;
        lorentz : T.t list;
        other : T.t list }

    type t = factor list

    let factor_stem token =
      { stem = token.T.stem;
        prefix = token.T.prefix;
        particle = [];
        color = [];
        flavor = [];
        lorentz = [];
        other = [] }

    let rev factor =
      { stem = factor.stem;
        prefix = List.rev factor.prefix;
        particle = List.rev factor.particle;
        color = List.rev factor.color;
        flavor = List.rev factor.flavor;
        lorentz = List.rev factor.lorentz;
        other = List.rev factor.other }

    let factor_add_prefix factor token =
      { factor with prefix = T.prefix_of_string token :: factor.prefix }

    let factor_add_particle factor token =
      { factor with particle = token :: factor.particle }

    let factor_add_color_index t factor token =
      { factor with color = token :: factor.color }

    let factor_add_lorentz_index t factor token =
      (* diagnostics: [Printf.eprintf "[L:[%s]]\n" (T.to_string token);] *)
      { factor with lorentz = token :: factor.lorentz }

    let factor_add_flavor_index t factor token =
      { factor with flavor = token :: factor.flavor }

    let factor_add_other_index factor token =
      { factor with other = token :: factor.other }

    let factor_add_kind factor token = function
      | S.Neutral | S.Charged | S.Anti -> factor_add_particle factor token
      | S.Index (S.Color (rep, group)) ->
	 factor_add_color_index (rep, group) factor token
      | S.Index (S.Flavor (rep, group)) ->
	 factor_add_flavor_index (rep, group) factor token
      | S.Index (S.Lorentz t) -> factor_add_lorentz_index t factor token
      | S.Tensor _ -> invalid_arg "factor_add_index: \\tensor"
      | S.Parameter -> invalid_arg "factor_add_index: \\parameter"
      | S.Derived -> invalid_arg "factor_add_index: \\derived"

    let factor_add_index symbol_table factor = function
      | T.Token "," -> factor
      | T.Token ("*" | "\\ast" as star) -> factor_add_prefix factor star
      | token ->
         begin
	   match S.kind_of_symbol_or_stem symbol_table token with
           | Some kind -> factor_add_kind factor token kind
           | None -> factor_add_other_index factor token
	 end

    let factor_of_token symbol_table token =
      let token = T.wrap_scripted token in
      rev (List.fold_left
             (factor_add_index symbol_table)
             (factor_stem token)
             (token.T.super @ token.T.sub))

    let list_to_string tag = function
      | [] -> ""
      | l -> "; " ^ tag ^ "=" ^ String.concat "," (List.map T.to_string l)

    let factor_to_string factor =
       "[" ^ T.to_string factor.stem ^
         (match factor.prefix with
         | [] -> ""
         | l -> "; prefix=" ^
		  String.concat "," (List.map T.prefix_to_string l)) ^
           list_to_string "particle" factor.particle ^
           list_to_string "color" factor.color ^
           list_to_string "flavor" factor.flavor ^
           list_to_string "lorentz" factor.lorentz ^
           list_to_string "other" factor.other ^ "]"

    let count_indices factors =
      ThoList.classify
	(ThoList.flatmap (fun f -> f.color @ f.flavor @ f.lorentz) factors)

    let format_mismatch (n, index) =
      Printf.sprintf "index %s appears %d times" (T.to_string index) n

    let indices_ok factors =
      match List.filter (fun (n, _) -> n <> 2) (count_indices factors) with
      | [] -> ()
      | mismatches ->
	 invalid_arg (String.concat ", " (List.map format_mismatch mismatches))
      
    let of_string s =
      let decls = parse_string s in
      let symbol_table = Symbol.load decls in
      (* diagnostics: [Symbol.dump stderr symbol_table;] *)
      let tokens =
        List.fold_left
          (fun acc -> function
          | Vertex_syntax.File.Vertex (_, v) -> T.wrap_list v :: acc
          | _ -> acc)
          [] decls in
      let vlist = List.map (List.map (factor_of_token symbol_table)) tokens in
      List.iter indices_ok vlist;
      vlist

    let to_string decls =
      String.concat "; "
        (List.map
	   (fun v -> String.concat " * " (List.map factor_to_string v))
	   decls)

    let of_string_and_back s =
      to_string (of_string s)

    type field =
      { name : T.t list }

  end

(* \thocwmodulesubsection{Complete Models} *)

module Modelfile =
  struct

  end

module Modelfile_Test =
  struct

    let example () =
      ()

    open OUnit

    let index_mismatches =
      "index mismatches" >:::
	[ "1" >::
	    (fun () ->
	     assert_raises
	       (Invalid_argument "index a_1 appears 1 times, \
				  index a_2 appears 1 times")
	       (fun () -> Declaration.of_string_and_back
			    "\\index{a}\\color{3}\
			     \\vertex{\\bar\\psi_{a_1}\\psi_{a_2}}"));
	  "3" >::
	    (fun () ->
	     assert_raises
	       (Invalid_argument "index a appears 3 times")
	       (fun () -> Declaration.of_string_and_back
			    "\\index{a}\\color{3}\
			     \\vertex{\\bar\\psi_a\\psi_a\\phi_a}")) ]

    let kind_conflicts =
      "kind conflictings" >:::
	[ "lorentz / color" >::
	    (fun () ->
	     assert_raises
	       (Invalid_argument
		  "conflicting stem kind: a_2 -> a -> \
		   Lorentz index vs color:SU(3):3 index")
	       (fun () -> Declaration.of_string_and_back
			    "\\index{a_1}\\color{3}\
			     \\index{a_2}\\lorentz{X}"));
	  "color / color" >::
	    (fun () ->
	     assert_raises
	       (Invalid_argument
		  "conflicting stem kind: a_2 -> a -> \
		   color:SU(3):8 index vs color:SU(3):3 index")
	       (fun () -> Declaration.of_string_and_back
			    "\\index{a_1}\\color{3}\
			     \\index{a_2}\\color{8}"));
	  "neutral / charged" >::
	    (fun () ->
	     assert_raises
	       (Invalid_argument
		  "conflicting stem kind: H^- -> H -> \
		   charged anti particle vs neutral particle")
	       (fun () -> Declaration.of_string_and_back
			    "\\neutral{H}\
			     \\charged{H^+}{H^-}")) ]

    let suite =
      "Modelfile_Test" >:::
        [ "ok" >::
            (fun () ->
              assert_equal ~printer:(fun s -> s)
                "[\\psi; prefix=\\bar; \
                  particle=e; color=a; lorentz=\\alpha_1] * \
                 [\\gamma; lorentz=\\mu,\\alpha_1,\\alpha_2] * \
                 [\\psi; particle=e; color=a; lorentz=\\alpha_2] * \
                 [A; lorentz=\\mu]"
                (Declaration.of_string_and_back
                   "\\charged{e^-}{e^+}\
                    \\index{a}\\color{\\bar3}\
                    \\index{b}\\color[SU(3)]{8}\
                    \\index{\\mu}\\lorentz{X}\
                    \\index{\\alpha}\\lorentz{X}\
                    \\vertex{\\bar{\\psi_e}_{a,\\alpha_1}\
                             \\gamma^\\mu_{\\alpha_1\\alpha_2}\
                             {\\psi_e}_{a,\\alpha_2}A_\\mu}"));
	  index_mismatches;
	  kind_conflicts;
          "QCD.omf" >::
            (fun () ->
              dump_file "QCD" (parse_file "QCD.omf"));
          "SM.omf" >::
            (fun () ->
              dump_file "SM" (parse_file "SM.omf"));
          "SM-error.omf" >::
            (fun () ->
	     assert_raises
	       (Invalid_argument
		  "SM-error.omf:32.22-32.27: syntax error (syntax error)")
	       (fun () -> parse_file "SM-error.omf"));
          "cyclic.omf" >::
            (fun () ->
	     assert_raises
	       (Invalid_argument "cyclic \\include{cyclic.omf}")
	       (fun () -> parse_file "cyclic.omf")) ]

  end

(* \thocwmodulesection{New Implementation: Obsolete Version~1} *)

(* Start of version 1 of the new implementation.  The old syntax
   will not be used in the real implementation, but the library
   for dealing with indices and permutations will remail important. *)

(* Note that [arity = length lorentz_reps = length color_reps].  Do we
   need to enforce this by an abstract type constructor?

   A cleaner approach would be
   [type context = (Coupling.lorentz, Color.t) array], but it would also
   require more tedious deconstruction of the pairs.  Well, an abstract
   type with accessors might be the way to go after all \ldots *)

type context =
    { arity : int;
      lorentz_reps : Coupling.lorentz array;
      color_reps : Color.t array }

let distinct2 i j =
  i <> j

let distinct3 i j k =
  i <> j && j <> k && k <> i

let distinct ilist =
  List.length (ThoList.uniq (List.sort compare ilist)) =
  List.length ilist

(* An abstract type that allows us to distinguish offsets
   in the field array from color and Lorentz indices in
   different representations. *)

module type Index =
  sig
    type t
    val of_int : int -> t
    val to_int : t -> int
  end

(* While the number of allowed indices is unlimited, the
   allowed offsets into the field arrays are of course
   restricted to the fields in the current [context]. *)

module type Field =
  sig
    type t
    exception Out_of_range of int
    val of_int : context -> int -> t
    val to_int : t -> int
    val get : 'a array -> t -> 'a
  end

module Field : Field =
  struct
    type t = int
    exception Out_of_range of int
    let of_int context i =
      if 0 <= i && i < context.arity then
        i
      else
        raise (Out_of_range i)
    let to_int i = 0
    let get = Array.get
  end

type field = Field.t

module type Lorentz =
  sig

    (* We combine indices~[I] and offsets~[F] into the field array
       into a single type so that we can unify vectors with vector
       components.  *)

    type index = I of int | F of field

    type vector = Vector of index

    type spinor = Spinor of index

    type conjspinor = ConjSpinor of index

    (* These are all the primitive ways to construct Lorentz tensors,
       a.\,k.\,a.~objects with Lorentz indices, from momenta, other
       Lorentz tensors and Dirac spinors: *)

    type primitive =
      | G of vector * vector                       (* $g_{\mu_1\mu_2}$ *)
      | E of vector * vector * vector * vector     (* $\epsilon_{\mu_1\mu_2\mu_3\mu_4}$ *)
      | K of vector * field                        (* $k_{2}^{\mu_1}$ *)
      | S of conjspinor * spinor                   (* $\bar\psi_1\psi_2$ *)
      | V of vector * conjspinor * spinor          (* $\bar\psi_1\gamma_{\mu_2}\psi_3$ *)
      | T of vector * vector * conjspinor * spinor (* $\bar\psi_1\sigma_{\mu_2\mu_3}\psi_4$ *)
      | A of vector * conjspinor * spinor          (* $\bar\psi_1\gamma_{\mu_2}\gamma_5\psi_3$ *)
      | P of conjspinor * spinor                   (* $\bar\psi_1\gamma_5\psi_2$ *)

    type tensor = int * primitive list

(* Below, we will need to permute fields.  For this purpose, we
   introduce the function
   [map_primitive v_idx v_fld s_idx s_fld c_idx c_fld tensor]
   that returns a structurally identical tensor, with
   [v_idx : int -> int] applied to all vector indices,
   [v_fld : field -> field] to all vector fields,
   [s_idx] and [c_idx] to all (conj)spinor indices and
   [s_fld] and [c_fld] to all (conj)spinor fields.

   Note we must treat spinors and vectors differently,
   even for simple permuations, in order to handle the
   statistics properly.  *)

    val map_tensor :
      (int -> int) -> (field -> field) -> (int -> int) -> (field -> field) ->
      (int -> int) -> (field -> field) -> tensor -> tensor

(* Check whether the [tensor] is well formed in the [context]. *)

    val tensor_ok : context -> tensor -> bool

(* The lattice $\mathbf{N}+\mathrm{i}\mathbf{N}\subset\mathbf{C}$, which
   suffices for representing the matrix elements of Dirac matrices.
   We hope to be able to avoid the lattice
   $\mathbf{Q}+\mathrm{i}\mathbf{Q}\subset\mathbf{C}$ or
   $\mathbf{C}$ itself down the road. *)

    module Complex :
      sig
        type t = int * int
        type t' =
          | Z (* $0$ *)
          | O (* $1$ *)
          | M (* $-1$ *)
          | I (* $\mathrm{i}$ *)
          | J (* $-\mathrm{i}$ *)
          | C of int * int (* $x+\mathrm{i}y$ *)
        val to_fortran : t' -> string
      end

    (* Sparse Dirac matrices as maps from Lorentz and Spinor indices
       to complex numbers.  This is supposed to be independent of
       the representation. *)

    module type Dirac =
      sig
        val scalar : int -> int -> Complex.t'
        val vector : int -> int -> int -> Complex.t'
        val tensor : int -> int -> int -> int -> Complex.t'
        val axial : int -> int -> int -> Complex.t'
        val pseudo : int -> int -> Complex.t'
      end

    (* Dirac matrices as tables of nonzero entries.  There will
       be one concrete Module per realization. *)

    module type Dirac_Matrices =
      sig
        type t = (int * int * Complex.t') list
        val scalar : t
        val vector : (int * t) list
        val tensor : (int * int * t) list
        val axial : (int * t) list
        val pseudo : t
      end

    (* E.\,g.~the chiral representation: *)

    module Chiral : Dirac_Matrices

    (* Here's the functor to create the maps corresponding to
       a given realization. *)

    module Dirac : functor (M : Dirac_Matrices) -> Dirac

  end

module Lorentz : Lorentz =
  struct

    
    type index =
      | I of int (* $\mu_0,\mu_1,\ldots$, not $0,1,2,3$ *)
      | F of field

    let map_index fi ff = function
      | I i -> I (fi i)
      | F i -> F (ff i)

    let indices = function
      | I i -> [i]
      | F _ -> []

    (* Is the following level of type checks useful or redundant? *)

    (* TODO: should we also support a [tensor] like $F_{\mu_1\mu_2}$? *)

    type vector = Vector of index
    type spinor = Spinor of index
    type conjspinor = ConjSpinor of index

    let map_vector fi ff (Vector i) = Vector (map_index fi ff i)
    let map_spinor fi ff (Spinor i) = Spinor (map_index fi ff i)
    let map_conjspinor fi ff (ConjSpinor i) = ConjSpinor (map_index fi ff i)

    let vector_ok context = function
      | Vector (I _) ->
        (* we could perform additional checks! *)
        true
      | Vector (F i) ->
          begin
            match Field.get context.lorentz_reps i with
            | Coupling.Vector -> true
            | Coupling.Vectorspinor ->
                failwith "Lorentz.vector_ok: incomplete"
            | _ -> false
          end
      
    let spinor_ok context = function
      | Spinor (I _) ->
        (* we could perfrom additional checks! *)
        true
      | Spinor (F i) ->
          begin
            match Field.get context.lorentz_reps i with
            | Coupling.Spinor -> true
            | Coupling.Vectorspinor | Coupling.Majorana ->
                failwith "Lorentz.spinor_ok: incomplete"
            | _ -> false
          end

    let conjspinor_ok context = function
      | ConjSpinor (I _) ->
        (* we could perform additional checks! *)
        true
      | ConjSpinor (F i) ->
          begin
            match Field.get context.lorentz_reps i with
            | Coupling.ConjSpinor -> true
            | Coupling.Vectorspinor | Coupling.Majorana ->
                failwith "Lorentz.conjspinor_ok: incomplete"
            | _ -> false
          end

    (* Note that [distinct2 i j] is automatically guaranteed
       for Dirac spinors, because the $\bar\psi$ and $\psi$ can
       not appear in the same slot.  This is however not the
       case for Weyl and Majorana spinors. *)

    let spinor_sandwitch_ok context i j =
      conjspinor_ok context i && spinor_ok context j

    type primitive =
      | G of vector * vector
      | E of vector * vector * vector * vector
      | K of vector * field
      | S of conjspinor * spinor
      | V of vector * conjspinor * spinor
      | T of vector * vector * conjspinor * spinor
      | A of vector * conjspinor * spinor
      | P of conjspinor * spinor

    let map_primitive fvi fvf fsi fsf fci fcf = function
      | G (mu, nu) ->
          G (map_vector fvi fvf mu, map_vector fvi fvf nu)
      | E (mu, nu, rho, sigma) ->
          E (map_vector fvi fvf mu,
             map_vector fvi fvf nu,
             map_vector fvi fvf rho,
             map_vector fvi fvf sigma)
      | K (mu, i) ->
          K (map_vector fvi fvf mu, fvf i)
      | S (i, j) ->
          S (map_conjspinor fci fcf i, map_spinor fsi fsf j)
      | V (mu, i, j) ->
          V (map_vector fvi fvf mu,
             map_conjspinor fci fcf i,
             map_spinor fsi fsf j)
      | T (mu, nu, i, j) ->
          T (map_vector fvi fvf mu,
             map_vector fvi fvf nu,
             map_conjspinor fci fcf i,
             map_spinor fsi fsf j)
      | A (mu, i, j) ->
          A (map_vector fvi fvf mu,
             map_conjspinor fci fcf i,
             map_spinor fsi fsf j)
      | P (i, j) ->
          P (map_conjspinor fci fcf i, map_spinor fsi fsf j)

    let primitive_ok context =
      function
        | G (mu, nu) ->
            distinct2 mu nu &&
            vector_ok context mu && vector_ok context nu
        | E (mu, nu, rho, sigma) ->
            let i = [mu; nu; rho; sigma] in
            distinct i && List.for_all (vector_ok context) i
        | K (mu, i) ->
            vector_ok context mu
        | S (i, j) | P (i, j) ->
            spinor_sandwitch_ok context i j
        | V (mu, i, j) | A (mu, i, j) ->
            vector_ok context mu && spinor_sandwitch_ok context i j
        | T (mu, nu, i, j) ->
            vector_ok context mu && vector_ok context nu &&
            spinor_sandwitch_ok context i j

    let primitive_vector_indices = function
      | G (Vector mu, Vector nu) | T (Vector mu, Vector nu, _, _) ->
          indices mu @ indices nu
      | E (Vector mu, Vector nu, Vector rho, Vector sigma) ->
          indices mu @ indices nu @ indices rho @ indices sigma
      | K (Vector mu, _)
      | V (Vector mu, _, _)
      | A (Vector mu, _, _) -> indices mu
      | S (_, _) | P (_, _) -> []

    let vector_indices p =
      ThoList.flatmap primitive_vector_indices p

    let primitive_spinor_indices = function
      | G (_, _) | E (_, _, _, _) | K (_, _) -> []
      | S (_, Spinor alpha) | V (_, _, Spinor alpha)
      | T (_, _, _, Spinor alpha)
      | A (_, _, Spinor alpha) | P (_, Spinor alpha) -> indices alpha

    let spinor_indices p =
      ThoList.flatmap primitive_spinor_indices p

    let primitive_conjspinor_indices = function
      | G (_, _) | E (_, _, _, _) | K (_, _) -> []
      | S (ConjSpinor alpha, _) | V (_, ConjSpinor alpha, _)
      | T (_, _, ConjSpinor alpha, _)
      | A (_, ConjSpinor alpha, _) | P (ConjSpinor alpha, _) -> indices alpha

    let conjspinor_indices p =
      ThoList.flatmap primitive_conjspinor_indices p

    let vector_contraction_ok p =
      let c = ThoList.classify (vector_indices p) in
      print_endline
        (String.concat ", "
           (List.map
              (fun (n, i) -> string_of_int n ^ " * " ^ string_of_int i)
              c));
      flush stdout;
      let res = List.for_all (fun (n, _) -> n = 2) c in
      res

    let two_of_each indices p =
      List.for_all (fun (n, _) -> n = 2) (ThoList.classify (indices p))

    let vector_contraction_ok = two_of_each vector_indices
    let spinor_contraction_ok = two_of_each spinor_indices
    let conjspinor_contraction_ok = two_of_each conjspinor_indices

    let contraction_ok p =
      vector_contraction_ok p &&
      spinor_contraction_ok p && conjspinor_contraction_ok p

    type tensor = int * primitive list

    let map_tensor fvi fvf fsi fsf fci fcf (factor, primitives) =
      (factor, List.map (map_primitive fvi fvf fsi fsf fci fcf ) primitives)

    let tensor_ok context (_, primitives) =
      List.for_all (primitive_ok context) primitives &&
      contraction_ok primitives

    module Complex =
      struct

        type t = int * int

        type t' = Z | O | M | I | J | C of int * int

        let to_fortran = function
          | Z -> "(0,0)"
          | O -> "(1,0)"
          | M -> "(-1,0)"
          | I -> "(0,1)"
          | J -> "(0,-1)"
          | C (r, i) -> "(" ^ string_of_int r ^ "," ^ string_of_int i ^ ")"

      end

    module type Dirac =
      sig
        val scalar : int -> int -> Complex.t'
        val vector : int -> int -> int -> Complex.t'
        val tensor : int -> int -> int -> int -> Complex.t'
        val axial : int -> int -> int -> Complex.t'
        val pseudo : int -> int -> Complex.t'
      end

    module type Dirac_Matrices =
      sig
        type t = (int * int * Complex.t') list
        val scalar : t
        val vector : (int * t) list
        val tensor : (int * int * t) list
        val axial : (int * t) list
        val pseudo : t
      end

    module Chiral : Dirac_Matrices =
      struct

        type t = (int * int * Complex.t') list

        let scalar =
          [ (1, 1, Complex.O);
            (2, 2, Complex.O);
            (3, 3, Complex.O);
            (4, 4, Complex.O) ]

        let vector =
          [ (0, [ (1, 4, Complex.O);
                  (4, 1, Complex.O);
                  (2, 3, Complex.M);
                  (3, 2, Complex.M) ]);
            (1, [ (1, 3, Complex.O);
                  (3, 1, Complex.O);
                  (2, 4, Complex.M);
                  (4, 2, Complex.M) ]);
            (2, [ (1, 3, Complex.I);
                  (3, 1, Complex.I);
                  (2, 4, Complex.I);
                  (4, 2, Complex.I) ]);
            (3, [ (1, 4, Complex.M);
                  (4, 1, Complex.M);
                  (2, 3, Complex.M);
                  (3, 2, Complex.M) ]) ]

        let tensor =
          [ (* TODO!!! *) ]

        let axial =
          [ (0, [ (1, 4, Complex.M);
                  (4, 1, Complex.O);
                  (2, 3, Complex.O);
                  (3, 2, Complex.M) ]);
            (1, [ (1, 3, Complex.M);
                  (3, 1, Complex.O);
                  (2, 4, Complex.O);
                  (4, 2, Complex.M) ]);
            (2, [ (1, 3, Complex.J);
                  (3, 1, Complex.I);
                  (2, 4, Complex.J);
                  (4, 2, Complex.I) ]);
            (3, [ (1, 4, Complex.O);
                  (4, 1, Complex.M);
                  (2, 3, Complex.O);
                  (3, 2, Complex.M) ]) ]

        let pseudo =
          [ (1, 1, Complex.M);
            (2, 2, Complex.M);
            (3, 3, Complex.O);
            (4, 4, Complex.O) ]

      end

    module Dirac (M : Dirac_Matrices) : Dirac =
      struct

        module Map2 =
          Map.Make
            (struct
              type t = int * int
              let compare = pcompare
            end)
            
        let init2 triples =
          List.fold_left
            (fun acc (i, j, e) -> Map2.add (i, j) e acc)
            Map2.empty triples

        let bounds_check2 i j =
          if i < 1 || i > 4 || j < 0 || j > 4 then
            invalid_arg "Chiral.bounds_check2"

        let lookup2 map i j =
          bounds_check2 i j;
          try Map2.find (i, j) map with Not_found -> Complex.Z

        module Map3 =
          Map.Make
            (struct
              type t = int * (int * int)
              let compare = pcompare
            end)
            
        let init3 quadruples =
          List.fold_left
            (fun acc (mu, gamma) ->
             List.fold_right
               (fun (i, j, e) -> Map3.add (mu, (i, j)) e)
               gamma acc)
            Map3.empty quadruples

        let bounds_check3 mu i j =
          bounds_check2 i j;
          if mu < 0 || mu > 3 then
            invalid_arg "Chiral.bounds_check3"

        let lookup3 map mu i j =
          bounds_check3 mu i j;
          try Map3.find (mu, (i, j)) map with Not_found -> Complex.Z

        module Map4 =
          Map.Make
            (struct
              type t = int * int * (int * int)
              let compare = pcompare
            end)
            
        let init4 quadruples =
          List.fold_left
            (fun acc (mu, nu, gamma) ->
             List.fold_right
               (fun (i, j, e) -> Map4.add (mu, nu, (i, j)) e)
               gamma acc)
            Map4.empty quadruples

        let bounds_check4 mu nu i j =
          bounds_check3 nu i j;
          if mu < 0 || mu > 3 then
            invalid_arg "Chiral.bounds_check4"

        let lookup4 map mu nu i j =
          bounds_check4 mu nu i j;
          try Map4.find (mu, nu, (i, j)) map with Not_found -> Complex.Z

        let scalar_map = init2 M.scalar
        let vector_map = init3 M.vector
        let tensor_map = init4 M.tensor
        let axial_map = init3 M.axial
        let pseudo_map = init2 M.pseudo

        let scalar = lookup2 scalar_map
        let vector = lookup3 vector_map
        let tensor mu nu i j =
          lookup4 tensor_map mu nu i j
        let tensor mu nu i j =
          failwith "tensor: incomplete"
        let axial = lookup3 axial_map
        let pseudo = lookup2 pseudo_map

      end

  end

module type Color =
  sig
    module Index : Index
    type index = Index.t
    type color_rep = F of field | C of field | A of field
    type primitive =
      | D of field * field
      | E of field * field * field  (* only for $SU(3)$ *)
      | T of field * field * field
      | F of field * field * field
    val map_primitive : (field -> field) -> primitive -> primitive
    val primitive_indices : primitive -> field list
    val indices : primitive list -> field list
    type tensor = int * primitive list
    val map_tensor :
      (field -> field) -> 'a * primitive list -> 'a * primitive list
    val tensor_ok : context -> 'a * primitive list -> bool
  end

module Color : Color = 
  struct

    module Index : Index =
      struct
        type t = int
        let of_int i = i
        let to_int i = i
      end

    (* $a_0,a_1,\ldots$, not $0,1,\ldots$ *)
    type index = Index.t

    type color_rep =
      | F of field
      | C of field
      | A of field

    type primitive =
      | D of field * field
      | E of field * field * field
      | T of field * field * field
      | F of field * field * field

    let map_primitive f = function
      | D (i, j) -> D (f i, f j)
      | E (i, j, k) -> E (f i, f j, f k)
      | T (a, i, j) -> T (f a, f i, f j)
      | F (a, b, c) -> F (f a, f b, f c)

    let primitive_ok ctx =
      function
        | D (i, j) ->
            distinct2 i j &&
            (match Field.get ctx.color_reps i, Field.get ctx.color_reps j with
            | Color.SUN (n1), Color.SUN (n2) ->
                n1 = - n2 && n2 > 0
            | _, _ -> false)
        | E (i, j, k) ->
            distinct3 i j k &&
            (match Field.get ctx.color_reps i,
              Field.get ctx.color_reps j, Field.get ctx.color_reps k with
            | Color.SUN (n1), Color.SUN (n2), Color.SUN (n3) ->
                n1 = 3 && n2 = 3 && n3 = 3 ||
                n1 = -3 && n2 = -3 && n3 = -3
              | _, _, _ -> false)
        | T (a, i, j) ->
            distinct3 a i j &&
            (match Field.get ctx.color_reps a,
              Field.get ctx.color_reps i, Field.get ctx.color_reps j with
            | Color.AdjSUN(n1), Color.SUN (n2), Color.SUN (n3) ->
                n1 = n3 && n2 = - n3 && n3 > 0
            | _, _, _ -> false)
        | F (a, b, c) ->
            distinct3 a b c &&
            (match Field.get ctx.color_reps a,
              Field.get ctx.color_reps b, Field.get ctx.color_reps c with
            | Color.AdjSUN(n1), Color.AdjSUN (n2), Color.AdjSUN (n3) ->
                n1 = n2 && n2 = n3 && n1 > 0
            | _, _, _ -> false)

    let primitive_indices = function
      | D (_, _) -> []
      | E (_, _, _) -> []
      | T (a, _, _) -> [a]
      | F (a, b, c) -> [a; b; c]

    let indices p =
      ThoList.flatmap primitive_indices p

    let contraction_ok p =
      List.for_all
        (fun (n, _) -> n = 2)
        (ThoList.classify (indices p))

    type tensor = int * primitive list

    let map_tensor f (factor, primitives) =
      (factor, List.map (map_primitive f) primitives)

    let tensor_ok context (_, primitives) =
      List.for_all (primitive_ok context) primitives

  end

type t =
    { fields : string array;
      lorentz : Lorentz.tensor list;
      color : Color.tensor list }

module Test (M : Model.T) : Test =
  struct

    module Permutation = Permutation.Default

    let context_of_flavors flavors =
      { arity = Array.length flavors;
        lorentz_reps = Array.map M.lorentz flavors;
        color_reps = Array.map M.color flavors }

    let context_of_flavor_names names =
      context_of_flavors (Array.map M.flavor_of_string names)

    let context_of_vertex v =
      context_of_flavor_names v.fields

    let ok v =
      let context = context_of_vertex v in
      List.for_all (Lorentz.tensor_ok context) v.lorentz &&
        List.for_all (Color.tensor_ok context) v.color

    module PM =
      Partial.Make (struct type t = field let compare = compare end)

    let id x = x

    let permute v p =
      let context = context_of_vertex v in
      let sorted =
        List.map
          (Field.of_int context)
          (ThoList.range 0 (Array.length v.fields - 1)) in
      let permute =
        PM.apply (PM.of_lists sorted (List.map (Field.of_int context) p)) in
      { fields = Permutation.array (Permutation.of_list p) v.fields;
        lorentz = List.map
          (Lorentz.map_tensor id permute id permute id permute) v.lorentz;
        color = List.map (Color.map_tensor permute) v.color }

    let permutations v =
      List.map (permute v)
        (Combinatorics.permute (ThoList.range 0 (Array.length v.fields - 1)))

    let wf_declaration flavor =
      match M.lorentz (M.flavor_of_string flavor) with
      | Coupling.Vector -> "vector"
      | Coupling.Spinor -> "spinor"
      | Coupling.ConjSpinor -> "conjspinor"
      | _ -> failwith "wf_declaration: incomplete"

    module Chiral = Lorentz.Dirac(Lorentz.Chiral)

    let write_fusion v =
      match Array.to_list v.fields with
      | lhs :: rhs ->
          let name = lhs ^ "_of_" ^ String.concat "_" rhs in
          let momenta = List.map (fun n -> "k_" ^ n) rhs in
          Printf.printf "pure function %s (%s) result (%s)\n"
            name (String.concat ", "
                    (List.flatten
                       (List.map2 (fun wf p -> [wf; p]) rhs momenta)))
            lhs;
          Printf.printf "  type(%s) :: %s\n" (wf_declaration lhs) lhs;
          List.iter
            (fun wf ->
              Printf.printf "  type(%s), intent(in) :: %s\n"
                (wf_declaration wf) wf)
            rhs;
          List.iter
            (Printf.printf "  type(momentum), intent(in) :: %s\n")
            momenta;
          let rhs1 = List.hd rhs
          and rhs2 = List.hd (List.tl rhs) in
          begin match M.lorentz (M.flavor_of_string lhs) with
          | Coupling.Vector ->
              begin
                for mu = 0 to 3 do
                  Printf.printf "  %s(%d) =" lhs mu;
                  for i = 1 to 4 do
                    for j = 1 to 4 do
                      match Chiral.vector mu i j with
                      | Lorentz.Complex.Z -> ()
                      | c ->
                          Printf.printf " + %s*%s(%d)*%s(%d)"
                            (Lorentz.Complex.to_fortran c) rhs1 i rhs2 j
                    done
                  done;
                  Printf.printf "\n"
                done
              end;
          | Coupling.Spinor | Coupling.ConjSpinor ->
              begin
                for i = 1 to 4 do
                  Printf.printf "  %s(%d) =" lhs i;
                  for mu = 0 to 3 do
                    for j = 1 to 4 do
                      match Chiral.vector mu i j with
                      | Lorentz.Complex.Z -> ()
                      | c ->
                          Printf.printf " + %s*%s(%d)*%s(%d)"
                            (Lorentz.Complex.to_fortran c) rhs1 mu rhs2 j
                    done
                  done;
                  Printf.printf "\n"
                done
              end;
          | _ -> failwith "write_fusion: incomplete"
          end;
          Printf.printf "end function %s\n" name;
          ()
      | [] -> ()

    let write_fusions v =
      List.iter write_fusion (permutations v)

(* Testing: *)

    let vector_field context i =
      Lorentz.Vector (Lorentz.F (Field.of_int context i))

    let spinor_field context i =
      Lorentz.Spinor (Lorentz.F (Field.of_int context i))

    let conjspinor_field context i =
      Lorentz.ConjSpinor (Lorentz.F (Field.of_int context i))

    let mu = Lorentz.Vector (Lorentz.I 0)
    and nu = Lorentz.Vector (Lorentz.I 1)

    let tbar_gl_t = [| "tbar"; "gl"; "t" |]
    let context = context_of_flavor_names tbar_gl_t
     
    let vector_current_ok =
      { fields = tbar_gl_t;
        lorentz = [ (1, [Lorentz.V (vector_field context 1,
                                    conjspinor_field context 0,
                                    spinor_field context 2)]) ];
        color = [ (1, [Color.T (Field.of_int context 1,
                                Field.of_int context 0,
                                Field.of_int context 2)])] }

    let vector_current_vector_misplaced =
      { fields = tbar_gl_t;
        lorentz = [ (1, [Lorentz.V (vector_field context 2,
                                    conjspinor_field context 0,
                                    spinor_field context 2)]) ];
        color = [ (1, [Color.T (Field.of_int context 1,
                                Field.of_int context 0,
                                Field.of_int context 2)])] }

    let vector_current_spinor_misplaced =
      { fields = tbar_gl_t;
        lorentz = [ (1, [Lorentz.V (vector_field context 1,
                                    conjspinor_field context 0,
                                    spinor_field context 1)]) ];
        color = [ (1, [Color.T (Field.of_int context 1,
                                Field.of_int context 0,
                                Field.of_int context 2)])] }

    let vector_current_conjspinor_misplaced =
      { fields = tbar_gl_t;
        lorentz = [ (1, [Lorentz.V (vector_field context 1,
                                    conjspinor_field context 1,
                                    spinor_field context 2)]) ];
        color = [ (1, [Color.T (Field.of_int context 1,
                                Field.of_int context 0,
                                Field.of_int context 2)])] }

    let vector_current_out_of_bounds () =
      { fields = tbar_gl_t;
        lorentz = [ (1, [Lorentz.V (mu,
                                    conjspinor_field context 3,
                                    spinor_field context 2)]) ];
        color = [ (1, [Color.T (Field.of_int context 1,
                                Field.of_int context 0,
                                Field.of_int context 2)])] }

    let vector_current_color_mismatch =
      let names = [| "t"; "gl"; "t" |] in
      let context = context_of_flavor_names names in
      { fields = names;
        lorentz = [ (1, [Lorentz.V (mu,
                                    conjspinor_field context 0,
                                    spinor_field context 2)]) ];
        color = [ (1, [Color.T (Field.of_int context 1,
                                Field.of_int context 0,
                                Field.of_int context 2)])] }

    let wwzz = [| "W+"; "W-"; "Z"; "Z" |]
    let context = context_of_flavor_names wwzz

    let anomalous_couplings =
      { fields = wwzz;
        lorentz = [ (1, [ Lorentz.K (mu, Field.of_int context 0);
                          Lorentz.K (mu, Field.of_int context 1) ]) ];
        color = [ ] }
      
    let anomalous_couplings_index_mismatch =
      { fields = wwzz;
        lorentz = [ (1, [ Lorentz.K (mu, Field.of_int context 0);
                          Lorentz.K (nu, Field.of_int context 1) ]) ];
        color = [ ] }
      
    exception Inconsistent_vertex

    let example () =
      if not (ok vector_current_ok) then begin
        raise Inconsistent_vertex
      end;
      write_fusions vector_current_ok

    open OUnit

    let vertex_indices_ok =
      "indices/ok" >::
        (fun () ->
          List.iter
            (fun v ->
              assert_bool "vector_current" (ok v))
            (permutations vector_current_ok))
                
    let vertex_indices_broken =
      "indices/broken" >::
        (fun () ->
          assert_bool "vector misplaced"
            (not (ok vector_current_vector_misplaced));
          assert_bool "conjugate spinor misplaced"
            (not (ok vector_current_spinor_misplaced));
          assert_bool "conjugate spinor misplaced"
            (not (ok vector_current_conjspinor_misplaced));
          assert_raises (Field.Out_of_range 3)
            vector_current_out_of_bounds;
          assert_bool "color mismatch"
            (not (ok vector_current_color_mismatch)))
                
    let anomalous_couplings_ok =
      "anomalous_couplings/ok" >::
        (fun () ->
          assert_bool "anomalous couplings"
            (ok anomalous_couplings))
                
    let anomalous_couplings_broken =
      "anomalous_couplings/broken" >::
        (fun () ->
          assert_bool "anomalous couplings"
            (not (ok anomalous_couplings_index_mismatch)))
                
    let suite =
      "Vertex" >:::
        [vertex_indices_ok;
         vertex_indices_broken;
         anomalous_couplings_ok;
         anomalous_couplings_broken]
      
  end

