(* vertex_syntax.ml --

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

(* \thocwmodulesection{Abstract Syntax} *)

exception Syntax_Error of string * Lexing.position * Lexing.position

module Token =
  struct

    type t =
    | Digit of int
    | Token of string
    | Scripted of scripted
    | List of t list

    and scripted = 
      { stem : t;
	prefix : prefix list;
	super : t list;
	sub : t list }

    and prefix =
    | Bar | Hat | Tilde
    | Dagger | Star
    | Prime

    let prefix_of_string = function
      | "\\bar" | "\\overline" -> Bar
      | "\\hat" | "\\widehat" -> Hat
      | "\\tilde" | "\\widetilde" -> Tilde
      | "\\dagger" -> Dagger
      | "*" | "\\ast" -> Star
      | "\\prime" -> Prime
      | _ -> invalid_arg "Vertex_Syntax.Token.string_to_prefix"

    let prefix_to_string = function
      | Bar -> "\\bar"
      | Hat -> "\\hat"
      | Tilde -> "\\tilde"
      | Dagger -> "\\dagger"
      | Star -> "*"
      | Prime -> "\\prime"

    let wrap_scripted = function
      | Scripted st -> st
      | t ->  { stem = t; prefix = []; super = []; sub = [] }

    let wrap_list = function
      | List tl -> tl
      | _ as t -> [t]

    let digit i = 
      if i >= 0 && i <= 9 then
	Digit i
      else
	invalid_arg ("Vertex_Syntax.Token.digit: " ^ string_of_int i)

    let token s =
      Token s

    let list = function
      | [] -> List []
      | [Scripted {stem = t; prefix = []; super = []; sub = []}] -> t
      | [t] -> t
      | tl ->  List tl

    let optional = function
      | None -> []
      | Some t -> wrap_list t

    let scripted prefix token (super, sub) =
      match token, prefix, super, sub with
      | _, [], None, None -> token
      | (Digit _ | Token _ | List _) as t, _, _, _ ->
	Scripted { stem = t;
		   prefix =  List.map prefix_of_string prefix;
		   super = optional super;
		   sub = optional sub }
      | Scripted st, _, _, _ ->
	Scripted { stem = st.stem;
		   prefix =  List.map prefix_of_string prefix @ st.prefix;
		   super = st.super @ optional super;
		   sub = st.sub @ optional sub }

    let rec stem = function
      | Digit _ | Token _ as t -> t
      | Scripted { stem = t } -> stem t
      | List tl ->
	begin match List.rev tl with
	| [] -> List []
	| t :: _ -> stem t
	end

    (* Strip superfluous [List] and [Scripted] constructors. *)
    (* NB: This might be unnecessary, if we used smart constructors. *)

    let rec strip = function
      | Digit _ | Token _ as t -> t
      | Scripted { stem = t; prefix = []; super = []; sub = [] } -> strip t
      | Scripted { stem = t; prefix = prefix; super = super; sub = sub } ->
	Scripted { stem = strip t;
		   prefix = prefix;
		   super = List.map strip super;
		   sub = List.map strip sub }
      | List tl ->
	begin match List.map strip tl with
	| [] -> List []
	| [t] -> t
	| tl ->  List tl
	end

    (* Recursively merge nested [List] and [Scripted] constructors. *)
    (* NB: This might be unnecessary, if we used smart constructors. *)

    let rec flatten = function
      | Digit _ | Token _ as t -> t
      | List tl -> flatten_list tl
      | Scripted st -> flatten_scripted st

    and flatten_list tl =
      match List.map flatten tl with
      | [] -> List []
      | [t] -> t
      | tl ->  List tl

    and flatten_scripted = function
      | { stem = t; prefix = []; super = []; sub = [] } -> t
      | { stem = t; prefix = prefix; super = super; sub = sub } ->
	let super = List.map flatten super
	and sub = List.map flatten sub in
	begin match flatten t with
	| Digit _ | Token _ | List _ as t ->
	  Scripted { stem = t;
		     prefix = prefix;
		     super = super;
		     sub = sub }
	| Scripted st ->
	  Scripted { stem = st.stem;
		     prefix = prefix @ st.prefix;
		     super = st.super @ super;
		     sub = st.sub @ sub }
	end

    let ascii_A = Char.code 'A'
    let ascii_Z = Char.code 'Z'
    let ascii_a = Char.code 'a'
    let ascii_z = Char.code 'z'

    let is_char c =
      let a = Char.code c in
      (ascii_A <= a && a <= ascii_Z) || (ascii_a <= a && a <= ascii_z)

    let is_backslash c =
      c = '\\'

    let first_char s =
      s.[0]

    let last_char s =
      s.[String.length s - 1]

    let rec to_string = function
      | Digit i -> string_of_int i
      | Token s -> s
      | Scripted t -> scripted_to_string t
      | List tl -> "{" ^ list_to_string tl ^ "}"

    and list_to_string = function
      | [] -> ""
      | [Scripted { stem = t; super = []; sub = [] }] -> to_string t
      | [Scripted _ as t] -> "{" ^ to_string t ^ "}"
      | [t] -> to_string t
      | tl -> "{" ^ concat_tokens tl ^ "}"

    and scripted_to_string t =
      let super =
	match t.super with
	| [] -> ""
	| tl -> "^" ^ list_to_string tl
      and sub =
	match t.sub with
	| [] -> ""
	| tl -> "_" ^ list_to_string tl in
      String.concat "" (List.map prefix_to_string t.prefix) ^
	to_string t.stem ^ super ^ sub

    and required_space t1 t2 =
      let required_space' s1 s2 =
	if is_backslash (first_char s2) then
	  []
	else if is_backslash (first_char s1) && is_char (last_char s1) then
	  [Token " "]
	else
	  [] in
      match t1, t2 with
      | Token s1, Token s2 -> required_space' s1 s2
      | Scripted s1, Token s2 -> required_space' (scripted_to_string s1) s2
      | Token s1, Scripted s2 -> required_space' s1 (scripted_to_string s2)
      | Scripted s1, Scripted s2 ->
	required_space' (scripted_to_string s1) (scripted_to_string s2)
      | List _, _ | _, List _ | _, Digit _ | Digit _, _ -> []

    and interleave_spaces tl =
      ThoList.interleave_nearest required_space tl

    and concat_tokens tl =
      String.concat "" (List.map to_string (interleave_spaces tl)) 

    let	compare t1 t2 =
      pcompare t1 t2

  end

module Expr =
  struct

    type t =
    | Integer of int
    | Sum of t list | Diff of t * t
    | Product of t list | Ratio of t * t
    | Function of Token.t * t list

    let integer i = Integer i

    let rec add a b =
      match a, b with
      | Integer a, Integer b -> Integer (a + b)
      | Sum a, Sum b -> Sum (a @ b)
      | Sum a, b -> Sum (a @ [b])
      | a, Sum b -> Sum (a :: b)
      | a, b -> Sum ([a; b])

    (* (a1 - a2) - (b1 - b2) = (a1 + b2) - (a2 + b1) *)
    (* (a1 - a2) - b = a1 - (a2 + b) *)
    (* a - (b1 - b2) = (a + b2) - b1 *)

    and sub a b =
      match a, b with
      | Integer a, Integer b -> Integer (a - b)
      | Diff (a1, a2), Diff (b1, b2) -> Diff (add a1 b2, add a2 b1)
      | Diff (a1, a2), b -> Diff (a1, add a2 b)
      | a, Diff (b1, b2) -> Diff (add a b2, b1)	
      | a, b -> Diff (a, b)	

    and mult a b =
      match a, b with
      | Integer a, Integer b -> Integer (a * b)
      | Product a, Product b -> Product (a @ b)
      | Product a, b -> Product (a @ [b])
      | a, Product b -> Product (a :: b)
      | a, b -> Product ([a; b])

    and div a b =
      match a, b with
      | Ratio (a1, a2), Ratio (b1, b2) -> Ratio (mult a1 b2, mult a2 b1)
      | Ratio (a1, a2), b -> Ratio (a1, mult a2 b)
      | a, Ratio (b1, b2) -> Ratio (mult a b2, b1)	
      | a, b -> Ratio (a, b)	

    let apply f args =
      Function (f, args)

    let rec to_string = function
      | Integer i -> string_of_int i
      | Sum ts -> String.concat "+" (List.map to_string ts)
      | Diff (t1, t2) -> to_string t1 ^ "-" ^ to_string t2
      | Product ts -> String.concat "*" (List.map to_string ts)
      | Ratio (t1, t2) -> to_string t1 ^ "/" ^ to_string t2
      | Function (f, args) ->
	Token.to_string f ^
	  String.concat ""
	  (List.map (fun arg -> "{" ^ to_string arg ^ "}") args)

  end

(*i module TLSet = Set.Make (struct type t = Token.t list let compare = compare end) i*)

module Particle =
  struct

    type name =
    | Neutral of Token.t
    | Charged of Token.t * Token.t

    type attr =
    | TeX of Token.t list | TeX_Anti of Token.t list
    | Alias of Token.t list | Alias_Anti of Token.t list
    | Fortran of Token.t list | Fortran_Anti of Token.t list
    | Spin of Expr.t | Charge of Expr.t
    | Color of Token.t list * Token.t list
    | Mass of Token.t list | Width of Token.t list

(*i
    (* Combine the sets of aliases and use the
       rightmost version of the other attributes.  *)
    let rec cons_attr a = function
      | [] -> [a]
      | a' :: alist ->
	match a, a' with
	| TeX tl, TeX tl' -> a' :: alist
	| TeX_Anti tl, TeX_Anti tl' -> a' :: alist
	| Aliases tl, Aliases tl' ->
	  Aliases (TLSet.union tl tl') :: alist
	| Aliases_Anti tl, Aliases_Anti tl' ->
	  Aliases_Anti (TLSet.union tl tl') :: alist
	| Fortran tl, Fortran tl' -> a' :: alist
	| Fortran_Anti tl, Fortran_Anti tl' -> a' :: alist
	| Spin tl, Spin tl' -> a' :: alist
	| Color tl, Color tl' -> a' :: alist
	| Charge tl, Charge tl' -> a' :: alist
	| Mass tl, Mass tl' -> a' :: alist
	| Width tl, Width tl' -> a' :: alist
	| _, _ -> a' :: cons_attr a alist
i*)

    type t =
      { name : name;
	attr : attr list }

    let name_to_string = function
      | Neutral p ->
	 "\\neutral{" ^ Token.to_string p ^ "}"
      | Charged (p, ap) ->
	"\\charged{" ^ Token.to_string p ^ "}{" ^ Token.to_string ap ^ "}"

    let attr_to_string = function
      | TeX tl -> "\\tex{" ^ Token.list_to_string tl ^ "}"
      | TeX_Anti tl -> "\\anti\\tex{" ^ Token.list_to_string tl ^ "}"
      | Alias tl -> "\\alias{" ^ Token.list_to_string tl ^ "}"
      | Alias_Anti tl -> "\\anti\\alias{" ^ Token.list_to_string tl ^ "}"
      | Fortran tl -> "\\fortran{" ^ Token.list_to_string tl ^ "}"
      | Fortran_Anti tl -> "\\anti\\fortran{" ^ Token.list_to_string tl ^ "}"
      | Spin e -> "\\spin{" ^ Expr.to_string e ^ "}"
      | Color ([], rep) -> "\\color{" ^ Token.list_to_string rep ^ "}"
      | Color (group, rep) ->
	 "\\color[" ^ Token.list_to_string group ^ "]{"	 ^
	   Token.list_to_string rep ^ "}"
      | Charge e -> "\\charge{" ^ Expr.to_string e ^ "}"
      | Mass tl -> "\\mass{" ^ Token.list_to_string tl ^ "}"
      | Width tl -> "\\width{" ^ Token.list_to_string tl ^ "}"

    let to_string p =
      name_to_string p.name ^
	String.concat "" (List.map attr_to_string (List.sort compare p.attr))
	
  end

module Parameter =
  struct

    type attr =
    | TeX of Token.t list
    | Alias of Token.t list
    | Fortran of Token.t list

    type t' =
      { name : Token.t;
	value : Expr.t;
	attr : attr list}

(*i
    let rec cons_attr a = function
      | [] -> [a]
      | a' :: alist ->
	match a, a' with
	| TeX tl, TeX tl' -> a' :: alist
	| Aliases tl, Aliases tl' ->
	  Aliases (TLSet.union tl tl') :: alist
	| Fortran tl, Fortran tl' -> a' :: alist
	| _, _ -> a' :: cons_attr a alist
i*)

    type t =
    | Parameter of t'
    | Derived of t'

    let attr_to_string = function
      | TeX tl -> "\\tex{" ^ Token.list_to_string tl ^ "}"
      | Alias tl -> "\\alias{" ^ Token.list_to_string tl ^ "}"
      | Fortran tl -> "\\fortran{" ^ Token.list_to_string tl ^ "}"

    let to_string' p =
      "{" ^ Token.to_string p.name ^ "}{" ^ Expr.to_string p.value ^ "}" ^
	String.concat "" (List.map attr_to_string p.attr)

    let to_string = function
      | Parameter p -> "\\parameter" ^ to_string' p
      | Derived p -> "\\derived" ^ to_string' p

  end

module Lie =
  struct

    type group =
    | SU of int | U of int
    | SO of int | O of int
    | Sp of int
    | E6 | E7 | E8 | F4 | G2

    module T = Token

    let default_group = SU 3

    let invalid_group s =
      invalid_arg ("Vertex.Lie.group_of_string: " ^ s)

    let series s name n =
      match name, n with
      | "SU", n when n > 1 -> SU n
      | "U", n when n >= 1  -> U n
      | "SO", n when n > 1  -> SO n
      | "O", n when n >= 1  -> O n
      | "Sp", n when n >= 2  -> Sp n
      | _ -> invalid_group s

    let exceptional s name n =
      match name, n with
      | "E", 6 -> E6
      | "E", 7 -> E7
      | "E", 8 -> E8
      | "F", 4 -> F4
      | "G", 2 -> G2
      | _ -> invalid_group s

    let group_of_string s =
      try
	Scanf.sscanf s "%_[{]%[SUOp](%d)%_[}]%!" (series s)
      with
      | _ ->
	 try
	   Scanf.sscanf s "%_[{]%[EFG]_%d%_[}]%!" (exceptional s)
	 with
	 | _ -> invalid_group s

    let group_to_string = function
      | SU n -> "SU(" ^ string_of_int n ^ ")"
      | U n -> "U(" ^ string_of_int n ^ ")"
      | SO n -> "SO(" ^ string_of_int n ^ ")"
      | O n -> "O(" ^ string_of_int n ^ ")"
      | Sp n -> "Sp(" ^ string_of_int n ^ ")"
      | E6 -> "E6"
      | E7 -> "E7"
      | E8 -> "E8"
      | F4 -> "F4"
      | G2 -> "G2"

    type rep = int

    let rep_of_string group rep =
      match group with
      | SU 3 ->
	 begin
	   match rep with
	   | "3" -> 3
	   | "\\bar 3" -> -3
	   | "8" -> 8
	   | _ ->
	      invalid_arg ("Vertex.Lie.rep_of_string:" ^
			     " unsupported representation " ^ rep ^
			     " of " ^ group_to_string group)
	 end
      | _ -> invalid_arg ("Vertex.Lie.rep_of_string:" ^
			    " unsupported group " ^ group_to_string group)

    let rep_to_string r =
      string_of_int r

    type t = group * rep

  end

module Lorentz =
  struct

    type rep =
    | Scalar | Vector
    | Dirac | ConjDirac | Majorana
    | Weyl | ConjWeyl

  end

module Index =
  struct

    type attr =
    | Color of Token.t list * Token.t list
    | Flavor of Token.t list * Token.t list
    | Lorentz of Token.t list

    type t =
      { name : Token.t;
	attr : attr list }

    let attr_to_string = function
      | Color ([], rep) -> "\\color{" ^ Token.list_to_string rep ^ "}"
      | Color (group, rep) ->
	 "\\color[" ^ Token.list_to_string group ^ "]{"	 ^
	   Token.list_to_string rep ^ "}"
      | Flavor ([], rep) -> "\\flavor{" ^ Token.list_to_string rep ^ "}"
      | Flavor (group, rep) ->
	 "\\flavor[" ^ Token.list_to_string group ^ "]{"	 ^
	   Token.list_to_string rep ^ "}"
      | Lorentz tl -> "\\lorentz{" ^ Token.list_to_string tl ^ "}"

    let to_string i =
      "\\index{" ^ Token.to_string i.name ^ "}" ^
	String.concat "" (List.map attr_to_string i.attr)
  end

module Tensor =
  struct

    type attr =
    | Color of Token.t list * Token.t list
    | Flavor of Token.t list * Token.t list
    | Lorentz of Token.t list

    type t =
      { name : Token.t;
	attr : attr list }

    let attr_to_string = function
      | Color ([], rep) -> "\\color{" ^ Token.list_to_string rep ^ "}"
      | Color (group, rep) ->
	 "\\color[" ^ Token.list_to_string group ^ "]{"	 ^
	   Token.list_to_string rep ^ "}"
      | Flavor ([], rep) -> "\\flavor{" ^ Token.list_to_string rep ^ "}"
      | Flavor (group, rep) ->
	 "\\flavor[" ^ Token.list_to_string group ^ "]{"	 ^
	   Token.list_to_string rep ^ "}"
      | Lorentz tl -> "\\lorentz{" ^ Token.list_to_string tl ^ "}"

    let to_string t =
      "\\tensor{" ^ Token.to_string t.name ^ "}" ^
	String.concat "" (List.map attr_to_string t.attr)
  end

module File_Tree =
  struct

    type declaration =
    | Particle of Particle.t
    | Parameter of Parameter.t
    | Index of Index.t
    | Tensor of Tensor.t
    | Vertex of Expr.t * Token.t
    | Include of string

    type t = declaration list

    let empty = []

  end

module File =
  struct

    type declaration =
    | Particle of Particle.t
    | Parameter of Parameter.t
    | Index of Index.t
    | Tensor of Tensor.t
    | Vertex of Expr.t * Token.t

    type t = declaration list

    let empty = []

    (* We allow to include a file more than once, but we don't
       optimize by memoization, because we assume that this will
       be rare.  However to avoid infinite loops when including
       a child, we make sure that it has not yet been included as
       a parent.  *)

    let expand_includes parser unexpanded =
      let rec expand_includes' parents unexpanded expanded =
	List.fold_right (fun decl decls ->
	  match decl with
	  | File_Tree.Particle p -> Particle p :: decls
	  | File_Tree.Parameter p -> Parameter p :: decls
	  | File_Tree.Index i -> Index i :: decls
	  | File_Tree.Tensor t -> Tensor t :: decls
	  | File_Tree.Vertex (e, v) -> Vertex (e, v) :: decls
	  | File_Tree.Include f ->
	     if List.mem f parents then
	       invalid_arg ("cyclic \\include{" ^ f ^ "}")
	     else
	       expand_includes' (f:: parents) (parser f) decls)
	  unexpanded expanded in
      expand_includes' [] unexpanded []

    let to_strings decls =
      List.map
	(function
	| Particle p -> Particle.to_string p
	| Parameter p -> Parameter.to_string p
	| Index i -> Index.to_string i
	| Tensor t -> Tensor.to_string t
	| Vertex (Expr.Integer 1, t) -> 
	  "\\vertex{" ^ Token.to_string t ^ "}"
	| Vertex (e, t) ->
	  "\\vertex[" ^ Expr.to_string e ^ "]{" ^
	    Token.to_string t ^ "}")
	decls

  end
