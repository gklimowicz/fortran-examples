(* vertex_syntax.mli --

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

(* The concrete syntax described below is modelled on \LaTeX{}
   and correct model descriptions should be correct \LaTeX-input
   (provided a few simple macros have been loaded. *)


(* \thocwmodulesection{Abstract Syntax} *)

exception Syntax_Error of string * Lexing.position * Lexing.position

(* \thocwmodulesubsection{Tokens} *)

(* Tokenization follows \TeX's rules. *)
       
module Token :
  sig

    (* Single-character tokens other than digits
       are stored as one character strings.
       Multi-character tokens like \verb+\psi+ are stored
       as a string \emph{including} the leading \verb+\+.
       Since \verb+a_12+
       is interpreted by \TeX{} as \verb+{a_1}2+, we can not
       use the lexer to construct integers, but interpret them
       as lists of digits.  Below, in [Expr], the parser can
       interpret then as integers.  *)
    type t = private
    | Digit of int
    | Token of string
    | Scripted of scripted
    | List of t list

    (* TODO: investigate if it is possible to introduce [stem]
       as a separate type to allow more fine-grained compile-time
       checks. *)

    (* In addition to super- and subscripts, there are prefixes
       such as \verb+\bar+, \verb+\hat+, etc.  *)
    and scripted = private
      { stem : t;
	prefix : prefix list;
	super : t list;
	sub : t list }

    and prefix =
    | Bar | Hat | Tilde
    | Dagger | Star
    | Prime

    val prefix_of_string : string -> prefix
    val prefix_to_string : prefix -> string

    (* Smart constructors that avoid redundant nestings of lists
       and scripted tokens with empty scripts. *)
    val digit : int -> t
    val token : string -> t
    val scripted : string list -> t -> t option * t option -> t
    val list : t list -> t

    (* If it's [Scripted], return unchanged, else as a scripted
       token with empty prefix, super- and subscripts. *)
    val wrap_scripted : t -> scripted

    (* If it's a [List], return the list itself,
       otherwise a singleton list. *)
    val wrap_list : t -> t list

    (* Recursively strip all prefixes, super- and subscripts and
       return only the LAST token in a list.
       I.e. [stem "\\bar\\psi_i"] and [stem "\\bar{\\phi\\psi}'"]
       both yield ["\\psi"]. *)
    val stem : t -> t

    (* Unparse the abstract syntax.   Since the smart constructors
       perform some normalization and minimize nested braces, the
       result is not guaranteed to be identical to the string that
       has been parsed, just equivalent. *)
    val to_string : t -> string
    val scripted_to_string : scripted -> string
    val list_to_string : t list -> string

    val compare : t -> t -> int

  end

(* \thocwmodulesubsection{Expressions} *)

(* A straightforward type for recursive expressions.  Note
   that values (a.\,k.\,a.~variables) are represented as
   functions with an empty argument list. *)

module Expr :
  sig

    type t =
    | Integer of int
    | Sum of t list | Diff of t * t
    | Product of t list | Ratio of t * t
    | Function of Token.t * t list

    val integer : int -> t
    val add : t -> t -> t
    val sub : t -> t -> t
    val mult : t -> t -> t
    val div : t -> t -> t
    val apply : Token.t -> t list -> t

    val to_string : t -> string

  end

(* \thocwmodulesubsection{Particle Declarations} *)

module Particle :
  sig

    (* Neutral particles are known by a single name,
       charged particles also by the name of the
       anti-particle,  \ldots *)
    type name =
    | Neutral of Token.t
    | Charged of Token.t * Token.t

    (* \ldots{} and a list of attributes: aliases, external
       representations for \LaTeX{} and Fortran, quantum
       numbers and symbols for mass and width. *)
    type attr =
    | TeX of Token.t list | TeX_Anti of Token.t list
    | Alias of Token.t list | Alias_Anti of Token.t list
    | Fortran of Token.t list | Fortran_Anti of Token.t list
    | Spin of Expr.t | Charge of Expr.t
    | Color of Token.t list * Token.t list
    | Mass of Token.t list | Width of Token.t list

    type t =
      { name : name;
	attr : attr list }

    (* Unparsing: *)
    val to_string : t -> string

  end

(* \thocwmodulesubsection{Parameter Declarations} *)

module Parameter :
  sig

    type attr =
    | TeX of Token.t list
    | Alias of Token.t list
    | Fortran of Token.t list

    type t' =
      { name : Token.t;
	value : Expr.t;
	attr : attr list}

    type t =
    | Parameter of t'
    | Derived of t'

    (*i val cons_attr : attr -> attr list -> attr list i*)
    val to_string : t -> string

  end

(* \thocwmodulesubsection{Lie Groups and Algebras} *)

module Lie :
  sig

    (* The full list [SU of int | U of int | SO of int | O of int
       | Sp of int | E6 | E7 | E8 | F4 | G2] is not realistic.
       In practice, we will concentrate on SU(3) for now. *)

    type group

    val default_group : group (* SU(3), of course *)
    val group_of_string : string -> group
    val group_to_string : group -> string

    (* For now, we only support the~$\mathbf{3}$, $\bar{\mathbf{3}}$
       and~$\mathbf{8}$ of SU(3). *)

    type rep

    val rep_of_string : group -> string -> rep
    val rep_to_string : rep -> string

    type t = group * rep

  end

(* \thocwmodulesubsection{Lorentz Representations} *)

module Lorentz :
  sig

    type rep =
    | Scalar | Vector
    | Dirac | ConjDirac | Majorana
    | Weyl | ConjWeyl

  end

(* \thocwmodulesubsection{Indices} *)

module Index :
  sig

    type attr =
    | Color of Token.t list * Token.t list
    | Flavor of Token.t list * Token.t list
    | Lorentz of Token.t list

    type t =
      { name : Token.t;
	attr : attr list }

    val to_string : t -> string

  end

(* \thocwmodulesubsection{Tensors} *)

module Tensor :
  sig

    type attr =
    | Color of Token.t list * Token.t list
    | Flavor of Token.t list * Token.t list
    | Lorentz of Token.t list

    type t =
      { name : Token.t;
	attr : attr list }

    val to_string : t -> string

  end

(* \thocwmodulesubsection{Files} *)

(* The abstract representation of a file, immediately after lexical
   and syntactical analysis and before any type checking
   or semantic analysis, is a list of declarations.  *)

(* There is one version with unexpanded \verb+\include+
   statements. *)

module File_Tree :
  sig

    type declaration =
    | Particle of Particle.t
    | Parameter of Parameter.t
    | Index of Index.t
    | Tensor of Tensor.t
    | Vertex of Expr.t * Token.t
    | Include of string

    type t = declaration list

    val empty : t

  end

(* A linear file, just like [File_Tree], but with all
   the \verb+\include+ statements expanded. *)

module File :
  sig

    type declaration =
    | Particle of Particle.t
    | Parameter of Parameter.t
    | Index of Index.t
    | Tensor of Tensor.t
    | Vertex of Expr.t * Token.t

    type t = declaration list

    val empty : t

    (* [expand_includes parser file_tree] recursively
       expands all include statemens in [file_tree], using
       [parser] to map a filename to a [File_Tree.t]. *)
    val expand_includes : (string -> File_Tree.t) -> File_Tree.t -> t

    val to_strings : t -> string list

  end

