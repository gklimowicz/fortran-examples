(* vertex.mli --

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

val parse_string : string -> Vertex_syntax.File.t
val parse_file : string -> Vertex_syntax.File.t

module type Test =
  sig
    val example : unit -> unit
    val suite : OUnit.test
  end

module Test (M : Model.T) : Test

module Parser_Test : Test
module Modelfile_Test : Test

(*i
module Symbol :
  sig
    type table
    val load : Vertex_syntax.File.t -> table
  end

module Vertex :
  sig

    type factor =
      { stem : Vertex_syntax.Token.t;
	prefix : string list;
	particle : Vertex_syntax.Token.t list;
	color : Vertex_syntax.Token.t list;
	lorentz : Vertex_syntax.Token.t list;
	flavor : Vertex_syntax.Token.t list;
	other : Vertex_syntax.Token.t list }

    val factor_of_token : Symbol.table -> Vertex_syntax.Token.scripted -> factor

  end
i*)

