(* omega_unit.ml --

   Copyright (C) 1999-2014 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>

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

module M = Modellib_SM.SM(Modellib_SM.SM_no_anomalous)
module List_Test = Permutation.Test(Permutation.Using_Lists)
module Array_Test = Permutation.Test(Permutation.Using_Arrays)
module Vertex_Test = Vertex.Test(M)
module Parser_Test = Vertex.Parser_Test
module Model_Test = Vertex.Modelfile_Test

let _ =
  let my_name = Sys.argv.(0) in
  let skip_tests = ref false
  and skip_example = ref false
  and timing = ref false
  and verbose = ref false
  and usage = "usage: " ^ my_name ^ " ..." in
  Arg.parse
    (Arg.align 
       [ ("-skip-tests", Arg.Set skip_tests, " skip the tests");
	 ("-skip-example", Arg.Set skip_example, " skip the example");
	 ("-timing", Arg.Set timing, " provide timing information");
	 ("-v", Arg.Set verbose, " be more verbose");
	 ("-verbose", Arg.Set verbose, " be more verbose") ])
    (fun s -> raise (Arg.Bad s))
    usage;
  if not !skip_tests then begin
    let suite =
      OUnit.(>:::) "All" 
	[Vertex_Test.suite;
	 Parser_Test.suite;
	 Model_Test.suite] in
    ignore (OUnit.run_test_tt ~verbose:!verbose suite)
  end;
  if !timing then begin
    print_endline "List based:";
    List_Test.time ();
    print_endline "Array based:";
    Array_Test.time ()
  end;
  if not !skip_example then begin
    Vertex_Test.example ();
    Parser_Test.example ();
    Model_Test.example ()
  end;
  exit 0
