(* circe2/circe2_tool.ml --  *)
(* Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
   Circe2 is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by 
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
   Circe2 is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)  

(* \subsubsection{Large Numeric File I/O} *)

type input_file =
  | ASCII_ic of in_channel
  | ASCII_inf of string
  | Binary_inf of string

type output_file =
  | ASCII_oc of out_channel
  | ASCII_outf of string
  | Binary_outf of string

let read columns = function
  | ASCII_ic ic -> Events.of_ascii_channel columns ic
  | ASCII_inf inf -> Events.of_ascii_file columns inf
  | Binary_inf inf -> Events.of_binary_file columns inf

let write output array =
  match output with
  | ASCII_oc oc -> Events.to_ascii_channel oc array
  | ASCII_outf outf -> Events.to_ascii_file outf array
  | Binary_outf outf -> Events.to_binary_file outf array

(* The special case of writing a binary file with mapped I/O
   can be treated most efficiently: *)
let cat columns input output =
  match input, output with
  | ASCII_ic ic, Binary_outf outf ->
      ignore (Events.of_ascii_channel ~file:outf columns ic)
  | _, _ -> write output (read columns input)

let map_xy fx fy columns input output =
  let a = read columns input in
  for i2 = 1 to Bigarray.Array2.dim2 a do
    Bigarray.Array2.set a 1 i2 (fx (Bigarray.Array2.get a 1 i2));
    Bigarray.Array2.set a 2 i2 (fy (Bigarray.Array2.get a 2 i2))
  done;
  write output a

let log10_xy = map_xy log10 log10
let exp10_xy = map_xy (fun x -> 10.0 ** x) (fun y -> 10.0 ** y)

(* \subsubsection{Histogramming} *)

let scan_string s =
  let tokens = Lexing.from_string s in
  let t1 = Events.next_float tokens in
  let t2 = Events.next_float tokens in
  let t3 = Events.next_float tokens in
  (t1, t2, t3)

let histogram_ascii name histograms =
  let ic = open_in name
  and histos =
    List.map (fun (tag, f, n, x_min, x_max) ->
      (tag, f, Histogram.create n x_min x_max)) histograms in
  begin try
    while true do
      let x, y, w = scan_string (input_line ic) in
      List.iter (fun (_, f, h) -> Histogram.record h (f x y) w) histos
    done
  with
  | End_of_file -> ()
  end;
  close_in ic;
  List.map (fun (t, _, h) -> (t, h)) histos

let histogram_binary_channel ic histograms =
  let histos =
    List.map (fun (tag, f, n, x_min, x_max) ->
      (tag, f, Histogram.create n x_min x_max)) histograms in
  begin try
    while true do
      let x = Float.Double.input_binary_float ic
      and y = Float.Double.input_binary_float ic
      and w = Float.Double.input_binary_float ic in
      List.iter (fun (_, f, h) -> Histogram.record h (f x y) w) histos
    done
  with
  | End_of_file -> ()
  end;
  List.map (fun (t, _, h) -> (t, h)) histos

let histogram_binary name histograms =
  let a = Events.of_binary_file 3 name
  and histos =
    List.map (fun (tag, f, n, x_min, x_max) ->
      (tag, f, Histogram.create n x_min x_max)) histograms in
  for i2 = 1 to Bigarray.Array2.dim2 a do
    let x = Bigarray.Array2.get a 1 i2
    and y = Bigarray.Array2.get a 2 i2
    and w = Bigarray.Array2.get a 3 i2 in
    List.iter (fun (_, f, h) -> Histogram.record h (f x y) w) histos
  done;
  List.map (fun (t, _, h) -> (t, h)) histos

(*i
let histogram_binary name histograms =
  let a = Events.of_binary_file 3 name
  and histos =
    List.map (fun (tag, f, n, x_min, x_max) ->
      (tag, f, Histogram.create n x_min x_max)) histograms in
  for i2 = 1 to Bigarray.Array2.dim2 a do
    let x = a.{1,i2}
    and y = a.{2,i2}
    and w = a.{3,i2} in
    List.iter (fun (_, f, h) -> Histogram.record h (f x y) w) histos
  done;
  List.map (fun (t, _, h) -> (t, h)) histos
i*)

let histogram_data to_file n reader suffix =
  let histograms = reader
      [ ("x", (fun x y -> x), n, 0.0, 1.0);
        ("x_low", (fun x y -> x), n, 0.0, 1.0e-4);
        ("1-x_low", (fun x y -> 1.0 -. x), n, 0.0, 1.0e-2);
        ("1-x_low2", (fun x y -> 1.0 -. x), n, 1.0e-10, 1.0e-2);
        ("y", (fun x y -> y), n, 0.0, 1.0);
        ("y_low", (fun x y -> y), n, 0.0, 1.0e-4);
        ("1-y_low", (fun x y -> 1.0 -. y), n, 0.0, 1.0e-2);
        ("1-y_low2", (fun x y -> 1.0 -. y), n, 1.0e-10, 1.0e-2);
        ("xy", (fun x y -> x *. y), n, 0.0, 1.0);
        ("xy_low", (fun x y -> x *. y), n, 0.0, 1.0e-8);
        ("z", (fun x y -> sqrt (x *. y)), n, 0.0, 1.0);
        ("z_low", (fun x y -> sqrt (x *. y)), n, 0.0, 1.0e-4);
        ("x-y", (fun x y -> x -. y), n, -1.0, 1.0);
        ("x_fine", (fun x y -> x), n, 0.75, 0.85);
        ("y_fine", (fun x y -> y), n, 0.75, 0.85);
        ("xy_fine", (fun x y -> x *. y), n, 0.5, 0.7);
        ("x-y_fine", (fun x y -> x -. y), n, -0.1, 0.1) ] in
  List.iter (fun (tag, h) ->
    to_file (tag ^ suffix) (Histogram.normalize h))
    histograms

(* \subsubsection{Moments} *)

let moments_ascii name moments =
  let ic = open_in name
  and f = Array.of_list (List.map (fun (tag, f) -> f) moments)
  and m = Array.of_list (List.map (fun (tag, f) -> 0.0) moments)
  and sum_w = ref 0.0 in
  begin try
    while true do
      let x, y, w = scan_string (input_line ic) in
      sum_w := !sum_w +. w;
      for i = 0 to Array.length f - 1 do
        m.(i) <- m.(i) +. w *. (f.(i) x y)
      done
    done
  with
  | End_of_file -> ()
  end;
  close_in ic;
  List.map2 (fun (tag, f) m -> (tag, m /. !sum_w)) moments (Array.to_list m)

let moments_binary name moments =
  let a = Events.of_binary_file 3 name in
  let f = Array.of_list (List.map (fun (tag, f) -> f) moments)
  and m = Array.of_list (List.map (fun (tag, f) -> 0.0) moments)
  and sum_w = ref 0.0 in
  for i2 = 1 to Bigarray.Array2.dim2 a do
    let x = Bigarray.Array2.get a 1 i2
    and y = Bigarray.Array2.get a 2 i2
    and w = Bigarray.Array2.get a 3 i2 in
    sum_w := !sum_w +. w;
    for i = 0 to Array.length f - 1 do
      m.(i) <- m.(i) +. w *. (f.(i) x y)
    done
  done;
  List.map2 (fun (tag, f) m -> (tag, m /. !sum_w)) moments (Array.to_list m)

let fmt var = function
  | 0 -> ""
  | 1 -> var
  | n -> var ^ "^" ^ string_of_int n

let moment nx ny =
  (fmt "x" nx ^ fmt "y" ny, (fun x y -> x ** (float nx) *. y ** (float ny)))
  
let diff_moment n =
  (fmt "|x-y|" n, (fun x y -> (abs_float (x -. y)) ** (float n)))
  
let moments_data reader =
  let moments = reader
      (List.map (moment 0) [1; 2; 3; 4; 5; 6] @
       List.map (moment 1) [0; 1; 2; 3; 4; 5] @
       List.map (moment 2) [0; 1; 2; 3; 4] @
       List.map (moment 3) [0; 1; 2; 3] @
       List.map (moment 4) [0; 1; 2] @
       List.map (moment 5) [0; 1] @
       List.map (moment 6) [0] @
       List.map diff_moment [1; 2; 3; 4; 5; 6]) in
  List.iter (fun (tag, m) -> Printf.printf "%s = %g\n" tag m) moments

(* \subsubsection{Regression} *)

let regression_interval (tag, h) (log_min, log_max) =
  let a, b =
    Histogram.regression h
      (fun x -> x >= log_min && x <= log_max) (fun x -> x) (fun x -> log x) in
  Printf.printf "%g<%s<%g: a = %g, b = %g\n" log_min tag log_max a b

let intervals =
  [ (-7.0, -6.0);
    (-6.0, -5.0);
    (-5.0, -4.0);
    (-4.0, -3.0);
    (-3.0, -2.0);
    (-7.0, -5.0);
    (-6.0, -4.0);
    (-5.0, -3.0);
    (-4.0, -2.0);
    (-7.0, -4.0);
    (-6.0, -3.0);
    (-5.0, -2.0);
    (-7.0, -3.0);
    (-6.0, -2.0) ]

let intervals =
  [ (-7.0, -4.0);
    (-6.0, -3.0);
    (-7.0, -3.0);
    (-6.0, -2.0) ]

let regression_data n reader =
  let histograms = reader
      [ ("log(x1)", (fun x1 x2 -> log x1), n, -8.0, 0.0);
        ("log(x2)", (fun x1 x2 -> log x2), n, -8.0, 0.0) ] in
  List.iter (fun (tag, h) ->
    List.iter (regression_interval (tag, h)) intervals) histograms

(* \subsubsection{Visually Adapting Powermaps} *)

let power_map beta eta =
  Diffmap.Power.create ~alpha:(1.0 /. (1.0 +. beta)) ~eta 0.0 1.0

let power_data to_file n center resolution reader suffix =
  let histograms = reader
      (List.flatten
         (List.map (fun p ->
           let pm = power_map p 0.0 in
           let ihp = Diffmap.Power.ihp pm in
           [((Printf.sprintf "1-x_low.%.2f" p), (fun x1 x2 -> ihp (1.0-.x1)), n, 0.0, ihp 1.0e-4);
            ((Printf.sprintf "1-y_low.%.2f" p), (fun x1 x2 -> ihp (1.0-.x2)), n, 0.0, ihp 1.0e-4);
            ((Printf.sprintf "x_low.%.2f" p), (fun x1 x2 -> ihp x1), n, 0.0, ihp 1.0e-4);
            ((Printf.sprintf "y_low.%.2f" p), (fun x1 x2 -> ihp x2), n, 0.0, ihp 1.0e-4)])
            [center -. 2.0 *. resolution;
             center -. resolution; center; center +. resolution;
             center +. 2.0 *. resolution])) in
  List.iter (fun (tag, h) ->
    to_file (tag ^ suffix) (Histogram.normalize h)) histograms

(* \subsubsection{Testing} *)

let make_test_data n (x_min, x_max) (y_min, y_max) f =
  let delta_x = x_max -. x_min
  and delta_y = y_max -. y_min in
  let array =
    Bigarray.Array2.create Bigarray.float64 Bigarray.fortran_layout 3 n in
  for i = 1 to n do
    let x = x_min +. Random.float delta_x
    and y = y_min +. Random.float delta_y in
    Bigarray.Array2.set array 1 i x;
    Bigarray.Array2.set array 2 i y;
    Bigarray.Array2.set array 3 i (f x y)
  done;
  array

(*i
let make_test_data n (x_min, x_max) (y_min, y_max) f =
  let delta_x = x_max -. x_min
  and delta_y = y_max -. y_min in
  let array =
    Bigarray.Array2.create Bigarray.float64 Bigarray.fortran_layout 3 n in
  for i = 1 to n do
    let x = x_min +. Random.float delta_x
    and y = y_min +. Random.float delta_y in
    array.{1,i} <- x;
    array.{2,i} <- y;
    array.{3,i} <- f x y
  done;
  array
i*)

module Div = Division.Mono
module Grid = Grid.Make (Div)

let test_design grid =
  let channel =
    { Grid.pid1 = 22; Grid.pol1 = 0;
      Grid.pid2 = 22; Grid.pol2 = 0;
      Grid.lumi = 0.0; Grid.g = grid } in
  { Grid.name = "TEST";
    Grid.roots = 500.0;
    Grid.channels = [ channel ];
    Grid.comments = [ "unphysical test" ]}

let test verbose triangle shrink nbins name f =
  let data = make_test_data 100000 (0.4, 0.9) (0.2, 0.7) f in
  let initial_grid =
    Grid.create ~triangle
      (Div.create nbins 0.0 1.0)
      (Div.create nbins 0.0 1.0) in
  let grid =
    Grid.of_bigarray ~verbose
      ~fixed_x1_min:(not shrink) ~fixed_x1_max:(not shrink)
      ~fixed_x2_min:(not shrink) ~fixed_x2_max:(not shrink)
      data initial_grid in
  Grid.designs_to_file name [test_design grid]

let random_interval () =
  let x1 = Random.float 1.0
  and x2 = Random.float 1.0 in
  (min x1 x2, max x1 x2)

module Test_Power = Diffmap.Make_Test (Diffmap.Power)
module Test_Resonance = Diffmap.Make_Test (Diffmap.Resonance)

let test_maps seed =
  Random.init seed;
  let x_min, x_max = random_interval ()
  and y_min, y_max = random_interval () in
  let alpha = 1.0 +. Random.float 4.0
  and eta =
    if Random.float 1.0 > 0.5 then
      y_max +. Random.float 5.0
    else
      y_min -. Random.float 5.0 in
  Test_Power.all
    (Diffmap.Power.create ~alpha ~eta ~x_min ~x_max y_min y_max);
  let a = Random.float 1.0
  and eta = y_min +. Random.float (y_max -. y_min) in
  Test_Resonance.all
    (Diffmap.Resonance.create ~eta ~a ~x_min ~x_max y_min y_max)
  
(* \subsubsection{Main Program} *)

type format = ASCII | Binary

type action =
  | Nothing
  | Command_file of string
  | Commands of string
  | Cat
  | Histo of format * string
  | Moments of format * string
  | Regression of format * string
  | Test of string * (float -> float -> float)
  | Test_Diffmaps of int
  | Unit_Tests
  | Log10
  | Exp10
  | Power of format * string

let rec passed = function
  | [] -> true
  | (OUnit.RFailure _ | OUnit.RError _ | OUnit.RTodo _ ) :: _ -> false
  | (OUnit.RSuccess _ | OUnit.RSkip _) :: tail -> passed tail

let _ =
  let usage = "usage: " ^ Sys.argv.(0) ^ " [options]" in
  let nbins = ref 100
  and triangle = ref false
  and shrink = ref false
  and verbose = ref false
  and action = ref Nothing
  and suffix = ref ".histo"
  and input = ref (ASCII_ic stdin)
  and output = ref (ASCII_oc stdout)
  and columns = ref 3
  and histogram_to_file = ref Histogram.to_file
  and center = ref 0.0
  and resolution = ref 0.01 in
  Arg.parse
    [("-c", Arg.String (fun s -> action := Commands s), "commands");
     ("-f", Arg.String (fun f -> action := Command_file f), "command file");
     ("-ia", Arg.String (fun n -> input := ASCII_inf n),
      "ASCII input file");
     ("-ib", Arg.String (fun n -> input := Binary_inf n),
      "Binary input file");
     ("-oa", Arg.String (fun n -> output := ASCII_outf n),
      "ASCII output file");
     ("-ob", Arg.String (fun n -> output := Binary_outf n),
      "Binary output file");
     ("-cat", Arg.Unit (fun () ->
       input := ASCII_ic stdin; output := ASCII_oc stdout;
       action := Cat), "copy stdin to stdout");
     ("-log10", Arg.Unit (fun () ->
       input := ASCII_ic stdin; output := ASCII_oc stdout;
       action := Log10), "");
     ("-exp10", Arg.Unit (fun () ->
       input := ASCII_ic stdin; output := ASCII_oc stdout;
       action := Exp10), "");
     ("-ha", Arg.String (fun s -> action := Histo (ASCII, s)),
      "ASCII histogramming tests");
     ("-hb", Arg.String (fun s -> action := Histo (Binary, s)),
      "binary histogramming tests");
     ("-ma", Arg.String (fun s -> action := Moments (ASCII, s)),
      "ASCII moments  tests");
     ("-mb", Arg.String (fun s -> action := Moments (Binary, s)),
      "binary moments tests");
     ("-pa", Arg.String (fun s -> action := Power (ASCII, s)), "");
     ("-pb", Arg.String (fun s -> action := Power (Binary, s)), "");
     ("-C", Arg.Float (fun c -> center := c), "");
     ("-R", Arg.Float (fun r -> resolution := r), "");
     ("-Pa", Arg.String (fun s -> action := Regression (ASCII, s)), "");
     ("-Pb", Arg.String (fun s -> action := Regression (Binary, s)), "");
     ("-p", Arg.String (fun s -> suffix := s), "histogram name suffix");
     ("-h", Arg.Unit (fun () ->
       histogram_to_file := Histogram.as_bins_to_file), "");
     ("-b", Arg.Int (fun n -> nbins := n), "#bins");
     ("-s", Arg.Set shrink, "shrinkwrap interval");
     ("-S", Arg.Clear shrink, "don't shrinkwrap interval [default]");
     ("-t", Arg.Set triangle,
      "project symmetrical distribution onto triangle");
     ("-v", Arg.Set verbose, "verbose");
     ("-test", Arg.Unit (fun () -> action := Unit_Tests),
      "run unit test suite");
     ("-test1", Arg.String (fun s ->
       action := Test (s, fun x y -> 1.0)), "testing");
     ("-test2", Arg.String (fun s ->
       action := Test (s, fun x y -> x *. y)), "testing");
     ("-test3", Arg.String (fun s ->
       action := Test (s, fun x y -> 1.0 /. x +. 1.0 /. y)), "testing");
     ("-testm", Arg.Int (fun seed -> action := Test_Diffmaps seed),
      "testing maps") ]
    (fun names -> prerr_endline usage; exit 2)
    usage;
  begin try
    match !action with
    | Nothing -> ()
    | Commands name -> Commands.execute (Commands.parse_string name)
    | Command_file name -> Commands.execute (Commands.parse_file name)
    | Histo (ASCII, name) ->
        histogram_data !histogram_to_file !nbins
          (histogram_ascii name) !suffix
    | Histo (Binary, "-") ->
        histogram_data !histogram_to_file !nbins
          (histogram_binary_channel stdin) !suffix
    | Histo (Binary, name) ->
        histogram_data !histogram_to_file !nbins
          (histogram_binary name) !suffix
    | Moments (ASCII, name) -> moments_data (moments_ascii name)
    | Moments (Binary, name) -> moments_data (moments_binary name)
    | Power (ASCII, name) ->
        power_data !histogram_to_file !nbins !center !resolution
          (histogram_ascii name) !suffix
    | Power (Binary, name) ->
        power_data !histogram_to_file !nbins !center !resolution
          (histogram_binary name) !suffix
    | Regression (ASCII, name) -> regression_data !nbins (histogram_ascii name)
    | Regression (Binary, name) -> regression_data !nbins (histogram_binary name)
    | Cat -> cat !columns !input !output 
    | Log10 -> log10_xy !columns !input !output 
    | Exp10 -> exp10_xy !columns !input !output 
    | Test (name, f) -> test !verbose !triangle !shrink !nbins name f
    | Test_Diffmaps seed -> test_maps seed
    | Unit_Tests ->
        let suite =
          OUnit.(>:::) "All"
            [ThoArray.suite;
             ThoMatrix.suite;
             Filter.suite] in
        if passed (OUnit.run_test_tt ~verbose:!verbose suite) then
          exit 0
        else
          exit 1

  with
  | Syntax.Syntax_Error (msg, _, _) ->
      Printf.eprintf "%s: parse error: %s\n" Sys.argv.(0) msg;
      exit 1
  end;
  exit 0

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
