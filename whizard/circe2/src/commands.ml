(* circe2/commands.ml --  *)
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

exception Invalid_interval of float * float

type t = Syntax.t

open Printf

module Maps = Diffmaps.Default
module Div = Division.Make_Poly (Maps)
module Grid = Grid.Make (Div)

(* \subsubsection{Processing} *)

let smooth_grid channel grid =    
  List.fold_left
    (fun acc (width, area) -> Grid.smooth width area acc)
    grid channel.Syntax.smooth

let report msg =
  prerr_string msg;
  flush stderr

let process_channel ch =
  report ("reading: " ^ ch.Syntax.events ^ " ...");
  let data =
    if ch.Syntax.binary then
      Events.of_binary_file ch.Syntax.columns ch.Syntax.events
    else
      Events.of_ascii_file ch.Syntax.columns ch.Syntax.events in
  report " done.\n";
  begin match ch.Syntax.scale1, ch.Syntax.scale2 with
  | None, None -> ()
  | Some scale1, None -> Events.rescale scale1 1.0 data
  | None, Some scale2 -> Events.rescale 1.0 scale2 data
  | Some scale1, Some scale2 -> Events.rescale scale1 scale2 data
  end;
  let initial_grid =
    Grid.create ~triangle:ch.Syntax.triangle
      (Div.create ch.Syntax.intervals1 ch.Syntax.bins1
         ch.Syntax.x1_min ch.Syntax.x1_max)
      (Div.create ch.Syntax.intervals2 ch.Syntax.bins2
         ch.Syntax.x2_min ch.Syntax.x2_max) in
  let grid =
    Grid.of_bigarray ~verbose:true
      ~iterations:ch.Syntax.iterations
      ~fixed_x1_min:ch.Syntax.fixed_x1_min
      ~fixed_x1_max:ch.Syntax.fixed_x1_max
      ~fixed_x2_min:ch.Syntax.fixed_x2_min
      ~fixed_x2_max:ch.Syntax.fixed_x2_max
      ~areas:(List.map snd ch.Syntax.smooth)
      data initial_grid in
  let smoothed_grid = smooth_grid ch grid in
  begin match ch.Syntax.histogram with
  | Some name ->
      let oc = open_out name in
      Grid.to_channel_2d oc smoothed_grid;
(*i   List.iter
        (fun (_, area) ->
          report (Printf.sprintf " %g>%g>%g"
                    (Grid.variance_area area initial_grid)
                    (Grid.variance_area area grid)
                    (Grid.variance_area area smoothed_grid)))
        ch.Syntax.smooth;
      report " ";  *)
      close_out oc
  | None -> ()
  end;
  { Grid.pid1 = ch.Syntax.pid1;
    Grid.pol1 = ch.Syntax.pol1;
    Grid.pid2 = ch.Syntax.pid2;
    Grid.pol2 = ch.Syntax.pol2;
    Grid.lumi = ch.Syntax.lumi;
    Grid.g = smoothed_grid }

module S = Set.Make (struct type t = string let compare = compare end)

let channel_prerequisites acc ch =
  S.add ch.Syntax.events acc

let process_design oc name design =
  let channels = List.rev_map process_channel design.Syntax.channels
  and comments = List.rev design.Syntax.comments in
  let acc =
    { Grid.name = design.Syntax.design;
      Grid.roots = design.Syntax.roots;
      Grid.channels = channels;
      Grid.comments = comments } in
  report ("writing: " ^ name ^ " ...");
  Grid.design_to_channel oc acc;
  report " done.\n"

  
let design_prerequisites acc design =
  List.fold_left (channel_prerequisites) acc design.Syntax.channels

let write_file file =
  let oc = open_out file.Syntax.name in
  List.iter (process_design oc file.Syntax.name) file.Syntax.designs;
  close_out oc

let file_prerequisites acc file =
  List.fold_left (design_prerequisites) acc file.Syntax.designs

let prerequisites files =
  List.fold_left (file_prerequisites) S.empty files

let unreadable name =
  try
    Unix.access name [Unix.R_OK];
    false
  with
  | Unix.Unix_error (_, _, _) -> true

let execute files =
  let missing = S.filter unreadable (prerequisites files) in
  if S.is_empty missing then
    List.iter write_file files
  else
    eprintf "circe2_tool: unreadable input files: %s!\n"
      (String.concat ", " (S.elements missing))

(* \subsubsection{Translate.} *)

let rec update_fix acc = function
  | b, Syntax.X12, s -> update_fix (update_fix acc (b, Syntax.X2, s)) (b, Syntax.X1, s)
  | b, c, Syntax.Minmax -> update_fix (update_fix acc (b, c, Syntax.Max)) (b, c, Syntax.Min)
  | b, Syntax.X1, Syntax.Min -> { acc with Syntax.fixed_x1_min = b }
  | b, Syntax.X1, Syntax.Max -> { acc with Syntax.fixed_x1_max = b }
  | b, Syntax.X2, Syntax.Min -> { acc with Syntax.fixed_x2_min = b }
  | b, Syntax.X2, Syntax.Max -> { acc with Syntax.fixed_x2_max = b }

let rec update_pid acc = function
  | n, Syntax.X12 -> update_pid (update_pid acc (n, Syntax.X2)) (n, Syntax.X1)
  | n, Syntax.X1 -> { acc with Syntax.pid1 = n }
  | n, Syntax.X2 -> { acc with Syntax.pid2 = n }

let rec update_pol acc = function
  | n, Syntax.X12 -> update_pol (update_pol acc (n, Syntax.X2)) (n, Syntax.X1)
  | n, Syntax.X1 -> { acc with Syntax.pol1 = n }
  | n, Syntax.X2 -> { acc with Syntax.pol2 = n }

let rec update_bins acc = function
  | n, Syntax.X12 -> update_bins (update_bins acc (n, Syntax.X2)) (n, Syntax.X1)
  | n, Syntax.X1 -> { acc with Syntax.bins1 = n }
  | n, Syntax.X2 -> { acc with Syntax.bins2 = n }

let rec update_scale acc = function
  | x, Syntax.X12 -> update_scale (update_scale acc (x, Syntax.X2)) (x, Syntax.X1)
  | x, Syntax.X1 -> { acc with Syntax.scale1 = Some x }
  | x, Syntax.X2 -> { acc with Syntax.scale2 = Some x }

let rec update_x_min acc = function
  | x, Syntax.X12 -> update_x_min (update_x_min acc (x, Syntax.X2)) (x, Syntax.X1)
  | x, Syntax.X1 -> { acc with Syntax.x1_min = x }
  | x, Syntax.X2 -> { acc with Syntax.x2_min = x }

let rec update_x_max acc = function
  | x, Syntax.X12 -> update_x_max (update_x_max acc (x, Syntax.X2)) (x, Syntax.X1)
  | x, Syntax.X1 -> { acc with Syntax.x1_max = x }
  | x, Syntax.X2 -> { acc with Syntax.x2_max = x }

let rec update_map acc = function
  | m, Syntax.X12 -> update_map (update_map acc (m, Syntax.X2)) (m, Syntax.X1)
  | m, Syntax.X1 -> { acc with Syntax.intervals1 = m :: acc.Syntax.intervals1 }
  | m, Syntax.X2 -> { acc with Syntax.intervals2 = m :: acc.Syntax.intervals2 }

let update_smooth acc width area =
  { acc with Syntax.smooth = (width, area) :: acc.Syntax.smooth }

let channel design cmds =
  List.fold_left
    (fun acc -> function
      | Syntax.Pid (n, c) -> update_pid acc (n, c)
      | Syntax.Pol (p, c) -> update_pol acc (p, c)
      | Syntax.Lumi l -> { acc with Syntax.lumi = l }
      | Syntax.Diffmap (m, c) -> update_map acc (m, c)
      | Syntax.Bins (n, c) -> update_bins acc (n, c)
      | Syntax.Scale (x, c) -> update_scale acc (x, c)
      | Syntax.Xmin (x, c) -> update_x_min acc (x, c)
      | Syntax.Xmax (x, c) -> update_x_max acc (x, c)
      | Syntax.Fix (b, c, s) -> update_fix acc (b, c, s)
      | Syntax.Smooth (w, a) -> update_smooth acc w a
      | Syntax.Triangle b -> { acc with Syntax.triangle = b }
      | Syntax.Iterations i -> { acc with Syntax.iterations = i }
      | Syntax.Events s -> { acc with Syntax.events = s }
      | Syntax.Histogram s -> { acc with Syntax.histogram = Some s }
      | Syntax.Columns n ->
          if n < 2 then
            invalid_arg "#columns < 2"
          else
            { acc with Syntax.columns = n }
      | Syntax.Binary b -> { acc with Syntax.binary = b })
    (Syntax.default_channel design) cmds

let rec update_design_bins acc = function
  | n, Syntax.X12 ->
      update_design_bins (update_design_bins acc (n, Syntax.X2)) (n, Syntax.X1)
  | n, Syntax.X1 ->
      { acc with Syntax.design_bins1 = n }
  | n, Syntax.X2 ->
      { acc with Syntax.design_bins2 = n }

let rec update_design_scale acc = function
  | x, Syntax.X12 ->
      update_design_scale (update_design_scale acc (x, Syntax.X2)) (x, Syntax.X1)
  | x, Syntax.X1 ->
      { acc with Syntax.design_scale1 = Some x }
  | x, Syntax.X2 ->
      { acc with Syntax.design_scale2 = Some x }

let design cmds =
  List.fold_left
    (fun acc -> function
      | Syntax.Design s ->
          { acc with Syntax.design = s }
      | Syntax.Roots r ->
          { acc with Syntax.roots = r }
      | Syntax.Design_Bins (n, c) ->
          update_design_bins acc (n, c)
      | Syntax.Design_Scale (x, c) ->
          update_design_scale acc (x, c)
      | Syntax.Channels cmds ->
	  { acc with Syntax.channels = channel acc cmds :: acc.Syntax.channels }
      | Syntax.Comment c ->
          { acc with Syntax.comments = c :: acc.Syntax.comments })
    Syntax.default_design cmds

let file cmds =
  List.fold_right
    (fun cmd acc ->
      match cmd with
      | Syntax.File s ->
          { acc with Syntax.name = s }
      | Syntax.Designs cmds ->
          { acc with Syntax.designs = design cmds :: acc.Syntax.designs })
    cmds Syntax.default_file 

(* \subsubsection{API} *)

let parse_file name =
  let ic = open_in name in
  let cmds =
    List.map file (Parser.main Lexer.token (Lexing.from_channel ic)) in
  close_in ic;
  cmds

let parse_string s =
  List.map file (Parser.main Lexer.token (Lexing.from_string s))

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
