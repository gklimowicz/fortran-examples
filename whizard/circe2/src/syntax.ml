(* circe2/syntax.ml --  *)
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

exception Syntax_Error of string * int * int

let epsilon = 100. *. epsilon_float

type boundary =
| Closed of float
| Open of float
| Bin of int

type point =
| Delta of float
| Box of int

type interval = boundary * boundary

type area =
| Rect of interval * interval
| Slice1 of interval * point
| Slice2 of point * interval

type channel =
    { pid1 : int;
      pol1 : int;
      pid2 : int;
      pol2 : int;
      lumi : float;
      bins1 : int;
      scale1 : float option;
      x1_min : float;
      x1_max : float;
      fixed_x1_min : bool;
      fixed_x1_max : bool;
      intervals1 : (int * Diffmaps.Default.t) list;
      bins2 : int;
      scale2 : float option;
      x2_min : float;
      x2_max : float;
      fixed_x2_min : bool;
      fixed_x2_max : bool;
      intervals2 : (int * Diffmaps.Default.t) list;
      smooth : (float * area) list;
      triangle : bool;
      iterations : int;
      events : string;
      histogram : string option;
      binary : bool;
      columns : int } 

 
type design =
    { design : string;
      roots : float;
      design_bins1 : int;
      design_bins2 : int;
      design_scale1 : float option;
      design_scale2 : float option;
      channels : channel list;
      comments : string list }

let default_design =
    { design = "TESLA";
      roots = 500.0;
      design_bins1 = 20;
      design_bins2 = 20;
      design_scale1 = None;
      design_scale2 = None;
      channels = [];
      comments = [] }

let default_channel design =
  { pid1 = 11 (* $e^-$ *);
    pol1 = 0;
    pid2 = -11 (* $e^+$ *);
    pol2 = 0;
    lumi = 0.0;
    bins1 = design.design_bins1;
    scale1 = design.design_scale1;
    x1_min = 0.0;
    x1_max = 1.0;
    fixed_x1_min = false;
    fixed_x1_max = false;
    intervals1 = [];
    bins2 = design.design_bins2;
    scale2 = design.design_scale2;
    x2_min = 0.0;
    x2_max = 1.0;
    fixed_x2_min = false;
    fixed_x2_max = false;
    intervals2 = [];
    smooth = [];
    triangle = false;
    iterations = 1000;
    events = "circe2.events";
    histogram = None;
    binary = false;
    columns = 3 } 


type file = { name : string; designs : design list }

let default_file = { name = "circe2_tool.out"; designs = [] }

type t = file list

type coord = X1 | X2 | X12
type side = Min | Max | Minmax

type channel_cmd =
  | Pid of int * coord
  | Pol of int * coord
  | Lumi of float
  | Xmin of float * coord
  | Xmax of float * coord
  | Bins of int * coord
  | Scale of float * coord
  | Diffmap of (int * Diffmaps.Default.t) * coord
  | Smooth of float * area
  | Triangle of bool
  | Iterations of int
  | Events of string
  | Histogram of string
  | Binary of bool
  | Columns of int
  | Fix of bool * coord * side

type design_cmd =
  | Design of string
  | Roots of float
  | Design_Bins of int * coord
  | Design_Scale of float * coord
  | Channels of channel_cmd list
  | Comment of string

type file_cmd =
  | File of string
  | Designs of design_cmd list
    
type file_cmds = file_cmd list

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
