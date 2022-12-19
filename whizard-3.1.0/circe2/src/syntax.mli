(* circe2/syntax.mli --  *)
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

(* \subsection{Abstract Syntax and Default Values} *)

val epsilon : float

(* A channel is uniquely specified by PDG particle ids and
   polarizations $\{-1,0,+1\}$, which must match the `events'
   in the given file; as should the luminosity.  The options
   are for tuning the grid. *)

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

(* A parameter set is uniquely specified by PDG particle ids
   (\emph{par abus de langage}), polarizations (now a floating
   point number for the effective polarization of the beam), and
   center of mass energy.  This must match the `events' in the
   files given for the channels.  The other options are for tuning
   the grid. *)

type design =
    { design : string;
      roots : float;
      design_bins1 : int;
      design_bins2 : int;
      design_scale1 : float option;
      design_scale2 : float option;
      channels : channel list;
      comments : string list }

val default_design : design

val default_channel : design -> channel

(* One file can hold more than one grid. *)

type file = { name : string; designs : design list }

val default_file : file

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
