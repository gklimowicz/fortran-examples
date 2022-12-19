(* circe2/grid.mli --  *)
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

exception Out_of_range of string * float * (float * float)

module type T =
  sig
    module D : Division.T

    type t
    val copy : t -> t

    (* Create an initial grid. *)
    val create : ?triangle:bool -> D.t -> D.t -> t

    (* [record grid x1 x2 w] records the value~[w] in the
       bin corresponding to coordinates~[x1] and~[x2]. *)
    val record : t -> float -> float -> float -> unit

    (* VEGAS style rebinning. *)
    val rebin : ?power:float ->
      ?fixed_x1_min:bool -> ?fixed_x1_max:bool ->
        ?fixed_x2_min:bool -> ?fixed_x2_max:bool -> t -> t

    (* The sum of all the weights shall be one. *)
    val normalize : t -> t

    (* Adapt an initial grid to data. The [power] controls speed
       vs.~stability of adaption and is passed on to [Division.rebin].
       [iterations] provides a hard cutoff for the number of
       iterations (default: 1000), while  [margin] and [cutoff] control
       the soft cutoff of the adaption.  If the variance grows to the
       best value multiplied by [margin] of if there are no improvements
       for [cutoff] steps, the adaption is stopped (defaults:~$1.5$
       and~$20$).  The remaining options control if the boundaries
       are fixed or allowed to move towards the limits of the dataset.
       The defaults are all [false], meaning that the boundaries are
       allowed to move. *)
    val of_bigarray : ?verbose:bool -> ?power:float ->
      ?iterations:int -> ?margin:float -> ?cutoff:int ->
        ?fixed_x1_min:bool -> ?fixed_x1_max:bool ->
          ?fixed_x2_min:bool -> ?fixed_x2_max:bool ->
            ?areas:Syntax.area list ->
              (float, Bigarray.float64_elt,
               Bigarray.fortran_layout) Bigarray.Array2.t -> t -> t

    val smooth : float -> Syntax.area -> t -> t

    val variance_area : Syntax.area -> t -> float

    val to_channel_2d : out_channel -> t -> unit

    (* Write output that \KirkeTwo/ can read: *)
    type channel =
        { pid1 : int;
          pol1 : int;
          pid2 : int;
          pol2 : int;
          lumi : float;
          g : t }

    val to_channel : out_channel -> channel -> unit

    type design =
        { name : string;
          roots : float;
          channels : channel list;
          comments : string list }

    val design_to_channel : out_channel -> design -> unit
    val designs_to_channel : out_channel ->
      ?comments:string list-> design list -> unit
    val designs_to_file : string ->
      ?comments:string list -> design list -> unit

    val variance : t -> float

  end

module Make (D : Division.T) : T with module D = D

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
