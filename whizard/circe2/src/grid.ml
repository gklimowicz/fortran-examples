(* circe2/grid.ml --  *)
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

open Printf

module type T =
  sig
    module D : Division.T

    type t
    val copy : t -> t
    val create : ?triangle:bool -> D.t -> D.t -> t

    val record : t -> float -> float -> float -> unit

    val rebin : ?power:float ->
      ?fixed_x1_min:bool -> ?fixed_x1_max:bool ->
        ?fixed_x2_min:bool -> ?fixed_x2_max:bool -> t -> t

    val normalize : t -> t

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

module Make (D : Division.T) =
  struct
    module D = D

    type t =
        { d1 : D.t;
          d2 : D.t;
          w : float array array;
          var : float array array;
          triangle : bool }

    let copy grid =
      { d1 = D.copy grid.d1;
        d2 = D.copy grid.d2;
        w = ThoMatrix.copy grid.w;
        var = ThoMatrix.copy grid.var;
        triangle = grid.triangle }

    let create ?(triangle = false) d1 d2 =
      let n1 = D.n_bins d1
      and n2 = D.n_bins d2 in
      { d1 = d1;
        d2 = d2;
        w = Array.make_matrix n1 n2 0.0;
        var = Array.make_matrix n1 n2 0.0;
        triangle = triangle }

    (* We need
       \begin{subequations}
       \begin{align}
         \textit{upper}\; x\rbrack &= \textit{lower}\; \lbrack x \\
         \textit{upper}\; x\rbrack &= \textit{lower}\; (x - 1 \\
         \textit{upper}\; x) &= \textit{lower}\; \lbrack x - 1 \\
         \textit{upper}\; x) &= \textit{lower}\; (x - 2
       \end{align}
       \end{subequations}
       and
       \begin{subequations}
       \begin{align}
         \textit{upper}\; x\rbrack &= \textit{upper}\; x) + 1 \\
         \textit{lower}\; \lbrack x &= \textit{lower}\; (x - 1
       \end{align}
       \end{subequations} *)

    (* [lower_bin] had [Open] and [Closed] mixed up! (tho:2014-12-09) *)

    let lower_bin div limit =
      try
        begin match limit with
        | Syntax.Closed x -> D.find div x
        | Syntax.Open x -> D.find div x + 1
        | Syntax.Bin n -> n
        end
      with
      | Division.Below_min (_, _, n) -> n
      | Division.Above_max (x, range, _) ->
          raise (Out_of_range ("Grid.lower_bin", x, range))

    let upper_bin div limit =
      try
        begin match limit with
        | Syntax.Closed x -> D.find div x
        | Syntax.Open x -> D.find div x - 1
        | Syntax.Bin n -> n
        end
      with
      | Division.Above_max (_, _, n) -> n
      | Division.Below_min (x, range, _) ->
          raise (Out_of_range ("Grid.upper_bin", x, range))

    let enclosed_bins div (x1, x2) =
      (lower_bin div x1, upper_bin div x2)

    let enclosing_bin div = function
      | Syntax.Delta x -> D.find div x
      | Syntax.Box n -> n

    let smooth width area grid =
      let gaussian = Filter.gaussian width in
      let w =
        begin match area with
        | Syntax.Rect (i1, i2) ->
            let nx1, nx2 = enclosed_bins grid.d1 i1
            and ny1, ny2 = enclosed_bins grid.d2 i2 in
            Filter.apply12
              ~inf1:nx1 ~sup1:nx2 ~inf2:ny1 ~sup2:ny2
              gaussian gaussian grid.w
        | Syntax.Slice1 (i1, y) ->
            let nx1, nx2 = enclosed_bins grid.d1 i1
            and ny = enclosing_bin grid.d2 y in
            Filter.apply1
                ~inf1:nx1 ~sup1:nx2 ~inf2:ny ~sup2:ny
                gaussian grid.w
        | Syntax.Slice2 (x, i2) ->
            let nx = enclosing_bin grid.d1 x
            and ny1, ny2 = enclosed_bins grid.d2 i2 in
            Filter.apply2
              ~inf1:nx ~sup1:nx ~inf2:ny1 ~sup2:ny2
              gaussian grid.w
        end in
      { grid with w }

    let to_channel_2d oc grid =
      for i = 0 to D.n_bins grid.d1 - 1 do
        Printf.fprintf oc "%g" grid.w.(i).(0);
        for j = 1 to D.n_bins grid.d2 - 1 do
          Printf.fprintf oc " %g" grid.w.(i).(j)
        done;
        Printf.fprintf oc "\n"
      done
      
    let project_triangle triangle x y =
      if triangle then begin
        if x >= y then begin
          (x, y /. x)
        end else begin
          (y, x /. y)
        end
      end else
        (x, y)

    (* Note that there is \emph{no} jacobian here.  It is
       applied later by the Fortran program interpreting the grid as
       a distribution.  It is not needed for the event generator
       anyway. *)
    let record grid x y f =
      let x', y' = project_triangle grid.triangle x y in
      D.record grid.d1 x' f;
      D.record grid.d2 y' f;
      let n1 = D.find grid.d1 x'
      and n2 = D.find grid.d2 y' in
      grid.w.(n1).(n2) <- grid.w.(n1).(n2) +. f;
      grid.var.(n1).(n2) <- grid.var.(n1).(n2)
          +. f /. D.caj grid.d1 x' /. D.caj grid.d2 y'

    let rebin ?power ?fixed_x1_min ?fixed_x1_max
        ?fixed_x2_min ?fixed_x2_max grid =
      let n1 = D.n_bins grid.d1
      and n2 = D.n_bins grid.d2 in
      { d1 = D.rebin ?power
          ?fixed_min:fixed_x1_min ?fixed_max:fixed_x1_max grid.d1;
        d2 = D.rebin ?power
          ?fixed_min:fixed_x2_min ?fixed_max:fixed_x2_max grid.d2;
        w = Array.make_matrix n1 n2 0.0;
        var = Array.make_matrix n1 n2 0.0;
        triangle = grid.triangle }

    let normalize grid =
      let sum_w = ThoMatrix.sum_float grid.w in
      { d1 = D.copy grid.d1;
        d2 = D.copy grid.d2;
        w = ThoMatrix.map (fun w -> w /. sum_w) grid.w;
        var = ThoMatrix.copy grid.var;
        triangle = grid.triangle }

    (* Monitoring the variance in each cell is \emph{not} a good idea for
       approximating distributions of unweighted events: it always vanishes
       for unweighted events, even if they are distributed very unevenly.
       Therefore, we monitor the \emph{global} variance instead: *)

    let variance_area area grid =
      let (nx1, nx2), (ny1, ny2) =
        begin match area with
        | Syntax.Rect (i1, i2) ->
            (enclosed_bins grid.d1 i1, enclosed_bins grid.d2 i2)
        | Syntax.Slice1 (i1, y) ->
            let ny = enclosing_bin grid.d2 y in
            (enclosed_bins grid.d1 i1, (ny, ny))
        | Syntax.Slice2 (x, i2) ->
            let nx = enclosing_bin grid.d1 x in
            ((nx, nx), enclosed_bins grid.d2 i2)
        end in
      let n = float ((nx2 - nx1 + 1) * (ny2 - ny1 + 1)) in
      let w =
        ThoMatrix.sum_float
          ~inf1:nx1 ~sup1:nx2 ~inf2:ny1 ~sup2:ny2 grid.w /. n
      and w2 =
        ThoMatrix.fold_left
          ~inf1:nx1 ~sup1:nx2 ~inf2:ny1 ~sup2:ny2
          (fun acc w -> acc +. w *. w) 0.0 grid.w /. n in
      w2 -. w *. w

    let variance grid =
      let n = float (D.n_bins grid.d1 * D.n_bins grid.d2) in
      let w = ThoMatrix.sum_float grid.w /. n
      and w2 =
        ThoMatrix.fold_left (fun acc w -> acc +. w *. w) 0.0 grid.w /. n in
      w2 -. w *. w

    (* Find the grid with the lowest variance.  Allow local fluctuations and
       stop only after moving to twice the lowest value. *)

    let start_progress_report verbose var =
      if verbose then begin
        eprintf "adapting variance: %g" var;
        flush stderr
      end

    let progress_report verbose soft_limit best_var var =
      if verbose then begin
        if var < best_var then begin
          eprintf ", %g" var;
          flush stderr
        end else begin
          eprintf " [%d]" soft_limit;
          flush stderr
        end
      end

    let stop_progress_report verbose =
      if verbose then begin
        eprintf " done.\n";
        flush stderr
      end

(* Scan a bigarray.  Assume a uniform weight,
   if it has only 2 columns. *)

    let record_data data grid =
      let columns = Bigarray.Array2.dim1 data in
      if columns < 2 then
        eprintf "error: not enough columns"
      else
        for i2 = 1 to Bigarray.Array2.dim2 data do
          let x = Bigarray.Array2.get data 1 i2
          and y = Bigarray.Array2.get data 2 i2
          and w =
            if columns > 2 then
              Bigarray.Array2.get data 3 i2
            else
              1.0 in
          try
            record grid x y w
          with
          | Division.Out_of_range (x, (x_min, x_max)) ->
              eprintf "internal error: %g not in [%g,%g]\n" x x_min x_max
      done

(* The main routine constructing an adapted grid. *)

    let of_bigarray ?(verbose = false)
        ?power ?(iterations = 1000) ?(margin = 1.5) ?(cutoff = 10)
        ?fixed_x1_min ?fixed_x1_max ?fixed_x2_min ?fixed_x2_max
        ?areas data initial =

      let rebinner grid =
        rebin ?power
          ?fixed_x1_min ?fixed_x1_max ?fixed_x2_min ?fixed_x2_max grid in

      let rec improve_bigarray hard_limit soft_limit best_var best_grid grid =
        if soft_limit <= 0 || hard_limit <= 0 then
          normalize best_grid
        else begin
          record_data data grid;
          let var = variance grid in
          begin match areas with
          | None | Some [] -> ()
          | Some areas ->
              let normalized_grid = normalize grid in
              let variances =
                List.map
                  (fun area ->
                    variance_area area normalized_grid) areas in
              let msg =
                " (" ^ Printf.sprintf "%g" (variance normalized_grid) ^ ": " ^
                String.concat "; "
                  (List.map (fun x -> Printf.sprintf "%g" x) variances) ^
                ")" in
              prerr_string msg;
              flush stderr
          end;
          progress_report verbose soft_limit best_var var;
          if var >= margin *. best_var then
            normalize best_grid
          else
            let best_var, best_grid, soft_limit =
              if var < best_var then
                (var, grid, cutoff)
              else
                (best_var, best_grid, pred soft_limit) in

            (* Continuation passing makes recursion with exception handling
               tail recursive.  This is not really needed, because the
               data structures are not to big and recursion is not expected
               to be too deep. It doesn't hurt either, since the idiom is
               sufficiently transparent. *)
            let continue =
              try
                let grid' = rebinner grid in
                fun () -> improve_bigarray
                    (pred hard_limit) soft_limit best_var best_grid grid'
              with
              | Division.Rebinning_failure msg ->
                  eprintf "circe2: rebinning failed: %s!\n" msg;
                  fun () -> best_grid in
            continue ()
        end in

      record_data data initial;
      let var = variance initial in
      start_progress_report verbose var;

      let result =
        improve_bigarray iterations cutoff var initial (rebinner initial) in
      stop_progress_report verbose;
      result

    type channel =
        { pid1 : int;
          pol1 : int;
          pid2 : int;
          pol2 : int;
          lumi : float;
          g : t }

    (* NB: we need to transpose the weight matrix to
       get from our row major to \texttt{Fortran}'s
       column major array format expected by
       \texttt{circe2}! *)
    let to_channel oc ch =
      fprintf oc "pid1, pol1, pid2, pol2, lumi\n";
      fprintf oc " %d %d %d %d %G\n"
        ch.pid1 ch.pol1 ch.pid2 ch.pol2 ch.lumi;
      fprintf oc "#bins1, #bins2, triangle?\n";
      fprintf oc " %d %d %s\n"
        (D.n_bins ch.g.d1) (D.n_bins ch.g.d2)
        (if ch.g.triangle then "T" else "F");
      fprintf oc "x1, map1, alpha1, xi1, eta1, a1, b1\n";
      D.to_channel oc ch.g.d1;
      fprintf oc "x2, map2, alpha2, xi2, eta2, a2, b2\n";
      D.to_channel oc ch.g.d2;
      fprintf oc "weights\n";
      ThoMatrix.iter
        (fun x -> fprintf oc " %s\n" (Float.Double.to_string x))
        (ThoMatrix.transpose ch.g.w)

    type design =
        { name : string;
          roots : float;
          channels : channel list;
          comments : string list }

    type polarization_support =
      | Averaged
      | Helicities
      | Density_Matrices

    let polarization_support design =
      if List.for_all (fun ch -> ch.pol1 = 0 && ch.pol2 = 0)
          design.channels then
        Averaged
      else if List.for_all (fun ch -> ch.pol1 <> 0 && ch.pol2 <> 0)
          design.channels then
        Helicities
      else
        invalid_arg
          "Grid.polarization_support: mixed polarization support!"

    let format_polarization_support = function
      | Averaged -> "averaged"
      | Helicities -> "helicities"
      | Density_Matrices -> "density matrices"

    let getlogin () =
      (Unix.getpwuid (Unix.getuid ())).Unix.pw_name

    let design_to_channel oc design =
      let utc = Unix.gmtime (Unix.time ()) in
      List.iter (fun s -> fprintf oc "! %s\n" s) design.comments;
      fprintf oc "! generated with %s by %s@%s, "
        (Sys.argv.(0)) (getlogin ()) (Unix.gethostname ());
      fprintf oc "%4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d GMT\n"
        (utc.Unix.tm_year + 1900) (utc.Unix.tm_mon + 1) utc.Unix.tm_mday
        utc.Unix.tm_hour utc.Unix.tm_min utc.Unix.tm_sec;
      fprintf oc "CIRCE2 FORMAT#1\n";
      fprintf oc "design, roots\n";
      fprintf oc " '%s' %G\n" design.name design.roots;
      fprintf oc "#channels, pol.support\n";
      fprintf oc " %d '%s'\n"
        (List.length design.channels)
        (format_polarization_support (polarization_support design));
      List.iter (to_channel oc) design.channels;
      fprintf oc "ECRIC2\n"

    let designs_to_channel oc ?(comments = []) designs =
      List.iter (fun c -> fprintf oc "! %s\n" c) comments;
      List.iter (design_to_channel oc) designs

    let designs_to_file name ?comments designs =
      let oc = open_out name in
      designs_to_channel oc ?comments designs;
      close_out oc

  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
