(* progress.ml --

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

type channel =
  | Channel of out_channel
  | File of string
  | Open_File of string * out_channel

type state =
    { channel : channel;
      mutable steps : int;
      mutable digits : int;
      mutable step : int;
      created : float;
      mutable last_reset : float;
      mutable last_begin : float; }

type t = state option

let digits n =
  if n > 0 then
    succ (truncate (log10 (float n)))
  else
    invalid_arg "Progress.digits: non-positive argument"

let mod_float2 a b =
  let modulus = mod_float a b in
  ((a -. modulus) /. b, modulus)

let time_to_string seconds =
  let minutes, seconds = mod_float2 seconds 60. in
  if minutes > 0.0 then
    let hours, minutes = mod_float2 minutes 60. in
    if hours > 0.0 then
      let days, hours = mod_float2 hours 24. in
      if days > 0.0 then
        Printf.sprintf "%.0f:%02.0f days" days hours
      else
        Printf.sprintf "%.0f:%02.0f hrs" hours minutes
    else
      Printf.sprintf "%.0f:%02.0f mins" minutes seconds
  else
    Printf.sprintf "%.2f secs" seconds

let create channel steps = 
  let now = Sys.time () in
  Some { channel = channel;
         steps = steps;
         digits = digits steps;
         step = 0;
         created = now;
         last_reset = now;
         last_begin = now }

let dummy =
  None 

let channel oc =
  create (Channel oc)

let file name = 
  let oc = open_out name in
  close_out oc;
  create (File name)

let open_file name = 
  let oc = open_out name in
  create (Open_File (name, oc))

let close_channel state =
  match state.channel with
  | Channel oc ->
      flush oc
  | File _ -> ()
  | Open_File (_, oc) ->
      flush oc;
      close_out oc

let use_channel state f =
  match state.channel with
  | Channel oc | Open_File (_, oc) ->
      f oc;
      flush oc
  | File name ->
      let oc = open_out_gen [Open_append; Open_creat] 0o644 name in
      f oc;
      flush oc;
      close_out oc
  
let reset state steps msg =
  match state with
  | None -> ()
  | Some state ->
      let now = Sys.time () in
      state.steps <- steps;
      state.digits <- digits steps;
      state.step <- 0;
      state.last_reset <- now;
      state.last_begin <- now

let begin_step state msg =
  match state with
  | None -> ()
  | Some state ->
      let now = Sys.time () in
      state.step <- succ state.step;
      state.last_begin <- now;
      use_channel state (fun oc ->
        Printf.fprintf oc "[%0*d/%0*d] %s ..." state.digits state.step state.digits state.steps msg)

let end_step state msg =
  match state with
  | None -> ()
  | Some state ->
      let now = Sys.time () in
      let last = now -. state.last_begin in
      let elapsed = now -. state.last_reset in
      let estimated = float state.steps *. elapsed /. float state.step in
      let remaining = estimated -. elapsed in
      use_channel state (fun oc ->
        Printf.fprintf oc " %s. [time: %s, total: %s, remaining: %s]\n" msg
          (time_to_string last) (time_to_string estimated) (time_to_string remaining))

let summary state msg =
  match state with
  | None -> ()
  | Some state ->
      let now = Sys.time () in
      use_channel state (fun oc ->
        Printf.fprintf oc "%s. [total time: %s]\n" msg
          (time_to_string (now -. state.created)));
      close_channel state
        
(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)





