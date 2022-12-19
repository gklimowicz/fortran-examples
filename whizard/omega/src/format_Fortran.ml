(* format_Fortran.ml -- Fortran90+ continuation lines etc.

   Copyright (C) 2019-2022 by

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

let default_width = 80

let max_clines = ref (-1) (* 255 *)
exception Continuation_Lines of int

(* Fortran style line continuation: *)

type formatter =
  { formatter : Format.formatter;
    mutable current_cline : int;
    mutable width : int }

let formatter_of_formatter ?(width=default_width) ff =
  { formatter = ff;
    current_cline = 1;
    width = width }

(* Default function to output new lines. *)
let pp_output_function ff =
  fst (Format.pp_get_formatter_output_functions ff.formatter ())

(* Default function to output spaces (copied from \texttt{format.ml}). *)
let blank_line = String.make 80 ' '
let rec pp_display_blanks ff n =
  if n > 0 then
    if n <= 80 then
      pp_output_function ff blank_line 0 n
    else begin
        pp_output_function ff blank_line 0 80;
        pp_display_blanks ff (n - 80)
      end

let pp_display_newline ff =
  pp_output_function ff "\n" 0  1

(* [ff.current_cline]
   \begin{itemize}
     \item $\le0$: not continuing: print a straight newline,
     \item $>0$: continuing: append [" &"] until we run up to [!max_clines].
       NB: [!max_clines < 0] means \emph{unlimited} continuation lines.
   \end{itemize} *)

let pp_switch_line_continuation ff = function
  | false -> ff.current_cline <- 0
  | true -> ff.current_cline <- 1

let pp_fortran_newline ff () =
  if ff.current_cline > 0 then
    begin
      if !max_clines >= 0 && ff.current_cline > !max_clines then
        raise (Continuation_Lines ff.current_cline)
      else
        begin
          pp_output_function ff " &" 0 2;
          ff.current_cline <- succ ff.current_cline
        end
    end;
  pp_display_newline ff

let pp_newline ff () =
  pp_switch_line_continuation ff false;
  Format.pp_print_newline ff.formatter ();
  pp_switch_line_continuation ff true

(* Make a formatter with default functions to output spaces and new lines. *)

(*i
let unsafe_output oc s i j =
  try
    output oc s i j
  with
  | _ -> Printf.eprintf "unsafe_output: '%s'\n" s
i*)

let pp_setup ff =
  let formatter_out_functions =
    Format.pp_get_formatter_out_functions ff.formatter () in
  Format.pp_set_formatter_out_functions
    ff.formatter
    { formatter_out_functions with
      Format.out_newline = pp_fortran_newline ff;
      Format.out_spaces = pp_display_blanks ff };
  Format.pp_set_margin ff.formatter (ff.width - 2)

let std_formatter =
  let ff = formatter_of_formatter Format.std_formatter in
  pp_setup ff;
  ff

let formatter_of_out_channel ?(width=default_width) oc =
  let ff = formatter_of_formatter ~width (Format.formatter_of_out_channel oc) in
  pp_setup ff;
  ff

let formatter_of_buffer ?(width=default_width) b =
  let ff =
    { formatter = Format.formatter_of_buffer b;
      current_cline = 1;
      width = width } in
  pp_setup ff;
  ff

let pp_set_formatter_out_channel ff ?(width=default_width) oc =
  Format.pp_set_formatter_out_channel ff.formatter oc;
  ff.width <- width;
  pp_setup ff

let set_formatter_out_channel ?(width=default_width) oc =
  Format.pp_set_formatter_out_channel std_formatter.formatter oc;
  std_formatter.width <- width;
  pp_setup std_formatter

let fprintf ff fmt = Format.fprintf ff.formatter fmt
let pp_flush ff = Format.pp_print_flush ff.formatter

let printf fmt = fprintf std_formatter fmt
let newline = pp_newline std_formatter
let flush = pp_flush std_formatter
let switch_line_continuation = pp_switch_line_continuation std_formatter

module Test =
  struct

    open OUnit

    let input_line_opt ic =
      try
        Some (input_line ic)
      with
      | End_of_file -> None
 
    let read_lines ic =
      let rec read_lines' acc =
        match input_line_opt ic with
        | Some line -> read_lines' (line :: acc)
        | None -> List.rev acc
      in
      read_lines' []
 
    let lines_of_file filename =
      let ic = open_in filename in
      let lines = read_lines ic in
      close_in ic;
      lines

    let equal_or_dump_lines lhs rhs =
      if lhs = rhs then
        true
      else
        begin
          Printf.printf "Unexpected output:\n";
          List.iter (Printf.printf "< %s\n") lhs;
          List.iter (Printf.printf "> %s\n") rhs;
          false
        end

    let format_and_compare f expected () =
      bracket_tmpfile
        ~prefix:"omega-" ~suffix:".f90"
        (fun (name, oc) ->
          (* There can be something left in the queue from [OUnit]! *)
          Format.print_flush ();
          f oc;
          close_out oc;
          (* [OUnit] uses [Format.printf]! *)
          Format.set_formatter_out_channel stdout;
          assert_bool "" (equal_or_dump_lines expected (lines_of_file name)))
        ()

    let suite =
      "Format_Fortran" >:::
        [ "formatter_of_out_channel" >::
            format_and_compare
              (fun oc ->
                let ff = formatter_of_out_channel ~width:20 oc in
                let nl = pp_newline ff in
                List.iter
                  (fprintf ff)
                  ["@[<2>lhs = rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"];
                nl ())
              [ "lhs = rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs" ];

          "formatter_of_buffer" >::
            format_and_compare
              (fun oc ->
                let buffer = Buffer.create 1024 in
                let ff = formatter_of_buffer ~width:20 buffer in
                let nl = pp_newline ff in
                List.iter
                  (fprintf ff)
                  ["  @[<2>lhs = rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"];
                nl ();
                pp_flush ff ();
                let ff' = formatter_of_out_channel ~width:20 oc in
                fprintf ff' "do mu = 0, 3"; pp_newline ff' ();
                fprintf ff' "%s" (Buffer.contents buffer);
                fprintf ff' "end do";
                pp_newline ff' ())
              [ "do mu = 0, 3";
                "  lhs = rhs + rhs &";
                "    + rhs + rhs &";
                "    + rhs + rhs &";
                "    + rhs + rhs &";
                "    + rhs + rhs &";
                "    + rhs";
                "end do" ];

          "formatter_of_out_channel+indentation" >::
            format_and_compare
              (fun oc ->
                let ff = formatter_of_out_channel ~width:20 oc in
                let nl = pp_newline ff in
                List.iter
                  (fprintf ff)
                  ["  @[<4>lhs = rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"];
                nl ())
              [ "  lhs = rhs + rhs &";
                "      + rhs + rhs &";
                "      + rhs + rhs &";
                "      + rhs + rhs &";
                "      + rhs + rhs &";
                "      + rhs" ];

          "set_formatter_out_channel" >::
            format_and_compare
              (fun oc ->
                let nl = newline in
                set_formatter_out_channel ~width:20 oc;
                List.iter
                  printf
                  ["@[<2>lhs = rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs";
                   "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"; "@ + rhs"];
                nl ())
              [ "lhs = rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs + rhs &";
                "  + rhs" ]; ]

  end
