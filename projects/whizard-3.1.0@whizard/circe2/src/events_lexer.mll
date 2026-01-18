(* events_lexer.mll -- *)
(* Copyright (C) 2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
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

{
  (* Fortran allows ['d'] and ['D'] as exponent starter, but
     O'Caml's [float_of_string] doesn't accept it.  *)

  let normalize_ascii_floats orig =
    let normalized = Bytes.of_string orig in
    for i = 0 to Bytes.length normalized - 1 do
      let c = Bytes.get normalized i in
      if c = 'd' || c = 'D' then
        Bytes.set normalized i 'E'
    done;
    Bytes.to_string normalized
}

let digit = [ '0' - '9' ]
let exp_e = [ 'e' 'E' ]
let exp_d = [ 'd' 'D' ]
let white = [ ' ' '\t' '\r' ]

rule token = parse
    white        { token lexbuf }     (* skip blanks *)
  | ['+''-']? digit+
    ( '.' digit* ( exp_e digit+ )? | exp_e digit+ )
                 { Some (float_of_string (Lexing.lexeme lexbuf)) }
  | ['+''-']? digit+
    ( '.' digit* ( exp_d digit+ )? | exp_d digit+ )
                 { Some (float_of_string (normalize_ascii_floats (Lexing.lexeme lexbuf))) }
  | ['+''-']? digit+
                 { Some (float_of_string (Lexing.lexeme lexbuf)) }
  | _            { failwith (Lexing.lexeme lexbuf) }
  | eof          { None }

{
  (* Not used by circe2, just for illustration and testing. *)

  let float_list_of_string s =
    let lexbuf = Lexing.from_string s in
    let rec collect xs_rev =
      match token lexbuf with
      | None -> List.rev xs_rev
      | Some x -> collect (x :: xs_rev) in
    try
      collect []
    with
    | Failure c ->
       invalid_arg ("invalid token '" ^ c ^ "' in \"" ^ s ^ "\"")

  let float_array_of_string a s =
    let lexbuf = Lexing.from_string s in
    try
      for i = 0 to Array.length a - 1 do
        match token lexbuf with
        | None -> invalid_arg ("not enough floats in \"" ^ s ^ "\"")
        | Some x -> a.(i) <- x
      done
    with
    | Failure c ->
       invalid_arg ("invalid token '" ^ c ^ "' in \"" ^ s ^ "\"")
}
