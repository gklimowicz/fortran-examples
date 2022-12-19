type token =
  | STRING of ( string )
  | INT of ( int )
  | LPAREN
  | RPAREN
  | LBRACKET
  | RBRACKET
  | AND
  | PLUS
  | COLON
  | COMMA
  | NOT
  | HAT
  | ONSHELL
  | OFFSHELL
  | GAUSS
  | END

val main :
  (Lexing.lexbuf  -> token) -> Lexing.lexbuf ->  (string, int list, string) Cascade_syntax.t 
