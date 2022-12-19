type token =
  | INT of ( int )
  | FLOAT of ( float )
  | ID of ( string )
  | QUOTED of ( string )
  | PLUS
  | MINUS
  | TIMES
  | POWER
  | DIV
  | LPAREN
  | RPAREN
  | COMMA
  | DOT
  | END

val input :
  (Lexing.lexbuf  -> token) -> Lexing.lexbuf ->  UFOx_syntax.expr 
