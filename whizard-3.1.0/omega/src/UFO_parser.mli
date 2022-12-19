type token =
  | INT of ( int )
  | FLOAT of ( float )
  | STRING of ( string )
  | ID of ( string )
  | DOT
  | COMMA
  | COLON
  | EQUAL
  | PLUS
  | MINUS
  | DIV
  | LPAREN
  | RPAREN
  | LBRACE
  | RBRACE
  | LBRACKET
  | RBRACKET
  | END

val file :
  (Lexing.lexbuf  -> token) -> Lexing.lexbuf ->  UFO_syntax.t 
