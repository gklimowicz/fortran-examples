type token =
  | DIGIT of ( int )
  | CHAR of ( string )
  | PREFIX of ( string )
  | TOKEN of ( string )
  | SUPER
  | SUB
  | PRIME
  | LBRACE
  | RBRACE
  | LBRACKET
  | RBRACKET
  | LPAREN
  | RPAREN
  | COMMA
  | PLUS
  | MINUS
  | TIMES
  | DIV
  | EQUAL
  | INCLUDE of ( string )
  | END
  | NEUTRAL
  | CHARGED
  | ANTI
  | ALIAS
  | TEX
  | FORTRAN
  | SPIN
  | COLOR
  | CHARGE
  | MASS
  | WIDTH
  | PARAMETER
  | DERIVED
  | TENSOR
  | INDEX
  | FLAVOR
  | LORENTZ
  | VERTEX

val file :
  (Lexing.lexbuf  -> token) -> Lexing.lexbuf ->  Vertex_syntax.File_Tree.t 
