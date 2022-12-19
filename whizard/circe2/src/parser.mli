type token =
  | INT of ( int )
  | FLOAT of ( float )
  | STRING of ( string )
  | SLASH
  | EQUALS
  | STAR
  | PLUS
  | MINUS
  | LBRACKET
  | LPAREN
  | LANGLE
  | COMMA
  | RBRACKET
  | RPAREN
  | RANGLE
  | LBRACE
  | RBRACE
  | Ascii
  | Binary
  | Beta
  | Eta
  | Bins
  | Scale
  | Center
  | Columns
  | Comment
  | Design
  | Electron
  | Positron
  | Photon
  | Events
  | Histogram
  | File
  | Fix
  | Free
  | Id
  | Iterations
  | Lumi
  | Roots
  | Map
  | Min
  | Max
  | Notriangle
  | Pid
  | Pol
  | Unpol
  | Power
  | Resonance
  | Smooth
  | Triangle
  | Width
  | END

val main :
  (Lexing.lexbuf  -> token) -> Lexing.lexbuf ->  Syntax.file_cmds list 
