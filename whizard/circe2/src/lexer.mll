(* lexer.mll -- 
 *)

{
open Parser
let unquote s =
  String.sub s 1 (String.length s - 2)
}

let digit = ['0'-'9']
let upper = ['A'-'Z']
let lower = ['a'-'z']
let char = upper | lower
let white = [' ' '\t' '\n']

rule token = parse
    white        { token lexbuf }     (* skip blanks *)
  | '#' [^'\n']* '\n'
                 { token lexbuf }     (* skip comments *)
  | ['+''-']? digit+
    ( '.' digit* ( ['e''E'] digit+ )? | ['e''E'] digit+ )
                 { FLOAT (float_of_string (Lexing.lexeme lexbuf)) }
  | ['+''-']? digit+
                 { INT (int_of_string (Lexing.lexeme lexbuf)) }
  | '"' [^'"']* '"'
                 { STRING (unquote (Lexing.lexeme lexbuf)) }
  | '/'          { SLASH }
  | '['          { LBRACKET }
  | '('          { LPAREN }
  | '<'          { LANGLE }
  | ','          { COMMA }
  | ']'          { RBRACKET }
  | ')'          { RPAREN }
  | '>'          { RANGLE }
  | '{'          { LBRACE }
  | '}'          { RBRACE }
  | '='          { EQUALS }
  | '*'          { STAR }
  | '+'          { PLUS }
  | '-'          { MINUS }
  | "ascii"      { Ascii }
  | "beta"       { Beta }
  | "binary"     { Binary }
  | "bins"       { Bins }
  | "center"     { Center }
  | "columns"    { Columns }
  | "comment"    { Comment }
  | "design"     { Design }
  | "electron"   { Electron }
  | "eta"        { Eta }
  | "events"     { Events }
  | "file"       { File }
  | "fix"        { Fix }
  | "free"       { Free }
  | "histogram"  { Histogram }
  | "id"         { Id }
  | "iterations" { Iterations }
  | "lumi"       { Lumi }
  | "map"        { Map }
  | "max"        { Max }
  | "min"        { Min }
  | "notriangle" { Notriangle }
  | "photon"     { Photon }
  | "gamma"      { Photon }
  | "pid"        { Pid }
  | "pol"        { Pol }
  | "positron"   { Positron }
  | "power"      { Power }
  | "resonance"  { Resonance }
  | "roots"      { Roots }
  | "scale"      { Scale }
  | "smooth"     { Smooth }
  | "triangle"   { Triangle }
  | "unpol"      { Unpol }
  | "width"      { Width }
  | eof          { END }
