/* parser.mly --
 */

%{
open Syntax
module Maps = Diffmaps.Default
let parse_error msg =
  raise (Syntax_Error (msg, symbol_start (), symbol_end ()))
%}

%token < int > INT
%token < float > FLOAT
%token < string > STRING
%token SLASH EQUALS STAR PLUS MINUS
%token LBRACKET LPAREN LANGLE COMMA RBRACKET RPAREN RANGLE
%token LBRACE RBRACE
%token Ascii Binary
%token Beta Eta
%token Bins Scale
%token Center
%token Columns
%token Comment
%token Design
%token Electron Positron Photon
%token Events Histogram File
%token Fix
%token Free
%token Id
%token Iterations
%token Lumi Roots
%token Map
%token Min Max
%token Notriangle
%token Pid
%token Pol Unpol
%token Power Resonance
%token Smooth
%token Triangle
%token Width
%token END

%start main
%type < Syntax.file_cmds list > main

%%

main:
    files END                       { $1 }
;

files:
                                    { [] }
  | file files                      { $1 :: $2 }
;

file:
  | LBRACE file_cmds RBRACE         { $2 }
;

file_cmds:
                                    { [] }
  | file_cmd file_cmds              { $1 :: $2 }
;

file_cmd:
    File EQUALS STRING              { Syntax.File $3 }
  | LBRACE design_cmds RBRACE       { Syntax.Designs $2 }
;

design_cmds:
                                    { [] }
  | design_cmd design_cmds          { $1 :: $2 }
;

design_cmd:
    Bins coord EQUALS INT           { Syntax.Design_Bins ($4, $2) }
  | Scale coord EQUALS float        { Syntax.Design_Scale ($4, $2) }
  | Design EQUALS STRING            { Syntax.Design $3 }
  | Roots EQUALS float              { Syntax.Roots $3 }
  | LBRACE channel_cmds RBRACE      { Syntax.Channels $2 }
  | Comment EQUALS STRING           { Syntax.Comment $3 }
;

channel_cmds:
                                    { [] }
  | channel_cmd channel_cmds        { $1 :: $2 }
;

channel_cmd:
    Pid coord EQUALS particle       { Syntax.Pid ($4, $2) }
  | Pol coord EQUALS polarization   { Syntax.Pol ($4, $2) }
  | Fix coord EQUALS side           { Syntax.Fix (true, $2, $4) }
  | Free coord EQUALS side          { Syntax.Fix (false, $2, $4) }
  | Bins coord EQUALS INT           { Syntax.Bins ($4, $2) }
  | Scale coord EQUALS float        { Syntax.Scale ($4, $2) }
  | Min coord EQUALS float          { Syntax.Xmin ($4, $2) }
  | Max coord EQUALS float          { Syntax.Xmax ($4, $2) }
  | Map coord EQUALS map            { Syntax.Diffmap ($4, $2) }
  | Lumi EQUALS float               { Syntax.Lumi $3 }
  | Columns EQUALS INT              { Syntax.Columns $3 }
  | Iterations EQUALS INT           { Syntax.Iterations $3 }
  | Events EQUALS STRING            { Syntax.Events $3 }
  | Histogram EQUALS STRING         { Syntax.Histogram $3 }
  | Binary                          { Syntax.Binary true }
  | Ascii                           { Syntax.Binary false }
  | Smooth EQUALS float area        { Syntax.Smooth ($3, $4) }
  | Triangle                        { Syntax.Triangle true }
  | Notriangle                      { Syntax.Triangle false }
;

particle:
    INT                             { $1 }
  | Electron                        { 11 }
  | Positron                        { -11 }
  | Photon                          { 22 }
;

polarization:
    INT                             { $1 }
  | Unpol                           { 0 }
;

coord:
                                    { Syntax.X12 }
  | SLASH STAR                      { Syntax.X12 }
  | SLASH INT                       {
      match $2 with
      | 1 -> Syntax.X1
      | 2 -> Syntax.X2
      | n ->
          Printf.eprintf "circe2: ignoring dimension %d (not 1, 2, or *)\n" n;
          Syntax.X12 }
;

side:
    Min                             { Syntax.Min }
  | Max                             { Syntax.Max }
  | STAR                            { Syntax.Minmax }
;

map:
   Id LBRACE id RBRACE              { $3 }
 | Power LBRACE power RBRACE        { $3 }
 | Resonance LBRACE resonance RBRACE{ $3 }
;

area:
   interval interval                { Syntax.Rect ($1, $2) }
 | interval point                   { Syntax.Slice1 ($1, $2) }
 | point interval                   { Syntax.Slice2 ($1, $2) }
;

point:
 | LBRACKET float RBRACKET          { Syntax.Delta $2 }
 | LANGLE INT RANGLE                { Syntax.Box $2 }
;

id:
   INT real_interval                {
     let x_min, x_max = $2 in
     ($1, Maps.id x_min x_max) }
;

real_interval:
    left float COMMA float right    { ($2, $4) }
;

left:
    LBRACKET                        { }
  | LPAREN                          { }
;

right:
    RBRACKET                        { }
  | RPAREN                          { }
;

interval:
    lower COMMA upper               { ($1, $3) }
;

lower:
    LBRACKET float                  { Syntax.Closed $2 }
  | LPAREN float                    { Syntax.Open $2 }
  | LANGLE INT                      { Syntax.Bin $2 }
;

upper:
    float RBRACKET                  { Syntax.Closed $1 }
  | float RPAREN                    { Syntax.Open $1 }
  | INT RANGLE                      { Syntax.Bin $1 }
;

power:
   INT real_interval power_params   {
     let x_min, x_max = $2
     and beta, eta = $3 in
     if beta <= -1.0 then begin
       Printf.eprintf "circe2: ignoring invalid beta: %g <= -1\n" beta;
       flush stderr;
       ($1, Maps.id x_min x_max)
     end else
       let alpha = 1.0 /. (1.0 +. beta) in
       ($1, Maps.power ~alpha ~eta x_min x_max) }
;

power_params:
    beta eta                        { ($1, $2) }
  | eta beta                        { ($2, $1) }
;

beta:
    Beta EQUALS float               { $3 }
;

eta:
    Eta EQUALS float                { $3 }
;

resonance:
   INT real_interval resonance_params {
     let x_min, x_max = $2
     and eta, a = $3 in
     ($1, Maps.resonance ~eta ~a x_min x_max) }
;

resonance_params:
    center width                    { ($1, $2) }
  | width center                    { ($2, $1) }
;

center:
    Center EQUALS float             { $3 }
;

width:
    Width EQUALS float              { $3 }
;

float:
    float_or_int                    { $1 }
  | float_or_int PLUS               { $1 +. Syntax.epsilon }
  | float_or_int MINUS              { $1 -. Syntax.epsilon }
;

float_or_int:
    INT                             { float $1 }
  | FLOAT                           { $1 }
;
