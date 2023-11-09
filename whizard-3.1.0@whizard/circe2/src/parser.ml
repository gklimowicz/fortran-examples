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

open Parsing;;
let _ = parse_error;;
# 5 "../../../circe2/src/parser.mly"
open Syntax
module Maps = Diffmaps.Default
let parse_error msg =
  raise (Syntax_Error (msg, symbol_start (), symbol_end ()))
# 63 "parser.ml"
let yytransl_const = [|
  260 (* SLASH *);
  261 (* EQUALS *);
  262 (* STAR *);
  263 (* PLUS *);
  264 (* MINUS *);
  265 (* LBRACKET *);
  266 (* LPAREN *);
  267 (* LANGLE *);
  268 (* COMMA *);
  269 (* RBRACKET *);
  270 (* RPAREN *);
  271 (* RANGLE *);
  272 (* LBRACE *);
  273 (* RBRACE *);
  274 (* Ascii *);
  275 (* Binary *);
  276 (* Beta *);
  277 (* Eta *);
  278 (* Bins *);
  279 (* Scale *);
  280 (* Center *);
  281 (* Columns *);
  282 (* Comment *);
  283 (* Design *);
  284 (* Electron *);
  285 (* Positron *);
  286 (* Photon *);
  287 (* Events *);
  288 (* Histogram *);
  289 (* File *);
  290 (* Fix *);
  291 (* Free *);
  292 (* Id *);
  293 (* Iterations *);
  294 (* Lumi *);
  295 (* Roots *);
  296 (* Map *);
  297 (* Min *);
  298 (* Max *);
  299 (* Notriangle *);
  300 (* Pid *);
  301 (* Pol *);
  302 (* Unpol *);
  303 (* Power *);
  304 (* Resonance *);
  305 (* Smooth *);
  306 (* Triangle *);
  307 (* Width *);
  308 (* END *);
    0|]

let yytransl_block = [|
  257 (* INT *);
  258 (* FLOAT *);
  259 (* STRING *);
    0|]

let yylhs = "\255\255\
\001\000\002\000\002\000\003\000\004\000\004\000\005\000\005\000\
\006\000\006\000\007\000\007\000\007\000\007\000\007\000\007\000\
\010\000\010\000\011\000\011\000\011\000\011\000\011\000\011\000\
\011\000\011\000\011\000\011\000\011\000\011\000\011\000\011\000\
\011\000\011\000\011\000\011\000\011\000\012\000\012\000\012\000\
\012\000\013\000\013\000\008\000\008\000\008\000\014\000\014\000\
\014\000\015\000\015\000\015\000\016\000\016\000\016\000\021\000\
\021\000\017\000\022\000\023\000\023\000\024\000\024\000\020\000\
\025\000\025\000\025\000\026\000\026\000\026\000\018\000\027\000\
\027\000\028\000\029\000\019\000\030\000\030\000\031\000\032\000\
\009\000\009\000\009\000\033\000\033\000\000\000"

let yylen = "\002\000\
\002\000\000\000\002\000\003\000\000\000\002\000\003\000\003\000\
\000\000\002\000\004\000\004\000\003\000\003\000\003\000\003\000\
\000\000\002\000\004\000\004\000\004\000\004\000\004\000\004\000\
\004\000\004\000\004\000\003\000\003\000\003\000\003\000\003\000\
\001\000\001\000\004\000\001\000\001\000\001\000\001\000\001\000\
\001\000\001\000\001\000\000\000\002\000\002\000\001\000\001\000\
\001\000\004\000\004\000\004\000\002\000\002\000\002\000\003\000\
\003\000\002\000\005\000\001\000\001\000\001\000\001\000\003\000\
\002\000\002\000\002\000\002\000\002\000\002\000\003\000\002\000\
\002\000\003\000\003\000\003\000\002\000\002\000\003\000\003\000\
\001\000\002\000\002\000\001\000\001\000\002\000"

let yydefred = "\000\000\
\000\000\000\000\000\000\086\000\000\000\000\000\000\000\000\000\
\000\000\000\000\001\000\003\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\004\000\006\000\034\000\
\033\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\037\000\000\000\000\000\
\000\000\036\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\008\000\010\000\007\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\015\000\018\000\046\000\045\000\
\000\000\000\000\016\000\013\000\084\000\085\000\014\000\000\000\
\000\000\000\000\029\000\031\000\032\000\000\000\000\000\030\000\
\028\000\000\000\000\000\000\000\000\000\000\000\000\000\011\000\
\012\000\082\000\083\000\023\000\024\000\049\000\047\000\048\000\
\021\000\022\000\000\000\000\000\000\000\027\000\025\000\026\000\
\038\000\039\000\040\000\041\000\019\000\042\000\043\000\020\000\
\000\000\000\000\000\000\035\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\066\000\000\000\053\000\054\000\000\000\
\000\000\055\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\056\000\057\000\065\000\067\000\000\000\000\000\064\000\
\060\000\061\000\058\000\000\000\050\000\000\000\051\000\000\000\
\052\000\070\000\068\000\069\000\000\000\000\000\000\000\071\000\
\000\000\000\000\000\000\000\000\076\000\000\000\000\000\000\000\
\000\000\000\000\072\000\073\000\000\000\000\000\077\000\078\000\
\000\000\074\000\075\000\079\000\080\000\062\000\063\000\059\000"

let yydgoto = "\002\000\
\004\000\005\000\006\000\009\000\010\000\019\000\020\000\046\000\
\079\000\043\000\044\000\117\000\120\000\105\000\110\000\124\000\
\141\000\143\000\145\000\125\000\126\000\155\000\156\000\192\000\
\127\000\152\000\168\000\169\000\170\000\173\000\174\000\175\000\
\080\000"

let yysindex = "\007\000\
\255\254\000\000\245\254\000\000\229\254\255\254\243\254\022\255\
\030\255\245\254\000\000\000\000\047\255\050\255\050\255\046\255\
\068\255\075\255\054\255\243\254\083\255\000\000\000\000\000\000\
\000\000\050\255\050\255\078\255\088\255\090\255\050\255\050\255\
\101\255\102\255\050\255\050\255\050\255\000\000\050\255\050\255\
\103\255\000\000\092\255\047\255\010\255\105\255\106\255\109\255\
\110\255\056\255\000\000\000\000\000\000\113\255\122\255\116\255\
\137\255\138\255\140\255\143\255\141\255\056\255\153\255\154\255\
\155\255\159\255\160\255\056\255\000\000\000\000\000\000\000\000\
\165\255\056\255\000\000\000\000\000\000\000\000\000\000\060\255\
\166\255\056\255\000\000\000\000\000\000\252\254\252\254\000\000\
\000\000\244\254\056\255\056\255\003\255\000\255\008\255\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\152\255\156\255\157\255\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\056\255\056\255\168\255\000\000\008\255\032\255\158\255\170\255\
\173\255\174\255\163\255\000\000\162\255\000\000\000\000\056\255\
\177\255\000\000\062\255\066\255\164\255\066\255\167\255\066\255\
\169\255\000\000\000\000\000\000\000\000\172\255\085\255\000\000\
\000\000\000\000\000\000\056\255\000\000\080\255\000\000\239\254\
\000\000\000\000\000\000\000\000\171\255\175\255\180\255\000\000\
\161\255\176\255\183\255\184\255\000\000\128\255\178\255\056\255\
\056\255\056\255\000\000\000\000\056\255\056\255\000\000\000\000\
\089\255\000\000\000\000\000\000\000\000\000\000\000\000\000\000"

let yyrindex = "\000\000\
\139\255\000\000\181\255\000\000\000\000\139\255\182\255\000\000\
\000\000\181\255\000\000\000\000\186\255\185\255\185\255\000\000\
\000\000\000\000\000\000\182\255\000\000\000\000\000\000\000\000\
\000\000\185\255\185\255\000\000\000\000\000\000\185\255\185\255\
\000\000\000\000\185\255\185\255\185\255\000\000\185\255\185\255\
\000\000\000\000\000\000\186\255\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\112\255\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\188\255\000\000\189\255\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\048\255\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000"

let yygindex = "\000\000\
\000\000\186\000\000\000\183\000\000\000\174\000\000\000\013\000\
\194\255\151\000\000\000\000\000\000\000\110\000\000\000\000\000\
\000\000\000\000\000\000\235\255\079\000\135\255\000\000\000\000\
\000\000\000\000\000\000\035\000\037\000\000\000\032\000\034\000\
\000\000"

let yytablesize = 208
let yytable = "\089\000\
\118\000\102\000\013\000\113\000\007\000\095\000\171\000\001\000\
\014\000\015\000\071\000\097\000\016\000\017\000\003\000\072\000\
\121\000\122\000\123\000\101\000\158\000\008\000\160\000\107\000\
\011\000\018\000\021\000\047\000\111\000\112\000\114\000\115\000\
\116\000\172\000\108\000\109\000\103\000\104\000\054\000\055\000\
\136\000\122\000\137\000\059\000\060\000\119\000\022\000\063\000\
\064\000\065\000\048\000\066\000\067\000\045\000\084\000\084\000\
\077\000\078\000\131\000\132\000\084\000\084\000\150\000\078\000\
\024\000\025\000\098\000\099\000\026\000\027\000\051\000\028\000\
\049\000\148\000\153\000\154\000\151\000\029\000\030\000\050\000\
\031\000\032\000\056\000\033\000\034\000\053\000\035\000\036\000\
\037\000\038\000\039\000\040\000\057\000\165\000\058\000\041\000\
\042\000\163\000\164\000\166\000\167\000\190\000\191\000\134\000\
\138\000\061\000\062\000\068\000\069\000\073\000\074\000\075\000\
\076\000\185\000\186\000\187\000\083\000\081\000\188\000\189\000\
\081\000\081\000\081\000\081\000\081\000\081\000\082\000\081\000\
\081\000\081\000\081\000\081\000\081\000\081\000\081\000\081\000\
\081\000\081\000\081\000\084\000\085\000\088\000\081\000\081\000\
\086\000\081\000\081\000\087\000\081\000\081\000\081\000\081\000\
\081\000\081\000\081\000\081\000\081\000\090\000\091\000\092\000\
\081\000\081\000\081\000\093\000\094\000\096\000\100\000\128\000\
\133\000\139\000\140\000\129\000\130\000\142\000\144\000\146\000\
\147\000\149\000\172\000\177\000\157\000\167\000\176\000\159\000\
\178\000\161\000\162\000\181\000\182\000\044\000\002\000\012\000\
\023\000\052\000\070\000\166\000\106\000\005\000\009\000\065\000\
\067\000\171\000\017\000\135\000\180\000\179\000\184\000\183\000"

let yycheck = "\062\000\
\001\001\006\001\016\001\001\001\016\001\068\000\024\001\001\000\
\022\001\023\001\001\001\074\000\026\001\027\001\016\001\006\001\
\009\001\010\001\011\001\082\000\142\000\033\001\144\000\036\001\
\052\001\039\001\005\001\015\000\091\000\092\000\028\001\029\001\
\030\001\051\001\047\001\048\001\041\001\042\001\026\000\027\000\
\009\001\010\001\011\001\031\000\032\000\046\001\017\001\035\000\
\036\000\037\000\005\001\039\000\040\000\004\001\007\001\008\001\
\001\001\002\001\121\000\122\000\013\001\014\001\001\001\002\001\
\018\001\019\001\007\001\008\001\022\001\023\001\017\001\025\001\
\005\001\136\000\009\001\010\001\139\000\031\001\032\001\005\001\
\034\001\035\001\005\001\037\001\038\001\003\001\040\001\041\001\
\042\001\043\001\044\001\045\001\005\001\156\000\005\001\049\001\
\050\001\013\001\014\001\020\001\021\001\013\001\014\001\125\000\
\126\000\005\001\005\001\005\001\017\001\005\001\005\001\003\001\
\003\001\176\000\177\000\178\000\001\001\005\001\181\000\182\000\
\009\001\010\001\011\001\012\001\013\001\014\001\005\001\016\001\
\017\001\018\001\019\001\020\001\021\001\022\001\023\001\024\001\
\025\001\026\001\027\001\003\001\003\001\001\001\031\001\032\001\
\005\001\034\001\035\001\005\001\037\001\038\001\039\001\040\001\
\041\001\042\001\043\001\044\001\045\001\005\001\005\001\005\001\
\049\001\050\001\051\001\005\001\005\001\001\001\001\001\016\001\
\001\001\012\001\001\001\016\001\016\001\001\001\001\001\013\001\
\015\001\001\001\051\001\005\001\017\001\021\001\012\001\017\001\
\005\001\017\001\015\001\005\001\005\001\005\001\052\001\006\000\
\010\000\020\000\044\000\020\001\087\000\017\001\017\001\012\001\
\012\001\024\001\017\001\125\000\170\000\169\000\175\000\174\000"

let yynames_const = "\
  SLASH\000\
  EQUALS\000\
  STAR\000\
  PLUS\000\
  MINUS\000\
  LBRACKET\000\
  LPAREN\000\
  LANGLE\000\
  COMMA\000\
  RBRACKET\000\
  RPAREN\000\
  RANGLE\000\
  LBRACE\000\
  RBRACE\000\
  Ascii\000\
  Binary\000\
  Beta\000\
  Eta\000\
  Bins\000\
  Scale\000\
  Center\000\
  Columns\000\
  Comment\000\
  Design\000\
  Electron\000\
  Positron\000\
  Photon\000\
  Events\000\
  Histogram\000\
  File\000\
  Fix\000\
  Free\000\
  Id\000\
  Iterations\000\
  Lumi\000\
  Roots\000\
  Map\000\
  Min\000\
  Max\000\
  Notriangle\000\
  Pid\000\
  Pol\000\
  Unpol\000\
  Power\000\
  Resonance\000\
  Smooth\000\
  Triangle\000\
  Width\000\
  END\000\
  "

let yynames_block = "\
  INT\000\
  FLOAT\000\
  STRING\000\
  "

let yyact = [|
  (fun _ -> failwith "parser")
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'files) in
    Obj.repr(
# 48 "../../../circe2/src/parser.mly"
                                    ( _1 )
# 362 "parser.ml"
               :  Syntax.file_cmds list ))
; (fun __caml_parser_env ->
    Obj.repr(
# 52 "../../../circe2/src/parser.mly"
                                    ( [] )
# 368 "parser.ml"
               : 'files))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'file) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'files) in
    Obj.repr(
# 53 "../../../circe2/src/parser.mly"
                                    ( _1 :: _2 )
# 376 "parser.ml"
               : 'files))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'file_cmds) in
    Obj.repr(
# 57 "../../../circe2/src/parser.mly"
                                    ( _2 )
# 383 "parser.ml"
               : 'file))
; (fun __caml_parser_env ->
    Obj.repr(
# 61 "../../../circe2/src/parser.mly"
                                    ( [] )
# 389 "parser.ml"
               : 'file_cmds))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'file_cmd) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'file_cmds) in
    Obj.repr(
# 62 "../../../circe2/src/parser.mly"
                                    ( _1 :: _2 )
# 397 "parser.ml"
               : 'file_cmds))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 66 "../../../circe2/src/parser.mly"
                                    ( Syntax.File _3 )
# 404 "parser.ml"
               : 'file_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'design_cmds) in
    Obj.repr(
# 67 "../../../circe2/src/parser.mly"
                                    ( Syntax.Designs _2 )
# 411 "parser.ml"
               : 'file_cmd))
; (fun __caml_parser_env ->
    Obj.repr(
# 71 "../../../circe2/src/parser.mly"
                                    ( [] )
# 417 "parser.ml"
               : 'design_cmds))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'design_cmd) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'design_cmds) in
    Obj.repr(
# 72 "../../../circe2/src/parser.mly"
                                    ( _1 :: _2 )
# 425 "parser.ml"
               : 'design_cmds))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 76 "../../../circe2/src/parser.mly"
                                    ( Syntax.Design_Bins (_4, _2) )
# 433 "parser.ml"
               : 'design_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 77 "../../../circe2/src/parser.mly"
                                    ( Syntax.Design_Scale (_4, _2) )
# 441 "parser.ml"
               : 'design_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 78 "../../../circe2/src/parser.mly"
                                    ( Syntax.Design _3 )
# 448 "parser.ml"
               : 'design_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 79 "../../../circe2/src/parser.mly"
                                    ( Syntax.Roots _3 )
# 455 "parser.ml"
               : 'design_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'channel_cmds) in
    Obj.repr(
# 80 "../../../circe2/src/parser.mly"
                                    ( Syntax.Channels _2 )
# 462 "parser.ml"
               : 'design_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 81 "../../../circe2/src/parser.mly"
                                    ( Syntax.Comment _3 )
# 469 "parser.ml"
               : 'design_cmd))
; (fun __caml_parser_env ->
    Obj.repr(
# 85 "../../../circe2/src/parser.mly"
                                    ( [] )
# 475 "parser.ml"
               : 'channel_cmds))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'channel_cmd) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'channel_cmds) in
    Obj.repr(
# 86 "../../../circe2/src/parser.mly"
                                    ( _1 :: _2 )
# 483 "parser.ml"
               : 'channel_cmds))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'particle) in
    Obj.repr(
# 90 "../../../circe2/src/parser.mly"
                                    ( Syntax.Pid (_4, _2) )
# 491 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'polarization) in
    Obj.repr(
# 91 "../../../circe2/src/parser.mly"
                                    ( Syntax.Pol (_4, _2) )
# 499 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'side) in
    Obj.repr(
# 92 "../../../circe2/src/parser.mly"
                                    ( Syntax.Fix (true, _2, _4) )
# 507 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'side) in
    Obj.repr(
# 93 "../../../circe2/src/parser.mly"
                                    ( Syntax.Fix (false, _2, _4) )
# 515 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 94 "../../../circe2/src/parser.mly"
                                    ( Syntax.Bins (_4, _2) )
# 523 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 95 "../../../circe2/src/parser.mly"
                                    ( Syntax.Scale (_4, _2) )
# 531 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 96 "../../../circe2/src/parser.mly"
                                    ( Syntax.Xmin (_4, _2) )
# 539 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 97 "../../../circe2/src/parser.mly"
                                    ( Syntax.Xmax (_4, _2) )
# 547 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'coord) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'map) in
    Obj.repr(
# 98 "../../../circe2/src/parser.mly"
                                    ( Syntax.Diffmap (_4, _2) )
# 555 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 99 "../../../circe2/src/parser.mly"
                                    ( Syntax.Lumi _3 )
# 562 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 100 "../../../circe2/src/parser.mly"
                                    ( Syntax.Columns _3 )
# 569 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 101 "../../../circe2/src/parser.mly"
                                    ( Syntax.Iterations _3 )
# 576 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 102 "../../../circe2/src/parser.mly"
                                    ( Syntax.Events _3 )
# 583 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 103 "../../../circe2/src/parser.mly"
                                    ( Syntax.Histogram _3 )
# 590 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    Obj.repr(
# 104 "../../../circe2/src/parser.mly"
                                    ( Syntax.Binary true )
# 596 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    Obj.repr(
# 105 "../../../circe2/src/parser.mly"
                                    ( Syntax.Binary false )
# 602 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'float) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'area) in
    Obj.repr(
# 106 "../../../circe2/src/parser.mly"
                                    ( Syntax.Smooth (_3, _4) )
# 610 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    Obj.repr(
# 107 "../../../circe2/src/parser.mly"
                                    ( Syntax.Triangle true )
# 616 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    Obj.repr(
# 108 "../../../circe2/src/parser.mly"
                                    ( Syntax.Triangle false )
# 622 "parser.ml"
               : 'channel_cmd))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 112 "../../../circe2/src/parser.mly"
                                    ( _1 )
# 629 "parser.ml"
               : 'particle))
; (fun __caml_parser_env ->
    Obj.repr(
# 113 "../../../circe2/src/parser.mly"
                                    ( 11 )
# 635 "parser.ml"
               : 'particle))
; (fun __caml_parser_env ->
    Obj.repr(
# 114 "../../../circe2/src/parser.mly"
                                    ( -11 )
# 641 "parser.ml"
               : 'particle))
; (fun __caml_parser_env ->
    Obj.repr(
# 115 "../../../circe2/src/parser.mly"
                                    ( 22 )
# 647 "parser.ml"
               : 'particle))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 119 "../../../circe2/src/parser.mly"
                                    ( _1 )
# 654 "parser.ml"
               : 'polarization))
; (fun __caml_parser_env ->
    Obj.repr(
# 120 "../../../circe2/src/parser.mly"
                                    ( 0 )
# 660 "parser.ml"
               : 'polarization))
; (fun __caml_parser_env ->
    Obj.repr(
# 124 "../../../circe2/src/parser.mly"
                                    ( Syntax.X12 )
# 666 "parser.ml"
               : 'coord))
; (fun __caml_parser_env ->
    Obj.repr(
# 125 "../../../circe2/src/parser.mly"
                                    ( Syntax.X12 )
# 672 "parser.ml"
               : 'coord))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 126 "../../../circe2/src/parser.mly"
                                    (
      match _2 with
      | 1 -> Syntax.X1
      | 2 -> Syntax.X2
      | n ->
          Printf.eprintf "circe2: ignoring dimension %d (not 1, 2, or *)\n" n;
          Syntax.X12 )
# 685 "parser.ml"
               : 'coord))
; (fun __caml_parser_env ->
    Obj.repr(
# 136 "../../../circe2/src/parser.mly"
                                    ( Syntax.Min )
# 691 "parser.ml"
               : 'side))
; (fun __caml_parser_env ->
    Obj.repr(
# 137 "../../../circe2/src/parser.mly"
                                    ( Syntax.Max )
# 697 "parser.ml"
               : 'side))
; (fun __caml_parser_env ->
    Obj.repr(
# 138 "../../../circe2/src/parser.mly"
                                    ( Syntax.Minmax )
# 703 "parser.ml"
               : 'side))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'id) in
    Obj.repr(
# 142 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 710 "parser.ml"
               : 'map))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'power) in
    Obj.repr(
# 143 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 717 "parser.ml"
               : 'map))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'resonance) in
    Obj.repr(
# 144 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 724 "parser.ml"
               : 'map))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'interval) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'interval) in
    Obj.repr(
# 148 "../../../circe2/src/parser.mly"
                                    ( Syntax.Rect (_1, _2) )
# 732 "parser.ml"
               : 'area))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'interval) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'point) in
    Obj.repr(
# 149 "../../../circe2/src/parser.mly"
                                    ( Syntax.Slice1 (_1, _2) )
# 740 "parser.ml"
               : 'area))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'point) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'interval) in
    Obj.repr(
# 150 "../../../circe2/src/parser.mly"
                                    ( Syntax.Slice2 (_1, _2) )
# 748 "parser.ml"
               : 'area))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'float) in
    Obj.repr(
# 154 "../../../circe2/src/parser.mly"
                                    ( Syntax.Delta _2 )
# 755 "parser.ml"
               : 'point))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 :  int ) in
    Obj.repr(
# 155 "../../../circe2/src/parser.mly"
                                    ( Syntax.Box _2 )
# 762 "parser.ml"
               : 'point))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 :  int ) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'real_interval) in
    Obj.repr(
# 159 "../../../circe2/src/parser.mly"
                                    (
     let x_min, x_max = _2 in
     (_1, Maps.id x_min x_max) )
# 772 "parser.ml"
               : 'id))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 4 : 'left) in
    let _2 = (Parsing.peek_val __caml_parser_env 3 : 'float) in
    let _4 = (Parsing.peek_val __caml_parser_env 1 : 'float) in
    let _5 = (Parsing.peek_val __caml_parser_env 0 : 'right) in
    Obj.repr(
# 165 "../../../circe2/src/parser.mly"
                                    ( (_2, _4) )
# 782 "parser.ml"
               : 'real_interval))
; (fun __caml_parser_env ->
    Obj.repr(
# 169 "../../../circe2/src/parser.mly"
                                    ( )
# 788 "parser.ml"
               : 'left))
; (fun __caml_parser_env ->
    Obj.repr(
# 170 "../../../circe2/src/parser.mly"
                                    ( )
# 794 "parser.ml"
               : 'left))
; (fun __caml_parser_env ->
    Obj.repr(
# 174 "../../../circe2/src/parser.mly"
                                    ( )
# 800 "parser.ml"
               : 'right))
; (fun __caml_parser_env ->
    Obj.repr(
# 175 "../../../circe2/src/parser.mly"
                                    ( )
# 806 "parser.ml"
               : 'right))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'lower) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'upper) in
    Obj.repr(
# 179 "../../../circe2/src/parser.mly"
                                    ( (_1, _3) )
# 814 "parser.ml"
               : 'interval))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 183 "../../../circe2/src/parser.mly"
                                    ( Syntax.Closed _2 )
# 821 "parser.ml"
               : 'lower))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 184 "../../../circe2/src/parser.mly"
                                    ( Syntax.Open _2 )
# 828 "parser.ml"
               : 'lower))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 185 "../../../circe2/src/parser.mly"
                                    ( Syntax.Bin _2 )
# 835 "parser.ml"
               : 'lower))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'float) in
    Obj.repr(
# 189 "../../../circe2/src/parser.mly"
                                    ( Syntax.Closed _1 )
# 842 "parser.ml"
               : 'upper))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'float) in
    Obj.repr(
# 190 "../../../circe2/src/parser.mly"
                                    ( Syntax.Open _1 )
# 849 "parser.ml"
               : 'upper))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 :  int ) in
    Obj.repr(
# 191 "../../../circe2/src/parser.mly"
                                    ( Syntax.Bin _1 )
# 856 "parser.ml"
               : 'upper))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  int ) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'real_interval) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'power_params) in
    Obj.repr(
# 195 "../../../circe2/src/parser.mly"
                                    (
     let x_min, x_max = _2
     and beta, eta = _3 in
     if beta <= -1.0 then begin
       Printf.eprintf "circe2: ignoring invalid beta: %g <= -1\n" beta;
       flush stderr;
       (_1, Maps.id x_min x_max)
     end else
       let alpha = 1.0 /. (1.0 +. beta) in
       (_1, Maps.power ~alpha ~eta x_min x_max) )
# 874 "parser.ml"
               : 'power))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'beta) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'eta) in
    Obj.repr(
# 208 "../../../circe2/src/parser.mly"
                                    ( (_1, _2) )
# 882 "parser.ml"
               : 'power_params))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'eta) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'beta) in
    Obj.repr(
# 209 "../../../circe2/src/parser.mly"
                                    ( (_2, _1) )
# 890 "parser.ml"
               : 'power_params))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 213 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 897 "parser.ml"
               : 'beta))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 217 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 904 "parser.ml"
               : 'eta))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  int ) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'real_interval) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'resonance_params) in
    Obj.repr(
# 221 "../../../circe2/src/parser.mly"
                                      (
     let x_min, x_max = _2
     and eta, a = _3 in
     (_1, Maps.resonance ~eta ~a x_min x_max) )
# 916 "parser.ml"
               : 'resonance))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'center) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'width) in
    Obj.repr(
# 228 "../../../circe2/src/parser.mly"
                                    ( (_1, _2) )
# 924 "parser.ml"
               : 'resonance_params))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'width) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'center) in
    Obj.repr(
# 229 "../../../circe2/src/parser.mly"
                                    ( (_2, _1) )
# 932 "parser.ml"
               : 'resonance_params))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 233 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 939 "parser.ml"
               : 'center))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'float) in
    Obj.repr(
# 237 "../../../circe2/src/parser.mly"
                                    ( _3 )
# 946 "parser.ml"
               : 'width))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'float_or_int) in
    Obj.repr(
# 241 "../../../circe2/src/parser.mly"
                                    ( _1 )
# 953 "parser.ml"
               : 'float))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'float_or_int) in
    Obj.repr(
# 242 "../../../circe2/src/parser.mly"
                                    ( _1 +. Syntax.epsilon )
# 960 "parser.ml"
               : 'float))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'float_or_int) in
    Obj.repr(
# 243 "../../../circe2/src/parser.mly"
                                    ( _1 -. Syntax.epsilon )
# 967 "parser.ml"
               : 'float))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 247 "../../../circe2/src/parser.mly"
                                    ( float _1 )
# 974 "parser.ml"
               : 'float_or_int))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  float ) in
    Obj.repr(
# 248 "../../../circe2/src/parser.mly"
                                    ( _1 )
# 981 "parser.ml"
               : 'float_or_int))
(* Entry main *)
; (fun __caml_parser_env -> raise (Parsing.YYexit (Parsing.peek_val __caml_parser_env 0)))
|]
let yytables =
  { Parsing.actions=yyact;
    Parsing.transl_const=yytransl_const;
    Parsing.transl_block=yytransl_block;
    Parsing.lhs=yylhs;
    Parsing.len=yylen;
    Parsing.defred=yydefred;
    Parsing.dgoto=yydgoto;
    Parsing.sindex=yysindex;
    Parsing.rindex=yyrindex;
    Parsing.gindex=yygindex;
    Parsing.tablesize=yytablesize;
    Parsing.table=yytable;
    Parsing.check=yycheck;
    Parsing.error_function=parse_error;
    Parsing.names_const=yynames_const;
    Parsing.names_block=yynames_block }
let main (lexfun : Lexing.lexbuf -> token) (lexbuf : Lexing.lexbuf) =
   (Parsing.yyparse yytables 1 lexfun lexbuf :  Syntax.file_cmds list )
