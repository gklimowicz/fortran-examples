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

open Parsing;;
let _ = parse_error;;
# 31 "../../../omega/src/UFO_parser.mly"
module U = UFO_syntax

let parse_error msg =
  raise (UFO_syntax.Syntax_Error
	   (msg, symbol_start_pos (), symbol_end_pos ()))

let invalid_parameter_attr () =
  parse_error "invalid parameter attribute"

# 34 "UFO_parser.ml"
let yytransl_const = [|
  261 (* DOT *);
  262 (* COMMA *);
  263 (* COLON *);
  264 (* EQUAL *);
  265 (* PLUS *);
  266 (* MINUS *);
  267 (* DIV *);
  268 (* LPAREN *);
  269 (* RPAREN *);
  270 (* LBRACE *);
  271 (* RBRACE *);
  272 (* LBRACKET *);
  273 (* RBRACKET *);
  274 (* END *);
    0|]

let yytransl_block = [|
  257 (* INT *);
  258 (* FLOAT *);
  259 (* STRING *);
  260 (* ID *);
    0|]

let yylhs = "\255\255\
\001\000\002\000\002\000\003\000\003\000\003\000\003\000\004\000\
\004\000\005\000\005\000\007\000\007\000\007\000\008\000\008\000\
\008\000\008\000\008\000\008\000\009\000\009\000\009\000\009\000\
\010\000\010\000\010\000\012\000\012\000\014\000\014\000\006\000\
\006\000\018\000\018\000\019\000\019\000\013\000\013\000\011\000\
\011\000\015\000\015\000\020\000\016\000\016\000\021\000\017\000\
\017\000\022\000\000\000"

let yylen = "\002\000\
\002\000\000\000\002\000\005\000\006\000\003\000\003\000\001\000\
\003\000\001\000\003\000\003\000\003\000\003\000\001\000\003\000\
\001\000\001\000\001\000\001\000\002\000\003\000\003\000\003\000\
\003\000\003\000\003\000\001\000\003\000\001\000\003\000\001\000\
\001\000\003\000\003\000\003\000\003\000\001\000\003\000\001\000\
\003\000\001\000\003\000\003\000\001\000\003\000\007\000\001\000\
\003\000\005\000\002\000"

let yydefred = "\000\000\
\000\000\000\000\000\000\051\000\000\000\000\000\000\000\001\000\
\003\000\000\000\008\000\000\000\007\000\000\000\032\000\033\000\
\000\000\000\000\000\000\000\000\009\000\040\000\000\000\037\000\
\000\000\000\000\004\000\000\000\000\000\041\000\000\000\035\000\
\000\000\005\000\000\000\000\000\017\000\000\000\000\000\000\000\
\019\000\012\000\013\000\014\000\000\000\011\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\021\000\000\000\000\000\000\000\000\000\000\000\016\000\000\000\
\000\000\000\000\025\000\026\000\027\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\022\000\023\000\024\000\044\000\
\000\000\000\000\043\000\000\000\046\000\000\000\049\000\031\000\
\029\000\039\000\000\000\000\000\000\000\050\000\000\000\000\000"

let yydgoto = "\002\000\
\004\000\005\000\006\000\058\000\028\000\013\000\029\000\042\000\
\043\000\044\000\059\000\066\000\061\000\062\000\050\000\051\000\
\052\000\015\000\016\000\053\000\054\000\055\000"

let yysindex = "\011\000\
\018\255\000\000\044\255\000\000\047\255\018\255\035\255\000\000\
\000\000\000\000\000\000\041\255\000\000\050\255\000\000\000\000\
\064\255\058\255\011\255\060\255\000\000\000\000\049\255\000\000\
\050\255\062\255\000\000\056\255\065\255\000\000\049\255\000\000\
\005\255\000\000\068\255\063\255\000\000\031\255\001\255\049\255\
\000\000\000\000\000\000\000\000\050\255\000\000\072\255\069\255\
\016\255\066\255\067\255\070\255\071\255\073\255\074\255\077\255\
\000\000\061\255\051\255\075\255\076\255\078\255\000\000\083\255\
\080\255\081\255\000\000\000\000\000\000\084\255\079\255\085\255\
\088\255\086\255\093\255\095\255\000\000\000\000\000\000\000\000\
\098\255\094\255\000\000\099\255\000\000\086\255\000\000\000\000\
\000\000\000\000\089\255\100\255\097\255\000\000\086\255\101\255"

let yyrindex = "\000\000\
\057\255\000\000\000\000\000\000\000\000\057\255\000\000\000\000\
\000\000\023\255\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\024\255\000\000\000\000\000\000\092\255\000\000\027\255\000\000\
\000\000\000\000\000\000\010\255\000\000\000\000\000\000\042\255\
\000\000\000\000\000\000\000\000\043\255\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\096\255\102\255\103\255\090\255\
\000\000\034\255\091\255\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\029\255"

let yygindex = "\000\000\
\000\000\072\000\000\000\249\255\074\000\248\255\000\000\000\000\
\000\000\000\000\252\255\218\255\035\000\039\000\043\000\044\000\
\042\000\000\000\096\000\000\000\000\000\000\000"

let yytablesize = 118
let yytable = "\012\000\
\060\000\056\000\014\000\022\000\011\000\036\000\037\000\022\000\
\011\000\024\000\023\000\001\000\031\000\025\000\026\000\015\000\
\065\000\057\000\038\000\011\000\039\000\003\000\015\000\027\000\
\041\000\040\000\006\000\036\000\045\000\036\000\034\000\040\000\
\034\000\048\000\047\000\089\000\036\000\010\000\011\000\034\000\
\006\000\036\000\049\000\047\000\034\000\017\000\028\000\020\000\
\018\000\018\000\028\000\007\000\019\000\017\000\020\000\018\000\
\075\000\018\000\020\000\076\000\022\000\011\000\030\000\011\000\
\008\000\017\000\074\000\021\000\034\000\033\000\035\000\026\000\
\063\000\047\000\002\000\064\000\070\000\009\000\071\000\072\000\
\067\000\068\000\073\000\080\000\069\000\081\000\048\000\096\000\
\056\000\011\000\084\000\077\000\078\000\082\000\079\000\022\000\
\086\000\030\000\091\000\065\000\092\000\093\000\094\000\095\000\
\010\000\017\000\030\000\038\000\046\000\090\000\042\000\088\000\
\083\000\087\000\085\000\032\000\045\000\048\000"

let yycheck = "\007\000\
\039\000\001\001\007\000\003\001\004\001\001\001\002\001\003\001\
\004\001\018\000\018\000\001\000\020\000\018\000\004\001\006\001\
\001\001\017\001\014\001\004\001\016\001\004\001\013\001\013\001\
\033\000\033\000\004\001\004\001\033\000\006\001\004\001\009\001\
\006\001\003\001\006\001\074\000\013\001\003\001\004\001\013\001\
\018\001\018\001\012\001\015\001\018\001\005\001\013\001\006\001\
\006\001\009\001\017\001\008\001\012\001\005\001\013\001\013\001\
\006\001\009\001\009\001\009\001\003\001\004\001\003\001\004\001\
\018\001\005\001\006\001\004\001\013\001\008\001\006\001\004\001\
\001\001\011\001\018\001\007\001\006\001\006\000\006\001\006\001\
\015\001\015\001\006\001\001\001\015\001\006\001\003\001\095\000\
\001\001\004\001\012\001\017\001\017\001\013\001\017\001\003\001\
\012\001\003\001\001\001\001\001\007\001\013\001\003\001\007\001\
\013\001\005\001\017\001\017\001\035\000\075\000\015\001\073\000\
\070\000\072\000\071\000\020\000\015\001\015\001"

let yynames_const = "\
  DOT\000\
  COMMA\000\
  COLON\000\
  EQUAL\000\
  PLUS\000\
  MINUS\000\
  DIV\000\
  LPAREN\000\
  RPAREN\000\
  LBRACE\000\
  RBRACE\000\
  LBRACKET\000\
  RBRACKET\000\
  END\000\
  "

let yynames_block = "\
  INT\000\
  FLOAT\000\
  STRING\000\
  ID\000\
  "

let yyact = [|
  (fun _ -> failwith "parser")
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'declarations) in
    Obj.repr(
# 59 "../../../omega/src/UFO_parser.mly"
                    ( _1 )
# 195 "UFO_parser.ml"
               :  UFO_syntax.t ))
; (fun __caml_parser_env ->
    Obj.repr(
# 63 "../../../omega/src/UFO_parser.mly"
                            ( [] )
# 201 "UFO_parser.ml"
               : 'declarations))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'declaration) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'declarations) in
    Obj.repr(
# 64 "../../../omega/src/UFO_parser.mly"
                            ( _1 :: _2 )
# 209 "UFO_parser.ml"
               : 'declarations))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 4 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 2 : 'name) in
    Obj.repr(
# 68 "../../../omega/src/UFO_parser.mly"
                                          ( { U.name = _1;
					      U.kind = _3;
					      U.attribs = [] } )
# 219 "UFO_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 5 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 3 : 'name) in
    let _5 = (Parsing.peek_val __caml_parser_env 1 : 'attributes) in
    Obj.repr(
# 71 "../../../omega/src/UFO_parser.mly"
                                          ( { U.name = _1;
					      U.kind = _3;
					      U.attribs = _5 } )
# 230 "UFO_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 74 "../../../omega/src/UFO_parser.mly"
                                          ( U.macro _1 (U.String _3) )
# 238 "UFO_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string_expr) in
    Obj.repr(
# 75 "../../../omega/src/UFO_parser.mly"
                                          ( U.macro _1 (U.String_Expr _3) )
# 246 "UFO_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 79 "../../../omega/src/UFO_parser.mly"
               ( [_1] )
# 253 "UFO_parser.ml"
               : 'name))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'name) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 80 "../../../omega/src/UFO_parser.mly"
               ( _3 :: _1 )
# 261 "UFO_parser.ml"
               : 'name))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'attribute) in
    Obj.repr(
# 84 "../../../omega/src/UFO_parser.mly"
                              ( [_1] )
# 268 "UFO_parser.ml"
               : 'attributes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'attribute) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'attributes) in
    Obj.repr(
# 85 "../../../omega/src/UFO_parser.mly"
                              ( _1 :: _3 )
# 276 "UFO_parser.ml"
               : 'attributes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'value) in
    Obj.repr(
# 89 "../../../omega/src/UFO_parser.mly"
                       ( { U.a_name = _1; U.a_value = _3 } )
# 284 "UFO_parser.ml"
               : 'attribute))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'list) in
    Obj.repr(
# 90 "../../../omega/src/UFO_parser.mly"
                       ( { U.a_name = _1; U.a_value = _3 } )
# 292 "UFO_parser.ml"
               : 'attribute))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'dictionary) in
    Obj.repr(
# 91 "../../../omega/src/UFO_parser.mly"
                       ( { U.a_name = _1; U.a_value = _3 } )
# 300 "UFO_parser.ml"
               : 'attribute))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 95 "../../../omega/src/UFO_parser.mly"
               ( U.Integer _1 )
# 307 "UFO_parser.ml"
               : 'value))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  int ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 96 "../../../omega/src/UFO_parser.mly"
               ( U.Fraction (_1, _3) )
# 315 "UFO_parser.ml"
               : 'value))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  float ) in
    Obj.repr(
# 97 "../../../omega/src/UFO_parser.mly"
               ( U.Float _1 )
# 322 "UFO_parser.ml"
               : 'value))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'string) in
    Obj.repr(
# 98 "../../../omega/src/UFO_parser.mly"
               ( U.String _1 )
# 329 "UFO_parser.ml"
               : 'value))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'string_expr) in
    Obj.repr(
# 99 "../../../omega/src/UFO_parser.mly"
               ( U.String_Expr _1 )
# 336 "UFO_parser.ml"
               : 'value))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'name) in
    Obj.repr(
# 100 "../../../omega/src/UFO_parser.mly"
               ( U.Name _1 )
# 343 "UFO_parser.ml"
               : 'value))
; (fun __caml_parser_env ->
    Obj.repr(
# 104 "../../../omega/src/UFO_parser.mly"
                            ( U.Empty_List )
# 349 "UFO_parser.ml"
               : 'list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'names) in
    Obj.repr(
# 105 "../../../omega/src/UFO_parser.mly"
                              ( U.Name_List _2 )
# 356 "UFO_parser.ml"
               : 'list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'strings) in
    Obj.repr(
# 106 "../../../omega/src/UFO_parser.mly"
                              ( U.String_List _2 )
# 363 "UFO_parser.ml"
               : 'list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'integers) in
    Obj.repr(
# 107 "../../../omega/src/UFO_parser.mly"
                              ( U.Integer_List _2 )
# 370 "UFO_parser.ml"
               : 'list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'orders) in
    Obj.repr(
# 111 "../../../omega/src/UFO_parser.mly"
                           ( U.Order_Dictionary _2 )
# 377 "UFO_parser.ml"
               : 'dictionary))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'couplings) in
    Obj.repr(
# 112 "../../../omega/src/UFO_parser.mly"
                           ( U.Coupling_Dictionary _2 )
# 384 "UFO_parser.ml"
               : 'dictionary))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'decays) in
    Obj.repr(
# 113 "../../../omega/src/UFO_parser.mly"
                           ( U.Decay_Dictionary _2 )
# 391 "UFO_parser.ml"
               : 'dictionary))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'name) in
    Obj.repr(
# 117 "../../../omega/src/UFO_parser.mly"
                    ( [_1] )
# 398 "UFO_parser.ml"
               : 'names))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'name) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'names) in
    Obj.repr(
# 118 "../../../omega/src/UFO_parser.mly"
                    ( _1 :: _3 )
# 406 "UFO_parser.ml"
               : 'names))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 122 "../../../omega/src/UFO_parser.mly"
                      ( [_1] )
# 413 "UFO_parser.ml"
               : 'integers))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  int ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'integers) in
    Obj.repr(
# 123 "../../../omega/src/UFO_parser.mly"
                      ( _1 :: _3 )
# 421 "UFO_parser.ml"
               : 'integers))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'literal_string_expr) in
    Obj.repr(
# 132 "../../../omega/src/UFO_parser.mly"
                       ( _1 )
# 428 "UFO_parser.ml"
               : 'string_expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'macro_string_expr) in
    Obj.repr(
# 133 "../../../omega/src/UFO_parser.mly"
                       ( _1 )
# 435 "UFO_parser.ml"
               : 'string_expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'name) in
    Obj.repr(
# 137 "../../../omega/src/UFO_parser.mly"
                                 ( [U.Literal _1; U.Macro _3] )
# 443 "UFO_parser.ml"
               : 'literal_string_expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'macro_string_expr) in
    Obj.repr(
# 138 "../../../omega/src/UFO_parser.mly"
                                 ( U.Literal _1 :: _3 )
# 451 "UFO_parser.ml"
               : 'literal_string_expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'name) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string) in
    Obj.repr(
# 142 "../../../omega/src/UFO_parser.mly"
                                 ( [U.Macro _1; U.Literal _3] )
# 459 "UFO_parser.ml"
               : 'macro_string_expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'name) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string_expr) in
    Obj.repr(
# 143 "../../../omega/src/UFO_parser.mly"
                                 ( U.Macro _1 :: _3 )
# 467 "UFO_parser.ml"
               : 'macro_string_expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'string) in
    Obj.repr(
# 147 "../../../omega/src/UFO_parser.mly"
                        ( [_1] )
# 474 "UFO_parser.ml"
               : 'strings))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'strings) in
    Obj.repr(
# 148 "../../../omega/src/UFO_parser.mly"
                        ( _1 :: _3 )
# 482 "UFO_parser.ml"
               : 'strings))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 152 "../../../omega/src/UFO_parser.mly"
                      ( _1 )
# 489 "UFO_parser.ml"
               : 'string))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'string) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 153 "../../../omega/src/UFO_parser.mly"
                      ( _1 ^ _3 )
# 497 "UFO_parser.ml"
               : 'string))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'order) in
    Obj.repr(
# 157 "../../../omega/src/UFO_parser.mly"
                      ( [_1] )
# 504 "UFO_parser.ml"
               : 'orders))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'order) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'orders) in
    Obj.repr(
# 158 "../../../omega/src/UFO_parser.mly"
                      ( _1 :: _3 )
# 512 "UFO_parser.ml"
               : 'orders))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 162 "../../../omega/src/UFO_parser.mly"
                    ( (_1, _3) )
# 520 "UFO_parser.ml"
               : 'order))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'coupling) in
    Obj.repr(
# 166 "../../../omega/src/UFO_parser.mly"
                            ( [_1] )
# 527 "UFO_parser.ml"
               : 'couplings))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'coupling) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'couplings) in
    Obj.repr(
# 167 "../../../omega/src/UFO_parser.mly"
                            ( _1 :: _3 )
# 535 "UFO_parser.ml"
               : 'couplings))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 5 :  int ) in
    let _4 = (Parsing.peek_val __caml_parser_env 3 :  int ) in
    let _7 = (Parsing.peek_val __caml_parser_env 0 : 'name) in
    Obj.repr(
# 171 "../../../omega/src/UFO_parser.mly"
                                          ( (_2, _4, _7) )
# 544 "UFO_parser.ml"
               : 'coupling))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'decay) in
    Obj.repr(
# 175 "../../../omega/src/UFO_parser.mly"
                      ( [_1] )
# 551 "UFO_parser.ml"
               : 'decays))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'decay) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'decays) in
    Obj.repr(
# 176 "../../../omega/src/UFO_parser.mly"
                      ( _1 :: _3 )
# 559 "UFO_parser.ml"
               : 'decays))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 3 : 'names) in
    let _5 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 180 "../../../omega/src/UFO_parser.mly"
                                    ( (_2, _5) )
# 567 "UFO_parser.ml"
               : 'decay))
(* Entry file *)
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
let file (lexfun : Lexing.lexbuf -> token) (lexbuf : Lexing.lexbuf) =
   (Parsing.yyparse yytables 1 lexfun lexbuf :  UFO_syntax.t )
