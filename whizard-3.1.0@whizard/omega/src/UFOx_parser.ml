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

open Parsing;;
let _ = parse_error;;
# 31 "../../../omega/src/UFOx_parser.mly"
module X = UFOx_syntax

let parse_error msg =
  raise (UFOx_syntax.Syntax_Error
	   (msg, symbol_start_pos (), symbol_end_pos ()))

let invalid_parameter_attr () =
  parse_error "invalid parameter attribute"

# 30 "UFOx_parser.ml"
let yytransl_const = [|
  261 (* PLUS *);
  262 (* MINUS *);
  263 (* TIMES *);
  264 (* POWER *);
  265 (* DIV *);
  266 (* LPAREN *);
  267 (* RPAREN *);
  268 (* COMMA *);
  269 (* DOT *);
  270 (* END *);
    0|]

let yytransl_block = [|
  257 (* INT *);
  258 (* FLOAT *);
  259 (* ID *);
  260 (* QUOTED *);
    0|]

let yylhs = "\255\255\
\001\000\002\000\002\000\002\000\002\000\002\000\002\000\002\000\
\002\000\002\000\002\000\002\000\002\000\002\000\002\000\002\000\
\002\000\003\000\003\000\000\000"

let yylen = "\002\000\
\002\000\002\000\002\000\001\000\001\000\001\000\001\000\003\000\
\003\000\003\000\003\000\002\000\002\000\003\000\003\000\003\000\
\004\000\001\000\003\000\002\000"

let yydefred = "\000\000\
\000\000\000\000\004\000\005\000\000\000\007\000\000\000\000\000\
\000\000\020\000\000\000\000\000\012\000\002\000\003\000\013\000\
\000\000\000\000\000\000\000\000\000\000\000\000\001\000\016\000\
\000\000\000\000\015\000\000\000\000\000\000\000\014\000\000\000\
\000\000\017\000\019\000"

let yydgoto = "\002\000\
\010\000\025\000\026\000"

let yysindex = "\008\000\
\042\255\000\000\000\000\000\000\005\255\000\000\042\255\052\255\
\042\255\000\000\252\254\021\255\000\000\000\000\000\000\000\000\
\096\255\042\255\042\255\042\255\042\255\042\255\000\000\000\000\
\088\255\003\255\000\000\004\255\004\255\008\255\000\000\008\255\
\042\255\000\000\000\000"

let yyrindex = "\000\000\
\000\000\000\000\000\000\000\000\028\255\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\017\255\000\000\000\000\074\255\078\255\054\255\000\000\064\255\
\000\000\000\000\000\000"

let yygindex = "\000\000\
\000\000\255\255\252\255"

let yytablesize = 107
let yytable = "\011\000\
\018\000\019\000\020\000\021\000\022\000\013\000\016\000\017\000\
\001\000\023\000\020\000\021\000\022\000\034\000\012\000\021\000\
\028\000\029\000\030\000\031\000\032\000\003\000\004\000\005\000\
\006\000\007\000\008\000\018\000\035\000\000\000\009\000\024\000\
\006\000\006\000\006\000\006\000\006\000\000\000\006\000\006\000\
\000\000\006\000\003\000\004\000\005\000\006\000\007\000\008\000\
\000\000\000\000\000\000\009\000\014\000\015\000\005\000\006\000\
\007\000\008\000\010\000\010\000\010\000\009\000\010\000\000\000\
\010\000\010\000\000\000\010\000\011\000\011\000\011\000\000\000\
\011\000\000\000\011\000\011\000\000\000\011\000\008\000\008\000\
\000\000\000\000\009\000\009\000\008\000\008\000\000\000\008\000\
\009\000\009\000\000\000\009\000\018\000\019\000\020\000\021\000\
\022\000\000\000\000\000\033\000\018\000\019\000\020\000\021\000\
\022\000\000\000\027\000"

let yycheck = "\001\000\
\005\001\006\001\007\001\008\001\009\001\007\000\008\000\009\000\
\001\000\014\001\007\001\008\001\009\001\011\001\010\001\008\001\
\018\000\019\000\020\000\021\000\022\000\001\001\002\001\003\001\
\004\001\005\001\006\001\011\001\033\000\255\255\010\001\011\001\
\005\001\006\001\007\001\008\001\009\001\255\255\011\001\012\001\
\255\255\014\001\001\001\002\001\003\001\004\001\005\001\006\001\
\255\255\255\255\255\255\010\001\001\001\002\001\003\001\004\001\
\005\001\006\001\005\001\006\001\007\001\010\001\009\001\255\255\
\011\001\012\001\255\255\014\001\005\001\006\001\007\001\255\255\
\009\001\255\255\011\001\012\001\255\255\014\001\005\001\006\001\
\255\255\255\255\005\001\006\001\011\001\012\001\255\255\014\001\
\011\001\012\001\255\255\014\001\005\001\006\001\007\001\008\001\
\009\001\255\255\255\255\012\001\005\001\006\001\007\001\008\001\
\009\001\255\255\011\001"

let yynames_const = "\
  PLUS\000\
  MINUS\000\
  TIMES\000\
  POWER\000\
  DIV\000\
  LPAREN\000\
  RPAREN\000\
  COMMA\000\
  DOT\000\
  END\000\
  "

let yynames_block = "\
  INT\000\
  FLOAT\000\
  ID\000\
  QUOTED\000\
  "

let yyact = [|
  (fun _ -> failwith "parser")
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 61 "../../../omega/src/UFOx_parser.mly"
            ( _1 )
# 148 "UFOx_parser.ml"
               :  UFOx_syntax.expr ))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 65 "../../../omega/src/UFOx_parser.mly"
                          ( X.integer (- _2) )
# 155 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 :  float ) in
    Obj.repr(
# 66 "../../../omega/src/UFOx_parser.mly"
                          ( X.float (-. _2) )
# 162 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 67 "../../../omega/src/UFOx_parser.mly"
                      ( X.integer _1 )
# 169 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  float ) in
    Obj.repr(
# 68 "../../../omega/src/UFOx_parser.mly"
                      ( X.float _1 )
# 176 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 69 "../../../omega/src/UFOx_parser.mly"
                      ( X.variable _1 )
# 183 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 70 "../../../omega/src/UFOx_parser.mly"
                         ( X.quoted _1 )
# 190 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 71 "../../../omega/src/UFOx_parser.mly"
                      ( X.add _1 _3 )
# 198 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 72 "../../../omega/src/UFOx_parser.mly"
                      ( X.subtract _1 _3 )
# 206 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 73 "../../../omega/src/UFOx_parser.mly"
                      ( X.multiply _1 _3 )
# 214 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 74 "../../../omega/src/UFOx_parser.mly"
                      ( X.divide _1 _3 )
# 222 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 75 "../../../omega/src/UFOx_parser.mly"
                          ( _2 )
# 229 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 76 "../../../omega/src/UFOx_parser.mly"
                          ( X.multiply (X.integer (-1)) _2 )
# 236 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 77 "../../../omega/src/UFOx_parser.mly"
                       ( X.power _1 _3 )
# 244 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 78 "../../../omega/src/UFOx_parser.mly"
                          ( _2 )
# 251 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 :  string ) in
    Obj.repr(
# 79 "../../../omega/src/UFOx_parser.mly"
                          ( X.apply _1 [] )
# 258 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 :  string ) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'args) in
    Obj.repr(
# 80 "../../../omega/src/UFOx_parser.mly"
                          ( X.apply _1 _3 )
# 266 "UFOx_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 84 "../../../omega/src/UFOx_parser.mly"
                   ( [_1] )
# 273 "UFOx_parser.ml"
               : 'args))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'args) in
    Obj.repr(
# 85 "../../../omega/src/UFOx_parser.mly"
                   ( _1 :: _3 )
# 281 "UFOx_parser.ml"
               : 'args))
(* Entry input *)
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
let input (lexfun : Lexing.lexbuf -> token) (lexbuf : Lexing.lexbuf) =
   (Parsing.yyparse yytables 1 lexfun lexbuf :  UFOx_syntax.expr )
