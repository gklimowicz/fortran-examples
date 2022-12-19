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

open Parsing;;
let _ = parse_error;;
# 26 "../../../omega/src/cascade_parser.mly"
open Cascade_syntax
let parse_error msg =
  raise (Syntax_Error (msg, symbol_start (), symbol_end ()))
# 26 "cascade_parser.ml"
let yytransl_const = [|
  259 (* LPAREN *);
  260 (* RPAREN *);
  261 (* LBRACKET *);
  262 (* RBRACKET *);
  263 (* AND *);
  264 (* PLUS *);
  265 (* COLON *);
  266 (* COMMA *);
  267 (* NOT *);
  268 (* HAT *);
  269 (* ONSHELL *);
  270 (* OFFSHELL *);
  271 (* GAUSS *);
  272 (* END *);
    0|]

let yytransl_block = [|
  257 (* STRING *);
  258 (* INT *);
    0|]

let yylhs = "\255\255\
\001\000\001\000\002\000\002\000\002\000\002\000\002\000\003\000\
\004\000\004\000\004\000\004\000\005\000\005\000\005\000\005\000\
\005\000\005\000\005\000\008\000\008\000\009\000\006\000\006\000\
\007\000\007\000\000\000"

let yylen = "\002\000\
\001\000\002\000\001\000\001\000\001\000\003\000\003\000\002\000\
\002\000\004\000\004\000\005\000\001\000\003\000\004\000\003\000\
\004\000\003\000\004\000\001\000\003\000\001\000\001\000\003\000\
\001\000\003\000\002\000"

let yydefred = "\000\000\
\000\000\000\000\022\000\000\000\000\000\000\000\001\000\027\000\
\000\000\003\000\004\000\005\000\000\000\020\000\000\000\023\000\
\000\000\000\000\000\000\000\000\002\000\000\000\000\000\000\000\
\000\000\006\000\000\000\000\000\000\000\000\000\007\000\021\000\
\000\000\000\000\000\000\000\000\000\000\000\000\024\000\011\000\
\000\000\010\000\000\000\000\000\000\000\000\000\000\000\012\000"

let yydgoto = "\002\000\
\008\000\009\000\010\000\011\000\012\000\028\000\029\000\013\000\
\014\000"

let yysindex = "\005\000\
\000\255\000\000\000\000\056\255\007\255\061\255\000\000\000\000\
\047\255\000\000\000\000\000\000\057\255\000\000\020\255\000\000\
\254\254\007\255\064\255\056\255\000\000\058\255\014\255\022\255\
\045\255\000\000\028\255\254\254\068\255\003\255\000\000\000\000\
\007\255\254\254\007\255\254\254\007\255\254\254\000\000\000\000\
\007\255\000\000\069\255\254\254\254\254\254\254\254\254\000\000"

let yyrindex = "\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\006\255\000\000\000\000\000\000\
\010\255\000\000\027\255\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\070\255\000\000\000\000\000\000\000\000\
\000\000\031\255\000\000\033\255\000\000\035\255\000\000\000\000\
\000\000\000\000\000\000\037\255\041\255\048\255\071\255\000\000"

let yygindex = "\000\000\
\000\000\001\000\000\000\000\000\000\000\251\255\020\000\000\000\
\039\000"

let yytablesize = 81
let yytable = "\017\000\
\019\000\003\000\004\000\016\000\015\000\001\000\027\000\016\000\
\042\000\013\000\005\000\006\000\013\000\008\000\016\000\007\000\
\008\000\034\000\036\000\038\000\031\000\013\000\016\000\026\000\
\033\000\008\000\020\000\044\000\039\000\045\000\009\000\046\000\
\035\000\009\000\014\000\047\000\016\000\014\000\018\000\016\000\
\015\000\018\000\009\000\015\000\017\000\016\000\014\000\017\000\
\016\000\043\000\018\000\019\000\015\000\020\000\019\000\037\000\
\017\000\003\000\004\000\003\000\032\000\016\000\021\000\019\000\
\022\000\018\000\005\000\006\000\030\000\023\000\024\000\025\000\
\027\000\040\000\048\000\025\000\026\000\041\000\041\000\025\000\
\026\000"

let yycheck = "\005\000\
\006\000\002\001\003\001\001\001\004\000\001\000\009\001\001\001\
\006\001\004\001\011\001\012\001\007\001\004\001\001\001\016\001\
\007\001\023\000\024\000\025\000\020\000\016\001\001\001\004\001\
\011\001\016\001\007\001\033\000\001\001\035\000\004\001\037\000\
\011\001\007\001\004\001\041\000\004\001\007\001\004\001\007\001\
\004\001\007\001\016\001\007\001\004\001\001\001\016\001\007\001\
\016\001\030\000\016\001\004\001\016\001\007\001\007\001\011\001\
\016\001\002\001\003\001\002\001\022\000\001\001\016\001\016\001\
\008\001\005\001\011\001\012\001\005\001\013\001\014\001\015\001\
\009\001\006\001\006\001\006\001\006\001\010\001\010\001\010\001\
\010\001"

let yynames_const = "\
  LPAREN\000\
  RPAREN\000\
  LBRACKET\000\
  RBRACKET\000\
  AND\000\
  PLUS\000\
  COLON\000\
  COMMA\000\
  NOT\000\
  HAT\000\
  ONSHELL\000\
  OFFSHELL\000\
  GAUSS\000\
  END\000\
  "

let yynames_block = "\
  STRING\000\
  INT\000\
  "

let yyact = [|
  (fun _ -> failwith "parser")
; (fun __caml_parser_env ->
    Obj.repr(
# 47 "../../../omega/src/cascade_parser.mly"
                                    ( mk_true () )
# 148 "cascade_parser.ml"
               :  (string, int list, string) Cascade_syntax.t ))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'cascades) in
    Obj.repr(
# 48 "../../../omega/src/cascade_parser.mly"
                                    ( _1 )
# 155 "cascade_parser.ml"
               :  (string, int list, string) Cascade_syntax.t ))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'exclusion) in
    Obj.repr(
# 52 "../../../omega/src/cascade_parser.mly"
                                    ( _1 )
# 162 "cascade_parser.ml"
               : 'cascades))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'vertex) in
    Obj.repr(
# 53 "../../../omega/src/cascade_parser.mly"
                                    ( _1 )
# 169 "cascade_parser.ml"
               : 'cascades))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'cascade) in
    Obj.repr(
# 54 "../../../omega/src/cascade_parser.mly"
                                    ( _1 )
# 176 "cascade_parser.ml"
               : 'cascades))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'cascades) in
    Obj.repr(
# 55 "../../../omega/src/cascade_parser.mly"
                                    ( _2 )
# 183 "cascade_parser.ml"
               : 'cascades))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'cascades) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'cascades) in
    Obj.repr(
# 56 "../../../omega/src/cascade_parser.mly"
                                    ( mk_and _1 _3 )
# 191 "cascade_parser.ml"
               : 'cascades))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 60 "../../../omega/src/cascade_parser.mly"
                                    ( mk_x_flavor _2 )
# 198 "cascade_parser.ml"
               : 'exclusion))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 64 "../../../omega/src/cascade_parser.mly"
                                    ( mk_x_vertex _2 [] )
# 205 "cascade_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'string_list) in
    Obj.repr(
# 66 "../../../omega/src/cascade_parser.mly"
                                    ( mk_x_vertex _2 [] )
# 212 "cascade_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'string_lists) in
    Obj.repr(
# 68 "../../../omega/src/cascade_parser.mly"
                                    ( mk_x_vertex [] _3 )
# 219 "cascade_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 3 : 'string_list) in
    let _4 = (Parsing.peek_val __caml_parser_env 1 : 'string_lists) in
    Obj.repr(
# 70 "../../../omega/src/cascade_parser.mly"
                                    ( mk_x_vertex _2 _4 )
# 227 "cascade_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'momentum_list) in
    Obj.repr(
# 74 "../../../omega/src/cascade_parser.mly"
                                    ( mk_any_flavor _1 )
# 234 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'momentum_list) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 76 "../../../omega/src/cascade_parser.mly"
                                    ( mk_on_shell _3 _1 )
# 242 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'momentum_list) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 78 "../../../omega/src/cascade_parser.mly"
                                    ( mk_on_shell_not _4 _1 )
# 250 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'momentum_list) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 80 "../../../omega/src/cascade_parser.mly"
                                    ( mk_off_shell _3 _1 )
# 258 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'momentum_list) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 82 "../../../omega/src/cascade_parser.mly"
                                    ( mk_off_shell_not _4 _1 )
# 266 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'momentum_list) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 83 "../../../omega/src/cascade_parser.mly"
                                    ( mk_gauss _3 _1 )
# 274 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 3 : 'momentum_list) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 85 "../../../omega/src/cascade_parser.mly"
                                    ( mk_gauss_not _4 _1 )
# 282 "cascade_parser.ml"
               : 'cascade))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'momentum) in
    Obj.repr(
# 89 "../../../omega/src/cascade_parser.mly"
                                    ( [_1] )
# 289 "cascade_parser.ml"
               : 'momentum_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'momentum_list) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'momentum) in
    Obj.repr(
# 90 "../../../omega/src/cascade_parser.mly"
                                    ( _3 :: _1 )
# 297 "cascade_parser.ml"
               : 'momentum_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 94 "../../../omega/src/cascade_parser.mly"
                                    ( _1 )
# 304 "cascade_parser.ml"
               : 'momentum))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 98 "../../../omega/src/cascade_parser.mly"
                                    ( [_1] )
# 311 "cascade_parser.ml"
               : 'string_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'string_list) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 99 "../../../omega/src/cascade_parser.mly"
                                    ( _3 :: _1 )
# 319 "cascade_parser.ml"
               : 'string_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 103 "../../../omega/src/cascade_parser.mly"
                                    ( [_1] )
# 326 "cascade_parser.ml"
               : 'string_lists))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'string_lists) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'string_list) in
    Obj.repr(
# 104 "../../../omega/src/cascade_parser.mly"
                                    ( _3 :: _1 )
# 334 "cascade_parser.ml"
               : 'string_lists))
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
   (Parsing.yyparse yytables 1 lexfun lexbuf :  (string, int list, string) Cascade_syntax.t )
