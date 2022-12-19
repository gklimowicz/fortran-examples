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

open Parsing;;
let _ = parse_error;;
# 31 "../../../omega/src/vertex_parser.mly"
module T = Vertex_syntax.Token
module E = Vertex_syntax.Expr
module P = Vertex_syntax.Particle
module V = Vertex_syntax.Parameter
module I = Vertex_syntax.Index
module X = Vertex_syntax.Tensor
module F = Vertex_syntax.File_Tree

let parse_error msg =
  raise (Vertex_syntax.Syntax_Error
	   (msg, symbol_start_pos (), symbol_end_pos ()))

let invalid_parameter_attr () =
  parse_error "invalid parameter attribute"

# 61 "vertex_parser.ml"
let yytransl_const = [|
  261 (* SUPER *);
  262 (* SUB *);
  263 (* PRIME *);
  264 (* LBRACE *);
  265 (* RBRACE *);
  266 (* LBRACKET *);
  267 (* RBRACKET *);
  268 (* LPAREN *);
  269 (* RPAREN *);
  270 (* COMMA *);
  271 (* PLUS *);
  272 (* MINUS *);
  273 (* TIMES *);
  274 (* DIV *);
  275 (* EQUAL *);
  277 (* END *);
  278 (* NEUTRAL *);
  279 (* CHARGED *);
  280 (* ANTI *);
  281 (* ALIAS *);
  282 (* TEX *);
  283 (* FORTRAN *);
  284 (* SPIN *);
  285 (* COLOR *);
  286 (* CHARGE *);
  287 (* MASS *);
  288 (* WIDTH *);
  289 (* PARAMETER *);
  290 (* DERIVED *);
  291 (* TENSOR *);
  292 (* INDEX *);
  293 (* FLAVOR *);
  294 (* LORENTZ *);
  295 (* VERTEX *);
    0|]

let yytransl_block = [|
  257 (* DIGIT *);
  258 (* CHAR *);
  259 (* PREFIX *);
  260 (* TOKEN *);
  276 (* INCLUDE *);
    0|]

let yylhs = "\255\255\
\001\000\002\000\002\000\003\000\003\000\003\000\003\000\003\000\
\003\000\004\000\004\000\012\000\012\000\012\000\009\000\009\000\
\011\000\015\000\015\000\017\000\017\000\010\000\010\000\018\000\
\018\000\018\000\018\000\018\000\018\000\018\000\018\000\018\000\
\018\000\018\000\018\000\005\000\005\000\020\000\020\000\021\000\
\021\000\021\000\021\000\021\000\021\000\021\000\021\000\021\000\
\006\000\022\000\022\000\023\000\023\000\023\000\023\000\023\000\
\007\000\024\000\024\000\025\000\025\000\025\000\025\000\025\000\
\008\000\008\000\008\000\008\000\008\000\013\000\013\000\013\000\
\013\000\013\000\013\000\013\000\013\000\013\000\013\000\029\000\
\029\000\019\000\019\000\019\000\027\000\027\000\030\000\030\000\
\030\000\030\000\030\000\016\000\016\000\014\000\028\000\033\000\
\033\000\033\000\033\000\033\000\033\000\033\000\033\000\035\000\
\035\000\036\000\036\000\032\000\032\000\037\000\038\000\038\000\
\034\000\034\000\031\000\031\000\031\000\031\000\031\000\031\000\
\031\000\031\000\031\000\031\000\026\000\026\000\026\000\026\000\
\026\000\026\000\026\000\026\000\026\000\026\000\026\000\000\000"

let yylen = "\002\000\
\002\000\000\000\002\000\001\000\001\000\001\000\001\000\001\000\
\001\000\003\000\003\000\003\000\003\000\003\000\003\000\003\000\
\002\000\003\000\003\000\003\000\003\000\000\000\002\000\002\000\
\003\000\002\000\003\000\002\000\003\000\002\000\002\000\003\000\
\002\000\002\000\002\000\004\000\004\000\000\000\002\000\002\000\
\002\000\002\000\001\000\001\000\001\000\001\000\001\000\001\000\
\003\000\000\000\002\000\002\000\003\000\002\000\003\000\002\000\
\003\000\000\000\002\000\002\000\003\000\002\000\003\000\002\000\
\002\000\003\000\004\000\004\000\002\000\001\000\003\000\003\000\
\003\000\003\000\003\000\003\000\003\000\003\000\002\000\000\000\
\002\000\003\000\003\000\003\000\001\000\002\000\001\000\003\000\
\003\000\004\000\004\000\001\000\002\000\003\000\003\000\000\000\
\001\000\001\000\002\000\002\000\001\000\002\000\002\000\002\000\
\002\000\002\000\002\000\000\000\002\000\001\000\001\000\002\000\
\001\000\001\000\001\000\001\000\001\000\001\000\001\000\001\000\
\001\000\001\000\001\000\001\000\001\000\001\000\001\000\001\000\
\001\000\001\000\001\000\001\000\001\000\001\000\001\000\002\000"

let yydefred = "\000\000\
\000\000\000\000\009\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\136\000\000\000\000\000\004\000\005\000\006\000\
\007\000\008\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\125\000\126\000\127\000\000\000\135\000\000\000\
\134\000\133\000\132\000\128\000\129\000\130\000\131\000\000\000\
\065\000\069\000\001\000\003\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\010\000\000\000\017\000\011\000\000\000\000\000\000\000\000\000\
\000\000\000\000\057\000\000\000\000\000\000\000\000\000\049\000\
\000\000\000\000\000\000\085\000\000\000\000\000\000\000\000\000\
\000\000\000\000\066\000\109\000\015\000\016\000\115\000\116\000\
\117\000\000\000\123\000\124\000\122\000\118\000\119\000\120\000\
\121\000\000\000\087\000\000\000\000\000\000\000\024\000\026\000\
\028\000\030\000\000\000\031\000\000\000\033\000\034\000\035\000\
\023\000\000\000\043\000\000\000\000\000\000\000\044\000\045\000\
\046\000\047\000\048\000\036\000\000\000\037\000\060\000\000\000\
\062\000\000\000\064\000\059\000\052\000\000\000\054\000\000\000\
\056\000\051\000\093\000\018\000\019\000\000\000\013\000\012\000\
\000\000\000\000\000\000\000\000\014\000\086\000\000\000\079\000\
\113\000\114\000\000\000\067\000\068\000\000\000\000\000\000\000\
\000\000\094\000\000\000\000\000\000\000\110\000\025\000\027\000\
\029\000\000\000\032\000\082\000\083\000\084\000\040\000\041\000\
\042\000\039\000\061\000\063\000\053\000\055\000\073\000\072\000\
\071\000\074\000\000\000\000\000\077\000\078\000\081\000\095\000\
\088\000\089\000\000\000\105\000\104\000\107\000\106\000\112\000\
\099\000\100\000\103\000\102\000\020\000\021\000\090\000\091\000"

let yydgoto = "\002\000\
\011\000\012\000\013\000\014\000\015\000\016\000\017\000\018\000\
\020\000\057\000\022\000\040\000\078\000\074\000\041\000\075\000\
\109\000\058\000\151\000\124\000\125\000\072\000\073\000\067\000\
\068\000\042\000\079\000\080\000\152\000\098\000\099\000\047\000\
\162\000\155\000\163\000\164\000\165\000\166\000"

let yysindex = "\005\000\
\020\255\000\000\000\000\036\255\036\255\036\255\036\255\036\255\
\036\255\034\000\000\000\002\255\020\255\000\000\000\000\000\000\
\000\000\000\000\025\255\190\000\036\255\190\000\040\255\040\255\
\031\255\056\255\000\000\000\000\000\000\025\255\000\000\046\255\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\049\255\
\000\000\000\000\000\000\000\000\025\255\254\254\088\000\073\255\
\059\255\059\255\059\255\040\255\080\255\040\255\059\255\059\255\
\000\000\190\000\000\000\000\000\046\255\199\000\199\000\080\255\
\080\255\059\255\000\000\031\255\080\255\080\255\059\255\000\000\
\056\255\025\255\015\255\000\000\046\255\161\000\050\255\040\255\
\018\255\012\255\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\025\255\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\066\255\000\000\059\255\059\255\059\255\000\000\000\000\
\000\000\000\000\025\255\000\000\059\255\000\000\000\000\000\000\
\000\000\172\000\000\000\059\255\059\255\059\255\000\000\000\000\
\000\000\000\000\000\000\000\000\199\000\000\000\000\000\059\255\
\000\000\059\255\000\000\000\000\000\000\059\255\000\000\059\255\
\000\000\000\000\000\000\000\000\000\000\128\000\000\000\000\000\
\046\255\046\255\046\255\046\255\000\000\000\000\040\255\000\000\
\000\000\000\000\066\255\000\000\000\000\043\255\052\000\070\000\
\069\255\000\000\081\255\103\255\081\255\000\000\000\000\000\000\
\000\000\054\255\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\118\255\118\255\000\000\000\000\000\000\000\000\
\000\000\000\000\029\255\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000"

let yyrindex = "\000\000\
\071\255\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\071\255\000\000\000\000\000\000\
\000\000\000\000\106\000\220\255\000\000\220\255\000\000\000\000\
\227\255\247\255\000\000\000\000\000\000\106\000\000\000\119\255\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\106\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\220\255\000\000\000\000\119\255\254\255\254\255\000\000\
\000\000\000\000\000\000\227\255\000\000\000\000\000\000\000\000\
\247\255\218\255\000\000\000\000\119\255\000\000\139\000\150\000\
\000\000\106\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\106\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\137\255\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\106\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\254\255\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\119\255\119\255\119\255\119\255\000\000\000\000\150\000\000\000\
\000\000\000\000\117\000\000\000\000\000\106\000\000\000\000\000\
\116\255\000\000\158\255\179\255\200\255\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\183\000\192\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\
\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000"

let yygindex = "\000\000\
\000\000\096\000\000\000\000\000\000\000\000\000\000\000\000\000\
\075\000\239\255\000\000\000\000\213\255\240\255\217\255\184\255\
\042\000\000\000\241\255\197\255\000\000\028\000\000\000\058\000\
\000\000\000\000\000\000\000\000\248\255\253\255\000\000\224\255\
\008\000\000\000\001\000\235\255\004\000\242\255"

let yytablesize = 487
let yytable = "\081\000\
\083\000\139\000\046\000\126\000\060\000\001\000\085\000\062\000\
\063\000\103\000\104\000\105\000\084\000\108\000\045\000\111\000\
\112\000\114\000\086\000\153\000\156\000\154\000\043\000\140\000\
\127\000\129\000\131\000\045\000\081\000\133\000\135\000\137\000\
\157\000\142\000\170\000\141\000\106\000\207\000\110\000\003\000\
\113\000\004\000\005\000\019\000\081\000\045\000\076\000\061\000\
\045\000\208\000\150\000\193\000\006\000\007\000\008\000\009\000\
\082\000\077\000\010\000\064\000\167\000\168\000\169\000\194\000\
\205\000\178\000\030\000\065\000\066\000\171\000\159\000\160\000\
\161\000\158\000\206\000\161\000\175\000\176\000\177\000\021\000\
\023\000\024\000\025\000\026\000\069\000\195\000\160\000\030\000\
\179\000\107\000\180\000\002\000\070\000\071\000\181\000\059\000\
\182\000\100\000\101\000\102\000\138\000\187\000\188\000\189\000\
\190\000\128\000\130\000\159\000\044\000\161\000\134\000\136\000\
\081\000\081\000\081\000\081\000\111\000\111\000\111\000\111\000\
\108\000\111\000\108\000\111\000\111\000\132\000\111\000\111\000\
\111\000\111\000\111\000\111\000\111\000\111\000\147\000\148\000\
\111\000\096\000\096\000\096\000\096\000\201\000\191\000\204\000\
\096\000\096\000\200\000\096\000\096\000\096\000\096\000\096\000\
\096\000\096\000\096\000\197\000\199\000\096\000\097\000\097\000\
\097\000\097\000\192\000\000\000\202\000\097\000\097\000\203\000\
\097\000\097\000\097\000\097\000\097\000\097\000\097\000\097\000\
\000\000\000\000\097\000\098\000\098\000\098\000\098\000\000\000\
\000\000\000\000\098\000\098\000\000\000\098\000\098\000\098\000\
\098\000\098\000\098\000\098\000\098\000\000\000\000\000\098\000\
\101\000\101\000\101\000\101\000\000\000\000\000\000\000\101\000\
\101\000\000\000\101\000\101\000\101\000\101\000\101\000\101\000\
\101\000\101\000\108\000\108\000\101\000\108\000\000\000\000\000\
\000\000\108\000\092\000\000\000\092\000\108\000\108\000\108\000\
\108\000\108\000\108\000\108\000\000\000\000\000\092\000\022\000\
\022\000\022\000\022\000\000\000\000\000\000\000\058\000\058\000\
\058\000\058\000\000\000\000\000\022\000\022\000\022\000\022\000\
\000\000\000\000\022\000\058\000\058\000\058\000\058\000\000\000\
\000\000\058\000\050\000\050\000\050\000\050\000\000\000\000\000\
\000\000\038\000\038\000\038\000\038\000\000\000\000\000\050\000\
\050\000\050\000\050\000\000\000\000\000\050\000\038\000\038\000\
\038\000\038\000\027\000\028\000\038\000\029\000\000\000\000\000\
\000\000\030\000\031\000\032\000\033\000\000\000\034\000\035\000\
\036\000\037\000\038\000\039\000\087\000\088\000\000\000\089\000\
\000\000\000\000\000\000\090\000\196\000\000\000\000\000\091\000\
\092\000\093\000\094\000\095\000\096\000\097\000\087\000\088\000\
\000\000\089\000\000\000\000\000\000\000\090\000\198\000\000\000\
\000\000\091\000\092\000\093\000\094\000\095\000\096\000\097\000\
\087\000\088\000\000\000\089\000\000\000\000\000\000\000\090\000\
\000\000\000\000\000\000\091\000\092\000\093\000\094\000\095\000\
\096\000\097\000\108\000\108\000\000\000\108\000\000\000\000\000\
\000\000\108\000\000\000\000\000\000\000\108\000\108\000\108\000\
\108\000\108\000\108\000\108\000\096\000\096\000\000\000\096\000\
\000\000\096\000\000\000\096\000\096\000\096\000\096\000\000\000\
\183\000\096\000\184\000\000\000\185\000\000\000\145\000\146\000\
\147\000\148\000\000\000\070\000\186\000\070\000\000\000\070\000\
\000\000\070\000\070\000\070\000\070\000\000\000\080\000\070\000\
\080\000\000\000\080\000\000\000\080\000\080\000\080\000\080\000\
\000\000\143\000\080\000\144\000\000\000\000\000\000\000\145\000\
\146\000\147\000\148\000\000\000\172\000\149\000\173\000\000\000\
\000\000\000\000\145\000\146\000\147\000\148\000\000\000\075\000\
\174\000\075\000\000\000\075\000\000\000\075\000\075\000\000\000\
\076\000\000\000\076\000\075\000\076\000\000\000\076\000\076\000\
\000\000\000\000\000\000\000\000\076\000\048\000\049\000\050\000\
\051\000\052\000\053\000\054\000\055\000\056\000\115\000\116\000\
\117\000\118\000\119\000\120\000\121\000\122\000\123\000"

let yycheck = "\032\000\
\040\000\074\000\019\000\063\000\022\000\001\000\009\001\023\000\
\024\000\049\000\050\000\051\000\045\000\053\000\003\001\055\000\
\056\000\061\000\021\001\002\001\009\001\004\001\021\001\009\001\
\064\000\065\000\066\000\003\001\061\000\069\000\070\000\071\000\
\021\001\077\000\107\000\021\001\052\000\009\001\054\000\020\001\
\058\000\022\001\023\001\008\001\077\000\003\001\001\001\008\001\
\003\001\021\001\001\001\009\001\033\001\034\001\035\001\036\001\
\008\001\012\001\039\001\029\001\100\000\101\000\102\000\021\001\
\011\001\125\000\008\001\037\001\038\001\109\000\005\001\006\001\
\007\001\090\000\021\001\007\001\116\000\117\000\118\000\005\000\
\006\000\007\000\008\000\009\000\029\001\158\000\006\001\008\001\
\128\000\010\001\130\000\021\001\037\001\038\001\134\000\021\000\
\136\000\025\001\026\001\027\001\073\000\145\000\146\000\147\000\
\148\000\064\000\065\000\005\001\013\000\007\001\069\000\070\000\
\145\000\146\000\147\000\148\000\001\001\002\001\003\001\004\001\
\002\001\006\001\004\001\008\001\009\001\068\000\011\001\012\001\
\013\001\014\001\015\001\016\001\017\001\018\001\017\001\018\001\
\021\001\001\001\002\001\003\001\004\001\163\000\151\000\165\000\
\008\001\009\001\161\000\011\001\012\001\013\001\014\001\015\001\
\016\001\017\001\018\001\159\000\160\000\021\001\001\001\002\001\
\003\001\004\001\155\000\255\255\164\000\008\001\009\001\164\000\
\011\001\012\001\013\001\014\001\015\001\016\001\017\001\018\001\
\255\255\255\255\021\001\001\001\002\001\003\001\004\001\255\255\
\255\255\255\255\008\001\009\001\255\255\011\001\012\001\013\001\
\014\001\015\001\016\001\017\001\018\001\255\255\255\255\021\001\
\001\001\002\001\003\001\004\001\255\255\255\255\255\255\008\001\
\009\001\255\255\011\001\012\001\013\001\014\001\015\001\016\001\
\017\001\018\001\001\001\002\001\021\001\004\001\255\255\255\255\
\255\255\008\001\009\001\255\255\011\001\012\001\013\001\014\001\
\015\001\016\001\017\001\018\001\255\255\255\255\021\001\020\001\
\021\001\022\001\023\001\255\255\255\255\255\255\020\001\021\001\
\022\001\023\001\255\255\255\255\033\001\034\001\035\001\036\001\
\255\255\255\255\039\001\033\001\034\001\035\001\036\001\255\255\
\255\255\039\001\020\001\021\001\022\001\023\001\255\255\255\255\
\255\255\020\001\021\001\022\001\023\001\255\255\255\255\033\001\
\034\001\035\001\036\001\255\255\255\255\039\001\033\001\034\001\
\035\001\036\001\001\001\002\001\039\001\004\001\255\255\255\255\
\255\255\008\001\009\001\010\001\011\001\255\255\013\001\014\001\
\015\001\016\001\017\001\018\001\001\001\002\001\255\255\004\001\
\255\255\255\255\255\255\008\001\009\001\255\255\255\255\012\001\
\013\001\014\001\015\001\016\001\017\001\018\001\001\001\002\001\
\255\255\004\001\255\255\255\255\255\255\008\001\009\001\255\255\
\255\255\012\001\013\001\014\001\015\001\016\001\017\001\018\001\
\001\001\002\001\255\255\004\001\255\255\255\255\255\255\008\001\
\255\255\255\255\255\255\012\001\013\001\014\001\015\001\016\001\
\017\001\018\001\001\001\002\001\255\255\004\001\255\255\255\255\
\255\255\008\001\255\255\255\255\255\255\012\001\013\001\014\001\
\015\001\016\001\017\001\018\001\008\001\009\001\255\255\011\001\
\255\255\013\001\255\255\015\001\016\001\017\001\018\001\255\255\
\009\001\021\001\011\001\255\255\013\001\255\255\015\001\016\001\
\017\001\018\001\255\255\009\001\021\001\011\001\255\255\013\001\
\255\255\015\001\016\001\017\001\018\001\255\255\009\001\021\001\
\011\001\255\255\013\001\255\255\015\001\016\001\017\001\018\001\
\255\255\009\001\021\001\011\001\255\255\255\255\255\255\015\001\
\016\001\017\001\018\001\255\255\009\001\021\001\011\001\255\255\
\255\255\255\255\015\001\016\001\017\001\018\001\255\255\009\001\
\021\001\011\001\255\255\013\001\255\255\015\001\016\001\255\255\
\009\001\255\255\011\001\021\001\013\001\255\255\015\001\016\001\
\255\255\255\255\255\255\255\255\021\001\024\001\025\001\026\001\
\027\001\028\001\029\001\030\001\031\001\032\001\024\001\025\001\
\026\001\027\001\028\001\029\001\030\001\031\001\032\001"

let yynames_const = "\
  SUPER\000\
  SUB\000\
  PRIME\000\
  LBRACE\000\
  RBRACE\000\
  LBRACKET\000\
  RBRACKET\000\
  LPAREN\000\
  RPAREN\000\
  COMMA\000\
  PLUS\000\
  MINUS\000\
  TIMES\000\
  DIV\000\
  EQUAL\000\
  END\000\
  NEUTRAL\000\
  CHARGED\000\
  ANTI\000\
  ALIAS\000\
  TEX\000\
  FORTRAN\000\
  SPIN\000\
  COLOR\000\
  CHARGE\000\
  MASS\000\
  WIDTH\000\
  PARAMETER\000\
  DERIVED\000\
  TENSOR\000\
  INDEX\000\
  FLAVOR\000\
  LORENTZ\000\
  VERTEX\000\
  "

let yynames_block = "\
  DIGIT\000\
  CHAR\000\
  PREFIX\000\
  TOKEN\000\
  INCLUDE\000\
  "

let yyact = [|
  (fun _ -> failwith "parser")
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'declarations) in
    Obj.repr(
# 75 "../../../omega/src/vertex_parser.mly"
                    ( _1 )
# 422 "vertex_parser.ml"
               :  Vertex_syntax.File_Tree.t ))
; (fun __caml_parser_env ->
    Obj.repr(
# 79 "../../../omega/src/vertex_parser.mly"
                                   ( [] )
# 428 "vertex_parser.ml"
               : 'declarations))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'declaration) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'declarations) in
    Obj.repr(
# 80 "../../../omega/src/vertex_parser.mly"
                                   ( _1 :: _2 )
# 436 "vertex_parser.ml"
               : 'declarations))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'particle) in
    Obj.repr(
# 84 "../../../omega/src/vertex_parser.mly"
                      ( F.Particle _1 )
# 443 "vertex_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'parameter) in
    Obj.repr(
# 85 "../../../omega/src/vertex_parser.mly"
                      ( F.Parameter _1 )
# 450 "vertex_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'index) in
    Obj.repr(
# 86 "../../../omega/src/vertex_parser.mly"
                      ( F.Index _1 )
# 457 "vertex_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'tensor) in
    Obj.repr(
# 87 "../../../omega/src/vertex_parser.mly"
                      ( F.Tensor _1 )
# 464 "vertex_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'vertex) in
    Obj.repr(
# 88 "../../../omega/src/vertex_parser.mly"
                      ( let e, t = _1 in
			F.Vertex (e, t) )
# 472 "vertex_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 90 "../../../omega/src/vertex_parser.mly"
                      ( F.Include _1 )
# 479 "vertex_parser.ml"
               : 'declaration))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'particle_attributes) in
    Obj.repr(
# 95 "../../../omega/src/vertex_parser.mly"
     ( { P.name = P.Neutral _2; P.attr = _3 } )
# 487 "vertex_parser.ml"
               : 'particle))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_arg_pair) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'particle_attributes) in
    Obj.repr(
# 97 "../../../omega/src/vertex_parser.mly"
     ( let p, ap = _2 in
       { P.name = P.Charged (p, ap); P.attr = _3 } )
# 496 "vertex_parser.ml"
               : 'particle))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 102 "../../../omega/src/vertex_parser.mly"
                          ( _2 )
# 503 "vertex_parser.ml"
               : 'expr_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 103 "../../../omega/src/vertex_parser.mly"
                          ( parse_error "expected `]', found `}'" )
# 510 "vertex_parser.ml"
               : 'expr_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 104 "../../../omega/src/vertex_parser.mly"
                          ( parse_error "missing `]'" )
# 517 "vertex_parser.ml"
               : 'expr_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'scripted_token) in
    Obj.repr(
# 108 "../../../omega/src/vertex_parser.mly"
                                  ( _2 )
# 524 "vertex_parser.ml"
               : 'token_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'scripted_token) in
    Obj.repr(
# 109 "../../../omega/src/vertex_parser.mly"
                                  ( parse_error "missing `}'" )
# 531 "vertex_parser.ml"
               : 'token_arg))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'token_arg) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_arg) in
    Obj.repr(
# 113 "../../../omega/src/vertex_parser.mly"
                       ( (_1, _2) )
# 539 "vertex_parser.ml"
               : 'token_arg_pair))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list) in
    Obj.repr(
# 117 "../../../omega/src/vertex_parser.mly"
                              ( _2 )
# 546 "vertex_parser.ml"
               : 'token_list_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list) in
    Obj.repr(
# 118 "../../../omega/src/vertex_parser.mly"
                              ( parse_error "missing `}'" )
# 553 "vertex_parser.ml"
               : 'token_list_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list) in
    Obj.repr(
# 124 "../../../omega/src/vertex_parser.mly"
                                  ( _2 )
# 560 "vertex_parser.ml"
               : 'token_list_opt_arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list) in
    Obj.repr(
# 125 "../../../omega/src/vertex_parser.mly"
                                  ( parse_error "missing `}'" )
# 567 "vertex_parser.ml"
               : 'token_list_opt_arg))
; (fun __caml_parser_env ->
    Obj.repr(
# 129 "../../../omega/src/vertex_parser.mly"
                                          ( [ ] )
# 573 "vertex_parser.ml"
               : 'particle_attributes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'particle_attribute) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'particle_attributes) in
    Obj.repr(
# 130 "../../../omega/src/vertex_parser.mly"
                                          ( _1 :: _2 )
# 581 "vertex_parser.ml"
               : 'particle_attributes))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 134 "../../../omega/src/vertex_parser.mly"
                                             ( P.Alias _2 )
# 588 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 135 "../../../omega/src/vertex_parser.mly"
                                             ( P.Alias _3 )
# 595 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 136 "../../../omega/src/vertex_parser.mly"
                                             ( P.TeX _2 )
# 602 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 137 "../../../omega/src/vertex_parser.mly"
                                             ( P.TeX_Anti _3 )
# 609 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 138 "../../../omega/src/vertex_parser.mly"
                                             ( P.Fortran _2 )
# 616 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 139 "../../../omega/src/vertex_parser.mly"
                                             ( P.Fortran_Anti _3 )
# 623 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'arg) in
    Obj.repr(
# 140 "../../../omega/src/vertex_parser.mly"
                                             ( P.Spin _2 )
# 630 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 141 "../../../omega/src/vertex_parser.mly"
                                                       ( P.Color ([], _2) )
# 637 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list_opt_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 142 "../../../omega/src/vertex_parser.mly"
                                                       ( P.Color (_2, _3) )
# 645 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'arg) in
    Obj.repr(
# 143 "../../../omega/src/vertex_parser.mly"
                                             ( P.Charge _2 )
# 652 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 144 "../../../omega/src/vertex_parser.mly"
                                             ( P.Mass _2 )
# 659 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 145 "../../../omega/src/vertex_parser.mly"
                                             ( P.Width _2 )
# 666 "vertex_parser.ml"
               : 'particle_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'token_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'arg) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_attributes) in
    Obj.repr(
# 150 "../../../omega/src/vertex_parser.mly"
     ( V.Parameter { V.name = _2; V.value = _3; V.attr = _4 } )
# 675 "vertex_parser.ml"
               : 'parameter))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'token_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'arg) in
    let _4 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_attributes) in
    Obj.repr(
# 152 "../../../omega/src/vertex_parser.mly"
     ( V.Derived { V.name = _2; V.value = _3; V.attr = _4 } )
# 684 "vertex_parser.ml"
               : 'parameter))
; (fun __caml_parser_env ->
    Obj.repr(
# 156 "../../../omega/src/vertex_parser.mly"
                                            ( [ ] )
# 690 "vertex_parser.ml"
               : 'parameter_attributes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'parameter_attribute) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'parameter_attributes) in
    Obj.repr(
# 157 "../../../omega/src/vertex_parser.mly"
                                            ( _1 :: _2 )
# 698 "vertex_parser.ml"
               : 'parameter_attributes))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 161 "../../../omega/src/vertex_parser.mly"
                          ( V.Alias _2 )
# 705 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 162 "../../../omega/src/vertex_parser.mly"
                          ( V.TeX _2 )
# 712 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 163 "../../../omega/src/vertex_parser.mly"
                          ( V.Fortran _2 )
# 719 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    Obj.repr(
# 164 "../../../omega/src/vertex_parser.mly"
                          ( invalid_parameter_attr () )
# 725 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    Obj.repr(
# 165 "../../../omega/src/vertex_parser.mly"
                          ( invalid_parameter_attr () )
# 731 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    Obj.repr(
# 166 "../../../omega/src/vertex_parser.mly"
                          ( invalid_parameter_attr () )
# 737 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    Obj.repr(
# 167 "../../../omega/src/vertex_parser.mly"
                          ( invalid_parameter_attr () )
# 743 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    Obj.repr(
# 168 "../../../omega/src/vertex_parser.mly"
                          ( invalid_parameter_attr () )
# 749 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    Obj.repr(
# 169 "../../../omega/src/vertex_parser.mly"
                          ( invalid_parameter_attr () )
# 755 "vertex_parser.ml"
               : 'parameter_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'index_attributes) in
    Obj.repr(
# 173 "../../../omega/src/vertex_parser.mly"
                                    ( { I.name = _2; I.attr = _3 } )
# 763 "vertex_parser.ml"
               : 'index))
; (fun __caml_parser_env ->
    Obj.repr(
# 177 "../../../omega/src/vertex_parser.mly"
                                    ( [ ] )
# 769 "vertex_parser.ml"
               : 'index_attributes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'index_attribute) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'index_attributes) in
    Obj.repr(
# 178 "../../../omega/src/vertex_parser.mly"
                                    ( _1 :: _2 )
# 777 "vertex_parser.ml"
               : 'index_attributes))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 182 "../../../omega/src/vertex_parser.mly"
                                             ( I.Color ([], _2) )
# 784 "vertex_parser.ml"
               : 'index_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list_opt_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 183 "../../../omega/src/vertex_parser.mly"
                                             ( I.Color (_2, _3) )
# 792 "vertex_parser.ml"
               : 'index_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 184 "../../../omega/src/vertex_parser.mly"
                                             ( I.Flavor ([], _2) )
# 799 "vertex_parser.ml"
               : 'index_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list_opt_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 185 "../../../omega/src/vertex_parser.mly"
                                             ( I.Flavor (_2, _3) )
# 807 "vertex_parser.ml"
               : 'index_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 186 "../../../omega/src/vertex_parser.mly"
                                             ( I.Lorentz _2 )
# 814 "vertex_parser.ml"
               : 'index_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'tensor_attributes) in
    Obj.repr(
# 190 "../../../omega/src/vertex_parser.mly"
                                      ( { X.name = _2; X.attr = _3 } )
# 822 "vertex_parser.ml"
               : 'tensor))
; (fun __caml_parser_env ->
    Obj.repr(
# 194 "../../../omega/src/vertex_parser.mly"
                                      ( [ ] )
# 828 "vertex_parser.ml"
               : 'tensor_attributes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'tensor_attribute) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'tensor_attributes) in
    Obj.repr(
# 195 "../../../omega/src/vertex_parser.mly"
                                      ( _1 :: _2 )
# 836 "vertex_parser.ml"
               : 'tensor_attributes))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 199 "../../../omega/src/vertex_parser.mly"
                                             ( X.Color ([], _2) )
# 843 "vertex_parser.ml"
               : 'tensor_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list_opt_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 200 "../../../omega/src/vertex_parser.mly"
                                             ( X.Color (_2, _3) )
# 851 "vertex_parser.ml"
               : 'tensor_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 201 "../../../omega/src/vertex_parser.mly"
                                             ( X.Flavor ([], _2) )
# 858 "vertex_parser.ml"
               : 'tensor_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token_list_opt_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 202 "../../../omega/src/vertex_parser.mly"
                                             ( X.Flavor (_2, _3) )
# 866 "vertex_parser.ml"
               : 'tensor_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 203 "../../../omega/src/vertex_parser.mly"
                                             ( X.Lorentz _2 )
# 873 "vertex_parser.ml"
               : 'tensor_attribute))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 207 "../../../omega/src/vertex_parser.mly"
                                   ( (E.integer 1, T.list _2) )
# 880 "vertex_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr_arg) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'token_list_arg) in
    Obj.repr(
# 208 "../../../omega/src/vertex_parser.mly"
                                   ( (_2, T.list _3) )
# 888 "vertex_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'expr_arg) in
    Obj.repr(
# 209 "../../../omega/src/vertex_parser.mly"
                                   ( (_2, T.list []) )
# 895 "vertex_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'expr_arg) in
    Obj.repr(
# 210 "../../../omega/src/vertex_parser.mly"
                                   ( parse_error "missing `}'" )
# 902 "vertex_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'not_arg_or_token_list) in
    Obj.repr(
# 211 "../../../omega/src/vertex_parser.mly"
                                   ( parse_error "expected `[' or `{'" )
# 909 "vertex_parser.ml"
               : 'vertex))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'integer) in
    Obj.repr(
# 217 "../../../omega/src/vertex_parser.mly"
                            ( E.integer _1 )
# 916 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 218 "../../../omega/src/vertex_parser.mly"
                            ( _2 )
# 923 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 219 "../../../omega/src/vertex_parser.mly"
                              ( parse_error "expected `)', found `]'" )
# 930 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 220 "../../../omega/src/vertex_parser.mly"
                            ( parse_error "expected `)', found `}'" )
# 937 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 221 "../../../omega/src/vertex_parser.mly"
                         ( parse_error "missing `)'" )
# 944 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 222 "../../../omega/src/vertex_parser.mly"
                            ( E.add _1 _3 )
# 952 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 223 "../../../omega/src/vertex_parser.mly"
                            ( E.sub _1 _3 )
# 960 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 224 "../../../omega/src/vertex_parser.mly"
                            ( E.mult _1 _3 )
# 968 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'expr) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'expr) in
    Obj.repr(
# 225 "../../../omega/src/vertex_parser.mly"
                            ( E.div _1 _3 )
# 976 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'bare_scripted_token) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'arg_list) in
    Obj.repr(
# 226 "../../../omega/src/vertex_parser.mly"
                                ( E.apply _1 _2 )
# 984 "vertex_parser.ml"
               : 'expr))
; (fun __caml_parser_env ->
    Obj.repr(
# 233 "../../../omega/src/vertex_parser.mly"
                           ( [] )
# 990 "vertex_parser.ml"
               : 'arg_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'arg) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'arg_list) in
    Obj.repr(
# 234 "../../../omega/src/vertex_parser.mly"
                           ( _1 :: _2 )
# 998 "vertex_parser.ml"
               : 'arg_list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 238 "../../../omega/src/vertex_parser.mly"
                        ( _2 )
# 1005 "vertex_parser.ml"
               : 'arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 239 "../../../omega/src/vertex_parser.mly"
                        ( parse_error "expected `}', found `]'" )
# 1012 "vertex_parser.ml"
               : 'arg))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'expr) in
    Obj.repr(
# 240 "../../../omega/src/vertex_parser.mly"
                        ( parse_error "missing `}'" )
# 1019 "vertex_parser.ml"
               : 'arg))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 244 "../../../omega/src/vertex_parser.mly"
                   ( _1 )
# 1026 "vertex_parser.ml"
               : 'integer))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'integer) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 245 "../../../omega/src/vertex_parser.mly"
                   ( 10 * _1 + _2 )
# 1034 "vertex_parser.ml"
               : 'integer))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'bare_token) in
    Obj.repr(
# 249 "../../../omega/src/vertex_parser.mly"
                                           ( _1 )
# 1041 "vertex_parser.ml"
               : 'token))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'scripted_token) in
    Obj.repr(
# 250 "../../../omega/src/vertex_parser.mly"
                                           ( _2 )
# 1048 "vertex_parser.ml"
               : 'token))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'scripted_token) in
    Obj.repr(
# 251 "../../../omega/src/vertex_parser.mly"
                                           ( parse_error "missing `}'" )
# 1055 "vertex_parser.ml"
               : 'token))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'scripted_token) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'token_list) in
    Obj.repr(
# 252 "../../../omega/src/vertex_parser.mly"
                                           ( T.list (_2 :: _3) )
# 1063 "vertex_parser.ml"
               : 'token))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 2 : 'scripted_token) in
    let _3 = (Parsing.peek_val __caml_parser_env 1 : 'token_list) in
    Obj.repr(
# 253 "../../../omega/src/vertex_parser.mly"
                                           ( parse_error "missing `}'" )
# 1071 "vertex_parser.ml"
               : 'token))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'scripted_token) in
    Obj.repr(
# 260 "../../../omega/src/vertex_parser.mly"
                             ( [_1] )
# 1078 "vertex_parser.ml"
               : 'token_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'scripted_token) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token_list) in
    Obj.repr(
# 261 "../../../omega/src/vertex_parser.mly"
                             ( _1 :: _2 )
# 1086 "vertex_parser.ml"
               : 'token_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'prefixes) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'token) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'optional_scripts) in
    Obj.repr(
# 265 "../../../omega/src/vertex_parser.mly"
                                   ( T.scripted _1 _2 _3 )
# 1095 "vertex_parser.ml"
               : 'scripted_token))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 2 : 'prefixes) in
    let _2 = (Parsing.peek_val __caml_parser_env 1 : 'name) in
    let _3 = (Parsing.peek_val __caml_parser_env 0 : 'optional_scripts) in
    Obj.repr(
# 269 "../../../omega/src/vertex_parser.mly"
                                   ( T.scripted _1 _2 _3 )
# 1104 "vertex_parser.ml"
               : 'bare_scripted_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 273 "../../../omega/src/vertex_parser.mly"
                  ( (None, None) )
# 1110 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'super) in
    Obj.repr(
# 274 "../../../omega/src/vertex_parser.mly"
                  ( (_1, None) )
# 1117 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'sub) in
    Obj.repr(
# 275 "../../../omega/src/vertex_parser.mly"
                  ( (None, _1) )
# 1124 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'super) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'sub) in
    Obj.repr(
# 276 "../../../omega/src/vertex_parser.mly"
                  ( (_1, _2) )
# 1132 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'sub) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'super) in
    Obj.repr(
# 277 "../../../omega/src/vertex_parser.mly"
                  ( (_2, _1) )
# 1140 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'primes) in
    Obj.repr(
# 278 "../../../omega/src/vertex_parser.mly"
                  ( (_1, None) )
# 1147 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'primes) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'sub) in
    Obj.repr(
# 279 "../../../omega/src/vertex_parser.mly"
                  ( (_1, _2) )
# 1155 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 : 'sub) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'primes) in
    Obj.repr(
# 280 "../../../omega/src/vertex_parser.mly"
                  ( (_2, _1) )
# 1163 "vertex_parser.ml"
               : 'optional_scripts))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token) in
    Obj.repr(
# 284 "../../../omega/src/vertex_parser.mly"
                ( Some _2 )
# 1170 "vertex_parser.ml"
               : 'super))
; (fun __caml_parser_env ->
    Obj.repr(
# 285 "../../../omega/src/vertex_parser.mly"
                ( parse_error "superscript can't start with `}'" )
# 1176 "vertex_parser.ml"
               : 'super))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'token) in
    Obj.repr(
# 291 "../../../omega/src/vertex_parser.mly"
                ( Some _2 )
# 1183 "vertex_parser.ml"
               : 'sub))
; (fun __caml_parser_env ->
    Obj.repr(
# 292 "../../../omega/src/vertex_parser.mly"
                ( parse_error "subscript can't start with `}'" )
# 1189 "vertex_parser.ml"
               : 'sub))
; (fun __caml_parser_env ->
    Obj.repr(
# 298 "../../../omega/src/vertex_parser.mly"
                    ( [] )
# 1195 "vertex_parser.ml"
               : 'prefixes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 1 :  string ) in
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'prefixes) in
    Obj.repr(
# 299 "../../../omega/src/vertex_parser.mly"
                    ( _1 :: _2 )
# 1203 "vertex_parser.ml"
               : 'prefixes))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 : 'prime_list) in
    Obj.repr(
# 303 "../../../omega/src/vertex_parser.mly"
                ( Some (T.list _1) )
# 1210 "vertex_parser.ml"
               : 'primes))
; (fun __caml_parser_env ->
    Obj.repr(
# 307 "../../../omega/src/vertex_parser.mly"
                    ( [T.token "\\prime"] )
# 1216 "vertex_parser.ml"
               : 'prime_list))
; (fun __caml_parser_env ->
    let _2 = (Parsing.peek_val __caml_parser_env 0 : 'prime_list) in
    Obj.repr(
# 308 "../../../omega/src/vertex_parser.mly"
                    ( T.token "\\prime" :: _2 )
# 1223 "vertex_parser.ml"
               : 'prime_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 312 "../../../omega/src/vertex_parser.mly"
            ( T.token _1 )
# 1230 "vertex_parser.ml"
               : 'name))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 313 "../../../omega/src/vertex_parser.mly"
            ( T.token _1 )
# 1237 "vertex_parser.ml"
               : 'name))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 317 "../../../omega/src/vertex_parser.mly"
            ( T.digit _1 )
# 1244 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 318 "../../../omega/src/vertex_parser.mly"
            ( T.token _1 )
# 1251 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 319 "../../../omega/src/vertex_parser.mly"
            ( T.token _1 )
# 1258 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 320 "../../../omega/src/vertex_parser.mly"
            ( T.token "+" )
# 1264 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 321 "../../../omega/src/vertex_parser.mly"
            ( T.token "-" )
# 1270 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 322 "../../../omega/src/vertex_parser.mly"
            ( T.token "*" )
# 1276 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 323 "../../../omega/src/vertex_parser.mly"
            ( T.token "/" )
# 1282 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 324 "../../../omega/src/vertex_parser.mly"
            ( T.token "," )
# 1288 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 325 "../../../omega/src/vertex_parser.mly"
            ( T.token "(" )
# 1294 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    Obj.repr(
# 326 "../../../omega/src/vertex_parser.mly"
            ( T.token ")" )
# 1300 "vertex_parser.ml"
               : 'bare_token))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  int ) in
    Obj.repr(
# 330 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1307 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 331 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1314 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    let _1 = (Parsing.peek_val __caml_parser_env 0 :  string ) in
    Obj.repr(
# 332 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1321 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 333 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1327 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 334 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1333 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 335 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1339 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 336 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1345 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 337 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1351 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 338 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1357 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 339 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1363 "vertex_parser.ml"
               : 'not_arg_or_token_list))
; (fun __caml_parser_env ->
    Obj.repr(
# 340 "../../../omega/src/vertex_parser.mly"
            ( () )
# 1369 "vertex_parser.ml"
               : 'not_arg_or_token_list))
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
   (Parsing.yyparse yytables 1 lexfun lexbuf :  Vertex_syntax.File_Tree.t )
