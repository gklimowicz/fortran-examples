define([foreach], [pushdef([$1])_foreach($@)popdef([$1])])
define([_arg1], [$1])
define([_foreach], [ifelse([$2], [()], [],
   [define([$1], _arg1$2)$3[]$0([$1], (shift$2), [$3])])])

define([forloop], [ifelse(eval([($2) <= ($3)]), [1],
  [pushdef([$1])_$0([$1], eval([$2]),
    eval([$3]), [$4])popdef([$1])])])
define([_forloop],
  [define([$1], [$2])$4[]ifelse([$2], [$3], [],
    [$0([$1], incr([$2]), [$3], [$4])])])


define([foreachq],
  [ifelse([$2], [], [], [_$0([$1], [$3], $2)])])
define([_foreachq],
  [pushdef([$1], forloop([$1], [3], [$#],
  [$0_([1], [2], indir([$1]))])[popdef(
    [$1])])indir([$1], $@)])
define([_foreachq_],
  [[define([$$1], [$$3])$$2[]]])


define([_cat],[$1$2])
