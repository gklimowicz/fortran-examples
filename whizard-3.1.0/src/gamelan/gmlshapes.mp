% gmlshapes.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
vardef polygon primary n =
  if n > 2:
    for i = 1 upto n:
      (.5up rotated (360i/n)) --
    endfor
    cycle
  else:
    fullcircle
  fi
enddef;
vardef polygram primary n =
  if n > 2:
    for i = 1 upto n:
      (.5up rotated (360i/n)) --
      (.2up rotated (360(i+.5)/n)) --
    endfor
    cycle
  else:
    fullcircle
  fi
enddef;
vardef polycross primary n =
  save i;
  for i = 1 upto n:
    origin -- .5 dir (360(i-.5)/n) --
  endfor
  cycle
enddef;
def circle = fullcircle enddef;
def square = (unitsquare shifted -(.5,.5)) enddef;
def triagon = (polygon 3) enddef;
let triangle = triagon;
def tetragon = (polygon 4) enddef;
let diamond = tetragon;
def pentagon = (polygon 5) enddef;
def hexagon = (polygon 6) enddef;
def triagram = (polygram 3) enddef;
def tetragram = (polygram 4) enddef;
def pentagram = (polygram 5) enddef;
def hexagram = (polygram 6) enddef;
def triacross = (polycross 3) enddef;
def tetracross = (polycross 4) enddef;
let cross = tetracross;
def pentacross = (polycross 5) enddef;
def hexacross = (polycross 6) enddef;
vardef arc primary u_ =
  save a_,i_; pair a_;  a_ = 8/360 if pair u_: u_ else (0,u_) fi;
  subpath a_ of fullcircle
enddef;
vardef wedge primary u_ =
  arc u_ -- origin -- cycle
enddef;
vardef doublearc(expr f_) primary u_ =
  (arc u_ scaled f_) -- (reverse arc u_) -- cycle
enddef;
vardef dshape(expr s_)(text d_) = image(draw s_ d_) enddef;
vardef fshape(expr s_)(text f_) = image(fill s_ f_) enddef;
vardef fdshape(expr s_)(text f_)(text d_) =
  image(fill s_ f_; draw s_ d_)
enddef;
vardef ufdshape(expr s_)(text d_) =
  image(unfill s_; draw s_ d_)
enddef;
