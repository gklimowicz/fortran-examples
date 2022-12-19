% graphx.mp --
% Copyright (C) 1995-2011 by ohl@physik.uni-wuerzburg.de
%
% Graphx is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% Graphx is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input graph;
input sarith;
numeric symbol_width, error_width;
symbol_width = 2mm;
error_width = 3mm;
def gMdraw = Gaddto_ doublepath GMdraw_ enddef;
def gMfill = Gaddto_ contour GMfill_ enddef;
def GMdraw_ expr p =
  if path p or pair p:
    GMcvi_ (true) p
  else:
    Gerr_ (p, "gMdraw argument should be a pair or a path")
    origin
  fi
  withpen currentpen Gwithlist_ _op_
enddef;
def GMfill_ expr p =
  if cycle p:
    GMcvi_ (true) p
  else:
    Gerr_ (p, "gMfill argument should be a cyclic path")
    origin .. cycle
  fi
  Gwithlist_ _op_
enddef;
vardef glabel_err (expr p, f) (suffix s) text w =
  if path p or picture p:
    glabel_err_symbol (p, f, s) w;
  elseif string p or numeric p:
    glabel_err_area (p, f, s);
  else:
    Gerr_ (p, "glabel_err's first argument must be a path or numeric")
  fi
enddef;
vardef glabel_err_symbol (expr p, f) (suffix s) text w =
  interim warningcheck := 0;
  save i, Mx, My, Mep, Mem, Mpep, Mpem, Sep, Sem;
  string Sep, Sem;
  glabel (if picture p:
            p
          else:
            image (if f: fill else: draw fi p)
          fi, s1, s2) w;
  (Mx, My) = (Mlog_Str s1, Mlog_Str s2);
  Sep := "0";
  Sem := "0";
  (Mpep, Mpem) = (Mzero, Mzero);
  i := 3;
  forever:
    exitunless iserror s[i];
    if (strip s[i]) Sneq 0:
      if (substring (0,1) of s[i]) = "+":
        Sep := Sep Sadd s[i];
        Mep := Mlog_Str Sep;
        glabel (image (draw bar), s1, s2 Sadd Sep) w;
        gMdraw (Mx, My Madd Mpep) -- (Mx, My Madd Mep) w;
        Mpep := Mep;
      elseif (substring (0,1) of s[i]) = "-":
        Sem := Sem Sadd (Sabs s[i]);
        Mem := Mlog_Str Sem;
        glabel (image (draw bar), s1, s2 Ssub Sem) w;
        gMdraw (Mx, My Msub Mpem) -- (Mx, My Msub Mem) w;
        Mpem := Mem;
      elseif isnumeric (substring (0,1) of s[i]):
        Sem := Sem Sadd (Sabs s[i]);
        Mem := Mlog_Str Sem;
        glabel (image (draw bar), s1, s2 Ssub Sem) w;
        gMdraw (Mx, My Msub Mpem) -- (Mx, My Msub Mem) w;
        Mpem := Mem;
        Sep := Sep Sadd s[i];
        Mep := Mlog_Str Sep;
        glabel (image (draw bar), s1, s2 Sadd Sep) w;
        gMdraw (Mx, My Madd Mpep) -- (Mx, My Madd Mep) w;
        Mpep := Mep;
      elseif (substring (0,1) of s[i]) = "^":
        Sep := Sep Sadd (substring (1,infinity) of s[i]);
        Mep := Mlog_Str Sep;
        glabel (image (fill uparr), s1, s2 Sadd Sep) w;
        gMdraw (Mx, My Madd Mpep) -- (Mx, My Madd Mep) w;
        Mpep := Mep;
      elseif (substring (0,1) of s[i]) = "_":
        Sem := Sem Sadd (Sabs (substring (1,infinity) of s[i]));
        Mem := Mlog_Str Sem;
        glabel (image (fill downarr), s1, s2 Ssub Sem) w;
        gMdraw (Mx, My Msub Mpem) -- (Mx, My Msub Mem) w;
        Mpem := Mem;
      fi
    fi
    i := i + 1;
  endfor
enddef;
vardef issymbolic expr s =
  (s = "+") or (s = "-") or (s = "^") or (s = "_")
enddef;
vardef strip expr s =
  if issymbolic (substring (0,1) of s):
    substring (1,infinity) of s
  else:
    s
  fi
enddef;
vardef isnumeric expr s =
  if known s:
    (isdigit (substring (0,1) of s))
    or (substring (0,1) of s = ".")
  else:
    false
  fi
enddef;
vardef iserror expr s =
  if known s:
    (isnumeric (substring (0,1) of s))
    or (issymbolic (substring (0,1) of s))
  else:
    false
  fi
enddef;
vardef glabel_err_area (expr w, col) (suffix s) =
  interim warningcheck := 0;
  save n, i, Mw, Mx, My, Mep, Mem, Mpep, Mpem, Mxp, Mxm;
  Mw := Mlog_Str w;
  (Mx, My) = (Mlog_Str s1, Mlog_Str s2);
  (Mxm, Mxp) = (Mx Msub Mw, Mx Madd Mw);
  n := 1;
  forever:
    exitunless iserror s[n];
    n := n + 1;
  endfor
  n := n - 1;
  (Mpep, Mpem) = (Mzero, Mzero);
  for i = 3 upto n:
    if s[i] Sneq 0:
      if (substring (0,1) of s[i]) = "+":
        Mep := Mpep Madd Mlog_Str s[i];
        rectangle (Mxm, Mxp, My Madd Mpep, My Madd Mep, col, (i-3)/(n-2));
        Mpep := Mep;
      elseif (substring (0,1) of s[i]) = "-":
        Mem := Mpem Madd (Mabs Mlog_Str s[i]);
        rectangle (Mxm, Mxp, My Msub Mpem, My Msub Mem, col, (i-3)/(n-2));
        Mpem := Mem;
      elseif isnumeric (substring (0,1) of s[i]):
        Mep := Mpep Madd Mlog_Str s[i];
        rectangle (Mxm, Mxp, My Madd Mpep, My Madd Mep, col, (i-3)/(n-2));
        Mpep := Mep;
        Mem := Mpem Madd (Mabs Mlog_Str s[i]);
        rectangle (Mxm, Mxp, My Msub Mpem, My Msub Mem, col, (i-3)/(n-2));
        Mpem := Mem;
      fi
    fi
  endfor
enddef;
vardef rectangle (expr xm, xp, ym, yp, col, f) =
  gMfill (xm, ym) -- (xm, yp) -- (xp, yp) -- (xp, ym) -- cycle
    withcolor
      f[if color col:
          col
        elseif numeric col:
          col*white
        else:
          .5white
        fi, white]
enddef;
vardef gdraw_err_scaled (expr file, scale, p, f, off) text w =
  gdata (file, s,
         gscale (s, off, scale);
         glabel_err (p, f, s) w;)
enddef;
vardef gscale (suffix s) (expr off, scale) =
  save i, pfx;
  string pfx;
  if known off and known s1:
    s[1] := s[1] Sadd off;
  fi
  if known scale:
    i := 2;
    forever:
      exitunless iserror s[i];
      pfx := substring (0,1) of s[i];
      if issymbolic pfx:
        s[i] := (substring (1,infinity) of s[i]) Smul scale;
        s[i] := pfx & s[i];
      else:
        s[i] := s[i] Smul scale;
      fi
      i := i + 1;
    endfor
  fi     
enddef;
vardef gdraw_err (expr file, p, f, off) text w =
  gdraw_err_scaled (file, 1, p, f, off) w
enddef;
vardef bar =
  (-error_width/2,0) -- (error_width/2,0)
enddef;
vardef uparr =
  (-error_width/2,-error_width/2)
  -- origin
  -- (error_width/2,-error_width/2)
  -- (0,-error_width/4)
  -- cycle
enddef;
vardef downarr =
  (-error_width/2,error_width/2)
  -- origin
  -- (error_width/2,error_width/2)
  -- (0,error_width/4)
  -- cycle
enddef;
vardef square =
  unitsquare shifted (-.5,-.5) scaled symbol_width
enddef;
vardef circle =
  fullcircle scaled symbol_width
enddef;
vardef diamond =
  unitsquare shifted (-.5,-.5) rotated 45 scaled symbol_width
enddef;
vardef polygon expr n =
  if n > 2:
    (for i = 1 upto n:
      (.5up rotated (360i/n)) --
    endfor
    cycle)
  else:
    fullcircle
  fi
  scaled symbol_width
enddef;
vardef polygram expr n =
  if n > 2:
    (for i = 1 upto n:
      (.5up rotated (360i/n)) --
      (.2up rotated (360(i+.5)/n)) --
    endfor
    cycle)
  else:
    fullcircle
  fi
  scaled symbol_width
enddef;
vardef polycross expr n =
  save i;
  (for i = 1 upto n:
    origin -- .5 dir (360(i-.5)/n) --
  endfor
  cycle) scaled symbol_width
enddef;
secondarydef a Madd b =
  if a = Mzero:
    b
  elseif b = Mzero:
    a
  elseif a >= b:
    (Mlog(1579 + Mexp(b Mmul (1885.32404-a))) + a-1885.32404)
  else:
    b Madd a
  fi
enddef;
vardef Garwhd_ (expr p) (text w) =
  addto Gcur_ also
    image (fill arrowhead p w; Gsetb_ point infinity of p..cycle);
enddef;
endinput;
% The End.
