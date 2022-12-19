% gmlaux.mp -- supplements to plain.mp (including string.mp)
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
vardef take expr i of s = 
  if string s:      substring(0,i) of 
  elseif path s:    subpath(0,i) of
  elseif picture s: subpicture(0,i) of
  fi s
enddef;
vardef drop expr i of s = 
  if string s:      substring(i,infinity) of 
  elseif path s:    subpath(i,infinity) of
  elseif picture s: subpicture(i,infinity) of
  fi s 
enddef;
vardef subpicture expr u of p =
  image(save i,x,y; i=0; z=(xpart u, ypart u);
        for q within p: 
          if (x<(incr i)) and (y>=i): addto currentpicture also q; fi 
        endfor)
enddef;
def number_from expr s =
  scantokens
    if substring(0,1) of s = "[":  substring(1, length s - 1) of fi  s
enddef;
vardef whatever_str = save ?; string ?; ? enddef;
def ?? = whatever enddef;
string TAB,CR; TAB=char 9; CR=char 13;
vardef thesymbol(expr s) =
  if picture s: s
  elseif path s: 
    image(addto currentpicture doublepath s withpen currentpen _op_)
  else:  s infont defaultfont scaled defaultscale
  fi
enddef;
vardef thelabel @#(expr s,z) =  % Position s near z
  save p; picture p; p=thesymbol(s);
  p shifted (z + labs @# *
       (labeloffset*laboff @#
        - (labxf @#*lrcorner p + labyf @#*ulcorner p 
           + (1-labxf @#-labyf @#)*llcorner p)))
enddef;
pair laboff.o; laboff.o=origin; labxf.o=labyf.o=0;
labs.o=0; 
labs=labs.lft=labs.rt=labs.bot=labs.top=
labs.ulft=labs.urt=labs.llft=labs.lrt=1;
def dotlabel = _dotlabel enddef;
vardef _dotlabel@#(expr s,z) =
  interim linecap:=rounded;
  draw z withpen pencircle scaled dotscale;
  label@#(s,z);
enddef;
newinternal dotscale; dotscale:=3bp;
vardef thetransform(expr p) =
  save t; transform t;
  if picture p:
    point 0 of ttempl_ transformed t = llcorner p;
    point 1 of ttempl_ transformed t = lrcorner p;
    point 2 of ttempl_ transformed t = urcorner p;
  elseif path p:
    for i=0,1,2:
      point i of ttempl_ transformed t = point i of p;
    endfor
  fi
  t
enddef;
path ttempl_; ttempl_=(0,0)--(1,0)--(1,1);
def penscale = pickup pencircle scaled enddef;
def withpenscale = withpen pencircle scaled enddef;
vardef colorpart primary p =
  (redpart p, greenpart p, bluepart p)
enddef;
def withoptions_from primary p =
  withpen penpart p withcolor colorpart p dashed dashpart p
enddef;
primarydef p colored c =
  image(addto currentpicture also p withcolor c)
enddef;
def savecurrentpen =
  pen _cpen; _cpen=currentpen;
  save currentpen; pen currentpen; pickup _cpen;
enddef;
numeric thick, thin; thin=1pt; thick=2thin;
color yellow, magenta, cyan;
yellow = red+green; magenta = red+blue; cyan = blue+green;
vardef gray primary nn = (nn/100)*white enddef;
def grey = gray enddef;
vardef newcolor @# expr s_ = 
  if not color @#: scantokens("color " & templateof @# & ";"); fi
  if color s_: @#:=s_; elseif string s_: @#:=color_from_string s_; fi
enddef;
vardef color_from_string expr s =
  save l_; l_=length s;
  if l_=6: 
    (hex substring(0,2) of s,
     hex substring(2,4) of s,
     hex substring(4,6) of s)/255
  elseif l_=11:
    (scantokens substring(0,3) of s,
     scantokens substring(4,7) of s,
     scantokens substring(8,11)of s)/255
  else: (??,??,??)
  fi
enddef;
vardef pencolor expr s_ = 
  color thepencolor;
  if color s_: thepencolor:=s_; 
  elseif string s_: thepencolor:=color_from_string s_; 
  fi
enddef;
pencolor black;
def savepencolor = 
  color _fgnd; _fgnd=thepencolor;
  save thepencolor; pencolor _fgnd;
enddef;
def setpencolor =
  drawoptions if thepencolor=black: () else: (withcolor thepencolor) fi;
enddef;
picture withdashdots;
withdashdots = dashpattern(on 3 off 3 on 0 off 3);
vardef typeof expr x =
  if numeric x: "numeric"
  elseif boolean x: "boolean"
  elseif pair x: "pair"
  elseif color x: "color"
  elseif transform x: "transform"
  elseif string x: "string"
  elseif path x: "path"
  elseif pen x: "pen"
  elseif picture x: "picture"
  else: "UNKNOWN"
  fi
enddef;
vardef templateof @# =
  save i_,p_,q_,s_,t_,skip_; i_=0; skip_=0;
  string s_,t_; s_=str @#;
  forever: t_:=substring(i_,incr i_) of s_;
    if skip_=0:
      if t_="[": skip_:=1; p_:=i_; 
      elseif isdigit t_: skip_:=2; p_:=i_-1;
      fi
    elseif skip_=1:
      if t_="]": skip_:=0; q_:=i_-1; 
        s_:= (substring(0,p_) of s_)&(substring(q_,infinity) of s_);
        i_:=i_-q_+p_;
      fi
    else:
      if not ((isdigit t_) or (t_=".")): 
        skip_:=0; q_:=i_-1;
        s_:= (substring(0,p_) of s_)&"[]"&(substring(q_,infinity) of s_);
        i_:=i_-q_+p_+1; 
      fi
    fi
    exitif t_="";
  endfor;
  s_
enddef;
def eatcolon(text t) =
  gobble begingroup
    save:; def : = endgroup t enddef;
enddef;
vardef compare_filedate(expr sfd,fd,v) =
  message("Package: `gamelan' "&v&" <"&fd&">");
  if sfd<>fd: 
    message "Warning: You are using gamelan.sty from "&sfd;
    message "         together with gamelan.mp  from "&fd;
  fi
  message "";
enddef;
vardef isdigit primary d =
  ("0"<=d)and(d<="9") 
enddef;
vardef cspan(expr s)(text c) =
  0
  for i=1 upto length s:
    exitunless c substring (i-1,i) of s;
    + 1
  endfor
enddef;
vardef loptok suffix s =
  save t_,j_,k_; string t_;
  j_ = cspan(s," ">=);
  k_ = j_ for i=j_ upto length s: exitif theseparator=substring(i,i+1) of s;
          +1 endfor;
  t_ = substring (j_,k_) of s;
  s := substring (k_+1,infinity) of s;
  t_
enddef;
