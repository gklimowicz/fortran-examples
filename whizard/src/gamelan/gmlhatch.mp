% gmlhatch.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
def fill expr p =
  addto currentpicture
  if picture p:
    also p
  else:
    contour p
  fi
  _op_
enddef;
tertiarydef p_ hatched q_ = 
  if known q_:
    begingroup
    save q; picture q; q = q_;
    if hatchable p_: clip q to p_; q
    elseif picture p_:
      save pic,pa; picture pic; path pa; 
      pic=p_; pa=pathpart pic;
      if hatchable pa: clip q to pa; addto pic also q;  fi 
      pic
    else: p_
    fi 
    endgroup
  else: p_
  fi
enddef;
def withbackground = hatched enddef;
def hatchpattern(text t_) =
  image(interim linecap:=butt;
        addto currentpicture doublepath (0,-hatchdim)..(0,hatchdim) 
        withpen (pencircle scaled hatchdim)
        dashed dashpattern(t_))
enddef;
def withstripes = hatchpattern(on 3 off 3) enddef;
def withlines = hatchpattern(on .5 off 5) enddef;
newinternal hatchdim; hatchdim := 20cm;
def hatchable primary p =
  cycle p and (llcorner p<>ulcorner p) and (llcorner p<>lrcorner p)
enddef;
def hatchpart expr p =
  image(
    if clipped p: 
      for q within p: addto currentpicture also q; endfor
    else: picture currentpicture;
    fi )
enddef;
