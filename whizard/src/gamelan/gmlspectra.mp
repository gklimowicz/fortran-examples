% gmlspectra.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
if unknown hexacross: input gmlshapes fi
vardef create_array @#(text typ_)(text val_) =
  typ_ @#[]; numeric @#?; @#?=-1;
  augment_array @#(val_);
enddef;
vardef augment_array @#(text val_) =
  if unknown @#?:
    Gundef_array_err_(@#);
  else:
    for v_=val_: @#[incr @#?]:=v_; endfor
  fi
enddef;  
vardef show_array @# =
  if unknown @#?:
    Gundef_array_err_(@#);
  elseif @#?<0:
    message "Array "& str @# &" is empty.";
  else:
    message "Contents of array "& str @# &":";
    for i=0 upto @#?: show @#[i]; endfor
  fi
enddef;  
vardef interpolate_array @#(expr x_) =
  if unknown @#?:  Gundef_array_err_(@#); 0
  elseif @#?<0:  Gempty_array_err_(@#); 0
  else:
    save i_,u_; u_=x_;    
    if (u_<0): u_:=0; elseif (u_>1): u_:=1; fi
    if (u_=1): @#[@#?]
    else:
      u_:=u_*@#?;  i_=floor u_; u_:=u_-i_; u_[@#[i_],@#[i_+1]]
    fi
  fi  
enddef;
def Gundef_array_err_(suffix @#) =
  errhelp("The array "& str @# &" has not been defined.  Use `create_array' "&
          "first.                   "&
          "(Maybe you misspelled the name.)");
  show str @#;
  errmessage("Unknown array");
enddef;
def Gempty_array_err_(suffix @#) =
  errmessage("The array "& str @# &" is empty, so I can't interpolate. "&
             "Use `augment_array' to fill it.");
  show str @#;
  errmessage("Empty array"); 0
enddef;
def spectrum(expr ncol)(text cols) =
  image(
    save ca,x,dx,u; create_array ca(color)(cols);
    dx = 1/ncol;
    for i=0 upto ncol-1: 
      x:=i/ncol;  u:=i/(ncol-1);
      fill unitsquare xscaled dx shifted (x,0)
          withcolor interpolate_array ca(u);
    endfor
  )
enddef;
def rayspectrum(expr ncol, angle)(text cols) =
  image(
    save ca,da; create_array ca(color)(cols);
    da = angle/ncol;
    for i=0 upto ncol-1: 
      fill wedge (i*da, (i+1)*da) 
        withcolor interpolate_array ca(i/(ncol-1));
    endfor
  )
enddef;    
def spotspectrum(expr ncol, angle)(text cols) =
  image(
    save ca; create_array ca(color)(cols);
    for i=0 upto ncol-1: 
      fill doublearc(i/(i+1)) (0,angle) scaled (i+1)/ncol
        withcolor interpolate_array ca(i/(ncol-1));
    endfor
  )
enddef;    
