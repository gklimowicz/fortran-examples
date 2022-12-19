% gmlbox.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
def fbox(expr p_)(text tf_)(text tl_) =
  image(save q_; picture q_; q_=thesymbol(p_);
        fill bbox q_ tf_; draw q_ tl_)
enddef;
def ufbox(expr p_)(text tl_) =
  image(save q_; picture q_; q_=thesymbol(p_);
        unfill bbox q_; draw q_ tl_)
enddef;
vardef fillbox @#(expr p_, sc_)(text tf_)(text tl_) =
  image(fill unitsquare xscaled xpart sc_ yscaled ypart sc_ tf_; 
        _drilabel_ @#(p_, (1-labxf @#,1-labyf @#) 
          xscaled xpart sc_ yscaled ypart sc_) tl_)
enddef;
vardef unfillbox @#(expr p_, sc_)(text tl_) =
  image(unfill unitsquare xscaled xpart sc_ yscaled ypart sc_; 
        _drilabel_ @#(p_, (1-labxf @#,1-labyf @#) 
          xscaled xpart sc_ yscaled ypart sc_) tl_)
enddef;
def dbox(expr p_)(text td_)(text tl_) =
  image(save q_; picture q_; q_=thesymbol(p_);
        draw bbox q_ td_; draw q_ tl_)
enddef;
vardef drawbox @#(expr p_, sc_)(text td_)(text tl_) =
  image(draw unitsquare xscaled xpart sc_ yscaled ypart sc_  td_; 
        _drilabel_ @#(p_, (1-labxf @#,1-labyf @#)
          xscaled xpart sc_ yscaled ypart sc_) tl_)
enddef;
def fdbox(expr p_)(text tf_)(text td_)(text tl_) =
  image(save q_; picture q_; q_=thesymbol(p_);
        fill bbox q_ tf_; draw bbox q_ td_; draw q_ tl_)
enddef;
vardef filldrawbox @#(expr p_, sc_)(text tf_)(text td_)(text tl_) =
  image(fill unitsquare xscaled xpart sc_ yscaled ypart sc_ tf_;
        draw unitsquare xscaled xpart sc_ yscaled ypart sc_ td_; 
        _drilabel_ @#(p_, (1-labxf @#,1-labyf @#)
          xscaled xpart sc_ yscaled ypart sc_) tl_)
enddef;
def ufdbox(expr p_)(text td_)(text tl_) =
  image(save q_; picture q_; q_=thesymbol(p_);
        unfill bbox q_; draw bbox q_  td_; draw q_ tl_)
enddef;
vardef unfilldrawbox @#(expr p_, sc_)(text td_)(text tl_) =
  image(unfill unitsquare xscaled xpart sc_ yscaled ypart sc_;
        draw unitsquare xscaled xpart sc_ yscaled ypart sc_ td_; 
        _drilabel_ @#(p_, (1-labxf @#,1-labyf @#)
          xscaled xpart sc_ yscaled ypart sc_) tl_)
enddef;
vardef _drilabel_ @#(expr s_,z_) text t_ =
  save p_; picture p_; p_=thesymbol(s_)
  addto currentpicture also p_ shifted (z_ +
    labeloffset*(-laboff @#) -
       ((1-labxf @#)*lrcorner p_ + (1-labyf @#)*ulcorner p_
         - (1-labxf @#-labyf @#)*llcorner p_)) _op_ t_;
enddef;
