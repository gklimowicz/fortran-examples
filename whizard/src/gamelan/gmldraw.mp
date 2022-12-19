% gmldraw.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
let plaindraw=draw;
let plainfill=fill;
let plainfilldraw=filldraw;
let plainlabel=label;
let plaindotlabel=dotlabel;
def enable_enhanced_drawing =
  save draw,fill,filldraw,phantom,label,dotlabel,enhanced_drawing;
  def draw     = Gdraw_(false,true ) enddef;
  def fill     = Gdraw_(true, false) enddef;
  def filldraw = Gdraw_(true ,true ) enddef;
  def phantom  = Gdraw_(false,false) enddef;
  def label    = Glabel_ enddef;
  def dotlabel = Gdotlabel_ enddef;
  enhanced_drawing=1;
enddef;
def disable_enhanced_drawing =
  save draw,fill,filldraw,phantom,label,dotlabel,enhanced_drawing;
  let draw=plaindraw;
  let fill=plainfill;
  let filldraw=plainfilldraw;
  let label=plainlabel;
  let dotlabel=plaindotlabel;
enddef;
enable_enhanced_drawing;
vardef image(text t) =
  save currentpicture,Gfixed_;
  picture currentpicture;
  currentpicture := nullpicture;
  t;
  currentpicture
enddef;
def Gdraw_(expr fl,dr) =
  begingroup
  initdata d_; d_0_data0=nullpicture;
  d_0_fill=fl; d_0_draw=dr; d_0_cycl=false; d_0_immd=true;
  d_0_lnjn=linejoin; d_0_lncp=linecap; d_0_lnml=miterlimit;
  Gdef_dsuf_; Gredef_plot_;
  save labn_,ltyp_,ldsc_,lloc_; labn_=0; pair lloc_;
  save leg_; picture leg_;
  Gdrawbody_
enddef;
def Gdrawbody_ expr p =
  interim linejoin:=d_0_lnjn; interim linecap:=d_0_lncp;
  interim miterlimit:=d_0_lnml;
  addto d_0_data0
    if picture p: 
      if cycle pathpart p: cyclic fi
      also p
    elseif path p:
      if cycle p: cyclic fi
      if d_0_fill: contour else: doublepath fi  p  
      if d_0_draw: withpen currentpen fi  _op_
    else: also nullpicture
    fi
  Gdrawend_
enddef;
def Gdrawend_ text t_ =
  t_;
  if (not d_0_immd) and (unknown d_0_link): d_0_link:=defaultlink; fi
  Gupdate_legend_;
  drawdata d_0;
  endgroup
enddef;
vardef Gupdate_legend_ =
  if known leg_:
    save p,dt; path p; picture dt; dt=nullpicture;
    if d_0_draw or d_0_fill:
      p = if d_0_cycl: unitsquare shifted (-.5,-.5) 
          else: ((-.5,0)--(.5,0))
          fi xscaled legendbarwidth yscaled legendbarheight;
      Gfindopts_(d_0_data0);
      Gaddl_(dt)(d_0_fill) p Gopt_;
      if (d_0_cycl):
        if known d_0_outl: Gfindopts_(d_0_outl);
          addto dt doublepath p Gopt_;
        fi
      else:
        if known d_0_endp0:
          addto dt also d_0_endp0 rotated 180 shifted point 0 of p Gopt_;
          addto dt also d_0_endp[infinity]    shifted point 1 of p Gopt_;
        fi
      fi
    fi
    if known d_0_symb: addto dt also d_0_symb; fi
    legendentry(dt,leg_);    
  fi
enddef;
newinternal legendbarwidth; legendbarwidth:=10mm;
newinternal legendbarheight; legendbarheight:=3mm;
def Gaddl_(suffix pic)(expr fl) expr p =
  addto pic
    if picture p: also
    elseif fl: contour
    else: doublepath
    fi p
enddef;
def Gdef_dsuf_ =
  save cyclic,piecewise,pointwise,pc_;  pc_=0; 
  let cyclic=Gdcycl_; let piecewise=Gdpc_; let pointwise=Gdpt_;
enddef;
def Gdcycl_ = hide(d_0_cycl:=true) enddef;
def Gdpc_ = hide(pc_:=1) enddef;
def Gdpt_ = hide(pc_:=2) enddef;
vardef from(text chains_) =
  if known Gfixed_:
    d_0_immd:=false;
    save i_; i_=1; d_0_data1=d_0_symp=nullpicture;
    Gscan_init_(chains_);
    Gadd_
      if pc_=0: Gscan_ elseif pc_=1: Gscan_pc_ else: Gscan_pt_ fi
      (Gwsymb_($!) Gpint_ z 
       if $~=2: -- 
       elseif $~=1: ; Gadd_
       else:  ; d_0_data[incr i_]=nullpicture;  Gadd_
       fi) 0;
  else:
    errhelp("Data drawing commands are possible only within a graph"
           &"environment, and outside of `image' expressions."
           &"Sorry, I have to discard this one.");
    errmessage("`from' not allowed here");
    d_0_draw:=d_0_fill:=false;
  fi
  ttempl_ if d_0_cycl or d_0_fill: --cycle fi
enddef;
vardef Gwsymb_(expr f) primary z =
  if f: addto d_0_symp doublepath z; fi 
  z
enddef;
def Gadd_ expr p_ =
  if path p_: addto d_0_data[i_] doublepath p_; fi
enddef;
def Gredef_plot_ =
  let plot_=plot; let hist_=hist; let table_=table;
  save tmp_,plot,hist,table,tables;
  let plot=Gdplot_; let hist=Gdhist_; let table=Gdtab_;
enddef;
def Gdplot_(text zz_) =
  hide(plot_(tmp_) zz_) from(tmp_)
enddef;
def Gdhist_(text zz_) =
  hide(hist_(tmp_) zz_) from(tmp_)
enddef;
def Gdtab_ =
  hide(let plot=Gtplot_; let hist=Gthist_)
enddef;
def Gtplot_(text xy_) =
  hide(table_ plot_(tmp_)(xy_)) from(tmp_)
enddef;
def Gthist_(text xy_) =
  hide(table_ hist_(tmp_)(xy_)) from(tmp_)
enddef;
def Glabel_ = phantom 0 withlabel enddef;
def Gdotlabel_ = phantom 0 withdotlabel enddef;
def withlabel =
  hide(ltyp_:=0) hide(ldsc_:=0) gobble Gwithlab_
enddef;
def withdotlabel =
  hide(ltyp_:=0) hide(ldsc_:=dotscale) gobble Gwithlab_
enddef;
vardef Gwithlab_ @#(expr t, loc) =
  save p; picture p; p=nullpicture;
  addto p also thelabel @#(t,if ltyp_>0: loc else: origin fi);
  if ltyp_=4: outerlabel @#(p);
  else:
    d_0_lab[incr labn_]=p;
    d_0_lloc[labn_]=if ltyp_>0: lloc_;  pair lloc_ else: loc fi;
    d_0_ldsc[labn_]=ldsc_;
    d_0_ltyp[labn_]=ltyp_;
  fi
enddef;
vardef on graph at primary loc = 
  if known Gfixed_:
    ltyp_:=1; lloc_=loc;
  else:
    errhelp("You can't say `on graph' if there is no graph. "
           &"I'll put the label at the coordinate origin.");
    errmessage("`on graph at' outside graph or inside `image' expression");
  fi
  origin   
enddef;
vardef on curve at primary loc = 
  ltyp_:= if pair loc: 3; loc  else: 2; (loc,0) fi
    = lloc_; origin
enddef;
vardef out = ltyp_:=4; origin enddef;
vardef margin @#(text u_) =
  (width*(1-labxf @#), height*(1-labyf @#))
  for u=u_: + u*(2labxf @#-1, 2labyf @#-1) endfor
enddef;
def withlinejoin primary j = hide(d_0_lnjn := j) enddef;
def withlinecap  primary c = hide(d_0_lncp := c) enddef;
def withmiterlimit primary c = hide(d_0_lnml := c) enddef;
def linked primary l =
  hide(d_0_link := if numeric l: (l,l) else: l fi)
enddef;
newinternal none,straight,smoothly,tight;
none:=0; straight:=1; smoothly:=2; tight:=3;
pair defaultlink; defaultlink:=(straight,straight);
def withsymbol primary s =
  hide(d_0_symb:=thesymbol(s))
enddef;
def withlegend primary s =
  hide(leg_ :=thesymbol(s))
enddef;
def outlined text t_ =
  hide(d_0_outl:=nullpicture)
  hide(addto d_0_outl doublepath origin withpen currentpen t_)
enddef;
def withshadow text t_ =
  hide(d_0_shad:=nullpicture)
  hide(Gaddl_(d_0_shad)(d_0_cycl) ttempl_ if d_0_cycl: --cycle fi t_)
enddef;
def withends(text t_) =
  for p_=t_: 
    hide(d_0_endp if known d_0_endp0: [infinity]:= else: [0]= fi
         if known p_: thesymbol(p_) else: nullpicture fi)
  endfor  
  hide(if unknown d_0_endp[infinity]: d_0_endp[infinity]=d_0_endp0 fi)
enddef;
def withend(expr p_) = withends(??,p_) enddef;
def witharrows = withends(thearrow scaled arrowscale) enddef;
def witharrow = withends(??,thearrow scaled arrowscale) enddef;
def withticks = withends(thetick scaled (ticklength*tickscale)) enddef;
def withtick = withends(??,thetick scaled (ticklength*tickscale)) enddef;
newinternal arrowscale; arrowscale:=1;
def thearrow = 
  image(plainfilldraw arrowhead (100,0)--(0,0)) rotated 180 
enddef;
newinternal ticklength; ticklength:=7bp;
newinternal tickscale;  tickscale:=4/7;
def thetick = (down--up) enddef;
def theotick = (down--origin) enddef;
def theitick = (origin--up) enddef;
