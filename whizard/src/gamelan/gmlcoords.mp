% gmlcoords.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
Gzz_:=zero/2;
vardef Gpint_ primary p = 
  (if abs Gtype_.x=linear: Gcint_ fi xpart p, 
   if abs Gtype_.y=linear: Gcint_ fi ypart p)
enddef;
def Gcint_ expr x =
  if x<>Mabs x: - fi (.5x-Gzz_)
enddef;
vardef Gpext_ primary p =
  (if abs Gtype_.x=linear: Gcext_ fi xpart p, 
   if abs Gtype_.y=linear: Gcext_ fi ypart p)
enddef;
def Gcext_ expr x =
  if x>=0: 2(x+Gzz_)
  else:    2(-x+Gzz_) + epsilon
  fi
enddef;
def begingraph(expr w, h) =
  begingroup
  save width, height; width=w; height=h;  
  save theframe; path theframe; theframe=unitsquare xscaled w yscaled h;
  if known hatchdim: hatchdim:=2(w++h); fi
  savecurrentpen;
  save Gneedgr_, Gneedfr_, Gneedol_, Gb_, Glb_;  
  boolean Gneedgr_, Gneedfr_, Gneedol_;  Gneedgr_=Gneedfr_=Gneedol_=true;
  pair Gb_.low, Gb_.high, Glb_.low, Glb_.high;
  forsuffixes $$=currentpicture, theouterlabels, thelegend:
    save $$; picture $$; $$=nullpicture;
  endfor
  enable_enhanced_drawing;
  initdata Gdata_;
  setup(linear,linear);
enddef;
def setup(expr xtype,ytype) =
  if Gdata_.count>0: freeze; fi
  save Gtype_; (Gtype_.x,Gtype_.y)=(xtype,ytype);
  save Gfixed_, Gr_, Glsc_, Glmax_, Gmr_;  pair Gr_.low, Gr_.high;
  pair Gmr_.low, Gmr_.high;
  boolean Gfixed_; Gfixed_=false;
  save convert;
enddef;
newinternal linear, log; linear:=1; log:=2;
vardef graphrange text t_ = 
  Grdrect_(Gr_) t_;
  if (known Gr_.high) and (known Gr_.low): 
    Gfixed_:=true;
    Gfindconv_(Gfixed_,Gr_,Glsc_,Glmax_); 
  fi
enddef; 
vardef labelbounds text t_ =
  Grdrect_(Glb_) t_;
enddef;
def framelabels =
  labelbounds llcorner theframe, urcorner theframe;
enddef;
vardef graphbounds text t_ = 
  Grdrect_(Gb_) t_;
enddef;
def framebounds = 
  graphbounds llcorner theframe, urcorner theframe; 
enddef;
vardef Grdrect_(suffix rr) text t_ =
  save i_,r_; i_=0;
  for x_=
      for p_=t_: if pair p_: xpart p_,ypart fi p_, endfor:
    r_[incr i_] = x_; exitif i_=4;
  endfor
  i_:=0;
  forsuffixes $$=low, high:  pair rr$$;
    rr$$ = (r_[incr i_], r_[incr i_]);
  endfor
enddef;
def initdata suffix d_ =
  save d_; d_.count:=0;
  boolean d_[]_immd;
  picture d_[]_data[];  pair d_[]_link;
  boolean d_[]_draw, d_[]_fill, d_[]_cycl;
  numeric d_[]_lnjn, d_[]_lncp, d_[]_lnml;
  picture d_[]_shad, d_[]_outl, d_[]_endp[], d_[]_symb, d_[]_symp;
  picture d_[]_lab[];  pair d_[]_lloc[];  
  numeric d_[]_ltyp[], d_[]_ldsc;
enddef;
vardef Gcleardata_ @# =
  forsuffixes p = @#_shad, @#_outl, @#_endp[0], @#_endp[infinity],
      @#_symb, @#_symp:
    p:=nullpicture;
  endfor
  save i; i=0;
  forever: exitif unknown @#_data[incr i];
    @#_data[i]:=nullpicture;
  endfor
enddef;
def drawdata suffix d_ =
  if if known Gfixed_: Gfixed_ else: true fi:
    addto currentpicture also Gimage_ d_;
  else:
    Gstoredata_ d_;
  fi
enddef;
vardef Gstoredata_ suffix d_ =
  save i_,j_;  i_ := incr Gdata_.count;
  forsuffixes @@=_immd, _link, _draw, _fill, _cycl, _lnjn, _lncp, _lnml,
    _shad, _outl, _symb, _symp, _endp[0], _endp[infinity]:
    if known d_@@: Gdata_[i_]@@=d_@@; fi
  endfor
  if d_._immd: Gdata_[i_]_data0 = d_._data0; 
  else: j_:=-1;
    forever: exitif unknown d_._data[incr j_];
      Gdata_[i_]_data[j_] = d_._data[j_];
    endfor
  fi
  j_:=0;
  forever: exitif unknown d_._lab[incr j_];
    Gdata_[i_]_lab[j_] = d_._lab[j_];
    Gdata_[i_]_lloc[j_]= d_._lloc[j_];
    Gdata_[i_]_ldsc[j_]= d_._ldsc[j_];
    Gdata_[i_]_ltyp[j_]= d_._ltyp[j_];
  endfor
enddef;
def freeze =
  if not Gfixed_: 
    Gfindconv_(Gfixed_,Gr_,Glsc_,Glmax_); 
    Gfixed_:=true;
    for i=1 upto Gdata_.count:
      addto currentpicture also Gimage_ Gdata_[i];
      Gcleardata_ Gdata_[i];
    endfor
    initdata Gdata_;
  fi
enddef;
def endgraph =
  freeze;
  if Gneedgr_: defaultgrid(); fi
  if Gneedfr_: frame; fi
  if Gneedol_: outerlabels(Gbbox_(Glb_,currentpicture)); fi
  interim truecorners:=1;
  setbounds currentpicture to Gbbox_(Gb_,currentpicture);
  currentpicture
  endgroup
enddef;
vardef Gbbox_(suffix b,pic) =
  interim bboxmargin:=0;
  Gfixb_(b,pic,low )(llcorner);
  Gfixb_(b,pic,high)(urcorner);
  bbox(b.low--b.high)
enddef;
def Gfixb_(suffix b,pic,hl)(text xxcorner) =
  (if unknown xpart b.hl: xpart b.hl else: ?? fi,
   if unknown ypart b.hl: ypart b.hl else: ?? fi)
  = xxcorner pic;
enddef;
vardef outerlabel @#(expr pic) =
  addto theouterlabels also
    image(currentpicture:=thelabel.@#(pic,origin);
          setbounds currentpicture to origin..cycle)
      shifted ((1,1)-(labxf @#, labyf @#));
enddef;
vardef outerlabels(expr p) =
  Gneedol_:=false;
  if length theouterlabels>0:
    save t; transform t; t=thetransform(p);
    Gadd_dummy_(theouterlabels);
    for q within theouterlabels:
      addto currentpicture also q
        shifted ((-llcorner q) + ((llcorner q) transformed t));
    endfor
    theouterlabels:=nullpicture;
  fi
enddef;
def legendentry(expr dt,dd) =
  interim labeloffset:=legendoffset;
  save p; picture p; p=nullpicture;
  if known dt: addto p also thelabel.lft(dt, origin); fi
  if known dd: addto p also thelabel.rt (dd, origin); fi
  addto thelegend also p
    shifted (0, (ypart llcorner thelegend)
               -(ypart llcorner p)-legendlineskip);
enddef;
newinternal legendoffset; legendoffset:=labeloffset;
newinternal legendlineskip; legendlineskip:=5mm;
def legend =
  addto currentpicture also thelegend
    hide(thelegend:=nullpicture)
enddef;
vardef Gfindconv_(expr f)(suffix r,s,m) =
  if not f: 
    save b; picture b; b=nullpicture;
    Gfindbbox_ b;
    Gfindrange_(b) r;
  fi
  Gfindscale_(urcorner theframe - llcorner theframe, r) s;
  Gfindmax_(s) m;
  Gfindtransf_(r) t;
  Gdefconv_(t);
enddef;
def Gfindbbox_ suffix b =
  save j;
  interim bboxmargin:=0;
  for i=1 upto Gdata_.count: j:=0;
    forever: 
      exitif unknown Gdata_[i]_data[incr j];
      addto b doublepath llcorner Gdata_[i]_data[j] .. 
        urcorner Gdata_[i]_data[j];
    endfor
    b := image(addto currentpicture doublepath llcorner b .. urcorner b);
  endfor
  if length b=0:  message "Warning: graph may be empty.";
    setbounds b to  Gpint_(Gorig_.x, Gorig_.y) .. Gpint_(ten,ten) -- cycle;
  else:
    if xpart llcorner b = xpart lrcorner b:
      message "Warning: graph degenerate in horizontal direction.";
      addto b also Gdummyp_(x)(Gpext_ llcorner b);
    fi
    if ypart llcorner b = ypart ulcorner b:
      message "Warning: graph degenerate in vertical direction.";
      addto b also Gdummyp_(y)(Gpext_ llcorner b);
    fi
  fi
enddef;  
def Gorig_ suffix xy =
  if abs Gtype_.xy=linear: zero else: one fi
enddef;
def Gdummyp_(suffix xy)(expr p) =
  image(setbounds currentpicture to 
        (Gshift_(xy)(p, #1, #10) .. Gshift_(xy)(p, #-1, #1/10) -- cycle))
enddef;
def Gshift_(suffix xy)(expr p, linshift, logshift) =
  Gshiftt_(xy)(p)
    (if abs Gtype_.xy=linear: plus linshift else: times logshift fi)
enddef;
def Gshiftt_(suffix xy)(expr p)(text mod) =
  Gpint_
  if str xy="x": (xpart p mod, ypart p)
  else:          (xpart p, ypart p mod)
  fi
enddef;
def Gfindrange_(suffix b) suffix r =
  save bl,bh; pair bl,bh;  
  bl = Gpext_ llcorner b;  bh = Gpext_ urcorner b;
  Gfix_(bl,bh) r.low;
  Gfix_(bh,bl) r.high;
enddef;
def Gfix_(expr ba,bb) suffix r =
  Gfixx_(xpart)(x)(ba,bb) r;
  Gfixx_(ypart)(y)(ba,bb) r;
enddef;
def Gfixx_(text xypart)(suffix xy)(expr ba,bb) suffix r =
  if unknown xypart r:
    xypart r =
      if abs Gtype_.xy=linear:
        #(1+marginfraction) times xypart ba 
        plus #(-marginfraction) times xypart bb
      else:
        (1+marginfraction)*(xypart ba) 
        + (-marginfraction)*xypart bb
      fi;
  fi
enddef;
newinternal marginfraction; marginfraction:=.07;
def Gfindscale_(expr d)(suffix r) suffix s =
  s.x := Mabs(#(xpart d) over max(xpart r.low, xpart r.high));
  s.y := Mabs(#(ypart d) over max(ypart r.low, ypart r.high));
enddef;
def Gfindmax_(suffix s) suffix m =
  m.x := Mlog 2infinity - s.x;
  m.y := Mlog 2infinity - s.y;
enddef;
def Gfindtransf_(suffix r) suffix t =
  transform t;  xypart t = yxpart t = 0;
  Gconv_ r.low  transformed t = llcorner theframe;
  Gconv_ r.high transformed t = urcorner theframe;
  if Gtype_.x<0: 
    t:=t reflectedabout (.5[llcorner theframe, lrcorner theframe], 
                         .5[ulcorner theframe, urcorner theframe]);
  fi
  if Gtype_.y<0: 
    t:=t reflectedabout (.5[llcorner theframe, ulcorner theframe], 
                         .5[lrcorner theframe, urcorner theframe]);
  fi
enddef;
vardef Gconv_ primary p =
  (if abs Gtype_.x=linear: Gtolin_.x fi xpart p,
   if abs Gtype_.y=linear: Gtolin_.y fi ypart p)
enddef;
vardef Gtolin_ @# primary x =
  if known x:
    if Mabs x < Glmax_ @#:  Mexp (x + Glsc_ @#)
    elseif x=Mabs x: 2infinity
    else:           -2infinity
    fi
  else: ??
  fi
enddef;
def Gdefconv_(expr t) =
  vardef convert primary p_ = Gconv_ p_ transformed t enddef;
enddef;
vardef Gimage_ @# =
  save Gopt_,cp,tmp,pp,tf; picture cp,tmp; cp=tmp=nullpicture; path pp;
  if not @#_immd: transform tf; tf=thetransform(pathpart @#_data0); fi
  interim linejoin:=@#_lnjn; interim linecap:=@#_lncp;
  interim miterlimit:=@#_lnml;
  Gdraw_data_(@#_data, @#_immd, @#_fill, @#_link, @#_cycl) tmp;
  if @#_draw or @#_fill:
    Gdraw_shadows_(@#_fill, @#_shad, tmp) cp;
    addto cp also tmp;
    if @#_cycl: Gdraw_outlines_(@#_outl, tmp) cp;
    else:       Gdraw_ends_(@#_endp, colorpart @#_data0, tmp) cp;
    fi
  fi
  Gdraw_symbols_(@#_immd, @#_symp, @#_symb) cp;
  cp := cp if known tf: transformed tf fi;
  if known Gfixed_: clip cp to theframe; fi
  Gdraw_labels_(@#_lab, @#_lloc, @#_ldsc, @#_ltyp, @#_immd, tf, tmp) cp;
  cp
enddef;
vardef Gdraw_data_(suffix data)(expr immd, fl, link, cycl) suffix pic =
  if immd: addto pic also data[0];
  else:
    Gdefmodes_(fl,link,cycl)(pic);
    Gfindopts_(data[0]);
    save j,k,p; path p; j:=0;
    forever: hide(k:=length data[incr j]) exitif k=0;
      Gaddc_  
        for d within data[j]: 
          hide(p:=pathpart d)
          convert Gpext_ point 0 of p
            for i=1 upto length p:
              Gccc_ convert Gpext_ point i of p
            endfor
            if (decr k>0): Gcc_ else: Gcyc_ Gopt_; fi
        endfor
    endfor
  fi
  Gadd_dummy_(pic);
enddef;
picture Gdummy_; Gdummy_=image(setbounds currentpicture to origin..cycle);
def Gadd_dummy_(suffix pic) =
  if length pic>0: addto pic also Gdummy_; fi
enddef;
def Gdefmodes_(expr fl,link,cycl)(suffix pic) =
  Gdef_(Gaddc_ expr p)
    if fl: (addto pic if picture p: also else: contour fi p _op_)
    else:  (addto pic if picture p: also else: doublepath fi p _op_)
    fi;
  Gcdef_(Gccc_)(xpart link);
  Gcdef_(Gcc_)(ypart link);
  Gdef_(Gcyc_) if cycl and (ypart link<>none): (Gcc_ cycle) else: () fi;
enddef;
def Gcdef_(text m)(expr ln) =
  Gdef_(m) if ln=straight: (--)
           elseif ln=smoothly: (..)
           elseif ln=tight: (...)
           else: (Gcyc_ Gopt_; Gaddc_)
           fi
enddef;
vardef Gfindopts_(expr pic) =
  save gc,gp,gd,gh; color gc; pen gp; picture gd,gh;
  gh = hatchpart pic; 
  if known gh: gc=colorpart gh;
  else: gc = colorpart pic; gp = penpart pic; gd = dashpart pic; 
  fi
  Gdefopts_(gc,gp,gd,gh);
enddef;
def Gdefopts_(expr gc,gp,gd,gh) =
  Gdef_(Gopt_)
    if known gh: (hatched gh) 
    else: (withpen gp dashed gd withcolor gc) 
    fi;
enddef;
def Gdef_(text m)(text b) =
  def m= b enddef
enddef;
def Gdraw_shadows_(expr fl, shad, data) suffix pic =
  if known shad:  
    save shtf; transform shtf; shtf=thetransform(pathpart shad);
    Gfindopts_(shad);
    Gdef_(Gaddc_ expr p)
      if fl: (addto pic if picture p: also else: contour fi p _op_)
      else:  (addto pic if picture p: also else: doublepath fi p _op_)
      fi;
    for q within data:
      if (stroked q) or (filled q) or (clipped q):
        Gaddc_ pathpart q transformed shtf Gopt_;
      fi
    endfor
  fi
enddef;
def Gdraw_outlines_(expr outl, data) suffix pic =
  if known outl:
    Gfindopts_(outl);
    for q within data:
      if (stroked q) or (filled q) or (clipped q):
        addto pic doublepath pathpart q Gopt_;
      fi
    endfor
  fi
enddef;
vardef Gdraw_ends_(suffix endp)(expr col, data) suffix pic =
  save p; path p;
  for q within data:
    if (stroked q) or (filled q) or (clipped q):
      p:=pathpart q;
      for j=0,infinity:
        if known endp[j]:
          addto pic also endp[j]
            rotated angle if j=0: - fi direction j of p
            shifted point j of p withcolor col;
        fi
      endfor
    fi
  endfor
enddef;
def Gdraw_symbols_(expr immd, symp, symb) suffix pic =
  if (known symb) and (known symp):
    for q within symp:
      addto pic also symb
        shifted if not immd: convert Gpext_ fi (llcorner q);
    endfor
  fi
enddef;
vardef Gdraw_labels_(suffix lab, loc, dsc, typ)
    (expr immd, traf, data) suffix pic =
  if known lab[1]:
    save i,p; i=0; path p;  p=Gconnected_path_(data);
    interim linecap:=rounded; 
    forever: exitif unknown lab[incr i];
      Gfindlloc_(immd, typ[i], p, traf) loc[i];
      if dsc[i]>0: addto pic doublepath loc[i] withpenscale dsc[i]; fi
      addto pic also lab[i] shifted loc[i];
    endfor
  fi
enddef;
def Gfindlloc_(expr immd, typ, p, traf) suffix loc =
  if typ=1:
    loc := convert loc;
    (if unknown xpart loc: xpart loc else: ?? fi, 
     if unknown ypart loc: ypart loc else: ?? fi)
    = llcorner theframe;
  elseif typ=2:  loc := point (xpart loc) of p;
  elseif typ=3:
    if not immd: loc := convert loc; fi
    if known loc: relax
    elseif known xpart loc:
      loc:=((xpart loc,-2infinity)--(xpart loc,2infinity)) intersectionpoint p;
    elseif known ypart loc:
      loc:=((-2infinity,ypart loc)--(2infinity,ypart loc)) intersectionpoint p;
    else: loc:=point 0 of p;
    fi
  fi
  if known traf: loc := loc transformed traf; fi
enddef;
vardef Gconnected_path_(expr data) =
  if length data=0: origin
  else:
    save c; boolean c; c=false;
    for d within data:
      if (stroked d) or (filled d) or (clipped d):
        if c: -- else: hide(c:=true) fi (pathpart d)
      fi
    endfor
  fi
enddef;
