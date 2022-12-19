% gmlgrid.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
newinternal ticknumber; ticknumber := 4;
newinternal tickfactor; tickfactor := 1.1;
newinternal tickminlog; tickminlog := 3.0;
def withlength primary x = hide(Gtlen_:=x) enddef;
def withscale primary x = hide(Gtscl_:=x) enddef;
def withnumber primary x = hide(Gtnum_:=x) enddef;
def withfactor primary x = hide(Gtfac_:=x) enddef;
def withminlog primary x = hide(Gtmlg_:=x) enddef;
def Greset_ticks_ =
  numeric Gtlen_,Gtscl_,Gtnum_,Gtfac_,Gtmlg_;
enddef;
vardef notick @#(expr lab_,loc_) text w_ = 
  Gtick_ @#(??,lab_)(loc_)(false)(w_)
enddef;
vardef dot @#(expr lab_,loc_) text w_ = 
  Gtick_ @#(origin,lab_)(loc_)(false)(w_)
enddef;
vardef tick @#(expr lab_,loc_) text w_ = 
  Gtick_ @#(thetick,lab_)(loc_)(false)(w_)
enddef;
vardef itick @#(expr lab_,loc_) text w_ = 
  Gtick_ @#(theitick,lab_)(loc_)(false)(w_)
enddef;
vardef otick @#(expr lab_,loc_) text w_ = 
  Gtick_ @#(theotick,lab_)(loc_)(false)(w_)
enddef;
vardef line @#(expr lab_,loc_) text w_ =
  Gtick_ @#(theitick,lab_)(loc_)(true)(w_)
enddef;
vardef Gtick_ @#(expr tick_,lab_,loc_,f_)(text w_) =
  save p_,q_;  path p_; picture q_;
  q_=image(plaindraw ttempl_ w_);
  if unknown Gtraf_:
    save Gtraf_,Gopt_; transform Gtraf_;
    Gfindopts_(q_); Gtraf_=thetransform(pathpart q_); 
  fi
  q_:=nullpicture;
  if (known loc_) and (abs(laboff @# dotprod (1,1))=1):
    if known tick_:
      p_ = tick_ scaled
            if f_: abs(laboff @# dotprod (width,height))
            elseif known Gtlen_:  (Gtlen_ if known Gtsc_: *Gtsc_ fi)
            else: ticklength
            fi
           rotated 90 zscaled laboff @#;
      addto q_ doublepath p_ transformed Gtraf_ withpen currentpen _op_ Gopt_;
    fi
    if known lab_:
      addto q_ also 
        thelabel.@#
          (if string lab_: Mformat(lab_, loc_) else: lab_ fi, 
           if known p_: point 0 of p_ else: origin transformed Gtraf_ fi);
    fi
    if unknown Gtsc_: phantom 0 fi 
    withlabel.o(q_, Gloc_.@#(loc_));
  else:
    errhelp("Either the tick mark position is unknown "
           &"or you specified a wrong location     "
           &"suffix (for tick marks only "
           &"lft,rt,top,bot is allowed).                        "
           &"Don't be surprised if this error message "
           &"occurs many times, once for each tick mark.");
    errmessage("Illegal location for tick mark");
  fi
  if unknown Gtsc_: Greset_ticks_; fi
enddef;
vardef Gloc_.bot(expr loc_) = on graph at (loc_,??) enddef;
vardef Gloc_.lft(expr loc_) = on graph at (??,loc_) enddef;
vardef Gloc_.top(expr loc_) = on graph at (loc_,??) + (0,height) enddef;
vardef Gloc_.rt (expr loc_) = on graph at (??,loc_) + (width,0) enddef;
vardef grid @#(text markdef) text opts =
enddef;
vardef grid @#(text mtyp_) text w_ =
  save Gopt_, Gcount_, Gmks_, Gttyp_, Gtnum_, Gtform_, Gtloc_, Gtsc_;
  save suf_,Gtraf_,p_; string suf_; transform Gtraf_; picture p_;
  Gtlen_:=ticklength; Gtscale_:=tickscale; Gtnum_:=ticknumber;
  Gtfac_:=tickfactor; Gtmlg_:=tickminlog;
  p_ = image(plaindraw ttempl_ w_);
  Gfindopts_(p_); Gtraf_=thetransform(pathpart p_);
  string Gttyp_[], Gtform_[], Gtloc_[];
  Gcount_:=1; Ggdef_(mtyp_);
  forsuffixes $$=scantokens(Glocs_ @#):  
    suf_:=str $$; numeric Gmks_; 
    for c_=1 upto Gcount_-1:
      if false forsuffixes s_=scantokens Gtloc_[c_]: or (suf_=str s_) endfor:
        Gmarks:=if known Gtnum_[c_]: Gtnum_[c_]
                elseif known Gmks_: Gmks_*Gtfac_
                else: Gtnum_
                fi;
        Gmks_:=0;
        phantom 0
          for loc_= auto if labxf$$=.5: x else: y fi:
            hide(Gmks_:=Gmks_+1)  hide(Gtsc_:=Gtsc_[c_])
            hide(scantokens(Gttyp_[c_])$$(Gtform_[c_], loc_))
          endfor;
      fi
    endfor
  endfor
  Greset_ticks_;
  Gneedgr_:=false;
enddef;
string Glocs_, Glocs_.lft, Glocs_.rt, Glocs_.bot, Glocs_.top;
string Glocs_.llft, Glocs_.ulft, Glocs_.lrt, Glocs_.urt;
string Glocs_.lr, Glocs_.bt;
Glocs_ = "lft, rt, bot, top";
Glocs_.lft = "lft";
Glocs_.rt  = "rt" ;
Glocs_.bot = "bot";
Glocs_.top = "top";
Glocs_.llft = "lft, bot";
Glocs_.ulft = "lft, top";
Glocs_.lrt  = "rt, bot";
Glocs_.urt  = "rt, top";
Glocs_.lr = "lft, rt" ;
Glocs_.bt = "bot, top";
vardef fullgrid @#(text tnums_) text w_ =
  save n_,i_; i_=0; n_[1]=n_[2];
  for m_=tnums_: n_[incr i_]:=m_; endfor
  grid @#(n_[1] labeled lines llft, n_[1] noticks urt,
         n_[2] iticks) w_;
enddef;
vardef standardgrid @#(text tnums_) text w_ =
  save n_,i_; i_=0; n_[1]=n_[2]=0;
  for m_=tnums_: n_[incr i_]:=m_; endfor
  grid @#(n_[1] labeled iticks llft, n_[1] iticks urt,
         n_[2] small iticks) w_;
enddef;
vardef defaultgrid @#(text tnums_) text w_ =
  save n_,i_; i_=0; n_=0; 
  for m_=tnums_: n_:=m_; endfor
  grid @#(n_ labeled oticks llft) w_;
enddef;
vardef nogrid =
  Gneedgr_:=false;
enddef;
vardef Ggdef_(text mtyp_) =
  save c_; c_=0;
  for tn_=mtyp_: c_:=c_+1;
    if known tn_:
      Gtnum_[c_]=if tn_>0: tn_ else: ?? fi; 
    fi
  endfor
enddef;
def labeled = hide(Gtform_[Gcount_]:=autoform) enddef;
def small = 
  hide(Gtsc_[Gcount_]:=Gtscale_ if known Gtsc_[Gcount_]: *Gtsc_[Gcount_] fi)
enddef;
def noticks = Ggtdef_("notick") enddef;
def dots = Ggtdef_("dot") enddef;
def ticks = Ggtdef_("tick") enddef;
def iticks = Ggtdef_("itick") enddef;
def oticks = Ggtdef_("otick") enddef;
def lines = Ggtdef_("line") enddef;
def Ggtdef_(expr type) suffix loc =
  +0 hide(Gtsc_[Gcount_]:=
            if known Gtsc_[Gcount_]: Gtsc_[Gcount_] else: 1 fi;
          Gttyp_[Gcount_]:=type; Gtloc_[Gcount_]:=Glocs_ loc;
          Gcount_:=Gcount_+1)
enddef;
string Gemarks;      
Gemarks="20,10,5,2,1";
def Gfor_(text t) = for i=t endfor enddef;  
string Glmarks[];
Glmarks[1]="1,2,5";
Glmarks[2]="1,2,3,4,5,6,7,8,9";
Glmarks[3]="1Gfor_(6upto10:,i/5)Gfor_(5upto10:,i/2)Gfor_(6upto9:,i)";
Glmarks[4]="1Gfor_(11upto20:,i/10)Gfor_(11upto25:,i/5)Gfor_(11upto19:,i/2)";
Glmarks[5]="1Gfor_(21upto40:,i/20)Gfor_(21upto50:,i/10)Gfor_(26upto49:,i/5)";
string Gumarks;  
Gumarks="10,5,2";       % start with 10 and go down; a final `,1' is appended
vardef Gpick_ @# =
  save b; 
  forsuffixes $$=low,high:
    b $$ = Gpart_ @# Gr_ $$;
    Gmr_ $$ := Meform(b $$)+Gbias_;
  endfor
  if abs Gtype_ @#=log: if b.high-b.low >=mlog Gtmlg_: log else: linear fi
  else: linear
  fi
enddef;
pair Gbias_; Gbias_=(0,3);
vardef Gpart_.x expr p = xpart p enddef;
vardef Gpart_.y expr p = ypart p enddef;
def Gmsc_(expr k, l, h, f)(text t) =
  for m=if f: scantokens Glmarks[k] 
        elseif string k: scantokens k 
        else: for i=1 until 10-eps step 1/k: Gcma_ i endfor 
        fi:
    exitif m>h;
    if m>=l: t fi
  endfor
enddef;
def Gmscan_(expr k, f)(text t) =
  for e=ypart Gmr_.low upto ypart Gmr_.high:
    Gmsc_(k, if e>ypart Gmr_.low: 1 else: xpart Gmr_.low/1000 fi,
          if e<ypart Gmr_.high: 10 else: xpart Gmr_.high/1000 fi, f,  t)
  endfor
enddef;
vardef Gkpick_ =
  save k;
  k = 0;
  forever:
    exitif unknown Glmarks[k+1];
    exitif 0 Gmscan_(incr k, true, +1) >= Gmarks;
  endfor
  k
enddef;
vardef Gempick_ @# =
  numeric @#;
  for e=scantokens Gemarks:
    @# = e;
    exitif floor(ypart Gmr_.high/e)-floor(Gey_(Gmr_.low)/e) >= Gmarks;
    numeric @#;
  endfor
  known @#
enddef;
vardef Gey_(expr p) = ypart p  if xpart p=1000: -1 fi  enddef;
vardef Gipick_ =
  save m, n, d;
  m = Gmarks;
  n = 1 for i=1 upto 
    mlog(xpart Gmr_.high-xpart Gmr_.low)/ten - mlog m/(ten-epsilon):
        *10 endfor;
  if n<=1000:
    for x=scantokens Gumarks:
      d = n*x;
      exitif 0 Gigen_(d,+1)>=m;
      numeric d;
    endfor
  fi
  if known d: d else: n fi
enddef;
def Gigen_(expr d)(text t) =
  if numeric d:
    for m = d*ceiling(xpart Gmr_.low/d) step d until xpart Gmr_.high:
      t
    endfor
  else:
    for m = scantokens d:
      if (m>=xpart Gmr_.low) and (m<=xpart Gmr_.high): t fi
    endfor
  fi
enddef;
def Gemgen_(expr d)(text t) =
  for e = d*floor(Gey_(Gmr_.low)/d+1)
      step d until d*floor(ypart Gmr_.high/d):  t
  endfor
enddef;
vardef Gesame_ =
  save e;
  e+3 = if Gmr_.low=Gbias_: ypart Gmr_.high
        elseif Gmr_.high=Gbias_: ypart Gmr_.low
        else: max(ypart Gmr_.low, ypart Gmr_.high) 
        fi;
  forsuffixes $$=low, high:
    Gmr_ $$ := (xpart Gmr_ $$ for i=ypart Gmr_ $$ upto e+2: /(10) endfor, e);
  endfor
  e
enddef;
vardef Gpack_(expr m, e) =
  if m=0: zero else: Mlog m + e*ten fi
enddef;
string autoform; autoform = "%g";
def auto suffix $ =
  hide(freeze)
  hide(Fe_base:=Fen_base)
  hide(def Gcma_= hide(def Gcma_=,enddef) enddef) 
  if Gmarks>0:
    if Gpick_$=log:
      if Gempick_.Gme_:  hide(Fe_base:=Fee_base) 
        Gemgen_(Gme_, Gcma_ Gpack_(1,e))
      else: 
        Gmscan_(Gkpick_, true, Gcma_ Gpack_(m,e))
      fi
    else:
      hide(Gme_:=Gesame_)
      Gigen_(Gipick_, Gcma_ Gpack_(m,Gme_))
    fi
  fi
enddef;
vardef subframe @# =
  if str @#<>"":  subpath round(angle laboff @#*Gfra_+Gfrb_) of  fi theframe
enddef;
pair Gfra_, Gfrb_;  Gfra_=(1,1)/90;  Gfrb_=(.75,2.25);
vardef frame @# text w_ =
  if     str @#="lr": frame.lft w_; frame.rt w_; 
  elseif str @#="bt": frame.bot w_; frame.top w_;
  else:  addto currentpicture doublepath subframe @# withpen currentpen w_;
  fi  noframe;
enddef;
vardef noframe = Gneedfr_:=false; enddef;
