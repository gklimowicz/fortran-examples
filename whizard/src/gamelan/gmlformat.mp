picture Fmneg_, Femarker_, Feneg_, Fe_base, Fen_base, Fee_base, Fe_plus;
string Fmfont_, Fefont_;
numeric Fmscale_, Fescale_, Feraise_;
vardef init_numbers(expr s, m, x, sn, e) =
  Fmneg_ := s;
  for p within m: 
    Fmfont_ := fontpart p;
    Fmscale_ := xxpart p;
    exitif true;
  endfor
  Femarker_ := x;
  Feneg_ := sn;
  for p within e:
    Fefont_ := fontpart p;
    Fescale_ := xxpart p;
    Feraise_ := ypart llcorner p;
    exitif true;
  endfor
  Fen_base := Fline_up_("1" infont Fmfont_ scaled Fmscale_, Femarker_);
  Fee_base := "10" infont Fmfont_ scaled Fmscale_;
  Fe_base := Fen_base;
  Fe_plus := nullpicture;
enddef;
vardef Fmant_(expr x) =
  (decimal abs x infont Fmfont_ scaled Fmscale_)
enddef;
vardef Fexp_(expr x) =
  (decimal x infont Fefont_ scaled Fescale_ shifted (0,Feraise_))
enddef;
vardef Fline_up_(text t_) =
  save p_, c_;
  picture p_; p_=nullpicture;
  pair c_; c_=(0,0);
  for q_=t_:
    addto p_ also q_ if string q_: infont defaultfont scaled defaultscale fi
        shifted c_;
    c_ := (xpart lrcorner p_, 0);
  endfor
  p_
enddef;
vardef Fdec_o_(expr x) =
  if x<0: Fline_up_(Fmneg_, Fmant_(x))
  else: Fmant_(x)
  fi
enddef;
vardef Fsci_o_(expr x, e) =
  Fline_up_(if x<0: Fmneg_,fi
      if abs x =1:  Fe_base  else:  Fmant_(x), Femarker_  fi,
      if e<0:  Feneg_  else:  Fe_plus  fi,
      Fexp_(abs e))
enddef;
vardef Feform_(expr x) =
  if string x:
    if unknown zero: scantokens "input marith"; fi
    Meform(Mlog_str x)
  else:
    save b, e;
    b=x; e=0;
    if abs b>=10000: (b/10, 1)
    elseif b=0: origin
    else: forever:
      exitif abs b>=1000;
      b:=b*10; e:=e-1;
      endfor
     (b, e)
    fi
  fi
enddef;
vardef Feadj_(expr x, y) =
  if abs x>=10000: (x/10, y+1)   else: (x,y)  fi
enddef;
vardef roundd(expr x, d) =
  if abs d >4:
    if d>0: x else: 0 fi
  elseif d>0:
    save i; i=floor x;
    i + round(Ten_to[d]*(x-i))/Ten_to[d]
  else:
    round(x/Ten_to[-d])*Ten_to[-d]
  fi
enddef;
Ten_to0=1; Ten_to1=10; Ten_to2=100; Ten_to3=1000; Ten_to4=10000;
primarydef p Fprec_ k =
  Feadj_ (roundd(xpart p,k-4), ypart p)  
enddef;
primarydef p Fdigs_ k =
  Feadj_ (roundd(xpart p,k+ypart p), ypart p)  
enddef;
vardef Formsci_(expr z) = Fsci_o_(xpart z/1000, ypart z +3)  enddef;
vardef Formdec_(expr z) =
  if ypart z>0: Formsci_(z)
  else:
    Fdec_o_(xpart z
      if ypart z>=-4: /Ten_to[-ypart z]
      else: for i=ypart z upto -5: /(10) endfor/10000  
      fi)
  fi
enddef;
vardef Formgen_(expr q) =
  clearxy; (x,y)=q;
  if x=0: Formdec_
  elseif y>=-6: Formdec_
  else: Formsci_
  fi (q)
enddef;
def Fset_item_(expr s) =
  if s<>"":  s infont defaultfont scaled defaultscale,  fi
enddef;
string fmt_[];
fmt_[ASCII "e"] = "Formsci_(z Fprec_ p)"; % scientific, p significant figures
fmt_[ASCII "f"] = "Formdec_(z Fdigs_ p)"; % decimal, p digits right of the point
fmt_[ASCII "g"] = "Formgen_(z Fprec_ p)"; % decimal or scientific, p sig. figs.
fmt_[ASCII "G"] = "Formgen_(z Fdigs_ p)"; % decimal or scientific, p digits
vardef dofmt_ @#(expr f, x) =
  if (unknown f) or (f=""): nullpicture
  else:
    save k, l, s, p, z;
    pair z; z=@#(x);
    k = 1 + cspan(f, "%"<>);
    l-k = cspan(substring(k,infinity)of f, isdigit);
    p = if l>k:  scantokens substring(k,l)of f  else: 3 fi;
    string s; s=fmt_[ASCII substring (l,l+1) of f];
    if unknown s:
      if k<=length f: errmessage("No valid format letter found in "&f); fi
      s = "nullpicture";
    fi
    Fline_up_(
      Fset_item_(substring (0,k-1) of f)
      scantokens s,
      Fset_item_(substring (l+1,infinity) of f)
      nullpicture)
  fi
enddef;
vardef format(expr f, x) = dofmt_.Feform_(f,x) enddef;
vardef Mformat(expr f, x) = dofmt_.Meform(f,x) enddef;
