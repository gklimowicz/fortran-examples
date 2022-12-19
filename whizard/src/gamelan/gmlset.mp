% gmlset.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
vardef Gdefset_ text @# =
  pair @#[], @#.@@;  boolean @#|;  numeric @#[]?, @#.$;
enddef;
vardef defined @# = if boolean @#| : true else: false fi enddef;
vardef cleared @# = known @#| enddef;
vardef Gset_ @#(expr typestr) =
  save v_,t_; string v_,t_; v_=templateof @#; t_=take 1 of v_;
  if (t_="") or (t_="["): Gnumtagerr_ t_;
  else:
    if tracingsets>0: message "Defining "&typestr&" dataset: "&v_; fi
    scantokens(typestr&" "&v_&"[]!;");
    scantokens("Gdefset_ "&v_);  
  fi
enddef;
newinternal tracingsets;
def Gnumtagerr_ expr t =
  errhelp("The first character in a dataset name should be a letter.");
  errmessage("Dataset name `"&t&"' empty or beginning with numeric suffix");
enddef;
def Gnotagerr_ =
  errhelp("Dataset tags must consist entirely of letters "
         &"(including underscores)          "
         &"and numeric suffixes.");
  errmessage("Illegal character in dataset tag");
enddef;
vardef clear(expr x) text t_ =
  forsuffixes @#=t_: 
    if not defined @#: Gset_ @#(typeof x); fi
    @#.$:=0;  @#[0]?:=0;  @#| :=true;  @#.@@:=(??,??);
    if tracingsets>0: message "Cleared dataset: "&str @#; fi
  endfor
enddef;
vardef clearn(expr n)(text tmp) text t_ =
  forsuffixes @#=t_:
    for i=0 upto n: clear(tmp) @#[i]; endfor
  endfor
enddef;
vardef empty @# =
  if cleared @#: if @#.$=0: true else: false fi
  else: Gnrerr_ @#; true
  fi
enddef;
vardef Gnrerr_ @# =
  errhelp("The dataset you referred to is not set up correctly. "&
          "Maybe you forgot to       "&
          "initialize it with clear.");
  show str @#; errmessage("Irregular dataset");
enddef;
vardef Gemptwarning_ @# =
  message("Warning: dataset "& str @# &" is empty.")
enddef;
vardef length_of @# = @#[@#.$]? enddef;
vardef partpoint @#(expr i_,j_) =
  save ii_,gi_; ii_=i_ mod @#.$; gi_=@#[ii_]?; 
  gi_ + j_ mod (@#[ii_+1]?-gi_)
enddef;
vardef endpoint @#(expr i_, r_) =
  partpoint @#(i_, if r_: -1 else: 0 fi)
enddef;
vardef partindex @#(expr j_) =
  save jj_,k_; jj_=j_ mod @#[@#.$]?; k_=0;
  forever: exitif @#[incr k_]?>jj_; endfor k_-1
enddef;
vardef augment @#(expr !!) text pp_ =
  if not cleared @#: clear(!!)@#; fi
  save l_,s_; l_=length_of @#; s_=@#.$;
  for p_=pp_:
    if known p_:
      @#[l_] := p_;
      @#[l_]! := !!;
      @#[if @#| : incr fi @#.$]? := @#[s_]? + 1;
      @#| :=false;
    else:
      @#| :=true;
    fi
  endfor
enddef;
vardef cut @# = augment @#(true) ?? enddef;
vardef checkdata(expr showerr_) text ss_ =
  save e_; boolean e_; e_=false;
  forsuffixes s_=ss_:
    if str s_ = "":
      if showerr_:
        errhelp "The empty string can't be a dataset tag.";
        errmessage "Dataset tag empty";
      fi
      e_:=true;
    elseif (known s_.$) and (numeric s_.$) and (known s_|) and (boolean s_|)
      and (pair s_[0]):
      if (s_.$>0) and (known length_of s_) and (numeric length_of s_):
        if (s_[0]?=0) and (length_of s_>0):
          for i_=1 upto s_.$: 
            if known s_[i_]?:  if s_[i_]?<=s_[i_-1]?: e_:=true; fi 
            else: e_:=true;
            fi
          endfor
          for i_=0 upto length_of s_ -1: if unknown s_[i_]: e_:=true; fi endfor
        else: e_:=true;
        fi
        if e_:
          errhelp "Sorry, but this dataset is somewhat screwed up, "
                 &"so I can't use the contents.";
          errmessage "Dataset "&str s_ &" is unreadable";
        fi
      else:
        if showerr_:
          errhelp "This dataset is empty.  "
                 &"                                                       "
                 &"Possibly a file was missing "
                 &"or empty when it should have been filled.";
          errmessage "Dataset "&str s_ &" is empty";
        fi
        e_:=true;
      fi
    else:
      if showerr_:
        errhelp "This dataset has not been defined. (Look for spelling errors.)";
        errmessage "Undefined dataset "&str s_;
      fi
      e_:=true;
    fi
  endfor
  (not e_)
enddef;
vardef showdata text ss_ =
  forsuffixes s_=ss_:
    save j_,l_,ip_; l_=length_of s_; ip_=0;
    if l_=0: message("Dataset "&str s_&" is empty.");
    else: scan(s_)(message("Contents of dataset "&str s_&":");j_=-1;)
      (message("Point "&if boolean $!: if $!: "*"& fi fi str[incr j_]&":" 
              &" ("&($#$ x)&", "&($#$ y)&")");
       if not boolean $!: show $!; fi
       if $~<2: message("End of part "&str[ip_]&"."); ip_:=ip_+1; fi)
      (message("End of dataset "&str s_&".  "&
               "Number of points: "&str[length_of s_]););
    fi
  endfor
enddef;
def plot(text s_) text t_ = 
  Gplot_(Gplpt_)(s_) t_
enddef;
vardef Gplot_(text cmd_)(text ss_) text pp_ =
  save n_,i_,j_,k_;  n_=0 forsuffixes #=ss_: +1 endfor; i_=-1;
  if n_>0: 
    for p_=pp_:  Gplset_(cmd_)(incr i_, n_)(ss_)(p_);  endfor
  else:
    errhelp("'plot' and 'hist' require a list of datasets as first "
             &"argument.");
    errmessage("Missing set tag");
  fi
enddef;
def Gplset_(text cmd_)(expr i_, n_)(text ss_)(expr p_) =
  j_:=i_ mod n_; k_:=-1;
  forsuffixes @#=ss_: exitif (incr k_)>j_; 
    if k_=j_: cmd_ @#(p_); fi
  endfor
enddef;
vardef Gplpt_ @#(expr p_) =
  if (numeric p_): cut @#;
  elseif unknown (xpart p_, ypart p_): cut @#;
  elseif pair p_:  augment @#(true) p_;
  else: 
    if known xxpart p_: cut @#;
      augment @#(false)(xxpart p_, ypart p_);
      augment @#(true)(xpart p_, ypart p_);
      augment @#(false)(xypart p_, ypart p_);
    fi
    if known yxpart p_: cut @#;
      augment @#(false)(xpart p_, yxpart p_);
      if unknown xxpart p_: augment @#(true)(xpart p_, ypart p_); fi
      augment @#(false)(xpart p_, yypart p_);
    fi
  fi
enddef;
tertiarydef p hbar b =
  begingroup save t; transform t;
    xpart t = xpart p;  ypart t = ypart p;
    xxpart t = (xpart t) minus if pair b: xpart fi b;
    xypart t = (xpart t) plus  if pair b: ypart fi b;
    if transform p:  yxpart t = yxpart p;  yypart t = yypart p;  fi
  t endgroup
enddef;
tertiarydef p vbar b =
  begingroup save t; transform t;
    xpart t = xpart p;  ypart t = ypart p;
    yxpart t = (ypart t) minus if pair b: xpart fi b;
    yypart t = (ypart t) plus  if pair b: ypart fi b;
    if transform p:  xxpart t = xxpart p;  xypart t = xypart p;  fi
  t endgroup
enddef;
tertiarydef p hrange b =
  begingroup save t; transform t;
    xpart t = xpart p;  ypart t = ypart p;
    xxpart t = if pair b: xpart fi b;
    xypart t = if pair b: ypart fi b;
    if transform p:  yxpart t = yxpart p;  yypart t = yypart p;  fi
  t endgroup
enddef;
tertiarydef p vrange b =
  begingroup save t; transform t;
    xpart t = xpart p;  ypart t = ypart p;
    yxpart t = if pair b: xpart fi b;
    yypart t = if pair b: ypart fi b;
    if transform p:  xxpart t = xxpart p;  xypart t = xypart p;  fi
  t endgroup
enddef;
def hist (text s_) text t_ = 
  Gplot_(Ghspt_)(s_) t_
enddef;
vardef Ghspt_ @#(expr p_) =
  if (unknown @#.@@) or (unknown xpart p_): cut @#;
  elseif pair p_: cut @#;
    augment @#(false) @#.@@; 
    augment @#(true)((xpart @#.@@ plus xpart p_) over two, ypart @#.@@);
    augment @#(false)(xpart p_, ypart @#.@@);
  else: 
    errhelp("For histogramming only pair expressions "
           &"are meaningful.                        "
           &"So, if you wanted to plot error bars, "
           &"you should use 'plot' instead.");
    errmessage("Wrong argument type for histogramming");
  fi
  @#.@@:=p_;
enddef;
def table = Gtable_ enddef;
vardef Gtable_.Gplot_(suffix cmd_)(suffix s_) text t_ =
  save x,y,xy_; Gdefxy_ t_;
  if known currentfile:
    vardef x = #($d1) enddef;
    vardef y @# = #($d[if str @#="": 2 else: @#+1 fi]) enddef;
    for l_ withinblock: $l:=l_; Gplot_(cmd_)(s_) xy_; endfor
  elseif Gfromf_:
    for x_ within Gfrange_: x:=x_; y:=function(x); Gplot_(cmd_)(s_) xy_; endfor
  else: Gtaberr_;
  fi
enddef;
def Gdefxy_(text t_) =
  save i_;  i_=0;
  for d_=t_: i_:=i_+1; endfor
  if i_=0:     def xy_ = (x,y) enddef;
  elseif i_=1: def xy_ = (x,t_) enddef;
  else:        def xy_ = (t_) enddef;
  fi
enddef;
def Gtaberr_ =
  errhelp("You must be inside of either a `fromfile' or `fromfunction' "&
          " block to use this macro.");
  errmessage("Can't use `table' here");
enddef;
def tables = Gtables_ enddef;
vardef Gtables_.Gplot_(suffix cmd_)(text ss_) text t_ =
  save x,y,xy_; Gdefxy_ t_;
  if known currentfile:
    vardef x = #($d1) enddef;
    vardef y = #($d[ii_]) enddef;
    for l_ withinblock: ii_:=1; $l:=l_;
      forsuffixes @# = ss_:  ii_:=ii_+1; Gplot_(cmd_)(@#)(xy_); endfor
    endfor
  else: Gmkploterr_;
  fi
enddef;
vardef calculate @#(text chains_) text xy_ =
  save Gs_;
  begingroup interim tracingsets:=0;
    scan pointwise(chains_)
      ()(x[$i]:=x; y[$i]:=y;
         if $i=$n: augment Gs_($!)(xy_); if $~=1: cut Gs_; fi fi)();
  endgroup;
  copy(Gs_) @#;
enddef;
vardef transpose(text chains_) text setlist =
  save Gs_;
  begingroup interim tracingsets:=0;
    scan pointwise(chains_)
      ()(augment Gs_($!)(x,y); if $i=$n: cut Gs_; fi)();
  endgroup;
  copy(Gs_) setlist;
enddef;
vardef copy(suffix s_) text setlist =
  forsuffixes @#=setlist:
    clear(s_[0]!)@#;
    for i=0 upto length_of.s_-1:  @#[i]:=s_[i]; @#[i]!:=s_[i]!; endfor
    for i=0 upto s_$: @#[i]?:=s_[i]?; endfor
    @#.$:=s_$;
  endfor
enddef;
def Gsample_(expr xmin,xmax,n) =
  if known xmin: 
    xmin
    if known xmax:
      if unknown n: hide(n=defaultsample) fi
      for i=1 upto n:
        ,(xmax times #(i/n)) plus (xmin times #(1-i/n)) 
      endfor
    fi
  fi
enddef;
newinternal defaultsample; defaultsample:=30;
def fromfunction suffix f = 
  def function=f enddef; Gfromfun_
enddef;
def Gfromfun_(text t) =
  def Gfrange_=(t) enddef;
  eatcolon(Gtext_ beginfrom save Gfromf_; boolean Gfromf_; Gfromf_=true;)
enddef;
let withinpicture=within;
def within expr p =
  if picture p: withinpicture p
  elseif color p:
    = Gsample_(redpart p, greenpart p, bluepart p)
  elseif pair p:
    = Gsample_(xpart p, ypart p, ??)
  elseif numeric p:
    = p
  elseif string p:
    = if length p>0: 
        substring(0,1) of p
        for i=2 upto length p: ,substring(i-1,i) of p endfor
      fi
  else:
    =
    hide(errhelp("I'm expecting a picture, string, numeric, pair, "&
                 "or tripel here.                So I will ignore this loop.");
         errmessage("Wrong type argument to within"))
  fi
enddef;
vardef Glexij_(suffix i_,j_)(expr c) =
  i_:=j_:=0;
  if c<>"":
    save found_; boolean found_; found_:=false;
    forever: i_:=i_+1; j_:=0;
      exitif unknown Glexc_[i_][1];
      forever: j_:=j_+1;
        exitif unknown Glexc_[i_][j_];
        if c=Glexc_[i_][j_]: found_:=true; fi
        exitif found_;
      endfor
      exitif found_;
    endfor
    if not found_: i_:=-1; j_:=0; fi
  fi
enddef;
