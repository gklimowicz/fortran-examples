% gmlscan.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
vardef scan @#(text chains_)(text ta_)(text ti_)(text te_) =
  Gscan_init_(chains_);
  if str @#="": 
    ta_ Gscan_(ti_) te_
  else:
    save piecewise,pointwise;
    let piecewise=Gscan_pc_;  let pointwise=Gscan_pt_;  
    ta_ @#(ti_) te_
  fi
enddef;
def Gscan_init_(text chains_) =
  save x,y,$; boolean $!; numeric $~,$i,$n;  
  begingroup save \,\\,{,};  Gparse_chain_ chains_; endgroup;
  Gq_init_;
  save f_; boolean f_; f_=false;
enddef;
def Gscan_(text ti_) =
  hide($i:=0)
  forever: exitif (incr $i)>$n;
    forever: hide(Gq_next_($i)) exitif Gqerr_; ti_ exitif $~=0; endfor
  endfor
enddef;
def Gscan_pc_(text ti_) =
  forever: hide($i:=0)
    forever: exitif (incr $i)>$n; 
      hide(Gdeftag_($i))
      forever: hide(Gq_next_($i)) exitif Gqerr_;
        if $~=0:     if $i<$n: hide($~:=1) fi  hide(f_:=true)
        elseif $~=1: if $i=$n: hide($~:=0) fi
        fi ti_
        exitif $~<2; 
      endfor
      exitif Gqerr_;
    endfor
    exitif f_ or Gqerr_;
  endfor
enddef;
def Gscan_pt_(text ti_) =
  forever: hide($i:=0)
    forever: exitif (incr $i)>$n; 
      hide(Gdeftag_($i))
      hide(Gq_next_($i)) exitif Gqerr_;
      if $~=0: hide(f_:=true) fi
      if $i<$n: hide($~:=2) fi
      ti_
    endfor
    exitif f_ or Gqerr_;
  endfor
enddef;
def Glexinit_ expr s =
  save lexs_,lexp_,lexq_,tok_,tcl_,tix_;
  string lexs_,tok_; lexs_=s; lexq_:=0;
enddef;
vardef Glex_ =
  save cl,ix;
  lexp_:=lexq_; lexq_:=lexp_+1;
  Glexij_(tcl_,tix_)(substring(lexp_,lexq_) of lexs_);
  if tcl_<>0:
    forever:
      Glexij_(cl,ix)(substring(lexq_, lexq_+1) of lexs_);
      exitif tcl_<>cl;
      lexq_:=lexq_+1;
    endfor
  fi
  tok_:=substring(lexp_,lexq_) of lexs_;
enddef;
string Glexc_[][];
Glexc_[1][1]="|";  Glexc_[1][2]="~";
Glexc_[2][1]="/";  Glexc_[2][2]="\";
Glexc_[3][1]="$"; Glexc_[4][1]="'";
Glexc_[5][1]="{"; Glexc_[6][1]="}";
vardef Gparse_chain_ text ss_ =
  save ql_,q_; q_=0;
  Gparseinit_;
  forsuffixes s_ = ss_: q_:=q_+1;
    Glexinit_ str s_;
    ql_:=0; Gqc_[q_][ql_]=0; Gqd_[q_][ql_]=false;
    forever: Gparse_;
      exitif Gqc_[q_][ql_]=0; 
    endfor
    ql_[q_]=ql_;
  endfor
  $n=q_;
enddef;
def Gparseinit_ =
  string Gqt_[][];
  numeric Gqc_[][],Gqp_[][]; boolean Gqd_[][],Gqsf_[][],Gqpf_[][];
enddef;
vardef Gparse_ =
  save qll_,perr_,sf_,pf_; boolean perr_,sf_,pf_;
  qll_=ql_;
  Glex_;
  if tcl_<0:
    forsuffixes @#=scantokens tok_:
      perr_:=not checkdata(true) @#;
      if not perr_:
        Gqt_[q_][ql_+1]=tok_; 
        Glex_;
        if (tcl_=1) or (tcl_=0):  
          Gparse_subchain_ @#(Ggen_subchain_(true,true));
        elseif tcl_=2: 
          Gparse_dir_(sf_,pf_); 
          Gparse_subchain_ @#(Ggen_subchain_(sf_,pf_)); Glex_;
        elseif tcl_=3:           
          Gparse_subchain_ @#(Ggen_pchain_);
        elseif tcl_=5:           
          Gparse_subchain_ @#(Gsubchain_); Glex_; 
        else: perr_:=true;
        fi
      fi
      if not perr_:
        if tcl_=1:     Gqc_[q_][ql_]:=tix_;
        elseif tcl_=0: Gqc_[q_][ql_]:=0;
        else: perr_:=true;
        fi
      fi
    endfor
  else: perr_:=true;
  fi
  if perr_: 
    errhelp "The chain shown above is strange: "
           &"Maybe it contains a syntax error, "
           &"or it      addresses an empty dataset.  "
           &"I'll try to recover, but all or part of the chain may be lost.";
    show lexs_;
    errmessage "Invalid chain"; 
    ql_:=qll_; Gqc_[q_][ql_]:=0; 
  fi
enddef;
vardef Gparse_dir_(suffix sf_,pf_) =
  save t_; string t_; t_=substring(1,2) of tok_;
  if tix_=1:     sf_:=true;
  elseif tix_=2: sf_:=false;
  else: perr_:=true;
  fi
  if     t_=Glexc_[2][1]: pf_:=true;
  elseif t_=Glexc_[2][2]: pf_:=false;
  else:  pf_:=sf_;
  fi
enddef;
vardef Ggen_subchain_(expr sf_, pf_) =
  Glexc_[3][1]& if sf_: "0" else: "[-1]" fi
 &Glexc_[2][ if sf_: 1 else: 2 fi ]
 &Glexc_[2][ if pf_: 1 else: 2 fi ]
 &Glexc_[3][1]& if sf_: "[-1]" else: "0" fi
enddef;
vardef Ggen_pchain_ =
  save p_; p_=lexp_;
  forever:  Glex_; exitif (tcl_=0) or (tcl_=1);  endfor;
  substring(p_,lexp_) of lexs_
enddef;
vardef Gsubchain_ =
  save p_; p_=lexq_;
  forever:  Glex_; exitif (tcl_=0) or (tcl_=6);  endfor;
  substring(p_,lexp_) of lexs_
enddef;
vardef Gparse_subchain_ @#(expr s_) =
  Glexinit_ s_;
  Glex_; Gchain_next_ @#;
  forever: exitif perr_;
    exitif Gqc_[q_][ql_]=0;
    Gchain_next_ @#;
  endfor
enddef;
vardef Gchain_next_ @# =
  ql_:=ql_+1;
  if tcl_<0:
    Gqp_[q_][ql_]=Gget_cpt_ @#; Glex_;
    Gset_con_; Gset_dir_;
  elseif tcl_=3: Glex_;
    Gsubchain_part_ @#;
  else: perr_:=true;
  fi
enddef;
vardef Gsubchain_part_ @# =
  if tcl_<0:
      save s_; s_=Gget_cpart_ @#;  Glex_;
% ..$a'b..
    if tcl_=4: Glex_;
      if tcl_<0: 
        Gqp_[q_][ql_]=partpoint @#(s_, number_from tok_); Glex_;
        Gset_con_; Gset_dir_;
      else: perr_:=true;
      fi
% ..$a~
    elseif (tcl_=1) or (tcl_=0):
      Gset_con_; Gset_dir_;
      if Gqd_[q_][ql_] and Gqd_[q_][ql_-1]: perr_:=true;   % //$a~// forbidden
      elseif Gqd_[q_][ql_]:
        Gqp_[q_][ql_]=endpoint @#(s_, not Gqpf_[q_][ql_]); % ..$a~//
      elseif Gqd_[q_][ql_-1]:
        Gqp_[q_][ql_]=endpoint @#(s_, Gqpf_[q_][ql_-1]);   % ..//$a~..
      else:
        Gqp_[q_][ql_]=endpoint @#(s_, false);              % ..$a~..
        Gqd_[q_][ql_]:=Gqsf_[q_][ql_]:=Gqpf_[q_][ql_]:=true;
        Gqd_[q_][ql_+1]=false;
        Gqc_[q_][ql_+1]=Gqc_[q_][ql_];  Gqc_[q_][ql_]:=2;
        Gqp_[q_][ql_+1]=endpoint @#(s_, true);
        ql_:=ql_+1;
      fi
% ..$a/..
    elseif tcl_=2:
      if Gqd_[q_][ql_-1]: perr_:=true;                     % //$a/.. forbidden
      else:
        Gset_dir_;
        if (tcl_<0) or (tcl_=3): Gset_con_;
          Gqp_[q_][ql_]=endpoint @#(s_, not Gqpf_[q_][ql_]); % $a//..
        elseif (tcl_=1) or (tcl_=0): Gset_con_;
          Gqp_[q_][ql_]=endpoint @#(s_, not Gqpf_[q_][ql_]); % $a/~.. 
          Gqd_[q_][ql_+1]=false;
          Gqc_[q_][ql_+1]=Gqc_[q_][ql_];  Gqc_[q_][ql_]:=2;
          Gqp_[q_][ql_+1]=endpoint @#(s_, Gqpf_[q_][ql_]);
          ql_:=ql_+1;
        else: perr_:=true;                                 % $a/~//  forbidden
        fi
      fi
    else: perr_:=true;
    fi
  else: perr_:=true;
  fi
enddef;        
vardef Gget_cpt_ @# =
  (number_from tok_) mod length_of @#
enddef;
vardef Gget_cpart_ @# =
  (number_from tok_) mod @#.$
enddef;
vardef Gset_con_ =
  if tcl_=1:     Gqc_[q_][ql_]=tix_; Glex_;
  elseif tcl_=0: Gqc_[q_][ql_]=0;
  else:          Gqc_[q_][ql_]=-1;
  fi
enddef;
vardef Gset_dir_ =
  if tcl_=2: Gqd_[q_][ql_]=true; 
    Gparse_dir_(Gqsf_[q_][ql_],Gqpf_[q_][ql_]); Glex_;
  else:      Gqd_[q_][ql_]=false;
  fi
enddef;
def Gq_init_ =
  numeric Gtc_[],Gtp_[],Gts_[],Gqpos_[];  
  boolean Gtn_[],Gtu_[];  string Gqtag_[];
  for q_ = 1 upto $n:
    Gqpos_[q_]=1; Gtn_[q_]=true; Gtu_[q_]=false; Gqtag_[q_]="";
  endfor
  boolean Gqerr_; Gqerr_=false;
enddef;
vardef Gq_next_(expr q_) =        
  save qi_,tp_,ts_,tc_; 
  qi_=Gqpos_[q_]; tp_=Gtp_[q_]; ts_=Gts_[q_]; tc_=Gtc_[q_];
  Gq_nextx_(q_);
  if not Gqerr_:
    Gqpos_[q_]:=qi_;
    Gtc_[q_]:=tc_; Gtp_[q_]:=tp_; Gts_[q_]:=ts_;
    numeric x,y; z=Gtag_[tp_]; $!:=Gtag_[tp_]!; $~:=tc_;
  else:
    errhelp "There are no points to scan here, possibly a result "
           &"of a previous error.       "
           &"Just type return, and I will ignore this.";
    errmessage "Empty scanning queue"; Gtc_[q_]:=0; 
  fi
enddef;
def Gq_nextx_(expr q_) =
  if Gtn_[q_]: Gtn_[q_]:=false;
    if known Gqt_[q_][qi_]: 
      Gqtag_[q_]:=Gqt_[q_][qi_]; Gdeftag_(q_);
    fi
    tp_:=Gqp_[q_][qi_]; tc_:=Gqc_[q_][qi_];
    if (known tp_) and (known tc_):
      ts_:=partindex Gtag_(tp_);
      if Gqd_[q_][qi_]: Gqset_con_(q_, Gqpf_[q_][qi_], Gqc_[q_][qi_]);
      else: qi_:=qi_+1; Gtn_[q_]:=true; 
      fi
    else: Gqerr_:=true;
    fi
  else:
    if tp_<>Gqp_[q_][qi_+1]: 
      Gqinc_tp_(q_, Gqpf_[q_][qi_], Gqsf_[q_][qi_]); 
      Gqset_con_(q_, Gqpf_[q_][qi_], Gqc_[q_][qi_]); 
    fi
    if tp_=Gqp_[q_][qi_+1]:
      qi_:=qi_+1; Gtn_[q_]:=true; Gq_nextx_(q_);
    fi
  fi
enddef;
def Gdeftag_(expr q_) =
  scantokens("def Gtag_ = "&Gqtag_[q_]&" enddef;");
enddef;
def Gqset_con_(expr q_,pf_,c_) =
  Gtu_[q_]:=
    if pf_:  if tp_+1 = Gtag_[ts_+1]?: true else: false fi
    else:    if tp_   = Gtag_[ts_]?:   true else: false fi
    fi;  
  if c_<0:  tc_:= if Gtu_[q_]: 1 else: 2 fi;  fi
enddef;
def Gqinc_tp_(expr q_, pf_,sf_) =
  if pf_:  
    tp_:=tp_+1;  
    if Gtu_[q_]: Gqinc_ts_(q_, sf_); tp_:=Gtag_[ts_]?; fi
  else:    
    if Gtu_[q_]: Gqinc_ts_(q_, sf_); tp_:=Gtag_[ts_+1]?; fi  
    tp_:=tp_-1;
  fi
enddef;
def Gqinc_ts_(expr q_, sf_) =
  if sf_:  ts_:=ts_+1;  if ts_=Gtag_.$: ts_:=0; fi
  else:    if ts_=0: ts_:=Gtag_.$; fi  ts_:=ts_-1;
  fi
enddef;
