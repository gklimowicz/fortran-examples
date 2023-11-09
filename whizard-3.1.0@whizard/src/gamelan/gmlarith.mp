% gmlarith.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
newinternal zero,one,two,ten,pi; 
zero:=-16384; one:=0; two:=177.44568; ten:=589.46179; pi:=293.07803;
vardef Mabs primary x = x*.5*2 enddef;
vardef Mlog primary x =
  if     x>0: Mabs mlog x
  elseif x<0: epsilon + Mabs mlog(-x)
  else: zero
  fi
enddef;
vardef Mexp primary x =
  if x=Mabs x: mexp x
  else: -mexp x
  fi
enddef;
vardef # primary u =
  if pair u: (# xpart u, # ypart u)
  elseif unknown u: ??
  elseif numeric u: Mlog u
  elseif string u:
    for s=take 1 of u:
      if s="-": -epsilon + Gnum_ drop 1 of
      elseif s="+": Gnum_ drop 1 of
      else: Gnum_
      fi u
    endfor
  else: 
    begingroup show u;
      errhelp("The # operator can handle only strings, numerics, and pairs.");
      errmessage("Not a number");
      ??
     endgroup
  fi
enddef;
def ##=# enddef;
vardef Gism_ primary m =
  (isdigit m)or(m=".")
enddef;
vardef Gnum_ primary u =
  save n,p; p=cspan(u,isdigit)+1;   
  if substring(p-1,p) of u=".": 
    n= p + cspan(drop p of u, isdigit)+1;
    Gmm_((take p-1 of u)&(substring(p,n-1) of u), p, drop n of u)
  else: Gmm_(take p-1 of u, p, drop p of u)
  fi
enddef;
vardef Gmm_(expr m,p,e) =
  save l,n,s; l=length m; n=cspan(m,"0"=); string s; 
  s= if n>0: drop n if n=l: -1 fi of fi m; 
  if s="0": zero
  elseif s="": ??
  else:
    Gee_(Mabs(mlog(scantokens((take 4 of (s&"000"))&"."&(drop 4 of s))))
         + (p-l-5+length s)*ten + 2epsilon)(e)
  fi
enddef;
vardef Gee_(expr m,u) =
  if u="": m
  else: save n,e; string e;
    e = for s=take 1 of u:
          if (s="-")or(s="+"): drop 1 of fi u
        endfor;
    if cspan(e,isdigit)<>length e: ??
    else: n=m+(scantokens u)*ten;
      if n<zero: zero else: n fi
    fi
  fi
enddef;
vardef $#$ primary u =
  if unknown u: begingroup save ?; string ?; ? endgroup
  else:
    save m,n,e,i,a,b,f,s; (m,e)=Meform(u);
    if m=0: "0"
    else:  
      n=abs m; i=floor n; f=n-i; i:=i/1000; a=floor i; b=i-a; 
      string s; s="00"&decimal abs (e+3);
      if m<0: "-"& fi
      decimal a & "." & (substring(2,5) of (decimal b &"0000"))
                      & (substring(2,7) of (decimal f &"000000"))
                & if e>-4: "E+" else: "E-" fi 
                & (drop (length s) - 2 of s)
    fi
  fi
enddef;
def $##$=$#$ enddef;
vardef Meform(expr q) =
  if q<=zero: (0,0)
  else:
    save e; e=floor((q-1768.38985)/ten);
    (Mexp(q-e*ten), e)
  fi
enddef;
def #$=Mexp enddef;
def ##$=Mexp enddef;
vardef showfloat text t_ = 
  for xf_=t_: 
    if numeric xf_: show $#$ xf_; 
    elseif pair xf_: showfloat xpart xf_; showfloat ypart xf_;
    else: show xf_; 
    fi 
  endfor 
enddef;
def floatabs=Mabs enddef;
vardef neg primary a =
  if unknown a: whatever
  elseif a=zero: zero
  elseif a=floatabs a: a-epsilon
  else: a+epsilon
  fi
enddef;
vardef sign primary u =
  if u=zero: zero
  elseif u=floatabs u: 0
  else: epsilon
  fi
enddef; 
tertiarydef a less b =
  if a=floatabs a:  if b=floatabs b: a<b else: false fi
  elseif b=floatabs b: true
  else: b<a
  fi
enddef;
tertiarydef a lesseq b =
  if a=b: true
  else: a less b
  fi
enddef;
tertiarydef a greatereq b = (b lesseq a) enddef;
tertiarydef a greater b = (b less a) enddef;
let equal = =; let unequal = <>;
primarydef a times b =
  if (unknown a) or (unknown b): whatever
  elseif (a<=zero) or (b<=zero): zero
  else: (a+b)
  fi
enddef;
primarydef a over b =
  if (unknown a) or (unknown b): whatever
  elseif b=zero: begingroup
                    errhelp "I can't divide by zero or by very "
                           &"small numbers.  I'll insert 'infinity'.";
                    errmessage "Division by zero";
                    -2epsilon-zero endgroup
  elseif a<=zero: zero
  else: (a-b)
  fi
enddef;
secondarydef a plus b =
  if (unknown a) or (unknown b): whatever
  elseif a = zero: b
  elseif b = zero: a
  elseif a >= b: (Mlog(1579 + Mexp(b times (1885.32404-a))) + a-1885.32404)
  else:  b plus a
  fi
enddef;
secondarydef a minus b = a plus (b-epsilon) enddef;
primarydef a ^ b = if a=zero: zero else: b*a fi enddef;
primarydef a to_the b = 
  if (unknown a) or (unknown b): whatever
  elseif a=zero: 
    if b=floatabs b: zero
    else:
      begingroup interim showstopping:=0; showfloat a;
        errhelp "#0 to_the b is illegal if b is "
               &"negative.  I'll insert 'infinity'.";
        errmessage "Negative power of zero is undefined";
        -2epsilon-zero 
      endgroup
    fi
  elseif a=floatabs a: (Mexp b)*a
  else:
    begingroup interim showstopping:=0; showfloat a;
      errhelp "a to_the b is undefined if a is "
             &"negative.  I'll use the absolute value of a "
             &"     instead.  Nevertheless, if b is an "
             &"integer, you can use a^b.";
      errmessage "Power of negative argument is undefined";
      (floatabs a) to_the b 
    endgroup
  fi
enddef;
vardef logn primary y =
  if unknown y: whatever
  elseif y=zero: errhelp "I can't take the logarithm of zero or very "
                       &"small numbers.  I'll insert         '-infinity'.";
                  errmessage "Logarithm of 0 is undefined";
                  -epsilon-zero
  elseif y<>floatabs y: interim showstopping:=0; showfloat y; 
                  errhelp "I can't take the logarithm of negative "
                         &"numbers.  I'll take log(abs(x)) instead.";
                  errmessage "Logarithm of negative argument is undefined";
                  logn(floatabs y)
  else: Mlog(y/256)
  fi
enddef;
newinternal logten; logten:=logn ten;
vardef logd primary y = (logn y over logten) enddef;
vardef exp primary y =
  if unknown y: whatever
  else: 256 Mexp y
  fi
enddef;
vardef root primary y = 
  if unknown y: whatever
  elseif y=zero: zero
  elseif y<>floatabs y: interim showstopping:=0; showfloat y;
                  errhelp "I can't take the square root of negative "
                         &"numbers.  I'll take sqrt(abs(x))      instead.";
                  errmessage "Square root of negative argument is undefined";
                  floatabs(.5 y)
  else: floatabs(.5 y)
  fi
enddef;
vardef inverse primary y =
  if unknown y: whatever
  elseif y=zero: errhelp "I can't take the inverse of zero or very "
                         &"small numbers. I'll insert 'infinity'.";
                  errmessage "Division by zero";
                  -2epsilon-zero 
  else: -y
  fi
enddef;
vardef sinh primary y = (exp y minus exp neg y) over two enddef;
vardef cosh primary y = (exp y plus exp neg y) over two enddef;
vardef tanh primary y = (sinh y) over (cosh y) enddef;
vardef coth primary y = (cosh y) over (sinh y) enddef;
vardef arsinh primary y = logn(y plus root(y^2 plus 0)) enddef;
vardef arcosh primary y = logn(y plus root(y^2 minus 0)) enddef;
vardef artanh primary y = logn root((0 plus y) over (0 minus y)) enddef;
vardef arcoth primary y = logn root((y plus 0) over (y minus 0)) enddef;
def arcsinh=arsinh enddef;
def arccosh=arcosh enddef;
def arctanh=artanh enddef;
def arccoth=arcoth enddef;
vardef sin primary y = Mlog sind(57.29578 Mexp y) enddef;
vardef cos primary y = Mlog cosd(57.29578 Mexp y) enddef;
vardef tan primary y = (sin y) over (cos y) enddef;
vardef cot primary y = (cos y) over (sin y) enddef;
vardef arcsin primary y = arctan(y over root(0 minus y^2)) enddef;  
vardef arccos primary y = arccot(y over root(0 minus y^2)) enddef;
vardef arctan primary y = Mlog angle(1,Mexp y) over 1036.3461 enddef;
vardef arccot primary y = Mlog angle(Mexp y,1) over 1036.3461 enddef;
