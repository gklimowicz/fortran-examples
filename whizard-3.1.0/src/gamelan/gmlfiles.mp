% gmlfiles.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
def fromfile expr f =
  eatcolon(Gtext_ beginfrom
    save currentfile;
    string currentfile; currentfile:=thedirectory & f & theextension;
    if known defaultkey: readkey defaultkey; fi
    if known defaultskip: skip defaultskip; fi   
    )
enddef;
def Gtext_(text t) = gobble begingroup t endgroup; enddef;
delimiters beginfrom endfrom;
numeric defaultskip;
def fileprefix=directory enddef;
def fileextension=extension enddef;
vardef key expr k =
  save l,$; l=length k; string $;
  forever: 
    $:=readfrom currentfile;  exitif k=take l of $;
    if $=EOF:  Gkeyerr_(k); fi  exitif $=EOF;
  endfor
enddef;
vardef skip expr sk =
  save $; string $;
  for i=1 upto sk: $:=readfrom currentfile;
    if $=EOF: Gskiperr_; fi exitif $=EOF;
  endfor
enddef;
def directory expr d =
  save thedirectory; string thedirectory; thedirectory:=d;
enddef;
def extension expr e =
  save theextension; string theextension; theextension:=e;
enddef;
def separator expr s = 
  save theseparator; string theseparator; theseparator=s; 
enddef;
def endkey expr k = 
  save theendkey; string theendkey; theendkey=k; 
enddef;
def prefix expr p = 
  save theprefix,Gpl_; string theprefix; theprefix=p; Gpl_=length p; 
enddef;
def comment expr c = 
  save thecomment,Gcl_; string thecomment; thecomment=c; Gcl_=length c; 
enddef;
directory ""; extension "";
separator " "; endkey ""; comment "#";
vardef endfile =
  save $; string $;
  forever:  $:=readfrom currentfile; exitif $=EOF;  endfor
enddef;
def readline =
  string $d; numeric $i; 
  forever:
    if known theprefix:
      forever: 
        $d:=readfrom currentfile; exitif theprefix=take Gpl_ of $d; 
        exitif $d=EOF;
      endfor
      if $d<>EOF:  $d := drop Gpl_ of $d; fi
    else: $d:=readfrom currentfile;
    fi
    exitif $d=EOF;
    string $d[]; $d0=$d; $d1=loptok $d0;
    exitif unknown thecomment;
    exitif thecomment<>take Gcl_ of $d1;
  endfor
  if $d=EOF: relax
  elseif $d1=theendkey: string $d;
  else: $i=1; 
    forever:  
      $d[incr $i]=loptok $d0; 
      if known thecomment: 
        if thecomment=take Gcl_ of $d[$i]: $d[$i]:=""; fi 
      fi
      exitif $d[$i]="";
    endfor 
    $d[$i]:=whatever_str;
  fi
  $i:=0;
enddef;
let \\ = readline;
vardef get text vars =
  forsuffixes  v=vars:  v:=#($d[incr $i]);  endfor
enddef;
def getn(expr n_) suffix v =
  for i_=1 upto n_: get v[i_]; endfor
enddef;
vardef getstring text vars =
  forsuffixes v=vars:  v:=$d[incr $i];  endfor
enddef;
def getnstring(expr n_) suffix v =
  for i_=1 upto n_: getstring v[i_]; endfor
enddef;
def withinblock =
  =  eatcolon(1 upto infinity: hide(readline)
              exitif unknown $d; exitif $d=EOF;
  )
enddef;
def Gkeyerr_(expr key) =
  errhelp("You should check the data file:  "
         &"Either it doesn't exist,                      "
         &"or it doesn't contain the key I was looking for.");
  show currentfile,key;
  errmessage("Key not found");
enddef;
def Gskiperr_ =
  errhelp("While executing a 'skip' command I reached EOF.  "
         &"Maybe the file doesn't exist.");
  show currentfile;
  errmessage("File ended while skipping lines");
enddef;
def Gfileerr_ =
  errhelp("The data file I was reading "
         &"either doesn't exist or doesn't contain            "
         &"any more data.  (Maybe you specified a wrong prefix.)");
  show currentfile;
  errmessage("File not found or empty");
enddef;
