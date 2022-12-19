% gmleps.mp --
% Copyright (C) 1997,1998 by Wolfgang Kilian <kilian@x4u2.desy.de>
%   v0.40 1998/11/02
% This file lacks all documentation.  See the source gamelan.nw
% Note: This is a test version and not released code.  
%
% This file is called by the main program  gamelan.mp
picture gmlepspicture; 
string extra_begingmleps; extra_begingmleps = "";
string extra_endgmleps; extra_endgmleps = "";
boolean extepsverbose; extepsverbose = true;

%% String handling tool
string string_split[];
def splitstring expr S =
  begingroup
    save __splitctr; numeric __splitctr; __splitctr = 0;
    save __prevchar; string __prevchar, __currentchar; __prevchar = " ";
    for i = 0 upto infinity:
      __currentchar := substring(i, i+1) of S;
      if (__currentchar = " ") and (__prevchar = " "):
        relax;
      elseif (__currentchar <> " ") and (__prevchar = " "):
        string_split[__splitctr] := __currentchar;
      elseif (__currentchar <> " ") and (__prevchar <> " "):
        string_split[__splitctr] := string_split[__splitctr] & __currentchar;
      elseif (__currentchar = " ") and (__prevchar <> " "):
        __splitctr := __splitctr+1;
      fi
      __prevchar := __currentchar;
    endfor
  endgroup;
enddef;
%% End string handling tool

def begingmleps text F =
  begingroup;
    save file;                string file;        file = F;    
    save angle;               numeric angle;      angle = 0;
    save clipping;            boolean clipping;   clipping = false;
    save scale;               pair scale;         scale = (1,1);
    save base;                pair base;          base = origin;
    save __bbxfound;          boolean __bbxfound; __bbxfound = false;
    save grid;                boolean grid;       grid = false;
    save gridstep;            numeric gridstep;   gridstep = 10;
    save __base;              pair __base;
    save  __eps__currentline; string __eps__currentline;
    save __bbxline;           string __bbxline;
    save llx, lly, urx, ury;  numeric llx, lly, urx, ury;
    save pct;                 numeric pct;
    save width;               numeric width;
    save height;              numeric height;
    save clippath;            path clippath;
    save largefile;           boolean largefile;  largefile = false;
    save gridllx;             numeric gridllx;
    save gridlly;             numeric gridlly;
    save gridurx;             numeric gridurx;
    save gridury;             numeric gridury;
%% Finding the bounding box
    forever:
      __eps__currentline := readfrom F;
      if substring(0,14) of __eps__currentline = "%%BoundingBox:":
        exitif __eps__currentline = EOF; % PATCH D. Roegel 23-sep-2006
        __bbxline := substring(14, infinity) of __eps__currentline;
        __bbxfound := true;
        splitstring __bbxline;
        llx = scantokens string_split[0];
        lly = scantokens string_split[1];
        urx = scantokens string_split[2];
        ury = scantokens string_split[3];
      fi
      exitif __bbxfound;
    endfor
    if not __bbxfound:
      message "Warning:  No bounding box found.";
      message "          Setting bounding box to 0 0 1 1";
      llx = lly = 0;
      urx = ury = 1;
    fi
    closefrom F;
    scantokens extra_begingmleps;
    __base = -(llx,lly);
    pct = (urx - llx)/100;
%% To ensure the right bounding box of the output file
%% a picture with the same size as the eps figure is added.
    gmlepspicture := nullpicture;
    clippath = (0,0)--(0,ury-lly)--(urx-llx,ury-lly)--(urx-llx,0)--cycle;
    setbounds gmlepspicture to clippath;
enddef;
  
def endgmleps =
%% Calculating scale if width and/or height is known
  if (known width) and (known height):
    scale := (width/(urx - llx),height/(ury - lly));
  elseif known width:
    scale := (width/(urx - llx),width/(urx - llx));
  elseif known height:
    scale := (height/(ury - lly),height/(ury - lly));
  fi
%% The graphics inclusion commands
  special "gsave";
  if base <> origin:
    special decimal.xpart.base & " " & decimal.ypart.base & " translate";
  fi
  if angle <> 0:
    special decimal angle & " rotate";
  fi
  if __base <> origin:
    special decimal.xpart.__base & " " & decimal.ypart.__base & " translate";
  fi
  if scale <> (1,1):
    gmlepspicture := gmlepspicture scaled xpart.scale
    if xpart.scale <> ypart.scale:
      yscaled (ypart.scale/xpart.scale)
    fi;
    special decimal xpart.scale & " " & decimal ypart.scale & " scale";
  fi
  if angle <> 0:
    gmlepspicture := gmlepspicture rotatedaround(origin)(angle);
  fi
%% Drawing the grid !! After the scaling :-)
  if unknown gridllx:
    gridllx = 0;
  fi
  if unknown gridlly:
    gridlly = 0;
  fi
  if unknown gridurx:
    gridurx = (urx - llx)/pct;
  fi
  if unknown gridury:
    gridury = (ury - lly)/pct;
  fi
  if grid:
    save __gridpicture; picture __gridpicture; __gridpicture := nullpicture;
    for i = gridllx*pct step gridstep*pct until (epsilon + gridurx*pct):
      addto __gridpicture doublepath (i*xpart.scale,gridlly*pct*ypart.scale)--(i*xpart.scale,gridury*pct*ypart.scale) withpen currentpen;
      addto __gridpicture also thelabel.bot(((decimal.(i/pct) & "%") infont defaultfont) rotated -90, (i*xpart.scale,gridlly*pct*ypart.scale));
    endfor
    for i = gridlly*pct step gridstep*pct until (epsilon + gridury*pct):
      addto __gridpicture doublepath (gridllx*pct*xpart.scale,i*ypart.scale)--(gridurx*pct*xpart.scale,i*ypart.scale) withpen currentpen;
      addto __gridpicture also thelabel.lft(((decimal.(i/pct) & "%") infont defaultfont), (gridllx*pct*xpart.scale,i*ypart.scale));
    endfor
    if angle <> 0:
      __gridpicture := __gridpicture rotatedaround(origin)(angle);
    fi
    addto gmlepspicture also __gridpicture;
  fi
  if clipping:
    save __clippath; path __clippath;
    __clippath=clippath shifted (llx,lly);
    special "newpath " & decimal.xpart.point 0 of __clippath
    & " " & decimal.ypart.point 0 of __clippath & " moveto";
    for i = 0 upto length.__clippath-1:
      special decimal.xpart.postcontrol i of __clippath & " " &
      decimal.ypart.postcontrol i of __clippath & " " &
      decimal.xpart.precontrol (i+1) of __clippath & " " &
      decimal.ypart.precontrol (i+1) of __clippath & " " &
      decimal.xpart.point (i+1) of __clippath & " " &
      decimal.ypart.point (i+1) of __clippath & " curveto";
    endfor;
    special "closepath clip";
  fi
  special "save";
  special "userdict begin";
  special "/showpage {} def";
    special "%%BeginDocument: " & file;
    if largefile:
      special "%% MetaPost exteps large file->" & file;
      if extepsverbose:
        message "exteps notification:  File " & file & " not inserted into " & jobname & "." & decimal.charcode;
        message "                      Run 'delfin " & jobname & "." & decimal.charcode & "' to insert " & file;
        message "                      This is caused by setting 'largefile:=true'";
        message "";
      fi
    else:
      if extepsverbose:
        message "Inserting " & file & " into " & jobname & "." & decimal.charcode;
      fi
      forever:
        __eps__currentline := readfrom file;
        exitunless __eps__currentline <> EOF;      
        special __eps__currentline;
      endfor
    fi
    special "%%EndDocument: " & file;
    scantokens extra_endgmleps;
    special "end restore";
    special "grestore";
    closefrom file;
    if base <> (0,0):
      gmlepspicture :=  gmlepspicture shifted base;
    fi
    addto currentpicture also gmlepspicture;
  endgroup;
enddef;

%% Special drawing commands
def epsfill expr c = addto gmlepspicture contour c _op_ enddef;

def epsdraw expr p =
  addto gmlepspicture
  if picture p:
    also p
  else:
    doublepath p withpen currentpen
  fi
  _op_
enddef;

def epsfilldraw expr c =
  addto gmlepspicture contour c withpen currentpen
  _op_ enddef;

def epsdrawdot expr z =
  addto gmlepspicture contour makepath currentpen shifted z
  _op_ enddef;

def epslabel = epsdraw thelabel enddef;

endinput
