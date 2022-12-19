%%
%% This is file `flex.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% flex.dtx  (with options: `class')
%% 
%% Copyright (C) 1997 by Thorsten.Ohl@Physik.TH-Darmstadt.de
%% 
%% This file is NOT the source for flex, because almost all comments
%% have been stripped from it.  It is NOT the preferred form of flex
%% for making modifications to it.
%% 
%% Therefore you can NOT redistribute and/or modify THIS file.  You can
%% however redistribute the complete source (flex.dtx and flex.ins)
%% and/or modify it under the terms of the GNU General Public License as
%% published by  the Free Software Foundation; either version 2, or (at
%% your option) any later version.
%% 
%% As a special exception, you can redistribute parts of this file for
%% the electronic distribution of scientific papers, provided that you
%% include a short note pointing to the complete source.
%% 
%% Flex is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% \CheckSum{1603}
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\fileversion{v0.00}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
{\def\RCS#1#2\endRCS{%
  \ifx$#1%
    \@RCS $#2 \endRCS
  \else
    \@RCS $*: #1#2$ \endRCS
  \fi}%
 \def\@RCS $#1: #2,v #3 #4 #5 #6 #7$ \endRCS{%
   \gdef\filename{#2}%
   \gdef\filerevision{#3}%
   \gdef\filedate{#4}%
   \gdef\filemaintainer{#6}}%
\RCS $Id: flex.dtx,v 1.2 1997/06/22 17:50:51 ohl Exp $ \endRCS}%
\ProvidesPackage{flex}[\filedate\space\fileversion\space
  FLEXible Document Class for LaTeX2e (\filemaintainer)]
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagetrue
\newif\if@openright
\DeclareOption{a4paper}
   {\setlength{\paperheight}{297mm}\setlength{\paperwidth}{210mm}}
\DeclareOption{a5paper}
   {\setlength{\paperheight}{210mm}\setlength{\paperwidth}{148mm}}
\DeclareOption{b5paper}
   {\setlength{\paperheight}{250mm}\setlength{\paperwidth}{176mm}}
\DeclareOption{letterpaper}
   {\setlength{\paperheight}{11in}\setlength{\paperwidth}{8.5in}}
\DeclareOption{legalpaper}
   {\setlength{\paperheight}{14in}\setlength{\paperwidth}{8.5in}}
\DeclareOption{executivepaper}
   {\setlength{\paperheight}{10.5in}\setlength{\paperwidth}{7.25in}}
\DeclareOption{landscape}
   {\setlength{\@tempdima}{\paperheight}%
    \setlength{\paperheight}{\paperwidth}%
    \setlength{\paperwidth}{\@tempdima}}
\DeclareOption{10pt}{\renewcommand{\@ptsize}{0}}
\DeclareOption{11pt}{\renewcommand{\@ptsize}{1}}
\DeclareOption{12pt}{\renewcommand{\@ptsize}{2}}
\DeclareOption{oneside}{\@twosidefalse\@mparswitchfalse}
\DeclareOption{twoside}{\@twosidetrue\@mparswitchtrue}
\DeclareOption{draft}{\setlength{\overfullrule}{5pt}}
\DeclareOption{final}{\setlength{\overfullrule}{0pt}}
\DeclareOption{titlepage}{\@titlepagetrue}
\DeclareOption{notitlepage}{\@titlepagefalse}
\DeclareOption{openright}{\@openrighttrue}
\DeclareOption{openany}{\@openrightfalse}
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
    \renewcommand{\@openbib@code}{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep\z@}%
    \renewcommand{\newblock}{\par}}}
\DeclareOption{chapters}{\@flx@chapterstrue}
\DeclareOption{nochapters}{\@flx@chaptersfalse}
\newif\if@flx@chapters
\@flx@chaptersfalse
\ExecuteOptions{letterpaper,10pt,oneside,onecolumn,final,openany}
\ProcessOptions
\input{size1\@ptsize.clo}
\setlength{\lineskip}{1\p@}
\setlength{\normallineskip}{1\p@}
\renewcommand{\baselinestretch}{}
\setlength{\parskip}{0\p@\@plus\p@}
\@lowpenalty=51
\@medpenalty=151
\@highpenalty=301
\setcounter{topnumber}{2}
\renewcommand{\topfraction}{.7}
\setcounter{bottomnumber}{1}
\renewcommand{\bottomfraction}{.3}
\setcounter{totalnumber}{3}
\renewcommand{\textfraction}{.2}
\renewcommand{\floatpagefraction}{.5}
\setcounter{dbltopnumber}{2}
\renewcommand{\dbltopfraction}{.7}
\renewcommand{\dblfloatpagefraction}{.5}
\if@twoside
  \def\ps@headings{%
    \let\@oddfoot\@empty
    \let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\markboth
    \def\chaptermark##1{%
      \markboth{\MakeUppercase{%
        \ifnum\c@secnumdepth>\m@ne
           \@chapapp\ \thechapter. \ %
        \fi
        ##1}}{}}%
    \def\sectionmark##1{%
      \markright{\MakeUppercase{%
        \ifnum\c@secnumdepth>\z@
          \thesection. \ %
        \fi
        ##1}}}}
\else
  \def\ps@headings{%
    \let\@oddfoot\@empty
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\markboth
    \def\chaptermark##1{%
      \markright{\MakeUppercase{%
        \ifnum\c@secnumdepth>\m@ne
           \@chapapp\ \thechapter. \ %
        \fi
        ##1}}}}
\fi
\def\ps@myheadings{%
  \let\@oddfoot\@empty
  \let\@evenfoot\@empty
  \def\@evenhead{\thepage\hfil\slshape\leftmark}%
  \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
  \let\@mkboth\@gobbletwo
  \let\chaptermark\@gobble
  \let\sectionmark\@gobble}
\if@titlepage
  \newcommand{\maketitle}{%
    \begin{titlepage}%
      \let\footnotesize\small
      \let\footnoterule\relax
      \let\footnote\thanks
      \null\vfil
      \vskip 60\p@
      \begin{center}%
        {\LARGE \@title \par}%
        \vskip 3em%
        {\large
         \lineskip .75em%
          \begin{tabular}[t]{c}%
            \@author
          \end{tabular}\par}%
          \vskip 1.5em%
        {\large \@date \par}%
      \end{center}\par
      \@thanks
      \vfil\null
    \end{titlepage}%
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax}
\else
  \newcommand{\maketitle}{\par
    \begingroup
      \renewcommand{\thefootnote}{\@fnsymbol\c@footnote}%
      \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
      \long\def\@makefntext##1{%
        \parindent 1em\noindent
        \hb@xt@1.8em{\hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
      \if@twocolumn
        \ifnum\col@number=\@ne
          \@maketitle
        \else
          \twocolumn[\@maketitle]%
        \fi
      \else
        \newpage
        \global\@topnum\z@   % Prevents figures from going at top of page.
        \@maketitle
      \fi
      \thispagestyle{plain}\@thanks
    \endgroup
    \setcounter{footnote}{0}%
    \global\let\thanks\relax
    \global\let\maketitle\relax
    \global\let\@maketitle\relax
    \global\let\@thanks\@empty
    \global\let\@author\@empty
    \global\let\@date\@empty
    \global\let\@title\@empty
    \global\let\title\relax
    \global\let\author\relax
    \global\let\date\relax
    \global\let\and\relax}
  \def\@maketitle{%
    \newpage
    \null
    \vskip 2em%
    \begin{center}%
      \let\footnote\thanks
      {\LARGE \@title \par}%
      \vskip 1.5em%
      {\large
        \lineskip .5em%
        \begin{tabular}[t]{c}%
          \@author
        \end{tabular}\par}%
      \vskip 1em%
      {\large \@date}%
    \end{center}%
    \par
    \vskip 1.5em}
\fi
\newcommand{\flxDefFont}[2]{\@namedef{flx:font:#1}{#2}}
\newcommand{\flxFont}[1]{\@nameuse{flx:font:#1}}
\newcommand{\flxDefSize}[2]{\@namedef{flx:Size:#1}{#2}}
\newcommand{\flxSize}[1]{\@nameuse{flx:Size:#1}}
\newcommand{\flxDefAlign}[2]{\@namedef{flx:align:#1}{#2}}
\newcommand{\flxAlign}[1]{\@nameuse{flx:align:#1}}
\newcommand{\flxDefFontAll}[1]{%
  \flxDefFont{part}{#1}%
  \flxDefFont{chapter}{#1}%
  \flxDefFont{section}{#1}%
  \flxDefFont{subsection}{#1}%
  \flxDefFont{subsubsection}{#1}%
  \flxDefFont{paragraph}{#1}%
  \flxDefFont{subparagraph}{#1}}
\newcommand{\flxDefSizeAll}[1]{%
  \flxDefSize{part}{#1}%
  \flxDefSize{chapter}{#1}%
  \flxDefSize{section}{#1}%
  \flxDefSize{subsection}{#1}%
  \flxDefSize{subsubsection}{#1}%
  \flxDefSize{paragraph}{#1}%
  \flxDefSize{subparagraph}{#1}}
\newcommand{\flxDefAlignAll}[1]{%
  \flxDefAlign{part}{#1}%
  \flxDefAlign{chapter}{#1}%
  \flxDefAlign{section}{#1}%
  \flxDefAlign{subsection}{#1}%
  \flxDefAlign{subsubsection}{#1}%
  \flxDefAlign{paragraph}{#1}%
  \flxDefAlign{subparagraph}{#1}}
\newcommand{\flxDefNumber}[2]{\@namedef{flx:num:#1}{#2}}
\newcommand{\flxNumber}[1]{\@nameuse{flx:num:#1}}
\def\flxChapter{\@chapapp}
\newcommand{\flx@startsection}[5]{%
  \@startsection{#1}{#2}{#3}{#4}{#5}%
    {\normalfont\flxSize{#1}\flxFont{#1}\flxAlign{#1}}}
\newcommand*{\chaptermark}[1]{}
\setcounter{secnumdepth}{2}
\newcounter{part}
\renewcommand{\thepart}{\@Roman\c@part}
\newcommand{\part}{%
  \if@openright
    \cleardoublepage
  \else
    \clearpage
  \fi
  \thispagestyle{plain}%
  \if@twocolumn
    \onecolumn
    \@tempswatrue
  \else
    \@tempswafalse
  \fi
  \null\vfil
  \secdef\@part\@spart}
\def\@part[#1]#2{%
    \ifnum\c@secnumdepth>-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    \markboth{}{}%
    {\interlinepenalty\@M
     \normalfont\flxSize{chapter}\flxFont{chapter}\flxAlign{chapter}%
     \ifnum\c@secnumdepth>-2\relax
       {\flxNumber{part}}%
     \fi
     #2\par}%
    \@endpart}
\def\@spart#1{%
    {\interlinepenalty\@M
     \normalfont\flxSize{chapter}\flxFont{chapter}\flxAlign{chapter}%
     #1\par}%
    \@endpart}
\def\@endpart{%
  \vfil\newpage
  \if@twoside
    \null
    \thispagestyle{empty}%
    \newpage
  \fi
  \if@tempswa
    \twocolumn
  \fi}
\newcommand*{\l@part}[2]{%
  \ifnum\c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength{\@tempdima}{3em}%
      \parindent\z@
      \rightskip\@pnumwidth
      \parfillskip-\@pnumwidth
      {\leavevmode
       \large\flxFont{part}%
       #1\hfil\hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \global\@nobreaktrue
       \everypar{\global\@nobreakfalse\everypar{}}%
    \endgroup
  \fi}
\if@flx@chapters
  \newcommand{\@chapapp}{\chaptername}
  \newcounter{chapter}
  \renewcommand{\thechapter}{\@arabic\c@chapter}
  \newcommand{\chapter}{%
    \if@openright
      \cleardoublepage
    \else
      \clearpage
    \fi
    \thispagestyle{plain}%
    \global\@topnum\z@
    \@afterindentfalse
    \secdef\@chapter\@schapter}
  \def\@chapter[#1]#2{%
    \ifnum\c@secnumdepth>\m@ne
      \refstepcounter{chapter}%
      \typeout{\@chapapp\space\thechapter.}%
      \addcontentsline{toc}{chapter}{\protect\numberline{\thechapter}#1}%
    \else
      \addcontentsline{toc}{chapter}{#1}%
    \fi
    \chaptermark{#1}%
    \addtocontents{lof}{\protect\addvspace{10\p@}}%
    \addtocontents{lot}{\protect\addvspace{10\p@}}%
    \if@twocolumn
      \@topnewpage[\@makechapterhead{#2}]%
    \else
      \@makechapterhead{#2}%
      \@afterheading
    \fi}
  \def\@makechapterhead#1{%
    \vspace*{50\p@}%
    {\parindent\z@
     \normalfont\flxSize{chapter}\flxFont{chapter}\flxAlign{chapter}%
     \ifnum\c@secnumdepth>\m@ne
       {\flxNumber{chapter}}%
     \fi
     \interlinepenalty\@M
     #1\par\nobreak
     \vskip 40\p@}}
  \def\@schapter#1{%
    \if@twocolumn
      \@topnewpage[\@makeschapterhead{#1}]%
    \else
      \@makeschapterhead{#1}%
      \@afterheading
    \fi}
  \def\@makeschapterhead#1{%
    \vspace*{50\p@}%
    {\parindent\z@
     \interlinepenalty\@M
     \normalfont\flxSize{chapter}\flxFont{chapter}\flxAlign{chapter}%
     #1\par\nobreak
     \vskip 40\p@}}
  \newcommand*\l@chapter[2]{%
    \ifnum \c@tocdepth >\m@ne
      \addpenalty{-\@highpenalty}%
      \vskip 1.0em \@plus\p@
      \setlength\@tempdima{1.5em}%
      \begingroup
        \parindent \z@ \rightskip \@pnumwidth
        \parfillskip -\@pnumwidth
        \leavevmode
        \flxFont{chapter}%
        \advance\leftskip\@tempdima
        \hskip -\leftskip
        #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
        \penalty\@highpenalty
      \endgroup
    \fi}
\fi
\if@flx@chapters
  \newcounter{section}[chapter]
  \renewcommand{\thesection}{\thechapter.\@arabic\c@section}
\else
  \newcounter{section}
  \renewcommand{\thesection}{\@arabic\c@section}
\fi
\newcommand{\section}{%
  \flx@startsection{section}{1}%
    {\z@}%
    {-3.5ex \@plus -1ex \@minus -.2ex}%
    {2.3ex \@plus .2ex}}
\newcommand*{\l@section}{%
  \flxFont{section}\@dottedtocline{1}{1.5em}{2.3em}}
\newcounter{subsection}[section]
\renewcommand{\thesubsection}{\thesection.\@arabic\c@subsection}
\newcommand{\subsection}{%
  \flx@startsection{subsection}{2}%
    {\z@}%
    {-3.25ex\@plus -1ex \@minus -.2ex}%
    {1.5ex \@plus .2ex}}
\newcommand*{\l@subsection}{%
  \flxFont{subsection}\@dottedtocline{2}{3.8em}{3.2em}}
\newcounter{subsubsection}[subsection]
\renewcommand{\thesubsubsection}{\thesubsection .\@arabic\c@subsubsection}
\newcommand{\subsubsection}{%
  \flx@startsection{subsubsection}{3}%
    {\z@}%
    {-3.25ex\@plus -1ex \@minus -.2ex}%
    {1.5ex \@plus .2ex}}
\newcommand*{\l@subsubsection}{%
  \flxFont{subsubsection}\@dottedtocline{3}{7.0em}{4.1em}}
\newcounter{paragraph}[subsubsection]
\renewcommand{\theparagraph}{\thesubsubsection.\@arabic\c@paragraph}
\newcommand{\paragraph}{%
  \flx@startsection{paragraph}{4}%
    {\z@}%
    {3.25ex \@plus1ex \@minus.2ex}%
    {-1em}}
\newcommand*{\l@paragraph}{%
  \flxFont{paragraph}\@dottedtocline{4}{10em}{5em}}
\newcounter{subparagraph}[paragraph]
\renewcommand{\thesubparagraph}{\theparagraph.\@arabic\c@subparagraph}
\newcommand{\subparagraph}{%
  \flx@startsection{subparagraph}{5}%
    {\parindent}%
    {3.25ex \@plus1ex \@minus .2ex}%
    {-1em}}
\newcommand*{\l@subparagraph}{%
  \flxFont{subparagraph}\@dottedtocline{5}{12em}{6em}}
\if@twocolumn
  \setlength{\leftmargini}{2em}
\else
  \setlength{\leftmargini}{2.5em}
\fi
\leftmargin\leftmargini
\setlength{\leftmarginii}{2.2em}
\setlength{\leftmarginiii}{1.87em}
\setlength{\leftmarginiv}{1.7em}
\if@twocolumn
  \setlength{\leftmarginv}{.5em}
  \setlength{\leftmarginvi}{.5em}
\else
  \setlength{\leftmarginv}{1em}
  \setlength{\leftmarginvi}{1em}
\fi
\setlength{\labelsep}{.5em}
\setlength{\labelwidth}{\leftmargini}
\addtolength{\labelwidth}{-\labelsep}
\@beginparpenalty=-\@lowpenalty
\@endparpenalty=-\@lowpenalty
\@itempenalty=-\@lowpenalty
\renewcommand{\theenumi}{\@arabic\c@enumi}
\renewcommand{\theenumii}{\@alph\c@enumii}
\renewcommand{\theenumiii}{\@roman\c@enumiii}
\renewcommand{\theenumiv}{\@Alph\c@enumiv}
\newcommand{\labelenumi}{\theenumi.}
\newcommand{\labelenumii}{(\theenumii)}
\newcommand{\labelenumiii}{\theenumiii.}
\newcommand{\labelenumiv}{\theenumiv.}
\renewcommand{\p@enumii}{\theenumi}
\renewcommand{\p@enumiii}{\theenumi(\theenumii)}
\renewcommand{\p@enumiv}{\p@enumiii\theenumiii}
\newcommand{\labelitemi}{\textbullet}
\newcommand{\labelitemii}{\normalfont\bfseries \textendash}
\newcommand{\labelitemiii}{\textasteriskcentered}
\newcommand{\labelitemiv}{\textperiodcentered}
\newenvironment{description}
  {\list{}{%
     \labelwidth\z@
     \itemindent-\leftmargin
     \let\makelabel\descriptionlabel}}
  {\endlist}
\newcommand*{\descriptionlabel}[1]{%
  \hspace\labelsep\normalfont\bfseries #1}
\if@titlepage
  \newenvironment{abstract}
    {\titlepage
     \null\vfil
     \@beginparpenalty\@lowpenalty
     \begin{center}%
       \bfseries\abstractname
       \@endparpenalty\@M
     \end{center}}%
    {\par\vfil\null\endtitlepage}
\else
  \newenvironment{abstract}
    {\if@twocolumn
       \section*{\abstractname}%
     \else
       \small
       \begin{center}%
         {\bfseries\abstractname\vspace{-.5em}\vspace{\z@}}%
       \end{center}%
       \quotation
     \fi}
     {\if@twocolumn\else\endquotation\fi}
\fi
\newenvironment{verse}
  {\let\\\@centercr
   \list{}{%
     \itemsep=\z@
     \itemindent=-1.5em%
     \listparindent=\itemindent
     \rightmargin=\leftmargin
     \advance\leftmargin by 1.5em}%
   \item\relax}
  {\endlist}
\newenvironment{quotation}
  {\list{}{%
     \listparindent=1.5em%
     \itemindent=\listparindent
     \rightmargin=\leftmargin
     \parsep=\z@\@plus\p@}%
   \item\relax}
  {\endlist}
\newenvironment{quote}
  {\list{}{\rightmargin\leftmargin}\item\relax}
  {\endlist}
\newenvironment{titlepage}
  {\if@twocolumn
     \@restonecoltrue\onecolumn
   \else
     \@restonecolfalse\newpage
   \fi
   \thispagestyle{empty}%
   \setcounter{page}\@ne}%
  {\if@restonecol
     \twocolumn
   \else
     \newpage
   \fi
   \if@twoside\else
     \setcounter{page}\@ne
   \fi}
\newcommand{\appendix}{%
  \par
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \renewcommand\@chapapp{\appendixname}%
  \renewcommand\thechapter{\@Alph\c@chapter}}
\setlength{\arraycolsep}{5\p@}
\setlength{\tabcolsep}{6\p@}
\setlength{\arrayrulewidth}{.4\p@}
\setlength{\doublerulesep}{2\p@}
\setlength{\tabbingsep}{\labelsep}
\skip\@mpfootins=\skip\footins
\setlength{\fboxsep}{3\p@}
\setlength{\fboxrule}{.4\p@}
\if@flx@chapters
  \@addtoreset{equation}{chapter}
  \renewcommand{\theequation}{%
     \ifnum\c@chapter>\z@\thechapter.\fi\@arabic\c@equation}
\else
  \renewcommand{\theequation}{\@arabic\c@equation}
\fi
\if@flx@chapters
  \newcounter{figure}[chapter]
  \renewcommand{\thefigure}{%
     \ifnum\c@chapter>\z@\thechapter.\fi\@arabic\c@figure}
\else
  \newcounter{figure}
  \renewcommand{\thefigure}{\@arabic\c@figure}
\fi
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
  {\@float{figure}}
  {\end@float}
\newenvironment{figure*}
  {\@dblfloat{figure}}
  {\end@dblfloat}
\if@flx@chapters
  \newcounter{table}[chapter]
  \renewcommand{\thetable}{%
    \ifnum\c@chapter>\z@\thechapter.\fi\@arabic\c@table}
\else
  \newcounter{table}
  \renewcommand{\thetable}{\@arabic\c@table}
\fi
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
  {\@float{table}}
  {\end@float}
\newenvironment{table*}
  {\@dblfloat{table}}
  {\end@dblfloat}
\newlength{\abovecaptionskip}
\newlength{\belowcaptionskip}
\setlength{\abovecaptionskip}{10\p@}
\setlength{\belowcaptionskip}{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim\wd\@tempboxa>\hsize
    #1: #2\par
  \else
    \global\@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*{\cal}{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*{\mit}{\@fontswitch\relax\mathnormal}
\newcommand{\@pnumwidth}{1.55em}
\newcommand{\@tocrmarg}{2.55em}
\newcommand{\@dotsep}{4.5}
\setcounter{tocdepth}{2}
\newcommand{\tableofcontents}{%
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \chapter*{%
    \contentsname
    \@mkboth{\MakeUppercase\contentsname}%
            {\MakeUppercase\contentsname}}%
  \@starttoc{toc}%
  \if@restonecol
    \twocolumn
  \fi}
\newcommand{\listoffigures}{%
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \chapter*{%
    \listfigurename
    \@mkboth{\MakeUppercase\listfigurename}%
            {\MakeUppercase\listfigurename}}%
  \@starttoc{lof}%
  \if@restonecol
    \twocolumn
  \fi}
\newcommand*{\l@figure}{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \chapter*{%
    \listtablename
    \@mkboth{\MakeUppercase\listtablename}%
            {\MakeUppercase\listtablename}}%
  \@starttoc{lot}%
  \if@restonecol
    \twocolumn
  \fi}
\let\l@table\l@figure
\newdimen{\bibindent}
\setlength{\bibindent}{1.5em}
\newenvironment{thebibliography}[1]
  {\chapter*{%
     \bibname
     \@mkboth{\MakeUppercase\bibname}%
             {\MakeUppercase\bibname}}%
   \list{\@biblabel{\@arabic\c@enumiv}}{%
     \settowidth{\labelwidth}{\@biblabel{#1}}%
     \leftmargin=\labelwidth
     \advance\leftmargin\labelsep
     \@openbib@code
     \usecounter{enumiv}%
     \let\p@enumiv\@empty
     \renewcommand{\theenumiv}{\@arabic\c@enumiv}}%
   \sloppy
   \clubpenalty4000
   \@clubpenalty \clubpenalty
   \widowpenalty4000%
   \sfcode`\.\@m}
  {\def\@noitemerr{%
     \@latex@warning{Empty `thebibliography' environment}}%
   \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
  {\if@twocolumn
     \@restonecolfalse
   \else
     \@restonecoltrue
   \fi
   \columnseprule \z@
   \columnsep 35\p@
   \twocolumn[\@makeschapterhead{\indexname}]%
   \@mkboth{\MakeUppercase\indexname}%
           {\MakeUppercase\indexname}%
   \thispagestyle{plain}\parindent\z@
   \parskip\z@ \@plus .3\p@\relax
   \let\item\@idxitem}
  {\if@restonecol
     \onecolumn
   \else
     \clearpage
   \fi}
\newcommand{\@idxitem}{\par\hangindent 40\p@}
\newcommand{\subitem}{\@idxitem \hspace*{20\p@}}
\newcommand{\subsubitem}{\@idxitem \hspace*{30\p@}}
\newcommand{\indexspace}{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand{\footnoterule}{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\@addtoreset{footnote}{chapter}
\newcommand{\@makefntext}[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\newcommand{\contentsname}{Contents}
\newcommand{\listfigurename}{List of Figures}
\newcommand{\listtablename}{List of Tables}
\newcommand{\bibname}{Bibliography}
\newcommand{\indexname}{Index}
\newcommand{\figurename}{Figure}
\newcommand{\tablename}{Table}
\newcommand{\partname}{Part}
\newcommand{\chaptername}{Chapter}
\newcommand{\appendixname}{Appendix}
\newcommand{\abstractname}{Abstract}
\newcommand{\today}{}
\edef\today{%
  \ifcase\month
  \or January%
  \or February%
  \or March%
  \or April%
  \or May%
  \or June%
  \or July%
  \or August%
  \or September%
  \or October%
  \or November%
  \or December%
  \fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi
\flxDefAlignAll{\centering}
\flxDefAlign{paragraph}{\relax}
\flxDefAlign{subparagraph}{\relax}
\flxDefSizeAll{\normalsize}
\flxDefSize{part}{\LARGE}
\flxDefSize{chapter}{\LARGE}
\flxDefSize{section}{\Large}
\flxDefSize{subsection}{\large}
\flxDefFontAll{\itshape}
\flxDefFont{part}{\scshape}
\flxDefFont{chapter}{\scshape}
\flxDefNumber{part}{%
  \huge\bfseries\partname~\thepart\par
  \nobreak\vskip20pt\relax}
\flxDefNumber{chapter}{%
  \huge\bfseries\flxChapter\space\thechapter\par
  \nobreak\vskip20pt\relax}
\flxDefNumber{part}{\Huge---\thepart---\par\nobreak}
\flxDefNumber{chapter}{\Huge---\thechapter---\par\nobreak}
\endinput
%%
%% End of file `flex.cls'.
