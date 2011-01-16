\NeedsTeXFormat{LaTeX2e}%
\ProvidesClass{myumnStatThesis}[2008/08/05 v1.0 University of Minnesota School of Statistics Doctoral Thesis]
\LoadClass[12pt, oneside]{book}

% v0.1  initial version
% v0.2  has headings and some on intro pages in small caps
% v0.3  2008/05/09 makes page nums bold and chapter numbers (eg Chapter # text) smaller
% v0.4  2008/05/14 revised title page font sizes, also spacing of title
% v0.41 2008/05/20 page nums bold on plain pages too
% v0.5  2008/06/13 add coadviser capability
% v0.51 2008/06/13 check for twoside when making headings
%        also move chaptermark before setting pagestyle to fancy
% v0.6  2008/08/05 allow front pages (acknowledgements, dedication, abstract) to move up/down
% v1.0  2009/03/10 extra features were separate but are now included

% ----- Useful AMS and other packages -----
\usepackage{amsbsy,amsmath,amsthm,amssymb,graphicx,verbatim}
\usepackage{subfigure}
\usepackage{multirow}
\usepackage[all]{xy}

\def\phantomsection{%
 \Hy@GlobalStepCount\Hy@linkcounter
 \xdef\@currentHref{section*.\the\Hy@linkcounter}%
 \Hy@raisedlink{\hyper@anchorstart{\@currentHref}\hyper@anchorend}%
}


% ----- Set Margins -----
\RequirePackage[total={6in,9in},top=1in,inner=1.5in,includeheadfoot,headheight=15pt]{geometry}

% ----- For Double Spacing, or any Variant -----
\RequirePackage{setspace}
\setstretch{1.5}
% Make the chapter headings look right,  given double spacing
%\let\origmakechapterhead\@makechapterhead
\def\origmakechapterhead#1{%
  \vspace*{50\p@}%
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
%        \huge\bfseries \@chapapp\space \thechapter
        \LARGE\bfseries \@chapapp\space \thechapter
%        \LARGE\textsc{\@chapapp\space \thechapter}
        \par\nobreak
        \vskip 20\p@
      \fi
    \fi
    \interlinepenalty\@M
    \Huge \bfseries #1\par\nobreak
    \vskip 40\p@
  }}

\let\origmakeschapterhead\@makeschapterhead
\renewcommand{\@makechapterhead}[1]{\begin{singlespace}\origmakechapterhead{#1}\end{singlespace}}
\renewcommand{\@makeschapterhead}[1]{\begin{singlespace}\origmakeschapterhead{#1}\end{singlespace}}

% ----- ADDED TITLING COMMANDS -----
\def\adviser#1{\gdef\@adviser{#1}}
\def\@adviser{\@latex@warning@no@line{No \noexpand\adviser given}}
\def\month#1{\gdef\@month{#1}}
\def\@month{\@latex@warning@no@line{No \noexpand\month given}}
\def\year#1{\gdef\@year{#1}}
\def\@year{\@latex@warning@no@line{No \noexpand\year given}}
\global\let\@coadviser\@empty
\def\coadviser#1{\gdef\@coadviser{#1}}

% ----- FOR BUILDING PRE MATTER PAGES -----
\newcommand{\clearasneeded}{\if@openright\cleardoublepage\else\clearpage\fi}
\newenvironment{prepage}
{\clearasneeded\thispagestyle{empty}\begin{singlespace}\begin{center}}
{\end{center}\end{singlespace}\clearasneeded}
%\newcommand\thesisEmph[1]{\MakeUppercase{#1}}
\newcommand\thesisEmph[1]{\textsc{\MakeLowercase{#1}}}


% ----- SIGNATURE PAGE -----
\newcommand\twoadvisersignatures[2]{%
\begin{minipage}{2.8in}\begin{center}%
{\large#1\par}%
\vskip 2pt\hrule\par\vskip 0.5em%
Name of Faculty Co-Adviser\par%
\vskip 5em\hrule\par\vskip 0.5em%
Signature of Faculty Co-Adviser\par%
\end{center}\end{minipage}%
\hskip 0.4in%
\begin{minipage}{2.8in}\begin{center}%
{\large#2\par}%
\vskip 2pt\hrule\par\vskip 0.5em%
Name of Faculty Co-Adviser\par%
\vskip 5em\hrule\par\vskip 0.5em%
Signature of Faculty Co-Adviser\par%
\end{center}\end{minipage}%
\vskip 5em\hrule\vskip 0.5em%
Date\par%
}

\newcommand\oneadvisersignature[1]{%
\begin{minipage}{5in}\begin{center}%
{\large#1\par}%
\vskip 2pt\hrule\par\vskip 0.5em%
Name of Faculty Adviser\par%
\vskip 5em\hrule\par\vskip 0.5em%
Signature of Faculty Adviser\par%
\vskip 5em\hrule\vskip 0.5em%
Date\par%
\end{center}\end{minipage}%
}

\newcommand\makesignaturepage{%
\begin{prepage}%
{\Large\thesisEmph{University of Minnesota}\par}%
\vfill%
{This is to certify that I have examined\\this copy of a doctoral dissertation by\par}%
\vskip 1em%
{\large\@author\par}%
\vskip 1em%
{and have found that it is complete and satisfactory in all respects,\\
and that any and all revisions required by the final \\
examining committee have been made.\par}
\vfill%
\ifx\@coadviser\@empty\oneadvisersignature{\@adviser}%
\else\twoadvisersignatures{\@adviser}{\@coadviser}\fi%
\vfill%
{\Large\thesisEmph{Graduate School}}%
\end{prepage}%
}

% ----- TITLE PAGE -----
\newcommand\maketitlepage{%
\begin{prepage}%
%\vspace*{\fill}
\vspace*{0.3in}
{\Large\begin{spacing}{1.4}\@title\end{spacing}\par}%
\vfill%
{\large\thesisEmph{A Dissertation}\\
\thesisEmph{Submitted to the Faculty of the Graduate School}\\
\thesisEmph{of the University of Minnesota}\\
\thesisEmph{by}\par}%
\vspace{3em}%
{\large\@author\par}%
\vfill%
{\large\thesisEmph{In Partial Fulfillment of the Requirements}\\
\thesisEmph{for the Degree of}\par
\thesisEmph{Doctor of Philosophy}\par}%
\vspace{3em}
\ifx\@coadviser\@empty{\large\@adviser, Adviser\par}%
\else{\large\@adviser\ and\\\@coadviser, Advisers\par}\fi%
\vspace{3em}
{\large\@month\ \@year}%
\end{prepage}%
}

% ----- COPYRIGHT PAGE -----
\newcommand\makecopyrightpage{%
\begin{prepage}%
\vspace*{\fill}%
\copyright\ \@author\ \@year%
\end{prepage}%
}

% ----- OTHER FRONT PAGES -----
\newenvironment{frontpage}[2][0in]
{\clearasneeded\thispagestyle{plain}\vspace*{#1}%
{\centering\large\textsc{#2}\par}%
\vspace{3em}\noindent\ignorespaces}
{\clearasneeded}

\newenvironment{acknowledgementspage}[1][0in]{\begin{frontpage}[#1]{Acknowledgements}}{\end{frontpage}}
\newenvironment{dedicationpage}[1][0in]{\begin{frontpage}[#1]{Dedication}}{\end{frontpage}}
\newenvironment{abstract}[1][0in]{\begin{frontpage}[#1]{Abstract}}{\end{frontpage}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----- The rest isn't strictly necessary, but I find it nice -----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----- Make it ragged bottom -----
% this prevents funny breaks in the middle of pages with little text
\raggedbottom

% ----- Add new Chapter mechanism -----
\newcommand{\chapterX}[1]{\chapter*{#1
        \markboth{\MakeUppercase{#1}}{\MakeUppercase{#1}}}%
        \addcontentsline{toc}{chapter}{#1}
}
  
% ----- Bibliography Stuff -----
% Don't give bibliography a chapter number
\RequirePackage{natbib}
\renewcommand\bibsection{\chapterX{\bibname}}
% Make Bibliography called "References"
\renewcommand\bibname{References}

% ----- Make headings look nicer -----
\RequirePackage{fancyhdr}
\fancyhf{}
\if@twoside
\fancyhead[LE]{\textbf{\thepage}}
\fancyhead[RE]{\nouppercase{\sc\leftmark}}
\fi
\fancyhead[RO]{\textbf{\thepage}}
\fancyhead[LO]{\nouppercase{\sc\rightmark}}
\renewcommand{\headrulewidth}{0pt}
% make chapters mark both sides to do the right thing 
% when a chapter has no sections
\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{\markboth
	{\MakeUppercase\@chapapp\ \thechapter. \ #1}
	{\MakeUppercase\@chapapp\ \thechapter. \ #1} }

% redefine the plain style to have bold page numbers
\fancypagestyle{plain}{% 
\fancyhf{} % clear all header and footer fields 
\fancyfoot[C]{\bfseries \thepage} % except the center 
\renewcommand{\headrulewidth}{0pt} 
\renewcommand{\footrulewidth}{0pt}}

% ----- Make blank pages really blank if using twoside -----
\let\origdoublepage\cleardoublepage
\newcommand{\clearemptydoublepage}{%
  \clearpage%
  {\pagestyle{empty}\origdoublepage}%
}
\let\cleardoublepage\clearemptydoublepage



% ---- Don't let figures float away from their own section ----
\RequirePackage[section,above,below]{placeins}



% ----- Example Environment -----
% textsquare from amsthm
%%\DeclareRobustCommand{\textsquare}{\begingroup\usefont{U}{msa}{m}{n}\thr@@\endgroup}
%%\RequirePackage[thmmarks,amsmath,fleqn]{ntheorem}

\newtheorem{theorem}{Theorem}[chapter]
\newtheorem{corollary}{Corollary}[chapter]
\newtheorem{lemma}{Lemma}[chapter]
\newtheorem{proposition}{Proposition}[chapter]
\newtheorem{claim}{Claim}

\theoremstyle{remark}
\newtheorem{example}{Example}[chapter]
\newtheorem{definition}{Definition}[chapter]
\newtheorem{remark}{Remark}[chapter]

\def\baro{\vskip  .2truecm\hfill \hrule height.5pt \vskip  .2truecm}
\def\barba{\vskip -.1truecm\hfill \hrule height.5pt \vskip .4truecm}

% ----- Other Commands -----
\newcommand\lla{
 \left\langle
}
\newcommand\rra{
 \right\rangle
}
