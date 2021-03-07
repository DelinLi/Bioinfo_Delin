% This class file is provided for authors submitting to the
% Genetics Society of America's G3: Genes | Genomes |
% Genetics Journal
%
% Version: 1.2, 2017/02/10
%
% Produced by Overleaf (www.overleaf.com)
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{report}[05/04/2017, v1.6]
\AtEndOfClass{\RequirePackage[tracking=smallcaps]{microtype}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{extarticle}}
\ProcessOptions*
\LoadClass{extarticle}

% Fonts and language
\RequirePackage[utf8]{inputenc}
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{mathpazo}
\RequirePackage[scaled]{helvet}
\RequirePackage[T1]{fontenc}
\RequirePackage{url}
\RequirePackage[colorlinks=black, allcolors=black]{hyperref}
\RequirePackage{pifont}

\newcommand{\watermarkfont}{\color{black50}\normalfont\fontsize{10}{12}\selectfont}
\newcommand{\headerfont}{\normalfont\sffamily\fontsize{8}{10}\selectfont}
\newcommand{\titlefont}{\color{color3}\normalfont\sffamily\bfseries\fontsize{21}{23}\selectfont}
\newcommand{\datesfont}{\normalfont\sffamily\itshape\fontsize{8}{10}\selectfont}
\renewcommand\UrlFont{\color{black}\sffamily} %blue
\newcommand{\absfont}{\normalfont\sffamily\linespread{1.2}\fontsize{10}{12}\selectfont}
\newcommand{\keywordstitlefont}{\color{color4}\normalfont\sffamily\linespread{1.2}\fontsize{8}{10}\selectfont}
\newcommand{\keywordsfont}{\normalfont\sffamily\linespread{1.2}\fontsize{10}{12}\selectfont}
\newcommand{\copyrightfont}{\normalfont\sffamily\fontsize{7}{9}\selectfont}
\newcommand{\footerfont}{\normalfont\sffamily\fontsize{8}{10}\selectfont}

\RequirePackage[nodayofweek]{datetime}
\newdateformat{mydate}{\twodigit{\THEDAY}{ }\shortmonthname[\THEMONTH], \THEYEAR}
\newdateformat{dashdate}{\twodigit{\THEDAY}{-}\twodigit{\THEMONTH}{-}\THEYEAR}

\RequirePackage{enumitem} 
\setlist[itemize]{noitemsep}


% Author and affiliation
\RequirePackage{authblk}
\setlength{\affilsep}{0.1em}
\newcommand{\clientfont}{\color{black}\normalfont\sffamily\bfseries\fontsize{12}{14}\selectfont}
\newcommand{\intfont}{\color{black}\normalfont\sffamily\fontsize{10}{12}\selectfont}
%\renewcommand\Authfont{\color{color0}\normalfont\sffamily\bfseries\fontsize{9}{13}\selectfont}
%\renewcommand\Affilfont{\color{color0}\normalfont\sffamily\fontsize{9}{13}\selectfont}
\makeatletter
\renewcommand\AB@affilsepx{, \protect\Affilfont}
\makeatother
\renewcommand\Authands{ and }

%Frame
 \usepackage[framemethod=TikZ]{mdframed}
 \mdfdefinestyle{MyFrame}{%
    linecolor=gray!50!white, %red
    outerlinewidth=2pt,
    roundcorner=8pt,
    innertopmargin=4pt,%\baselineskip,
    innerbottommargin=4pt,%\baselineskip,
    innerrightmargin=4pt,
    innerleftmargin=4pt,
    backgroundcolor=gray!10!white}

% Choose article type
\newcommand*{\articletype}[1]{%
  \RequirePackage{styles/#1}}
  
% Options for element switching
\RequirePackage{xifthen}
\newboolean{displaycopyright}

% Copyright statement
\RequirePackage{textcomp} % For copyright symbol styling
%\makeatletter
%\newcommand{\copyrightstatement}{Copyright  \textcopyright{} \the\year{} by the Genetics Society of America}
%\makeatother

% Graphics, tables and other formatting
\RequirePackage{graphicx,xcolor}
\RequirePackage{colortbl}
\RequirePackage{booktabs}
\RequirePackage{tikz}
\RequirePackage{algorithm}
\RequirePackage[noend]{algpseudocode}
\RequirePackage{changepage}
\RequirePackage[left=54pt,%
                right=36pt,%
                top=48pt,%
                bottom=60pt,%
                headheight=40pt,%
                headsep=10pt,%
                a4paper]{geometry}% letterpaper,twoside]{geometry}
\RequirePackage[labelfont={bf,sf},%
                labelsep=space,%
                figurename=Figure,%
                singlelinecheck=off,%
                justification=RaggedRight]{caption}
\setlength{\columnsep}{24pt} % Distance between the two columns of text
\setlength{\parindent}{12pt} % Paragraph indent

%additional packages for table
\usepackage{siunitx,booktabs}
\sisetup{
  input-decimal-markers = .,input-ignore = {,},table-number-alignment = right,
  group-separator={,}, group-four-digits = true
}
\newcommand{\splitcell}[1]{{\begin{tabular}{@{}c@{}}#1\end{tabular}}}
\renewcommand{\arraystretch}{1.2} 
\usepackage{multirow} 

% Set document color scheme
\definecolor{black50}{gray}{0.5} % 50% black for hrules
\definecolor{color0}{RGB}{0,0,0} % Base color
\definecolor{color1}{RGB}{246, 245, 209} % Table column heads background
% color2 defined by the specific article type style file
\definecolor{color3}{RGB}{0,55,123} % Article Titles
\definecolor{color4}{RGB}{0,65,123}%{134,0,55} % Subheadings: A heads and “Abstract”
\definecolor{tabledingbat}{RGB}{184,108,108} % Table dingbat icon

% Bibliography
\usepackage{natbib} \bibpunct{(}{)}{;}{author-year}{}{,}
%\bibliographystyle{genetics}
%\addto{\captionsenglish}{\renewcommand{\refname}{Reference}}
\setlength{\bibsep}{0.0pt}



% Table style
\arrayrulecolor{color0} % Set the color of the horizontal rules
\renewcommand{\arraystretch}{1.5} % Increased line spacing
\captionsetup[table]{position=below} % Reduce spacing below caption
\colorlet{tableheadcolor}{color1} % Set table header colour
\newcommand{\header}{\rowcolor{tableheadcolor}} % create the command for the header
\captionsetup[table]{name={\hspace{0.8em}Table}} % Add the icon to the left of the table name
\captionsetup{justification=centering}
\usepackage{threeparttable} % for table notes May 9, 2017
% Algorithm style
\makeatletter
\newcommand\fs@coloruled{\def\@fs@cfont{\bfseries}\let\@fs@capt\floatc@coloruled
  \def\@fs@pre{\kern2pt}%
  \def\@fs@post{\kern2pt{\color{color2}\hrule}\relax}%
  \def\@fs@mid{\kern2pt{\color{color2}\hrule}\kern2pt}%
  \let\@fs@iftopcapt\iftrue}
\makeatother

\floatstyle{coloruled}
\restylefloat{algorithm}

% Headers and footers
\RequirePackage{fancyhdr}  % custom headers/footers
\RequirePackage{lastpage}  % Number of pages in the document
\pagestyle{fancy}          % Enables the custom headers/footers

\fancypagestyle{firststyle}{
   \fancyhead[L]{}
   \fancyhead[R]{}
    \fancyfoot[L]{\textbf{\small Bioinformatics Team}}%{\includegraphics[height=2em]{NGS_horiz}}
   \fancyfoot[R]{{\footerfont \hspace{1em}\mydate\today\hspace{1em}|\hspace{1em}\thepage}}
}

% Headers
\fancyhead[LO,RE,LE,RO]{} % Blank headers for non-title pages.
% Footers
\lfoot{}
\cfoot{}
\rfoot{}
\makeatletter
\fancyfoot[LO]{\textbf{\small Bioinformatics Team}} %\hspace{0.2em} \textit{et al.}
%\fancyfoot[LE]{{\footerfont \thepage\hspace{1em}|\hspace{1em}For Dr. \@correspondingauthor}} %\hspace{0.2em} \textit{et al.}
\fancyfoot[RO]{{\footerfont \hspace{1em}\mydate\today\hspace{1em}|\hspace{1em}\thepage\ }}
\makeatother
%\fancyfoot[LO,RE]{\includegraphics[height=1.8em]{Bioinformatics Team}}
\renewcommand{\headrulewidth}{0pt}% % No header rule
\renewcommand{\footrulewidth}{0pt}% % No footer rule

% Article type header box
\newcommand{\logomark}{
\begin{textblock}{180}(56,26) %200 36,36
\begin{flushleft}
\includegraphics[width=2.129in]{latex} %1.7129in
\end{flushleft}
\end{textblock}
\begin{textblock}{160}(400,36) %200 36,36
\begin{flushright}
\centering
%\textbf{\sffamily\bfseries\LARGE Report Date: \dashdate\today}
\end{flushright}
\end{textblock}
}

% Article type header box
\newcommand{\articletypemark}{
\begin{textblock}{200}(378,36)
\begin{flushright}
{\raggedleft \headerfont\color{black} \textls[96]{\articletypename}}
\end{flushright}
\end{textblock}
}

% Colored margin box
\newcommand{\marginmark}{
\begin{textblock}{72}(506,0) 
\noindent
\begin{tikzpicture}[x=1,y=1]
\fill[color=color2] (0,0) rectangle (72,36);
\end{tikzpicture}
\end{textblock}
}

%% 13 July 2016: make sure sections are unnumbered
\setcounter{secnumdepth}{0}
\setcounter{tocdepth}{0}

% Section/subsection/paragraph set-up
\RequirePackage[explicit]{titlesec}
\renewcommand{\thesubsection}{\Alph{subsection}}

\titleformat{\section}
  {\large\sffamily\bfseries\color{color4}}
  {\thesection.}
  {0.5em}
  {\MakeUppercase{#1}}
  []
\titleformat{name=\section,numberless}
  {\large\sffamily\bfseries\color{color4}}
  {}
  {0em}
  {\MakeUppercase{#1}}
  []
\titleformat{\subsection}
  {\sffamily\bfseries}
  {\thesubsection.}
  {0.5em}
  {#1}
  []
\titleformat{\subsubsection}[runin]
  {\sffamily\bfseries\itshape}
  {\thesubsubsection.}
  {0.5em}
  {#1}
  []
\titleformat{\paragraph}[runin]
  {\sffamily\small\bfseries}
  {}
  {0em}
  {#1}
\titlespacing*{\section}{0pc}{3ex \@plus4pt \@minus3pt}{5pt}
\titlespacing*{\subsection}{0pc}{2.5ex \@plus3pt \@minus2pt}{2pt}
\titlespacing*{\subsubsection}{0pc}{2ex \@plus2.5pt \@minus1.5pt}{2pt}
\titlespacing*{\paragraph}{0pc}{1.5ex \@plus2pt \@minus1pt}{12pt}

% Article meta data additional fields
\newcommand{\correspondingauthor}[1]{\def\@correspondingauthor{#1}}
\newcommand{\keywords}[1]{\def\@keywords{#1}}
\newcommand{\runningtitle}[1]{\def\@runningtitle{#1}}

\newcommand{\keywordname}{KEYWORDS} % Defines the keywords heading name

\def\xabstract{abstract}
\long\def\abstract#1\end#2{\def\two{#2}\ifx\two\xabstract
\long\gdef\theabstract{\ignorespaces#1}
\def\go{\end{abstract}}\else
\typeout{^^J^^J PLEASE DO NOT USE ANY \string\begin\space \string\end^^J
COMMANDS WITHIN ABSTRACT^^J^^J}#1\end{#2}
\gdef\theabstract{\vskip12pt BADLY FORMED ABSTRACT: PLEASE DO
NOT USE {\tt\string\begin...\string\end} COMMANDS WITHIN
THE ABSTRACT\vskip12pt}\let\go\relax\fi
\go}

% Set environment for abstract content and styling
\newcommand{\abscontent}{
\begin{minipage}[t]{505pt}
{\color{black50}\rule{\linewidth}{0.5pt}}
\vskip9pt%
{\absfont\bfseries\color{color4} SUMMARY} \absfont \theabstract%
\vskip1pt%
%{\color{black50}\rule{\linewidth}{0.5pt}}
\end{minipage}\hspace{10pt}
%\begin{minipage}[t]{72pt}
%\vskip18.3pt%
%{\keywordstitlefont\bfseries\keywordname} \\ \keywordsfont \@keywords%
%\end{minipage}
}%

% Custom title page
\renewcommand{\@maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{adjustwidth}{0pt}{24pt}
\begin{flushleft}
{%
\vskip60pt%
{\raggedright \titlefont \@title\par}%
\vskip11pt
{\raggedright \@author\par}
\vskip12pt%
}%
\end{flushleft}
\end{adjustwidth}
\egroup
{%
\abscontent
}%
\vskip25pt%
}%

% Other packages
\RequirePackage{enumitem} % For reducing bullet list item separation

\RequirePackage{tabularx}

\RequirePackage[flushmargin,ragged]{footmisc}

\newlength{\myFootnoteWidth}%
\newlength{\myFootnoteLabel}%
\setlength{\myFootnoteLabel}{0.5em}
\renewcommand{\@makefntext}[1]{%
  \setlength{\myFootnoteWidth}{\columnwidth}%
  \addtolength{\myFootnoteWidth}{-\myFootnoteLabel}%
  \noindent\makebox[\myFootnoteLabel][r]{\@makefnmark\ }%
  \parbox[t]{\myFootnoteWidth}{#1}
  }

% Commands to set the footnote for the front page.
\renewcommand{\footnoterule}{%
  \kern -3pt
  {\color{black50} \hrule width \columnwidth height 0.5pt}
  \kern 2.5pt
}

\newcommand\blfootnote[1]{%
  \begingroup
  \setlength{\myFootnoteLabel}{0em}
  \renewcommand\thefootnote{}\footnote{#1}%
  \addtocounter{footnote}{-1}%
  \addtocounter{Hfootnote}{-1}%
  \endgroup
}

\newcommand\correspondingauthoraffiliation[1]{%
  \begingroup
  \setlength{\myFootnoteLabel}{0em}
\renewcommand\thefootnote{}\footnote{\footnotelayout\copyrightfont\textsuperscript{1}#1}%
  \addtocounter{footnote}{-1}%
  \addtocounter{Hfootnote}{-1}%
  \endgroup
}

\newcommand{\firstpagefootnote}{\blfootnote{\copyrightfont %\ifthenelse{\boolean{displaycopyright}}%{
%\copyrightstatement}{}\\%
%Manuscript compiled: \today%
\begin{mdframed}[style=MyFrame]
{\large NOTE:} Results will be stored for up to 12 weeks after delivery of the report. If you would like your results purged sooner, please send your request via email to \underline{\textit{\href{mailto:delin.bio@xxx.com}{delin.bio@xxx.com}}}.
\end{mdframed}
}}

\newenvironment{tableminipage}[1]{\begin{minipage}{#1}\renewcommand\footnoterule{ \kern -1ex}%
\setlength{\myFootnoteLabel}{0.5em}%
}{\end{minipage}}

% Package to enable absolute positioning of figures
\RequirePackage[absolute,overlay]{textpos}
\setlength{\TPHorizModule}{1pt}
\setlength{\TPVertModule}{1pt}
\textblockorigin{0pt}{0pt} % start everything near the top-left corner

%default path
\graphicspath{ {figures/},{pictures/} }

%for text in margin
\usepackage[all]{background}
\SetBgContents{\textbf{\normalsize{Bioinformatics Team}}}
%\textbf{\normalsize{test\textsuperscript{\textregistered}:
\SetBgPosition{-1cm,-10.0cm}% Select location
\SetBgOpacity{1.0}% Select opacity
\SetBgAngle{90.0}% Select rotation of logo
\SetBgScale{1.4}% Select scale factor of logo
\SetBgColor{black!50!}%{green!35!blue!60!}
%for more color option, try http://latexcolor.com/

\usepackage{grffile} %to support file name with dots
\usepackage{seqsplit} %to support wrap long words
