"""
Build the supplementary material of the CMPB article (paper/cmpb_supplement.tex) from the full preprint paper/main.tex.

The supplement carries the full methods and results, the deviations log, the proof of Theorem 1 and the pre-registration
audit trail, word for word from main.tex, so that the two never disagree. Figures and the tables that appear in the
article are removed and referenced by their number in the article. Run:  python paper/build_cmpb_supplement.py
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(HERE, "main.tex"), encoding="utf-8").read()


def cut(a, b):
    return s[s.index(a):s.index(b)]


methods = cut("\\section{Setting}", "\\section{Results}")
res_rq1 = cut("\\subsection{Information-saturation stopping (RQ1)}", "\\subsection{Certified threshold decisions (RQ3)}")
res_rq34 = cut("\\subsection{Certified threshold decisions (RQ3)}", "\\section{Discussion}")
proof = cut("\\section{Proof of Theorem~\\ref{thm:minimax}}", "\\section{Deviations from the protocol")
devs = cut("\\section{Deviations from the protocol", "\\section{Pre-registration audit trail}")
audit = cut("\\section{Pre-registration audit trail}", "\\end{document}")

# demote section levels
methods = methods.replace("\\section{", "\\subsection{")
devs = devs.replace("\\section{", "\\subsection{")
res_rq1 = res_rq1.replace("\\subsection{Information-saturation stopping (RQ1)}\\label{sec:partI}",
                          "\\section{Extended results: information-saturation stopping (RQ1)}\\label{sec:partI}")
proof = proof.replace("\\section{Proof of Theorem~\\ref{thm:minimax}}", "\\section{Proof of Theorem~\\ref{thm:minimax}}")
body = ("\\section{Extended methods}\\label{s:methods}\n" + methods + devs + res_rq1 +
        "\\section{Extended results: certification (RQ3) and assumptions (RQ4)}\\label{s:rq34}\n" + res_rq34 + proof + audit)

# drop figures, and the tables that appear in the article
body = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}\n?", "", body, flags=re.S)
for lab in ("tab:mscale", "tab:decision", "tab:bench", "tab:hyp"):
    body = re.sub(r"\\begin\{table\}(?:(?!\\begin\{table\}).)*?\\label\{" + lab + r"\}.*?\\end\{table\}\n?", "", body, flags=re.S)

ART_FIG = {"fig:partI": "1", "fig:partII": "2", "fig:partIII": "3", "fig:rq4": "4"}
ART_TAB = {"tab:mscale": "1", "tab:decision": "2", "tab:bench": "3", "tab:hyp": "4"}
for lab, n in ART_FIG.items():
    body = re.sub(r"Figure~\\ref\{" + lab + r"\}([a-c]?)", r"Fig.~" + n + r"\1 of the article", body)
    body = re.sub(r"Fig\.~\\ref\{" + lab + r"\}([a-c]?)", r"Fig.~" + n + r"\1 of the article", body)
for lab, n in ART_TAB.items():
    body = body.replace("Table~\\ref{" + lab + "}", "Table~" + n + " of the article")
body = body.replace("Appendix~\\ref{", "Section~\\ref{")
assert "\\ref{fig:" not in body and "\\ref{tab:" not in body, "unmapped figure or table reference"

PRE = r"""\documentclass[11pt]{article}
\usepackage[a4paper,margin=2.5cm]{geometry}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{xurl}
\usepackage[numbers,sort&compress]{natbib}
\usepackage[colorlinks=true,linkcolor=blue!50!black,citecolor=blue!50!black,urlcolor=blue!50!black]{hyperref}
\emergencystretch=1.5em
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thetable}{S\arabic{table}}
\newtheorem{proposition}{Proposition}
\newtheorem{theorem}{Theorem}
\newtheorem{corollary}{Corollary}
\theoremstyle{remark}
\newtheorem{remark}{Remark}
\theoremstyle{definition}
\newtheorem{assumption}{Assumption}
\renewcommand{\theproposition}{S\arabic{proposition}}
\renewcommand{\thetheorem}{S\arabic{theorem}}
\renewcommand{\thecorollary}{S\arabic{corollary}}
\renewcommand{\theremark}{S\arabic{remark}}
\renewcommand{\theassumption}{S\arabic{assumption}}
\newcommand{\Prb}{\mathbb{P}}
\newcommand{\E}{\mathbb{E}}

\title{Supplementary material for\\[2pt] \large ``Certified adaptive examination of peripheral blood smears: when information-gain stopping fails, and what can be guaranteed instead''}
\author{Qudrat Ullah, Bal\'azs Harangi}
\date{}

\begin{document}
\maketitle

\noindent This supplement gives the complete methods and results of every experiment summarised in the article, the deviations from the protocol and from the reimplemented method, the proof of Theorem~1 and the pre-registration audit trail. Figures 1--4 and Tables 1--4 are those of the article. Formal statements are numbered separately: Proposition~S1, Assumption~S1, Theorem~S1 and Remark~S1 are Proposition~1, Assumption~1, Theorem~1 and Remark~1 of the article; Proposition~S2 (unequal field sizes) is Proposition~3, and Proposition~S3 (ranked search) is Proposition~2. Section~S1.4 corresponds to Sections~2.4--2.8 of the article. Code, outputs and the audit trail are at \url{https://github.com/Afridi-7/certified-adaptive-pbs}.

"""
out = PRE + body + "\n\\bibliographystyle{unsrtnat}\n\\bibliography{refs}\n\\end{document}\n"
open(os.path.join(HERE, "cmpb_supplement.tex"), "w", encoding="utf-8").write(out)
print("cmpb_supplement.tex written:", len(out.split()), "tokens")
