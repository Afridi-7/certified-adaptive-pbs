"""
Build the Computer Methods and Programs in Biomedicine (Elsevier) manuscript paper/cmpb.tex from paper/main.tex.

The body (Introduction to Conclusion, tables, figures, appendices) is taken from main.tex unchanged, so both versions
carry the same numbers. Only the front matter (elsarticle class, structured abstract, keywords) and the end matter
(declarations required by Elsevier) differ. Highlights are written to paper/cmpb_highlights.txt (Elsevier asks for them
as a separate file). Run from the repository root or from paper/:  python paper/build_cmpb.py
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(HERE, "main.tex"), encoding="utf-8").read()

body = s[s.index("\\section{Introduction}"):s.index("\\paragraph{Data and code availability.}")]
data = s[s.index("\\paragraph{Data and code availability.}"):s.index("\\paragraph{Acknowledgements.}")]
data = data.replace("\\paragraph{Data and code availability.}", "\\section*{Data and code availability}").replace("Appendix~\\ref{", "\\ref{")
app = s[s.index("\\appendix"):s.index("\\end{document}")]
theorems = s[s.index("\\newtheorem{proposition}"):s.index("\\title{")]

ABSTRACT = r"""\begin{abstract}
\noindent\textbf{Background and Objective:} Automated smear examination increasingly decides which fields to inspect and when to stop. We asked what such a procedure can guarantee about unexamined material, and how many expert-verified cells are needed to certify whether the blast proportion lies above a diagnostic threshold.

\noindent\textbf{Methods:} We reimplemented the information-saturation halting rule of InfoDPP-PAC and compared it, at matched cost, with uniform sampling, a fixed field order and confidence-based stopping on pre-registered composite specimens built from 28,280 expert-labelled cells (KU-Optofil) with controlled rarity and dilution. We derived the minimax number of fields needed to certify the absence of rare evidence. For threshold decisions at 10\% and 20\% blasts we evaluated design-based, anytime-valid tests (ALPHA and weighted betting without replacement) with uniform and classifier-driven sampling, and tested their assumptions under acquisition shift (RV-PBS, MLL23), spatial clustering, unequal field sizes and reader disagreement in annotated bone-marrow regions (BoneMarrowWSI).

\noindent\textbf{Results:} The saturation rule's stopping point barely changed with rarity. With one critical field among 200 it missed 76\% of positive specimens, against 85\% for equal-size uniform sampling, and on real marrow tiles it did worse than uniform sampling in low-blast specimens. Without location information, any strategy must examine 95\% of the fields to certify absence at a 5\% error level. The valid tests kept false certification at or below 0.031 wherever their assumptions held. Their cost in verified cells was governed mainly by the margin to the threshold; near the threshold a 100-cell count was wrong in up to 37\% of specimens, whereas valid procedures were mostly undetermined. Classifier-driven sampling needed 0.28--0.73 times the verified cells of uniform sampling in distribution, 0.76--0.90 times on MLL23 and 1.77 times on RV-PBS, while validity held. Analysing the cells of sampled fields as independent draws exceeded the nominal error in 27 of 75 real specimens, and with a single reader as verifier the tests certified that reader's count rather than the consensus.

\noindent\textbf{Conclusions:} Information saturation does not certify unexamined material. A classifier can reduce verification effort, but validity depends on the sampling design, the unit of analysis and its weighting, and the reliability of verification.
\end{abstract}"""

KEYWORDS = r"""\begin{keyword}
adaptive sampling \sep anytime-valid inference \sep risk-limiting audit \sep peripheral blood smear \sep blast percentage \sep rare-event certification
\end{keyword}"""

HIGHLIGHTS = [
    "Saturation stopping lost its advantage over random sampling as cells were pooled",
    "Certifying absence of rare evidence needs near-census examination in the worst case",
    "Anytime-valid design-based tests certified blast thresholds at the nominal error",
    "Classifier-guided sampling saved verified cells in distribution, not under shift",
    "Validity required the sampled unit, cell-count weights and bounded verifier error",
]
assert all(len(h) <= 85 for h in HIGHLIGHTS), [len(h) for h in HIGHLIGHTS]

END = r"""
\section*{CRediT authorship contribution statement}
% To be confirmed by the authors.
\textbf{Qudrat Ullah:} Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Visualization, Writing -- original draft. \textbf{Bal\'azs Harangi:} Supervision, Writing -- review \& editing.

\section*{Declaration of competing interest}
% To be confirmed by the authors.
The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

\section*{Funding}
[Funding statement to be added by the authors.]

\section*{Acknowledgements}
[To be added.]

% Elsevier policy: any use of generative AI or AI-assisted tools in preparing the manuscript (beyond basic grammar and
% spelling checks) must be disclosed in a section titled "Declaration of Generative AI and AI-assisted technologies in
% the writing process", placed directly above the references. To be completed by the authors before submission.
"""

PREAMBLE = r"""\documentclass[preprint,12pt]{elsarticle}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{multirow}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{xurl}
\usepackage[colorlinks=true,linkcolor=blue!50!black,citecolor=blue!50!black,urlcolor=blue!50!black]{hyperref}
\usepackage{lineno}
\usepackage{adjustbox}
\usepackage[section]{placeins}
\emergencystretch=1.5em
\biboptions{numbers,sort&compress}
\journal{Computer Methods and Programs in Biomedicine}

""" + theorems + r"""
\begin{document}
\begin{frontmatter}
\title{Certified adaptive examination of peripheral blood smears: when information-gain stopping fails, and what can be guaranteed instead}
% Author list to be confirmed with the supervisor before submission.
\author[ud]{Qudrat Ullah\corref{cor}}
\ead{qudratafridi@mailbox.unideb.hu}
\author[ud]{Bal\'azs Harangi}
\cortext[cor]{Corresponding author.}
\affiliation[ud]{organization={Faculty of Informatics, University of Debrecen},
                 city={Debrecen}, country={Hungary}}

""" + ABSTRACT + "\n\n" + KEYWORDS + r"""
\end{frontmatter}
\linenumbers

"""

body = body.replace("\\begin{tabular}", "\\begin{adjustbox}{max width=\\linewidth}\n\\begin{tabular}").replace(
    "\\end{tabular}", "\\end{tabular}\n\\end{adjustbox}")
app = app.replace("\\begin{tabular}", "\\begin{adjustbox}{max width=\\linewidth}\n\\begin{tabular}").replace(
    "\\end{tabular}", "\\end{tabular}\n\\end{adjustbox}")
for name in ("body", "app"):
    t = locals()[name]
    t = t.replace("\\begin{table}[t]", "\\begin{table}[!htbp]").replace("\\begin{figure}[t]", "\\begin{figure}[!htbp]")
    t = t.replace("Appendix~\\ref{", "\\ref{")          # elsarticle's \ref to an appendix already reads "Appendix X"
    if name == "body": body = t
    else: app = t
out = PREAMBLE + body + data + END + "\n\\bibliographystyle{elsarticle-num}\n\\bibliography{refs}\n\n" + app + "\\end{document}\n"
out = out.replace("The author thanks", "The authors thank")
open(os.path.join(HERE, "cmpb.tex"), "w", encoding="utf-8").write(out)
open(os.path.join(HERE, "cmpb_highlights.txt"), "w", encoding="utf-8").write(
    "Highlights\n\n" + "\n".join("- " + h for h in HIGHLIGHTS) + "\n")
words = len(re.sub(r"\\[a-zA-Z]+|[{}]", " ", ABSTRACT[ABSTRACT.index("Background"):ABSTRACT.index("\\end{abstract}")]).split())
assert words <= 350, words
print("cmpb.tex written; abstract words:", words, "; highlights:", [len(h) for h in HIGHLIGHTS])
