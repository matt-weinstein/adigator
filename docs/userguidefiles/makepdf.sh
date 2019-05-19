#!/bin/bash

pdflatex ADiGatorUserGuide.tex
bibtex ADiGatorUserGuide
pdflatex ADiGatorUserGuide.tex
pdflatex ADiGatorUserGuide.tex
rm *.log *.toc *.out *.bbl *.blg *.aux
mv ADiGatorUserGuide.pdf ../