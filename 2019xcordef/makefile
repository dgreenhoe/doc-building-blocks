#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = 2019xcor

#--------------------------------------
# Files
#--------------------------------------
FILE_APP  = metrics.tex srccode.tex               

FILE_SRC  = \
	../common/backmata.tex      \
	main.tex                    \
	abstract.tex                \
	case_studies.tex            \
	case_addnoise.tex           \
	case_estHx.tex              \
	case_operator_mnoise.tex    \
	case_parallel_ops.tex       \
	case_dual_mnoise.tex        \
	case_dual_noise.tex         \
	case_systemid_nonlinear.tex \
	dsp.tex                     \
	expectation.tex             \
	intro.tex                   \
	normalg.tex                 \
	randseq.tex                 \
	xcordef.tex                 \
	xcor_xyreal.tex

FILE_BIB  = \
	../common/bib/algebra.bib   \
	../common/bib/analysis.bib  \
	../common/bib/approx.bib    \
	../common/bib/found.bib     \
	../common/bib/mathmisc.bib  \
	../common/bib/order.bib     \
	../common/bib/wavelets.bib  

FILE_STY  = \
	../common/cls/djgart.cls   \
	../common/sty/article.sty  \
	../common/sty/math.sty     \
	../common/sty/ocs.sty      \
	../common/sty/dan.sty      \
	../common/sty/fonts.sty    \
	../common/sty/packages.sty \
	../common/sty/article.sty  \
	../common/sty/wavelets.sty 

#--------------------------------------
# Tools
#--------------------------------------
 PRG_TYPESET = xelatex.exe
 PRG_BIB     = bibtex.exe
 PRG_INDEX   = makeindex.exe
 PRG_XDV2PDF = xdvipdfmx.exe

#--------------------------------------
# Options
#--------------------------------------
  OPT_XDV2PDF = -p a4

#--------------------------------------
# Build Control
#--------------------------------------
$(TARGET).pdf: mainx.pdf 
  copy /y mainx.pdf $(TARGET).pdf

mainx.pdf: mainx.xdv 
  @echo start generating $(TARGET).pdf
  @time /T
  $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv
  @time /T

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_STY) $(FILE_BIB) $(FILE_INP)  makefile
   cd graphics
   make
   cd ..
   echo "\indexentry{{dummy index entry}}{1}" > mainx.idx # erase and/or initalize index file
   echo "\indexentry{{dummy index entry}}{1}" > xref.idx  # erase and/or initalize index file
   echo "\indexentry{{dummy index entry}}{1}" > xsym.idx  # erase and/or initalize index file
   $(PRG_INDEX) mainx.idx
   $(PRG_INDEX)  xref.idx
   $(PRG_INDEX)  xsym.idx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_BIB)      mainx
   $(PRG_INDEX) mainx.idx
   $(PRG_INDEX)  xref.idx
   $(PRG_INDEX)  xsym.idx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_TYPESET)  -no-pdf mainx

mainx.tex: main.tex
  copy /Y main.tex mainx.tex

#--------------------------------------
# Commands
#--------------------------------------
clean:
   -del mainx.*
   -del *.aux
   -del *.idx
   -del *.ind

new:
   make clean
   make

xdv:
  del mainx.xdv
  make mainx.xdv

pdf:
  del mainx.xdv
  make

print:
  #xdvipdfmx  -s 1-4 -p a4 -q -o 2015sphx_intro.pdf mainx
