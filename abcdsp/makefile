#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = abcdsp

#--------------------------------------
# Files
#--------------------------------------
FILE_SRC  =                                 \
            front.tex                       \
            frontcover.tex                  \
            ../common/symbols.tex           \
            ../common/backmat.tex           \
            ../common/calculus.tex          \
            ../common/dsp.tex               \
            ../common/dspcalculus.tex       \
            ../common/laplace.tex           \
            ../common/pwspace.tex           \
            ../common/pwrspec.tex           \
            ../common/sampling.tex          \
            ../common/src.tex               \
            ../common/transop.tex           \
            ../common/zcalc1.tex            \
            ../common/zcalc2.tex            \
            ../common/zcalc2p.tex           \
            ../common/zchar.tex             

FILE_APP  =  \
            mathtbls.tex            \
            src_code.tex            \
            ../common/binomial.tex  \
            ../common/dsp.tex       \
            ../common/frames.tex    \
            ../common/interpo.tex   \
            ../common/operator.tex  \
            ../common/polynom.tex   \
            ../common/spline.tex    \
            ../common/harTrig.tex   \
            ../common/harPoly.tex   \
            ../common/harFour.tex  

FILE_BIB  = ../common/bib/analysis.bib  \
            ../common/bib/approx.bib    \
            ../common/bib/algebra.bib   \
            ../common/bib/mathhist.bib  \
            ../common/bib/mathmisc.bib  \
            ../common/bib/order.bib     \
            ../common/bib/wavelets.bib  

FILE_STY  = ../common/sty/wavelets.sty \
            ../common/sty/math.sty     \
            ../common/sty/dan.sty      \
            ../common/sty/defaults.sty \
            ../common/sty/fonts.sty    \
            ../common/sty/packages.sty \
           ..\common\sty\xsd.sty       \
            ../common/sty/language.sty

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
  copy /Y mainx.log tmp.log
  @echo start generating $(TARGET).pdf
  $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_APP) $(FILE_STY) $(FILE_BIB)
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

fast: 
   copy /Y main.tex mainx.tex
   echo % >> mainx.idx # make sure general   index file exists
   echo % >> mainx.rdx # make sure reference index file exists
   echo % >> mainx.sdx # make sure symbol    index file exists
   $(PRG_INDEX) -q -o mainx.ind mainx.idx
   $(PRG_INDEX) -q -o mainx.rnd mainx.rdx
   $(PRG_INDEX) -q -o mainx.snd mainx.sdx
   $(PRG_BIB)      mainx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv
   copy /y mainx.pdf $(TARGET).pdf

mainx.tex: main.tex
  copy /Y main.tex mainx.tex

#--------------------------------------
# Commands
#--------------------------------------
x: mainx.tex
  lualatex mainx

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

dvi:
  del mainx.xdv
  make mainx.xdv

pdf:
  del mainx.xdv
  make

covers470x305:
  xelatex covers470x305
  rgbcheck.bat covers470x305
  #grep -i rgb covers470x305.ps
 
covers470x318:
  xelatex covers470x318
  rgbcheck.bat covers470x318
  #grep -i rgb covers470x318.ps

covers533x305:
  xelatex covers533x305
  rgbcheck.bat covers533x305 
  #grep -i rgb covers533x305.ps

frontcover203x254:
  xelatex frontcover203x254
  rgbcheck.bat frontcover203x254

backcover203x254:
  xelatex backcover203x254
  rgbcheck.bat backcover203x254

spine203x254:
  xelatex spine203x254
  rgbcheck.bat spine203x254

print:
  #xdvipdfmx  -s 11-20 -p a4 -q -o complete_closed.pdf mainx
  #xdvipdfmx  -s 1,5,19-22,34,38-39 -p a4 -q -o boxbug.pdf mainx
  #xdvipdfmx  -s 1-2,18-26,29-61 -p a4 -q -o bases.pdf mainx
  #xdvipdfmx  -s 127-133 -p a4 -q -o hconditions.pdf mainx
  #xdvipdfmx  -s 391-394 -p a4 -q -o fwt.pdf mainx


