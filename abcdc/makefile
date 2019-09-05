#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = abcdc

#--------------------------------------
# Files
#--------------------------------------
FILE_SRC  = main.tex                          \
            symbols.tex                       \
            front.tex                         \
            frontcover.tex                    \
            mathtbls.tex                      \
            ../common/cdma.tex                \
            ../common/estimate_phase.tex      \
            ../common/distorted_frequency.tex \
            ../common/isi.tex                 \
            ../common/multipath.tex           \
            ../common/modulation_coding.tex   \
            ../common/com_channel.tex         \
            ../common/mod_sin.tex             \
            ../common/narrowband.tex          \
            ../common/projection_statistics.tex \
            ../common/em.tex

FILE_APP  = ../common/archive.tex      \
            ../common/estimate.tex     \
            ../common/harmonic.tex     \
            ../common/info.tex         \
            ../common/operator.tex     \
            ../common/mat_calc.tex

FILE_STY  =                            \
            ../common/sty/dan.sty      \
            ../common/sty/packages.sty \
            ../common/sty/math.sty     \
            ../common/sty/physics.sty  \
            ../common/sty/wavelets.sty \
            ../common/sty/language.sty

FILE_BIB  = ../common/bib/algebra.bib  \
            ../common/bib/analysis.bib \
            ../common/bib/em.bib       \
            ../common/bib/found.bib    \
            ../common/bib/mathhist.bib \
            ../common/bib/mathmisc.bib \
            ../common/bib/digcom.bib

#--------------------------------------
# Tools
#--------------------------------------
 PRG_TYPESET = xelatex.exe
 PRG_INDEX   = makeindex.exe
 PRG_BIB     = bibtex
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
  $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_APP) $(FILE_STY) $(FILE_BIB) makefile
   #cd graphics
   #make
   #cd ..
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
   -del *.idx
   -del *.aux

new:
   make clean
   make


pdf:
  make dcsd.pdf






