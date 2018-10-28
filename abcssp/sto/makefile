#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = 2014sto

#--------------------------------------
# Files
#--------------------------------------
FILE_SRC  = main.tex                   \
            argmin.tex                 \
            graph.tex                  \
            lattice.tex                \
            metric.tex                 \
            ocs.tex                    \
            ocsrv.tex                  \
            ocsop.tex                  \
            order.tex                  \
            oms.tex                    \
            prbspace.tex               \
            relations.tex              \
            stats.tex                  \
            sums.tex                   
           

FILE_GRP  = ../common/graphics/watercraft/ghind_djg_blue.eps   \
            ../common/math/graphics/m2tol4phi.tex              \
            ../common/math/graphics/ftypes.tex                 \
            ../common/math/graphics/relfnctd.tex               \
            ../common/math/graphics/argmin/argmin_xx-2x2.tex   \
            ../common/math/graphics/metric/mpf_x1x2.tex        \
            ../common/math/graphics/ocs/ocs_fdie.tex           \
            ../common/math/graphics/ocs/ocs_rdie.tex           \
            ../common/math/graphics/ocs/ocs_wagon.tex          \
            ../common/math/graphics/oms/realline.tex           \
            ../common/math/graphics/oms/Zline.tex              \
            ../common/math/graphics/oms/Cplane.tex             \
            ../common/math/graphics/oms/oms_l2e3_235.tex       \
            ../common/math/graphics/oms/gsp_markov2.tex        \
            ../common/math/graphics/oms/pspinner.tex           \
            ../common/math/graphics/sto/archerytarget.tex      \
            ../common/math/graphics/sto/dqpsk.tex              \
            ../common/math/graphics/sto/fairdieXO6.tex         \
            ../common/math/graphics/sto/lcg7x1m9_XYZ.tex       \
            ../common/math/graphics/sto/lcg7x1m9_seqorder.tex  \
            ../common/math/graphics/sto/lcg7x1m9_dgraph.tex    \
            ../common/math/graphics/sto/realdieXRYZ.tex        \
            ../common/math/graphics/sto/rline_argminmax_11312.tex \
            ../common/math/graphics/sto/rline_argminmax_11312a.tex \
            ../common/math/graphics/ocsrv/ocsrv_spinner_xy.tex \
            ../common/math/graphics/ocs/ocs_wdie.tex           \
            ../common/math/graphics/sto/fairdieXL6.tex         \
            ../common/math/graphics/sto/realdieXL6.tex         \
            ../common/math/graphics/sto/dartstarget.tex        \
            ../common/math/graphics/sto/dna.tex                \
            ../common/math/graphics/sto/dmetricdnatol4.tex     \
            ../common/math/graphics/sto/dmetricdnatoc.tex      \
            ../common/math/graphics/sto/dnaXO4c.tex            \
            ../common/math/graphics/sto/ocsarchery.tex         \
            ../common/math/graphics/sto/ocsdarts.tex           \
            ../common/math/graphics/sto/pairdice.tex           \
            ../common/math/graphics/sto/pdice_xy.tex           \
            ../common/math/graphics/sto/ssea.tex               \
            ../common/math/graphics/sto/spinner6tol6.tex       \
            ../common/math/graphics/sto/spinnerXO6c.tex        \
            ../common/math/graphics/sto/wdie_xy.tex         


FILE_BIB  = ../common/bib/algebra.bib   \
            ../common/bib/analysis.bib  \
            ../common/bib/approx.bib    \
            ../common/bib/digcom.bib    \
            ../common/bib/found.bib     \
            ../common/bib/mathmisc.bib  \
            ../common/bib/order.bib     \
            ../common/bib/wavelets.bib  

FILE_STY  = ../common/sty/wavelets.sty \
            ../common/sty/math.sty     \
            ../common/sty/dan.sty      \
            ../common/sty/fonts.sty    \
            ../common/sty/packages.sty \
            ../common/sty/xsd.sty      \
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
  @echo start generating $(TARGET).pdf
  $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_STY) $(FILE_BIB) $(FILE_INP) $(FILE_GRP)
   echo % >> mainx.idx # make sure general   index file exists
   echo % >> mainx.rdx # make sure reference index file exists
   echo % >> mainx.sdx # make sure symbol    index file exists
   $(PRG_INDEX) -q -o mainx.ind mainx.idx
   $(PRG_INDEX) -q -o mainx.rnd mainx.rdx
   $(PRG_INDEX) -q -o mainx.snd mainx.sdx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_BIB)      mainx
   $(PRG_INDEX) -q -o mainx.ind mainx.idx
   $(PRG_INDEX) -q -o mainx.rnd mainx.rdx
   $(PRG_INDEX) -q -o mainx.snd mainx.sdx
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

nocite: 
   copy /Y main.tex mainx.tex
   echo % >> mainx.idx # make sure general   index file exists
   echo % >> mainx.rdx # make sure reference index file exists
   echo % >> mainx.sdx # make sure symbol    index file exists
   $(PRG_INDEX) -q -o mainx.ind mainx.idx
   $(PRG_INDEX) -q -o mainx.rnd mainx.rdx
   $(PRG_INDEX) -q -o mainx.snd mainx.sdx
   #$(PRG_BIB)      mainx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv
   copy /y mainx.pdf $(TARGET).pdf

mainx.tex: main.tex
  copy /Y main.tex mainx.tex

#--------------------------------------
# Commands
#--------------------------------------
gspdf:
  c:\p\gs902\bin\gswin32c -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(TARGET)_gs.pdf $(TARGET).pdf

x: mainx.tex
  lualatex mainx

clean:
   -del mainx.*
   -del *.aux
   -del ..\common\*.aux
   -del ..\common\wavelets\*.aux

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
  #xdvipdfmx  -s 1,7,8,11 -p a4 -q -o latp_ortho.pdf mainx


