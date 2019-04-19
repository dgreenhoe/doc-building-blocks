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

FILE_SRC  = ../common/backmata.tex      \
            main.tex                    \
            abstract.tex                \
            dtft.tex                    \
            estHx.tex                   \
            expectation.tex             \
            randseq.tex                 \
            xcordef.tex

FILE_BIB  = ../common/bib/algebra.bib   \
            ../common/bib/analysis.bib  \
            ../common/bib/approx.bib    \
            ../common/bib/found.bib     \
            ../common/bib/mathmisc.bib  \
            ../common/bib/order.bib     \
            ../common/bib/wavelets.bib  

FILE_STY  = ../common/cls/djgart.cls   \
            ../common/sty/article.sty  \
            ../common/sty/math.sty     \
            ../common/sty/ocs.sty      \
            ../common/sty/dan.sty      \
            ../common/sty/fonts.sty    \
            ../common/sty/packages.sty \
            ../common/sty/article.sty  \
            ../common/sty/wavelets.sty 

FILE_GRP  = ../common/graphics/watercraft/ghind_djg_blue.eps   \
            ../common/graphics/cc/ccby40djg.pdf                \
            ../common/math/graphics/2xyfnct1.tex               \
            ../common/math/graphics/2xyfnct2.tex               \
            ../common/math/graphics/2xynotfnct1.tex            \
            ../common/math/graphics/2xynotfnct2.tex            \
            ../common/math/graphics/function_x4y4.tex          \
            ../common/math/graphics/surjective_x4y4.tex        \
            ../common/math/graphics/injective_x4y4.tex         \
            ../common/math/graphics/bijective_x4y4.tex         \
            ../common/math/graphics/setstr/lattopxy.tex        \
            ../common/math/graphics/setstr/lattopxyz.tex       \
            ../common/math/graphics/setstr/latlattopxyz.tex    \
            ../common/math/graphics/metric/mpf_0121.tex        \
            ../common/math/graphics/metric/mpf_11x1.tex        \
            ../common/math/graphics/metric/mpf_ascaled.tex     \
            ../common/math/graphics/metric/mpf_bounded.tex     \
            ../common/math/graphics/metric/mpf_discrete.tex    \
            ../common/math/graphics/metric/mpf_snowflake.tex   \
            ../common/math/graphics/metric/mpf_truncated.tex   \
            ../common/math/graphics/metric/mpf_x11x1.tex       \
            ../common/math/graphics/metric/mpf_x1x2.tex        \
            ../common/math/graphics/sto/fairdieXO6.tex         \
            ../common/math/graphics/sto/fairdieXO6c.tex        \
            ../common/math/graphics/sto/fairdieXRYZ.tex        \
            ../common/math/graphics/ocs/ocs_fdie.tex           \
            ../common/math/graphics/ocs/ocs_rdie.tex           \
            ../common/math/graphics/ocs/ocs_wdie.tex           \
            ../common/math/graphics/ocs/ocs_spinner.tex        \
            ../common/math/graphics/sto/dqpsk.tex              \
            ../common/math/graphics/sto/fairdieXO6.tex         \
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
            ../common/symseq/graphics/earc_def.tex             \
            ../common/symseq/graphics/larc_def.tex             \
            ../common/symseq/graphics/larc_pq.tex              \
            ../common/symseq/graphics/larc_earc_compare.tex    \
            ../common/symseq/graphics/larc_euclid_compare.tex  \
            ../common/symseq/graphics/larc_metex.tex           \
            ../common/symseq/graphics/larc_hmg.tex             \
            ../common/symseq/graphics/larc_R2_balls.tex        \
            ../common/symseq/graphics/larc_ti.tex              \
            ../common/symseq/graphics/larc_trieq.tex           \
            ../common/symseq/graphics/limrp0.tex               \
            ../common/symseq/graphics/limrp0Rpq.tex   

#            ../common/symseq/graphics/spin_5eed_16ones_histo.tex         \


FILES_CPP = ..\2015ssp\c\larc.cpp     \
            ..\2015ssp\c\r2.cpp       \
            ..\2015ssp\c\r2op.cpp      

FILES_H   = ..\2015ssp\c\larc.h       \
            ..\2015ssp\c\r2.h         \
            ..\2015ssp\c\r2op.h        

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

gspdf:
  c:\p\gs902\bin\gswin32c -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(TARGET)gs.pdf $(TARGET).pdf

covers470x305:
  xelatex covers470x305
  rgbcheck.bat covers470x305
  #grep -i rgb covers470x305.ps
 

print:
  #xdvipdfmx  -s 1-4 -p a4 -q -o 2015sphx_intro.pdf mainx
  #xdvipdfmx  -s 7-9,61-63,65 -p a4 -q -o intbug.pdf mainx
  #xdvipdfmx  -s 1-2,18-26,29-61 -p a4 -q -o bases.pdf mainx
  #xdvipdfmx  -s 127-133 -p a4 -q -o hconditions.pdf mainx
  #xdvipdfmx  -s 1,7,8,11 -p a4 -q -o latp_ortho.pdf mainx


