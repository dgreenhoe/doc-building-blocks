#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = abcstat
#--------------------------------------
# Files
#--------------------------------------
FILE_SRC  =                       \
  front.tex                       \
  frontcover.tex                  \
  main.tex                        \
  src_code.tex                    \
  ../common/xsdfront.tex          \
  ../common/quotepagepair.tex     \
  ../common/binomial.tex          \
  ../common/calculus.tex          \
  ../common/com_channel.tex       \
  ../common/dtft.tex              \
  ../common/expectation.tex       \
  ../common/kl_basis.tex          \
  ../common/moment_estimation.tex \
  ../common/norm_minimization.tex \
  ../common/operator.tex          \
  ../common/pdfcdf.tex            \
  ../common/probspace.tex         \
  ../common/spline.tex            \
  ../common/randprocess.tex       \
  ../common/randseq.tex           \
  ../common/randseq_addnoise.tex  \
  ../common/randseq_dsp.tex       \
  ../common/randseq_eigen.tex     \
  ../common/randvar.tex           \
  ../common/info.tex              \
  ../common/estimate.tex          \
  ../common/netdet.tex            \
  ../common/setstrct.tex          \
  ../common/spectral.tex          \
  ../common/sums.tex              \
  ../common/measure.tex           \
  ../common/integrat.tex          \
  ../common/ortholat.tex          \
  ../common/mat_calc.tex          \
  ../common/subspace.tex          \
  ../common/pnseq.tex             \
  ../common/problatn.tex          \
  ../common/problatoc.tex         \
  ../common/symbols.tex           \
  ../common/systemid.tex          \
  ../common/backmat.tex           

FILE_APP  =  \
            mathtbls.tex                    \
            ../common/latc.tex

GRP_LAT = \
  graphics/lat4_m2_1xy0.pdf             \
  graphics/lat5_m3_xyz.pdf              \
  graphics/lat6_o6slash_123.pdf         \
  graphics/lat5_pentbslash_1xyz0.pdf    \
  graphics/lat5_l2onm2_1xyz0.pdf        \
  graphics/lat8_2e3_set235.pdf          \
  graphics/lat6_plat_P12345.pdf         \
  graphics/lat4_l4_1234.pdf             \
  graphics/lat6_n6_0abcp1.pdf           \
  graphics/lat5_n5_0acp1.pdf

FILE_GRP  = ../common/graphics/watercraft/ghind_djg_blue.eps  \
            ../common/graphics/cc/ccbyncnd40djg.pdf           \
            graphics/lat8_l2e3_prob_coin.pdf

FILE_BIB  = ../common/bib/analysis.bib  \
            ../common/bib/algebra.bib   \
            ../common/bib/approx.bib    \
            ../common/bib/mathmisc.bib  \
            ../common/bib/probstat.bib  \
            ../common/bib/wavelets.bib  

FILE_STY  = ../common/sty/wavelets.sty   \
            ../common/sty/colors_rgb.sty \
            ../common/sty/math.sty       \
            ../common/sty/dan.sty        \
            ../common/sty/fonts.sty      \
            ../common/sty/packages.sty   \
            ../common/sty/xsd.sty 

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

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_LAT) $(FILE_GRP) $(FILE_APP) $(FILE_STY) $(FILE_BIB) $(FILE_INP) makefile
   cd graphics
   make
   cd ..
   echo "\indexentry{{dummy index entry}}{1}" > mainx.idx # erase and/or initalize index file
   echo "\indexentry{{dummy index entry}}{1}" > xref.idx  # erase and/or initalize index file
   echo "\indexentry{{dummy index entry}}{1}" > xsym.idx  # erase and/or initalize index file
   $(PRG_INDEX) -q -o mainx.ind mainx.idx
   $(PRG_INDEX) -q -o  xref.ind  xref.idx
   $(PRG_INDEX) -q -o  xsym.ind  xsym.idx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_BIB)      mainx
   $(PRG_INDEX) -q -o mainx.ind mainx.idx
   $(PRG_INDEX) -q -o  xref.ind  xref.idx
   $(PRG_INDEX) -q -o  xsym.ind  xsym.idx
   $(PRG_TYPESET)  -no-pdf mainx
   $(PRG_TYPESET)  -no-pdf mainx

mainx.tex: main.tex
  copy /Y main.tex mainx.tex

#--------------------------------------
# Commands
#--------------------------------------
clean:
   -del mainx.*
   -del     *.aux
   -del     *.idx
   -del     *.ind

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
  xdvipdfmx  -s 57-63 -p a4 -q -o expectation.pdf mainx
  #xdvipdfmx  -s 1,5,19-22,34,38-39 -p a4 -q -o boxbug.pdf mainx
  #xdvipdfmx  -s 1-2,18-26,29-61 -p a4 -q -o bases.pdf mainx
  #xdvipdfmx  -s 127-133 -p a4 -q -o hconditions.pdf mainx


