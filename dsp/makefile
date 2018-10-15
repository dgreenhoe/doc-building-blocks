#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
TARGET = dsp

#--------------------------------------
# Files
#--------------------------------------
FILE_SRC  =                                 \
            front.tex                       \
            ../common/symbols.tex           \
            ../common/backmat.tex           \
            ../common/calculus.tex          \
            ../common/dsp.tex               \
            ../common/dspcalculus.tex       \
            ../common/pwrspec.tex           \
            ../common/src.tex               \
            ../common/transop.tex           \
            ../common/zcalc1.tex            \
            ../common/zcalc2.tex            \
            ../common/zcalc2p.tex           \
            ../common/zchar.tex             \
            ../common/ztricks.tex           

FILE_GRP  =  \
  ../common/graphics/watercraft/ghind_djg_blue.eps \
  ../common/graphics/cc/ccbyncnd40djg.pdf           \
  graphics/baslat_cheby.pdf                        \
  graphics/baslat_cose.pdf                         \
  graphics/baslat_cosh.pdf                         \
  graphics/baslat_d1.pdf                           \
  graphics/baslat_d2.pdf                           \
  graphics/baslat_h3.pdf                           \
  graphics/baslat_haar.pdf                         \
  graphics/bspline_phi_h.pdf                       \
  graphics/bspline_psi_g.pdf                       \
  graphics/bspline_pz.pdf                          \
  graphics/cencubs.pdf                             \
  graphics/cencubs_pun.pdf                         \
  graphics/coscos11.pdf                            \
  graphics/coscos11ft.pdf                          \
  graphics/d10_phi_h.pdf                           \
  graphics/d10_psi_g.pdf                           \
  graphics/d11_phi_h.pdf                           \
  graphics/d11_psi_g.pdf                           \
  graphics/d12_phi_h.pdf                           \
  graphics/d12_psi_g.pdf                           \
  graphics/D12_pz.pdf                              \
  graphics/d16_phi_h.pdf                           \
  graphics/d16_psi_g.pdf                           \
  graphics/D16_pz.pdf                              \
  graphics/d1_g.pdf                                \
  graphics/d1_h.pdf                                \
  graphics/d1_phi.pdf                              \
  graphics/d1_phi_h.pdf                            \
  graphics/d1_psi.pdf                              \
  graphics/d1_psi_g.pdf                            \
  graphics/D1_pz.pdf                               \
  graphics/d2_phi_h.pdf                            \
  graphics/d2_psi_g.pdf                            \
  graphics/D2_pz.pdf                               \
  graphics/d3_phi_h.pdf                            \
  graphics/d3_psi_g.pdf                            \
  graphics/D3_pz.pdf                               \
  graphics/d4_phi_h.pdf                            \
  graphics/d4_psi_g.pdf                            \
  graphics/D4_pz.pdf                               \
  graphics/d5_phi_h.pdf                            \
  graphics/d5_psi_g.pdf                            \
  graphics/d6_phi_h.pdf                            \
  graphics/d6_psi_g.pdf                            \
  graphics/d7_phi_h.pdf                            \
  graphics/d7_psi_g.pdf                            \
  graphics/d8_phi_h.pdf                            \
  graphics/d8_psi_g.pdf                            \
  graphics/D8_pz.pdf                               \
  graphics/d9_phi_h.pdf                            \
  graphics/d9_psi_g.pdf                            \
  graphics/Dp_pz.pdf                               \
  graphics/DT_tentfx.pdf                           \
  graphics/Fhw.pdf                                 \
  graphics/Fpulse.pdf                              \
  graphics/Ftent012.pdf                            \
  graphics/haar0_sin_a.pdf                         \
  graphics/haar0_sin_t.pdf                         \
  graphics/haar1_sin_a.pdf                         \
  graphics/haar1_sin_t.pdf                         \
  graphics/haar2_sin_a.pdf                         \
  graphics/haar2_sin_t.pdf                         \
  graphics/hn1.pdf                                 \
  graphics/lat6_plat_xyzv10.pdf                    \
  graphics/latmn.pdf                               \
  graphics/latmra.pdf                              \
  graphics/latp_1235711.pdf                        \
  graphics/lattopxy.pdf                            \
  graphics/latwav.pdf                              \
  graphics/mexhatwavelet.pdf                       \
  graphics/morletwavelet.pdf                       \
  graphics/opDil.pdf                               \
  graphics/opTrn.pdf                               \
  graphics/pounityz.pdf                            \
  graphics/pulse.pdf                               \
  graphics/pulse03.pdf                             \
  graphics/pulse2.pdf                              \
  graphics/pulse_pun.pdf                           \
  graphics/rcos.pdf                                \
  graphics/rcos_pun.pdf                            \
  graphics/rsin.pdf                                \
  graphics/rsin02hn.pdf                            \
  graphics/rsin_pun.pdf                            \
  graphics/S10_pz.pdf                              \
  graphics/S11_pz.pdf                              \
  graphics/s12_phi_h.pdf                           \
  graphics/s12_psi_g.pdf                           \
  graphics/S12_pz.pdf                              \
  graphics/S13_pz.pdf                              \
  graphics/S14_pz.pdf                              \
  graphics/S15_pz.pdf                              \
  graphics/s16_phi_h.pdf                           \
  graphics/s16_psi_g.pdf                           \
  graphics/S16_pz.pdf                              \
  graphics/S17_pz.pdf                              \
  graphics/s4_phi_h.pdf                            \
  graphics/s4_psi_g.pdf                            \
  graphics/S4_pz.pdf                               \
  graphics/S5_pz.pdf                               \
  graphics/S6_pz.pdf                               \
  graphics/S7_pz.pdf                               \
  graphics/s8_phi_h.pdf                            \
  graphics/s8_psi_g.pdf                            \
  graphics/S8_pz.pdf                               \
  graphics/S9_pz.pdf                               \
  graphics/sinsin02.pdf                            \
  graphics/sinsin02hn.pdf                          \
  graphics/sinsin02_pun.pdf                        \
  graphics/sin_pi_t.pdf                            \
  graphics/TD_tentfx.pdf                           \
  graphics/tent012.pdf                             \
  graphics/tentfx.pdf                              \
  graphics/tent_phi_h.pdf                          \
  graphics/tent_psi_g.pdf                          \
  graphics/tent_pun.pdf                            \
  graphics/unitcircle.pdf                          \
  graphics/vecres.pdf
  

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
  @echo start generating $(TARGET).pdf
  $(PRG_XDV2PDF) $(OPT_XDV2PDF) -o mainx.pdf mainx.xdv

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_GRP) $(FILE_APP) $(FILE_STY) $(FILE_BIB) makefile
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


