#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------
#--------------------------------------
# TARGET
#--------------------------------------
#.RECIPEPREFIX=>
#--------------------------------------
# Paths
#--------------------------------------
PATHCC     = .

#--------------------------------------
# Files
#--------------------------------------
#FILE_GRP  = ../common/watercraft/ghind_djg_blue.eps  \
#            ../../common/wavelets/morletwavelet.tex

#  latlogics.pdf                     \
#  latnegations.pdf                  \
#  latortholat.pdf                   \

GRP_CC = \
  ccby40djg.pdf                     \
  ccbyncnd40djg.pdf    


FILE_STY  = ../../common/sty/wavelets.sty \
            ../../common/sty/math.sty     \
            ../../common/sty/dan.sty      \
            ../../common/sty/defaults.sty \
            ../../common/sty/fonts.sty    \
            ../../common/sty/packages.sty \
            ../../common/sty/xsd.sty       \
            ../../common/sty/language.sty

#--------------------------------------
# Programs
#--------------------------------------
PRG_TYPESET = xelatex.exe
#PRG_XDV2PDF = xdvipdfmx.exe

#--------------------------------------
# master build control
#--------------------------------------
top: $(GRP_CC) 
	@echo Done ... sto/graphics pdfs are ready!

#--------------------------------------
# implicit build control
#--------------------------------------
#[../../common/math/graphics].tex.pdf:
#{[c:/dan/r/common/math/graphics/]}.tex.pdf:
#..\..\..\common\math\graphics\%.tex%.pdf:
#$&.pdf: ../../common/math/graphics/$&.tex
#[$(PATHSRC)].tex.pdf:
#{../}.tex.pdf:

#{$(PATHCC)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf


#--------------------------------------
# explicit build control
#--------------------------------------
ccby40djg.pdf: ccby40djg.tex shell_ccby40djg.tex
	$(PRG_TYPESET) shell_ccby40djg.tex
	copy /V/Y shell_ccby40djg.pdf ccby40djg.pdf
	del shell_ccby40djg.pdf

ccbyncnd40djg.pdf: ccbyncnd40djg.tex shell_ccbyncnd40djg.tex
	$(PRG_TYPESET) shell_ccbyncnd40djg.tex
	copy /V/Y shell_ccbyncnd40djg.pdf ccbyncnd40djg.pdf
	del shell_ccbyncnd40djg.pdf

#latmn.pdf: ../../common/wavelets/graphics/latmn.tex 
#	$(PRG_TYPESET) shell_latmn.tex
#	copy /V/Y shell_latmn.pdf latmn.pdf
#	del shell_latmn.pdf

#--------------------------------------
# commands
#--------------------------------------
clean:
  del *.aux
  del *.log
  del *.out

scrub:
  del *.pdf
  make clean

