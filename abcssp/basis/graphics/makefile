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
PATHMATH = ../../../common/math/graphics/
PATHWAVE = ../../../common/wavelets/graphics/
#PATHSTO  = ../../../common/math/graphics/sto/

#--------------------------------------
# Files
#--------------------------------------
#FILE_GRP  = ../common/watercraft/ghind_djg_blue.eps  \
#            ../../../common/wavelets/morletwavelet.tex

FILE_GRP  = \
            baslat_cheby.pdf  \
            baslat_cose.pdf   \
            baslat_cosh.pdf   \
            baslat_haar.pdf   \
            latwav.pdf        \
            mexhatwavelet.pdf \
            morletwavelet.pdf \
            opTrn.pdf         \
            opDil.pdf         \
            latmn.pdf


FILE_STY  = ../../../common/sty/wavelets.sty \
            ../../../common/sty/math.sty     \
            ../../../common/sty/dan.sty      \
            ../../../common/sty/defaults.sty \
            ../../../common/sty/fonts.sty    \
            ../../../common/sty/packages.sty \
            ../../../common/sty/xsd.sty       \
            ../../../common/sty/language.sty

#--------------------------------------
# Programs
#--------------------------------------
PRG_TYPESET = xelatex.exe
#PRG_XDV2PDF = xdvipdfmx.exe

#--------------------------------------
# master build control
#--------------------------------------
top: $(FILE_GRP) $(FILE_STY)
	@echo Done ... basis/graphics pdfs are ready!
#  echo % >> mainx.idx # make sure general   index file exists

#--------------------------------------
# implicit build control
#--------------------------------------
#[../../../common/math/graphics].tex.pdf:
#{[c:/dan/r/common/math/graphics/]}.tex.pdf:
#..\..\..\common\math\graphics\%.tex%.pdf:
#$&.pdf: ../../../common/math/graphics/$&.tex
#[$(PATHSRC)].tex.pdf:
#{../}.tex.pdf:

{$(PATHMATH)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHWAVE)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

#{$(PATHSTO)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf

.pdf.eps:
	pdftops -eps $&.pdf $&.eps
	#pdftops -eps -noembt1 -noembtt $&.pdf $&.eps
	#pdftops -eps -noembt1 -noembtt -noembcidps -noembcidtt $&.pdf $&.eps

#--------------------------------------
# explicit build control
#--------------------------------------
#opTrn.pdf: ../../../common/math/graphics/opTrn.tex 
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#latmn.pdf: ../../../common/wavelets/graphics/latmn.tex 
#	$(PRG_TYPESET) shell_latmn.tex
#	copy /V/Y shell_latmn.pdf latmn.pdf
#	del shell_latmn.pdf

#--------------------------------------
# commands
#--------------------------------------
clean:
 #del *.eps
  del *.aux
  del *.log
  del *.out

make scrub:
  make clean
  del *.pdf

