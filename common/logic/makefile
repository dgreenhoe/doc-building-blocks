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
PATHARGMIN = ../../common/math/graphics/argmin/
PATHBALLS  = ../../common/math/graphics/balls/
PATHMATH   = ../../common/math/graphics/
PATHMETRIC = ../../common/math/graphics/metric/
PATHSTO    = ../../common/math/graphics/sto/
PATHSYMSEQ = ../../common/symseq/graphics/
PATHOMS    = ../../common/math/graphics/oms/
PATHOCS    = ../../common/math/graphics/ocs/
PATHOCSRV  = ../../common/math/graphics/ocsrv/

#--------------------------------------
# Files
#--------------------------------------
#FILE_GRP  = ../common/watercraft/ghind_djg_blue.eps  \
#            ../../common/wavelets/morletwavelet.tex

GRP_SHAPES = \
           gateand.pdf        \
           gateandABC.pdf     \
           gatenand.pdf       \
           gatenandABC.pdf    \
           gateor.pdf         \
           gateorABC.pdf      \
           gatenor.pdf        \
           gatenorABC.pdf     \
           gatexor.pdf        \
           gatexorABC.pdf     \
           gatenxor.pdf       \
           gatenxorABC.pdf


#GRP_ARGMIN    = \
#           ring6.pdf                    

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
top: $(GRP_SHAPES) 
	@echo Done ... sto/graphics pdfs are ready!
#  echo % >> mainx.idx # make sure general   index file exists

#--------------------------------------
# implicit build control
#--------------------------------------
#[../../common/math/graphics].tex.pdf:
#{[c:/dan/r/common/math/graphics/]}.tex.pdf:
#..\..\..\common\math\graphics\%.tex%.pdf:
#$&.pdf: ../../common/math/graphics/$&.tex
#[$(PATHSRC)].tex.pdf:
#{../}.tex.pdf:



.tex.pdf:
	$(PRG_TYPESET) $&.tex

#{$(PATHMETRIC)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHARGMIN)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHBALLS)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHSTO)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHSYMSEQ)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHOCS)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHOCSRV)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
#{$(PATHOMS)}.tex.pdf:
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf

#--------------------------------------
# explicit build control
#--------------------------------------

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

