#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------
#--------------------------------------
# Paths
#--------------------------------------
PATHMATH   = ../../common/math/graphics/
PATHSYSTEM = ../../common/math/graphics/system/
#--------------------------------------
# Files
#--------------------------------------
GRP_MATH = \
	opH_mnoise.pdf

FILE_STY  = \
	../../common/sty/wavelets.sty  \
	../../common/sty/math.sty      \
	../../common/sty/dan.sty       \
	../../common/sty/defaults.sty  \
	../../common/sty/fonts.sty     \
	../../common/sty/packages.sty  \
	../../common/sty/xsd.sty       \
	../../common/sty/language.sty

#--------------------------------------
# Programs
#--------------------------------------
PRG_TYPESET = xelatex.exe

#--------------------------------------
# master build control
#--------------------------------------
top: $(GRP_MATH) $(FILE_STY)
	@echo Done ... graphics pdf(s) are ready!

#--------------------------------------
# implicit build control
#--------------------------------------
{$(PATHMATH)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSYSTEM)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf
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

