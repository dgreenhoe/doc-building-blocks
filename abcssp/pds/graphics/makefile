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
PATHMATH   = ../../common/math/graphics/
PATHMETRIC = ../../common/math/graphics/metric/

#--------------------------------------
# Files
#--------------------------------------
GRP_MATH = ballinopen.pdf cont_yx.pdf cont_yfy.pdf cont_ydy.pdf

GRP_METRIC = trirel_sigmap.pdf  \
             mpf_ascaled.pdf    \
             mpf_snowflake.pdf  \
             mpf_truncated.pdf  \
             mpf_bounded.pdf    \
             mpf_discrete.pdf   \
             mpf_x1x2.pdf  

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
PRG_PDFTOEPS= pdftops -eps -pagecrop
#PRG_XDV2PDF = xdvipdfmx.exe

#--------------------------------------
# master build control
#--------------------------------------
top: $(GRP_MATH) $(GRP_METRIC)
	@echo Done ... sto/graphics pdfs are ready!
#  echo % >> mainx.idx # make sure general   index file exists

#--------------------------------------
# implicit build control
#--------------------------------------
#.pdf.eps:
#	$(PRG_PDFTOEPS) $&.pdf $&.eps

{$(PATHMATH)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHMETRIC)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

#--------------------------------------
# explicit build control
#--------------------------------------
#ballinopen.pdf: shell_ballinopen.tex $(PATHMATH)/ballinopen.tex
#	$(PRG_TYPESET) shell_ballinopen.tex
#	copy /V/Y shell_ballinopen.pdf ballinopen.pdf
#	$(PRG_PDFTOEPS) ballinopen.pdf ballinopen.eps
#	del shell_ballinopen.pdf
#
#trirel_sigmap.pdf: shell_trirel_sigmap.tex $(PATHMETRIC)/trirel_sigmap.tex
#	$(PRG_TYPESET) shell_trirel_sigmap.tex
#	copy /V/Y shell_trirel_sigmap.pdf trirel_sigmap.pdf
#	del shell_trirel_sigmap.pdf

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
