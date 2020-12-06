#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------

#--------------------------------------
# TARGET
#--------------------------------------
#.RECIPEPREFIX=>

TARGET = abcssp

#--------------------------------------
# Files
#--------------------------------------
FILE_SRC  =                                 \
            mainx.tex                       \
            ../common/backmat.tex           \
            ../common/symbols.tex           \
            front/abstract.tex              \
            front/ackno.tex                 \
            front/front.tex                 \
            front/frontcover.tex            \
            front/preface.tex               \
            front/quotepagepair.tex         \
            front/symbols.tex               \
            srccode.tex                     \
            intro/chp_intro.tex             \
            intro/dft.tex                   \
            intro/fieldops.tex              \
            intro/filter.tex                \
            intro/order.tex                 \
            intro/overview.tex              \
            intro/relations.tex             \
            intro/sets.tex                  \
            intro/sequences.tex             \
            intro/stats.tex                 \
            intro/probability.tex           \
            chp_sto.tex                     \
            ../common/ocs.tex               \
            ../common/ocsrv.tex             \
            ../common/ocsop.tex             \
            ../common/oms.tex               \
            ssp/ocsseq.tex                  \
            ssp/sspdft.tex                  \
            ssp/chp_ssp.tex

FILE_GRP  = ../common/graphics/watercraft/ghind_djg_blue.eps  \
            ../common/graphics/cc/ccbyncnd40djg.pdf           \
            ../common/wavelets/graphics/morletwavelet.tex

#FILE_GRP  =   \
#            sto/graphics/rdie_wxyz.pdf

FILE_APP  =  \
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
            ../common/bib/biology.bib   \
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
            ../common/sty/xsd.sty       \
            ../common/sty/language.sty

#--------------------------------------
# Tools
#--------------------------------------
PRG_TYPESET = xelatex.exe
PRG_BIB     = bibtex.exe
PRG_INDEX   = makeindex.exe
PRG_XDV2PDF = xdvipdfmx.exe
#PRG_COPY    = cmd /c copy
PRG_COPY    = xcopy
#SHELL = cmd.exe

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

mainx.xdv: mainx.tex $(FILE_SRC) $(FILE_APP) $(FILE_STY) $(FILE_BIB)
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
	-del *.aux
	-del ..\common\*.aux
	-del ..\common\wavelets\*.aux
	-del *.ilg
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
	#xdvipdfmx  -s 7 -p a4 -q -o abstract_dan_20160625.pdf mainx
	#xdvipdfmx  -s 46,27,28 -p a4 -q -o dft_dan_20160625.pdf mainx
	#xdvipdfmx  -s 28,9,10 -p a4 -q -o dftex_dan_20160625.pdf mainx
	#xdvipdfmx  -s 7,11-16,20, -p a4 -q -o front_dan_20160626.pdf mainx
	xdvipdfmx  -s 1-142, -p a4 -q -o c:\dan\print\color.pdf mainx
	xdvipdfmx  -s 143-274, -p a4 -q -o c:\dan\print\bw.pdf mainx

#test:
#	$(PRG_COPY) /Y shell_morletwavelet.pdf a

