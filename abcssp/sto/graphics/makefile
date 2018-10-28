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
PATHARGMIN = ../../../common/math/graphics/argmin/
PATHMATH   = ../../../common/math/graphics/
PATHMETRIC = ../../../common/math/graphics/metric/
PATHSTO    = ../../../common/math/graphics/sto/
PATHOMS    = ../../../common/math/graphics/oms/
PATHOCS    = ../../../common/math/graphics/ocs/
PATHOCSRV  = ../../../common/math/graphics/ocsrv/

#--------------------------------------
# Files
#--------------------------------------
#FILE_GRP  = ../common/watercraft/ghind_djg_blue.eps  \
#            ../../common/wavelets/morletwavelet.tex

GRP_MATH = \
           2xyfnct1.pdf                 \
           2xyfnct2.pdf                 \
           2xynotfnct1.pdf              \
           2xynotfnct2.pdf              \
           function_x4y4.pdf            \
           surjective_x4y4.pdf          \
           injective_x4y4.pdf           \
           bijective_x4y4.pdf           \
           m2tol4phi.pdf                \
           mpf_ascaled.pdf              \
           mpf_snowflake.pdf            \
           mpf_truncated.pdf            \
           mpf_bounded.pdf              \
           mpf_discrete.pdf             \
           mpf_x1x2.pdf                 \
           mpf_11x1.pdf                 \
           mpf_x11x1.pdf                \
           mpf_0121.pdf

GRP_ARGMIN = \
	argmin_xx-2x2.pdf


GRP_STO  = \
            archerytarget.pdf           \
            dartstarget.pdf             \
            dicetol6.pdf                \
            die_realline.pdf            \
            dmetricdnatoc.pdf           \
            dmetricdnatol4.pdf          \
            dna.pdf                     \
            dnatoc.pdf                  \
            dnatol4.pdf                 \
            dnaXO4c.pdf                 \
            dna_Rmax.pdf                \
            dqpsk.pdf                   \
            fairdieXL6.pdf              \
            fairdieXO6.pdf              \
            fairdieXO6c.pdf             \
            fairdieXO6cseg.pdf          \
            fairdieXRYZ.pdf             \
            gsp_markov2.pdf             \
            gspRZ.pdf                   \
            lcg7x1m9_dgraph.pdf         \
            lcg7x1m9_max.pdf            \
            lcg7x1m9_seqorder.pdf       \
            lcg7x1m9_XYZ.pdf            \
            m2tol4phi.pdf               \
            network.pdf                 \
            ocsarchery.pdf              \
            ocsdarts.pdf                \
            ocsrv_spinner_xy.pdf        \
            pairdice.pdf                \
            pdice_xy.pdf                \
            pEocsE_max.pdf              \
            pspinner_xy.pdf             \
            realdietoO6.pdf             \
            rdie_wxyz.pdf               \
            realdietoO6.pdf             \
            realdietoO6c.pdf            \
            realdieXL6.pdf              \
            realdieXR.pdf               \
            realdieXRYZ.pdf             \
            rline_argminmax_11312.pdf   \
            rline_argminmax_11312a.pdf  \
            spinner6tol6.pdf            \
            spinnerXO6c.pdf             \
            spinnerXO6Ycircle.pdf       \
            spinner_x2x.pdf             \
            ssea.pdf                    \
            udicetol6.pdf               \
            vending.pdf                 \
            wdietoO6c.pdf               \
            wdieXO6YO6c.pdf             \
            wdieXR.pdf                  \
            wdie_max.pdf                \
            wdie_xy.pdf                 \
            weighteddieXL6.pdf          \
            wspinnermax.pdf             \
            wspinner_Rmax.pdf

GRP_OCS    = \
            l4_abcd.pdf                 \
            ocs_fdie.pdf                \
            ocs_rdie.pdf                \
            ocs_spinner.pdf             \
            ocs_spinner_121141.pdf      \
            ocs_wagon.pdf               \
            ocs_wagon_121141.pdf        \
            ocs_wdie.pdf                \
            ocs_wdie_121141.pdf         \
            square_1234.pdf             \
            wring5.pdf                  \
            wring5short.pdf             \
            wring5shortd.pdf

GRP_OMS    = \
           Cplane.pdf                   \
           discretemetric6.pdf          \
           gsp_markov2.pdf              \
           oms_dna.pdf                  \
           oms_fdie.pdf                 \
           oms_l2e3_235.pdf             \
           oms_rdie.pdf                 \
           oms_spinner.pdf              \
           oms_wring5shortd.pdf         \
           pspinner.pdf                 \
           realline.pdf                 \
           ring6.pdf                    \
           wspinner.pdf                 \
           Zline.pdf


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
top: $(GRP_MATH) $(GRP_STO) $(GRP_OCS) $(GRP_OMS)  $(GRP_ARGMIN)  #$(FILE_STY)
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



{$(PATHMATH)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHMETRIC)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHARGMIN)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSTO)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHOCS)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHOCSRV)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHOMS)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

#--------------------------------------
# explicit build control
#--------------------------------------
#opTrn.pdf: ../../common/math/graphics/opTrn.tex 
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf
#
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

