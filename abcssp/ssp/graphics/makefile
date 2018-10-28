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

GRP_MATH = \
           2xyfnct1.pdf                 \
           2xyfnct2.pdf                 \
           2xynotfnct1.pdf              \
           2xynotfnct2.pdf              \
           function_x4y4.pdf            \
           surjective_x4y4.pdf          \
           injective_x4y4.pdf           \
           bijective_x4y4.pdf           \
           ae_ball(00_00_05).pdf        \
           larc_ball(00_00_05).pdf      \
           lattopxy.pdf                 \
           lat2xyzdotted.pdf            \
           lattopxyz.pdf                \
           latlattopxyz.pdf             \
           m2tol4.pdf                   \
           m3ton5.pdf                   \
           n5ton5.pdf                   \
           mpf_0121.pdf

GRP_ARGMIN = \
	argmin_xx-2x2.pdf


GRP_BALLS  =r.pdf                       \
            rdie_3d_larcball.pdf

GRP_METRIC =                 \
            mpf_0121.pdf     \
            mpf_11x1.pdf     \
            mpf_ascaled.pdf  \
            mpf_bounded.pdf  \
            mpf_discrete.pdf \
            mpf_snowflake.pdf\
            mpf_truncated.pdf\
            mpf_x11x1.pdf    \
            mpf_x1x2.pdf     

GRP_STO  = \
            dmetricdnatoc.pdf             \
            gspRZ.pdf                     \
            ocsrv_spinner_xy.pdf          \
            rdie_wxyz.pdf                 

GRP_SYMSEQ=                                                    \
            dna_5eed_16000_histo.pdf                           \
            dna_ebola_histo.pdf                                \
            dna_mpbacterium_histo.pdf                          \
            dna_mpbacterium_auto.pdf                           \
            dna_papaya1446_histo.pdf                           \
            dna_papaya1446_auto.pdf                            \
            dna_sars_histo.pdf                                 \
            fairdieXO6.pdf                                     \
            fairdieXO6c.pdf                                    \
            fairdieXRYZ.pdf                                    \
            hanninglp50.pdf                                    \
            hanninghp50.pdf                                    \
            earc_def.pdf                                       \
            fdie_5eed_16002_histo.pdf                          \
            fdie_5eed_16ones_r6ae.pdf                          \
            fdie_5eed_16ones_R1_histo.pdf                      \
            fdie_5eed_16ones_R6_histo.pdf                      \
            larc_def.pdf                                       \
            larc_pq.pdf                                        \
            larc_earc_compare.pdf                              \
            larc_euclid_compare.pdf                            \
            larc_hmg.pdf                                       \
            larc_metex.pdf                                     \
            larc_trieq.pdf                                     \
            larc_ti.pdf                                        \
            larc_R2_balls.pdf                                  \
            mca_R2_balls.pdf                                   \
            rdie_16ones.pdf                                    \
            rdie_16hp.pdf                                      \
            rdie_50hp.pdf                                      \
            rdie_5eed_16002_R1_rect_histo.pdf                  \
            rdie_5eed_16002_R3_rect_larc_histo.pdf             \
            rdie_5eed_16000_16ones_R3_auto.pdf                 \
            rdie_5eed_16ones.pdf                               \
            rdie_5eed_16ones_histo.pdf                         \
            rdie_5eed_16ones_histo_euclid.pdf                  \
            rdie_5eed_16ones_histo_larc.pdf                    \
            rdie_5eed_16ones_histo_R1.pdf                      \
            rdie_5eed_16ones_r3euclid.pdf                      \
            rdie_5eed_16ones_r3larc.pdf                        \
            rdie_5eed_16000_histo.pdf                          \
            rdie_r3met.pdf                                     \
            seq_lprect16.pdf                                   \
            seq_hprect16.pdf                                   \
            seq_shaar16.pdf                                    \
            seq_whaar16.pdf                                    \
            seq_shaar8.pdf                                     \
            seq_whaar8.pdf                                     \
            seq_whaar9.pdf                                     \
            spinner_r2.pdf                                     \
            spin_5eed_16002_histo.pdf                          \
            spin_5eed_16002_auto.pdf                           \
            spin_5eed_16000_histo.pdf                          \
            spin_5eed_16ones_histo_larc.pdf                    \
            spin_5eed_16ones_histo_R1.pdf                      \
            spinner_5eed.pdf                                   \
            spinner_5eed_16ones.pdf                            \
            spin_5eed_16ones_r2larc.pdf                        \
            test.pdf                                           \
            wspin_5eed_16002_histo.pdf                         \
            wspin_5eed_16002_auto.pdf                          \
            wspin_5eed_16000_histo.pdf                         \
            wdie_5eed_16002_histo.pdf                          \
            wdie_5eed_160000_down50_histo.pdf                  \
            wdie_5eed_160000_hanning_down50_R1_histo.pdf       \ 
            wdie_5eed_160000_hanning_down50_R3larc_histo.pdf   \ 
            wdie_5eed_160000_recthp_down50_R3larc_histo.pdf    \
            wdie_5eed_16002_auto.pdf                           \
            wdie_5eed_hanninghp_R1.pdf                         \
            wdie_5eed_hanninghp_R3_euclid.pdf                  \
            wdie_5eed_hanninghp_R3_larc.pdf                    \
            wdie_5eed_hanninghp_R1_histo.pdf                   \
            wdie_5eed_hanninghp_R3_euclid_histo.pdf            \
            wdie_5eed_hanninghp_R3_larc_histo.pdf              \
            wdie_5eed_recthp_R1.pdf                            \
            wdie_5eed_recthp_R3_euclid.pdf                     \
            wdie_5eed_recthp_R3_larc.pdf                       \
            wdie_5eed_recthp_R1_histo.pdf                      \
            wdie_5eed_recthp_R3_euclid_histo.pdf               \
            wdie_5eed_recthp_R3_larc_histo.pdf                 \
            wdie_5eed_5555_75_5_16hp.pdf                       \
            wdie_5eed_5555_75_5_16hp_r3larc.pdf                \
            wdie_5eed_5555_75_5_16hp_r3larc_16down.pdf 

GRP_OCS    =                                            \
            ocs_dna.pdf                                 \
            ocs_dnan.pdf                                \
            ocs_fdie.pdf                                \
            ocs_rdie.pdf                                \
            ocs_wdie.pdf                                \
            ocs_wrdie.pdf                               \
            ocs_spinner.pdf                             \
            wring5shortd.pdf

GRP_OMS    = \
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
top: $(GRP_MATH) $(GRP_BALLS) $(GRP_METRIC) $(GRP_STO) $(GRP_SYMSEQ) $(GRP_OCS) $(GRP_OMS)  $(GRP_ARGMIN)   #$(FILE_STY)
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

{$(PATHBALLS)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSTO)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSYMSEQ)}.tex.pdf:
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
r.pdf: $(PATHBALLS)/r.tex $(PATHBALLS)/r.dat $(PATHBALLS)/r2.dat $(PATHBALLS)/r3.dat shell_r.tex
	$(PRG_TYPESET) shell_r.tex
	copy /V/Y shell_r.pdf r.pdf
	del shell_r.pdf

larc_R2_balls.pdf: $(PATHSYMSEQ)/larc_R2_balls.tex $(PATHSYMSEQ)/larc_ball(0_0).dat $(PATHSYMSEQ)/larc_ball(00_05).dat $(PATHSYMSEQ)/larc_ball(1_1).dat  $(PATHSYMSEQ)/larc_ball(07_07).dat $(PATHSYMSEQ)/larc_ball(pio2_pio2).dat $(PATHSYMSEQ)/larc_ball(2_2).dat
	$(PRG_TYPESET) shell_larc_R2_balls.tex
	copy /V/Y shell_larc_R2_balls.pdf larc_R2_balls.pdf
	del shell_larc_R2_balls.pdf

mca_R2_balls.pdf: $(PATHSYMSEQ)/mca_R2_balls.tex $(PATHSYMSEQ)/larc_ball(0_0).dat $(PATHSYMSEQ)/mca_ball(00_05).dat $(PATHSYMSEQ)/larc_ball(1_1).dat  $(PATHSYMSEQ)/larc_ball(07_07).dat $(PATHSYMSEQ)/larc_ball(pio2_pio2).dat $(PATHSYMSEQ)/larc_ball(2_2).dat
	$(PRG_TYPESET) shell_mca_R2_balls.tex
	copy /V/Y shell_mca_R2_balls.pdf mca_R2_balls.pdf
	del shell_mca_R2_balls.pdf


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

