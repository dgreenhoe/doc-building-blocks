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
PATHPLOTS  = ../plots/

#--------------------------------------
# Files
#--------------------------------------
PLOTS =                                                         \
  fdie_5eed_51_seq.pdf                                          \
  fdie_5eed_16002_histo.pdf                                     \
  fdie_5eed_16002_auto.pdf                                      \
  rdie_5eed_16002_histo.pdf                                     \
  rdie_5eed_16002_auto.pdf                                      \
  spin_5eed_51_seq.pdf                                          \
  spin_5eed_16002_histo.pdf                                     \
  spin_5eed_16002_auto.pdf                                      \
  wrdie_5eed_51_seq.pdf                                         \
  wrdie_5eed_16002_histo.pdf                                    \
  wrdie_5eed_16002_auto.pdf                                     \
  wdie_5eed_16002_histo.pdf                                     \
  wdie_5eed_16002_auto.pdf                                      \
  wspin_5eed_16002_histo.pdf                                    \
  wspin_5eed_16002_auto.pdf                                     \
  dna_5eed_51_seq.pdf                                           \
  dna_5eed_16000_histo.pdf                                      \
  dna_5eed_16000_auto.pdf                                       \
  dna_sars_51_seq.pdf                                           \
  dna_sars_histo.pdf                                            \
  dna_sars_auto.pdf                                             \
  dna_ebola_51_seq.pdf                                          \
  dna_ebola_histo.pdf                                           \
  dna_ebola_auto.pdf                                            \
  dna_mpbacterium_51_seq.pdf                                    \
  dna_mpbacterium_histo.pdf                                     \
  dna_mpbacterium_auto.pdf                                      \
  dna_papaya1446_51_seq.pdf                                     \
  dna_papaya1446_histo.pdf                                      \
  dna_papaya1446_auto.pdf                                       \
  diedftR1_1525_1200m120.pdf                                    \
  diedftC1_1525_1200m120.pdf                                    \
  diedftR6_1525_1200m120.pdf                                    \
  diedftR6_1620_12000m1200.pdf                                  \
  diedftR6_1620_12000m120.pdf                                   \
  dnadftR1.pdf                                                  \
  dnadftC1.pdf                                                  \
  dnadftR4.pdf                                                  \
  dna_AF086833_ebola_dftR1.pdf                                  \
  dna_AF086833_ebola_dftC1.pdf                                  \
  dna_AF086833_ebola_dftR4.pdf                                  \
  wrdie_hp_1200_R1_rect50_euclid_seq.pdf                       \
  wrdie_hp_1200_R1_rect50_euclid_histo.pdf                     \
  wrdie_hp_1200_R1_rect50_euclid_auto.pdf                      \
  wrdie_hp_1200_R3_rect50_larc_seq.pdf                         \
  wrdie_hp_1200_R3_rect50_larc_histo.pdf                       \
  wrdie_hp_1200_R3_rect50_larc_auto.pdf                        \
  wrdie_hp_1200_R3_hann50_larc_seq.pdf                         \
  wrdie_hp_1200_R3_hann50_larc_histo.pdf                       \
  wrdie_hp_1200_R3_hann50_larc_auto.pdf                        \
  wspin_hp_1200_R1_rect50_euclid_seq.pdf                       \
  wspin_hp_1200_R1_rect50_euclid_histo.pdf                     \
  wspin_hp_1200_R1_rect50_euclid_auto.pdf                      \
  wspin_hp_1200_R1_hann50_euclid_seq.pdf                       \
  wspin_hp_1200_R1_hann50_euclid_histo.pdf                     \
  wspin_hp_1200_R1_hann50_euclid_auto.pdf                      \
  wspin_hp_1200_R2_rect50_larc_seq.pdf                         \
  wspin_hp_1200_R2_rect50_larc_histo.pdf                       \
  wspin_hp_1200_R2_rect50_larc_auto.pdf                        \
  wspin_hp_1200_R2_hann50_larc_seq.pdf                         \
  wspin_hp_1200_R2_hann50_larc_histo.pdf                       \
  wspin_hp_1200_R2_hann50_larc_auto.pdf                        \
  wspin_hp_1200_R2_hann50_euclid_histo.pdf                     \
  wspin_hp_1200_R2_hann50_euclid_auto.pdf                      \
  wrdie_hp_1200_R1_rect50_euclid_seq.pdf                       \
  wdie_hp_1200_R1_rect16_euclid_histo.pdf                      \
  wdie_hp_1200_R1_rect16_euclid_auto.pdf                       \
  wdie_hp_1200_R6_rect16_euclid_seq.pdf                        \
  wdie_hp_1200_R6_rect16_euclid_histo.pdf                      \
  wdie_hp_1200_R6_rect16_euclid_auto.pdf                       \
  wdie_hp_1200_R6_hann16_euclid_seq.pdf                        \
  wdie_hp_1200_R6_hann16_euclid_histo.pdf                      \
  wdie_hp_1200_R6_hann16_euclid_auto.pdf                       \
  wdie_hp_1200_R6_hann50_euclid_seq.pdf                        \
  wdie_hp_1200_R6_hann50_euclid_histo.pdf                      \
  wdie_hp_1200_R6_hann50_euclid_auto.pdf                       \
  wdie_hp_1200_R6_rect50_euclid_seq.pdf                        \
  wdie_hp_1200_R6_rect50_euclid_histo.pdf                      \
  wdie_hp_1200_R6_rect50_euclid_auto.pdf                       \
  rdie_lp_12000_R1_rect16_euclid_seq.pdf                        \
  rdie_lp_12000_R1_rect16_euclid_histo.pdf                      \
  rdie_lp_12000_R1_rect16_euclid_auto.pdf                       \
  rdie_lp_12000_R3_rect16_larc_seq.pdf                          \
  rdie_lp_12000_R3_rect16_larc_histo.pdf                        \
  rdie_lp_12000_R3_rect16_larc_auto.pdf                         \
  rdie_lp_12000_R3_hann16_larc_seq.pdf                          \
  rdie_lp_12000_R3_hann16_larc_histo.pdf                        \
  rdie_lp_12000_R3_hann16_larc_auto.pdf                         \
  rdie_lp_12000_R3_hann50_larc_histo.pdf                        \
  rdie_lp_12000_R3_hann50_larc_auto.pdf                         \
  rdie_lp_12000_R3_rect50_larc_histo.pdf                        \
  rdie_lp_12000_R3_rect50_larc_auto.pdf                         \
  spin_lp_12000_R1_rect16_euclid_seq.pdf                        \
  spin_lp_12000_R1_rect16_euclid_histo.pdf                      \
  spin_lp_12000_R1_rect16_euclid_auto.pdf                       \
  spin_lp_12000_R2_rect16_larc_seq.pdf                          \
  spin_lp_12000_R2_rect16_larc_histo.pdf                        \
  spin_lp_12000_R2_rect16_larc_auto.pdf                         \
  spin_lp_12000_R2_rect16_euclid_histo.pdf                      \
  spin_lp_12000_R2_rect16_euclid_auto.pdf                       \
  spin_lp_12000_R2_hann16_larc_seq.pdf                          \
  spin_lp_12000_R2_hann16_larc_histo.pdf                        \
  spin_lp_12000_R2_hann16_larc_auto.pdf                         \
  spin_lp_12000_R2_hann16_euclid_histo.pdf                      \
  spin_lp_12000_R2_hann16_euclid_auto.pdf                       \
  spin_lp_12000_R2_hann50_larc_histo.pdf                        \
  spin_lp_12000_R2_hann50_larc_auto.pdf                         \
  spin_lp_12000_R2_rect50_larc_histo.pdf                        \
  spin_lp_12000_R2_rect50_larc_auto.pdf                         \
  fdie_lp_12000_R6_rect16_euclid_seq.pdf                        \
  fdie_lp_12000_R6_rect16_euclid_histo.pdf                      \
  fdie_lp_12000_R6_rect16_euclid_auto.pdf                       \
  fdie_lp_12000_R6_hann16_euclid_seq.pdf                        \
  fdie_lp_12000_R6_hann16_euclid_histo.pdf                      \
  fdie_lp_12000_R6_hann16_euclid_auto.pdf                       \
  fdie_lp_12000_R6_hann50_euclid_histo.pdf                      \
  fdie_lp_12000_R6_hann50_euclid_auto.pdf                       \
  fdie_lp_12000_R6_rect50_euclid_histo.pdf                      \
  fdie_lp_12000_R6_rect50_euclid_auto.pdf                       \
  diehaarR1_12000m4000_h200_1050_D10.pdf                        \
  diehaarC1_12000m4000_h200_1050_D10.pdf                        \
  diehaarR6_12000m4000_h200_1050_D10.pdf                        \
  dnahaarR1_12000m4000_h200_1749_D5.pdf                         \
  dnahaarC1_12000m4000_h200_1749_D5.pdf                         \
  dnahaarR4_12000m4000_h200_1749_D5.pdf                         \
  dna_NC001416_phagelambdaR1b_hs1600.pdf                        \
  dna_NC001416_phagelambdaR1_h1600.pdf                          \
  dna_NC001416_phagelambdaC1_h1600.pdf                          \
  dna_NC001416_phagelambdaR4_h1600.pdf                          \
  dna_NC001416_phagelambdaR4_h4000.pdf    

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
top: $(PLOTS)
	@echo Done ... plot pdfs are ready!

#--------------------------------------
# implicit build control
#--------------------------------------
#[../../common/math/graphics].tex.pdf:
#{[c:/dan/r/common/math/graphics/]}.tex.pdf:
#..\..\..\common\math\graphics\%.tex%.pdf:
#$&.pdf: ../../common/math/graphics/$&.tex
#[$(PATHSRC)].tex.pdf:
#{../}.tex.pdf:

{$(PATHPLOTS)}.tex.pdf:
	$(PRG_TYPESET) $&.tex


#--------------------------------------
# explicit build control
#--------------------------------------
#r.pdf: $(PATHBALLS)/r.tex $(PATHBALLS)/r.dat $(PATHBALLS)/r2.dat $(PATHBALLS)/r3.dat shell_r.tex
#	$(PRG_TYPESET) shell_r.tex
#	copy /V/Y shell_r.pdf r.pdf
#	del shell_r.pdf

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

