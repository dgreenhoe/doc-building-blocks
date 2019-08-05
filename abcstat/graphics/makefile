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
PATHDSP      = ../../common/math/graphics/dsp/
PATHFRAME      = ../../common/math/graphics/frames/
PATHMATH     = ../../common/math/graphics/
PATHBSPLINES = ../../common/math/graphics/bsplines/
PATHOSET     = ../../common/math/graphics/oset
PATHOSETU    = ../../common/math/graphics/osetu
PATHLAT      = ../../common/math/graphics/lat
PATHLATU     = ../../common/math/graphics/latu
PATHNEG      = ../../common/math/graphics/neg
PATHLATOC    = ../../common/math/graphics/latoc
PATHNEG      = ../../common/math/graphics/neg
PATHLOGIC    = ../../common/math/graphics/
PATHFUZZY    = ../../common/math/graphics/fuzzy
PATHPROB     = ../../common/math/graphics/prob
PATHVENN     = ../../common/math/graphics/venn
PATHSETOPS   = ../../common/math/graphics/setops
PATHSETSTR   = ../../common/math/graphics/setstr
PATHSYSTEM   = ../../common/math/graphics/system
PATHWAVE     = ../../common/wavelets/graphics/

#--------------------------------------
# Files
#--------------------------------------
GRP_DSP  =                     \
  sinc.pdf                     \
  sinsin02.pdf                 \
  sinsin02hn.pdf               \
  sinsin02_pun.pdf             \
  sin_pi_t.pdf                 \
  haar0_sin_a.pdf              \
  haar0_sin_t.pdf              \
  haar1_sin_a.pdf              \
  haar1_sin_t.pdf              \
  haar2_sin_a.pdf              \
  haar2_sin_t.pdf              \
  dfI_order1_abc.pdf           \
  dfI_order1_lpf2pi3.pdf       \
  dfI_order2_156.pdf           \
  dfI_order2_fir.pdf           \
  dfI_order2_156.pdf           \
  dfII_order1_abc.pdf          \
  D4_pz.pdf                    \
  S4_pz.pdf                    \
  d4_phi_h.pdf                 \
  s4_phi_h.pdf                 \
  hanninglp50.pdf              \
  hpfwc.pdf                    \
  lpfwc.pdf                    \
  iir2n.pdf                    \
  iir2rphi.pdf                 \
  iir1pole.pdf                 \
  pz_allpass.pdf               \
  pz_minphase.pdf              \
  pz_minphase_inv.pdf          \
  pz_stable.pdf                \
  pz_pole00.pdf                \
  pz_realcoefs.pdf             \
  pz_realcoefs_11.pdf          \
  pz_unstable.pdf              \
  pz_unstable2.pdf             \
  pz_unall.pdf                 \
  rectangular_1N.pdf           \
  samplat.pdf                  \
  sin_pi_t_sampled.pdf         \
  sinclat.pdf                  \
  unitcircle.pdf               \
  vecres.pdf                   

GRP_FRAMES  =                  \
  frameR2N.pdf                 \
  frameR2pq.pdf                \
  mframe.pdf                   \
  mframe01.pdf                 \
  mframer01.pdf                \
  mframepm.pdf                 

GRP_MATH = \
  latr3subspaces.pdf           \
  r3subspaces_xy.pdf           \
  orderingsxyz.pdf             \
  pulseabc.pdf                 \
  pulseabcft.pdf               \
  triangleabc.pdf              \
  triangleabcft.pdf            \
  lat2xyz_antichains.pdf       \
  m3m3_antichain_abcpqr.pdf

GRP_PROB =                     \
  latestimation.pdf            \
  lat3_L3_1x0_-0-x-1_prob.pdf  \
  lat4_m2_psE_coin.pdf         \
  lat4_m2_boa_prob.pdf         \
  lat4_m2_prob_die.pdf         \
  lat6_o6_prob_10adbc.pdf      \
  lat8_l2e3_prob_coin.pdf      \
  lat8_l2e3_prob_die.pdf       \
  pdf_dicesum_bluered.pdf      \
  pdf_fairdie_blue.pdf         \
  pdf_fairdie_red.pdf

GRP_LATU = \
  lat1.pdf                     \
  lat2_l2.pdf                  \
  lat4_l4.pdf                  \
  lat6_m2m2.pdf                \
  lat6_o6bslash.pdf            \
  lat7_m2inm2_alt.pdf          \
  lat7_m2inm2.pdf              \
  m2tol4.pdf                   \
  m3ton5.pdf                   \
  n5ton5.pdf                   \
  lat8_l2e3.pdf                \
  lat4_l4.pdf                  \
  lat4_m2.pdf                  \
  lat5_l2onm2.pdf              \
  lat5_m2onl2.pdf

GRP_LAT = \
  lat4_m2_1xy0.pdf             \
  lat5_m3_xyz.pdf              \
  lat6_o6slash_123.pdf         \
  lat5_pentbslash_1xyz0.pdf    \
  lat5_l2onm2_1xyz0.pdf        \
  lat8_2e3_set235.pdf          \
  lat6_plat_P12345.pdf         \
  lat4_l4_1234.pdf             \
  lat6_n6_0abcp1.pdf           \
  lat5_n5_0acp1.pdf

GRP_NEG = \
  lat2_L2_10_-0-1.pdf                 \
  lat3_L3_1-1_x-x_0-0.pdf             \
  lat3_L3_1-0_x-1_0-x.pdf             \
  lat3_L3_1-x_x-0_0-1.pdf             \
  lat6_o6_1-0_a-b_b-a_c-d_d-c_0-1.pdf \
  lat3_L3_1x0_-0-1-x.pdf              \
  lat3_L3_1x0_-1x0.pdf                \
  lat3_L3_1-0-a_a_0-1.pdf             \
  lat3_L3_1-a-0_a-1_0.pdf             \
  lat3_L3_1x0_-0-x-1.pdf              \
  lat3_L3_1_a-0_0-1-a.pdf             \
  lat3_L3_1-0_a_0-1-x.pdf             \
  lat4_m2_1-0_a_b-b_0-1-a.pdf         \
  lat4_m2_a-b_b-a.pdf                 \
  lat4_bn4_a-a_b-b.pdf                \
  lat4_m2_1-0_a_b-b_0-1-a.pdf         \
  lat3_L3_1-0_a_0-1-x.pdf             \
  lat3_L3_1x0_-0-x-1.pdf              \
  lat5_m3_negdm.pdf                   \
  lat5_m3_negkl.pdf                   \
  lat6_o6_negdm.pdf                   \
  lat6_o6_negor.pdf                   \
  lat6_o6slash_nonantitone.pdf        \
  lat6_o6slash_kleene.pdf             \
  lat8_l2e3_setxyz.pdf                \
  lat8_l2e3_abc_dm.pdf                \
  lat8_l2e3_abc_ortho.pdf

GRP_OSET = \
  oset4_Nref.pdf                   \
  oset3_wedge.pdf                  \
  oset4_Nref.pdf                   \
  oset3_wedge.pdf                  \
  oset4_Nref.pdf                   \
  oset3_wedge.pdf                  \
  oset7_wedge_oplus_Nref.pdf       \
  oset4_Nref.pdf                   \
  oset3_wedge.pdf                  \
  oset12_wedge_times_Nref.pdf      \
  oset4_Nref.pdf                   \
  oset3_wedge.pdf                  \
  oset12_wedge_otimes_Nref.pdf

GRP_ORTHOLAT =               \
  lat6_o6_latoc_pq.pdf       \
  lat6_o6_xy.pdf             \
  lat8_o8_xyp.pdf            \
  lat8_l3l4l4l3_xypoc.pdf    \
  lat8_2e3m2_xyzoc.pdf       \
  lat10_n5n5oc_xyp.pdf       \
  lat10_m4m4oc_wxyz.pdf      \
  lat10_o8m2_xypqoc.pdf      \
  lat10_hsum_pxyz.pdf        \
  latocdil7_abxyzcd.pdf      \
  lat6_m4_wxyz.pdf           \
  lat8_m6_xyzoc.pdf          \
  lat1_01.pdf                \
  lat2_l2_01.pdf             \
  lat4_m2_latoc_p.pdf        \
  lat8_l2e3_latoc_pqr.pdf    \
  lat16_l2e4_latoc_pqrs.pdf  \
  lat32_l2e5_latoc_pqrst.pdf \
  r3subspaces_xyzxcyczc.pdf  \
  lat6_o6_pqab01.pdf         \
  lat6_primo_pqra01.pdf      \
  lat6_o6.pdf                \
  lat8_l2e3.pdf              \
  lat6_mo2_xy.pdf            \
  lat2_l2.pdf                \
  lat12_mo2x2.pdf

GRP_POUNITY =                \
  coscos11.pdf               \
  sinsin02_pun.pdf           \
  coscos11ft.pdf             \
  rcos.pdf                   \
  rcos_pun.pdf               \
  pounityz.pdf               \
  Fhw.pdf                    \
  d1_phi_h.pdf               \
  d1_psi_g.pdf               \
  tent_phi_h.pdf             \
  tent_psi_g.pdf             \
  bspline_phi_h.pdf          \
  bspline_psi_g.pdf          \
  Snn_bspline.pdf            \
  sinsin02.pdf               \
  sinsin02_pun.pdf           \
  sinsin02hn.pdf             \
  rsin.pdf                   \
  rsin_pun.pdf               \
  rsin02hn.pdf               
                             
GRP_SETOPS =                 \
  latsetops.pdf              \
  setop_0000.pdf             \
  setop_0001.pdf             \
  setop_0010.pdf             \
  setop_0011.pdf             \
  setop_0100.pdf             \
  setop_0101.pdf             \
  setop_0110.pdf             \
  setop_0111.pdf             \
  setop_1000.pdf             \
  setop_1001.pdf             \
  setop_1010.pdf             \
  setop_1011.pdf             \
  setop_1100.pdf             \
  setop_1101.pdf             \
  setop_1110.pdf             \
  setop_1111.pdf             
                             
GRP_BOOLEAN =                \
  lat2_l2_setx.pdf           \
  lat4_l2e2_setxy.pdf        \
  lat4_m2_setxy.pdf          \
  lat8_2e3_setxyz.pdf        
                             
GRP_SETSTR = \               
  lat5_m3_nontop.pdf         \
  lat5_n5_nontop.pdf         \
  lat5_l2onm2_top.pdf        \
  lat5_m2onl2_top.pdf        \
  lat5_l5_top.pdf            \
  lat6_primo_ss_xyz.pdf      \
  lat5_m3_ss_xyz.pdf         \
  latwavxyz.pdf              \
  latlattopxyz.pdf           \
  lattopxy.pdf               \
  lat2xyzdotted.pdf          \
  latlatalgwxyz.pdf          \
  latlatringxyz.pdf          \
  latlatpartitionsxyz.pdf    \
  latlatringxy.pdf           \
  latlatalgxyz.pdf           \
  latlatpartwxyz.pdf         \
  latlattopxyz_covers.pdf

GRP_SPLINE =                 \
  fstep.pdf                  \
  tent_pun.pdf               \
  n0.pdf                     \
  n1.pdf                     \
  n2.pdf                     \
  n3.pdf                     \
  n3x.pdf                    \
  n4.pdf                     \
  n5.pdf                     \
  n6.pdf                     \
  n7.pdf                     \
  b0_g.pdf                   \
  b1_g.pdf                   \
  b2_g.pdf                   \
  n0_h.pdf                   \
  n1_h.pdf                   \
  n2_h.pdf                   \
  n3_h.pdf                   \
  n4_h.pdf                   \
  n0punity.pdf               \
  n1punity.pdf               \
  n2punity.pdf               \
  n0dyadic.pdf               \
  n1dyadic.pdf               \
  n2dyadic.pdf               \
  pulse_pun.pdf              \
  Fpulse.pdf                 \
  tent012.pdf                \
  tent_pun.pdf               \
  Ftent012.pdf               \
  FN2.pdf                    \
  cencubs.pdf                \
  cencubs_pun.pdf  

GRP_SYSTEM =                 \
  sysH_xy.pdf                \
  sysHw_xy.pdf               \
  sysHz_xy.pdf               \
  sysH_cnoise.pdf            \
  sysH_xy0v-q.pdf            \
  sysH_xyuv--.pdf            \
  sysH_xy0v--.pdf            \
  sysH_xyu0--.pdf            \
  sysH_xyvq.pdf              \
  sysH_xyup.pdf              \
  sysH_xypw.pdf              \
  sysH_estH.pdf              \
  sysT_xyvq.pdf              \
  sysT_xyup.pdf              \
  sysT_xypw.pdf              \
  sys_xvy.pdf                \
  transferHs_p.pdf           \
  transferHs_ddp.pdf         \
  xGy_xHw.pdf                \
  xGy_xHw_addnoise.pdf       \
  xTy_xHw.pdf                \
  opH_cnoise.pdf             \
  opH_mnoise.pdf             \
  opT_cnoise.pdf             \
  opT_mnoise.pdf             \
  opT_opH_cnoise.pdf         \
  opT_opH_mnoise.pdf         \
  opT_estH_cnoise.pdf        \
  opT_estH_mnoise.pdf        \
  opT_estG_estH_mnoise.pdf   \
  opT_estG_estH_cnoise.pdf   \
  dual_addnoise.pdf          \
  sysH_mnoise.pdf

GRP_TOPOLOGY =               \
  lat5_m3_nontop.pdf         \
  lat5_m3_ss_xyz.pdf         \
  lat5_n5_nontop.pdf         \
  lat5_l2onm2_top.pdf        \
  lat5_m2onl2_top.pdf        \
  lat5_l5_top.pdf            \
  latwavxyz.pdf              \
  lattopxy.pdf               \
  latlattopxyz.pdf           \
  lat2xyzdotted.pdf          \
  lattopxyz.pdf              \
  butterfly_setA.pdf         \
  butterfly_clsA.pdf         \
  butterfly_intA.pdf         \
  butterfly_bndA.pdf         \
  butterfly_drvA.pdf         \
  cont_yx.pdf                \
  cont_yfy.pdf               \
  cont_ydy.pdf               


GRP_VENN =                   \
  venn_emptyset.pdf          \
  venn_cmpA.pdf              \
  venn_setdiff.pdf           \
  venn_cmpB.pdf              \
  venn_symdiff.pdf           \
  venn_intersect.pdf         \
  venn_union.pdf             \
  venn_universal.pdf         \
  venn_ABc.pdf               \
  venn_AB.pdf                \
  venn_AcB.pdf               \
  venn_AcBc.pdf

FILE_STY  = ../../common/sty/wavelets.sty \
            ../../common/sty/math.sty     \
            ../../common/sty/dan.sty      \
            ../../common/sty/defaults.sty \
            ../../common/sty/fonts.sty    \
            ../../common/sty/packages.sty \
            ../../common/sty/xsd.sty      \
            ../../common/sty/language.sty

#--------------------------------------
# Programs
#--------------------------------------
PRG_TYPESET = xelatex.exe
#PRG_XDV2PDF = xdvipdfmx.exe

#--------------------------------------
# master build control
#--------------------------------------
top: $(GRP_MATH) $(GRP_SETOPS) $(GRP_DSP) $(GRP_FRAMES) $(GRP_SETSTR) $(GRP_LATU) $(GRP_LAT) $(GRP_LATTICE) $(GRP_LATA) $(GRP_LATM) $(GRP_LATD) $(GRP_LATC) $(GRP_NEG) $(GRP_OSETU) $(GRP_OSET) $(GRP_BOOLEAN) $(GRP_ORTHOLAT) $(GRP_POUNITY) $(GRP_PROB) $(GRP_SPLINE) $(GRP_SYSTEM) $(GRP_TOPOLOGY)
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

{$(PATHMATH)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHFRAME)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHDSP)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHBSPLINES)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHLAT)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHLATU)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHNEG)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHLATOC)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

%.pdf : shell_%.tex
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHNEG)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHLOGIC)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHPROB)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHVENN)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSETOPS)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSETSTR)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHOSET)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHOSETU)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHSYSTEM)}.tex.pdf:
	$(PRG_TYPESET) shell_$&.tex
	copy /V/Y shell_$&.pdf $&.pdf
	del shell_$&.pdf

{$(PATHWAVE)}.tex.pdf:
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

