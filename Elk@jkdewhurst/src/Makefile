
include ../make.inc

# modules
SRC_modules = modmain.f90 $(SRC_MPI) modmpi.f90 modramdisk.f90 $(SRC_LIBXC) \
 modxcifc.f90 modfxcifc.f90 moddftu.f90 modrdm.f90 modphonon.f90 modtest.f90 \
 modrandom.f90 modpw.f90 modvars.f90 modtddft.f90 modgw.f90 modulr.f90 \
 modbog.f90 moddelf.f90 $(SRC_OMP) modomp.f90 $(SRC_MKL) mkl_init.f90 \
 $(SRC_OBLAS) oblas_init.f90 $(SRC_BLIS) blis_init.f90 $(SRC_W90S) modw90.f90

# Elk program
SRC_elk = elk.f90

# main subroutines and functions
SRC_main0 = \
 factn.f90 factr.f90 factn2.f90 r3mm.f90 r3mtm.f90 r3mmt.f90 r3vo.f90 gcd.f90

SRC_main1 = \
 symdmat.f90 axangsu2.f90 gradzvcln.f90 gentpmae.f90 zfcmtwr.f90 cfcmtwr.f90 \
 projsbf.f90 gradzf.f90 timesec.f90 zfmtftoc.f90 ztorfmt.f90 zfmtctof.f90 \
 eveqnfvz.f90 energynn.f90 gaunt.f90 potxcmt.f90 sdelta_mp.f90 stheta_mp.f90 \
 sdelta_fd.f90 stheta_fd.f90 writelsj.f90 eveqnzh.f90 gensmatk.f90 \
 gensfacgp.f90 rhomag.f90 wsplintp.f90 checkmt.f90 reciplat.f90 gradrf.f90 \
 genzfrm.f90 mixadapt.f90 zftrf.f90 sbessel.f90 gengclgq.f90 gencfun.f90 \
 genstress.f90 force.f90 genvbmatk.f90 forcek.f90 writeepsinv.f90 olprad.f90 \
 polark.f90 zfmtinp.f90 hflocal.f90 mixbroyden.f90 grad2rfmt.f90 clebgor.f90 \
 rfirctof.f90 effmass.f90 zpotcoul.f90 sfacmag.f90 polynm.f90 epsinv.f90 \
 rfinpc.f90 writeiad.f90 genapwfr.f90 gengclq.f90 addbfsm.f90 genylmv.f90 \
 olpfv.f90 hmlfv.f90 genexpmt.f90 rfinp.f90 sphcover.f90 sctovec.f90 \
 genppts.f90 eveqnzg.f90 fderiv.f90 genpmat.f90 genexpmat.f90 writeforces.f90 \
 r3cross.f90 match.f90 genspchi0.f90 genzrho.f90 potefield.f90 rhoinit.f90 \
 eveqnsv.f90 roteuler.f90 genrmesh.f90 splint.f90 wavefcr.f90 lopzflm.f90 \
 brzint.f90 rotdmat.f90 genylmg.f90 gentaucr.f90 readinput.f90 numlist.f90 \
 hermite.f90 genfspecies.f90 atpstep.f90 dmatls.f90 writeengy.f90 writechg.f90 \
 rbsht.f90 rfsht.f90 genlofr.f90 rotrfmt.f90 zftzf.f90 genpmatk.f90 \
 rvfcross.f90 writeinfo.f90 mossbauer.f90 initoep.f90 i3minv.f90 r3mv.f90 \
 genidxlo.f90 nuclei.f90 writeefg.f90 writeevsp.f90 gencore.f90 r3mtv.f90 \
 r3minv.f90 r3frac.f90 genapwlofr.f90 gridsize.f90 symrfir.f90 erf.f90 \
 geomplot.f90 ylmroty.f90 ssfext.f90 writestate.f90 sdelta_sq.f90 potxcir.f90 \
 potxc.f90 rzfmtinp.f90 rhoplot.f90 occupy.f90 symveca.f90 jprplot.f90 \
 genvmatk.f90 geomopt.f90 eveqnss.f90 hartfock.f90 potplot.f90 genwfpw.f90 \
 hmlxbse.f90 hmlxbsek.f90 hmldbse.f90 hmldbsek.f90 gradwfcr2.f90 oepresk.f90 \
 addlorbcnd.f90 energy.f90 wfplot.f90 oepvcl.f90 oepvclk.f90 polar.f90 \
 dbxcplot.f90 sbesseldm.f90 writehmlbse.f90 dielectric_bse.f90 vecplot.f90 \
 nesting.f90 writeevbse.f90 emdplot3d.f90 emdplot2d.f90 emdplot1d.f90 \
 rfhkintp.f90 writewfpw.f90 vblocal.f90 writefsm.f90 getwfpw.f90 plotpt1d.f90 \
 init1.f90 init2.f90 findsymcrys.f90 plotpt3d.f90 elnes.f90 fermisurfbxsf.f90 \
 genspecies.f90 writeexpmat.f90 writekpts.f90 wfcrplot.f90 writedos.f90 \
 energykncr.f90 ggamt_1.f90 ggair_1.f90 ggamt_sp_1.f90 ggair_sp_1.f90 \
 ggamt_2a.f90 ggair_2a.f90 ggamt_2b.f90 ggair_2b.f90 ggamt_sp_2a.f90 \
 ggair_sp_2a.f90 ggamt_sp_2b.f90 ggair_sp_2b.f90 writeemd.f90 reademd.f90 \
 emdplot.f90 init3.f90 init4.f90 bfieldfsm.f90 testcheck.f90 pade.f90 \
 wsplint.f90 gensocfr.f90 pades.f90 genscss.f90 rfcmtwr.f90 nonlinopt.f90 \
 writelat.f90 rhocore.f90 rzfinp.f90 readspecies.f90 findband.f90 exxengy.f90 \
 gentauk.f90 gradrfmt.f90 checkfsm.f90 genrlmv.f90 mae.f90 writegvecrf.f90 \
 zcfmtwr.f90 wfmtsv.f90 olpistl.f90 hmlistl.f90 olpalo.f90 hmlalo.f90 \
 olplolo.f90 hmllolo.f90 olpaa.f90 hmlaa.f90 gengclg.f90 rhomagsh.f90 \
 charge.f90 moment.f90 rfint0.f90 r3mdet.f90 gencfrc.f90 checkstop.f90 \
 genbs.f90 rfirsm.f90 rfmtsm.f90 rfint.f90 gradwf2.f90 closefiles.f90 \
 writelinen.f90 rfmtlm.f90 putpmat.f90 getpmat.f90 rfmtpack.f90 zfmtpack.f90 \
 putevalfv.f90 getevalfv.f90 wfmtfv.f90 putkmat.f90 getkmat.f90 putepsinv.f90 \
 getcfgq.f90 putevalsv.f90 getevalsv.f90 getevecsv.f90 putevecsv.f90 \
 genshtmat.f90 allatoms.f90 putevecfv.f90 getevecfv.f90 writesym.f90 \
 putoccsv.f90 getoccsv.f90 torque.f90 findqpt.f90 findkpt.f90 genevfsv.f90 \
 genwfsv.f90 sdelta.f90 stheta.f90 sortidx.f90 symmat.f90 radnucl.f90 \
 symrvf.f90 writeqpts.f90 rndevsv.f90 gengkvec.f90 mixlinear.f90 writesf.f90 \
 gengqf.f90 genkmat.f90 genolpq.f90 genjlgprmt.f90 rfinterp.f90 fsmooth.f90 \
 axangrot.f90 rotaxang.f90 wspline.f90 splined.f90 findswidth.f90 \
 findsymlat.f90 writepmat.f90 genstrain.f90 symmetry.f90 linengy.f90 \
 writestress.f90 eulerrot.f90 cftwfir.f90 rschrodint.f90 writegclq.f90 \
 writegeom.f90 writemom.f90 gendmat.f90 curlrvf.f90 genzvclmt.f90 genjtot.f90 \
 symvec.f90 checkwrite.f90 genjprk.f90 proj2d.f90 genidxbse.f90 symrvfir.f90 \
 genjpr.f90 eveqnfvr.f90 gentau.f90 genjlgpr.f90 writestrain.f90 \
 writefermi.f90 rfmtctof.f90 findsym.f90 zfsht.f90 rhonorm.f90 genvchi0.f90 \
 zfmtconj.f90 stheta_sq.f90 sdelta_lr.f90 stheta_lr.f90 findsymsite.f90 \
 genlmirep.f90 writespecies.f90 rhomagk.f90 symrvfmt.f90 wigner3jf.f90 \
 genhvec.f90 ylmrot.f90 wxcplot.f90 massnucl.f90 gendmatk.f90 rfplot.f90 \
 straingkq.f90 eveqnit.f90 oepmain.f90 rdirac.f90 bandstr.f90 sfacinit.f90 \
 findscq.f90 plot1d.f90 plot2d.f90 plot3d.f90 init0.f90 readstate.f90 \
 eveqnhf.f90 zbsht.f90 i3mtv.f90 spline.f90 gauntyry.f90 findnjcmax.f90 \
 gndstate.f90 sstask.f90 fermisurf.f90 spiralsc.f90 atom.f90 maginit.f90 \
 exxengyk.f90 xc_c_tb09.f90 moke.f90 hmlrad.f90 dos.f90 findprimcell.f90 \
 plotpt2d.f90 sfacrho.f90 bdipole.f90 dielectric.f90 elfplot.f90 vclcore.f90 \
 zfmtwr.f90 zvcldisp.f90 vecfbz.f90 lopzflmn.f90 zpotclmt.f90 \
 readfermi.f90 gensdmat.f90 symrfmt.f90 sphcrd.f90 potcoul.f90 genzrm.f90 \
 rminv.f90 rtozfmt.f90 genwfsvp.f90 genwfsvp_sp.f90 rdiracint.f90 rotzflm.f90 \
 zfinp.f90 latvstep.f90 gradzfmt.f90 eveqnfv.f90 cbsht.f90 cfsht.f90 \
 wfmtsv_sp.f90 wfirsv.f90 wfirsv_sp.f90 genwfsv_sp.f90 ggamt_4.f90 ggair_4.f90 \
 trzhmm.f90 gengvec.f90 wigner3j.f90 gengvc.f90 potnucl.f90 zfmtint.f90 \
 rfmtint.f90 genvsig.f90 rfmtinp.f90 findngkmax.f90 potks.f90 mtdmin.f90 \
 rfcopy.f90 rhomagv.f90 eveqn.f90 genffacgp.f90 unitary.f90 zmdet.f90 \
 rfmtftoc.f90 zminv.f90 symrf.f90 genvmat.f90 writeeval.f90 mixerifc.f90 \
 splinew.f90 nfftifc.f90 batchdv.f90 readevalsv.f90 readoccsv.f90 z2mm.f90 \
 z2mctm.f90 z2mmct.f90 zmctmu.f90 zmctm.f90

SRC_main = $(SRC_main0) $(SRC_main1)

SRC_BOG0 = \
 dmatuv.f90 dmatwx.f90 eveqnuv.f90 eveqnwxy.f90

SRC_BOG1 = \
 writeevaluv.f90 occupyuv.f90 putevaluv.f90 getevaluv.f90 putevecuv.f90 \
 getevecuv.f90 putevecwxy.f90 getevecwxy.f90 putevalwx.f90 getevalwx.f90 \
 aceplot.f90 writeevalwx.f90

SRC_BOG = $(SRC_BOG0) $(SRC_BOG1)

SRC_phonon0 = \
 dolpistl.f90 dhmlistl.f90 olpaloq.f90 hmlaloq.f90 dolpalo.f90 dhmlalo.f90 \
 dhmllolo.f90 olpaaq.f90 dolpaa.f90 dhmlaa.f90 hmlaaq.f90 initeph.f90 \
 dengyeph.f90 gndsteph.f90 eveqneph.f90 hmlepha.f90 hmlephb.f90 hmlephde.f90 \
 getephmkq.f90 ephdos.f90

SRC_phonon1 = \
 gradrhomt.f90 gendcfun.f90 gendvsig.f90 drhomagsh.f90 doccupy.f90 dhmlrad.f90 \
 gendwfsv.f90 dwfmtfv.f90 gengqvec.f90 dmatch.f90 writedvs.f90 readdyn.f90 \
 dyntask.f90 dpotks.f90 dpotxc.f90 deveqnfv.f90 dpotcoul.f90 drhomagk.f90 \
 gengvsmt.f90 dforce.f90 dforcek.f90 dynqtor.f90 dynevs.f90 phonon.f90 \
 dynrtoq.f90 mcmillan.f90 genephmat.f90 gengkqvec.f90 genmcph.f90 \
 putdevecfv.f90 getdevecfv.f90 putdevecsv.f90 getdevecsv.f90 phononsc.f90 \
 genscph.f90 phdisp.f90 dynsymapp.f90 dynsym.f90 dynev.f90 phdos.f90 \
 sumrule.f90 writephn.f90 dynfext.f90 writegamma.f90 readgamma.f90 \
 phlwidth.f90 alpha2f.f90 writelambda.f90 ephcouple.f90 eliashberg.f90 \
 phscdvs.f90 readdvs.f90 readalpha2f.f90 gendsocfr.f90 bornechg.f90 \
 bectask.f90 becfext.f90 putephmat.f90 getephmat.f90

SRC_phonon = $(SRC_phonon0) $(SRC_phonon1)

SRC_ULR0 = \
 zfirctof.f90 writeevalu.f90 writekpa.f90 genhdbulr.f90 vblocalu.f90 \
 chargeu.f90 addbfsmu.f90 momentu.f90 occupyulr.f90 potksu.f90 writevclr.f90 \
 genzvmatk.f90 genzvbmatk.f90 rhomagq.f90 writeengyu.f90 rfzfftq.f90

SRC_ULR1 = \
 rhocoreu.f90 eveqnulr.f90 rhomaguk.f90 genhmlu.f90 potxcu.f90 writeinfou.f90 \
 zfplot.f90 vclqinit.f90 writechgrmt.f90 initulr.f90 writemomrmt.f90 \
 writestulr.f90 readstulr.f90 potefieldu.f90 puthdbulr.f90 gethdbulr.f90 \
 potuinit.f90 readvclr.f90 genkpakq.f90 plotu1d.f90 plotu2d.f90 plotu3d.f90 \
 gndstulr.f90 plotulr.f90 potuplot.f90 potcoulu.f90 rhouplot.f90 maguplot.f90

SRC_ULR = $(SRC_ULR0) $(SRC_ULR1)

SRC_DFTU0 = \
 gendmatmt.f90 genveedu.f90 tm2todm.f90 tm3todm.f90 dmtotm3.f90

SRC_DFTU1 = \
 writedftu.f90 sbesseli.f90 shankeli.f90 genfdu.f90 genfdufr.f90 engyfdu.f90 \
 findlambda.f90 vmatmtsc.f90 genvmatmt.f90 vmatmtdu.f90 vmatmtftm.f90 \
 fyukawa.f90 fyukawa0.f90 readdmatmt.f90 writetm.f90 writetm3.f90 gendmftm.f90 \
 writeftm.f90 writetm3td.f90 genwkpr0.f90 tm3rtoz.f90

SRC_DFTU = $(SRC_DFTU0) $(SRC_DFTU1)

SRC_RDMFT0 = \
 rdmft.f90 rdmwriteengy.f90 rdmminc.f90 rdmvaryc.f90 rdmdedc.f90 rdmengyxc.f90 \
 rdmminn.f90 rdmvaryn.f90

SRC_RDMFT1 = \
 rdmdkdc.f90 rdmenergy.f90 rdmwritededn.f90 rdmdedn.f90 rdmdexcdn.f90 \
 genvcl1223.f90 getvcl1221.f90 getvcl1223.f90 writevcl1221.f90 \
 writevcl1223.f90 rdmdexcdc.f90 rdmdtsdn.f90 rdmentropy.f90 rdmeval.f90 \
 genvcl1221.f90

SRC_RDMFT = $(SRC_RDMFT0) $(SRC_RDMFT1)

SRC_GW0 = \
 gwdmat.f90 gwlocal.f90 genwgw.f90 ksinvert.f90 gtwsum.f90 minf_nm.f90 \
 acpole.f90

SRC_GW1 = \
 gwrhomag.f90 gwdmatk.f90 gwefermi.f90 dysonr.f90 gwsefm.f90 putgwsefm.f90 \
 getgwsefm.f90 writeigw.f90 polefit.f90 gwchgk.f90 gwsefmk.f90 gwspecf.f90 \
 writegwsf.f90 acgwse.f90 fmin_nm.f90 zfpole.f90 gwbandstr.f90 gwscrho.f90

SRC_GW = $(SRC_GW0) $(SRC_GW1)

SRC_TDDFT0 = \
 atptstep.f90 potkst.f90 tdbackup.f90 writetimes.f90 readtimes.f90 \
 gentimes.f90 readafieldt.f90 genjr.f90 tddft.f90 timestep.f90

SRC_TDDFT1 = \
 tddftlr.f90 genhmlt.f90 writetddos.f90 writeafpdt.f90 readafindt.f90 \
 genvfxc.f90 genfxcr.f90 tdtemp.f90 writetdlsj.f90 genspfxcr.f90 genvfxcg.f90 \
 tddftsplr.f90 genspfxcg.f90 genafieldt.f90 writetddft.f90 dielectric_tdrt.f90 \
 readjtot.f90 afindtstep.f90 genefieldt.f90 energytd.f90 engyknk.f90 \
 readforcet.f90 tdrestart.f90 readdatposc.f90 writedatposc.f90 rhostatic.f90 \
 jtotk.f90 rhosplot.f90 bornecdyn.f90 becforce.f90 tdinit.f90 readrhos.f90 \
 gengvnsmt.f90

SRC_TDDFT = $(SRC_TDDFT0) $(SRC_TDDFT1)

SRC_XC = \
 xc_pzca.f90 xc_pwca.f90 xc_pbe.f90 xc_am05.f90 xc_xalpha.f90 xc_wc06.f90 \
 x_wc06.f90 x_pbe.f90 c_pbe.f90 c_pbe_gcor.f90 xc_vbh.f90 fxc_pwca.f90 \
 k_tf.f90 k_tf_sp.f90 k_tf1.f90 k_tfvw.f90 k_tfvw_sp.f90 k_tfvw1.f90 k_vwlb.f90

SRC_W90 = \
 initw90.f90 writew90.f90 writew90win.f90 writew90eig.f90 setupw90.f90 \
 writew90mmn.f90 writew90spn.f90

SRC = $(SRC_modules) $(SRC_FFT) $(SRC_elk) $(SRC_main) $(SRC_phonon) \
 $(SRC_XC) $(SRC_DFTU) $(SRC_RDMFT) $(SRC_GW) $(SRC_TDDFT) \
 $(SRC_ULR) $(SRC_W90) $(SRC_BOG)

OBJ = $(SRC:.f90=.o)
EXE = elk

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
elk:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(F90_LIB) $(LIB_LIBXC) $(LIB_W90)

all:	elk

clean:
	rm -f *.o *.mod *~ ifc* *.gcno gmon.out *.aux *.dvi *.log \
         elk.pdf *.tex *.toc *.out $(EXE)

cleanall:
	$(MAKE) clean

doc:
	rm -f elk.aux elk.bbl elk.blg elk.log elk.pdf elk.tex elk.toc elk.lst
	ls $(SRC_modules) $(SRC_main) $(SRC_XC) $(SRC_phonon) $(SRC_DFTU) \
         $(SRC_RDMFT) $(SRC_TDDFT) > elk.lst
	./protex -F -s $(SRC_elk) $$(cat elk.lst) > elk.tex
	pdflatex elk;pdflatex elk;pdflatex elk

backup:
	tar -czf elk.tgz *.f90 libfftv eos spacegroup junk Makefile notes.txt \
         ../docs ../examples ../tests ../tests-libxc protex ../COPYING \
         ../README elk_silhouette.pdf depend gentest genalltests checklist.txt \
         rmspaces grepelk vimelk

VERSION = $$(grep "version(3)" modmain.f90 | sed 's/.*\[\([^]]*\)\].*/\1/g' | sed 's/,/./g')
NAME = elk-$(VERSION)
RELEASE = ../release/$(NAME)

version:
	echo $(VERSION)

release:
	rm -rf $(RELEASE)
	mkdir $(RELEASE)
	cd eos; $(MAKE) clean
	cd spacegroup; $(MAKE) clean
	cd ..; find . -type f -name *~ -exec rm -f {} \;
	cd ..; find . -type f -name gmon.out -exec rm -f {} \;
	cd ..; find . -type f -name RUNNING -exec rm -f {} \;
	mkdir $(RELEASE)/src
	cp -r eos $(RELEASE)/src
	cp -r spacegroup $(RELEASE)/src
	cp $(SRC) Makefile protex elk_silhouette.pdf $(RELEASE)/src
	cp $(SRC) rmspaces grepelk vimelk $(RELEASE)/src
	cp mpi_stub.f90 omp_stub.f90 mkl_stub.f90 oblas_stub.f90 blis_stub.f90 w90_stub.f90 $(RELEASE)/src
	cp libxcf90.f90 libxcifc.f90 libxcifc_stub.f90 $(RELEASE)/src
	cp zfftifc_fftw.f90 cfftifc_fftw.f90 zfftifc_mkl.f90 cfftifc_mkl.f90 $(RELEASE)/src
	cp -r ../examples $(RELEASE)
	cp -r ../tests $(RELEASE)
	cp -r ../tests-libxc $(RELEASE)
	mkdir $(RELEASE)/species
	cp ../species/*.in $(RELEASE)/species
	cp -r ../utilities $(RELEASE)
	cp ../COPYING $(RELEASE)
	cp ../README $(RELEASE)
	cp ../release_notes.txt $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../make.inc $(RELEASE)
	cp ../elk.sh $(RELEASE)
	$(MAKE) doc
	cp elk.pdf ../docs
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../../docs;$(MAKE) clean
	cp -r ../docs $(RELEASE)
	mkdir $(RELEASE)/release
	tar -C ../release -czf ../release/$(NAME).tgz $(NAME)
	cp ../docs/elk.pdf ../release
	cp ../docs/spacegroup.pdf ../release
	echo $(NAME)

lines:
	cat $(SRC) | wc -l

spaces:
	./rmspaces $(SRC)

vim:
	./vimelk

