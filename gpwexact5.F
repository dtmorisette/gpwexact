      program gpwexact5
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Program "gpwexact" calculates a set of theoretical points
c     w(i),gpwd(i) for an n-type MOS capacitor, given parameters
c     defining the energy-dependent interface state density dit(E),
c     energy-dependent electron capture cross section sign(E), and
c     actual standard deviation of fixed charge sigqf.  The exact
c     us-vg relation is used in each elementary patch under the
c     gate to calculate the local surface potential as a function
c     of the local qf.  This avoids having to deal with a distorted
c     probability density function for surface potential caused
c     by the nonlinear us-vg relationship near flat band.  The
c     exact us-vg relationship includes the spreading effect of
c     charge in interface states qit(ufs) in each patch.
c
c     The w(i),gpwd(i) points produced by this program are then
c     regarded as experimental data, and are passed to a least-
c     squares fitting routine to extract the APPARENT interface
c     state parameters dit, sigus, and taubar.  The apparent
c     sig.n can be calculated from taubar.
c
c     Definition of energy variables:
c
c     us = position of Ei at surface relative to Ei in the bulk
c     uf = position of Ef relative to Ei in the bulk
c     ut = position of trap at energy Et relative to Ei in the bulk
c
c     uts = ut - us (position of the trap relative to midgap)
c     utf = ut - uf (position of trap relative to the Fermi level)
c     ufs = uf - us (position of the Fermi level relative to midgap
c                    at the surface)
c
c     us = the position of midgap at the surface relative to midgap
c     in the bulk
c
c     usbar = the AVERAGE value of us, as set by the gate voltage
c
c     Sign convention for energy variables:  us, uf, and ut are
c     POSITIVE DOWNWARD on a band diagram.  Hence, they are regarded
c     as potentials, not energies.  uf is negative for an n-type
c     sample.  us < 0 represents depletion (or inversion), while
c     us > 0 represents accumulation.
c
c     Note that utf = uts - ufs.
c
c        SCHEMATIC BAND DIAGRAM (Potentials positive downward)
c
c                                  |**
c ------o-------o--------------o-  M   ****
c       ^       ^     Utf < 0  |   |       ********************* Ec
c       |    Uts < 0       ----o---|-------------------o-------- Ef
c   Ut < 0      |     Ufs < 0  |   |                   ^
c       |   ----o--------------o---|**                 | Uf < 0
c       |              Us < 0  |   |   ****            |
c ------o----------------------o-- |       ************o******** Ei
c                                  |
c                                  |
c                                  |**
c                                      ****
c                                          *********************Ev
c
c
c     gpwexact5.F allows for a uniform distribution of states
c     into the oxide, yielding an additional time constant
c     dispersion (if the parameter wtunnel is greater than zero).
c     Reference: JAC PhD thesis (1973), Chapter 2 and Appendix G.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real nd,ktq,ni,ld
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      common /para2/qfbar,sigqf,qitfb,vgvfbbar,usmax,usmin
      common /dit1/dm,dv,dc,evd,ecd,ec,em,dlhalf,duhalf,da
      common /dit2/dmin,dc0,betac,dv0,betav,ucs,uvs
      common /sign/signmax,signcb0,gammacb,signvb0,gammavb
      common /usvg/vgvfb(600),usvg(600)
      dimension w(151),gpwd(11),ufsval(11)
      character outfile*20

c     define fundamental parameters

      q=1.6e-19
      e0=8.85e-14
      ktq=0.0259
      vt=1e7
      
c     define 4H-SiC and SiO2 material parameters
      
c=====================
      eg=3.23
      ni=6.66e-9
c=====================
      ug=eg/ktq
      es=10*e0
      eox=3.9*e0
      ld=sqrt(es*ktq/(2*q*ni))
      
c     define tunneling parameters

      snxfact=1
      
c     define specific device parameters
      
c=====================
      nd=2e17
      tox=17e-7
c=====================
      uf=-log(nd/ni)
      cox=eox/tox
      phims=0
      qfbar=q*1e13
      
c     Define energy-dependent dit parameters.  Dit is exponential
c     toward both bands and constant near midgap.  ditcb0 is the
c     extrapolated value at midgap of the dit equation for the
c     upper half of the bandgap, and betacb is the exponential
c     slope factor for the upper half.  ditvb0 and betavb are the
c     corresponding variables for the lower half of the bandgap.
c     ditmin is the value of dit in the middle of the bandgap.

c     (For now, set dit in the lower half to be a mirror image
c      of dit in the upper half.)

c     da = 1 means states are donor-like in the upper half and
c     acceptor-like in the lower half. da = 0 is the opposite.

      dmin=3e11
      dc0=1.36e7
      betac=0.207
      
c==================================
      dmin=5e11
      betac=0
c==================================
      
c     Setting betac = betav = 0 forces dit(E) = dmin = constant.
c     betac=0

      dv0=dc0
      betav=betac
      da=1
      
c     Define surface state density parameters in a format that is
c     compatible with subroutine qitsub. Log(Dit) decreases
c     linearly from "dc" at Ec to "dm" at energy "ecd", is constant
c     from "ecd" to "evd", and increases linearly to "dv" at Ev.
c     Ev is used as the reference energy in qitsub, so Ev = 0.

      dm=dmin
      ec=eg
      em=eg/2
      
c     Check that the exponential definition of dit in the lower
c     half of the BG actually intersects dm between ev and ei.
c     If not, set dit = dm in the lower half of the bandgap.

      dv=dv0*exp(betav*ug/2)
      if(dv.gt.dm.and.dv0.lt.dm)then
         evd=em-ktq/betav*log(dm/dv0)
      else
         dv=dm
         evd=em/2
      endif
      
c     Check that the exponential definition of dit in the upper
c     half of the BG actually intersects dm between ei and ec.
c     If not, set dit = dm in the upper half of the bandgap.

      dc=dc0*exp(betac*ug/2)
      if(dc.gt.dm.and.dc0.lt.dm)then
         ecd=em+ktq/betac*log(dm/dc0)
      else
         dc=dm
         ecd=3*em/4
      endif

c     Calculate the total number of states in the lower and upper
c     halves of the bandgap by integrating Dit(E) given above.

      if(dv.eq.dm)then
      dlhalf=dm*(em-0)
      else
      dlhalf=dm*((em-evd)+(evd-0)/log(dv/dm)*(dv/dm-1))
      endif
      if(dc.eq.dm)then
      duhalf=dm*(ec-em)
      else
      duhalf=dm*((ecd-em)+(ec-ecd)/log(dc/dm)*(dc/dm-1))
      endif
         
c     Define intersection points where the exponential segments
c     cross the constant middle segment.

      if(betac.ne.0)then
        ucs=-log(dmin/dc0)/betac
      else
        ucs=-ug/2
      endif
      if(ucs.lt.-ug/2)ucs=-ug/2
      
      if(betav.ne.0)then
        uvs=log(dmin/dv0)/betav
      else
        uvs=ug/2
      endif
      if(uvs.gt.ug/2)uvs=ug/2
      
c===============================================================

      print *,'dc, dm, dv =',dc,dm,dv
      print *,'ecd, evd =',ecd,evd
      print *,'Total Dit in upper/lower half =',duhalf,dlhalf
      
c===============================================================
      
c     Calculate the average flat-band voltage, vfbbar.

      call qitsub(0,qitfb)

      vfbbar=phims-qfbar/cox-qitfb/cox
      print *,'Qit(FB)/q =',qitfb/q,', Vfb(bar) =',vfbbar
      
c     Create a table of (vg-vfb) vs us, including spreading due
c     to interface states.  This table will be used during
c     execution of the program to look up the values for local
c     surface potential uslocal, based on the local fixed charge
c     density, qflocal.

      call usvgtable
      
c===============================================================
c     open(4,file='usvg.csv')
c     do 5 i=1,600
c   5 write(4,*)vgvfb(i),',',usvg(i)
c     stop
c===============================================================
      
c     Define energy-dependent sig.n parameters. sign is exponential
c     toward both bands and constant near midgap.  signcb0 is the
c     extrapolated value at midgap of the sign equation for the
c     upper half of the bandgap, and gammacb is the exponential
c     slope factor for the upper half.  signvb0 and gammavb are the
c     corresponding variables for the lower half of the bandgap.
c     signmax is the value of dit in the middle of the bandgap.

c     (For now, set sign = signmax in the lower half of bandgap.)

      signmax=1e12
      signcb0=1.867e4
      gammacb=0.785
 
c===================================
c     Use Ashish's parameters...
c===================================
c     signcb0=2.341e3
c     gammacb=0.85
c     signmax=1e12
      
c===============================================================
c     Set signcb0 and gammacb to match earlier run of gpwerr.F
c===============================================================
c     signcb0=1.84e3
c     gammacb=0.9
      
      signvb0=signmax
      gammavb=0
      
c     Read sigqfq.
      
      print *,'Enter Sigma(Qf/q) and tunneling width (cm)...'
      read *,sigqfq,wtunnel
      
      sigqf=q*sigqfq
      
c     Read output file name and open output file.
      
      print *,'Enter output file name...'
      read 901,outfile
  901 format(a20)
  
      open(4,file=outfile)
      
c     Set up frequencies at which to calculate gpwd.

      nfreq=151
      fstart=1e-3
      
      w(1)=1.00*fstart
      w(2)=1.26*fstart
      w(3)=1.58*fstart
      w(4)=2.00*fstart
      w(5)=2.51*fstart
      w(6)=3.16*fstart
      w(7)=3.98*fstart
      w(8)=5.01*fstart
      w(9)=6.31*fstart
      w(10)=7.94*fstart

      do 10 i=11,nfreq
   10 w(i)=10*w(i-10)
   
c     Set up ufs values at which to calculate gpwd.

      nufs=10
      
      ufsval(1)=-60
      
c=====================
      nufs=3
      ufsval(1)=-60
      delufs=10
c=====================
      
      do 15 j=2,nufs
   15 ufsval(j)=ufsval(j-1)+˙˙
   
      write(4,*)0.0,',',ufsval(1),',',ufsval(2),',',ufsval(3),',',
     #ufsval(4),',',ufsval(5),',',ufsval(6),',',ufsval(7),',',
     #ufsval(8),',',ufsval(9),',',ufsval(10)
     
      usmax=-5000
      usmin=5000

c---------------------------------------------------------------------   
c     Step frequency and calculate gpwd point for each ufs curve.
c---------------------------------------------------------------------

      do 50 i=1,nfreq

c---------------------------------------------------------------------      
c     Step through all ufs points to be calculated for this frequency.
c---------------------------------------------------------------------

      do 30 j=1,nufs
      
      ufsbar=ufsval(j)
      usbar=uf-ufsbar
      
c     Calculate the quantity vgvfbbar = (vg-vfbbar).

      us=usbar
      
c     Must write F(us,uf) this way to avoid exp(us) overflow.

      f=sqrt(exp(uf-us)+exp(uf)*us-exp(uf)+exp(us-uf)-exp(-uf)*us-
     #exp(-uf))

      sign=1
      if(us.lt.0)sign=-1

      call qitsub(us,qit)

      vgvfbbar=ktq*(us+sign*es/(cox*ld)*f)-(qit-qitfb)/cox
      
c     Calculate maximum range of uts and ufs within integrals.

c     ufsmax=ufsbar+3*sigus0
c     ufsmin=ufsbar-3*sigus0
c     utsmax=ufsmax+5
c     utsmin=ufsmin-5    
c     print *,'range of uts values within integrals'
c     print *,'utsmax =',utsmax,', band edge at uts =',ug/2
c     print *,'utsmin =',utsmin,', band edge at uts =',-ug/2
c     print *,'ufsmax =',ufsmax
c     print *,'ufsmin =',ufsmin
   
c     Calculate gpwd(i) for each w(i).

   30 call gpwdqf(w(i),gpwd(j))
      
      write(4,*)w(i),',',gpwd(1),',',gpwd(2),',',gpwd(3),',',
     #gpwd(4),',',gpwd(5),',',gpwd(6),',',gpwd(7),',',
     #gpwd(8),',',gpwd(9),',',gpwd(10)
   
   50 continue
   
      print *,'------------------------------------------'
      print *,'Maximum Us Encountered =',usmax
      print *,'Minimum Us Encountered =',usmin
   
      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine "gpwdqf" calculates a theoretical gp/w value at
c     frequency wi by integrating the gp/w expression over fixed
c     charge qf from qfbar-3*sigqf to qfbar+3*sigqf.  At each
c     value of qf in the integral, the gp/w expression is
c     integrated over energy from uf-5 to uf+5.
c     The energy dependence of dit(E) and the electron capture
c     cross section sign(E) are included explicitly.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gpwdqf(wi,gpwdi)
      
      real ktq,ni,ld
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      common /para2/qfbar,sigqf,qitfb,vgvfbbar,usmax,usmin
      
      nmax=100
      delqf=6*sigqf/nmax
      sigqf2=sigqf**2

c     Evaluate gpw integral from (qf-qfbar)=-3*sigqf to (qf-qfbar)=
c     +3*sigqf.  Variable of integration is dqf=(qf-qfbar).  For
c     each value of qf, calculate local vfb=vfb* from the equation
c     (vg0-vfb*)=(vg0-vfbbar)+(qf-qfbar)/cox.  Then look up the
c     local value of surface potential us* from a pre-computed
c     table of (vg0-vfb*) vs us*.  The quantity (vg0-vfbbar) is
c     calculated by calling function usvg(vgvfbbar,usbar), which
c     uses the exact us-vg relationship (Reference: Appendix B,
c     Semiconductor Device Fundamentals, R. F. Pierret), e.g.
c     (vg0-vfbbar)=ktq*[usbar+sign(usbar)*es/(cox*ld)*F(usbar,uf)]
c     -[qit(usbar)-qit(0)]/cox.  The second term corrects for the
c     charge in interface states when us=usbar.

      sum=0
      do 20 j=1,nmax+1
      dqf=-3.0*sigqf+(j-1)*delqf
      qflocal=dqf+qfbar
      
c     look up local value of us for this local value of qf

      vgvfb=vgvfbbar+(qflocal-qfbar)/cox
      call usvglook(uslocal,vgvfb)
      if(uslocal.gt.usmax)usmax=uslocal
      if(uslocal.lt.usmin)usmin=uslocal
      
c     call subroutine to integrate Gp/w function over energy utf
c     within this patch having local surface potential uslocal
      
      ufslocal=uf-uslocal
      
c     Call a different subroutine if tunneling into oxide traps
c     is included in the model (wtunnel > 0).
      
      if(wtunnel.gt.0) then
        call gpwx(ufslocal,wi,gpwlocal)
        else
        call gpwutf(ufslocal,wi,gpwlocal)
      endif
      
      if(gpwlocal.gt.0)sum=sum+gpwlocal*exp(-dqf**2/(2*sigqf2))
      
c     If gpwlocal is NaN, the comparison above will be FALSE,
c     (see p. 461 of Sun F77 reference), and NaN will not be
c     added to the sum.
      
c===============================================
c     print *,j,qflocal/q,vgvfb,uslocal,gpwlocal,sum
c===============================================

   20 continue

      gpwdi=sum*delqf/sqrt(2*3.14159*sigqf2)

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine "gpwx" is called if the tunneling parameter
c     wtunnel.gt.0.  "gpwx" assumes the states are uniformly
c     distributed with respect to distance into the oxide,
c     normalized so that the integral of Dit(x) w/r/t x is unity.
c     "gpwx" integrates gp/w over distance into the oxide,
c     calculatinga capture cross section reduction factor for
c     each x according to snxfact = exp(-2*k0*x).  Reference:
c     JAC PhD thesis (1973), Chapter 2 and Appendix G.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gpwx(ufslocal,wi,gpwlocal)
      
      real ns,ktq,ni,ld,k0
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      
      k0=5e7
      delx=wtunnel/10
      sum=0
      
      do 10 i=1,10
      
      x=delx*(i-1)
      snxfact=exp(-2*k0*x)
      
      call gpwutf(ufslocal,wi,gpwlocal)
 
c-------------------------------------------
c     Rectangular state distribution w/r/t x.
c-------------------------------------------
   10 sum=sum+gpwlocal

c-------------------------------------------
c     Triangular state distribution w/r/t x.
c-------------------------------------------
c  10 sum=sum+gpwlocal*2*(wtunnel-x)/wtunnel
   
      gpwlocal=sum*delx/wtunnel
   
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine "gpwutf" calculates gpw for a given angluar
c     frequency wi in a patch where the local surface potential
c     is uslocal.  This places the Fermi level at ufslocal =
c     uf-uslocal relative to midgap.
c     This subroutine integrates the gpw expression over bandgap
c     energy from Ef-5kT (utf=-5) to Ef+5kT (utf=5).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gpwutf(ufslocal,wi,gpwlocal)
      
      real ns,ktq,ni,ld
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      common /para2/qfbar,sigqf,qitfb,vgvfbbar,usmx,usmin
      
      deluts=0.1
      uts=ufslocal-5-deluts
      sum=0
      count=0
      
   10 uts=uts+deluts
      count=count+1
   
      utf=uts-ufslocal
      ft=1/(1+exp(-utf))
      ns=ni*exp(-ufslocal)
      
      cs=dit(uts)*ft*(1-ft)
      sn=sigman(uts)
      sn=sn*snxfact
      
c     Check to see if sn = 0 or sn = Inf or ns = 0 or ns = Inf.
c     If so, don't compute the integrand, to avoid result of
c     "not a number" (NaN). Reference on extreme numbers:
c     pp.460-461 of Sun F77 manual (on Mac).

      if(sn.eq.0.or.1/sn.eq.0.or.ns.eq.0.or.1/ns.eq.0)go to 20
      
      ts=ft/(sn*vt*ns)
      wt=wi*ts
      
      sum=sum+cs*wt/(1+wt**2)
      
   20 if(uts.lt.ufslocal+5.and.count.le.100)go to 10
      
      gpwlocal=q*sum*deluts
      
      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine usvgtable creates a table of us and corresponding
c     (vg-vfb) values.  The table will be used to look up the us
c     value corresponding to the local value of (vg - vfb) at
c     each (i,j) grid point.  Reference:  Appendix B, Semiconductor
c     Device Fundamentals, R. F. Pierret.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine usvgtable

      real ktq,ni,ld
      common /usvg/vgvfb(600),usvg(600)
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      common /para2/qfbar,sigqf,qitfb,vgvfbbar,usmax,usmin
      
c     Iterate us values from us=ustart to us=ustop, calculating
c     the corresponding values of vgvfb(i)=(vg-vfb).

c     ustart is negative for n-type.  ustart is the band bending
c     that places the Fermi level at Ev (strong inversion).
c
c     ustop is positive for n-type.  ustop is the band bending
c     that places the Fermi level at Ec (strong accumulation)
c
c     ustep = (total us range)/600

c     ustart and ustop are adjusted to cover the full range used
c     in Ashish's MatLab program

      ustart=-(ug/2+abs(uf))-8
      ustop=ug/2-abs(uf)+8
      ustep=abs(ustop-ustart)/600

      do 10 i=1,600

      us=ustart+(i-1)*ustep
      usvg(i)=us

c     Must write F(us,uf) this way to avoid exp(us) overflow.

      f=sqrt(exp(uf-us)+exp(uf)*us-exp(uf)+exp(us-uf)-exp(-uf)*us-
     #exp(-uf))

      sign=1
      if(us.lt.0)sign=-1

      call qitsub(us,qit)

   10 vgvfb(i)=ktq*(us+sign*es/(cox*ld)*f)-(qit-qitfb)/cox

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine usvglook interpolates in a precomputed table of
c     us and corresponding (vg-vfb) to return the local surface
c     potential coresponding to the local value of (vg-vfb), where
c     vfb is determined by the local value of fixed charge, qflocal.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine usvglook(uslocal,vgvfbx)
      
      common /usvg/vgvfb(600),usvg(600)
      
c     If vgvfbx is more negative than the most negative entry
c     in the table, set uslocal equal to the most negative
c     uslocal entry in the table and return.
      
      if(vgvfbx.gt.vgvfb(1))go to 10
      uslocal=usvg(1)
      return
      
c     Step through the table to find uslocal.
      
   10 i=1
      
   20 i=i+1
      if(i.gt.600)go to 40
      
      if(vgvfbx.gt.vgvfb(i))go to 20
      
      uslocal=usvg(i-1)+(usvg(i)-usvg(i-1))/(vgvfb(i)-vgvfb(i-1))*
     #(vgvfbx-vgvfb(i-1))
      return
      
   40 uslocal=usvg(600)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Function "dit(uts)" returns dit at the energy uts.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function dit(uts)
      
      common /dit2/dmin,dc0,betac,dv0,betav,ucs,uvs
      
      if(uts.lt.0)then
      dit=dc0*exp(-betac*uts)
      else
      dit=dv0*exp(betav*uts)
      endif
      
      if(dit.lt.dmin)dit=dmin
      
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Function "sigman(uts)" returns the majority carrier (electron)
c     capture cross section sign(uts) at the energy uts.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function sigman(uts)
      
      common /sign/signmax,signcb0,gammacb,signvb0,gammavb
      
      if(uts.lt.0)then
      sigman=signcb0*exp(gammacb*uts)
      else
      sigman=signvb0*exp(-gammavb*uts)
      endif
      
      if(sigman.gt.signmax)sigman=signmax
      
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine qitsub calculates the charge in interface states,
c     given the position of the Fermi level within the bandgap.
c     Dit(E) is modeled as a constant value "dm" in the middle
c     of the bandgap, increasing exponentially toward both band
c     edges. Log(Dit) decreases linearly from "dc" at Ec to "dm"
c     at energy "ecd", is constant at value "dm" from "ecd" to
c     "evd", and increases linearly to "dv" at Ev.  Ev is used
c     as the reference energy, so Ev = 0.  If variable "da" = 1,
c     states in the upper half of the BG are donor-like and in
c     the lower half are acceptor-like. If "da" = 0, the opposite.
c
c     duhalf is the total number of states in the upper half
c     fuhalf is the number of FILLED states in the upper half
c
c     dlhalf is the total number of states in the lower half
c     flhalf is the number of FILLED states in the lower half     
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine qitsub(us,qit)
      
      real ktq,ni,ld
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      common /para2/qfbar,sigqf,qitfb,vgvfbbar,usmax,usmin
      common /dit1/dm,dv,dc,evd,ecd,ec,em,dlhalf,duhalf,da
      
      ef=em+ktq*(us-uf)

      if(ef.gt.evd)go to 10

c     Ev < Ef < Evd

      if(dv.eq.dm)then
      flhalf=dm*(ef-0)
      else
      flhalf=dm*(evd-0)/log(dv/dm)*(dv/dm-(dv/dm)**((evd-ef)/(evd-0)))
      endif
      fuhalf=0
      go to 50

   10 if(ef.gt.em)go to 20

c     Evd < Ef < Em

      if(dv.eq.dm)then
      flhalf=dm*(ef-0)
      else
      flhalf=dm*(evd-0)/log(dv/dm)*(dv/dm-1)+dm*(ef-evd)
      endif
      fuhalf=0
      go to 50

   20 if(ef.gt.ecd)go to 30

c     Em < Ef < Ecd

      if(dv.eq.dm)then
      flhalf=dm*(em-0)
      else
      flhalf=dm*(evd-0)/log(dv/dm)*(dv/dm-1)+dm*(em-evd)
      endif
      fuhalf=dm*(ef-em)
      go to 50

c     Ecd < Ef < Ec

   30 if(dv.eq.dm)then
      flhalf=dm*(em-0)
      else
      flhalf=dm*(evd-0)/log(dv/dm)*(dv/dm-1)+dm*(em-evd)
      endif
      if(dc.eq.dm)then
      fuhalf=dm*(ef-em)
      else
      fuhalf=dm*(ecd-em)+dm*(ec-ecd)/log(dc/dm)*((dc/dm)**
     #((ef-ecd)/(ec-ecd))-1)
      endif

   50 if(da.eq.0)go to 60

c     Dit in upper half are donor-like (neutral when filled)
c     Dit in lower half are acceptor-like (neg when filled)

      quhalf=+q*(duhalf-fuhalf)
      qlhalf=-q*flhalf
      qit=quhalf+qlhalf
      return

c     Dit in upper half are acceptor-like (neg when filled)
c     Dit in lower half are donor-like (neutral when filled)

   60 quhalf=-q*fuhalf
      qlhalf=+q*(dlhalf-flhalf)
      qit=quhalf+qlhalf
      return 

      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine qitsub2 calculates the charge in interface states,
c     given the position of the Fermi level ufs within the bandgap.
c     Dit(uts) is modeled as a constant value "dm" in the middle
c     of the bandgap, and as an exponential toward both bands.
c     dc0 is the extrapolated value at midgap of the dit equation
c     for the upper half of the BG, and betac is the exponential
c     slope factor for the upper half.  dv0 and betav are the
c     corresponding variables for the lower half of the bandgap.
c     dmin is the value of dit in the middle of the bandgap.
c
c     If variable "da" = 1, states are donor-like in the upper
c     half of the BG and acceptor-like in the lower half.  If
c     "da" = 0, the opposite.
c
c     duhalf = total number of states in the upper half
c     fuhalf = number of FILLED states in the upper half
c     euhalf = number of EMPTY states in the upper half
c
c     dlhalf = total number of states in the lower half
c     flhalf = number of FILLED states in the lower half
c     elhalf = number of EMPTY states in the lower half
c 
c             ^  Dit(E-Ev)                            ^ Dit(uts)
c             |                   .                   |
c             |                   .                   |
c          dv |*       slope =    .  slope =         *| dc
c             |  * .....betav     .   betac ...... *  |
c             |    *              .              *    |
c             |------*-------------------------*------| dm=dmin
c             |        *          .          *        |
c             |          *        .        *          |
c             |            *      .      *            |
c             |              *    .    *              |
c             |                *  .  *                |
c             |           dv0 –>  *  <– dc0           |
c             |–––––––––––––––––––|–––––––––––––––––––|––> (E-Ev)
c           ev=0     ^            em           ^    ec=eg
c                    |                         |
c                   evd                       ecd
c       uts   <––––––|––––––––––––|––––––––––––|–––––––
c                  uvs>0          0          ucs<0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine qitsub2(us,qit)
      
      real ktq,ni,ld
      common /para1/q,ktq,uf,ug,ni,vt,es,cox,ld,wtunnel,snxfact
      common /para2/qfbar,sigqf,qitfb,vgvfbbar,usmax,usmin
      common /dit1/dm,dv,dc,evd,ecd,eg,em,dlhalf,duhalf,da
      common /dit2/dmin,dc0,betac,dv0,betav,ucs,uvs
      
      ufs=uf-us
      ug=eg/ktq
      
c     Calculate density of filled states in the lower half of BG.

      if(betav.gt.0)go to 10
      flhalf=dmin*ktq*(ug/2-ufs)
      if(ufs.lt.0)flhalf=dmin*ktq*ug/2
      go to 20
      
   10 if(ufs.gt.uvs)then
        flhalf=dv0*kTq/betav*(exp(betav*ug/2)-exp(betav*ufs))
      else if(ufs.le.0)then
        flhalf=dv0*ktq/betav*(exp(betav*ug/2)-exp(betav*uvs))+
     #  dmin*ktq*uvs
      else
        flhalf=dv0*ktq/betav*(exp(betav*ug/2)-exp(betav*uvs))-
     #  dmin*ktq*(ufs-uvs)
      endif
      
c     Calculate density of empty states in the upper half of BG.

   20 if(betac.gt.0)go to 30
      euhalf=-dmin*ktq*(-ug/2-ufs)
      if(ufs.gt.0)euhalf=dmin*ktq*ug/2
      go to 40

   30 if(ufs.lt.ucs)then
        euhalf=dc0*ktq/betac*(exp(betac*ug/2)-exp(-betac*ufs))
      else if(ufs.ge.0)then
        euhalf=dc0*ktq/betac*(exp(betac*ug/2)-exp(-betac*ucs))-
     #  dmin*ktq*ucs
      else
        euhalf=dc0*ktq/betac*(exp(betac*ug/2)-exp(-betac*ucs))-
     #  dm*ktq*(ucs-ufs)
      endif
      
   40 elhalf=dlhalf-flhalf
      fuhalf=duhalf-euhalf
      
      if(da.eq.1)then
        qit=q*(euhalf-flhalf)
      else
        qit=q*(elhalf-fuhalf)
      endif
      
      return 
      end
      