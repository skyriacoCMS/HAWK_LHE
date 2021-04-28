#define PRIVATE

***********************************************************************
*     file public.F                                                   *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     30.01.06  Ansgar Denner     last changed  15.02.12              *
***********************************************************************

* read input from file
#define READINPUT

#define RANLUX


#ifdef RECOLA
#ifdef RECOLA1
#include "recola_drivers/init_recola1.f90"
#endif
#ifdef RECOLA2
#include "recola_drivers/init_recola2.f90"
#endif
#endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     read input file                                             c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readinput
#ifdef RECOLA2
        use recola, only: get_modelname_rcl
#endif
      implicit none

#include "mcinput.h"
#include "mcparams.h"
#include "mccommon.h"

#ifdef RECOLA2
        character*100 :: modelname
#endif

c for initcarry
      integer x(0:24),carry,initsteps
      common/randomseed/x,carry,initsteps

c local variables
      integer i1,isep1,isep2,linputline,hemicut,higgsbjets
      integer printweightmax
      integer sonsborn,sonsqcd
      integer lqedfacscale,lqcdrenscale
      character*100 identifier,orgidentifier
      character*101 directive,orgdirective,storeinput
      character*100 inputline
      integer      autofreq

      logical has_photon

      data linputline /100/, autofreq /1/, printweightmax /1/
     &    lprintweightmax /.true./


c anomalous HVV couplings
      integer qhvv
      real*8 rsm,d,db,dt,dtb,lambdahvv,
     &             a1hww,a2hww,a3hww,a1haa,a2haa,a3haa,
     &             a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz
      common/hvv/rsm,d,db,dt,dtb,lambdahvv,
     &             a1hww,a2hww,a3hww,a1haa,a2haa,a3haa,
     &             a1haz,a2haz,a3haz,a1hzz,a2hzz,a3hzz,qhvv

#ifdef RANLUX
      integer ranluxseed
      common/ranluxseed/ranluxseed
#endif

      integer qalps
      common/qalps/qalps

      integer selprocess
      common/selprocess/selprocess

      integer num


#ifdef multpdf
#include "multpdf.h"
      integer nset
#endif

#ifdef RECOLA2
        call get_modelname_rcl(modelname)
#endif

c initialization
      lcmsgen=.true.
      lcmsgen=.false.
      do gen=0,maxg
        sgen(gen)=0
        sborn(gen)=0
        swborn(gen)=0
        svirt(gen)=0
        spbrems(gen)=0
        sgbrems(gen)=0
        spinc(gen)=0
        sginc(gen)=0
        ssub(gen)=0
        ssli(gen)=0
        seff(gen)=0
        sendp(gen)=0
        scuts(gen)=0
        srecomb(gen)=0
        sconv(gen)=0
        soaisr(gen)=0
        sincludecuts(gen)=1
        do icon=1,maxcon
          freqgencon(gen,icon)=1
          sgencon(gen,icon)=0
        enddo
        nextgen(gen)=5
      enddo
      nextmax=nextgen(0)

c test input

      outputfile=' '
      outputfile2='sigma.dat'
      outputvar=' '
      energy=13000d0

      nevents=100 000 000
       nevents=20



      outputfile=' '
      outputfile2='sigma2.dat'

c initsteps: initialization steps for initcarry
      initsteps = 100000

c lnoutmc=5 full output of montecarlo
c lnoutmc=3 medium output of montecarlo
c lnoutmc=1 standard output of montecarlo
      lnoutmc=1

c lnoutgen=1 minimal output of generator
c lnoutgen=0 no output of generator
      lnoutgen=0

c maximal number of error messages
      maxout=100
      maxout=0

c selprocess=0  jet jet H   final state
c selprocess=1  l+  nu  H   final state
c selprocess=2  nu~ l-  H   final state
c selprocess=3  l+  l-  H   final state
c selprocess=4  nu~ nu  H   final state
      selprocess=0

c sppbar=0  p p initial state
c sppbar=1  p pbar initial state
c sppbar=2  pbar pbar initial state
c sppbar=3  pbar p initial state
      sppbar=1
      sppbar=0

c salp=0   alpha(0) scheme
c salp=1   alpha(mz) scheme
c salp=2   alpha(gf) scheme
      salp =2

c sborn=0   Born not included
c sborn=1   Born included
      sborn(0)=1

c sconv=0   no convolution with pdf
c sconv=1   convolution with pdf
      sconv(0)=1

c swidth=0  zero width
c swidth=1  complex mass scheme
c swidth=2  naive fixed width scheme
      swidth=1

c swdiag=0  W-exchange diagrams excluded
c swdiag=1  W-exchange diagrams included
      swdiag=1

c szdiag=0  Z-exchange diagrams excluded
c szdiag=1  Z-exchange diagrams included
      szdiag=1

c scp=0  CP symmetry not used
c scp=1  CP symmetry used in lowest order matrix elements
      scp=0

c svirt=0   virtual corrections not included
c svirt=1   virtual corrections included
      svirt(0)=1

c sew=0  electroweak corrections not included
c sew=1  electroweak corrections included
      sew=1

c sweak=0  weak virtual corrections not included
c sweak=1  weak virtual (fermionic and bosonic) corrections included
c sweak=1  weak virtual fermionic corrections included
c sweak=1  weak virtual bosonic corrections included
      sweak=1

c sqcd=0  qcd corrections not included
c sqcd=1  qcd corrections included
      sqcd=1

c sqcddiag=0  diagonal qcd corrections not included
c sqcddiag=1  diagonal qcd corrections included
      sqcddiag=1

c the switches sqcdnondiag, sqcdgsplit, sqcdggfus are not relevant for WHZH
c sqcdnondiag=0  nondiagonal qcd corrections not included
c sqcdnondiag=1  nondiagonal qcd corrections included
      sqcdnondiag=1

c sqcdgsplit=0  qqg qcd corrections not included
c sqcdgsplit=1  qqg qcd corrections included
      sqcdgsplit=1

c sqcdggfus=0  gg fusion qcd corrections not included
c sqcdggfus=1  gg fusion qcd corrections included
      sqcdggfus=1

c The switches sscha, stcha, sch2 and schint must be set to 1 for WH/ZH
c sscha= 0:  switches off contributions connected to s-channel W/Z's in LO
c sscha= 1:  switches on contributions connected to s-channel W/Z's in LO
      sscha= 1

c stcha= 0:  switches off contributions connected to t-channel W/Z's in LO
c stcha= 1:  switches on contributions connected to t-channel W/Z's in LO
      stcha= 1
c     Note: qschan & qtchan =/= 0 do not work for photon-induced processes,
c           because this splitting is not gauge invariant there.

c sch2  = 0:  switches off contributions connected to squared LO diagrams
c sch2  = 1:  switches on contributions connected to squared LO diagrams
      sch2  = 1

c schint= 0:  switches off contributions connected to LO interferences
c schint= 1:  switches on contributions connected to LO interferences
      schint= 1

c sbini = 0:  switches off b-quarks in the initial state
c sbini = 1:  switches on b-quarks in the initial state
      sbini=1

c sbfin = 0:  switches off b-quarks in the final state
c sbfin = 1:  switches on b-quarks in the final state
      sbfin=1
c                Note: channels with external b-quarks included only in LO,
c                      and qbini & qbfin are only active in subroutines
c                      Mat2born, Mat2born1, Mat2born2.

c shh2  = 0:  switches off 2-loop heavy-Higgs effects of O(GF^2*MH^4)
c shh2  = 1:  switches on 2-loop heavy-Higgs effects of O(GF^2*MH^4)
c                        included for squared t-channel diagrams
c                        following Frink et al., Phys.Rev.D54 (1996) 4548
c                                                [hep-ph/9606310]
c does not apply to WH/ZH branch  (contributions not implemented)
#ifdef RECOLA2
      shh2=0
#else
      shh2=1
#endif

c sfactqed=0   no factorization of QED
c sfactqed=1   MSbar factorization of QED
c sfactqed=2   DIS factorization of QED
      sfactqed=2

c sfactqcd=0   no factorization for QCD
c sfactqcd=1   MSbar factorization of QCD
c sfactqcd=2   DIS factorization of QCD
      sfactqcd=1

c spbrems(1)=0 photon bremsstrahlung not included
c spbrems(1)=1 photon bremsstrahlung included
      spbrems(1)=1
      spbrems(0)=spbrems(1)

c sgbrems(2)=0 gluon bremsstrahlung not included
c sgbrems(2)=1 gluon bremsstrahlung included
      sgbrems(2)=1
      sgbrems(0)=sgbrems(2)

c spinc(3)=0 incoming photon contribution 1 not included
c spinc(3)=1 incoming photon contribution 1 included
      spinc(3)=1
      spinc(0)=spinc(3)

c spinc(4)=0 incoming photon contribution 2 not included
c spinc(4)=1 incoming photon contribution 2 included
      spinc(4)=1

c sginc(5)=0 incoming gluon contribution 1 not included
c sginc(5)=1 incoming gluon contribution 1 included
      sginc(5)=1
      sginc(0)=sginc(5)

c sginc(6)=0 incoming gluon contribution 2 not included
c sginc(6)=1 incoming gluon contribution 2 included
      sginc(6)=1

c ssub=1    massless subtraction terms of STD
c ssub=2    massless subtraction terms of MAR
      ssub(0)=1
      do gen=1,maxg
        ssub(gen)=ssub(0)
      enddo

c ssli=1    massless slicing terms of CHM
      ssli(0)=0
      do gen=1,maxg
        ssli(gen)=ssli(0)
      enddo

c seff=1    effective collinear factors
      seff(0)=0
      do gen=1,maxg
        seff(gen)=seff(0)
      enddo

c sendp(0)=0  no virtual endpoint contribution included
c sendp(0)=-1  virtual endpoint contribution included in M_virt
c sendp(0)=1  virtual endpoint contribution included in MC
      sendp(0)=1

c sincludecuts=1  cuts included in generator
c sincludecuts=0  cuts not included in generator
      sincludecuts(0)=1

c srecomb=0  no recombination
      srecomb(0)=0
c srecomb=1  kt-clustering   ycut=5
      srecomb(0)=1
c srecomb=2  kt-clustering   ycut=1d10
c     srecomb(0)=2

c sbarelep=1  no photon-lepton recombination
c sbarelep=0  photon-lepton recombination
      sbarelep=1

c colllep=0  no initial-photon to collinear lepton splitting
c colllep=1  with  initial-photon to collinear lepton splitting
c set internally if ptlcut=0
c      colllep=0

c scuts=0   no cuts applied
      scuts(0)=0
c scuts=0   no predefined cuts
      scuts(0)=-1
c scuts=1   standard ppqqh cuts with pt ordering
c      scuts(0)=1
c scuts=2   standard ppqqh cuts with E ordering
c      scuts(0)=2


c spdf=1   MRST2004QED [with O(alpha) effects]
c      spdf=1
c spdf=2   MRST2004QED [without O(alpha) effects]
c       spdf=2
c spdf=11  cteq6m pdfs [Standard MSbar scheme (alphas(mz)=0.118)]
c       spdf=11
c spdf=14  cteq6l1 pdfs [Leading Order (alphas(mz)=0.130)]
c       spdf=14
c spdf=99  LHApdf     (needs pdfname or pdfpath)
       spdf=99
       pdfpath=''
       pdfname='NNPDF31_nlo_as_0118_luxqed'

c pdfmap=1 flat mapping of pdfs
c       pdfmap=1
c pdfmap=4 1/x and 1/tau mapping of pdfs
       pdfmap=4

c treatment of Higgs boson
c shtr=0   on-shell Higgs boson
       shtr=0
c shtr=1   off-shell Higgs boson
c      shtr=1
c shtr=2   off-shell Higgs boson decaying into a pair of massless singlets
c      shtr=2

c treatment of Higgs-boson width
c sgh=0   Higgs-boson width as specified in the input
c      sgh=0
c sgh=1   Higgs-boson width from an interpolation of arXiv:1101.0593
c      sgh=1
c sgh=2   Higgs-boson width according to complex-mass scheme
       sgh=2

c Higgs boson Breit-Wigner
c shbw=0   standard Breit-Wigner with fixed width
c      shbw=0
c shbw=1   off-shell propagator according to arXiV:1107.0683
       shbw=1
c shbw=2   Breit-Wigner with s-dependent width
c      shbw=2
c shbw=8   Breit-Wigner with s-dependent width in numerator
c      shbw=8
c shbw=9   off-shell propagator with fixed width
c      shbw=9

c on-shell projection of off-shell Higgs boson for LO
c lonsborn=F  no on-shell projection
       lonsborn=.false.
c sonsborn=T  on-shell projection
c      lonsborn=.true.

c on-shell projection of off-shell Higgs boson for QCD corrections
c lonsqcd=F  no on-shell projection
       lonsqcd=.false.
c lonsqcd=T  on-shell projection
c      lonsqcd=.true.

c on-shell projection of off-shell Higgs boson for EW corrections required
c lonsew=T  on-shell projection
       lonsew=.true.

c on-shell projection of off-shell Higgs boson
c sonsproj=0  no on-shell projection
c      sonsproj=0
c sonsproj=1  on-shell projection via rescaling of outgoing momenta
       sonsproj=1
c sonsproj=2   on-shell projection via rescaling of all momenta
c      sonsproj=2

c Higgs branching ratio
       hbr=1d0

c scalings for factorization and renormalization scales
       qcdfacscalefac=1d0
       qedfacscalefac=1d0
       qcdrenscalefac=1d0

c shvv=0  anomalous HVV couplings not included
c shvv=1  anomalous HVV couplings included
      shvv  = 0
      qhvv  = shvv
      rsm   = 1d0
      d     = 0d0
      db    = 0d0
      dt    = 0d0
      dtb   = 0d0
      lambdahvv = -1d0
      a1hww = 0d0
      a2hww = 0d0
      a3hww = 0d0
      a1haa = 0d0
      a2haa = 0d0
      a3haa = 0d0
      a1haz = 0d0
      a2haz = 0d0
      a3haz = 0d0
      a1hzz = 0d0
      a2hzz = 0d0
      a3hzz = 0d0


c produce histograms
       shisto=1
       if(selprocess.eq.0) then
         nhisto=12
       elseif(selprocess.eq.1)then
         nhisto=11
       elseif(selprocess.eq.2)then
         nhisto=11
       elseif(selprocess.eq.3)then
         nhisto=18
       elseif(selprocess.eq.4)then
         nhisto=8
       endif
       nhisto=maxh
       do i1=0,maxh
         lhisto(i1)=1
       enddo

 201  continue

c initialize parameters
      ssoft=1            !  choice of subtraction terms
      sexcl=0

      slightfermions=1


      initrecola = 0
#ifdef RECOLA2
      if (modelname .eq. "THDM" .or. modelname .eq. "HS") then
        higgsflavor = "Hl"
      else
        higgsflavor = "H"
      end if
#else
      higgsflavor = "H"
#endif

      if(selprocess.eq.0)then
        if(srecomb(0).eq.1)then
          ktd=0.4d0
          partetacut=5d0
          ktpower=1
          sbarelep=0
        elseif(srecomb(0).eq.2)then
          ktd=0.4d0
          ktpower=1
          sbarelep=0
          partetacut=1d10
        endif
        if(scuts(0).eq.1.or.scuts(0).eq.2)then
          ptcutj=20d0
          ycutj=4.5d0
          dycutjj=4d0
          lhemicut=.true.
          ycuth=1d10
          ecuth=0d0
          ptcuth=0d0
          mlcuth=1d-2             ! lower cut on Higgs invariant mass
          mucuth=1d30             ! upper cut on Higgs invariant mass
          lhiggsbjets=.false.
          ptcutj1=20d0
          ptcutj2=20d0
          ycutj1=4.5d0
          ycutj2=4.5d0
          ptmaxj1=1d30
          ptmaxj2=1d30
          ptmaxj3=1d30
          yminj1=0d0
          yminj2=0d0
          yminj3=0d0
          ptcutvis=0d0
          ptmaxvis=1d30
          dphicutj1h=0d0
          dphicutj2h=0d0
          mlcutjj=0d0
          mucutjj=1d30
          mlcutjjh=0d0
          mucutjjh=1d30
          ptcuts=0d0
          ycuts=1d10
          drcutjs=0d0
        elseif(scuts(0).eq.-1)then
          mlcuth=1d-2             ! lower cut on Higgs invariant mass
          mucuth=1d30             ! upper cut on Higgs invariant mass
          mucutjj=1d30            ! upper cut on jet-jet invariant mass
          mucutjjh=1d30           ! upper cut on jet-jet-Higgs invariant mass
          ptmaxvis=1d30           ! upper cut on pT of jet-jet-Higgs system
        elseif(scuts(0).eq.0)then
          mlcuth=1d-2             ! lower cut on Higgs invariant mass
          mucuth=1d30             ! upper cut on Higgs invariant mass
        endif
      elseif(selprocess.ge.1.and.selprocess.le.4)then
         dgammaparameter=0.1d0

         ycutl=1d10
         ptcutl=0d0
         ptmaxl=1d10
         yminl=0d0
         if(scuts(0).eq.1)then
            ycutl=2.5d0
            ptcutl=20d0
            ptmaxl=1d10
            yminl=0d0
         elseif(scuts(0).eq.2)then
            ptmaxl=20d0
            yminl=2.5d0
         endif
         ecutl=0d0

         if(scuts(0).eq.0)then
            ptcutv=0d0
            ptcutm=0d0
            ptcuth=0d0
         else
            ptcutv=190d0
            ptcutm=25d0
            ptcuth=200d0
         endif
         mlcutll=0d0
         mucutll=1d10
         mlcutllh=0d0
         mucutllh=1d10
         mcutjl=0d0
         dRcutll=0d10
         dRcutjl=0d10
         dycutll=0d0
         dycutjl=0d0
      endif


c default slicing cut parameters
      mcslideltac=1d-6
      mcslideltas=1d-4

#ifdef RANLUX
c default ranlux initialization
      ranluxseed=314159265
      ranluxseed=0
#endif

#ifdef READINPUT
       nevents=100
       pdfmember=0
#ifdef multpdf
       pdfmemberfrom=1
       pdfmemberto=0
       pdfmembercentral=pdfmember
#endif

#ifdef mpiuse
       
       open(unit=5,file='inputfile',status='old')
 50    read(5,'((a))',end=51,err=52) inputline
#else
 50    read(*,'((a))',end=51,err=52) inputline
#endif

       if(inputline(1:1).eq.'!'.or.inputline(1:1).eq.'#') goto 50
       isep1=0
       isep2=linputline+1


       do i1=1,linputline
         if(inputline(i1:i1).eq.'=') then
           isep1=i1
           goto53
         elseif(inputline(i1:i1).eq.'!') then
           isep2=i1
           goto54
         endif
       enddo
 53    continue


       do i1=isep1+1,linputline
         if(inputline(i1:i1).eq.'!') then
           isep2=i1
           goto54
         endif
       enddo
 54    continue


       if(isep1.gt.0)then
         identifier= inputline(1:isep1-1)
         orgidentifier=identifier
         do i1=isep1-1,1,-1
           if (identifier(i1:i1).eq.' ') then
             identifier=identifier(1:i1-1)//identifier(i1+1:linputline)
     &           //' '
           elseif(identifier(i1:i1).ge.'A'.and.identifier(i1:i1).le.'Z')
     &         then
             if(identifier(i1:i1).eq.'A') then
               identifier(i1:i1)='a'
             elseif(identifier(i1:i1).eq.'B') then
               identifier(i1:i1)='b'
             elseif(identifier(i1:i1).eq.'C') then
               identifier(i1:i1)='c'
             elseif(identifier(i1:i1).eq.'D') then
               identifier(i1:i1)='d'
             elseif(identifier(i1:i1).eq.'E') then
               identifier(i1:i1)='e'
             elseif(identifier(i1:i1).eq.'F') then
               identifier(i1:i1)='f'
             elseif(identifier(i1:i1).eq.'G') then
               identifier(i1:i1)='g'
             elseif(identifier(i1:i1).eq.'H') then
               identifier(i1:i1)='h'
             elseif(identifier(i1:i1).eq.'I') then
               identifier(i1:i1)='i'
             elseif(identifier(i1:i1).eq.'J') then
               identifier(i1:i1)='j'
             elseif(identifier(i1:i1).eq.'K') then
               identifier(i1:i1)='k'
             elseif(identifier(i1:i1).eq.'L') then
               identifier(i1:i1)='l'
             elseif(identifier(i1:i1).eq.'M') then
               identifier(i1:i1)='m'
             elseif(identifier(i1:i1).eq.'N') then
               identifier(i1:i1)='n'
             elseif(identifier(i1:i1).eq.'O') then
               identifier(i1:i1)='o'
             elseif(identifier(i1:i1).eq.'P') then
               identifier(i1:i1)='p'
             elseif(identifier(i1:i1).eq.'Q') then
               identifier(i1:i1)='q'
             elseif(identifier(i1:i1).eq.'R') then
               identifier(i1:i1)='r'
             elseif(identifier(i1:i1).eq.'S') then
               identifier(i1:i1)='s'
             elseif(identifier(i1:i1).eq.'T') then
               identifier(i1:i1)='t'
             elseif(identifier(i1:i1).eq.'U') then
               identifier(i1:i1)='u'
             elseif(identifier(i1:i1).eq.'V') then
               identifier(i1:i1)='v'
             elseif(identifier(i1:i1).eq.'W') then
               identifier(i1:i1)='w'
             elseif(identifier(i1:i1).eq.'X') then
               identifier(i1:i1)='x'
             elseif(identifier(i1:i1).eq.'Y') then
               identifier(i1:i1)='y'
             elseif(identifier(i1:i1).eq.'Z') then
               identifier(i1:i1)='z'
             endif
           endif
         enddo
         write(storeinput,*) inputline(isep1+1:isep2-1)


         if(identifier.eq.'nevents') then
            read(storeinput,'(i100)',err=61) nevents
            read(storeinput,*) nevents


         elseif(identifier.eq.'selprocess') then
            read(storeinput,'(i100)',err=61) selprocess
            read(storeinput,*) selprocess

            if(selprocess.ge.1.and.selprocess.le.4)then
              dgammaparameter=0.1d0
              sbarelep=1

              ycutl=1d10
              ptcutl=0d0
              ptmaxl=1d10
              yminl=0d0
              if(scuts(0).eq.1)then
                ycutl=2.5d0
                ptcutl=20d0
              elseif(scuts(0).eq.2)then
                ptmaxl=20d0
                yminl=2.5d0
              endif
              ecutl=0d0
              if(scuts(0).eq.0.or.scuts(0).eq.-1)then
                 ptcutv=0d0
                 ptcutm=0d0
                 ptcuth=0d0
              else
                 ptcutv=190d0
                 ptcutm=25d0
                 ptcuth=200d0
              endif
              mlcutll=0d0
              mucutll=1d10
              mlcutllh=0d0
              mucutllh=1d10
              mcutjl=0d0
              dRcutll=0d0
              dRcutjl=0d0
              dycutll=0d0
              dycutjl=0d0
            elseif(selprocess.eq.0)then
              if(scuts(0).eq.0.or.scuts(0).eq.-1)then
                 ptcutj=0d0
                 ycutj=1d10
                 dycutjj=0d0
                 lhemicut=.false.
                 ptcutj1=0d0
                 ptcutj2=0d0
                 ycutj1=1d10
                 ycutj2=1d10
              else
                 ptcutj=20d0
                 ycutj=4.5d0
                 dycutjj=4d0
                 lhemicut=.true.
                 ptcutj1=20d0
                 ptcutj2=20d0
                 ycutj1=4.5d0
                 ycutj2=4.5d0
              endif
              ktd=0.4d0
              partetacut=5d0
              ktpower=1
              sbarelep=0
              ycuth=1d10
              ecuth=0d0
              ptcuth=0d0
              mlcuth=1d-2       ! lower cut on Higgs invariant mass
              mucuth=1d30       ! upper cut on Higgs invariant mass
              lhiggsbjets=.false.
              yminj1=0d0
              yminj2=0d0
              yminj3=0d0
              ptmaxj1=1d30
              ptmaxj2=1d30
              ptmaxj3=1d30
              ptcutvis=0d0
              ptmaxvis=1d30
              dphicutj1h=0d0
              dphicutj2h=0d0
              mlcutjj=0d0
              mucutjj=1d30
              mlcutjjh=0d0
              mucutjjh=1d30
              ptcuts=0d0
              ycuts=1d10
              drcutjs=0d0
           else
              write(*,'(/(a))') ' Error message: HAWK Stops here!'
              write(*,*) 'The input selprocess can only be set to' //
     &                   ' 0,1,2,3,4'
              write(*,*) 'You try to use selprocess=',selprocess
              stop
            endif


         elseif(identifier.eq.'energy') then
            read(storeinput,'(f50.50)',err=61) energy
            read(storeinput,*) energy


         elseif(identifier.eq.'outputfile') then
           read(storeinput,'(a100)',err=61) outputfile
           read(storeinput,*) outputfile


         elseif(identifier.eq.'outputfile2') then
           read(storeinput,'(a100)',err=61) outputfile2
           read(storeinput,*) outputfile2


        elseif(identifier.eq.'higgsflavor' .or.
     &   identifier.eq.'higgsflavour') then
           read(storeinput,'(a2)',err=61) higgsflavor
           read(storeinput,*) higgsflavor

           if (trim(adjustl(higgsflavor)) .ne. 'H'  .and.
     &         trim(adjustl(higgsflavor)) .ne. 'Hl' .and.
     &         trim(adjustl(higgsflavor)) .ne. 'Hh' .and.
     &         trim(adjust l(higgsflavor)) .ne. 'Ha') then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) trim(adjustl(higgsflavor))
             write(*,*) 'The input higgsflavor can only be set to '//
     &                   'H,Hl,Hh,Ha'
             stop
           end if

         elseif(identifier.eq.'outputvar') then
           read(storeinput,*) outputvar


         elseif(identifier.eq.'sppbar') then
           read(storeinput,'(i100)',err=61) sppbar
           read(storeinput,*) sppbar
           if(sppbar.lt.0.or.sppbar.ge.4) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sppbar can only be set to'//
     &                   ' 0,1,2,3'
             write(*,*) 'You try to use sppbar=',sppbar
             stop
           endif
         elseif(identifier.eq.'1/alpha0') then
           read(storeinput,'(f50.50)',err=61) alpha0in
           read(storeinput,*) alpha0in
           alpha0in=1d0/alpha0in
         elseif(identifier.eq.'1/alphaz') then
           read(storeinput,'(f50.50)',err=61) alphazin
           read(storeinput,*) alphazin
           alphazin=1d0/alphazin
         elseif(identifier.eq.'gf') then
           read(storeinput,'(f50.50)',err=61) gfin
           read(storeinput,*) gfin
         elseif(identifier.eq.'mz') then
           read(storeinput,'(f50.50)',err=61) mzin
           read(storeinput,*) mzin
         elseif(identifier.eq.'gz') then
           read(storeinput,'(f50.50)',err=61) gzin
           read(storeinput,*) gzin
         elseif(identifier.eq.'mw') then
           read(storeinput,'(f50.50)',err=61) mwin
           read(storeinput,*) mwin
         elseif(identifier.eq.'gw') then
           read(storeinput,'(f50.50)',err=61) gwin
           read(storeinput,*) gwin
         elseif(identifier.eq.'me') then
           read(storeinput,'(f50.50)',err=61) mein
           read(storeinput,*) mein
         elseif(identifier.eq.'mmu') then
           read(storeinput,'(f50.50)',err=61) mmuin
           read(storeinput,*) mmuin
         elseif(identifier.eq.'mtau') then
           read(storeinput,'(f50.50)',err=61) mtauin
           read(storeinput,*) mtauin
         elseif(identifier.eq.'mu') then
           read(storeinput,'(f50.50)',err=61) muin
           read(storeinput,*) muin
         elseif(identifier.eq.'md') then
           read(storeinput,'(f50.50)',err=61) mdin
           read(storeinput,*) mdin
         elseif(identifier.eq.'ms') then
           read(storeinput,'(f50.50)',err=61) msin
           read(storeinput,*) msin
         elseif(identifier.eq.'mc') then
           read(storeinput,'(f50.50)',err=61) mcin
           read(storeinput,*) mcin
         elseif(identifier.eq.'mb') then
           read(storeinput,'(f50.50)',err=61) mbin
           read(storeinput,*) mbin
         elseif(identifier.eq.'mt') then
           read(storeinput,'(f50.50)',err=61) mtin
           read(storeinput,*) mtin
         elseif(identifier.eq.'mh') then
           read(storeinput,'(f50.50)',err=61) mhin
           read(storeinput,*) mhin
         elseif(identifier.eq.'sinthetac') then
           read(storeinput,'(f50.50)',err=61) sinthetacin
           read(storeinput,*) sinthetacin
         elseif(identifier.eq.'gh') then
           read(storeinput,'(f50.50)',err=61) ghin
           read(storeinput,*) ghin
         elseif(identifier.eq.'mumsbsm') then
           read(storeinput,'(f50.50)',err=61) mumsbsmin
           read(storeinput,*) mumsbsmin
         elseif(identifier.eq.'tb') then
           read(storeinput,'(f50.50)',err=61) tbin
           read(storeinput,*) tbin
         elseif(identifier.eq.'cab') then
           read(storeinput,'(f50.50)',err=61) cabin
           read(storeinput,*) cabin
         elseif(identifier.eq.'l5') then
           read(storeinput,'(2f50.50)',err=61) l5in
           read(storeinput,*) l5in
         elseif(identifier.eq.'ytypein') then
           read(storeinput,'(i100)',err=61) ytypein
           read(storeinput,*) ytypein
         elseif(identifier.eq.'mhh') then
           read(storeinput,'(f50.50)',err=61) mhhin
           read(storeinput,*) mhhin
         elseif(identifier.eq.'mhc') then
           read(storeinput,'(f50.50)',err=61) mhcin
           read(storeinput,*) mhcin
         elseif(identifier.eq.'mha') then
           read(storeinput,'(f50.50)',err=61) mhain
           read(storeinput,*) mhain
         elseif(identifier.eq.'ghh') then
           read(storeinput,'(f50.50)',err=61) ghhin
           read(storeinput,*) ghhin
         elseif(identifier.eq.'ghc') then
           read(storeinput,'(f50.50)',err=61) ghcin
           read(storeinput,*) ghcin
         elseif(identifier.eq.'gha') then
           read(storeinput,'(f50.50)',err=61) ghain
           read(storeinput,*) ghain
         elseif(identifier.eq.'mixrs') then
           read(storeinput,'(A16)',err=61) mixrsin
           read(storeinput,*) mixrsin
         elseif(identifier.eq.'sa') then
           read(storeinput,'(f50.50)',err=61) sain
           read(storeinput,*) sain
         elseif(identifier.eq.'l3') then
           read(storeinput,'(f50.50)',err=61) l3in
           read(storeinput,*) l3in
         elseif(identifier.eq.'hbr') then
           read(storeinput,'(f50.50)',err=61) hbr
           read(storeinput,*) hbr
         elseif(identifier.eq.'lambda') then
           read(storeinput,'(f50.50)',err=61) lambdain
           read(storeinput,*) lambdain
         elseif(identifier.eq.'qcdfacscalefac') then
           read(storeinput,'(f50.50)',err=61) qcdfacscalefac
           read(storeinput,*) qcdfacscalefac
         elseif(identifier.eq.'qcdrenscalefac') then
           read(storeinput,'(f50.50)',err=61) qcdrenscalefac
           read(storeinput,*) qcdrenscalefac
         elseif(identifier.eq.'qedfacscalefac') then
           read(storeinput,'(f50.50)',err=61) qedfacscalefac
           read(storeinput,*) qedfacscalefac
         elseif(identifier.eq.'qcdfacscalefac2') then
           read(storeinput,'(f50.50)',err=61) qcdfacscalefac
           read(storeinput,*) qcdfacscalefac
           qcdfacscalefac=sqrt(qcdfacscalefac)
         elseif(identifier.eq.'qcdrenscalefac2') then
           read(storeinput,'(f50.50)',err=61) qcdrenscalefac
           read(storeinput,*) qcdrenscalefac
           qcdrenscalefac=sqrt(qcdrenscalefac)
         elseif(identifier.eq.'qedfacscalefac2') then
           read(storeinput,'(f50.50)',err=61) qedfacscalefac
           read(storeinput,*) qedfacscalefac
           qedfacscalefac=sqrt(qedfacscalefac)
         elseif(identifier.eq.'lqcdrenscale') then
           read(storeinput,'(i100)',err=61) lqcdrenscale
           read(storeinput,*) lqcdrenscale
           if (lqcdrenscale.eq.1) then
             qcdrenscalefac=qcdfacscalefac
           elseif (lqcdrenscale.eq.-1) then
             qcdrenscalefac=1d0/qcdfacscalefac
           else
             goto 61
           endif
         elseif(identifier.eq.'lqedfacscale') then
           read(storeinput,'(i100)',err=61) lqedfacscale
           read(storeinput,*) lqedfacscale
           if (lqedfacscale.eq.1) then
             qedfacscalefac=qcdfacscalefac
           else
             goto 61
           endif
         elseif(identifier.eq.'srecomb') then
           read(storeinput,'(i100)',err=61) srecomb(0)
           read(storeinput,*) srecomb(0)
           if(selprocess.eq.0)then
             if(srecomb(0).eq.1)then
               ktd=0.4d0
               partetacut=5d0
               ktpower=1
               sbarelep=0
             elseif(srecomb(0).eq.2)then
               ktd=0.4d0
               partetacut=1d10
               ktpower=1
               sbarelep=0
             elseif(srecomb(0).eq.0)then
               write(*,'(/(a))') ' Warning: srecomb=0, '//
     &                 'recombination switched off'
               write(*,'(a)') ' use with care!'
             else
               write(*,'(/(a))') ' Error message: HAWK Stops here!'
               write(*,*) 'The input srecomb can only be set to'//
     &                    ' 0,1,2'
               write(*,*) 'You try to use srecomb=',srecomb(0)
               stop
             endif
           endif
         elseif(identifier.eq.'scuts') then
           read(storeinput,'(i100)',err=61) scuts(0)
           read(storeinput,*) scuts(0)
           if(selprocess.eq.0)then
             if(scuts(0).eq.1.or.scuts(0).eq.2)then
               ptcutj=20d0
               ycutj=4.5d0
               dycutjj=4d0
               lhemicut=.true.
               ptcutj1=20d0
               ptcutj2=20d0
               ycutj1=4.5d0
               ycutj2=4.5d0
               yminj1=0d0
               yminj2=0d0
               yminj3=0d0
               ptmaxj1=1d30
               ptmaxj2=1d30
               ptmaxj3=1d30
               ptcutvis=0d0
               ptmaxvis=1d30
               dphicutj1h=0d0
               dphicutj2h=0d0
               mlcutjj=0d0
               mucutjj=1d30
               mlcutjjh=0d0
               mucutjjh=1d30
               ptcuth=0d0
               ycuth=1d10
               ecuth=0d0
               mlcuth=1d-2             ! lower cut on Higgs invariant mass
               mucuth=1d30             ! upper cut on Higgs invariant mass
             elseif(scuts(0).eq.0.or.scuts(0).eq.-1) then
               ptcutj=0d0
               ycutj=1d10
               dycutjj=0d0
               lhemicut=.false.
               ptcutj1=0d0
               ptcutj2=0d0
               ycutj1=1d10
               ycutj2=1d10
               yminj1=0d0
               yminj2=0d0
               yminj3=0d0
               ptmaxj1=1d30
               ptmaxj2=1d30
               ptmaxj3=1d30
               ptcutvis=0d0
               ptmaxvis=1d30
               dphicutj1h=0d0
               dphicutj2h=0d0
               mlcutjj=0d0
               mucutjj=1d30
               mlcutjjh=0d0
               mucutjjh=1d30
               ptcuth=0d0
               ycuth=1d10
               ecuth=0d0
               mlcuth=1d-2             ! lower cut on Higgs invariant mass
               mucuth=1d30             ! upper cut on Higgs invariant mass
             else
               write(*,'(/(a))') ' Error message: HAWK Stops here!'
               write(*,*) 'The input scuts can only be set to'//
     &                    ' 0,1,2'
               write(*,*) 'You try to use scuts',scuts(0)
               stop
             endif
           else if(selprocess.ge.1.and.selprocess.le.4)then
              if(scuts(0).eq.0.or.scuts(0).eq.-1)then
                 ycutl=1d10
                 ptcutl=0d0
                 ptmaxl=1d10
                 yminl=0d0
              elseif(scuts(0).eq.1)then
                 ycutl=2.5d0
                 ptcutl=20d0
                 ptmaxl=1d10
                 yminl=0d0
              elseif(scuts(0).eq.2)then
                 ycutl=1d10
                 ptcutl=0d0
                 ptmaxl=20d0
                 yminl=2.5d0
              else
                 write(*,'(/(a))') 'Error message: HAWK Stops here!'
                 write(*,*) 'The input scuts can only be set to'//
     &                ' 0,1,2'
                 write(*,*) 'You try to use scuts',scuts(0)
                 stop
              endif
              ecutl=0d0
              if(scuts(0).eq.0.or.scuts(0).eq.-1)then
                 ptcutv=0d0
                 ptcutm=0d0
                 ptcuth=0d0
              else
                 ptcutv=190d0
                 ptcutm=25d0
                 ptcuth=200d0
              endif
              ycuth=1d10
              ecuth=0d0
              mlcuth=1d-2       ! lower cut on Higgs invariant mass
              mucuth=1d30       ! upper cut on Higgs invariant mass
              mlcutll=0d0
              mucutll=1d10
              mlcutllh=0d0
              mucutllh=1d10
              mcutjl=0d0
              dRcutll=0d0
              dRcutjl=0d0
              dycutll=0d0
              dycutjl=0d0
           endif
         elseif(identifier.eq.'sincludecuts') then
           read(storeinput,'(i100)',err=61) sincludecuts(0)
           read(storeinput,*) sincludecuts(0)
           if(sincludecuts(0).lt.0.or.sincludecuts(0).gt.2) then
             write(*,*) 'The input sincludecuts can only be set to'//
     &                    ' 0,1,2'
             write(*,*) 'You try to use sincludecuts =',sincludecuts(0)
             write(*,*) 'Your input is replaced by the default: 1'
             sincludecuts(0)=1
           endif
         elseif(identifier.eq.'sirtr') then
           read(storeinput,'(i100)',err=61) sirtr
           read(storeinput,*) sirtr
           if(sirtr.eq.0)then
             ssub(0)=1
             ssli(0)=0
           elseif(sirtr.eq.1)then
             ssub(0)=0
             ssli(0)=1
           else
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sirtr can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sirtr',sirtr
             stop
           endif
         elseif(identifier.eq.'sendp') then
           read(storeinput,'(i100)',err=61) sendp(0)
         elseif(identifier.eq.'ssub') then
           read(storeinput,'(i100)',err=61) ssub(0)
           read(storeinput,*) ssub(0)
         elseif(identifier.eq.'ssli') then
           read(storeinput,'(i100)',err=61) ssli(0)
           read(storeinput,*) ssli(0)
         elseif(identifier.eq.'seff') then
           read(storeinput,'(i100)',err=61) seff(0)
           read(storeinput,*) seff(0)
         elseif(identifier.eq.'sborn') then
           read(storeinput,'(i100)',err=61) sborn(0)
           read(storeinput,*) sborn(0)
           if(sborn(0).ne.0.and.sborn(0).ne.1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sborn can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sborn',sborn(0)
             stop
           endif
         elseif(identifier.eq.'sew') then
           read(storeinput,'(i100)',err=61) sew
           read(storeinput,*) sew
           if(sew.eq.1)then
             svirt(0)=1
             sweak=1
#ifdef RECOLA
             shh2=0
#else
             shh2=1
#endif
             spbrems(0)=1
             spbrems(1)=1
             spinc(0)=1
             spinc(3)=1
             spinc(4)=1
           elseif(sew.eq.0)then
             sweak=0
             shh2=0
             spbrems(0)=0
             spbrems(1)=0
             spinc(0)=0
             spinc(3)=0
             spinc(4)=0
           else
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sew can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sew',sew
             stop
           endif
         elseif(identifier.eq.'sqcd') then
           read(storeinput,'(i100)',err=61) sqcd
           read(storeinput,*) sqcd
           if(sqcd.eq.1)then
             svirt(0)=1
             sqcddiag=1
             sqcdnondiag=1
             sqcdggfus=1
             sqcdgsplit=1
             sgbrems(0)=1
             sgbrems(2)=1
             sginc(0)=1
             sginc(5)=1
             sginc(6)=1
           elseif(sqcd.eq.0)then
             sqcddiag=0
             sqcdnondiag=0
             sqcdggfus=0
             sqcdgsplit=0
             sgbrems(0)=0
             sgbrems(2)=0
             sginc(0)=0
             sginc(5)=0
             sginc(6)=0
           else
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sqcd can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sqcd',sqcd
             stop
           endif
         elseif(identifier.eq.'shtr') then
           read(storeinput,'(i100)',err=61) shtr
           read(storeinput,*) shtr
#ifdef RECOLA
           if(shtr.ne.0)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'Off-shell Higgs-boson not supported by RECOLA.'
             stop
           end if
#endif
           if(shtr.lt.0.or.shtr.ge.3)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input shtr can only be set to'//
     &                    ' 0,1,2'
             write(*,*) 'You try to use shtr',shtr
             stop
           endif
         elseif(identifier.eq.'sgh') then
           read(storeinput,'(i100)',err=61) sgh
           read(storeinput,*) sgh
           if(sgh.lt.0.or.sgh.ge.3)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sgh can only be set to'//
     &                    ' 0,1,2'
             write(*,*) 'You try to use sgh',sgh
             stop
           endif
         elseif(identifier.eq.'shbw') then
           read(storeinput,'(i100)',err=61) shbw
           read(storeinput,*) shbw
           if (shbw.ne.0.and.shbw.ne.1.and.shbw.ne.2
     &         .and.shbw.ne.8.and.shbw.ne.9
     &         )then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input shbw can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use shbw',shbw
             stop
           endif
           if (shbw.eq.1) then
             sgh=2
           endif
         elseif(identifier.eq.'sonsborn') then
           read(storeinput,'(i100)',err=61) sonsborn
           read(storeinput,*) sonsborn
           if(sonsborn.eq.1) then
             lonsborn=.true.
           elseif(sonsborn.eq.0) then
             lonsborn=.false.
           else
             write(*,*)  'input not known: sonsborn = ',sonsborn
           endif
         elseif(identifier.eq.'sonsqcd') then
           read(storeinput,'(i100)',err=61) sonsqcd
           read(storeinput,*) sonsqcd
           if(sonsqcd.eq.1) then
             lonsqcd=.true.
           elseif(sonsqcd.eq.0) then
             lonsqcd=.false.
           else
             write(*,*)  'input not known: sonsqcd = ',sonsqcd
           endif
         elseif(identifier.eq.'sonsproj') then
           read(storeinput,'(i100)',err=61) sonsproj
           read(storeinput,*) sonsproj
           if(sonsproj.ne.0.and.sonsproj.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sonsproj can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sonsproj',sonsproj
             stop
           endif
         elseif(identifier.eq.'svirt') then
           read(storeinput,'(i100)',err=61) svirt(0)
           read(storeinput,*) svirt(0)
           if(svirt(0).ne.0.and.svirt(0).ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input svirt can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use svirt',svirt
             stop
           endif
         elseif(identifier.eq.'sweak') then
           read(storeinput,'(i100)',err=61) sweak
           read(storeinput,*) sweak
           if(sweak.ne.0.and.sweak.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sweak can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sweak',sweak
             stop
           endif
         elseif(identifier.eq.'sqcddiag') then
           read(storeinput,'(i100)',err=61) sqcddiag
           read(storeinput,*) sqcddiag
           if(sqcddiag.ne.0.and.sqcddiag.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sqcddiag can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sqcddiag',sqcddiag
             stop
           endif
           if(sqcddiag.eq.1.and.sqcd.eq.0) then
             sqcddiag=0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'sqcddiag=1 is only possible for sqcd=1'
             write(*,*) 'Your input sqcddiag=1 is ignored!'
           endif
         elseif(identifier.eq.'sqcdnondiag') then
           read(storeinput,'(i100)',err=61) sqcdnondiag
           read(storeinput,*) sqcdnondiag
           if(sqcdnondiag.ne.0.and.sqcdnondiag.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sqcdnondiag can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sqcdnondiag',sqcdnondiag

             stop
           endif
           if(sqcdnondiag.eq.1.and.sqcd.eq.0) then
             sqcdnondiag=0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'sqcdnondiag=1 is only possible for sqcd=1'
             write(*,*) 'Your input sqcdnondiag=1 is ignored!'
           endif
         elseif(identifier.eq.'sqcdggfus') then
           read(storeinput,'(i100)',err=61) sqcdggfus
           read(storeinput,*) sqcdggfus
           if(sqcdggfus.ne.0.and.sqcdggfus.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sqcdggfus can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sqcdggfus',sqcdggfus

             stop
           endif
           if(sqcdggfus.eq.1.and.sqcd.eq.0) then
             sqcdggfus=0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'sqcdggfus=1 is only possible for sqcd=1'
             write(*,*) 'Your input sqcdggfus=1 is ignored!'
           endif
         elseif(identifier.eq.'sqcdgsplit') then
           read(storeinput,'(i100)',err=61) sqcdgsplit
           read(storeinput,*) sqcdgsplit
           if(sqcdgsplit.ne.0.and.sqcdgsplit.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sqcdgsplit can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sqcdgsplit',sqcdgsplit
             stop
           endif
           if(sqcdgsplit.eq.1.and.sqcd.eq.0) then
             sqcdgsplit=0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'sqcdgsplit=1 is only possible for sqcd=1'
             write(*,*) 'Your input sqcdgsplit=1 is ignored!'
           endif
         elseif(identifier.eq.'sscha') then
           read(storeinput,'(i100)',err=61) sscha
           read(storeinput,*) sscha
#ifdef RECOLA
           if(sscha.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'sscha!=1 selection not supported by RECOLA.'
             stop
           end if
#endif
           if(sscha.ne.0.and.sscha.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sscha can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sscha',sscha
             stop
           endif
         elseif(identifier.eq.'stcha') then
           read(storeinput,'(i100)',err=61) stcha
           read(storeinput,*) stcha
#ifdef RECOLA
           if(stcha.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'stcha!=1 selection not supported by RECOLA.'
             stop
           end if
#endif
           if(stcha.ne.0.and.stcha.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input stcha can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use stcha',stcha
             stop
           endif
         elseif(identifier.eq.'sch2') then
           read(storeinput,'(i100)',err=61) sch2
           read(storeinput,*) sch2
#ifdef RECOLA
           if(sch2.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'sch2!=1 selection not supported by RECOLA.'
             stop
           end if
#endif
           if(sch2.ne.0.and.sch2.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sch2 can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sch2',sch2
             stop
           endif
         elseif(identifier.eq.'schint') then
           read(storeinput,'(i100)',err=61) schint
           read(storeinput,*) schint
#ifdef RECOLA
           if(schint.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'schint!=1 selection not supported by RECOLA.'
             stop
           end if
#endif
           if(schint.ne.0.and.schint.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input schint can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use schint',schint
             stop
           endif
         elseif(identifier.eq.'sbini') then
           read(storeinput,'(i100)',err=61) sbini
           read(storeinput,*) sbini
           if(sbini.ne.0.and.sbini.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sbini can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sbini',sbini
             stop
           endif
         elseif(identifier.eq.'sbfin') then
           read(storeinput,'(i100)',err=61) sbfin
           read(storeinput,*) sbfin
           if(sbfin.ne.0.and.sbfin.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sbfin can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sbfin',sbfin
             stop
           endif
         elseif(identifier.eq.'shh2') then
           read(storeinput,'(i100)',err=61) shh2
           read(storeinput,*) shh2
#ifdef RECOLA
           if (shh2 .eq. 1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*)
     &         '2-loop heavy Higgs corrections not supported in RECOLA.'
             stop
           end if
#endif
           if(shh2.ne.0.and.shh2.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input shh2 can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use shh2',shh2
             stop
           endif
         elseif(identifier.eq.'scp') then
           read(storeinput,'(i100)',err=61) scp
           read(storeinput,*) scp
           if(scp.ne.0.and.scp.ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input scp can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use scp',scp
             stop
           endif
         elseif(identifier.eq.'spdf') then
           read(storeinput,'(i100)',err=61) spdf
           read(storeinput,*) spdf
           if(spdf.ne.99.and.spdf.ne.1.and.
     &        spdf.ne.11.and.spdf.ne.14)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input spdf can only be set to 1,11,14,99'
             write(*,*) 'You try to use spdf=',spdf
             stop
           endif
           if(spdf.eq.1) then
             Open(unit=99,file='qed6-10gridp.dat',status='old',err=55)
             Close(99)
           else if(spdf.eq.11) then
             Open(unit=99,file='cteq6m.tbl',status='old',err=56)
             Close(99)
           else if(spdf.eq.14) then
             Open(unit=99,file='cteq6l1.tbl',status='old',err=57)
             Close(99)
           endif
         elseif(identifier.eq.'pdfmap') then
           read(storeinput,'(i100)',err=61) pdfmap
           read(storeinput,*) pdfmap
           if(pdfmap.ne.4.and.pdfmap.ne.3.and.pdfmap.ne.2
     &       .and.pdfmap.ne.1.and.pdfmap.ne.0)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input pdfmap can only be set to'//
     &                    ' 0,1,2,3,4'
             write(*,*) 'You try to use pdfmap',pdfmap
             stop
           endif
         elseif(identifier.eq.'pdfname') then
           read(storeinput,'(a100)',err=61) pdfname
           read(storeinput,*) pdfname
           pdfpath=''


         elseif(identifier.eq.'pdfpath') then
           read(storeinput,'(a100)',err=61) pdfpath
           read(storeinput,*) pdfpath
           pdfname=''


         elseif(identifier.eq.'pdfmember') then
           read(storeinput,'(i100)',err=61) pdfmember
           read(storeinput,*) pdfmember
         elseif(identifier.eq.'pdfmemberfrom') then
#ifdef multpdf
           read(storeinput,'(i100)',err=61) pdfmemberfrom
           read(storeinput,*) pdfmemberfrom
#else
           write(*,'(/(a))') ' Error message: HAWK Stops here!'
           write(*,*) 'The input pdfmemberfrom can only be used'//
     &                    ' if LHAPDF is available.'
           write(*,*) 'Recompile with LHAPDF or remove '//
     &                'pdfmemberfrom from input file.'
           stop
#endif
         elseif(identifier.eq.'pdfmemberto') then
#ifdef multpdf
           read(storeinput,'(i100)',err=61) pdfmemberto
           read(storeinput,*) pdfmemberto
#else
           write(*,'(/(a))') ' Error message: HAWK Stops here!'
           write(*,*) 'The input pdfmemberto can only be used'//
     &                    ' if LHAPDF is available.'
           write(*,*) 'Recompile with LHAPDF or remove '//
     &                'pdfmemberto from input file.'
           stop
#endif
         elseif(identifier.eq.'sfactqcd') then
           read(storeinput,'(i100)',err=61) sfactqcd
           read(storeinput,*) sfactqcd
           if(sfactqcd.ne.2.and.sfactqcd.ne.1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sfactqcd can only be set to'//
     &                    ' 1,2'
             write(*,*) 'You try to use sfactqcd',sfactqcd
             stop
           endif
         elseif(identifier.eq.'sfactqed') then
           read(storeinput,'(i100)',err=61) sfactqed
           read(storeinput,*) sfactqed
           if(sfactqed.ne.2.and.sfactqed.ne.1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sfactqed can only be set to'//
     &                    ' 1,2'
             write(*,*) 'You try to use sfactqed',sfactqed
             stop
           endif
         elseif(identifier.eq.'spbrems') then
           read(storeinput,'(i100)',err=61) spbrems(0)
           read(storeinput,*) spbrems(0)
           if(spbrems(0).ne.0.and.spbrems(0).ne.1)then
             write(*,*) 'input not known: spbrems = ',spbrems(0)
             stop
           endif
           spbrems(1)=spbrems(0)
         elseif(identifier.eq.'sgbrems') then
           read(storeinput,'(i100)',err=61) sgbrems(0)
           read(storeinput,*) sgbrems(0)
           if(sgbrems(0).ne.0.and.sgbrems(0).ne.1)then
             write(*,*) 'input not known: sgbrems(0) = ',sgbrems(0)
             stop
           endif
           sgbrems(2)=sgbrems(0)
         elseif(identifier.eq.'spbrems(1)') then
           read(storeinput,'(i100)',err=61) spbrems(1)
           read(storeinput,*) spbrems(1)
           if(spbrems(1).ne.0.and.spbrems(1).ne.1)then
             write(*,*) 'input not known: spbrems(1) = ',spbrems(1)
             stop
           endif
         elseif(identifier.eq.'sgbrems(2)') then
           read(storeinput,'(i100)',err=61) sgbrems(2)
           read(storeinput,*) sgbrems(2)
           if(sgbrems(2).ne.0.and.sgbrems(2).ne.1)then
             write(*,*) 'input not known: sgbrems(2) = ',sgbrems(2)
             stop
           endif
         elseif(identifier.eq.'spinc') then
           read(storeinput,'(i100)',err=61) spinc(0)
           read(storeinput,*) spinc(0)
           if(spinc(0).ne.0.and.spinc(0).ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input spinc can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use spinc',spinc(0)
             stop
           endif
           spinc(3)=spinc(0)
           spinc(4)=spinc(0)
         elseif(identifier.eq.'sginc') then
           read(storeinput,'(i100)',err=61) sginc(0)
           read(storeinput,*) sginc(0)
           if(sginc(0).ne.0.and.sginc(0).ne.1)then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sginc can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sginc',sginc(0)
             stop
           endif
           sginc(5)=sginc(0)
           sginc(6)=sginc(0)
         elseif(identifier.eq.'lnoutmc') then
           read(storeinput,'(i100)',err=61) lnoutmc
           read(storeinput,*) lnoutmc
         elseif(identifier.eq.'maxout') then
           read(storeinput,'(i100)',err=61) maxout
           read(storeinput,*) maxout
         elseif(identifier.eq.'lnoutgen') then
           read(storeinput,'(i100)',err=61) lnoutgen
           read(storeinput,*) lnoutgen
         elseif(identifier.eq.'nhisto') then
           read(storeinput,'(i100)',err=61) nhisto
           read(storeinput,*) nhisto
           if(shisto.lt.0.or.nhisto.gt.maxh) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input nhisto can only be set to'//
     &                    ' 0<nhisto<',maxh,''
             write(*,*) 'You try to use nhisto',nhisto
             stop
           endif
         elseif(identifier.eq.'shisto') then
           read(storeinput,'(i100)',err=61) shisto
           read(storeinput,*) shisto
           if(shisto.ne.0.and.shisto.ne.1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input shisto can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use shisto',shisto
             stop
           endif
         elseif(identifier.eq.'deltas') then
           read(storeinput,'(f50.50)',err=61) mcslideltas
           read(storeinput,*) mcslideltas
         elseif(identifier.eq.'deltac') then
           read(storeinput,'(f50.50)',err=61) mcslideltac
           read(storeinput,*) mcslideltac
         elseif(identifier.eq.'deltas2') then
           read(storeinput,'(f50.50)',err=61) mcslideltas
           read(storeinput,*) mcslideltas
           mcslideltas=sqrt(mcslideltas)
         elseif(identifier.eq.'deltac2') then
           read(storeinput,'(f50.50)',err=61) mcslideltac
           read(storeinput,*) mcslideltac
           mcslideltac=sqrt(mcslideltac)
         elseif(identifier.eq.'ptcut(jet)') then
           read(storeinput,'(f50.50)',err=61) ptcutj
           read(storeinput,*) ptcutj
           ptcutj1=ptcutj
           ptcutj2=ptcutj
           if(ptcutj.gt.0d0.and.scuts(0).eq.0) then
             ptcutj=0d0
             ptcutj1=0d0
             ptcutj2=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(jet) is ignored!'
           endif
         elseif(identifier.eq.'ptcut(jet1)') then
           read(storeinput,'(f50.50)',err=61) ptcutj1
           read(storeinput,*) ptcutj1
           if(ptcutj1.gt.0d0.and.scuts(0).eq.0) then
             ptcutj1=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(jet1) is ignored!'
           endif
         elseif(identifier.eq.'ptcut(jet2)') then
           read(storeinput,'(f50.50)',err=61) ptcutj2
           read(storeinput,*) ptcutj2
           if(ptcutj2.gt.0d0.and.scuts(0).eq.0) then
             ptcutj2=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(jet2) is ignored!'
           endif
         elseif(identifier.eq.'ptmax(jet1)') then
           read(storeinput,'(f50.50)',err=61) ptmaxj1
           read(storeinput,*) ptmaxj1
           if(ptmaxj1.gt.0d0.and.scuts(0).eq.0) then
             ptmaxj1=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptmax(jet1) is ignored!'
           endif
         elseif(identifier.eq.'ptmax(jet2)') then
           read(storeinput,'(f50.50)',err=61) ptmaxj2
           read(storeinput,*) ptmaxj2
           if(ptmaxj2.gt.0d0.and.scuts(0).eq.0) then
             ptmaxj2=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptmax(jet2) is ignored!'
           endif
         elseif(identifier.eq.'ptmax(jet3)') then
           read(storeinput,'(f50.50)',err=61) ptmaxj3
           read(storeinput,*) ptmaxj3
           if(ptmaxj3.gt.0d0.and.scuts(0).eq.0) then
             ptmaxj3=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptmax(jet3) is ignored!'
           endif
         elseif(identifier.eq.'ycut(jet)') then
           read(storeinput,'(f50.50)',err=61) ycutj
           read(storeinput,*) ycutj
           ycutj1=ycutj
           ycutj2=ycutj
           if(ycutj.lt.1d10.and.scuts(0).eq.0) then
             ycutj=1d10
             ycutj1=1d10
             ycutj2=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ycut(jet) is ignored!'
           endif
         elseif(identifier.eq.'ycut(jet1)') then
           read(storeinput,'(f50.50)',err=61) ycutj1
           read(storeinput,*) ycutj1
           if(ycutj1.lt.1d10.and.scuts(0).eq.0) then
             ycutj1=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ycut(jet1) is ignored!'
           endif
         elseif(identifier.eq.'ycut(jet2)') then
           read(storeinput,'(f50.50)',err=61) ycutj2
           read(storeinput,*) ycutj2
           if(ycutj2.lt.1d10.and.scuts(0).eq.0) then
             ycutj2=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ycut(jet2) is ignored!'
           endif
         elseif(identifier.eq.'ymin(jet1)') then
           read(storeinput,'(f50.50)',err=61) yminj1
           read(storeinput,*) yminj1
           if(yminj1.gt.0d0.and.scuts(0).eq.0) then
             yminj1=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ymin(jet1) is ignored!'
           endif
         elseif(identifier.eq.'ymin(jet2)') then
           read(storeinput,'(f50.50)',err=61) yminj2
           read(storeinput,*) yminj2
           if(yminj2.gt.0d0.and.scuts(0).eq.0) then
             yminj2=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ymin(jet2) is ignored!'
           endif
         elseif(identifier.eq.'ymin(jet3)') then
           read(storeinput,'(f50.50)',err=61) yminj3
           read(storeinput,*) yminj3
           if(yminj3.gt.0d0.and.scuts(0).eq.0) then
             yminj3=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ymin(jet3) is ignored!'
           endif
         elseif(identifier.eq.'dycut(jet,jet)') then
           read(storeinput,'(f50.50)',err=61) dycutjj
           read(storeinput,*) dycutjj
           if(dycutjj.gt.0d0.and.scuts(0).eq.0) then
             dycutjj=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for dycut(jet,jet) is ignored!'
           endif
         elseif(identifier.eq.'mcut(jet,jet)') then
           read(storeinput,'(f50.50)',err=61) mlcutjj
           read(storeinput,*) mlcutjj
           if(mlcutjj.gt.0d0.and.scuts(0).eq.0) then
             mlcutjj=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mcut(jet,jet) is ignored!'
           endif
         elseif(identifier.eq.'mlcut(jet,jet)') then
           read(storeinput,'(f50.50)',err=61) mlcutjj
           read(storeinput,*) mlcutjj
           if(mlcutjj.gt.0d0.and.scuts(0).eq.0) then
             mlcutjj=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mlcut(jet,jet) is ignored!'
           endif
         elseif(identifier.eq.'mucut(jet,jet)') then
           read(storeinput,'(f50.50)',err=61) mucutjj
           read(storeinput,*) mucutjj
           if(mucutjj.gt.0d0.and.scuts(0).eq.0) then
             mucutjj=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mucut(jet,jet) is ignored!'
           endif
         elseif(identifier.eq.'mlcut(jet,jet,higgs)'.or.
     &         identifier.eq.'mlcut(higgs,jet,jet)') then
           read(storeinput,'(f50.50)',err=61) mlcutjjh
           read(storeinput,*) mlcutjjh
           if(mlcutjjh.gt.0d0.and.scuts(0).eq.0) then
             mlcutjjh=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mlcut(jet,jet,higgs) is ignored!'
           endif
         elseif(identifier.eq.'mucut(jet,jet,higgs)'.or.
     &         identifier.eq.'mucut(higgs,jet,jet)') then
           read(storeinput,'(f50.50)',err=61) mucutjjh
           read(storeinput,*) mucutjjh
           if(mucutjjh.gt.0d0.and.scuts(0).eq.0) then
             mucutjjh=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mucut(jet,jet,higgs) is ignored!'
           endif
         elseif(identifier.eq.'ptcut(visible)') then
           read(storeinput,'(f50.50)',err=61) ptcutvis
           read(storeinput,*) ptcutvis
           if(ptcutvis.gt.0d0.and.scuts(0).eq.0) then
             ptcutvis=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(visible) is ignored!'
           endif
         elseif(identifier.eq.'ptmax(visible)') then
           read(storeinput,'(f50.50)',err=61) ptmaxvis
           read(storeinput,*) ptmaxvis
           if(ptmaxvis.lt.1d10.and.scuts(0).eq.0) then
             ptmaxvis=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptmax(visible) is ignored!'
           endif
         elseif(identifier.eq.'dphicut(jet1,higgs)'
     &         .or.identifier.eq.'dphicut(higgs,jet1)') then
           read(storeinput,'(f50.50)',err=61) dphicutj1h
           read(storeinput,*) dphicutj1h
           if(dphicutj1h.gt.0d0.and.scuts(0).eq.0) then
             dphicutj1h=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for dphicut(higgs,jet1) is ignored!'
           endif
         elseif(identifier.eq.'dphicut(jet2,higgs)'
     &         .or.identifier.eq.'dphicut(higgs,jet2)') then
           read(storeinput,'(f50.50)',err=61) dphicutj2h
           read(storeinput,*) dphicutj2h
           if(dphicutj2h.gt.0d0.and.scuts(0).eq.0) then
             dphicutj2h=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for dphicut(jet2,higgs) is ignored!'
           endif
         elseif(identifier.eq.'hemispherecut') then
           read(storeinput,'(i100)',err=61) hemicut
           read(storeinput,*) hemicut
           if(hemicut.eq.1) then
             lhemicut=.true.
           elseif(hemicut.eq.0) then
             lhemicut=.false.
           else
             goto 61
           endif
           if(hemicut.ne.0.and.scuts(0).eq.0) then
             hemicut=0
             lhemicut=.false.
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for hemispherecut is ignored!'
           endif
         elseif(identifier.eq.'higgsbetweenjets') then
           read(storeinput,'(i100)',err=61) higgsbjets
           read(storeinput,*) higgsbjets
           if(higgsbjets.eq.1) then
             lhiggsbjets=.true.
           elseif(higgsbjets.eq.0) then
             lhiggsbjets=.false.
           else
             goto 61
           endif
           if(higgsbjets.ne.0.and.scuts(0).eq.0) then
             higgsbjets=0
             lhiggsbjets=.false.
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for higgsbetweenjets is ignored!'
           endif
         elseif(identifier.eq.'ptcut(higgs)') then
           read(storeinput,'(f50.50)',err=61) ptcuth
           read(storeinput,*) ptcuth
           if(ptcuth.gt.0d0.and.scuts(0).eq.0) then
             ptcuth=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(higgs) is ignored!'
           endif
         elseif(identifier.eq.'ycut(higgs)') then
           read(storeinput,'(f50.50)',err=61) ycuth
           read(storeinput,*) ycuth
           if(ycuth.lt.1d10.and.scuts(0).eq.0) then
             ycuth=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ycut(higgs) is ignored!'
           endif
         elseif(identifier.eq.'ecut(higgs)') then
           read(storeinput,'(f50.50)',err=61) ecuth
           read(storeinput,*) ecuth
           if(ecuth.gt.0d0.and.scuts(0).eq.0) then
             ecuth=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ecut(higgs) is ignored!'
           endif
         elseif(identifier.eq.'mucut(higgs)') then
           read(storeinput,'(f50.50)',err=61) mucuth
           read(storeinput,*) mucuth
         elseif(identifier.eq.'mlcut(higgs)') then
           read(storeinput,'(f50.50)',err=61) mlcuth
           read(storeinput,*) mlcuth
         elseif(identifier.eq.'ptcut(decp)') then
           read(storeinput,'(f50.50)',err=61) ptcuts
           read(storeinput,*) ptcuts
           if(ptcuts.gt.0d0.and.scuts(0).eq.0) then
             ptcuts=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(decp) is ignored!'
           endif
         elseif(identifier.eq.'ycut(decp)') then
           read(storeinput,'(f50.50)',err=61) ycuts
           read(storeinput,*) ycuts
           if(ycuts.lt.1d10.and.scuts(0).eq.0) then
             ycuts=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptmax(visible) is ignored!'
           endif
         elseif(identifier.eq.'ktpower') then
           read(storeinput,'(i100)',err=61) ktpower
           read(storeinput,*) ktpower
         elseif(identifier.eq.'sbarelep') then
           read(storeinput,'(i100)',err=61) sbarelep
           read(storeinput,*) sbarelep
           if(sbarelep.ne.0.and.sbarelep.ne.1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input sbarelep can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use sbarelep',sbarelep
             stop
           endif
         elseif(identifier.eq.'colllep') then
           read(storeinput,'(i100)',err=61) colllep
           read(storeinput,*) colllep
           if(colllep.ne.0.and.colllep.ne.1) then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input colllep can only be set to'//
     &                    ' 0,1'
             write(*,*) 'You try to use colllep',colllep
             stop
           endif
         elseif(identifier.eq.'dparameter') then
           read(storeinput,'(f50.50)',err=61) ktd
           read(storeinput,*) ktd
         elseif(identifier.eq.'etacut(parton)') then
           read(storeinput,'(f50.50)',err=61) partetacut
           read(storeinput,*) partetacut

         elseif(identifier.eq.'ptcut(lep)') then
           read(storeinput,'(f50.50)',err=61) ptcutl
           read(storeinput,*) ptcutl
           if(ptcutl.gt.0d0.and.scuts(0).eq.0) then
             ptcutl=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(lep) is ignored!'
           endif
         elseif(identifier.eq.'ycut(lep)') then
           read(storeinput,'(f50.50)',err=61) ycutl
           read(storeinput,*) ycutl
           if(ycutl.lt.1d10.and.scuts(0).eq.0) then
             ycutl=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ycut(lep) is ignored!'
           endif
         elseif(identifier.eq.'ecut(lep)') then
           read(storeinput,'(f50.50)',err=61) ecutl
           read(storeinput,*) ecutl
           if(ecutl.gt.0d0.and.scuts(0).eq.0) then
             ecutl=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ecut(lep) is ignored!'
           endif
         elseif(identifier.eq.'mucut(lep,lep)') then
           read(storeinput,'(f50.50)',err=61) mucutll
           read(storeinput,*) mucutll
           if(mucutll.lt.1d10.and.scuts(0).eq.0) then
             mucutll=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mucut(lep,lep) is ignored!'
           endif
         elseif(identifier.eq.'mlcut(lep,lep)') then
           read(storeinput,'(f50.50)',err=61) mlcutll
           read(storeinput,*) mlcutll
           if(mlcutll.gt.0d0.and.scuts(0).eq.0) then
             mlcutll=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mlcut(lep,lep) is ignored!'
           endif
         elseif(identifier.eq.'mucut(lep,lep,higgs)'.or.
     &         identifier.eq.'mucut(higgs,lep,lep)') then
           read(storeinput,'(f50.50)',err=61) mucutllh
           read(storeinput,*) mucutllh
           if(mucutllh.lt.1d10.and.scuts(0).eq.0) then
             mucutllh=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mucut(lep,lep,higgs) is ignored!'
           endif
         elseif(identifier.eq.'mlcut(lep,lep,higgs)'.or.
     &         identifier.eq.'mlcut(higgs,lep,lep)') then
           read(storeinput,'(f50.50)',err=61) mlcutllh
           read(storeinput,*) mlcutllh
           if(mlcutllh.gt.0d0.and.scuts(0).eq.0) then
             mlcutllh=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mlcut(lep,lep,higgs) is ignored!'
           endif
         elseif(identifier.eq.'drcut(lep,lep)') then
           read(storeinput,'(f50.50)',err=61) dRcutll
           read(storeinput,*) dRcutll
           if(dRcutll.gt.0d0.and.scuts(0).eq.0) then
             dRcutll=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for drcut(lep,lep) is ignored!'
           endif
         elseif(identifier.eq.'drcut(jet,lep)') then
           read(storeinput,'(f50.50)',err=61) dRcutjl
           read(storeinput,*) dRcutjl
           if(dRcutjl.gt.0d0.and.scuts(0).eq.0) then
             dRcutjl=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for drcut(jet,lep) is ignored!'
           endif
         elseif(identifier.eq.'dgammaparameter') then
           read(storeinput,'(f50.50)',err=61) dgammaparameter
           read(storeinput,*) dgammaparameter
         elseif(identifier.eq.'ptmax(lep)') then
           read(storeinput,'(f50.50)',err=61) ptmaxl
           read(storeinput,*) ptmaxl
           if(ptmaxl.lt.1d10.and.scuts(0).eq.0) then
             ptmaxl=1d10
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptmax(lep) is ignored!'
           endif
         elseif(identifier.eq.'ymin(lep)') then
           read(storeinput,'(f50.50)',err=61) yminl
           read(storeinput,*) yminl
           if(yminl.gt.0d0.and.scuts(0).eq.0) then
             yminl=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ymin(lep) is ignored!'
           endif
         elseif(identifier.eq.'ptcut(miss)') then
           read(storeinput,'(f50.50)',err=61) ptcutm
           read(storeinput,*) ptcutm
           if(ptcutm.gt.0d0.and.scuts(0).eq.0) then
             ptcutm=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(miss) is ignored!'
           endif
         elseif(identifier.eq.'ptcut(v)') then
           read(storeinput,'(f50.50)',err=61) ptcutv
           read(storeinput,*) ptcutv
           if(ptcutv.gt.0d0.and.scuts(0).eq.0) then
             ptcutv=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(v) is ignored!'
           endif
         elseif(identifier.eq.'ptcut(V)') then
           read(storeinput,'(f50.50)',err=61) ptcutv
           read(storeinput,*) ptcutv
           if(ptcutv.gt.0d0.and.scuts(0).eq.0) then
             ptcutv=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for ptcut(V) is ignored!'
           endif
         elseif(identifier.eq.'dycut(jet,lep)') then
           read(storeinput,'(f50.50)',err=61) dycutjl
           read(storeinput,*) dycutjl
           if(dycutjl.gt.0d0.and.scuts(0).eq.0) then
             dycutjl=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for dycut(jet,lep) is ignored!'
           endif
         elseif(identifier.eq.'dycut(lep,lep)') then
           read(storeinput,'(f50.50)',err=61) dycutll
           read(storeinput,*) dycutll
           if(dycutll.gt.0d0.and.scuts(0).eq.0) then
             dycutll=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for dycut(lep,lep) is ignored!'
           endif
         elseif(identifier.eq.'mcut(jet,lep)') then
           read(storeinput,'(f50.50)',err=61) mcutjl
           read(storeinput,*) mcutjl
           if(mcutjl.gt.0d0.and.scuts(0).eq.0) then
             mcutjl=0d0
             write(*,'(/(a))') ' Warning message: HAWK input ignored!'
             write(*,*) 'scuts=0 calculates the inclusive cross section.'
             write(*,*) 'Your input for mcut(jet,lep) is ignored!'
           endif

         elseif(identifier.eq.'autofreq') then
           read(storeinput,'(i100)',err=61) autofreq
           read(storeinput,*) autofreq
         elseif(identifier.eq.'printweightmax') then
           read(storeinput,'(i100)',err=61) printweightmax
           read(storeinput,*) printweightmax
          if(printweightmax.eq.1) then
             lprintweightmax=.true.
           elseif(printweightmax.eq.0) then
             lprintweightmax=.false.
           else
             goto 61
           endif
         elseif(identifier.eq.'initsteps') then
           read(storeinput,'(i100)',err=61) initsteps
           read(storeinput,*) initsteps
         elseif(identifier.eq.'ranluxseed') then
#ifdef RANLUX
           read(storeinput,'(i100)',err=61) ranluxseed
           read(storeinput,*) ranluxseed
#endif
c input for anomalous HVV couplings
         elseif(identifier.eq.'qhvv') then
         read(storeinput,'(i100)',err=61) qhvv
         read(storeinput,*) qhvv
#ifdef RECOLA
         if (qhvv .ne. 0) then
           write(*,'(/(a))') ' Error message: HAWK Stops here!'
           write(*,*) 'anomalous HVV couplings (qhvv) not yet'//
     &                ' supported in RECOLA.'
           stop
         end if
#else
           shvv=qhvv
           if (qhvv.ne.0.and.qhvv.ne.1.and.qhvv.ne.2
     &         )then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input qhvv can only be set to'//
     &                    ' 0,1,2'
             write(*,*) 'You try to use qhvv',qhvv
             stop
           endif
#endif
         elseif(identifier.eq.'shvv') then
         read(storeinput,'(i100)',err=61) shvv
         read(storeinput,*) shvv
#ifdef RECOLA
         if (shvv .ne. 0) then
           write(*,'(/(a))') ' Error message: HAWK Stops here!'
           write(*,*) 'anomalous HVV couplings (shvv) not yet supported in RECOLA.'
           stop
         end if
#else
           qhvv=shvv
           if (shvv.ne.0.and.shvv.ne.1.and.shvv.ne.2
     &         )then
             write(*,'(/(a))') ' Error message: HAWK Stops here!'
             write(*,*) 'The input shvv can only be set to'//
     &                    ' 0,1,2'
             write(*,*) 'You try to use shvv',shvv
             stop
           endif
#endif
         elseif(identifier.eq.'rsm') then
           read(storeinput,'(f50.50)',err=61) rsm
           read(storeinput,*) rsm
         elseif(identifier.eq.'d') then
           read(storeinput,'(f50.50)',err=61) d
           read(storeinput,*) d
         elseif(identifier.eq.'db') then
           read(storeinput,'(f50.50)',err=61) db
           read(storeinput,*) db
         elseif(identifier.eq.'dt') then
           read(storeinput,'(f50.50)',err=61) dt
           read(storeinput,*) dt
         elseif(identifier.eq.'dtb') then
           read(storeinput,'(f50.50)',err=61) dtb
           read(storeinput,*) dtb

         elseif(identifier.eq.'lambdahvv') then
           read(storeinput,'(f50.50)',err=61) lambdahvv
           read(storeinput,*) lambdahvv
         elseif(identifier.eq.'a1hww') then
           read(storeinput,'(f50.50)',err=61) a1hww
           read(storeinput,*) a1hww
         elseif(identifier.eq.'a2hww') then
           read(storeinput,'(f50.50)',err=61) a2hww
           read(storeinput,*) a2hww
         elseif(identifier.eq.'a3hww') then
           read(storeinput,'(f50.50)',err=61) a3hww
           read(storeinput,*) a3hww
         elseif(identifier.eq.'a1haa') then
           read(storeinput,'(f50.50)',err=61) a1haa
           read(storeinput,*) a1haa
         elseif(identifier.eq.'a2haa') then
           read(storeinput,'(f50.50)',err=61) a2haa
           read(storeinput,*) a2haa
         elseif(identifier.eq.'a3haa') then
           read(storeinput,'(f50.50)',err=61) a3haa
           read(storeinput,*) a3haa
         elseif(identifier.eq.'a1haz') then
           read(storeinput,'(f50.50)',err=61) a1haz
           read(storeinput,*) a1haz
         elseif(identifier.eq.'a2haz') then
           read(storeinput,'(f50.50)',err=61) a2haz
           read(storeinput,*) a2haz
         elseif(identifier.eq.'a3haz') then
           read(storeinput,'(f50.50)',err=61) a3haz
           read(storeinput,*) a3haz
         elseif(identifier.eq.'a1hzz') then
           read(storeinput,'(f50.50)',err=61) a1hzz
           read(storeinput,*) a1hzz
         elseif(identifier.eq.'a2hzz') then
           read(storeinput,'(f50.50)',err=61) a2hzz
           read(storeinput,*) a2hzz
         elseif(identifier.eq.'a3hzz') then
           read(storeinput,'(f50.50)',err=61) a3hzz
           read(storeinput,*) a3hzz

         elseif(identifier.ne.' ') then
           write(*,'(/(a))') ' Error message: HAWK Stops here!'
           write(*,*) 'The following Identifier in the input'//
     &                ' file is not known:  ',identifier
           write(*,*) 'The following Identifier in the input'//
     &                ' file is not known:  ',orgidentifier
           stop
         endif

         goto 50
       else
         directive= inputline(1:isep2-1)
         orgdirective=directive
         do i1=isep2-1,1,-1
           if (directive(i1:i1).eq.' ') then
             directive=directive(1:i1-1)//directive(i1+1:linputline)
     &           //' '
           elseif(directive(i1:i1).ge.'A'.and.directive(i1:i1).le.'Z')
     &         then
             if(directive(i1:i1).eq.'A') then
               directive(i1:i1)='a'
             elseif(directive(i1:i1).eq.'B') then
               directive(i1:i1)='b'
             elseif(directive(i1:i1).eq.'C') then
               directive(i1:i1)='c'
             elseif(directive(i1:i1).eq.'D') then
               directive(i1:i1)='d'
             elseif(directive(i1:i1).eq.'E') then
               directive(i1:i1)='e'
             elseif(directive(i1:i1).eq.'F') then
               directive(i1:i1)='f'
             elseif(directive(i1:i1).eq.'G') then
               directive(i1:i1)='g'
             elseif(directive(i1:i1).eq.'H') then
               directive(i1:i1)='h'
             elseif(directive(i1:i1).eq.'I') then
               directive(i1:i1)='i'
             elseif(directive(i1:i1).eq.'J') then
               directive(i1:i1)='j'
             elseif(directive(i1:i1).eq.'K') then
               directive(i1:i1)='k'
             elseif(directive(i1:i1).eq.'L') then
               directive(i1:i1)='l'
             elseif(directive(i1:i1).eq.'M') then
               directive(i1:i1)='m'
             elseif(directive(i1:i1).eq.'N') then
               directive(i1:i1)='n'
             elseif(directive(i1:i1).eq.'O') then
               directive(i1:i1)='o'
             elseif(directive(i1:i1).eq.'P') then
               directive(i1:i1)='p'
             elseif(directive(i1:i1).eq.'Q') then
               directive(i1:i1)='q'
             elseif(directive(i1:i1).eq.'R') then
               directive(i1:i1)='r'
             elseif(directive(i1:i1).eq.'S') then
               directive(i1:i1)='s'
             elseif(directive(i1:i1).eq.'T') then
               directive(i1:i1)='t'
             elseif(directive(i1:i1).eq.'U') then
               directive(i1:i1)='u'
             elseif(directive(i1:i1).eq.'V') then
               directive(i1:i1)='v'
             elseif(directive(i1:i1).eq.'W') then
               directive(i1:i1)='w'
             elseif(directive(i1:i1).eq.'X') then
               directive(i1:i1)='x'
             elseif(directive(i1:i1).eq.'Y') then
               directive(i1:i1)='y'
             elseif(directive(i1:i1).eq.'Z') then
               directive(i1:i1)='z'
             endif
           endif
         enddo
         write(storeinput,*) inputline(1:isep2-1)


         if(directive.eq.'tev') then
           energy=1960d0
           sppbar=1
         elseif(directive.eq.'lhc7') then
           energy=7000d0
           sppbar=0
         elseif(directive.eq.'lhc14') then
           energy=14000d0
           sppbar=0
         elseif(directive.eq.'qcdonly') then
           sborn(0)=1
           svirt(0)=1
           sweak=0
           shh2=0
           sqcddiag=1
           sqcdnondiag=1
           sqcdggfus=1
           sqcdgsplit=1
           spbrems(0)=0
           sgbrems(0)=1
           spbrems(1)=0
           sgbrems(2)=1
           spinc(0)=0
           sginc(0)=1
           spinc(3)=0
           spinc(4)=0
           sginc(5)=1
           sginc(6)=1
         elseif(directive.eq.'bornonly') then
           sborn(0)=1
           svirt(0)=0
           sweak=0
           shh2=0
           sqcddiag=0
           sqcdnondiag=0
           sqcdggfus=0
           sqcdgsplit=0
           spbrems(0)=0
           sgbrems(0)=0
           spbrems(1)=0
           sgbrems(2)=0
           spinc(0)=0
           sginc(0)=0
           spinc(3)=0
           spinc(4)=0
           sginc(5)=0
           sginc(6)=0
           ssub(0)=0
           ssli(0)=0
           seff(0)=0
         elseif(directive.eq.'ewonly') then
           sborn(0)=1
           svirt(0)=1
           sweak=1
           sqcddiag=0
           sqcdnondiag=0
           sqcdggfus=0
           sqcdgsplit=0
           spbrems(0)=1
           sgbrems(0)=0
           spbrems(1)=1
           sgbrems(2)=0
           spinc(0)=1
           sginc(0)=0
           spinc(3)=1
           spinc(4)=1
           sginc(5)=0
           sginc(6)=0
         elseif(directive.eq.'endinput') then
           goto 51
         elseif(directive.ne.' ') then
           write(*,'(/(a))') ' Error message: HAWK Stops here!'
           write(*,*) 'The following directive in the input'//
     &                ' file is not known:  ',directive
           write(*,*) 'The following directive in the input'//
     &                ' file is not known:  ',orgdirective
           stop
         endif
         goto 50
       end if

 55    write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'You have to copy the file qed6-10gridp.dat '//
     &     'from the HAWK'
       write(*,*) 'subdirectory PDFs to the working directory!'
       write(*,*) 'Then you can restart HAWK.'
       stop
 56    write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'You have to copy the file cteq6m.tbl '//
     &     'from the HAWK'
       write(*,*) 'subdirectory PDFs to the working directory!'
       write(*,*) 'Then you can restart HAWK.'
       stop
 57    write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'You have to copy the file cteq6l1.tbl '//
     &     'from the HAWK'
       write(*,*) 'subdirectory PDFs to the working directory!'
       write(*,*) 'Then you can restart HAWK.'
       stop

 52    write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'The following input line in the input'//
     &                ' file is not appropriate:  ',inputline
       stop
 61    write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'In the input file the input is not '//
     &            'appropriate for identifier ''',
     &     identifier(1:isep1-1),''' : ',storeinput
       stop

 51    continue
#ifdef mpiuse
       close(5)
#endif

#endif


ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  put scuts(0)=1 for scuts(0)=-1 initialization  c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (scuts(0).eq.-1) scuts(0)=1

ccccccccccccccccccccccccccccccccccccccccc
c  calculate CKM matrix from sinthetac  c
ccccccccccccccccccccccccccccccccccccccccc
      vin(1,2) = sinthetacin
      vin(1,1) = sqrt(1d0-vin(1,2)*vin(1,2))
      vin(2,1) = -vin(1,2)
      vin(2,2) = vin(1,1)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  collinear lepton emission off incoming photon is         c
c  automatically switched on if ptcut(lep) is set to zero.  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccsk      if (ptcutl.eq.0d0) colllep=3
      if (ptcutl.eq.0d0) colllep=2

cccccccccccccccccccccccccccccccccccccccccc
c  colllep = 1: HAWK selects dipole      c
cccccccccccccccccccccccccccccccccccccccccc
c            if (colllep.eq.1) colllep=2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  redefinition of pT_miss and pT_V in case of missing lepton  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc       if(ptcutl.eq.0d0.and.(selprocess.eq.1.or.selprocess.eq.2)) then
       if ((ptcutl.eq.0d0).and.
     &      (selprocess.eq.1.or.selprocess.eq.2)) then
          if (ptcutm.lt.ptcutv) then
             ptcutm=ptcutv
          endif
          ptcutv=0d0
       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  set switches that don't have a meaning in WH/ZH to  c
c  adaequate values.                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (selprocess.ge.1) then
          sqcdnondiag = 0
          sqcdggfus = 0
          sqcdgsplit = 0
          if (sqcd.eq.1) then
             sqcddiag = 1
          else
             sqcddiag = 0
          endif
          sscha = 1
          stcha = 1
          sch2 = 1
          schint = 1
          shh2 = 0                 ! irrelevant
       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  set switches that don't have a meaning in VBF to    c
c  adaequate values.                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (selprocess.eq.0) then
          sbarelep=0
       endif



c check for inconsistent switches
       if(srecomb(0).eq.0.and.scuts(0).gt.0) then
         write(*,'(/(a))') ' Error message: HAWK Stops here!'
         write(*,'(a49,i4)')' Inconsistent input: '//
     &       'srecomb = 0 and 0 < scuts = ',scuts(0)
         stop
       endif

c check for inconsistent switches
       if(shtr.gt.0.and.sew.ne.0.and.(sonsproj.eq.0.or..not.lonsew))
     &     then
         write(*,'(/(a))') ' Error message: HAWK Stops here!'
         write(*,'(a44,i2,/a50)') ' Inconsistent input: '//
     &       'sonsproj = 0 and sew = ',sew,
     &       'electroweak corrections need on-shell projection'
         stop
       endif

       if(shtr.gt.0.and.shbw.eq.1.and.sgh.eq.0)
     &     then
         write(*,'(/(a))') ' Error message: HAWK Stops here!'
         write(*,'(a44,i2,/a70/a20)') ' Inconsistent input: '//
     &       'shbw = 1 and sgh = ',sgh,
     &       'input Higgs width for propagator with off-shell width '//
     &       'not allowed',
     &       'choose sgh=2!'
         stop
       endif

       if(shtr.gt.0.and.shbw.eq.1.and.sgh.eq.1)
     &     then
         write(*,'(/(a))') ' Error message: HAWK Stops here!'
         write(*,'(a44,i2,/a70/a50)') ' Warning: '//
     &       'shbw = 1 and sgh = ',sgh,
     &       'off-shell propagator with on-shell width! ',
     &       'use complex mass scheme! (sgh=2)'
         stop
       endif

       if(shtr.eq.2.and.selprocess.gt.0)
     &     then
         write(*,'(/(a))') ' Error message: HAWK Stops here!'
         write(*,'(a44,i2,/a70/a50)') ' Inconsistent input: '//
     &       'shtr = 2 and selprocess = ',selprocess,
     &       'Higgs decay into two singlets not supported for ',
     &       'associated production'
         stop
       endif

       if(shvv.ne.0.and.sqcdnondiag.ne.0)
     &     then
         write(*,'(a17,i2,a19,i2,a13/a60/a26/)') ' Warning: '//
     &       'shvv = ',shvv,' and sqcdnondiag = ',sqcdnondiag,
     &       ' not allowed!',
     &       'sqcdnondiag with anomalous higgs couplings'//
     &        ' not supported ',
     &       'sqcdnondiag set to zero'
         sqcdnondiag=0
       endif

       if(shvv.ne.0.and.sqcdggfus.ne.0)
     &     then
         write(*,'(a17,i2,a17,i2,a13/a58/a24/)') ' Warning: '//
     &       'shvv = ',shvv,' and sqcdggfus = ',sqcdggfus,
     &       ' not allowed!',
     &       'sqcdggfus with anomalous Higgs couplings'//
     &        ' not supported ',
     &       'sqcdggfus set to zero'
         sqcdggfus=0
       endif

       if(shvv.ne.0.and.sqcdgsplit.ne.0)
     &     then
         write(*,'(a17,i2,a18,i2,a13/a59/a25/)') ' Warning: '//
     &       'shvv = ',shvv,' and sqcdgsplit = ',sqcdgsplit,
     &       ' not allowed!',
     &       'sqcdgsplit with anomalous Higgs couplings'//
     &        ' not supported ',
     &       'sqcdgsplit set to zero'
         sqcdgsplit=0
       endif










       nsteps = max(nevents/100,10000)

c initialize pdfs
       if(sppbar.eq.0)then
         shadin1=1
         shadin2=1
       elseif(sppbar.eq.1)then
         shadin1=1
         shadin2=-1
       elseif(sppbar.eq.3)then
         shadin1=-1
         shadin2=1
       elseif(sppbar.eq.2)then
         shadin1=-1
         shadin2=-1
       else
         write(*,'(a44,i2,/a50)') ' Unknown input value: '//
     &       'sppbar = ',sppbar
         stop
       endif


      if (spdf.ge.10.and.spdf.lt.20) then
        spinc(0)=0
        spinc(3)=0
        spinc(4)=0
        Call SetCtq6(spdf-10)
      else if (spdf.eq.99) then
#ifdef includeLHAPDF
        if (pdfpath.ne.'') then
#ifndef mpiuse
          write(*,*) 'pdfpath = ',pdfpath
#endif
          call InitPDFset(pdfpath)
        elseif(pdfname.ne.'') then
#ifndef mpiuse
          write(*,*) 'pdfname = ',pdfname
#endif
          call InitPDFsetByName(pdfname)
        else
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'no pdf set specified for LHApdf in input file!'
          write(*,*) 'Use pdfname or pdfpath!'
          stop
        endif
        call numberPDF(num)
        If(pdfmember.LT.0.OR.pdfmember.GT.num) Then
          write(*,*) 'pdfmember out of range: set to         0'
          pdfmember=0
        Endif
#ifdef multpdf
        pdfmembercentral=pdfmember
#endif
        call InitPDF(pdfmember)

         if (has_photon().eqv..FALSE.) then
            spinc(0)=0
            spinc(3)=0
            spinc(4)=0
          write(*,*) 'Chosen pdf set does not contain photon pdfs:'
          write(*,*) 'Incoming-photon contributions are set to 0.'
          write(*,*)
         endif

#else
       write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'LHAPDF not included: recompile with LHAPDF to use '//
     &            'LHAPDF interface, i.e. spdf=99!'
       write(*,*) 'To run HAWK without LHAPDF (spdf=1,11,14),'//
     &            ' change spdf in the input file'
       Stop
#endif
      endif

      qalps=spdf

c contributions for generator 0 (Born generator)
      sgen(0)=1
c contribution 1: Born
      if (sborn(0).eq.1) sgencon(0,1)=1
      swborn(0)=1
c contribution 2: electroweak virtual - end-point part
      if ((svirt(0).eq.1.or.sendp(0).eq.1).and.
     &    (sweak.ne.0.or.shh2.ne.0))
     &    sgencon(0,2)=1
c contribution 3: singular parts of photonic real and virtual corrections
      if ((sweak.eq.1.or.sweak.eq.3).and.spbrems(0).ne.0.and.
     &    (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0))then
        sgencon(0,3)=1
      endif
c contribution 4: singular parts of gluonic real and virtual corrections
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.sgbrems(0).ne.0.and.
     &     (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0)) then
        sgencon(0,4)=1
      endif
c contribution 5: qcd virtual - end-point part
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.
     &     (svirt(0).eq.1.or.sendp(0).eq.1)) then
         sgencon(0,5)=1
      endif
c contribution 6: incoming photon contribution 1
      if (spinc(0).ne.0.and.
     &    (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0)) then
        sgencon(0,6)=1
      endif
c contribution 7: incoming photon contribution 2
      if (spinc(0).ne.0.and.
     &    (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0)) then
        sgencon(0,7)=1
      endif
c contribution 8: incoming gluon contribution 1
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.sginc(0).ne.0.and.
     &     (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0)) then
         sgencon(0,8)=1
      endif
c contribution 9: incoming gluon contribution 2
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.sginc(0).ne.0.and.
     &     (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0)) then
        sgencon(0,9)=1
      endif
c contribution 10: collinear-photon from lepton (muon) contribution
      if(sbarelep.ne.0.and.sew.ne.0.and.
     &    (ssub(0).gt.0.or.ssli(0).gt.0.or.seff(0).gt.0)) then
        sgencon(0,10)=1
      endif


c contributions for generator 1 (real photon radiation)
c contribution 1: real corrections - subtraction terms
      if (spbrems(1).ne.0) sgen(1)=1

c contributions for generator 2 (real gluon radiation)
c contribution 1: real corrections - subtraction terms
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.sgbrems(2).ne.0) then
         sgen(2)=1
      endif

c contributions for generator 3 (incoming photon contribution 1)
c contribution 1: real corrections - subtraction terms
      if (spinc(3).ne.0) sgen(3)=1

c contributions for generator 4 (incoming photon contribution 2)
c contribution 1: real corrections - subtraction terms
      if (spinc(4).ne.0) sgen(4)=1

c contributions for generator 5 (incoming gluon contribution 1)
c contribution 1: real corrections - subtraction terms
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.sginc(5).ne.0) then
         sgen(5)=1
      endif

c contributions for generator 6 (incoming gluon contribution 2)
c contribution 1: real corrections - subtraction terms
      if ((sqcddiag.eq.1.or.sqcdnondiag.eq.1).and.sginc(6).ne.0) then
         sgen(6)=1
      endif
c contributions for generator 7 (incoming l- contribution 1)
c contribution 1: real corrections - subtraction terms
      if (spinc(0).ne.0.and.ptcutl.eq.0d0.and.
     &     (selprocess.eq.1.or.selprocess.eq.3)) then
         sgen(7)=1
      endif

c contributions for generator 8 (incoming l- contribution 2)
c contribution 1: real corrections - subtraction terms
      if (spinc(0).ne.0.and.ptcutl.eq.0d0.and.
     &     (selprocess.eq.1.or.selprocess.eq.3)) then
         sgen(8)=1
      endif

c contributions for generator 9 (incoming l+ contribution 1)
c contribution 1: real corrections - subtraction terms
      if (spinc(0).ne.0.and.ptcutl.eq.0d0.and.
     &     (selprocess.eq.2.or.selprocess.eq.3)) then
         sgen(9)=1
      endif

c contributions for generator 10 (incoming l+ contribution 2)
c contribution 1: real corrections - subtraction terms
      if (spinc(0).ne.0.and.ptcutl.eq.0d0.and.
     &     (selprocess.eq.2.or.selprocess.eq.3)) then
         sgen(10)=1
      endif


      if (sincludecuts(0).gt.2.or.sincludecuts(0).lt.0) then
        write(*,*)
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) ' sincludecuts > 2 or < 0 not allowed'
        write(*,*) ' sincludecuts = ',sincludecuts(0)
        stop
      endif

      do gen=1,maxg
        sgencon(gen,1)=sgen(gen)
        ssub(gen)=ssub(0)
        ssli(gen)=ssli(0)
        seff(gen)=seff(0)
        sincludecuts(gen)=sincludecuts(0)
      enddo

      if (ssli(0).ne.0.and.ssub(0).ne.0.or.
     &    ssli(0).ne.0.and.seff(0).ne.0.or.
     &    seff(0).ne.0.and.ssub(0).ne.0) then
        write(*,*)
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) ' either slicing or subtraction or '
     &      ,'effective collinear factors must be used'
        write(*,*) ' ssli = ',ssli(1)
        write(*,*) ' ssub = ',ssub(1)
        write(*,*) ' seff = ',seff(1)
        stop
      endif


      if(autofreq.eq.1)then
        if(sppbar.eq.0.or.sppbar.eq.2)then
          if(sgencon(0,1).eq.1.or.sgencon(0,4).eq.1) freqgencon(0,2)=100
          if(sgencon(0,1).eq.1.or.sgencon(0,4).eq.1) freqgencon(0,5)=10
          if(sgencon(0,4).eq.1) freqgencon(0,3)=10
          if(sgencon(0,8).eq.1) freqgencon(0,6)=10
          if(sgencon(0,9).eq.1) freqgencon(0,7)=10
          if(sgencon(0,10).eq.1) freqgencon(0,10)=10
          if(sgencon(2,1).eq.1) freqgencon(1,1)=10
          if(sgencon(5,1).eq.1) freqgencon(3,1)=10
          if(sgencon(6,1).eq.1) freqgencon(4,1)=10
          if(sgencon(5,1).eq.1) freqgencon(7,1)=10
          if(sgencon(6,1).eq.1) freqgencon(8,1)=10
          if(sgencon(5,1).eq.1) freqgencon(9,1)=10
          if(sgencon(6,1).eq.1) freqgencon(10,1)=10
        elseif(sppbar.eq.1.or.sppbar.eq.3)then
          if(sgencon(0,1).eq.1.or.sgencon(0,4).eq.1) freqgencon(0,2)=200
          if(sgencon(0,1).eq.1.or.sgencon(0,4).eq.1) freqgencon(0,5)=10
          if(sgencon(0,4).eq.1) freqgencon(0,3)=20
          if(sgencon(0,8).eq.1) freqgencon(0,6)=50
          if(sgencon(0,9).eq.1) freqgencon(0,7)=50
          if(sgencon(0,10).eq.1) freqgencon(0,10)=20
          if(sgencon(2,1).eq.1) freqgencon(1,1)=20
          if(sgencon(5,1).eq.1) freqgencon(3,1)=10
          if(sgencon(6,1).eq.1) freqgencon(4,1)=10
          if(sgencon(5,1).eq.1) freqgencon(7,1)=10
          if(sgencon(6,1).eq.1) freqgencon(8,1)=10
          if(sgencon(5,1).eq.1) freqgencon(9,1)=10
          if(sgencon(6,1).eq.1) freqgencon(10,1)=10
       endif

       if (sqcd.eq.0.and.sew.eq.0) then
          freqgencon(0,6)=1
          freqgencon(0,7)=1
          freqgencon(3,1)=1
          freqgencon(4,1)=1
          freqgencon(7,1)=1
          freqgencon(8,1)=1
          freqgencon(9,1)=1
          freqgencon(10,1)=1
       endif

c for subtraction
      elseif(autofreq.eq.2)then
        if(sgencon(0,1).eq.1) then
          freqgencon(0,2)=1000
          freqgencon(0,3)=100
          freqgencon(0,4)=5
          freqgencon(0,5)=100
          freqgencon(0,6)=100
          freqgencon(0,7)=100
          freqgencon(0,8)=10
          freqgencon(0,9)=10
          freqgencon(0,10)=100
          freqgencon(1,1)=100
          freqgencon(2,1)=5
          freqgencon(3,1)=100
          freqgencon(4,1)=100
          freqgencon(5,1)=10
          freqgencon(6,1)=10
          freqgencon(7,1)=100
          freqgencon(8,1)=100
          freqgencon(9,1)=100
          freqgencon(10,1)=100
        else if(
     &      sgencon(2,1).eq.1.or.sgencon(0,4).eq.1)then
          freqgencon(0,2)=200
          freqgencon(0,3)=20
          freqgencon(0,4)=1
          freqgencon(0,5)=20
          freqgencon(0,6)=20
          freqgencon(0,7)=20
          freqgencon(0,8)=2
          freqgencon(0,9)=2
          freqgencon(0,10)=20
          freqgencon(1,1)=20
          freqgencon(2,1)=1
          freqgencon(3,1)=20
          freqgencon(4,1)=20
          freqgencon(5,1)=2
          freqgencon(6,1)=2
          freqgencon(7,1)=20
          freqgencon(8,1)=20
          freqgencon(9,1)=20
          freqgencon(10,1)=20
        elseif(sgencon(0,8).eq.1.or.sgencon(0,9).eq.1.or.
     &    sgencon(5,1).eq.1.or.sgencon(6,1).eq.1)then
          freqgencon(0,2)=100
          freqgencon(0,3)=10
          freqgencon(0,5)=10
          freqgencon(0,6)=10
          freqgencon(0,7)=10
          freqgencon(0,8)=1
          freqgencon(0,9)=1
          freqgencon(0,10)=10
          freqgencon(1,1)=10
          freqgencon(3,1)=10
          freqgencon(4,1)=10
          freqgencon(5,1)=1
          freqgencon(6,1)=1
          freqgencon(7,1)=10
          freqgencon(8,1)=10
          freqgencon(9,1)=10
          freqgencon(10,1)=10
        elseif(sgencon(0,2).eq.1.or.sgencon(0,3).eq.1.or.
     &        sgencon(0,5).eq.1.or.
     &        sgencon(0,8).eq.1.or.sgencon(0,9).eq.1.or.
     &        sgencon(3,1).eq.1.or.sgencon(4,1).eq.1) then
          freqgencon(0,2)=10
          freqgencon(0,3)=1
          freqgencon(0,5)=1
          freqgencon(0,6)=1
          freqgencon(0,7)=1
          freqgencon(0,10)=1
          freqgencon(1,1)=1
          freqgencon(3,1)=1
          freqgencon(4,1)=1
          freqgencon(7,1)=1
          freqgencon(8,1)=1
          freqgencon(9,1)=1
          freqgencon(10,1)=1
        endif
c for slicing
      elseif(autofreq.eq.3)then
        if(sgencon(0,1).eq.1.or.
     &      sgencon(2,1).eq.1.or.sgencon(0,4).eq.1)then
          freqgencon(0,2)=10000
          freqgencon(0,3)=100
          freqgencon(0,4)=1
          freqgencon(0,5)=1000
          freqgencon(0,6)=10000
          freqgencon(0,7)=10000
          freqgencon(0,8)=100
          freqgencon(0,9)=100
          freqgencon(0,10)=100
          freqgencon(1,1)=10
          freqgencon(2,1)=1
          freqgencon(3,1)=1000
          freqgencon(4,1)=1000
          freqgencon(5,1)=10
          freqgencon(6,1)=10
          freqgencon(7,1)=1000
          freqgencon(8,1)=1000
          freqgencon(9,1)=1000
          freqgencon(10,1)=1000
        elseif(sgencon(1,1).eq.1.or.
     &    sgencon(5,1).eq.1.or.sgencon(6,1).eq.1)then
          freqgencon(0,2)=10
          freqgencon(0,2)=1000
          freqgencon(0,3)=10
          freqgencon(0,5)=100
          freqgencon(0,6)=1000
          freqgencon(0,7)=1000
          freqgencon(0,8)=10
          freqgencon(0,9)=10
          freqgencon(0,10)=10
          freqgencon(1,1)=1
          freqgencon(3,1)=100
          freqgencon(4,1)=100
          freqgencon(5,1)=1
          freqgencon(6,1)=1
          freqgencon(7,1)=100
          freqgencon(8,1)=100
          freqgencon(9,1)=100
          freqgencon(10,1)=100
        elseif(sgencon(0,2).eq.1.or.sgencon(0,3).eq.1.or.
     &    sgencon(0,8).eq.1.or.sgencon(0,9).eq.1)then
          freqgencon(0,2)=1
          freqgencon(0,2)=100
          freqgencon(0,3)=1
          freqgencon(0,5)=10
          freqgencon(0,6)=100
          freqgencon(0,7)=100
          freqgencon(0,8)=1
          freqgencon(0,9)=1
          freqgencon(0,10)=1
          freqgencon(3,1)=10
          freqgencon(4,1)=10
          freqgencon(7,1)=10
          freqgencon(8,1)=10
          freqgencon(9,1)=10
          freqgencon(10,1)=10
        elseif(sgencon(0,5).eq.1.or.
     &    sgencon(3,1).eq.1.or.sgencon(4,1).eq.1) then
          freqgencon(0,2)=10
          freqgencon(0,5)=1
          freqgencon(0,6)=10
          freqgencon(0,7)=10
          freqgencon(3,1)=1
          freqgencon(4,1)=1
          freqgencon(7,1)=10
          freqgencon(8,1)=10
          freqgencon(9,1)=10
          freqgencon(10,1)=10
        endif
      endif



      if(lprintweightmax) then
        open(unit=20,file='dat.weightmax',status='unknown')
      endif

      if(lnoutmc.lt.3) lhisto(0)=0


      if(sqcddiag.eq.0.and.sqcdnondiag.eq.0.and.sgen(2).eq.0.and.
     &    sgen(5).eq.0.and.sgen(6).eq.0)then
        if(shtr.eq.0)then
          smodel=1
        else
          smodel=15
        endif
        if(spinc(4).eq.1)then
          ngenerator=4
        elseif(spinc(3).eq.1)then
          ngenerator=3
        else
          ngenerator=1
        endif
      elseif(sqcd.eq.0.and.sew.eq.0) then

        if(shtr.eq.0)then
          smodel=1
        else
          smodel=15
        endif
        if (selprocess.eq.0) then
           if(spinc(4).eq.1)then
              ngenerator=4
           elseif(spinc(3).eq.1)then
              ngenerator=3
           else
              ngenerator=1
           endif
        endif

      else
        if(shtr.eq.0)then
          smodel=12
        else
          smodel=125
        endif
        if (selprocess.eq.0) then
           if(sginc(6).eq.1)then
              ngenerator=6
           elseif(sginc(5).eq.1)then
              ngenerator=5
           else
              ngenerator=2
           endif
        endif
      endif

      if (selprocess.gt.0) ngenerator=10


c equal parameters for all generators
      do gen=0,ngenerator
        scuts(gen)=scuts(0)
        srecomb(gen)=srecomb(0)
        recparam(1,gen)=recparam(1,0)
        recparam(2,gen)=recparam(2,0)
        recparam(3,gen)=recparam(3,0)
        recparam(4,gen)=recparam(4,0)
        sconv(gen)=sconv(0)
        if(sppbar.eq.0)then
          name(1,gen)='proton'
        elseif(sppbar.eq.1)then
          name(1,gen)='anti-proton'
        endif
        name(2,gen)='proton'

        if (selprocess.eq.0)then
           name(3,gen)='jet'
           name(4,gen)='jet'
        elseif (selprocess.eq.1) then
           if(sbarelep.eq.0) then
              name(3,gen)='anti-e'
              name(4,gen)='nu_e'
           else
              name(3,gen)='anti-mu'
              name(4,gen)='nu_mu'
           endif
        elseif (selprocess.eq.2) then
           if(sbarelep.eq.0) then
              name(3,gen)='anti-nu_e'
              name(4,gen)='e'
           else
              name(3,gen)='anti-nu_mu'
              name(4,gen)='mu'
           endif
        elseif (selprocess.eq.3) then
           if(sbarelep.eq.0) then
              name(3,gen)='anti-e'
              name(4,gen)='e'
           else
              name(3,gen)='anti-mu'
              name(4,gen)='mu'
           endif
        elseif (selprocess.eq.4) then
           name(3,gen)='anti-nu_e'
           name(4,gen)='nu_e'
        endif

        name(5,gen)='h0'
        name(6,gen)='photon'
      enddo

c fixing external particles of bremsstrahlung generators
      if(sgen(1).eq.1)then
        gen=1
        nextgen(gen)=nextgen(0)+1
        name(nextgen(gen),gen)='photon'
        nextmax=max(nextmax,nextgen(gen))
      endif
      if(sgen(2).eq.1)then
        gen=2
        nextgen(gen)=nextgen(0)+1
        nextmax=max(nextmax,nextgen(gen))
        name(nextgen(gen),gen)='gluon'
      endif
c fixing external particles of incoming photon generators
      if(sgen(3).eq.1)then
        gen=3
        nextgen(gen)=nextgen(0)+1
        nextmax=max(nextmax,nextgen(gen))
        name(nextgen(gen),gen)='jet'
        name(1,gen)='photon'
      endif
      if(sgen(4).eq.1)then
        gen=4
        nextgen(gen)=nextgen(0)+1
        nextmax=max(nextmax,nextgen(gen))
        name(nextgen(gen),gen)='jet'
        name(2,gen)='photon'
      endif
c fixing external particles of incoming gluon generators
      if(sgen(5).eq.1)then
        gen=5
        nextgen(gen)=nextgen(0)+1
        nextmax=max(nextmax,nextgen(gen))
        name(nextgen(gen),gen)='jet'
        name(1,gen)='gluon'
      endif
      if(sgen(6).eq.1)then
        gen=6
        nextgen(gen)=nextgen(0)+1
        nextmax=max(nextmax,nextgen(gen))
        name(nextgen(gen),gen)='jet'
        name(2,gen)='gluon'
      endif
c fixing external particles of incoming lepton generators
cccsk maybe split between incoming l- (7,8) and incoming l+ (9,10) generators ???
      if(sgen(7).eq.1)then
        gen=7
        nextgen(gen)=nextgen(0)
        nextmax=max(nextmax,nextgen(gen))
        if (selprocess.eq.1) then
           if(sbarelep.eq.0) then
              name(1,gen)='e'
              name(3,gen)='nu_e'
           elseif(sbarelep.eq.1) then
              name(1,gen)='mu'
              name(3,gen)='nu_mu'
           endif
           name(4,gen)='jet'
        elseif (selprocess.eq.3) then
           if(sbarelep.eq.0) then
              name(1,gen)='e'
              name(3,gen)='e'
           elseif(sbarelep.eq.1) then
              name(1,gen)='mu'
              name(3,gen)='mu'
           endif
           name(4,gen)='jet'
        endif
      endif
      if(sgen(8).eq.1)then
        gen=8
        nextgen(gen)=nextgen(0)
        nextmax=max(nextmax,nextgen(gen))
        if (selprocess.eq.1) then
           if(sbarelep.eq.0) then
              name(2,gen)='e'
              name(4,gen)='nu_e'
           elseif(sbarelep.eq.1) then
              name(2,gen)='mu'
              name(4,gen)='nu_mu'
           endif
        elseif (selprocess.eq.3) then
           if(sbarelep.eq.0) then
              name(2,gen)='e'
              name(4,gen)='e'
           elseif(sbarelep.eq.1) then
              name(2,gen)='mu'
              name(4,gen)='mu'
           endif
        endif
        name(3,gen)='jet'
      endif
      if(sgen(9).eq.1)then
        gen=9
        nextgen(gen)=nextgen(0)
        nextmax=max(nextmax,nextgen(gen))
        if (selprocess.eq.2) then
           if(sbarelep.eq.0) then
              name(1,gen)='anti-e'
              name(3,gen)='anti-nu_e'
           elseif(sbarelep.eq.1) then
              name(1,gen)='anti-mu'
              name(3,gen)='anti-nu_mu'
           endif
        elseif (selprocess.eq.3) then
           if(sbarelep.eq.0) then
              name(1,gen)='anti-e'
              name(3,gen)='anti-e'
           elseif(sbarelep.eq.1) then
              name(1,gen)='anti-mu'
              name(3,gen)='anti-mu'
           endif
        endif
        name(4,gen)='jet'
      endif
      if(sgen(10).eq.1)then
        gen=10
        nextgen(gen)=nextgen(0)
        nextmax=max(nextmax,nextgen(gen))
        if (selprocess.eq.2) then
           if(sbarelep.eq.0) then
              name(2,gen)='anti-e'
              name(4,gen)='anti-nu_e'
           elseif(sbarelep.eq.1) then
              name(2,gen)='anti-mu'
              name(4,gen)='anti-nu_mu'
           endif
        elseif (selprocess.eq.3) then
           if(sbarelep.eq.0) then
              name(2,gen)='anti-e'
              name(4,gen)='anti-e'
           elseif(sbarelep.eq.1) then
              name(2,gen)='anti-mu'
              name(4,gen)='anti-mu'
           endif
        endif
        name(3,gen)='jet'
      endif

c process
#ifndef RECOLA2
      if(selprocess.eq.0) then
        finalstate='j j H (VBF) '
      elseif(selprocess.eq.1) then
        finalstate='l+ nu H     '
      elseif(selprocess.eq.2) then
        finalstate='nu~ l- H    '
      elseif(selprocess.eq.3) then
        finalstate='l+ l- H     '
      elseif(selprocess.eq.4) then
        finalstate='nu~ nu H    '
      endif
#else
      if(selprocess.eq.0) then
        finalstate='j j ' // trim(higgsflavor) //' (VBF) '
      elseif(selprocess.eq.1) then
        finalstate='l+ nu ' // trim(higgsflavor) // '     '
      elseif(selprocess.eq.2) then
        finalstate='nu~ l- ' // trim(higgsflavor) // '    '
      elseif(selprocess.eq.3) then
        finalstate='l+ l- ' // trim(higgsflavor) // '     '
      elseif(selprocess.eq.4) then
        finalstate='nu~ nu ' // trim(higgsflavor) // '    '
      endif
#endif

      if(sppbar.eq.0) then
        initialstate='p p  '
      elseif(sppbar.eq.1) then
        initialstate='p p~ '
      elseif(sppbar.eq.2) then
        initialstate='p~ p '
      elseif(sppbar.eq.3) then
        initialstate='p~ p~'
      endif


c nextsub: charged external particles:  p p , j j for born

      nextsub(0)=4
      do gen=0,ngenerator
        nextme(gen)=nextgen(gen)
        nextsub(gen)=nextsub(0)
      enddo

      if(shtr.gt.0)then
        do gen=0,ngenerator
          nextgen(gen)=nextgen(gen)+1
          if(gen.eq.0)then
            name(5,gen)='s0'
            name(nextgen(gen),gen)='s0'
          else
            name(nextgen(gen),gen)= name(nextgen(gen)-1,gen)
            name(5,gen)='s0'
            name(6,gen)='s0'
          endif
        enddo
        nextmax=nextmax+1
      endif




c switches for contributions of real generators
      do gen=1,ngenerator
        sgencon(gen,1)=sgen(gen)
      enddo


c equal recombination parameters for all generators
      do gen=0,ngenerator
        recparam(1,gen)=recparam(1,0)
        recparam(2,gen)=recparam(2,0)
        recparam(3,gen)=recparam(3,0)
        recparam(4,gen)=recparam(4,0)
      enddo


#ifdef RECOLA
      if (selprocess .eq. 0) then
        if (sinthetacin .ne. 0) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'Non-zero Cabibbo angle not supported in VBF + RECOLA.'
          stop
        end if
        if (sqcd .eq. 1 .and. sqcddiag .ne. 1) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcddiag selection not supported in VBF + RECOLA.'
          stop
        end if
        if (sqcd .eq. 1 .and. sqcdnondiag .ne. 1) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcdnondiag selection not supported in VBF + RECOLA.'
          stop
        end if
        if (sqcd .eq. 1 .and. sqcdggfus .ne. 1) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcdggfus selection not supported in VBF + RECOLA.'
          stop
        end if
        if (sqcd .eq. 1 .and. sqcdgsplit .ne. 1) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcdgsplit selection not supported in VBF + RECOLA.'
          stop
        end if
        if (sqcd .eq. 0 .and. sqcddiag .ne. 0) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcddiag selection not supported in VBF + RECOLA.'
          write(*,*) 'Enable all qcd corrections via sqcd=1.'
          stop
        end if
        if (sqcd .eq. 0 .and. sqcdnondiag .ne. 0) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcdnondiag selection not supported in VBF + RECOLA.'
          write(*,*) 'Enable all qcd corrections via sqcd=1.'
          stop
        end if
        if (sqcd .eq. 0 .and. sqcdggfus .ne. 0) then
          write(*,'(/(a))') ' Error message: HAWK input ignored!'
          write(*,*) 'sqcdggfus selection not supported in VBF + RECOLA.'
          write(*,*) 'Enable all qcd corrections via sqcd=1.'
          stop
        end if
        if (sqcd .eq. 0 .and. sqcdgsplit .ne. 0) then
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'sqcdgsplit selection not supported in VBF + RECOLA.'
          write(*,*) 'Enable all qcd corrections via sqcd=1.'
          stop
        end if
      end if
#endif

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     setting input parameters                                    c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
***********************************************************************
      subroutine inputparameter
***********************************************************************
*     setting defaults for global input parameters                    *
*---------------------------------------------------------------------*
*     24.07.07  Ansgar Denner     last changed  24.07.07              *
***********************************************************************
#ifdef RECOLA2
        use recola, only: get_modelname_rcl
#endif
      implicit none
#ifdef RECOLA2
        character*100 :: modelname
#endif
#include "mcinput.h"

#ifdef RECOLA2
        call get_modelname_rcl(modelname)
#endif

      alpha0in  = 1d0/137.035999139d0
      alphazin  = 1d0/128.936d0
      GFin      = 1.1663787D-5

      mzin      = 91.1876d0
      gzin      = 2.4952d0
      mwin      = 80.385d0
      gwin      = 2.085d0
      mhin      = 125d0
      ghin      = 0.004088d0

c      mein   = .51099892d-3
c      mmuin  = .1056583715d0
c      mtauin = 1.77699d0
c      mdin  = .066d0
c      muin  = .066d0
c      msin  = .150d0
c      mcin  = 1.2d0
c      mbin  = 4.3d0
      mtin  = 172.5d0

c     mein   = 1d-2
c     mdin   = 1d-2
c     muin   = 1d-3

      mein   = 1d-3
      mmuin  = .1056583715d0
      mtauin = 1d-3
      mdin  = 1d-3
      muin  = 1d-3
      msin  = 1d-3
      mcin  = 1d-3
      mbin  = 1d-3



      sinthetacin = 0.225d0

      vin(1,2) = sinthetacin
      vin(1,2) = sqrt(1d0-vin(1,2)*vin(1,2))
      vin(2,1) = -vin(1,2)
      vin(2,2) = vin(1,1)

      vin(1,3) = 0d0
      vin(2,3) = 0d0
      vin(3,1) = 0d0
      vin(3,2) = 0d0
      vin(3,3) = 1d0

#ifdef RECOLA2
      if (modelname .eq. "THDM") then
        cabin = 0.1d0
        tbin = 2d0
        l5in = -1.9d0
        ytypein = 2
        mhhin = 300d0
        mhcin = 460d0
        mhain = 460d0
        ghhin = 0d0
        ghcin = 0d0
        ghain = 0d0
        mixrsin = "OS12"
        mumsbsmin = (mhin + mhhin + 2*mhcin + mhain)/5
      else if (modelname .eq. "HS") then
        sain = 0.29d0
        l3in = 0.14d0
        mhhin = 200d0
        mixrsin = "OS"
        mumsbsmin = mhin
      end if
#endif

      lambdain=mwin
      mudim2in=mwin*mwin

      gevtofbin=389379660D3

      end
***********************************************************************
      subroutine fixparameter
***********************************************************************
*     setting dependent parameters                                    *
*---------------------------------------------------------------------*
*     24.07.07  Ansgar Denner     last changed  24.07.07              *
***********************************************************************
      implicit none
      real*8   wfac,zfac
      real*8   alphasf

      integer iNumProcs, iErr, iMyRank
      common/mpistuff/iNumProcs,iMyRank

#include "mcinput.h"
#include "mcparams.h"

c mcoutput        in generator and montecarlo
      integer nout,numout,maxout,noutgen,nsteps,lnoutmc
      common/mcoutput/nout,numout,maxout,noutgen,nsteps,lnoutmc

c qcdrenscale
      real*8 qcdrenscale,qcdrenscale2,qcdrenscalefac
      common /qcdrenscale/ qcdrenscale,qcdrenscale2,qcdrenscalefac

c qedfacscale
      real*8 qedfacscale,qedfacscale2,qedfacscalefac
      common /qedfacscale/ qedfacscale,qedfacscale2,qedfacscalefac

c qcdfacscale
      real*8 qcdfacscale,qcdfacscale2,qcdfacscalefac
      common /qcdfacscale/ qcdfacscale,qcdfacscale2,qcdfacscalefac

c mcpdf         in public and montecarlo
      real*8 pdf1(-6:6),pdf2(-6:6),pdf_phot1,pdf_phot2
      integer spdf,pdfmap
      common/mcpdf/pdf1,pdf2,pdf_phot1,pdf_phot2,spdf,pdfmap
c mcoptions  in public and montecarlo, initamps
      integer salp
      common/mcrs/salp

      integer nf
      real*8  lambdalo,lambdanlo
      common/qcd/lambdalo,lambdanlo,nf

c LHAPDF
      integer nfmax,num,order
      real*8  alphasPDF,muf,qcdl5

      integer selprocess
      common/selprocess/selprocess

#ifdef multpdf
#include "multpdf.h"
      integer nset
#endif

c transition to pole definition
      wfac = sqrt(1d0+(gwin/mwin)**2)
      mwin   = mwin/wfac
      gwin   = gwin/wfac
      zfac = sqrt(1d0+(gzin/mzin)**2)
      mzin   = mzin/zfac
      gzin   = gzin/zfac


C*** renormalization scheme
        if (salp.eq.0) then
          alphain = alpha0in
        elseif (salp.eq.1) then
          alphain = alphazin
        elseif (salp.eq.2) then
          alphain = sqrt(2d0)*gfin/pi*mwin*mwin
     &        *(1d0-mwin*mwin/(mzin*mzin))
        elseif (salp.eq.3) then
          alphain = alphazin
        endif

c initialize running alphas
      if (spdf.eq.1.or.spdf.eq.2) then
        alphaszin=0.1187d0
        call fixlambdas(mzin,alphaszin,2,5)
      elseif(spdf.eq.11) then
        alphaszin=0.118d0
        call fixlambdas(mzin,alphaszin,2,5)
      elseif(spdf.eq.14) then
        alphaszin=0.130d0
        call fixlambdas(mzin,alphaszin,1,5)
      elseif(spdf.eq.99) then
#ifdef includeLHAPDF
        alphaszin=alphasPDF(mzin)
        call numberPDF(num)
        if(lnoutmc.ge.3)then
          If(iMyrank.EQ.0) write(nout,*) 'members in pdf set:',num
Cam removed for compatability with LHAPDF 6.1.5 (functions no longer provided)
C          call GetLam5(1,qcdl5)
C          If(iMyrank.EQ.0) write(nout,*) 'lambda5_qcd:',qcdl5
C          call GetOrderPDF(order)
C          If(iMyrank.EQ.0) write(nout,*) 'order of pdf set:',order
C          call GetOrderAS(order)
C          If(iMyrank.EQ.0) write(nout,*)
C     &     'order of alphas evol:',order
C          call GetRenFac(muf)
C          If(iMyrank.EQ.0) write(nout,*)
C     &     'renormalization factor:',muf
          call GetNF(nfmax)
          If(iMyrank.EQ.0) write(nout,*) 'number of flavors:',nfmax
        endif
#else
       write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'LHAPDF not included: recompile with LHAPDF to use '//
     &            'LHAPDF interface, i.e. spdf=99!'
       write(*,*) 'To run HAWK without LHAPDF (spdf=1,11,14),'//
     &            ' change spdf in the input file'
       Stop
#endif
      endif





      if (selprocess.eq.0) then
         call getqcdrenscale(mwin)
         call getqcdfacscale(mwin)
         call getqedfacscale(mwin)
      elseif (selprocess.eq.1.or.selprocess.eq.2) then
         call getqcdrenscale(mwin+mhin)
         call getqcdfacscale(mwin+mhin)
         call getqedfacscale(mwin+mhin)
      elseif (selprocess.eq.3.or.selprocess.eq.4) then
         call getqcdrenscale(mzin+mhin)
         call getqcdfacscale(mzin+mhin)
         call getqedfacscale(mzin+mhin)
      endif


      if (spdf.eq.1.or.spdf.eq.2.or.spdf.eq.11) then
        alphasin=alphasf(qcdrenscale,2)
      elseif(spdf.eq.14) then
        alphasin=alphasf(qcdrenscale,1)
      elseif(spdf.eq.99) then
#ifdef includeLHAPDF
#ifdef multpdf
        If(pdfmemberfrom.LT.0.OR.pdfmemberfrom.GT.num) Then
          write(nout,*) 'pdfmemberfrom = ',pdfmemberfrom
          write(nout,*) 'pdfmemberfrom out of range: set to         0'
          pdfmemberfrom=0
        Endif
        If(pdfmemberto.LT.0.OR.pdfmemberto.GT.num) Then
          write(nout,*) 'pdfmemberto = ',pdfmemberto
          write(nout,*) 'pdfmemberto out of range: set to ',num
          pdfmemberto=num
        Endif
        if(lnoutmc.ge.2)then
        If(pdfmemberto.LT.pdfmemberfrom) Then
          write(nout,*) 'pdfmemberto smaller than pdfmemberfrom'
          write(nout,*) 'no additional PDF sets will be calculated'
        Endif
        endif
        If(pdfmemberto-pdfmemberfrom+1.GT.maxpdf) Then
          write(nout,*) 'Too many PDF sets requested'
          write(nout,*) 'current limit is',maxpdf,' PDFs at once'
          pdfmemberto=pdfmemberfrom+maxpdf-1
          write(nout,*) 'pdfmemberto set to ',pdfmemberto
        Endif
        call InitPDF(pdfmembercentral)
        alphasin=alphasPDF(qcdrenscale)
        do nset=1,pdfmemberto-pdfmemberfrom+1
          call InitPDF(pdfmemberfrom-1+nset)
          alphasratio(nset)=alphasPDF(qcdrenscale)/alphasin
        enddo
        call InitPDF(pdfmembercentral)
#else
        alphasin=alphasPDF(qcdrenscale)
#endif
#else
       write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'LHAPDF not included: recompile with LHAPDF to use '//
     &            'LHAPDF interface, i.e. spdf=99!'
       write(*,*) 'To run HAWK without LHAPDF (spdf=1,11,14),'//
     &            ' change spdf in the input file'
       Stop
#endif
      endif

#ifdef RECOLA
      call init_recola_param
#endif

      end
***********************************************************************
      subroutine mcparameter
***********************************************************************
*     setting input parameters for Monte Carlo generator              *
*---------------------------------------------------------------------*
*     30.01.06  Ansgar Denner     last changed  30.01.06              *
***********************************************************************
      implicit none

#include "mcinput.h"
#include "mcparams.h"
#include "mccommon.h"


c rcoptions
      integer qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &                   qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp
      common/rcoptions/qborn,qw,qz,qschan,qtchan,qch2,qchint,
     &                   qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2,
     &                   qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp

c rscheme
      complex*16 dalpz,drbos,drferm
      integer qalp
      common/rscheme/dalpz,drbos,drferm,qalp

      integer iNumProcs, iErr, iMyRank
      common/mpistuff/iNumProcs,iMyRank

      integer selprocess
      common/selprocess/selprocess


c conversion of units
      mcgevtofb=gevtofbin      ! conversion factor GeV^-2 to fbarns

c particle masses in mc (in generator: mcmass)
      do mcpart=1,maxv
         mcm(mcpart)=0
         mcw(mcpart)=0
      enddo
      mcm(1)=mdin          ! d-quark mass
      mcm(2)=muin          ! u-quark mass
      mcm(3)=msin          ! s-quark mass
      mcm(4)=mcin          ! c-quark mass
      mcm(5)=mbin          ! b-quark mass
      mcm(6)=mtin          ! t-quark mass
      mcw(6)=0d0           ! t-quark width
      mcm(11)=mein         ! electron mass
      mcm(13)=mmuin        ! muon mass
      mcm(15)=mtauin       ! tau mass
      mcm(23)=mzin         ! z-boson mass
      mcw(23)=gzin         ! z-boson width
      mcm(24)=mwin         ! w-boson boson
      mcw(24)=gwin         ! w-boson width

      if (higgsflavor .eq. 'H') then
        mcm(25)=mhin          ! higgs mass
      else if (higgsflavor .eq. 'Hl') then
        mcm(25)=mhin         ! heavy higgs mass
      else if (higgsflavor .eq. 'Hh') then
        mcm(25)=mhhin         ! heavy higgs mass
      else if (higgsflavor .eq. 'Ha') then
        mcm(25)=mhain         ! peudo scalar higgs mass
      end if

      if(shtr.eq.0) then
        mcw(25)=0d0        ! higgs width
      else if(sgh.eq.0)then
        mcw(25)=ghin            ! higgs width
        if(iMyRank.EQ.0) then
          write(nout,99) mcw(25)
        endif
      elseif (sgh.eq.1)then
        call gridHt(mcm(25),mcw(25))
        if(iMyRank.EQ.0) then
          write(nout,99) mcw(25)
        endif
      elseif(sgh.eq.2) then
        if(iMyRank.EQ.0) then
          write(nout,*) "The Higgs width is taken from:"
        endif
        call HAWK_HTO(mcm(25),mcm(6),mcw(25),nout,iMyRank)
      else
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) 'input not supported: sgh = ',sgh
        stop
      endif
 99   format('   GH = ',g12.7)
      mcm(29)=0d0          ! singlet mass

c for output only
      mcalphasz=alphaszin

c parameters for subtraction terms
      mcalpha=alphain
      mcalphas=alphasin
      mclambda=lambdain
c parton masses for subtraction terms
c     (only first generation masses enter in matrix element of SD!)
      nparton =4
      m2parton(1) = mdin*mdin   !  d
      m2parton(2) = muin*muin   !  u
      m2parton(3) = mdin*mdin   !  s
      m2parton(4) = muin*muin   !  c
      m2electron = mein*mein


      do iparton=1,nparton
        m2parton(-iparton)=m2parton(iparton)
      enddo

c number of evaluation of virtual corrections
      nvirt=100


c other redundant parameters
      mcalpha0=alpha0in
      mcalphaz=alphazin
      mcgf    = gfin
      mccw= mwin/mzin


c     parameters for mappings:
c     mapping of fermion splitting function denominators with power
      powersplitf=0.8d0
c     mapping of photon splitting function denominators with power
      powersplitp=0.5d0
c     mapping of x in phase-space transformations with power
      powerexmap=0.8d0
c     default mapping of massless propagators with power
      powerprop=0.9d0
c technical parameters in h function and subroutine process
      if(ssub(1).gt.0.or.ssub(2).gt.0)then
        techparam(1)=1d-4       ! small negative mass2 in mapinv-function
      else
        techparam(1)=1d-8
      endif
      techparam(2)=1d-10*energy**2  ! allowed uncertainty tmax in process
      techparam(3)=1d-12*energy     ! allowed uncertainty (q1+q2)_T in process
      techparam(4)=1d-10            ! technical cut for x integration

c mass for effective collinear factors
      mcollreg2=1d-6

c cut on collinear invariant, should be smaller than mcollreg2/shat

      mceffcut=1d-12

#ifndef READINPUT
c     cut parameters
c     slicing cuts
      mcslideltac=1d-6
      mcslideltas=1d-4

#endif



      if(ssli(0).ne.0)then
        omxmaxf=mcslideltas
      elseif(ssub(0).ne.0.or.seff(0).ne.0)then
        omxmaxf=techparam(4)
      endif


      omxmaxp=techparam(4)

      if(ssli(0).ne.0d0.or.ssub(0).ne.0d0)then
      if (powersplitf.ge.1d0.and.omxmaxf.eq.0) then
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) ' powersplitf >= 1d0 and omxmaxf = 0d0 not allowed'
        stop
      elseif (powersplitp.ge.1d0.and.omxmaxp.eq.0) then
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) ' powersplitp >= 1d0 and omxmaxp = 0d0 not allowed'
        stop
      elseif (powerexmap.ge.1d0) then
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) ' powerexmap >= 1d0 not allowed'
        stop
      else if (powerprop.ge.1d0.and.techparam(1).eq.0d0) then
        write(*,'(/(a))') ' Error message: HAWK Stops here!'
        write(*,*) ' powerprop >= 1d0 and',
     &      ' techparam(1) = 0d0 not allowed'
        stop
      end if
      endif

c     parameters for subtraction branch
c     technical cuts for subtraction
      tcdeltae=0d0
      tcdeltas=1d-4
      tcdeltac=1d0


      tcdeltas=1d-10

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     non-standard couplings                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function nonstandardcoup(p1,p2,p3,g1,g2,g3,schannel)
      implicit none
c local variables
      character*3 p1,p2,p3,g1,g2,g3
      logical nonstandardcoup,schannel
c nonsmcoupl
      real*8 dkg,lg,ac,a0,a0t,gaah
      common/nonsmcoupl/dkg,lg,ac,a0,a0t,gaah
      nonstandardcoup=.true.
      if(gaah.ne.0d0.and.g1.eq.'h0 '.and.g2.eq.'ph '.and.g3.eq.'ph ')
     *  return
      if(ac.ne.0d0.or.a0.ne.0d0.or.a0t.ne.0d0)then
        if(g1.eq.'ph '.and.g2.eq.'ph '.and.g3.eq.'16 ')return
        if(g1.eq.'Z0 '.and.g2.eq.'Z0 '.and.g3.eq.'16~')return
      endif
      nonstandardcoup=.false.
      end


***********************************************************************
      subroutine settings(kbeam,k,lweight,event,generator,step)
***********************************************************************
*     filling of individual histograms for individual processes       *
*                                                                     *
*     selprocess = 0     Hjj     p p -> jet jet H                     *
*     selprocess = 1     Hln     p p -> l+ n H                        *
*     selprocess = 2     Hnl     p p -> n~ l- H                       *
*     selprocess = 3     Hll     p p -> l+ l- H                       *
*     selprocess = 4     Hnn     p p -> n~ n H                        *
*---------------------------------------------------------------------*
*     09.06.10  Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),lweight(maxcon)
      integer event,generator,step
      integer selprocess

      common/selprocess/selprocess

      if (selprocess.eq.0) then
         call settings_Hjj(kbeam,k,lweight,event,generator,step)
      elseif (selprocess.eq.1) then
         call settings_Hln(kbeam,k,lweight,event,generator,step)
      elseif (selprocess.eq.2) then
         call settings_Hnl(kbeam,k,lweight,event,generator,step)
      elseif (selprocess.eq.3) then
         call settings_Hll(kbeam,k,lweight,event,generator,step)
      elseif (selprocess.eq.4) then
         call settings_Hnn(kbeam,k,lweight,event,generator,step)
      endif

      end

***********************************************************************
      subroutine getqcdrenscale(scale)
***********************************************************************
*     qcd renormalization scale                                       *
*---------------------------------------------------------------------*
*     20.03.06  Ansgar Denner     last changed  24.07.07              *
***********************************************************************
      implicit none
      real*8   scale

#include "mcinput.h"
#include "mcparams.h"
#include "mccommon.h"

      qcdrenscale=scale*qcdrenscalefac
      mumsin = qcdrenscale

      qcdrenscale2=qcdrenscale*qcdrenscale

      end

***********************************************************************
      subroutine getqcdfacscale(scale)
***********************************************************************
*     qcd factorization scale                                         *
*---------------------------------------------------------------------*
*     20.03.06  Ansgar Denner     last changed  24.07.07              *
***********************************************************************
      implicit none
      real*8   scale

#include "mcparams.h"
#include "mccommon.h"

      qcdfacscale=scale*qcdfacscalefac

      qcdfacscale2=qcdfacscale*qcdfacscale

      end

***********************************************************************
      subroutine getqedfacscale(scale)
***********************************************************************
*     factorization scale                                             *
*---------------------------------------------------------------------*
*     21.09.06  Ansgar Denner     last changed  21.09.06              *
***********************************************************************
      implicit none
      real*8   scale

#include "mcparams.h"
#include "mccommon.h"

* local variables

      qedfacscale=scale*qedfacscalefac

      qedfacscale2=qedfacscale*qedfacscale

      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     photon or gluon recombination                               c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine recombination(kbeam,k,krec,generator,lnext,cut,
     *  switch)
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8  kbeam(2,0:3),k(maxe,0:3),krec(maxe,0:3)
      integer generator,lnext,cut,switch

      integer selprocess
      integer phasespace, em
      real*8 z

      common/selprocess/selprocess
      if (em.eq.0) then
         em = 3
      endif

      if (selprocess.eq.0)then
         call vbfh_recombination(kbeam,k,krec,generator,lnext,cut,
     *        switch)
      elseif (selprocess.ge.1) then
         call whzh_recombination(kbeam,k,krec,generator,lnext,cut,
     *        switch, z, phasespace, em)
      endif

      end



***********************************************************************
      subroutine applyhiggscut(kbeam,k,generator,lnext,cut,switch)
***********************************************************************
*     application of separation cuts on Higgs boson                   *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     08.02.10  Ansgar Denner     last changed  08.02.10              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),p(0:3),kvec(maxe),kbeamvec(2)
      real*8 kth,yh,kh(0:3),kt(maxe),y(maxe)
      real*8 s,prod,ktiktj,dR,dy,dphi,sh
      integer lnext,cut,generator,switch,prevcuts

      if(switch.gt.0)return

      prevcuts=40

      if(shtr.eq.0)then
        do vcomp=0,3
          kh(vcomp)=k(5,vcomp)
        enddo
      else
        do vcomp=0,3
          kh(vcomp)=k(5,vcomp)+k(6,vcomp)
        enddo
      endif

      if(scuts(generator).gt.0.or.shtr.gt.0)then

        cut=prevcuts

        ext=5                   ! higgs

c apply transverse momentum cuts
        kth=sqrt(kh(1)**2+kh(2)**2)
        cut=cut+1

        if(kth.lt.ptcuth)return

c apply energy cuts


        cut=cut+1


        if(kh(0).lt.ecuth)return


c apply rapidity cuts
        cut=cut+1
        if(kh(0)-kh(3).ne.0d0)then
          yh=.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))
        else
          return
        endif

        if(abs(yh).gt.ycuth)return

c apply invariant mass cuts

        sh=kh(0)*kh(0)-kh(1)*kh(1)-kh(2)*kh(2)
     &      -kh(3)*kh(3)

        cut=cut+1
        if(sh.gt.sucuth.or.sh.lt.slcuth) return

c cuts on Higgs decay singlets
        if(shtr.eq.2)then

c apply transverse momentum cuts
          do ext=5,6
            kt(ext)=sqrt(k(ext,1)**2+k(ext,2)**2)

            cut=cut+1
            if(kt(ext).lt.ptcuts)return
          enddo

c apply rapidity cuts
          do ext=5,6
            if(k(ext,0)-k(ext,3).ne.0d0)then
              y(ext)=.5d0*log((k(ext,0)+k(ext,3))/(k(ext,0)-k(ext,3)))
            else
              return
            endif
            cut=cut+1
            if(abs(y(ext)).gt.ycuts)return
          enddo


c require singlets to lie in rapidity between tagging jets
          if(lhiggsbjets) then
          do ext=3,4
            if(k(ext,0)-k(ext,3).eq.0) then
              y(ext) = 1d10
            elseif(k(ext,0)+k(ext,3).eq.0) then
              y(ext) = -1d10
            else
              y(ext)=.5d0*log((k(ext,0)+k(ext,3))
     &            /(k(ext,0)-k(ext,3)))
            endif
          enddo

          do ext=5,6
            cut=cut+1
            if(y(ext).lt.min(y(3),y(4)).or.y(ext).gt.max(y(3),y(4))
     &        .or.y(ext).eq.1d10.or.y(ext).eq.-1d10)
     &          return
          enddo
          endif

        endif
        if(cut.gt.maxcut)then
          write(*,*) ' Error message: HAWK Stops here!'
          write(*,*) 'You have implemented too many cuts!'
          write(*,*) 'reset maxcut = ',cut
          stop
        endif

      endif

      cut=0

      end


***********************************************************************
      subroutine init_cuts(generator)
***********************************************************************
*     initialization of cuts for individual processes                 *
*                                                                     *
*     selprocess = 0     Hjj     p p -> jet jet H                     *
*     selprocess = 1     Hln     p p -> l+ n H                        *
*     selprocess = 2     Hnl     p p -> n~ l- H                       *
*     selprocess = 3     Hll     p p -> l+ l- H                       *
*     selprocess = 4     Hnn     p p -> n~ n H                        *
*---------------------------------------------------------------------*
*     09.06.10  Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

      integer generator
      integer selprocess

      common/selprocess/selprocess

      if (selprocess.eq.0)then
         call init_cuts_Hjj(generator)
      elseif (selprocess.eq.1) then
         call init_cuts_Hln(generator)
      elseif (selprocess.eq.2) then
         call init_cuts_Hnl(generator)
      elseif (selprocess.eq.3) then
         call init_cuts_Hll(generator)
      elseif (selprocess.eq.4) then
         call init_cuts_Hnn(generator)
      endif

      end


**********************************************************************
#ifdef multpdf
        subroutine wpdf(x,mfact,lhadin,pdf,mpdf,pdf_phot,mpdf_phot,spdf)
#else
        subroutine wpdf(x,mfact,lhadin,pdf,pdf_phot,spdf)
#endif
**********************************************************************
*       spdf=1:   MRST2004QEDx with    O(alpha) effects
*       spdf=2:   MRST2004QEDx without O(alpha) effects
*       spdf=11:  CTEQ6M
*       spdf=14:  CTEQ6L1
*       spdf=99:  LHApdf
*       lhadin=1:  p
*       lhadin=-1:  pbar
*---------------------------------------------------------------------
*       6.1.06 Stefan Dittmaier     last changed 6.12.10 A. Denner
**********************************************************************
        implicit none
#include "mcparams.h"
#ifdef multpdf
#include "multpdf.h"
        integer nset
        double precision mpdf(maxpdf,-6:6),mpdf_phot(maxpdf)
#endif
c*** boundary for MRST2004QED
        real*8 pdf(-6:6),x,pdf_phot,lpdf(-6:6)
        real*8 xmin,xmax,qsqmin,qsqmax,mfact
        real*8 upv,dnv,usea,dsea,str,chm,bot,xglu,phot
        real*8 ctq6pdf
        integer i,spdf,lhadin
        real*8 xphoton
        logical has_photon

        data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/


        if (spdf.eq.1) then
c*** MRST2004QED [with O(alpha) effects]
c    ===================================
          if ((x.lt.xmin).or.(x.gt.xmax)) then
            do i=-6,6
              pdf(i) = 0d0
            enddo
            pdf_phot = 0d0
            return
          endif
          if (mfact**2.lt.qsqmin) mfact = sqrt(qsqmin)
          if (mfact**2.gt.qsqmax) mfact = sqrt(qsqmax)
          call mrstqed(x,mfact,1,
     &                  upv,dnv,usea,dsea,str,chm,bot,xglu,phot)
          pdf(0)   = xglu / x
          pdf(+1*lhadin)  = (dnv+dsea) / x
          pdf(-1*lhadin)  =     dsea / x
          pdf(+2*lhadin)  = (upv+usea) / x
          pdf(-2*lhadin)  =     usea / x
          pdf(+3*lhadin)  = str / x
          pdf(-3*lhadin)  = str / x
          pdf(+4*lhadin)  = chm / x
          pdf(-4*lhadin)  = chm / x
          pdf(+5*lhadin)  = bot / x
          pdf(-5*lhadin)  = bot / x
          pdf_phot = phot / x
        elseif (spdf.eq.2) then

c*** MRST2004QEDx [without O(alpha) effects]
c    =======================================
          if ((x.lt.xmin).or.(x.gt.xmax)) then
            do i=-6,6
              pdf(i) = 0d0
            enddo
            pdf_phot = 0d0
            return
          endif
          if (mfact**2.lt.qsqmin) mfact = sqrt(qsqmin)
          if (mfact**2.gt.qsqmax) mfact = sqrt(qsqmax)
          call mrstqed(x,mfact,3,
     &                  upv,dnv,usea,dsea,str,chm,bot,xglu,phot)
          pdf(0)   = xglu / x
          pdf(+1*lhadin)  = (dnv+dsea) / x
          pdf(-1*lhadin)  =     dsea / x
          pdf(+2*lhadin)  = (upv+usea) / x
          pdf(-2*lhadin)  =     usea / x
          pdf(+3*lhadin)  = str / x
          pdf(-3*lhadin)  = str / x
          pdf(+4*lhadin)  = chm / x
          pdf(-4*lhadin)  = chm / x
          pdf(+5*lhadin)  = bot / x
          pdf(-5*lhadin)  = bot / x
          pdf_phot = 0d0


        elseif (spdf.ge.10.and.spdf.lt.20) then

          pdf(0)  = Ctq6Pdf( 0,x,mfact)
          pdf(+1*lhadin) = Ctq6Pdf( 2,x,mfact)
          pdf(-1*lhadin) = Ctq6Pdf(-2,x,mfact)
          pdf(+2*lhadin) = Ctq6Pdf( 1,x,mfact)
          pdf(-2*lhadin) = Ctq6Pdf(-1,x,mfact)
          pdf(+3*lhadin) = Ctq6Pdf( 3,x,mfact)
          pdf(-3*lhadin) = Ctq6Pdf(-3,x,mfact)
          pdf(+4*lhadin) = Ctq6Pdf( 4,x,mfact)
          pdf(-4*lhadin) = Ctq6Pdf(-4,x,mfact)
          pdf(+5*lhadin) = Ctq6Pdf( 5,x,mfact)
          pdf(-5*lhadin) = Ctq6Pdf(-5,x,mfact)




        elseif (spdf.eq.99) then
#ifdef includeLHAPDF
#ifdef multpdf
          do nset=1,pdfmemberto-pdfmemberfrom+1
            call InitPDF(pdfmemberfrom-1+nset)
            if (has_photon().eqv..TRUE.) then
               call evolvePDFphoton(x,mfact,lpdf,xphoton)
               mpdf_phot(nset) = xphoton/x
             else
               call evolvepdf(x,mfact,lpdf)
            endif
            mpdf(nset,0)  = lpdf(0)/x
            mpdf(nset,+1*lhadin) = lpdf(+1)/x
            mpdf(nset,-1*lhadin) = lpdf(-1)/x
            mpdf(nset,+2*lhadin) = lpdf(+2)/x
            mpdf(nset,-2*lhadin) = lpdf(-2)/x
            mpdf(nset,+3*lhadin) = lpdf(+3)/x
            mpdf(nset,-3*lhadin) = lpdf(-3)/x
            mpdf(nset,+4*lhadin) = lpdf(+4)/x
            mpdf(nset,-4*lhadin) = lpdf(-4)/x
            mpdf(nset,+5*lhadin) = lpdf(+5)/x
            mpdf(nset,-5*lhadin) = lpdf(-5)/x
          enddo
          If(pdfmemberto-pdfmemberfrom.GE.0) Then
            call InitPDF(pdfmembercentral)
          Endif
          if (has_photon().eqv..TRUE.) then
             call evolvePDFphoton(x,mfact,lpdf,xphoton)
             pdf_phot = xphoton/x
           else
             call evolvepdf(x,mfact,lpdf)
          endif
          pdf(0)  = lpdf(0)/x
          pdf(+1*lhadin) = lpdf(+1)/x
          pdf(-1*lhadin) = lpdf(-1)/x
          pdf(+2*lhadin) = lpdf(+2)/x
          pdf(-2*lhadin) = lpdf(-2)/x
          pdf(+3*lhadin) = lpdf(+3)/x
          pdf(-3*lhadin) = lpdf(-3)/x
          pdf(+4*lhadin) = lpdf(+4)/x
          pdf(-4*lhadin) = lpdf(-4)/x
          pdf(+5*lhadin) = lpdf(+5)/x
          pdf(-5*lhadin) = lpdf(-5)/x
#else
          if (has_photon().eqv..TRUE.) then
             call evolvePDFphoton(x,mfact,lpdf,xphoton)
             pdf_phot = xphoton/x
           else
             call evolvepdf(x,mfact,lpdf)
          endif

          pdf(0)  = lpdf(0)/x
          pdf(+1*lhadin) = lpdf(+1)/x
          pdf(-1*lhadin) = lpdf(-1)/x
          pdf(+2*lhadin) = lpdf(+2)/x
          pdf(-2*lhadin) = lpdf(-2)/x
          pdf(+3*lhadin) = lpdf(+3)/x
          pdf(-3*lhadin) = lpdf(-3)/x
          pdf(+4*lhadin) = lpdf(+4)/x
          pdf(-4*lhadin) = lpdf(-4)/x
          pdf(+5*lhadin) = lpdf(+5)/x
          pdf(-5*lhadin) = lpdf(-5)/x

#endif
#else
       write(*,'(/(a))') ' Error message: HAWK Stops here!'
       write(*,*) 'LHAPDF not included: recompile with LHAPDF to use '//
     &            'LHAPDF interface, i.e. spdf=99!'
       write(*,*) 'To run HAWK without LHAPDF (spdf=1,11,14),'//
     &            ' change spdf in the input file'
       Stop
#endif
        else
          write(*,'(/(a))') ' Error message: HAWK Stops here!'
          write(*,*) 'wpdf: not defined: spdf= ',spdf
          stop
        endif

        end
