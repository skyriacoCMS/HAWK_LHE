#define PRIVATE
***********************************************************************
*     file montecarlo.F                                               *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     30.01.06  Ansgar Denner     last changed  05.11.12              *
***********************************************************************


#define PRINTWEIGHTMAX

#define RANLUX

* switch subprocesses on
c#define SUBPROCESSES

c#define ANYTEST


***********************************************************************
* NOTES                                                               *
*                                                                     *
*   LOOPS OVER EMITTER AND SPECTATOR ASSUME UNCHARGED HIGGS!!!!       *
*   need to steer charge flow from composite particles!!!!!           *
*                                                                     *
*       fac(phasespace)=0.5d0/pi/gsum(phasespace)*wconv               *
*                                                                     *
*       factor(gen)=mcgevtofb/(2d0*pi)**(3*next(gen)-10)              *
*    &    /(2d0*sqrt(energy**4                                       *
*    &      -mcm(iabs(hepnum(1,gen)))**2*mcm(iabs(hepnum(2,gen)))**2))*
*                                                                     *
*     weight(1,phasespace)=m2pdfpdf_0(phasespace)*wconv               *
*                     *factor(generator)/gsum(phasespace)             *
*                                                                     *
*                                                                     *
*   exit 1:      generation: s < sppmin                               *
*   exit 2:      density:    sinv < 0                                 *
*   exit 3:      inv:        smin > smax                              *
*   exit 4:      inv:        smin < 0                                 *
*   exit 5:      inv:        1/density <= 0 (width>0)                 *
*   exit 6:      inv:        1/density <= 0 (width=0)                 *
*   exit 7:      process:    dabs(q1(1:2)+q2(1:2)) > acc              *
*   exit 8:      process:    lambda < 0                               *
*   exit 9:      process:    d < 0                                    *
*   exit 10:     process:    lambdas <= 0                             *
*   exit 11:     process:    lambdat <= 0                             *
*   exit 12:     process:    tmin >= tmax                             *
*   exit 13:     process:    tmax > 0                                 *
*   exit 14:     process:    1/density = 0 (width>0)                  *
*   exit 15:     process:    1/density = 0 (width=0)                  *
*   exit 16:     process:    s <= 0                                   *
*   exit 17:     decay:      lambda <= 0                              *
*   exit 18:     decay:      s <= 0                                   *
*   exit 19:     decay:      1/density <= 0                           *
*   exit 20:     mapinv:     xmin-m2, xmax-m2 < 0                     *
*   exit 21:     jacobianinv:xmin-m2, xmax-m2 < 0                     *
*   exit 22:     rotation:   |cost| > 1                               *
*   exit 23:     boost:      m2 <= 0                                  *
*   exit 24:     mapin:      1/g <= 0                                 *
*   exit 25:     mapout:     1/g <= 0                                 *
*   exit 26:     transverse: pvec = 0                                 *
*   exit 27:     integrand:  x = 1                                    *
*   exit 28:     integrand:  gsum = 0                                 *
*   exit 29:     integrand:  s(em,sp,...) < 0                         *
*   exit 30:     integrand:  s(em,sp,...) < 0                         *
*   exit 31:     integrand:  s(em,sp,1) <= 0                          *
*   exit 32:     conv:       x=1                                      *
*   exit 33:     dli2gen:    converges badly                          *
*   exit 34:     genbrems  : m2 <= 0                                  *
*   exit 35:     momentum conservation violated:                      *
*   exit 36:     onshellness violated:                                *
*   exit 40:     mapx:       xmin-m2, xmax-m2 < 0                     *
*   exit 41:     jacobianx:  xmin-m2, xmax-m2 < 0                     *
*   exit 42:     mapout:     z = 1 (eliminated)                       *
*   exit 43:     mapout:     z = 0 (eliminated)                       *
*   exit 44:     dets:       sinv = 0 (eliminated)                    *
*   exit 45:     ktcluster:  |cos(deltaphi)| > 0                      *
*   exit 46:     recombination:  |cos(deltaphi)| > 0                  *
*   exit 47:     integrand:   integrand: s(1,2,1) < 0                 *
*   exit 48:     integrand:   integrand: s(em,sp,1) < 0               *
*                                                                     *
***********************************************************************
* CHANGES with respect to MAR                                         *
*                                                                     *
*   nsteps                                                            *
*                                                                     *
***********************************************************************
* TO BE ADAPTED                                                       *
*                                                                     *
*    SWITCH weak                                                      *
*                                                                     *
*    histogram to weight(contr)                                       *
*                                                                     *
***********************************************************************
*                                                                     *
*     Monte Carlo generator                                           *
*                                                                     *
*     written by Markus Roth          2005                            *
*     adapted by Ansgar Denner                                        *
*                                                                     *
*     last changed  05.11.2012                                        *
*                                                                     *
***********************************************************************
      program montecarlo
      implicit none

#ifdef mpiuse
      include "mpif.h"
#endif

#include "mcparams.h"
#include "mccommon.h"

c local variables
      integer unweighting
      real*8  s,m2,ps3,ps4

      integer iNumProcs, iErr, iMyRank
      common/mpistuff/iNumProcs,iMyRank

#ifdef mpiuse
      call MPI_Init(iErr)
      call MPI_Comm_size(MPI_COMM_WORLD, iNumProcs, iErr)
      call MPI_Comm_rank(MPI_COMM_WORLD, iMyRank, iErr)
#else
      iMyRank=0
      iNumProcs=1
#endif

      if(iMyRank.EQ.0) then
	Print *,""
        print *,
     &      '======================================================='
        print *,' This is HAWK, version 3.0.0                       '
        print *,'  released May 3, 2019                             '
        print *,' authors: A. Denner, S. Dittmaier,                 '
        print *,'          S. Kallweit, J.-N. Lang, A. Muck         '
        print *,
     &      '======================================================='
#ifdef mpiuse
	Print *,"(running on ",iNumProcs," cores)"
#endif
	Print *,""
      endif


c initialization
      call initialization
c numerical integration
      unweighting=0
      call integration()
c write result
      call writeresult

#ifdef mpiuse
       call MPI_Finalize(iErr)
#endif

      stop

      s = energy**2
      m2 = mcm(25)**2

      if (m2.ne.0d0) then
        ps3 = 1d0/(2*pi)**5*(pi/2)**2/2/s
     &      *(s**2-m2**2+2*m2*s*log(m2/s))
        ps4 = 1d0/(2*pi)**8*(pi/2)**3/12/s
     &      *((s-m2)*(s**2+10*m2*s+m2**2)+6*m2*s*(s+m2)*log(m2/s))
      else
        ps3 = 1d0/(2*pi)**5*(pi/2)**2/2*s
        ps4 = 1d0/(2*pi)**8*(pi/2)**3/12*s**2
      endif

      write(*,*) 's   = ',s
      write(*,*) 'mh2 = ',m2
      write(*,*) 'ps3 = ',ps3*mcgevtofb/2/s
      write(*,*) 'ps4 = ',ps4*mcgevtofb/2/s

      end


***********************************************************************
      subroutine initialization
***********************************************************************
*     initialization of phase space integration                       *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     08.03.06  Ansgar Denner     last changed  22.02.10              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

#ifdef RANLUX
      integer ranluxseed
      common/ranluxseed/ranluxseed
#endif

#ifdef multpdf
#include "multpdf.h"
      integer nset
#endif
#ifdef mpiuse
      include "mpif.h"
#endif
      integer iNumProcs, iErr, iMyRank
      common/mpistuff/iNumProcs,iMyRank

c local variables
      real*8 sppmin,mtotcut
      logical endname
      integer i1
      integer switchinitrans
      character*3 madname(-maxv:maxv)

c date and time
      real*8     rcpu
      integer*4  int,icpu,icpuh,icpum,icpus
      character  datum*10,da1*2,da2*2,da3*4
      character  chtime*8
      integer iarray(3)

      real*8  qu,qd,ql,qn,qf,mu,mu2,md,md2,mlep,mlep2
      complex*16 guu,gdd,gnn,gll
      common/qf/qu,qd,ql,qn,qf(4),mu,mu2,md,md2,mlep,mlep2,
     &            guu(-1:1),gdd(-1:1),gnn(-1:1),gll(-1:1)

      real*8 constff,constfi,constif,constii
      common/constsub/constff,constfi,constif,constii

      integer oc
      common/ranluxoutput/oc

      save switchinitrans

c initialize
      av=0d0
      si=0d0
      average=0d0
      sigma=0d0
      sigmaa=0d0
      avborn=0d0
      siborn=0d0
#ifdef multpdf
          do nset=1,pdfmemberto-pdfmemberfrom+1
	    maverage(nset)=0d0
	    msigma(nset)=0d0
	    msigmaa(nset)=0d0
            mavborn(nset)=0d0
            msiborn(nset)=0d0
	  enddo
#endif
      wmax=0d0
      weighttotmax=0d0
      nunw=0
      nneg=0
      nevt=1
      nmaxtot=0
      do icut=1,maxcut
        mccutevt_0(icut)=0d0
      enddo
      do gen=0,maxg
        do icon=1,maxcon
          do icut=1,maxcut
            do iphsp=1,maxps
              mccutsubevt(icut,gen,icon,iphsp)=0d0
            enddo
            mccutevt(icut,gen,icon)=0d0
          enddo
        enddo
        do iexit=1,maxexit
          mcexitevt(iexit,gen)=0d0
        enddo
        do icon=1,maxcon
          processname(gen)=' '
          avgen(gen,icon)=0d0
          sigen(gen,icon)=0d0
#ifdef multpdf
          do nset=1,pdfmemberto-pdfmemberfrom+1
            mavgen(nset,gen,icon)=0d0
            msigen(nset,gen,icon)=0d0	
	  enddo
#endif
          nevtgen(gen,icon)=0d0
          weightmax(gen,icon)=0d0
          rejgen(gen,icon)=0
          nmax(gen,icon)=0
          channelmax(gen,icon)=0
          do cha=1,maxch
            avgencha(gen,icon,cha)=0d0
            sigencha(gen,icon,cha)=0d0
            weightmaxcha(gen,icon,cha)=0d0
            nevtgencha(gen,icon,cha)=0
            rejcha(gen,icon,cha)=0
            nmaxcha(gen,icon,cha)=0
          enddo
        enddo
      enddo

#ifdef SUBPROCESSES
      do parton1=-4,4
      do parton2=-4,4
        averageii(parton1,parton2)=0d0
        sigmaii(parton1,parton2)=0d0
        avbornii(parton1,parton2)=0d0
        sibornii(parton1,parton2)=0d0
        do gen=0,maxg
          do icon=1,maxcon
            avgenii(gen,icon,parton1,parton2)=0d0
            sigenii(gen,icon,parton1,parton2)=0d0
          enddo
        enddo
      enddo
      enddo
#endif

      do ipart=-maxip,maxip
        hepname(ipart)=' '
      enddo
      do ipart=-maxip,maxip
        class(ipart)=' '
      enddo
      do mcpart=1,maxv
        mcm(mcpart)=0d0
        mcw(mcpart)=0d0
      enddo
      do gen=0,maxg
      do ext=1,maxe
        hepnum(ext,gen)=0
        name(ext,gen)=' '
        ecut(ext,gen)=0d0
        ptcut(ext,gen)=0d0
        ycut(ext,gen)=1d10
        ecutgen(ext,gen)=0d0
        do ext2=1,maxe
          dycut(ext,ext2,gen)=0d0
          dRcut(ext,ext2,gen)=0d0
          scut(ext,ext2,gen)=0d0
          ccut(ext,ext2,gen)=0d0
          scutgen(ext,ext2,gen)=0d0
          ccutgen(ext,ext2,gen)=0d0
        enddo
      enddo
      enddo
      mtotcut=0d0

      xpmax=0
      xpmin=1
      x1max=0
      x1min=1
      x2max=0
      x2min=1
      x3max=0
      x3min=1

      ptcutj1=0d0        ! cut on pt of leading jet
      ptcutj2=0d0        ! cut on pt of subleading jet
      ycutj1=1d10        ! cut on rapidity of leading jet
      ycutj2=1d10        ! cut on rapidity of subleading jet
      ptmaxj1=1d30       ! upper cut on pt of leading jet
      ptmaxj2=1d30       ! upper cut on pt of subleading jet
      ptmaxj3=1d30       ! upper cut on pt of third jet
      yminj1=1d10        ! lower cut on rapidity of leading jet
      yminj2=1d10        ! lower cut on rapidity of subleading jet
      yminj3=1d10        ! lower cut on rapidity of third jet
      ptcutvis=0d0       ! lower cut on (p_jet1+p_jet2+p_h)_T
      ptmaxvis=1d30      ! upper cut on (p_jet1+p_jet2+p_h)_T
      dphicutj1h=0d0     ! cut on azimuth between leading jet and Higgs
      dphicutj2h=0d0     ! cut on azimuth between subleading jet and Higgs
      mucutjj=1d30       ! upper cut on invariant mass between two leading jets
      mucutjjh=1d30      ! upper cut on invariant mass between two leading jets and Higgs
      mlcutjj=0d0        ! lower cut on invariant mass between two leading jets
      dycutjj=0d0        ! rapidity cut between two leading jets

      ptcuth=0d0         ! Higgs pt cut
      ecuth=0d0          ! Higgs energy cut
      ycuth=1d10         ! Higgs rapidity cut
      mlcuth=1d-2        ! lower cut on Higgs invariant mass
      mucuth=1d30        ! upper cut on Higgs invariant mass

      ptcuts=0d0         ! singlet pt cut
      ycuts=1d10         ! singlet rapidity cut


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     setting parameters                                          c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c fixing SM/BSM input
      call inputparameter

c reading input
      call readinput
#ifdef RECOLA
         if (sendp(0) .ne. 1) then
           write(*,*) 'Error message: HAWK Stops here!'
           write(*,*) 'Only the input sendp=1 is compatible '//
     &                'when using RECOLA.'
           write(*,*) 'Please set sendp=1 in the input card.'
           stop
         end if
#endif

c     open outputfile


c     open lhe file                                                                                                                                                                                                             
      open(unit=1988,file="output.lhe",status='REPLACE')
      write(1988,'(A)') '<LesHouchesEvents version="1.0">'
      write(1988,'(A)') '<!--    '
      write(1988,'(A)') ' HAWK 3.0.0    '      
      write(1988,'(A)') ' HAWK -> lhe by  S.Kyriacou    '
      write(1988,'(A)') '--->    '
      write(1988,'(A)') '<header>'
      write(1988,'(A)') '</header>'
      write(1988,'(A)') '<init>'
      write(1988,'(A)') '</init>'
      
      


   
      endname=.false.
      do i1=1,50
        if(outputfile(i1:i1).eq.' ') endname=.true.
        if(endname)outputfile(i1:i1)=' '
      enddo
      if(outputfile(1:1).eq.' ')then
        nout=6
      else
        nout=22
	open(unit=nout,file=outputfile,status='unknown')
        if(iMyRank.EQ.0) then
          write(*,'(/A30,A60)') " Output is written into file: ",
     &    outputfile
	  write(nout,*) ""
	  write(nout,*)
     &        '======================================================='
	  write(nout,*)
     &        ' This is HAWK, version 3.0.0                       '
	  write(nout,*)
     &        '  released May 3, 2019                             '
	  write(nout,*)
     &        ' authors: A. Denner, S. Dittmaier,                 '
	  write(nout,*)
     &        '          S. Kallweit, J.-N. Lang, A. Mueck        '
	  write(nout,*)
     &        '======================================================='
#ifdef mpiuse
	  write(nout,*) "(running on ",iNumProcs," cores)"
#endif
	  write(nout,*) ""
        endif
      endif
      oc=nout


c fixing dependent parameters

      call fixparameter



 97   format(i2.2,'.',i2.2,'.20',i2.2)
 98   format(i2.2,'.',i2.2,'.',i2.2)
      call date_and_time(datum)
      da1=datum(7:8)
      da2=datum(5:6)
      da3=datum(1:4)
      datum=da1//'.'//da2//'.'//da3
      call itime(iarray)
      write(chtime,98)iarray(1),iarray(2),iarray(3)
      call cpu_time(rcpu)

      icpuh=int(rcpu/3600.)
      icpum=int(rcpu/60.-icpuh*60)
      icpus=int(rcpu-icpuh*3600-icpum*60)
      icpu =icpuh*60+icpum
      If(imyRank.EQ.0) Then
        write(nout,100) datum,chtime,icpuh,icpum,icpus
100     format(/' date: ',a10,' ,    time: ',a8,' hour,',
     &       '    cpu: ',i3,' h',i3,' m',i3,' s')
      Endif

      if(iMyRank.EQ.0) then
        write(nout,*)
        write(nout,*) 'Process: ',initialstate,'-> ',finalstate
      endif

c setting parameters for amplitudes
      call initqqqqhamps
      If(imyRank.EQ.0) Then
        call paramout
      Endif

c initialize names of possible external particles
      call init_particles()

c setting parmeters for monte carlo generator
      call mcparameter


      If(imyRank.EQ.0) Then
        write(nout,'(a)')' '
        write(nout,'(a)')'                       Monte Carlo'
        write(nout,'(a)')'                       ==========='
      Endif

      do gen=0,ngenerator
c setting external particle parameters
        processname(gen)=name(1,gen)//' '//name(2,gen)//' -> '
        do ext=3,nextgen(gen)
          processname(gen)=processname(gen)(1:54)//name(ext,gen)
          do i1=70,2,-1
            if(processname(gen)(i1-1:i1).eq.'  ')
     &          processname(gen)=processname(0)(1:i1-2)//
     &          processname(gen)(i1:70)
          enddo
        enddo

        do ext=1,nextgen(gen)
          do ipart=-maxip,maxip
            if(name(ext,gen).eq.hepname(ipart)) hepnum(ext,gen)=ipart
          enddo
          if(hepnum(ext,gen).eq.0)then
            write(*,'(a24,a16)')' Unknown particle name: ',name(ext,gen)
            write(*,'(a16,2x,i5)')' Model: smodel= ',smodel
            stop
          endif
      enddo
      enddo

      If(imyRank.EQ.0) Then
        write(nout,200)  processname(0)
 200    format(' Input process:  ',a60)

        if(lnoutmc.ge.3)then
          write(nout,300) hepnum(1,0),hepnum(2,0),hepnum(3,0),
     &        (hepnum(ext,0),ext=4,nextgen(0))
 300      format(' Input process:  ',i4,' + ',i4,' -> ',
     &        i3,(9(:,' + ',i3,:)))
        endif
      Endif

c input-parameter scheme
      if(salp.eq.0.or.salp.eq.3)then
	mcalpha=mcalpha0
      elseif(salp.eq.1)then
	mcalpha=mcalphaz
      elseif(salp.eq.2)then
	mcalpha=sqrt(2d0)*mcgf/pi*mcm(24)**2*(1d0-mccw**2)
      else
        write(*,'(a)')
     &      ' Unknown input for input-parameter scheme (salp)!'
        stop
      endif

c inititalize beam momenta
      pbeam(1,0)=0.5d0*energy
      pbeam(1,1)=0d0
      pbeam(1,2)=0d0
      pbeam(1,3)=0.5d0*energy
      pbeam(2,0)=0.5d0*energy
      pbeam(2,1)=0d0
      pbeam(2,2)=0d0
      pbeam(2,3)=(-0.5d0)*energy

c definition of end-point part of subtraction functions (const=Cij-1/2)
      do ext=1,maxe
      do ext2=1,maxe
        if(ssoft.eq.1)then
c racoonww paper eq(4.29)
          if(ext.ge.3.and.ext2.ge.3)const(ext,ext2)=1d0-pi2/3d0
          if(ext.ge.3.and.ext2.le.2)const(ext,ext2)=1d0-pi2/2d0
          if(ext.le.2.and.ext2.ge.3)const(ext,ext2)=-1.5d0+pi2/6d0
          if(ext.le.2.and.ext2.le.2)const(ext,ext2)=1.5d0-pi2/3d0
        elseif(ssoft.eq.2)then
          const(ext,ext2)=-0.5d0
c Cij = 0
        elseif(ssoft.eq.3)then
c YFS inspired
          if(ext.ge.3.and.ext2.ge.3)const(ext,ext2)=0.5d0-pi2*2d0/3d0
          if(ext.ge.3.and.ext2.le.2)const(ext,ext2)=0.5d0-pi2/6d0
          if(ext.le.2.and.ext2.ge.3)const(ext,ext2)=0.5d0-pi2/6d0
          if(ext.le.2.and.ext2.le.2)const(ext,ext2)=0.5d0-pi2*2d0/3d0
        endif
      enddo
      enddo

c output
      numout=0

c initialize phase-space generator
      mcmass(0)=0d0
      mcwidth(0)=0d0
      do mcpart=1,maxv
        mcmass(mcpart)=mcm(mcpart)
        mcwidth(mcpart)=mcw(mcpart)
      enddo


c neglect light-fermion masses
      if(slightfermions.eq.1)then
        do mcpart=1,4
          mcmass(mcpart)=0d0
        enddo
        do mcpart=11,13,2
          mcmass(mcpart)=0d0
        enddo
      elseif(slightfermions.eq.2)then
        do mcpart=1,5
          mcmass(mcpart)=0d0
        enddo
        do mcpart=11,15,2
          mcmass(mcpart)=0d0
        enddo
      endif

c initialize histograms


      call settings(pbeam,pevent,weight(1,1),0,0,1)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     initialize subcontributions                                 c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do gen=0,ngenerator

      if(sgen(gen).ne.0)then
c process name
        processname(gen)=name(1,gen)//' '//name(2,gen)//' -> '
        do ext=3,nextgen(gen)
          processname(gen)=processname(gen)(1:54)//name(ext,gen)
          do i1=70,2,-1
            if(processname(gen)(i1-1:i1).eq.'  ')
     &          processname(gen)=processname(gen)(1:i1-2)//
     &          processname(gen)(i1:70)
          enddo
        enddo


c initialize flux and phasespacefactor
        factor(gen)=mcgevtofb/(2d0*pi)**(3*nextgen(gen)-10)
     &    /(2d0*sqrt(energy**4
     &      -mcm(iabs(hepnum(1,gen)))**2*mcm(iabs(hepnum(2,gen)))**2))




c definitions for subtraction functions
        lambda2(gen)=mclambda*mclambda

c properties of external particles
        do ext=1,nextgen(gen)

c masses of external particles (for subtraction functions)
c does not work for composite particles!!!!
          mass2(ext,gen)=mcm(abs(hepnum(ext,gen)))**2


c direction of external particles

c charges of external fermions
          i1=abs(hepnum(ext,gen))


          mccharge(ext,gen)=0d0
          if(i1.eq.101.or.i1.eq.102) then
            mccharge(ext,gen)=1d20   ! only switch!
          else if(i1.eq.1.or.i1.eq.3.or.i1.eq.5) then
            mccharge(ext,gen)=-(1d0/3d0)
          else if(i1.eq.2.or.i1.eq.4.or.i1.eq.6) then
            mccharge(ext,gen)=2d0/3d0
          else if(i1.eq.11.or.i1.eq.13.or.i1.eq.15) then
            mccharge(ext,gen)=-1d0
          endif

        enddo

c photon recombination
        if(srecomb(gen).ne.0)then
          do ext=1,nextsub(gen)
            recparticle(ext,gen)=0
            if(mccharge(ext,gen).ne.0d0)recparticle(ext,gen)=1
          enddo
          if(srecomb(gen).eq.1.or.srecomb(gen).eq.2)then
            continue
          else
            write(*,*)' initialization: wrong input for srecomb:',
     &          srecomb(gen),gen
            stop
          endif
        endif

        call init_cuts(gen)
        sucuth=mucuth*mucuth
        slcuth=mlcuth*mlcuth

c initializing phase space generators
        if(lnoutgen.gt.0)then
          If(imyRank.EQ.0) Then
            write(nout,'(a)') ' '
            write(nout,'(a)') ' Initializing phase-space generators'
          Endif
        endif
        nchannel(gen)=0
        nphsp(gen)=maxps       ! to ensure correct initialization of
                               ! arrays in integrand for first event
        ncontr(gen)=maxcon     ! to ensure correct initialization of
                               ! arrays in integrand for first event

c determine partonic processes
        if(shtr.eq.0)then
          call init_subprocesses(energy,sppmin,gen,
     &        nextgen(gen),smodel,sincludecuts(gen),
     &        ssub(gen)+ssli(gen)+seff(gen))
        else
          call init_subprocesses(energy,sppmin,gen,
     &        nextgen(gen),smodel,sincludecuts(gen),
     &        ssub(gen)+ssli(gen)+seff(gen))
        endif


       if(sincludecuts(gen).ge.2) sppmin=max(sppmin,mtotcutgen**2)

       if(lnoutgen.gt.0)then
         If(imyRank.EQ.0) Then
           write(nout,*) 'initgenerator sppmin=',sppmin,sqrt(sppmin)
         endif
       endif

        xmin(gen)=sppmin/energy**2
        if(xmin(gen).gt.1d0)then
          write(*,'(a)')' error: energy too low for particle production'
          stop
        endif

c initialize adaptation
        noptran(gen)=2000*nchannel(gen)
        alminran=0.2d0
        do igrv=1,maxr
          beran(0,igrv,gen)=0d0
          do igrp=1,maxgr
            wran(igrp,igrv,gen)=0d0
            beran(igrp,igrv,gen)=dble(igrp)/dble(maxgr)
            alran(igrp,igrv,gen)=1d0/dble(maxgr)
          enddo
        enddo
      endif
      enddo

c initialize random number generator
      if(switchinitrans.eq.0)then
          If(imyRank.EQ.0) Then
            write(nout,'(a)')
            write(nout,'(a)') " Initializing Random Number Generator"
          Endif
#ifdef RANLUX
        call rluxgo(3,ranluxseed,0,0)
c       call rluxgo(0,ranluxseed,0,0)
#else
        call initcarry
#endif
        switchinitrans=1
      endif

c write output
      If(imyRank.EQ.0) Then
        call writeoutput
      Endif

      end

***********************************************************************
      subroutine writeoutput
***********************************************************************
*     writes output                                                   *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     21.02.06  Ansgar Denner     last changed  11.03.13              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"
#ifdef multpdf
#include "multpdf.h"
#endif

      integer selprocess
      common/selprocess/selprocess

c local variables
      integer i1
      character*70 contentname

      logical lwritecuts,lwriterec
      data lwritecuts /.true./, lwriterec /.true./


 1008 format(/" Number of unweighted events = ",i10,"   (",es9.3,")")
 1018 format(/" Number of weighted events   = ",i10,"   (",es9.3,")")
 1009 format( " Collider CMS energy         = ",f12.5, " TeV")
      write(nout,1018) nevents,dfloat(nevents)

      if (lnoutmc.le.1) then
        write(nout,1009) energy/1d3
      endif

c output
      if(lnoutmc.ge.2)then
        write(nout,'(a)')' '
        write(nout,'(a)')' Input parameters:',
     *   ' (masses and widths in GeV)'
        write(nout,'(a)')' '
 1000   format(9x,'energy = ',f18.12,',',/,
     *    7x,'alpha(0) = 1/', f16.12,',',2x,'alpha(MZ) = 1/',f16.12,','
     *      ,/,
     *    12x,'G_F = ',g18.8,',',6x,'alpha = 1/',f16.12,',',/,
     *    4x,'alpha_s(MZ) = ',f18.12,',',4x,'alpha_s = ',f18.12,',',/
     *    4x,'        MH  = ',f18.12,',',4x,'     GH = ',f18.12,'.',/)
 2000 format(12x,'mgl = ',f18.12,',',8x,'mmu = ',f18.12,',',/,
     *    12x,'ma0 = ',f18.12,',',9x,'tb = ',f18.12,',',/,
     *    13x,'m1 = ',f18.12,',',9x,'m2 = ',f18.12,',',/,
     *    10x,'msusy = ',f18.12,',',7x,'atau = ',f18.12,',',/,
     *    13x,'at = ',f18.12,',',9x,'ab = ',f18.12,'.',/)
 3000   format(9x,'msl(1) = ',f9.4,',  msl(2) = ',f9.4,
     *    ',  msl(3) = ',f9.4,',',/,
     *    9x,'mse(1) = ',f9.4,',  mse(2) = ',f9.4,
     *    ',  mse(3) = ',f9.4,',',/,
     *    9x,'msq(1) = ',f9.4,',  msq(2) = ',f9.4,
     *    ',  msq(3) = ',f9.4,',',/,
     *    9x,'msu(1) = ',f9.4,',  msu(2) = ',f9.4,
     *    ',  msu(3) = ',f9.4,',',/,
     *    9x,'msd(1) = ',f9.4,',  msd(2) = ',f9.4,
     *    ',  msd(3) = ',f9.4,'.',/)
        write(nout,1000)energy,1d0/mcalpha0,1d0/mcalphaz,mcgf
     *    ,1d0/mcalpha,mcalphasz,mcalphas,mcm(25),mcw(25)
        if(smodel.eq.3.or.smodel.eq.4.or.smodel.eq.34)then
          write(nout,2000)mgl,mue,ma0,tb,mgo1,mgo2,msusy,atau,at,ab
          write(nout,3000)msl(1),msl(2),msl(3),mse(1),mse(2),mse(3),
     &      msq(1),msq(2),msq(3),msu(1),msu(2),msu(3),msd(1),msd(2),
     &      msd(3)
        endif
      endif

      if(lnoutmc.ge.3)then
c output particles
        write(nout,'(a)')' '
        write(nout,'(a)')
     &      ' Particle           Mass                 Width'
        do ipart=1,ncomp
          if(hepname(ipart).ne.' ')then
            if(mcm(ipart).ge.1d0.and.mcw(ipart).ne.0d0)then
              write(nout,'(1x,a15,f15.10," GeV  ",f15.10," GeV")')
     &            hepname(ipart),mcm(ipart),mcw(ipart)
            elseif(mcm(ipart).gt.0d0.and.mcw(ipart).ne.0d0)then
              write(nout,'(1x,a15,f15.10," MeV  ",f15.10," MeV")')
     &            hepname(ipart),1d3*mcm(ipart),1d3*mcw(ipart)
            elseif(mcm(ipart).ge.1d0.and.mcw(ipart).eq.0d0)then
              write(nout,'(1x,a15,f15.10," GeV")')
     &            hepname(ipart),mcm(ipart)
            elseif(mcm(ipart).gt.0d0.and.mcw(ipart).eq.0d0)then
              write(nout,'(1x,a15,f15.10," MeV")')
     &            hepname(ipart),1d3*mcm(ipart)
            elseif(mcm(ipart).eq.0d0.and.mcw(ipart).eq.0d0)then
              write(nout,'(1x,a15,"    massless")')hepname(ipart)
            else
              write(*,'(a)')' writeoutput: width for massless particle!'
              stop
            endif
          endif
        enddo
        write(nout,'(a)')' '
        write(nout,'(a)')' Composite particle contains partons        '
        do ipart=ncomp,maxip
          if(hepname(ipart).ne.' ')then
            contentname = ' '
            do mcpart=-naux+1,naux-1
              if (partcont(ipart,mcpart)) then
                contentname=contentname(1:54)//hepname(mcpart)
                do i1=70,2,-1
                  if(contentname(i1-1:i1).eq.'  ')
     &                contentname=contentname(1:i1-2)//
     &                contentname(i1:70)
                enddo
              endif
            enddo
            write(nout,'(1x,a15,a50)')
     &          hepname(ipart),contentname(1:50)
          endif
        enddo
      endif

c output factorization scales
      write(nout,*)
      write(nout,'(a30,1g20.10)')'QCD factorization scale   =',
     &    qcdfacscale
      write(nout,'(a30,1g20.10)')'QED factorization scale   =',
     &    qedfacscale

c output renormalization scales
      write(nout,'(a30,1g20.10)')'QCD renormalization scale =',
     &    qcdrenscale











      if(lnoutmc.ge.5)then
        write(nout,'(a)')' '
        write(nout,'(a)')' Input parameters in full precision:'
        write(nout,'(a)')' '
 1101   format(9x,'energy = ',d25.16,',',/,
     *       7x,'alpha(0) = 1/', d25.16,',',
     &       2x,'alpha(MZ) = 1/',d25.16,',',/,
     *       12x,'G_F = ',g18.8,',',6x,'alpha = 1/',d25.16,',',/,
     *       4x,'alpha_s(MZ) = ',d25.16,',',
     &       4x,'alpha_s = ',d25.16,'.',/)
 2101   format(12x,'mgl = ',d25.16,',',8x,'mmu = ',d25.16,',',/,
     *       12x,'ma0 = ',d25.16,',',9x,'tb = ',d25.16,',',/,
     *       13x,'m1 = ',d25.16,',',9x,'m2 = ',d25.16,',',/,
     *       10x,'msusy = ',d25.16,',',7x,'atau = ',d25.16,',',/,
     *       13x,'at = ',d25.16,',',9x,'ab = ',d25.16,'.',/)
 3101   format(9x,'msl(1) = ',d25.16,',  msl(2) = ',d25.16,
     *       ',  msl(3) = ',d25.16,',',/,
     *       9x,'mse(1) = ',d25.16,',  mse(2) = ',d25.16,
     *       ',  mse(3) = ',d25.16,',',/,
     *       9x,'msq(1) = ',d25.16,',  msq(2) = ',d25.16,
     *       ',  msq(3) = ',d25.16,',',/,
     *       9x,'msu(1) = ',d25.16,',  msu(2) = ',d25.16,
     *       ',  msu(3) = ',d25.16,',',/,
     *       9x,'msd(1) = ',d25.16,',  msd(2) = ',d25.16,
     *       ',  msd(3) = ',d25.16,'.',/)
        write(nout,1101)energy,1d0/mcalpha0,1d0/mcalphaz,mcgf,
     *       1d0/mcalpha,mcalphasz,mcalphas

      if(smodel.eq.3.or.smodel.eq.4.or.smodel.eq.34)then
        write(nout,2101)mgl,mue,ma0,tb,mgo1,mgo2,msusy,atau,at,ab
        write(nout,3101)msl(1),msl(2),msl(3),mse(1),mse(2),mse(3),
     &    msq(1),msq(2),msq(3),msu(1),msu(2),msu(3),msd(1),msd(2),msd(3)
      endif

      if(lnoutmc.ge.3)then
c output particles
        write(nout,'(a)')' '
        write(nout,'(a)')
     &      ' Particle           Mass                 Width'
        do ipart=1,ncomp
          if(hepname(ipart).ne.' ')then
            if(mcm(ipart).ge.1d0.and.mcw(ipart).ne.0d0)then
              write(nout,'(1x,a15,d25.16," GeV  ",d25.16," GeV")')
     &            hepname(ipart),mcm(ipart),mcw(ipart)
            elseif(mcm(ipart).gt.0d0.and.mcw(ipart).ne.0d0)then
              write(nout,'(1x,a15,d25.16," MeV  ",d25.16," MeV")')
     &            hepname(ipart),1d3*mcm(ipart),1d3*mcw(ipart)
            elseif(mcm(ipart).ge.1d0.and.mcw(ipart).eq.0d0)then
              write(nout,'(1x,a15,d25.16," GeV")')
     &            hepname(ipart),mcm(ipart)
            elseif(mcm(ipart).gt.0d0.and.mcw(ipart).eq.0d0)then
              write(nout,'(1x,a15,d25.16," MeV")')
     &            hepname(ipart),1d3*mcm(ipart)
            elseif(mcm(ipart).eq.0d0.and.mcw(ipart).eq.0d0)then
              write(nout,'(1x,a15,"    massless")')hepname(ipart)
            else
              write(*,'(a)')' writeoutput: width for massless particle!'
              stop
            endif
          endif
        enddo
        write(nout,'(a)')' '
        write(nout,'(a)')' Composite particle contains partons        '
        do ipart=ncomp,maxip
          if(hepname(ipart).ne.' ')then
            contentname = ' '
            do mcpart=-naux+1,naux-1
              if (partcont(ipart,mcpart)) then
                contentname=contentname(1:54)//hepname(mcpart)
                do i1=70,2,-1
                  if(contentname(i1-1:i1).eq.'  ')
     &                contentname=contentname(1:i1-2)//
     &                contentname(i1:70)
                enddo
              endif
            enddo
            write(nout,'(1x,a15,a50)')
     &          hepname(ipart),contentname(1:50)
          endif
        enddo
      endif

c output factorization scales
      write(nout,*)
      write(nout,'(a28,d25.16)')'QCD factorization scale   =',
     &    qcdfacscale
      write(nout,'(a28,d25.16)')'QED factorization scale   =',
     &    qedfacscale

c output renormalization scales
      write(nout,'(a28,d25.16)')'QCD renormalization scale =',
     &    qcdrenscale


      endif

























      if(lnoutmc.ge.5) then
c output mapping parameters
      write(nout,*)	
      write(nout,'(a33,1g20.10)')'propagator mapping with power = ',
     &    powerprop
      write(nout,'(a33,1g20.10)')'propagator mapping with cut   = ',
     &    techparam(1)
      write(nout,'(a33,1g20.10)')'allowed uncertainty for tmax  = ',
     &    techparam(2)
      write(nout,'(a33,1g20.10)')'allowed uncertainty for P_T,in= ',
     &    techparam(3)
      write(nout,'(a33,1g20.10)')'mapping of x with cut         = ',
     &    techparam(4)
      write(nout,'(a33,1g20.10)')'1/(1-x) mapping with power = ',
     &    powersplitf
      write(nout,'(a33,1g20.10)')'ln(1-x) mapping with power = ',
     &    powersplitp
        write(nout,'(a33,1g20.10)')'phase-space mapping with power= ',
     &      powerexmap
      endif

c output technical parameters
      write(nout,*)
      if(ssli(0).ne.0)then
        write(nout,'(a33,1e20.10)')'slicing cut:           deltas = ',
     &      mcslideltas
        write(nout,'(a33,1e20.10)')'slicing cut:           deltac = ',
     &      mcslideltac
      endif
      if(seff(0).ne.0)then
        write(nout,'(a33,1e20.10)')'coll. mass param.:  mcollreg2 = ',
     &      mcollreg2
        write(nout,'(a33,1e20.10)')'technical cut:  deltaspk/shat = ',
     &      mceffcut
      endif
      if(lnoutmc.ge.5) then
      if(ssub(0).ne.0)then
        write(nout,'(a33,1e20.10)')'subtraction cut:       deltas = ',
     &      tcdeltas
        write(nout,'(a33,1e20.10)')'subtraction cut:       deltae = ',
     &      tcdeltae
        write(nout,'(a33,1e20.10)')'subtraction cut:       deltac = ',
     &      tcdeltac
      endif
      write(nout,'(a33,1e20.10)')'splittingfunction cut:omxmaxf = ',
     &    omxmaxf
      write(nout,'(a33,1e20.10)')'splittingfunction cut:omxmaxp = ',
     &    omxmaxp
      endif

      if(lnoutmc.ge.3) then
        write(nout,'(a29,i2)')'mapping of pdfs: pdfmap = ',
     &      pdfmap
      endif

c output options
c     write(nout,'(a)')' '

c initial state
      if(sppbar.eq.1.and.shadin1.eq.1.and.shadin2.eq.-1)then
        write(nout,'(a16,i2,(a))')'sppbar=',sppbar,
     &      ': p pbar initial state'
      elseif(sppbar.eq.0.and.shadin1.eq.1.and.shadin2.eq.1)then
        write(nout,'(a16,i2,(a))')'sppbar=',sppbar,
     &      ': p p initial state'
      elseif(sppbar.eq.3.and.shadin1.eq.-1.and.shadin2.eq.1)then
        write(nout,'(a16,i2,(a))')'sppbar=',sppbar,
     &      ': pbar p initial state'
      elseif(sppbar.eq.2.and.shadin1.eq.-1.and.shadin2.eq.-1)then
        write(nout,'(a16,i2,(a))')'sppbar=',sppbar,
     &      ': pbar pbar initial state'
      elseif(sppbar.ne.0.and.sppbar.ne.1.and.sppbar.ne.2
     &      .and.sppbar.ne.3) then
        write(nout,'(a16,i2,(a))')'sppbar=',sppbar,
     &      ': undefined initial state'
        stop
      else
        write(nout,'(a15,i2,(a))')'sppbar=',sppbar,
     &      ': inconsistent switches'
        write(nout,*) 'shadin1,shadin2=',shadin1,shadin2
        stop
      endif
c
      if(spdf.eq.1)then
        write(nout,'("           spdf=",i2,(a))') spdf,
     &      ': MRST MRST2004QED [with O(alpha) effects]'
      elseif(spdf.eq.2)then
        write(nout,'("           spdf=",i2,(a))') spdf,
     &      ': MRST MRST2004QED [without O(alpha) effects]'
      elseif(spdf.eq.11)then
        write(nout,'("           spdf=",i2,(a))') spdf,
     &      ': CTEQ6M [Standard MSbar scheme (alphas(mz)=0.118)]'
      elseif(spdf.eq.14)then
        write(nout,'("           spdf=",i2,(a))') spdf,
     &      ': CTEQ6L1 [Leading Order (alphas(mz)=0.130)]'
#ifdef includeLHAPDF
      elseif(spdf.eq.99)then
        if(pdfpath.ne.'')then
          write(nout,'("           spdf=",i2,(a),/20x,"pdfset ",i3)')
     &        spdf,': LHApdf: '//pdfpath,pdfmember
        else
          write(nout,'("           spdf=",i2,(a),/20x,"pdfset ",i3)')
     &        spdf,': LHApdf: '//pdfname,pdfmember
        endif
#endif
      elseif(sconv(gen).ne.0)then
        write(*,'(a23,i4)')' Unknown input: spdf = ',spdf
        stop
      endif

#ifdef multpdf
      If(pdfmemberto.GE.pdfmemberfrom) Then
        write(nout,*) ''
        write(nout,*) 'Usage of multiple PDF sets requested:'
	write(nout,'(" for sets ",i4," to ",i4,
     &              " results will be provided")')
     &                pdfmemberfrom,pdfmemberto
        write(nout,*) ''
      Endif
#endif
      if(lnoutmc.ge.3)then
c model
        if(smodel.eq.12)then
          write(nout,'(a15,i3,(a))')'smodel=',smodel,
     &        ': Standard Model with QCD'
        elseif(smodel.eq.1)then
          write(nout,'(a15,i3,(a))')'smodel=',smodel,
     &        ': Standard Model without QCD'
        elseif(smodel.eq.2)then
          write(nout,'(a15,i3,(a))')'smodel=',smodel,
     &        ': QCD'
        elseif(smodel.eq.15)then
          write(nout,'(a16,i3,(a))')'smodel=',smodel,
     &        ': Standard Model with Higgs singlet without: QCD'
        elseif(smodel.eq.125)then
          write(nout,'(a16,i3,(a))')'smodel=',smodel,
     &        ': Standard Model with Higgs singlet and QCD'
        elseif(smodel.eq.34)then
          write(nout,'(a15,i3,(a))')'smodel=',smodel,
     &        ': MSSM with SQCD'
        elseif(smodel.eq.3)then
          write(nout,'(a15,i3,(a))')'smodel=',smodel,
     &        ': MSSM without SQCD'
        elseif(smodel.eq.4)then
          write(nout,'(a15,i3,(a))')'smodel=',smodel,
     &        ': SQCD'
        else
          write(*,'(a30,i3)')' Unknown input: smodel = ',smodel
          stop
        endif

c input-parameter scheme
        if(salp.eq.0)then
          write(nout,'(a16,i2,(a))')'salp=',salp,
     &        ': alpha(0)-parametrization scheme'
        elseif(salp.eq.1)then
          write(nout,'(a16,i2,(a))')'salp=',salp,
     &        ': alpha(MZ)-parametrization scheme'
        elseif(salp.eq.2)then
          write(nout,'(a16,i2,(a))')'salp=',salp,
     &        ': GF-parametrization scheme'
        elseif(salp.eq.3)then
          write(nout,'(a16,i2,(a))')'salp=',salp,
     &        ': mixed alpha(0)-alpha_GF-parametrization scheme'
        else
          write(*,'(a)')' Unkown input: salp = ',salp
          stop
        endif

        if(scp.eq.0)then
          write(nout,'(a16,i2,(a))')'scp=',scp,
     &        ': CP symmetry not used'
        elseif(scp.eq.1)then
          write(nout,'(a16,i2,(a))')'scp=',scp,
     &        ': CP symmetry used in lowest order matrix elements'
        else
          write(*,'(a30,i3)')' Unknown input: scp = ',scp
          stop
        endif

        if(swdiag.eq.0)then
          write(nout,'(a16,i2,(a))')'swdiag=',swdiag,
     &        ': W-exchange diagrams excluded '
        elseif(swdiag.eq.1)then
          write(nout,'(a16,i2,(a))')'swdiag=',swdiag,
     &        ': W-exchange diagrams included '
        else
          write(*,'(a30,i3)')' Unknown input: swdiag = ',swdiag
          stop
        endif

        if(szdiag.eq.0)then
          write(nout,'(a16,i2,(a))')'szdiag=',szdiag,
     &        ': Z-exchange diagrams excluded '
        elseif(szdiag.eq.1)then
          write(nout,'(a16,i2,(a))')'szdiag=',szdiag,
     &        ': Z-exchange diagrams included '
        else
          write(*,'(a30,i3)')' Unknown input: szdiag = ',szdiag
          stop
        endif

        if(swidth.eq.0)then
          write(nout,'(a16,i2,(a))')'swidth=',swidth,
     &        ': zero width'
        elseif(swidth.eq.1)then
          write(nout,'(a16,i2,(a))')'swidth=',swidth,
     &        ': complex mass scheme'
        elseif(swidth.eq.2)then
          write(nout,'(a16,i2,(a))')'swidth=',swidth,
     &        ': naive fixed width scheme'
        else
          write(*,'(a30,i3)')' Unknown input: swidth = ',swidth
          stop
        endif

c ir-treatment
        if(sirtr.eq.0)then
           write(nout,'(a16,i2,(a))')'sirtr=',sirtr,
     &        ': subtraction method used'
        elseif(sirtr.eq.1)then
           write(nout,'(a16,i2,(a))')'sirtr=',sirtr,
     &          ': phase-space slicing  used'
        else
           write(*,'(a)')' Unknown input: sirtr = ',sirtr
           stop
        endif
      endif

c lepton recombination
      if(selprocess.ge.1.and.selprocess.le.3)then
         if(sbarelep.eq.1)then
            write(nout,'(a16,i2,(a))')'sbarelep=',sbarelep,
     &           ': no lepton recombination with photons:'//
     &           ' bare leptons'
         elseif(sbarelep.eq.0)then
            write(nout,'(a16,i2,(a))')'sbarelep=',sbarelep,
     &           ': leptons are recombined with photons'
         endif

         if(colllep.ge.1.and.colllep.le.3)then
            write(nout,'(a16,i2,(a))')'colllep=',colllep,
     &           ': initial gamma -> f fbar splitting included'
         elseif(colllep.eq.0) then
           if(lnoutmc.ge.3)then
              write(nout,'(a16,i2,(a))')'colllep=',colllep,
     &           ': no initial gamma -> f fbar splitting included'
           endif
         else
            write(*,'(a30,i3)')' Unknown input: colllep = ',colllep
            stop
         endif
      endif

c weak virtual corrections
      if(sew.eq.1)then
        write(nout,'(a16,i2,(a))')'sew=',sew,
     &      ': electroweak corrections included'
      elseif(sew.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sew = ',sew
        stop
      endif

      if(lnoutmc.ge.3)then
c weak virtual corrections
        if(sweak.eq.1)then
          write(nout,'(a16,i2,(a))')'sweak=',sweak,
     &        ': electroweak virtual corrections included'
        elseif(sweak.eq.2)then
          write(nout,'(a16,i2,(a))')'sweak=',sweak,
     &        ': only fermionic virtual corrections included'
        elseif(sweak.eq.3)then
          write(nout,'(a16,i2,(a))')'sweak=',sweak,
     &        ': only bosonic virtual corrections included'
        elseif(sweak.ne.0)then
          write(*,'(a30,i3)')' Unknown input: sweak = ',sweak
          stop
        endif
      endif

c incoming photon corrections
      if(spinc(0).eq.1)then
        write(nout,'(a16,i2,(a))')'spinc=',spinc(0),
     &      ': incoming photon corrections included'
      elseif(spinc(0).ne.0)then
        write(*,'(a30,i3)')' Unknown input: spinc = ',spinc(0)
        stop
      endif

c two-loop heavy Higgs corrections
      if(selprocess.eq.0)then
      if(shh2.eq.1)then
        write(nout,'(a16,i2,(a))')'shh2=',shh2,
     &    ': 2-loop heavy Higgs corrections included'
      elseif(shh2.ne.0)then
        write(*,'(a30,i3)')' Unknown input: shh2 = ',shh2
        stop
      endif
      endif

c qcd corrections
      if(sqcd.eq.1)then
        write(nout,'(a16,i2,(a))')'sqcd=',sqcd,
     &    ': QCD corrections included'
      elseif(sqcd.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sqcd = ',sqcd
        stop
      endif
      if(selprocess.eq.0)then
c diagonal qcd corrections
      if(sqcddiag.eq.1)then
        write(nout,'(a16,i2,(a))')'sqcddiag=',sqcddiag,
     &    ': diagonal QCD corrections included'
      elseif(sqcddiag.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sqcddiag = ',sqcddiag
        stop
      endif
c nondiagonal qcd corrections
      if(sqcdnondiag.eq.1)then
        write(nout,'(a16,i2,(a))')'sqcdnondiag=',sqcdnondiag,
     &    ': nondiagonal QCD corrections included'
      elseif(sqcdnondiag.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sqcdnondiag = ',sqcdnondiag
        stop
      endif
c qqg qcd corrections
      if(sqcdgsplit.eq.1)then
        write(nout,'(a16,i2,(a))')'sqcdgsplit=',sqcdgsplit,
     &    ': qqg QCD corrections included'
      elseif(sqcdgsplit.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sqcdgsplit = ',sqcdgsplit
        stop
      endif
c gg fusion qcd corrections
      if(sqcdggfus.eq.1)then
        write(nout,'(a16,i2,(a))')'sqcdggfus=',sqcdggfus,
     &    ': gg fusion QCD corrections included'
      elseif(sqcdggfus.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sqcdggfus = ',sqcdggfus
        stop
      endif
      endif

      if(selprocess.eq.0)then
c t-channel diagrams
      if(stcha.eq.0)then
        write(nout,'(a16,i2,(a))')'stcha=',stcha,
     &    ': t-channel diagrams not included'
      elseif(stcha.ne.1)then
        write(*,'(a30,i3)')' Unknown input: stcha = ',stcha
        stop
      endif
c s-channel diagrams
      if(sscha.eq.0)then
        write(nout,'(a16,i2,(a))')'sscha=',sscha,
     &    ': s-channel diagrams not included'
      elseif(sscha.ne.1)then
        write(*,'(a30,i3)')' Unknown input: sscha = ',sscha
        stop
      endif
c squared matrix elements
      if(sch2.eq.0)then
        write(nout,'(a16,i2,(a))')'sch2=',sch2,
     &    ': squared matrix elements not included'
      elseif(sch2.ne.1)then
        write(*,'(a30,i3)')' Unknown input: sch2 = ',sch2
        stop
      endif
c interferences
      if(schint.eq.0)then
        write(nout,'(a16,i2,(a))')'schint=',schint,
     &    ': interferences not included'
      elseif(schint.ne.1)then
        write(*,'(a30,i3)')' Unknown input: schint = ',schint
        stop
      endif
      endif

c bottom quark contributions in initial state
      if(sbini.eq.1)then
        write(nout,'(a16,i2,(a))')'sbini=',sbini,
     &    ': initial bottoms included in LO'
      elseif(sbini.eq.0)then
        write(nout,'(a16,i2,(a))')'sbini=',sbini,
     &    ': initial bottoms not included in LO'
      elseif(sbini.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sbini = ',sbini
        stop
      endif
c bottom quark contributions in final state
      if(sbfin.eq.1)then
        write(nout,'(a16,i2,(a))')'sbfin=',sbfin,
     &    ': final bottoms included in LO'
      elseif(sbfin.eq.0)then
        write(nout,'(a16,i2,(a))')'sbfin=',sbfin,
     &    ': final bottoms not included in LO'
      elseif(sbfin.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sbfin = ',sbfin
        stop
      endif

      if(lnoutmc.ge.3)then
c end-point part of subtraction terms
        if(ssoft.eq.1.)then
          write(nout,'(a16,i2,(a))')'ssoft=',ssoft,
     &        ': end-point part of subtraction functions of STD'
        elseif(ssoft.ne.0)then
          write(*,'(a30,i3)')' Unknown input: ssoft = ',ssoft
          stop
        endif

c light fermion masses neglected
        if(slightfermions.eq.1)then
          write(nout,'(a16,i2,(a))')' slightfermions=',slightfermions,
     &        ': small fermion masses neglected (el,mu,uq,dq,cq,sq)'
        elseif(slightfermions.eq.2)then
          write(nout,'(a16,i2,(a))')' slightfermions=',slightfermions,
     &        ': small fermion masses neglected'//
     &        ' (el,mu,ta,uq,dq,cq,sq,bq)'
        elseif(slightfermions.ne.0)then
          write(*,'(a30,i3)')' Unknown input: slightfermions = ',
     &        slightfermions
          stop
        endif

c real part of bremsstrahlung corrections
        if(sreal.eq.1)then
          write(*,'(a)')' sreal = 1 not implemented'
          stop
          write(nout,'(a16,i2,(a))')'sreal=',sreal,
     &        ': imaginary part of virtual corrections neglected'
        elseif(sreal.ne.0)then
          write(*,'(a30,i3)')' Unknown input: sreal = ',sreal
          stop
        endif
      endif

c qed factorization
      if(sfactqed.eq.1)then
        write(nout,'("       sfactqed=",i2,(a))') sfactqed,
     &      ': MSbar factorization for QED'
      else if(sfactqed.eq.2)then
        write(nout,'("       sfactqed=",i2,(a))') sfactqed,
     &      ': DIS factorization for QED  '
      elseif(sfactqed.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sfactqed = ',sfactqed
        stop
      endif
c qcd factorization
      if(sfactqcd.eq.1)then
        write(nout,'("       sfactqcd=",i2,(a))') sfactqcd,
     &      ': MSbar factorization for QCD'
      else if(sfactqcd.eq.2)then
        write(nout,'("       sfactqcd=",i2,(a))') sfactqcd,
     &      ': DIS factorization for QCD  '
      elseif(sfactqcd.ne.0)then
        write(*,'(a30,i3)')' Unknown input: sfactqcd = ',sfactqcd
        stop
      endif
      if(sfactqcd.eq.0.and.sfactqcd.ne.0.or.
     &    sfactqcd.ne.0.and.sfactqcd.eq.0) then
        write(*,'(a30,i3)')' Unknown input: sfactqed = ',sfactqed
        write(*,'(a30,i3)')' Unknown input: sfactqcd = ',sfactqcd
        stop
      endif

      if(lnoutmc.ge.3)then
c frame for momentum generation
        if(lcmsgen)then
          write(nout,'("        lcmsgen=",l2,(a))') lcmsgen,
     &        ': momenta generated in CMS'
        else
          write(nout,'("        lcmsgen=",l2,(a))') lcmsgen,
     &        ': momenta generated in LS'
        endif
      endif

      if(shtr.eq.0)then
        write(nout,'("           shtr=",i2,(a))') shtr,
     &      ': on-shell Higgs boson'
      elseif(shtr.eq.1)then
        write(nout,'("           shtr=",i2,(a))') shtr,
     &      ': off-shell Higgs boson'
      elseif(shtr.eq.2)then
        write(nout,'("           shtr=",i2,(a))') shtr,
     &      ': Higgs decaying to two singlets'
      else
        write(*,'(a30,i3)')' Unknown input: shtr = ',shtr
        stop
      endif

      if(shtr.gt.0)then
        if(sgh.eq.0)then
          write(nout,'("            sgh=",i2,(a))') sgh,
     &        ': input Higgs-boson width'
        elseif(sgh.eq.1)then
          write(nout,'("            sgh=",i2,(a))') sgh,
     &        ': on shell Higgs-boson width'
        elseif(sgh.eq.2)then
          write(nout,'("            sgh=",i2,(a))') sgh,
     &        ': complex-mass scheme Higgs-boson width'//
     &         ' according to cpHTO11.f by G. Passarino'
        else
          write(*,'(a30,i3)')' Unknown input: sgh = ',sgh
          stop
        endif
         if(shbw.eq.0)then
            write(nout,'("           shbw=",i2,(a))') shbw,
     &           ': standard Breit Wigner'
         elseif(shbw.eq.1)then
            write(nout,'("           shbw=",i2,(a))') shbw,
     &           ': off-shell propagator according to arXiV:1107.0683'
         elseif(shbw.eq.2)then
            write(nout,'("           shbw=",i2,(a))') shbw,
     &           ': Breit Wigner with s-dependent width'
         elseif(shbw.eq.8)then
            write(nout,'("           shbw=",i2,(a))') shbw,
     &           ': Breit Wigner with s-dependent width in numerator'
         elseif(shbw.eq.9)then
            write(nout,'("           shbw=",i2,(a))') shbw,
     &           ': off-shell propagator with fixed width'
         else
            write(*,'(a30,i3)')' Unknown input: shbw = ',shbw
            stop
         endif
         if(sonsproj.eq.0)then
            write(nout,'("       sonsproj=",i2,(a))') sonsproj,
     &           ': no on-shell projection'
         elseif(sonsproj.eq.1)then
            write(nout,'("       sonsproj=",i2,(a))') sonsproj,
     &           ': on-shell projection via rescaling '//
     &           'of final-state momenta'
         elseif(sonsproj.eq.2)then
            write(nout,'("       sonsproj=",i2,(a))') sonsproj,
     &           ': on-shell projection via rescaling '//
     &           'of all momenta'
         else
            write(*,'(a30,i3)')' Unknown input: sonsproj = ',sonsproj
            stop
         endif

         if(sborn(0).gt.0.and.lonsborn)then
            write(nout,'("       lonsborn=",l2,(a))') lonsborn,
     &           ': LO ME projected on-shell'
         endif

         if(sqcd.gt.0.and.lonsqcd)then
            write(nout,'("        lonsqcd=",l2,(a))') lonsqcd,
     &           ': NLO QCD ME projected on-shell'
         endif

         if(sew.gt.0.and.lonsew)then
            write(nout,'("         lonsew=",l2,(a))') lonsew,
     &           ': NLO EW ME projected on-shell'
         endif
      endif

c subcontributions

      do gen=0,ngenerator
      if(sgen(gen).ne.0)then
       if(lnoutmc.ge.2)then
        write(nout,*)' '
c subcontributions
        if(gen.eq.0)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for Born or virtual corrections'
        elseif(gen.eq.1)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for real photonic corrections'
        elseif(gen.eq.2)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for real gluonic corrections'
        elseif(gen.eq.3)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial photon contributions 1'
        elseif(gen.eq.4)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial photon contributions 2'
        elseif(gen.eq.5)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial gluon contributions 1'
        elseif(gen.eq.6)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial gluon contributions 2'
        elseif(gen.eq.7)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial lepton contributions 1'
        elseif(gen.eq.8)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial lepton contributions 2'
        elseif(gen.eq.9)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial l+ contributions 1'
        elseif(gen.eq.10)then
          write(nout,'("        sgen(",i1,")=",i2,(a))')gen,sgen(gen),
     &      ': branch for initial l+ contributions 2'
        else
          write(*,'(a)')' Unknown input: sgen = ',sgen(gen)
          stop
        endif
c event numbers
        if(freqgencon(0,6).ne.freqgencon(0,7))then
          write(*,'(a30,2i5)')' Unknown input: freqgencon 6=/=7 ',
     &        freqgencon(0,6),freqgencon(0,7)
          stop
        endif
        if(freqgencon(0,8).ne.freqgencon(0,9))then
          write(*,'(a30,2i5)')' Unknown input: freqgencon 8=/=9 ',
     &        freqgencon(0,8),freqgencon(0,9)
          stop
        endif
        if(freqgencon(0,2).lt.freqgencon(0,5))then
          write(*,*)' ew frequency freqgencon(0,2) > ',
     &        ' strong frequency freqgencon(0,5)',
     &        ' strong frequency reset to electroweak frequency',
     &        freqgencon(0,2),freqgencon(0,5)
          freqgencon(0,5)=freqgencon(0,2)
        endif

        if (sbarelep.eq.0) then
           ncontr(0)=9
        else
           ncontr(0)=10
        endif
        if(gen.eq.0)then
          do icon=1,ncontr(gen)
            if(sgencon(gen,icon).eq.1)then
              write(nout,'("   nevents(",i1,",",i2,")=",i10)')
     &            gen,icon,nevents/freqgencon(gen,icon)
            endif
          enddo
        else
          if(sgen(gen).eq.1)then
            write(nout,'("     nevents(",i1,")=",i10)')gen,
     &          nevents/freqgencon(gen,1)
          endif
        endif
c born amplitude
        if(sborn(gen).eq.1)then
          write(nout,'("       sborn(",i1,")=",i2,(a))')gen,sborn(gen),
     &      ': tree-level amplitude included'
        elseif(sborn(gen).ne.0)then
          write(*,'(a25,i3)')' Unknown input: sborn = ',sborn(gen)
          stop
        endif
c additional information for tree-level cross section
        if(swborn(gen).eq.1)then
          write(nout,'("      swborn(",i1,")=",i2,(a))')gen,swborn(gen),
     &      ': additional result for tree-level cross section'
        elseif(swborn(gen).ne.0)then
          write(*,'(a30,i3)')' Unknown input: swborn = ',swborn
          stop
        endif

        if(lnoutmc.ge.3)then
c convolution
          if(sconv(gen).eq.1)then
            write(nout,'("       sconv(",i1,")=",i2,(a))')
     &          gen,sconv(gen),
     &          ': convolution over parton distributions included'
            if(spdf.eq.1)then
              write(nout,'("           spdf=",i2,(a))') spdf,
     &            ': MRST MRST2004QED [with O(alpha) effects]'
            elseif(spdf.eq.2)then
              write(nout,'("           spdf=",i2,(a))') spdf,
     &            ': MRST MRST2004QED [without O(alpha) effects]'
            elseif(spdf.eq.11)then
              write(nout,'("           spdf=",i2,(a))') spdf,
     &            ': CTEQ6M [Standard MSbar scheme (alphas(mz)=0.118)]'
            elseif(spdf.eq.14)then
              write(nout,'("           spdf=",i2,(a))') spdf,
     &            ': CTEQ6L1 [Leading Order (alphas(mz)=0.130)]'
#ifdef includeLHAPDF
            elseif(spdf.eq.99)then
              if(pdfpath.ne.'')then
                write(nout,'("           spdf=",i2,(a),
     &              /20x,"pdfset ",i3)')
     &              spdf,': LHApdf: '//pdfpath,pdfmember
              else
                write(nout,'("           spdf=",i2,(a),
     &              /20x,"pdfset ",i3)')
     &              spdf,': LHApdf: '//pdfname,pdfmember
              endif
#endif
            elseif(sconv(gen).ne.0)then
              write(*,'(a)')' Unknown input: spdf = ',spdf
              stop
            endif
          elseif(sconv(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: sconv = ',sconv(gen)
            stop
          endif

c virtual corrections
          if(svirt(gen).eq.1)then
            write(nout,'("       svirt(",i1,")=",i2,(a))')
     &          gen,svirt(gen),
     &          ': virtual corrections included'
          elseif(svirt(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: svirt = ',svirt(gen)
            stop
          endif
c photonic bremsstrahlung corrections
          if(spbrems(gen).eq.1)then
            write(nout,'("     spbrems(",i1,")=",i2,(a))')
     &          gen,spbrems(gen),
     &          ': real photonic corrections included'
          elseif(spbrems(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: spbrems = ',spbrems(gen)
            stop
          endif
c gluonic bremsstrahlung corrections
          if(sgbrems(gen).eq.1)then
            write(nout,'("     sgbrems(",i1,")=",i2,(a))')
     &          gen,sgbrems(gen),
     &          ': real gluonic corrections included'
          elseif(sgbrems(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: sgbrems = ',sgbrems(gen)
            stop
          endif
c incoming photon corrections
          if(spinc(gen).eq.1)then
            write(nout,'("       spinc(",i1,")=",i2,(a))')
     &          gen,spinc(gen),
     &          ': incoming photon corrections included'
          elseif(spinc(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: spinc = ',spinc(gen)
            stop
          endif
c incoming gluon  corrections
          if(sginc(gen).eq.1)then
            write(nout,'("       sginc(",i1,")=",i2,(a))')
     &          gen,sginc(gen),
     &          ': incoming gluon corrections included'
          elseif(sginc(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: sginc = ',sginc(gen)
            stop
          endif
c slicing terms
          if(ssli(gen).eq.1)then
            write(nout,'("        ssli(",i1,")=",i2,(a))')gen,ssli(gen),
     &          ': slicing terms included'
          elseif(ssli(gen).ne.0)then
            write(*,'(a)')' Unknown input: ssli = ',ssli(gen)
            stop
          endif
c subtraction terms
          if(ssub(gen).eq.1.or.ssub(gen).eq.2)then
            write(nout,'("        ssub(",i1,")=",i2,(a))')gen,ssub(gen),
     &          ': subtraction terms included'
          elseif(ssub(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: ssub = ',ssub(gen)
            stop
          endif
c collinear effective factors
          if(seff(gen).eq.1)then
            write(nout,'("        seff(",i1,")=",i2,(a))')gen,seff(gen),
     &          ': collinear effective factors included'
          elseif(seff(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: sseff = ',seff(gen)
            stop
          endif
c end-point contributions
          if(sendp(gen).eq.1)then
            write(nout,'("       sendp(",i1,")=",i2,(a))')
     &          gen,sendp(gen),
     &          ': end-point contributions included in Monte Carlo'
          elseif(sendp(gen).eq.-1)then
            write(nout,'("       sendp(",i1,")=",i2,(a))')
     &          gen,sendp(gen),
     &          ': end-point contributions included in |M|^2 '
          elseif(sendp(gen).ne.0)then
            write(*,'(a30,i3)')' Unknown input: sendp = ',sendp(gen)
            stop
          endif
c include cuts in phase-space generator
          if(sincludecuts(gen).eq.1)then
            write(nout,'("sincludecuts(",i1,")=",i2,(a))')gen,
     &          sincludecuts(gen),
     &          ': threshold cuts taken into account'//
     &          ' in phase-space generation'
          elseif(sincludecuts(gen).eq.2)then
            write(nout,'("sincludecuts(",i1,")=",i2,(a))')gen,
     &          sincludecuts(gen),
     &          ': threshold cuts taken into account'//
     &          ' in phase-space and x generation'
          elseif(sincludecuts(gen).eq.3)then
            write(nout,'("sincludecuts(",i1,")=",i2,(a))')gen,
     &          sincludecuts(gen),
     &          ': cuts taken into account in phase-space generation'
          elseif(sincludecuts(gen).eq.0)then
            write(nout,'("sincludecuts(",i1,")=",i2,(a))')gen,
     &          sincludecuts(gen),
     &          ': cuts NOT taken into account in phase-space'//
     &          ' generation'
          else
            write(*,'(a30,i3)')' Unknown input: sincludecuts = ',
     &          sincludecuts(gen)
            stop
          endif
        endif
       endif

       if(lnoutmc.ge.2.or.lwriterec)then
c recombination cuts
        if(selprocess.eq.0)then
        if(srecomb(gen).eq.1.or.srecomb(gen).eq.2)then
          if(lnoutmc.ge.2.) then
            write(nout,'("     srecomb(",i1,")=",i2,(a))')
     &        gen,srecomb(0),': with kT jet clustering:'
          endif
          if(lwriterec) then
            lwriterec=.false.
          endif
          if(ktpower.eq.1)then
            write(nout,'(" kt algorithm")')
          elseif(ktpower.eq.-1)then
            write(nout,'(" anti-kt algorithm")')
          elseif(ktpower.eq.0)then
            write(nout,'(" Cambridge-Aachen algorithm")')
          else
            write(nout,'(" type of kt algorithm undefined:",i3)')
     &          ktpower
            stop
          endif
          write(nout,'(" resolution parameter D  = ",f11.5)')
     &      ktd
          write(nout,'(" partonic cut |y_parton| < ",f11.5)')
     &      partetacut
        else if(srecomb(gen).eq.3.or.srecomb(gen).eq.4)then
          write(nout,'("     srecomb(",i1,")=",i2,(a))')gen,srecomb(0),
     &      ': with photon recombination:'
          if(srecomb(gen).eq.4) then
            write(nout,'(" theta_p    = ",f11.5," deg")')
     &          dacos(recparam(1,gen))*180d0/pi
            write(nout,'(" E_p        = ",f11.5," GeV")')
     &          recparam(2,gen)
            write(nout,'(" theta_pf   = ",f11.5," deg")')
     &          dacos(recparam(3,gen))*180d0/pi
            write(nout,'(" M_pf       = ",f11.5," GeV")')
     &          sqrt(recparam(4,gen))
          elseif (srecomb(gen).eq.3)then
            write(nout,'(" y_p        = ",f11.5)') recparam(1,gen)
            write(nout,'(" E_p        = ",f11.5," GeV")')
     &          recparam(2,gen)
            write(nout,'(" dR_pj      = ",f11.5)')
     &          sqrt(recparam(3,gen))
            write(nout,'(" M_pj       = ",f11.5," GeV")')
     &          sqrt(recparam(4,gen))
          endif
        elseif(srecomb(gen).ne.0)then
          write(*,*)' Unknown input: srecomb= ',srecomb(gen)
          stop
        endif
        elseif(selprocess.ge.1.and.selprocess.le.3)then
          if(sbarelep.eq.1)then
             write(nout,'(a16,i2,(a))')'sbarelep=',sbarelep,
     &            ': no lepton recombination with photons:'//
     &            ' bare leptons'
          elseif(sbarelep.eq.0)then
             write(nout,'(a16,i2,(a))')'sbarelep=',sbarelep,
     &            ': leptons are recombined with photons'
             write(nout,'(" recombination parameter R = ",f11.5)')
     &       dgammaparameter
          endif
          if(lwriterec) then
            lwriterec=.false.
          endif
          if(lnoutmc.ge.2.) then
            if(colllep.eq.1)then
               write(nout,'(a16,i2,(a))')'colllep=',colllep,
     &              ': initial gamma -> f fbar splitting included'
            elseif(colllep.eq.0.and.lnoutmc.ge.3)then
               write(nout,'(a16,i2,(a))')'colllep=',colllep,
     &              ': no initial gamma -> f fbar splitting included'
            endif
          endif
        endif
       endif

       if(lnoutmc.ge.2.or.lwritecuts)then
c separation cuts
        if(selprocess.eq.0)then
        if(scuts(gen).ge.1.and.scuts(gen).le.2)then
          if(lnoutmc.ge.2) then
            write(nout,'("       scuts(",i1,")=",i2,(a))')
     &        gen,scuts(gen),': with separation cuts:'
          else
            write(nout,'("          scuts=",i2,(a))') scuts(gen),
     &        ': with separation cuts:'
          endif
          if(lwritecuts) then
            lwritecuts=.false.
          endif
          if(scuts(gen).eq.1)then
            write(nout,'((a))')
     &      ' tagging jets defined via pT ordering'
          elseif(scuts(gen).eq.2)then
            write(nout,'((a))')
     &      ' tagging jets defined via E ordering'
          endif
          if(ptcutj1.gt.0d0)then
            write(nout,'(" pT(jet1)   > ",f11.5," GeV")') ptcutj1
          endif
          if(ycutj1.lt.1d10)then
            write(nout,'(" y(jet1)    < ",f11.5)') ycutj1
          endif
          if(ptcutj2.gt.0d0)then
            write(nout,'(" pT(jet2)   > ",f11.5," GeV")') ptcutj2
          endif
          if(ycutj2.lt.1d10)then
            write(nout,'(" y(jet2)    < ",f11.5)') ycutj2
          endif
          if(mlcutjj.gt.0d0)then
            write(nout,'(" m(j1,j2)   > ",f11.5," GeV")') mlcutjj
          endif
          if(mucutjj.lt.1d10)then
            write(nout,'(" m(j1,j2)   < ",f11.5," GeV")') mucutjj
          endif
          if(dycutjj.gt.0d0)then
            write(nout,'(" dy(j1,j2)  > ",f11.5)') dycutjj
          endif
          if(ptcutvis.gt.0d0)then
            write(nout,'(" pt(j1+j2+H) > ",f11.5," GeV")')
     &          ptcutvis
          endif
          if(ptmaxvis.lt.1d10)then
            write(nout,'(" pt(j1+j2+H) < ",f11.5," GeV")')
     &          ptmaxvis
          endif
          if(dphicutj1h.gt.0d0)then
            write(nout,'(" phi(j1,H)  > ",f11.5," rad")') dphicutj1h
          endif
          if(dphicutj2h.gt.0d0)then
            write(nout,'(" phi(j2,H)  > ",f11.5," rad")') dphicutj2h
          endif
          if(ptmaxj1.lt.1d10.and.yminj1.gt.0d0)then
            write(nout,'(" pT(jet1)   < ",f11.5" GeV  or ",
     &          " y(jet1)    > ",f11.5,"  (no visible jets)")')
     &          ptmaxj1,yminj1
          endif
          if(ptmaxj2.lt.1d10.and.yminj2.gt.0d0)then
            write(nout,'(" pT(jet2)   < ",f11.5" GeV  or ",
     &          " y(jet2)    > ",f11.5,"  (no visible 2nd jet)")')
     &          ptmaxj2,yminj2
          endif
          if(ptmaxj3.lt.1d10.and.yminj3.gt.0d0)then
            write(nout,'(" pT(jet3)   < ",f11.5" GeV  or ",
     &          " y(jet3)    > ",f11.5,"  (no visible 3rd jet)")')
     &          ptmaxj3,yminj3
          endif
          if(mucutjjh.lt.1d10)then
            write(nout,'(" mass(j,j,h)< ",f11.5" GeV")') mucutjjh
          endif
          if(mlcutjjh.gt.0d0)then
            write(nout,'(" mass(j,j,h)> ",f11.5" GeV")') mlcutjjh
          endif
          if(ptcuth.gt.0d0)then
             write(nout,'(" pT(Higgs)  > ",f11.5" GeV")') ptcuth
          endif
          if(ycuth.lt.1d10)then
             write(nout,'(" y(Higgs)   < ",f11.5)') ycuth
          endif
          if(ecuth.gt.0d0)then
             write(nout,'(" E(Higgs)   > ",f11.5" GeV")') ecuth
          endif
          if((shtr.ge.1))then
            write(nout,'(" sqrt(s_H)  > ",f11.5, " GeV")') mlcuth
            if(mucuth.lt.1e29)
     &          write(nout,'(" sqrt(s_H)  < ",f11.5, " GeV")') mucuth
          endif
          if((scuts(gen).eq.1.or.scuts(gen).eq.2).and.lhemicut)then
            write(nout,'(" hemisphere cut:  y(jet1)*y(jet2) < 0")')
          endif
          if((shtr.ge.1))then
            if(ptcuts.gt.0d0)then
              write(nout,'(" pT(decp)   > ",f11.5," GeV")') ptcuts
            endif
            if(ycuts.lt.1d10)then
              write(nout,'(" y(decp)    < ",f11.5)') ycuts
            endif
          endif
          if(shtr.eq.2.and.(scuts(gen).eq.1.or.scuts(gen).eq.2)
     &        .and.lhiggsbjets)then
            write(nout,
     &        '(" rapidity cut:    y(jet_min) < y(decp) < y(jet_max)")')
          endif
        elseif(scuts(gen).eq.0)then
          if((shtr.ge.1))then
            write(nout,'(" sqrt(s_H)  > ",f11.5, " GeV")') mlcuth
            if(mucuth.lt.1e29)
     &          write(nout,'(" sqrt(s_H)  < ",f11.5, " GeV")') mucuth
          endif
        else
          write(*,'(a30,i3)')' Unknown input: scuts = ',scuts(gen)
          stop
        endif
        else  ! not jjh
          if(scuts(gen).ge.1.and.scuts(gen).le.2)then
            if(lnoutmc.ge.2) then
              if(selprocess.eq.1) then
                write(nout,'("       scuts(",i1,")=",i2,(a))')
     &            gen,scuts(gen),': with separation cuts: (3=l+) '
              elseif(selprocess.eq.2) then
                write(nout,'("       scuts(",i1,")=",i2,(a))')
     &            gen,scuts(gen),': with separation cuts: (4=l-) '
              elseif(selprocess.eq.3) then
                write(nout,'("       scuts(",i1,")=",i2,(a))')
     &            gen,scuts(gen),': with separation cuts: (3=l+,4=l-) '
              elseif(selprocess.eq.4) then
                write(nout,'("       scuts(",i1,")=",i2,(a))')
     &            gen,scuts(gen),': with separation cuts:'
              endif
            else
              if(selprocess.eq.1) then
                write(nout,'("          scuts=",i2,(a))') scuts(gen),
     &            ': with separation cuts: (3=l+) '
              elseif(selprocess.eq.2) then
                write(nout,'("          scuts=",i2,(a))') scuts(gen),
     &            ': with separation cuts: (4=l-) '
              elseif(selprocess.eq.3) then
                write(nout,'("          scuts=",i2,(a))') scuts(gen),
     &            ': with separation cuts: (3=l+, 4=l-) '
              elseif(selprocess.eq.4) then
                write(nout,'("          scuts=",i2,(a))') scuts(gen),
     &            ': with separation cuts: '
              endif
            endif
            if(lwritecuts) then
              lwritecuts=.false.
            endif
            do ext=1,nextgen(gen)
              if(ptcut(ext,gen).gt.0d0
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu')then
                write(nout,'(" pT(",i1,")      > ",f11.5," GeV")')
     &            ext,ptcut(ext,gen)
              endif
            enddo
            do ext=1,nextgen(gen)
              if(ecut(ext,gen).gt.0d0
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu')then
                write(nout,'(" energy(",i1,")  > ",f11.5," GeV")')
     &            ext,ecut(ext,gen)
              endif
            enddo
            do ext=1,nextgen(gen)
              if(ycut(ext,gen).lt.1d10
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu')then
                write(nout,'(" y(",i1,")       < ",f11.5)')
     &            ext,ycut(ext,gen)
              endif
            enddo
            do ext=1,nextgen(gen)
            do ext2=ext+1,nextgen(gen)
              if(ccut(ext,ext2,gen).lt.1d0
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu'
     &           .and.name(ext2,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext2,gen)(1:2).ne.'nu')then
                write(nout,'(" angle(",i1,",",i1,") > ",f11.5," deg")')
     &            ext,ext2,dacos(ccut(ext,ext2,gen))/pi*180d0
              endif
            enddo
            enddo
            do ext=1,nextgen(gen)
            do ext2=ext+1,nextgen(gen)
              if(scut(ext,ext2,gen).gt.0d0
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu'
     &           .and.name(ext2,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext2,gen)(1:2).ne.'nu')then
                write(nout,'(" mass(",i1,",",i1,")  > ",f11.5," GeV")')
     &            ext,ext2,sqrt(scut(ext,ext2,gen))
              endif
            enddo
            enddo
            do ext=3,nextgen(gen)
            do ext2=ext+1,nextgen(gen)
              if(dycut(ext,ext2,gen).gt.0d0
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu'
     &           .and.name(ext2,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext2,gen)(1:2).ne.'nu')then
                write(nout,'(" dy(",i1,",",i1,")    > ",f11.5)')
     &            ext,ext2,dycut(ext,ext2,gen)
              endif
            enddo
            enddo
            do ext=3,nextgen(gen)
            do ext2=ext+1,nextgen(gen)
              if(dRcut(ext,ext2,gen).gt.0d0
     &           .and.name(ext,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext,gen)(1:2).ne.'nu'
     &           .and.name(ext2,gen)(1:7).ne.'anti-nu'
     &           .and.name(ext2,gen)(1:2).ne.'nu')then
                write(nout,'(" dR(",i1,",",i1,")    > ",f11.5)')
     &            ext,ext2,dRcut(ext,ext2,gen)
              endif
            enddo
            enddo
            if(dRcutjl.gt.0d0)then
                 write(nout,'(" dR(j,l))   > ",f11.5)') dRcutjl
            endif
            if(mucutll.lt.1d10)then
                 write(nout,'(" mass(l,l)  < ",f11.5" GeV")') mucutll
            endif
            if(mucutllh.lt.1d10)then
                 write(nout,'(" mass(l,l,h)< ",f11.5" GeV")') mucutllh
            endif
            if(mlcutllh.gt.0d0)then
                 write(nout,'(" mass(l,l,h)> ",f11.5" GeV")') mlcutllh
            endif
            if(ptcutv.gt.0d0)then
                 write(nout,'(" pT(V)      > ",f11.5" GeV")') ptcutv
            endif
            if(ptcuth.gt.0d0)then
                 write(nout,'(" pT(Higgs)  > ",f11.5" GeV")') ptcuth
            endif
            if(ptcutm.gt.0d0.and.selprocess.ne.3)then
                 write(nout,'(" pT(miss)   > ",f11.5" GeV")') ptcutm
            endif
            if(ptmaxl.lt.1d10.and.yminl.gt.0d0)then
               write(nout,'(" pT(lep)    < ",f11.5" GeV or ",
     &          " y(lep)     > ",f11.5,"  (no visible lepton)")')
     &          ptmaxl,yminl
            endif
            if(ycuth.lt.1d10)then
              write(nout,'(" y(Higgs)   > ",f11.5" GeV")') ycuth
            endif
            if(ecuth.gt.0d0)then
              write(nout,'(" E(Higgs)   > ",f11.5" GeV")') ecuth
            endif
            if((shtr.ge.1))then
              write(nout,'(" sqrt(s_H)  > ",f11.5, " GeV")') mlcuth
              if(mucuth.lt.1e29)
     &            write(nout,'(" sqrt(s_H)  < ",f11.5, " GeV")') mucuth
            endif
          elseif(scuts(gen).eq.0)then
            if((shtr.ge.1))then
              write(nout,'(" sqrt(s_H)  > ",f11.5, " GeV")') mlcuth
              if(mucuth.lt.1e29)
     &            write(nout,'(" sqrt(s_H)  < ",f11.5, " GeV")') mucuth
            endif
          else
            write(*,'(a30,i3)')' Unknown input: scuts = ',scuts(gen)
            stop
          endif
        endif
       endif
      endif
      enddo
      write(nout,*)

      end

***********************************************************************
      subroutine writeresult
***********************************************************************
*     output of result                                                *
*     adapted from COFFERGG by Markus Roth                            *
*                                                                     *
*---------------------------------------------------------------------*
*     21.02.06  Ansgar Denner     last changed  18.06.07              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

#ifdef multpdf
#include "multpdf.h"
      integer nset
#endif
#ifdef mpiuse
      include "mpif.h"
      integer status(MPI_STATUS_SIZE)
      real*8 keep_sigmaa
      real*8 keep_avborn
      real*8 keep_siborn
#ifdef multpdf
      real*8 keep_mavborn(maxpdf)
      real*8 keep_msiborn(maxpdf)
#endif
      real*8 keep_avbornii(-4:4,-4:4)
      real*8 keep_sibornii(-4:4,-4:4)
      real*8 keep_averageii(-4:4,-4:4)
      real*8 keep_sigmaii(-4:4,-4:4)
      real*8 keep_avgenii(0:maxg,maxcon,-4:4,-4:4)
      real*8 keep_sigenii(0:maxg,maxcon,-4:4,-4:4)
      real*8 keep_avgencha(0:maxg,maxcon,maxch)
      real*8 keep_sigencha(0:maxg,maxcon,maxch)
      real*8 keep_avgen(0:maxg,maxcon)
      real*8 keep_sigen(0:maxg,maxcon)
#ifdef multpdf
      real*8 keep_mavgen(maxpdf,0:maxg,maxcon)
      real*8 keep_msigen(maxpdf,0:maxg,maxcon)
#endif
      real*8 keep_weighttotmax
      real*8 keep_weightmaxcha(0:maxg,maxcon,maxch)
      real*8 keep_weightmax(0:maxg,maxcon)
      integer keep_nmax(0:maxg,maxcon)
      integer keep_channelmax(0:maxg,maxcon)
      integer keep_rejcha(0:maxg,maxcon,maxch)
      integer keep_rejgen(0:maxg,maxcon)
#endif

      integer iNumProcs, iErr, iMyRank,ii
      common/mpistuff/iNumProcs,iMyRank

c local variables
       character*5 part(-4:4)
       integer i1,i2
      real*8 avgenchan(0:maxg,maxcon,maxch)
      real*8 sigenchan(0:maxg,maxcon,maxch)


c date and time
      real*8     rcpu
      integer*4  int,icpu,icpuh,icpum,icpus
      integer*4  count,nout2
      character  datum*10,da1*2,da2*2,da3*4
      character  chtime*8
      integer    iarray(3)

#ifdef mpiuse
      if(iMyRank.EQ.0.AND.iNumProcs.NE.1) Then
	do ii=1,iNumProcs-1
	  call MPI_Recv(keep_sigmaa, 1,
     &    MPI_DOUBLE_PRECISION, ii, 10, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_avborn, 1,
     &    MPI_DOUBLE_PRECISION, ii, 11, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_siborn, 1,
     &    MPI_DOUBLE_PRECISION, ii, 12, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_avbornii, 81,
     &    MPI_DOUBLE_PRECISION, ii, 13, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_sibornii, 81,
     &    MPI_DOUBLE_PRECISION, ii, 14, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_averageii, 81,
     &    MPI_DOUBLE_PRECISION, ii, 15, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_sigmaii, 81,
     &    MPI_DOUBLE_PRECISION, ii, 16, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_avgenii, 81*(maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 17, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_sigenii, 81*(maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 18, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_avgencha, (maxg+1)*(maxcon)*(maxch),
     &    MPI_DOUBLE_PRECISION, ii, 19, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_sigencha, (maxg+1)*(maxcon)*(maxch),
     &    MPI_DOUBLE_PRECISION, ii, 20, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_avgen, (maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 21, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_sigen, (maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 22, MPI_COMM_WORLD, status, iErr)
#ifdef multpdf
	  call MPI_Recv(keep_mavborn, maxpdf,
     &    MPI_DOUBLE_PRECISION, ii, 11, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_msiborn, maxpdf,
     &    MPI_DOUBLE_PRECISION, ii, 12, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_mavgen, maxpdf*(maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 30, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_msigen, maxpdf*(maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 31, MPI_COMM_WORLD, status, iErr)
#endif
	  call MPI_Recv(keep_weighttotmax, 1,
     &    MPI_DOUBLE_PRECISION, ii, 23, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_weightmaxcha, (maxg+1)*(maxcon)*(maxch),
     &    MPI_DOUBLE_PRECISION, ii, 24, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_weightmax, (maxg+1)*(maxcon),
     &    MPI_DOUBLE_PRECISION, ii, 25, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_nmax, (maxg+1)*(maxcon),
     &    MPI_INTEGER, ii, 26, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_channelmax, (maxg+1)*(maxcon),
     &    MPI_INTEGER, ii, 27, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_rejcha, (maxg+1)*(maxcon)*(maxch),
     &    MPI_INTEGER, ii, 28, MPI_COMM_WORLD, status, iErr)
	  call MPI_Recv(keep_rejgen, (maxg+1)*(maxcon),
     &    MPI_INTEGER, ii, 29, MPI_COMM_WORLD, status, iErr)
	  sigmaa=sigmaa+keep_sigmaa
	  avborn=avborn+keep_avborn
	  siborn=siborn+keep_siborn
#ifdef multpdf
	  do nset=1,pdfmemberto-pdfmemberfrom+1
	    mavborn(nset)=mavborn(nset)+keep_mavborn(nset)
	    msiborn(nset)=msiborn(nset)+keep_msiborn(nset)
	  enddo
#endif
	  do i1=-4,4
	   do i2=-4,4
	    avbornii(i1,i2)=avbornii(i1,i2)+keep_avbornii(i1,i2)
	    sibornii(i1,i2)=sibornii(i1,i2)+keep_sibornii(i1,i2)
	    averageii(i1,i2)=averageii(i1,i2)+keep_averageii(i1,i2)
	    sigmaii(i1,i2)  =sigmaii(i1,i2)  +keep_sigmaii(i1,i2)
            do gen=0,ngenerator
             do icon=1,ncontr(gen)
	      avgenii(gen,icon,i1,i2)=avgenii(gen,icon,i1,i2)+
     &                keep_avgenii(gen,icon,i1,i2)
	      sigenii(gen,icon,i1,i2)=sigenii(gen,icon,i1,i2)+
     &                keep_sigenii(gen,icon,i1,i2)
	     enddo
	    enddo
	   enddo
	  enddo	  	  	
          do gen=0,ngenerator
           do icon=1,ncontr(gen)
            do cha=1,nchannel(gen)
	     avgencha(gen,icon,cha)=avgencha(gen,icon,cha)+
     &                              keep_avgencha(gen,icon,cha)
	     sigencha(gen,icon,cha)=sigencha(gen,icon,cha)+
     &                              keep_sigencha(gen,icon,cha)
	    enddo
	    avgen(gen,icon)=avgen(gen,icon)+keep_avgen(gen,icon)
	    sigen(gen,icon)=sigen(gen,icon)+keep_sigen(gen,icon)
#ifdef multpdf
	    do nset=1,pdfmemberto-pdfmemberfrom+1
	      mavgen(nset,gen,icon)=mavgen(nset,gen,icon)+
     &                              keep_mavgen(nset,gen,icon)
	      msigen(nset,gen,icon)=msigen(nset,gen,icon)+
     &                              keep_msigen(nset,gen,icon)
            enddo
#endif
           enddo
	  enddo
	  If(keep_weighttotmax.GT.weighttotmax)
     &                           weighttotmax=keep_weighttotmax
          do gen=0,ngenerator
           do icon=1,ncontr(gen)
            do cha=1,nchannel(gen)
	      If(keep_weightmaxcha(gen,icon,cha).GT.
     &            weightmaxcha(gen,icon,cha)) Then
	        weightmaxcha(gen,icon,cha)=
     &                           keep_weightmaxcha(gen,icon,cha)
	      Endif	
	      rejcha(gen,icon,cha)=rejcha(gen,icon,cha)+
     &                             keep_rejcha(gen,icon,cha)	
	    enddo
	     If(keep_weightmax(gen,icon).GT.weightmax(gen,icon)) Then
	       weightmax(gen,icon)=keep_weightmax(gen,icon)
	       nmax(gen,icon)=keep_nmax(gen,icon)	
	       channelmax(gen,icon)=keep_channelmax(gen,icon)	
	     Endif	
	     rejgen(gen,icon)=rejgen(gen,icon)+keep_rejgen(gen,icon)
	   enddo
	  enddo
	enddo
      else if(iMyRank.NE.0.AND.iNumProcs.NE.1) Then
	call MPI_Send(sigmaa, 1,
     &  MPI_DOUBLE_PRECISION, 0, 10, MPI_COMM_WORLD, iErr)
	call MPI_Send(avborn, 1,
     &  MPI_DOUBLE_PRECISION, 0, 11, MPI_COMM_WORLD, iErr)
	call MPI_Send(siborn, 1,
     &  MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, iErr)
	call MPI_Send(avbornii, 81,
     &  MPI_DOUBLE_PRECISION, 0, 13, MPI_COMM_WORLD, iErr)
	call MPI_Send(sibornii, 81,
     &  MPI_DOUBLE_PRECISION, 0, 14, MPI_COMM_WORLD, iErr)
	call MPI_Send(averageii, 81,
     &  MPI_DOUBLE_PRECISION, 0, 15, MPI_COMM_WORLD, iErr)
	call MPI_Send(sigmaii, 81,
     &  MPI_DOUBLE_PRECISION, 0, 16, MPI_COMM_WORLD, iErr)
	call MPI_Send(avgenii, 81*(maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 17, MPI_COMM_WORLD, iErr)
	call MPI_Send(sigenii, 81*(maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 18, MPI_COMM_WORLD, iErr)
	call MPI_Send(avgencha, (maxg+1)*(maxcon)*(maxch),
     &  MPI_DOUBLE_PRECISION, 0, 19, MPI_COMM_WORLD, iErr)
	call MPI_Send(sigencha, (maxg+1)*(maxcon)*(maxch),
     &  MPI_DOUBLE_PRECISION, 0, 20, MPI_COMM_WORLD, iErr)
	call MPI_Send(avgen, (maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 21, MPI_COMM_WORLD, iErr)
	call MPI_Send(sigen, (maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 22, MPI_COMM_WORLD, iErr)
#ifdef multpdf
	call MPI_Send(mavborn, maxpdf,
     &  MPI_DOUBLE_PRECISION, 0, 11, MPI_COMM_WORLD, iErr)
	call MPI_Send(msiborn, maxpdf,
     &  MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, iErr)
	call MPI_Send(mavgen, maxpdf*(maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 30, MPI_COMM_WORLD, iErr)
	call MPI_Send(msigen, maxpdf*(maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 31, MPI_COMM_WORLD, iErr)
#endif
	call MPI_Send(weighttotmax, 1,
     &  MPI_DOUBLE_PRECISION, 0, 23, MPI_COMM_WORLD, iErr)
	call MPI_Send(weightmaxcha, (maxg+1)*(maxcon)*(maxch),
     &  MPI_DOUBLE_PRECISION, 0, 24, MPI_COMM_WORLD, iErr)
	call MPI_Send(weightmax, (maxg+1)*(maxcon),
     &  MPI_DOUBLE_PRECISION, 0, 25, MPI_COMM_WORLD, iErr)
	call MPI_Send(nmax, (maxg+1)*(maxcon),
     &  MPI_INTEGER, 0, 26, MPI_COMM_WORLD, iErr)
	call MPI_Send(channelmax, (maxg+1)*(maxcon),
     &  MPI_INTEGER, 0, 27, MPI_COMM_WORLD, iErr)
	call MPI_Send(rejcha, (maxg+1)*(maxcon)*(maxch),
     &  MPI_INTEGER, 0, 28, MPI_COMM_WORLD, iErr)
	call MPI_Send(rejgen, (maxg+1)*(maxcon),
     &  MPI_INTEGER, 0, 29, MPI_COMM_WORLD, iErr)
      endif
      if(iMyRank.EQ.0) Then
#endif

      part(-4) = 'cbar'
      part(-3) = 'sbar'
      part(-2) = 'ubar'
      part(-1) = 'dbar'
      part(1) = 'd'
      part(2) = 'u'
      part(3) = 's'
      part(4) = 'c'

 1052 format(" Cross section          = ",1g24.14)
 1053 format(" Standard deviation     = ",1g24.14)
      if(nunw.ne.0)then
        write(nout,'(a)')' '
        write(nout,'(a)')' '
        write(nout,'(a)')' Unweighted events'
        write(nout,1008)nunw
        write(nout,1009)nneg
      endif

      average=average/dble(nevt)
      sigma=sigma/dble(nevt)**2
      sigmaa=sigmaa/dble(nevt)**2
      sigma=sqrt(dabs(sigma-average**2/dble(nevt)))
      sigmaa=sqrt(dabs(sigmaa-average**2/dble(nevt)))
      avborn=avborn/dble(nevt)
      siborn=siborn/dble(nevt)**2
      siborn=sqrt(dabs(siborn-avborn**2/dble(nevt)))
#ifdef multpdf
      do nset=1,pdfmemberto-pdfmemberfrom+1
	maverage(nset)=maverage(nset)/dble(nevt)
	msigma(nset)=msigma(nset)/dble(nevt)**2
	msigmaa(nset)=msigmaa(nset)/dble(nevt)**2
	msigma(nset)=sqrt(dabs(msigma(nset)-maverage(nset)**2/
     &                     dble(nevt)))
	msigmaa(nset)=sqrt(dabs(msigmaa(nset)-maverage(nset)**2/
     &                     dble(nevt)))
	mavborn(nset)=mavborn(nset)/dble(nevt)
	msiborn(nset)=msiborn(nset)/dble(nevt)**2
	msiborn(nset)=sqrt(dabs(msiborn(nset)-
     &                     mavborn(nset)**2/dble(nevt)))
      enddo
#endif
#ifdef SUBPROCESSES
      do i1=-4,4
      do i2=-4,4
        avbornii(i1,i2)=avbornii(i1,i2)/dble(nevt)
        sibornii(i1,i2)=sibornii(i1,i2)/dble(nevt)**2
        sibornii(i1,i2)=sqrt(dabs(sibornii(i1,i2)
     &      -avbornii(i1,i2)**2/dble(nevt)))
        averageii(i1,i2)=averageii(i1,i2)/dble(nevt)
        sigmaii(i1,i2)=sigmaii(i1,i2)/dble(nevt)**2
        sigmaii(i1,i2)=sqrt(dabs(sigmaii(i1,i2)
     &      -averageii(i1,i2)**2/dble(nevt)))
        do gen=0,ngenerator
          do icon=1,ncontr(gen)
            avgenii(gen,icon,i1,i2)=avgenii(gen,icon,i1,i2)/dble(nevt)
            sigenii(gen,icon,i1,i2)=
     &          sigenii(gen,icon,i1,i2)/dble(nevt)**2
            sigenii(gen,icon,i1,i2)=sqrt(dabs(sigenii(gen,icon,i1,i2)
     &          -avgenii(gen,icon,i1,i2)**2/dble(nevt)))
          enddo
        enddo
      enddo
      enddo
#endif
      do gen=0,ngenerator
        do icon=1,ncontr(gen)
          do cha=1,nchannel(gen)
            if (nevtgencha(gen,icon,cha).gt.0) then
              avgencha(gen,icon,cha)=
     &            avgencha(gen,icon,cha)/nevtgencha(gen,icon,cha)
              sigencha(gen,icon,cha)=sigencha(gen,icon,cha)/
     &            dble(nevtgencha(gen,icon,cha))**2
              sigencha(gen,icon,cha)=sqrt(dabs(sigencha(gen,icon,cha)
     &            -avgencha(gen,icon,cha)**2/nevtgencha(gen,icon,cha)))
c     normalize to all events
              avgenchan(gen,icon,cha)=
     &          avgencha(gen,icon,cha)*nevtgencha(gen,icon,cha)/nevents
              sigenchan(gen,icon,cha)=sigencha(gen,icon,cha)*
     &            dble(nevtgencha(gen,icon,cha))/nevents
              weightmaxcha(gen,icon,cha)=weightmaxcha(gen,icon,cha)
            endif
          enddo
          avgen(gen,icon)=avgen(gen,icon)/dble(nevt)
          sigen(gen,icon)=sigen(gen,icon)/dble(nevt)**2
          sigen(gen,icon)=sqrt(dabs(sigen(gen,icon)
     &        -avgen(gen,icon)**2/dble(nevt)))

#ifdef multpdf
      do nset=1,pdfmemberto-pdfmemberfrom+1
          mavgen(nset,gen,icon)=mavgen(nset,gen,icon)/dble(nevt)
          msigen(nset,gen,icon)=msigen(nset,gen,icon)/dble(nevt)**2
          msigen(nset,gen,icon)=sqrt(dabs(msigen(nset,gen,icon)
     &        -mavgen(nset,gen,icon)**2/dble(nevt)))
      enddo
#endif
        enddo
      enddo




 1000 format(" Total cross section    = ",2g18.10)
 1001 format(" Number of events       = ",i10,:" of ",i10)
 1002 format(" Cross section          = ",2g18.10)
 1003 format(" Standard deviation     = ",2g18.10)
 1004 format(" Relative corrections   = ",1f18.10)
 1005 format(" Maximal weight         = ",1f18.10)
 1006 format(" Rejected events        = ",i10,1f17.4,"%")
 1007 format(" (event = ",i9,", channel = ",i8,")")
 1017 format(" (event = ",i9,")")
 1008 format(" Number of unweighted events          = ",i10)
 1009 format(" Number of negative unweighted events = ",i10)

 1021 format(1x,a5,a5,": Cross section          = ",1g18.10)
 1022 format(12x," Standard deviation     = ",1g18.10)

      if(lnoutmc.ge.2)then
      write(nout,'(a)')' '
      write(nout,'(a)')' Results:'
      write(nout,'(a)')' --------'
      write(nout,1001)nevt
      write(nout,1000)average
      write(nout,1003)sigma       ! from sum_n (sum_icon w_icon)^2
      write(nout,1003)sigmaa      ! from sum_n sum_icon (w_icon)^2
      write(nout,1005)weighttotmax/dble(nevt)
      if(average.ne.0d0.and.avborn.ne.0d0)then
        ngen=0
        do gen=0,ngenerator
          ngen=ngen+sborn(gen)
        enddo
        write(nout,1004)average/avborn-dble(ngen)
        write(nout,1003)dabs(average/avborn)
     &    *sqrt((sigma/average)**2+(siborn/avborn)**2)
      endif
      write(nout,'(a)')' '
      write(nout,'(a)')' '
      write(nout,'(a)')' Born contribution'
      write(nout,1000)avborn
      write(nout,1003)siborn
      write(nout,*)

#ifdef SUBPROCESSES
      do i1=-4,4
      do i2=-4,4
      if(i1*i2.ne.0d0) then
        write(nout,1021)part(i1),part(i2),avbornii(i1,i2)
        write(nout,1022)sibornii(i1,i2)
      endif
      enddo
      enddo
#endif

      write(nout,'(a)')' '
      write(nout,'(a)')' '
      write(nout,'(a)')' Complete results'
      write(nout,1000)average
      write(nout,1003)sigma
      write(nout,1003)sigmaa
      write(nout,*)

#ifdef SUBPROCESSES
      do i1=-4,4
      do i2=-4,4
      if(i1*i2.ne.0d0.and.averageii(i1,i2).ne.0d0) then
        write(nout,1021)part(i1),part(i2),averageii(i1,i2)
        write(nout,1022)sigmaii(i1,i2)
      endif
      enddo
      enddo
      i1=0
      part(i1)='p/gl'
      do i2=-4,4
      if(i2.ne.0d0.and.averageii(0,i2).ne.0d0) then
        write(nout,1021)part(i1),part(i2),averageii(i1,i2)
        write(nout,1022)sigmaii(i1,i2)
      endif
      enddo
      i2=0
      part(i2)='p/gl'
      do i1=-4,4
      if(i1.ne.0d0.and.averageii(i1,0).ne.0d0) then
        write(nout,1021)part(i1),part(i2),averageii(i1,i2)
        write(nout,1022)sigmaii(i1,i2)
      endif
      enddo
#endif


      write(nout,'(a)')' '
      do gen=0,ngenerator
        do icon=1,ncontr(gen)
          if(sgencon(gen,icon).ne.0)then
            write(nout,'(a)')' '
            write(nout,'(a23,i2)')' Phase-space generator ',gen
            write(nout,'(a23,i2," of ",i2)')
     &          ' contribution ',icon,ncontr(gen)
            write(nout,'(" Process: ",a60)')processname(gen)
            write(nout,1001)nevt/freqgencon(gen,icon),nevt
            write(nout,1000)avgen(gen,icon)
            write(nout,1003)sigen(gen,icon)
            write(nout,1006)rejgen(gen,icon),
     &          1d2*rejgen(gen,icon)/dble(nevt)*freqgencon(gen,icon)
            write(nout,1005)weightmax(gen,icon)/dble(nevt)
            write(nout,1007)nmax(gen,icon),channelmax(gen,icon)

#ifdef SUBPROCESSES
            if(gen.le.2.and.icon.le.5) then
              do i1=-4,4
                do i2=-4,4
                  if(i1*i2.ne.0d0.and.avgenii(gen,icon,i1,i2).ne.0d0)
     &                then
                    write(nout,1021)part(i1),part(i2),
     &                  avgenii(gen,icon,i1,i2)
                    write(nout,1022)sigenii(gen,icon,i1,i2)
                  endif
                enddo
              enddo
            endif

            if(gen.eq.3.or.gen.eq.0.and.icon.eq.6) then
              part(0) = 'phot'
              i1=0
              do i2=-4,4
                if(i2.ne.0d0) then
                  write(nout,1021)part(i1),part(i2),
     &                avgenii(gen,icon,i1,i2)
                  write(nout,1022)sigenii(gen,icon,i1,i2)
                endif
              enddo
            endif
            if(gen.eq.4.or.gen.eq.0.and.icon.eq.7) then
              part(0) = 'phot'
              i2=0
              do i1=-4,4
                if(i1.ne.0d0) then
                  write(nout,1021)part(i1),part(i2),
     &                avgenii(gen,icon,i1,i2)
                  write(nout,1022)sigenii(gen,icon,i1,i2)
                endif
              enddo
            endif
            if(gen.eq.5.or.gen.eq.0.and.icon.eq.8) then
              part(0) = 'gluo'
              i1=0
              do i2=-4,4
                if(i2.ne.0d0) then
                  write(nout,1021)part(i1),part(i2),
     &                avgenii(gen,icon,i1,i2)
                  write(nout,1022)sigenii(gen,icon,i1,i2)
                endif
              enddo
            endif
            if(gen.eq.6.or.gen.eq.0.and.icon.eq.9) then
              part(0) = 'gluo'
              i2=0
              do i1=-4,4
                if(i1.ne.0d0) then
                  write(nout,1021)part(i1),part(i2),
     &                avgenii(gen,icon,i1,i2)
                  write(nout,1022)sigenii(gen,icon,i1,i2)
                endif
              enddo
            endif
#endif

            if(lnoutmc.ge.5) then
              do cha=1,nchannel(gen)
                write(nout,'(a)')' '
                write(nout,'(a23,i3)')' channel        ',cha
                write(nout,'(a23,g23.8)')' alpha_channel  ',
     &              alphaopt(cha,opt(gen),gen)
                write(nout,1001)nevtgencha(gen,icon,cha)
     &              /freqgencon(gen,icon)
                write(nout,1002)avgenchan(gen,icon,cha)
                write(nout,1003)sigenchan(gen,icon,cha)
                write(nout,1006)rejcha(gen,icon,cha),
     &              ((1d2*rejcha(gen,icon,cha))/dble(nevt)
     &              *freqgencon(gen,icon))
                write(nout,1005)weightmaxcha(gen,icon,cha)/dble(nevt)
                write(nout,1017)nmaxcha(gen,icon,cha)
              enddo
            endif
          else
            write(nout,'(a)')' '
            write(nout,'(a23,i2)')' Phase-space generator ',gen
            write(nout,'(a23,i2," of ",i2)')
     &          ' contribution ',icon,ncontr(gen)
            write(nout,'(" Process: ",a60)')processname(gen)
            write(nout,'(a)')' no contribution'
          endif
        enddo
      enddo
      endif

 301  format(/' Relative corrections')
 400  format(a12,i2,a14,i2,a1,1g18.10,' +/- ',1g18.10)
 500  format(a31,1g18.10,' +/- ',1g18.10)
      write(nout,*)
      write(nout,*)
     & 'Summary of results for ',initialstate,'-> ',finalstate
      write(nout,500)' Born cross section:',
     &    avborn,siborn
      write(nout,500)' complete cross section:',
     &    average,sigma
      if(lnoutmc.ge.2)then
      if (sbarelep.eq.0) then
       write(nout,500)' complete cross section:',
     &    (avgen(0,1)+avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))
     &    +(avgen(0,2)+avgen(0,3)+avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1)),
     &    sqrt((sigen(0,1))**2+(sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2)
      else
       write(nout,500)' complete cross section:',
     &     (avgen(0,1)+avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))
     &    +(avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1)),
     &    sqrt((sigen(0,1))**2+(sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2)
      endif
      endif
      write(nout,500)' QCD corr. to cross section:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)),
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2)
      write(nout,500)' incoming gluon cross section:',
     &    (avgen(0,8)+avgen(0,9)+avgen(5,1)+avgen(6,1)),
     &    sqrt((sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(5,1))**2+(sigen(6,1))**2)
      write(nout,500)' virt+real QCD cross section:',
     &    (avgen(0,5)+avgen(0,4)+avgen(2,1)),
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2+(sigen(2,1))**2)
      write(nout,500)' virt. QCD cross section:',
     &    avgen(0,5),sigen(0,5)
      write(nout,500)' outgoing gluon cross section:',
     &    (avgen(0,4)+avgen(2,1)),
     &    sqrt((sigen(0,4))**2+(sigen(2,1))**2)
      if (sbarelep.eq.0) then
       write(nout,500) 'ew corr. to cross section:',
     &    (avgen(0,2)+avgen(0,3)+avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1)),
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2)
      else
       write(nout,500) 'ew corr. to cross section:',
     &    (avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1)),
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &        +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2)
      endif
      write(nout,500) 'incoming photon cross section:',
     &    (avgen(0,6)+avgen(0,7)
     &     +avgen(3,1)+avgen(4,1)
     &     +avgen(7,1)+avgen(8,1)
     &     +avgen(9,1)+avgen(10,1)),
     &     sqrt( (sigen(0,6))**2+(sigen(0,7))**2
     &     +(sigen(3,1))**2+(sigen(4,1))**2
     &     +(sigen(7,1))**2+(sigen(8,1))**2
     &     +(sigen(9,1))**2+(sigen(10,1))**2)
      if (sbarelep.eq.0) then
       write(nout,500) ' virt+real ew cross section:',
     &    (avgen(0,2)+avgen(0,3)+avgen(1,1)),
     &    sqrt( (sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(1,1))**2)
      else
       write(nout,500) ' virt+real ew cross section:',
     &    (avgen(0,2)+avgen(0,3)+avgen(0,10)+avgen(1,1)),
     &    sqrt( (sigen(0,2))**2+(sigen(0,3))**2
     &        +(sigen(0,10))**2
     &        +(sigen(1,1))**2)
      endif
      write(nout,500) ' virt. ew cross section:',
     &    avgen(0,2),sigen(0,2)
      if (sbarelep.eq.0) then
       write(nout,500) 'outgoing photon cross section:',
     &    (avgen(0,3)+avgen(1,1)),
     &    sqrt((sigen(0,3))**2+(sigen(1,1))**2)
      else
       write(nout,500) 'outgoing photon cross section:',
     &    (avgen(0,3)+avgen(1,1)+avgen(0,10)),
     &    sqrt((sigen(0,3))**2+(sigen(1,1))**2+(sigen(0,10))**2)
      endif

      if(lnoutmc.ge.2)then
      write(nout,500)
      do gen=0,ngenerator
        do icon=1,ncontr(gen)
          if(avgen(gen,icon).ne.0d0)then
            write(nout,400)' generator',gen,
     &    ', contribution',icon,':',avgen(gen,icon),sigen(gen,icon)
          endif
        enddo
      enddo
      endif


      write(nout,301)
      write(nout,500) 'total correction:',
     &    average/avborn-1d0,
     &    average/avborn*sqrt((sigma/average)**2+(siborn/avborn)**2)
      if(lnoutmc.ge.2)then
      if (sbarelep.eq.0) then
      write(nout,500)' total correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &    +avgen(0,2)+avgen(0,3)+avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt( (sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &         +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &         +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &         +avgen(0,2)+avgen(0,3)+avgen(0,6)+avgen(0,7)
     &         +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &         +avgen(7,1)+avgen(8,1)
     &         +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)
     &       /avborn
      else
      write(nout,500)' total correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &    +avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt( (sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &         +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &         +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &         +avgen(0,2)+avgen(0,3)
     &         +avgen(0,10)
     &         +avgen(0,6)+avgen(0,7)
     &         +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &         +avgen(7,1)+avgen(8,1)
     &         +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)
     &       /avborn
      endif
      endif
      write(nout,500) 'total QCD correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))/avborn,
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &        +(sigen(0,8))**2+(sigen(0,9))**2
     &        +(sigen(2,1))**2
     &        +(sigen(5,1))**2+(sigen(6,1))**2
     &        +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &        +avgen(2,1)+avgen(5,1)+avgen(6,1))*siborn/avborn)**2)
     &        /avborn
      write(nout,500) 'virt+real QCD correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(2,1))/avborn,
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &        +(sigen(2,1))**2
     &        +((avgen(0,5)+avgen(0,4)+avgen(2,1))*siborn/avborn)**2)
     &        /avborn
      write(nout,500) ' incoming gluon correction:',
     &    (avgen(0,8)+avgen(0,9)+avgen(5,1)+avgen(6,1))/avborn,
     &    sqrt((sigen(0,8))**2+(sigen(0,9))**2
     &        +(sigen(5,1))**2+(sigen(6,1))**2
     &        +((avgen(0,8)+avgen(0,9)
     &          +avgen(5,1)+avgen(6,1))*siborn/avborn)**2)
     &        /avborn
      if (sbarelep.eq.0) then
      write(nout,500) 'total ew correction:',
     &    (avgen(0,2)+avgen(0,3)+avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &        +((avgen(0,2)+avgen(0,3)+avgen(0,6)+avgen(0,7)
     &        +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &        +avgen(7,1)+avgen(8,1)
     &        +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn
      else
      write(nout,500) 'total ew correction:',
     &    (avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &        +((avgen(0,2)+avgen(0,3)
     &        +avgen(0,10)
     &        +avgen(0,6)+avgen(0,7)
     &        +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &        +avgen(7,1)+avgen(8,1)
     &        +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn
      endif
      if (sbarelep.eq.0) then
      write(nout,500) 'virt+real ew correction:',
     &    (avgen(0,2)+avgen(0,3)+avgen(1,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(1,1))**2
     &    +((avgen(0,2)+avgen(0,3)+avgen(1,1))*siborn/avborn)**2)/avborn
      else
      write(nout,500) 'virt+real ew correction:',
     &    (avgen(0,2)+avgen(0,3)+avgen(0,10)+avgen(1,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &        +(sigen(0,10))**2
     &        +(sigen(1,1))**2
     &    +((avgen(0,2)+avgen(0,3)+avgen(0,10)+avgen(1,1))
     &        *siborn/avborn)**2)/avborn
      endif
      write(nout,500) 'incoming photon correction:',
     &    (avgen(0,6)+avgen(0,7)+avgen(3,1)+avgen(4,1)
     &        +avgen(7,1)+avgen(8,1)
     &        +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &     +((avgen(0,6)+avgen(0,7)
     &     +avgen(3,1)+avgen(4,1)
     &     +avgen(7,1)+avgen(8,1)
     &     +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn

#ifdef multpdf
C include varying alpha_s from PDF-set to PDF-set
      do nset=1,pdfmemberto-pdfmemberfrom+1
	mavgen(nset,0,5)=mavgen(nset,0,5)*alphasratio(nset)
	mavgen(nset,0,4)=mavgen(nset,0,4)*alphasratio(nset)
	mavgen(nset,0,8)=mavgen(nset,0,8)*alphasratio(nset)
	mavgen(nset,0,9)=mavgen(nset,0,9)*alphasratio(nset)
	mavgen(nset,2,1)=mavgen(nset,2,1)*alphasratio(nset)
	mavgen(nset,5,1)=mavgen(nset,5,1)*alphasratio(nset)
	mavgen(nset,6,1)=mavgen(nset,6,1)*alphasratio(nset)
	msigen(nset,0,5)=msigen(nset,0,5)*alphasratio(nset)
	msigen(nset,0,4)=msigen(nset,0,4)*alphasratio(nset)
	msigen(nset,0,8)=msigen(nset,0,8)*alphasratio(nset)
	msigen(nset,0,9)=msigen(nset,0,9)*alphasratio(nset)
	msigen(nset,2,1)=msigen(nset,2,1)*alphasratio(nset)
	msigen(nset,5,1)=msigen(nset,5,1)*alphasratio(nset)
	msigen(nset,6,1)=msigen(nset,6,1)*alphasratio(nset)
C	mavgen(nset,0,5)=(mavgen(nset,0,5)+avgen(0,5))*
C     &                    alphasratio(nset)-avgen(0,5)
C	mavgen(nset,0,4)=(mavgen(nset,0,4)+avgen(0,4))*
C     &                    alphasratio(nset)-avgen(0,4)
C	mavgen(nset,0,8)=(mavgen(nset,0,8)+avgen(0,8))*
C     &                    alphasratio(nset)-avgen(0,8)
C	mavgen(nset,0,9)=(mavgen(nset,0,9)+avgen(0,9))*
C     &                    alphasratio(nset)-avgen(0,9)
C	mavgen(nset,2,1)=(mavgen(nset,2,1)+avgen(2,1))*
C     &                    alphasratio(nset)-avgen(2,1)
C	mavgen(nset,5,1)=(mavgen(nset,5,1)+avgen(5,1))*
C     &                    alphasratio(nset)-avgen(5,1)
C	mavgen(nset,6,1)=(mavgen(nset,6,1)+avgen(6,1))*
C     &                    alphasratio(nset)-avgen(6,1)
C	msigen(nset,0,5)=sqrt(msigen(nset,0,5)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)-1D0)**2*sigen(0,5)**2)
C	msigen(nset,0,4)=sqrt(msigen(nset,0,4)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)**2-1D0)**2*sigen(0,4)**2)
C	msigen(nset,0,8)=sqrt(msigen(nset,0,8)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)**2-1D0)**2*sigen(0,8)**2)
C	msigen(nset,0,9)=sqrt(msigen(nset,0,9)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)-1D0)**2*sigen(0,9)**2)
C	msigen(nset,2,1)=sqrt(msigen(nset,2,1)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)-1D0)**2*sigen(2,1)**2)
C	msigen(nset,5,1)=sqrt(msigen(nset,5,1)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)-1D0)**2*sigen(5,1)**2)
C	msigen(nset,6,1)=sqrt(msigen(nset,6,1)**2*alphasratio(nset)**2+
C     &                    (alphasratio(nset)-1D0)**2*sigen(6,1)**2)
      enddo
      if(pdfmemberto.ge.pdfmemberfrom)then
        write(nout,*) ''
        write(nout,*) 'Results for different PDF sets:'
      endif
      do nset=1,pdfmemberto-pdfmemberfrom+1
      write(nout,*) ''
      write(nout,*) 'PDF set:',pdfmemberfrom-1+nset
      write(nout,500)' Born cross section:',
     &    mavborn(nset),msiborn(nset)
C      write(nout,500)' complete cross section:',
C     &    maverage(nset),msigma(nset)
C      if(lnoutmc.ge.2)then
      if (sbarelep.eq.0) then
         write(nout,500)' complete cross section:',
     &        (mavgen(nset,0,1)+mavgen(nset,0,5)+mavgen(nset,0,4)
     &        +mavgen(nset,0,8)+mavgen(nset,0,9)
     &        +mavgen(nset,2,1)+mavgen(nset,5,1)+mavgen(nset,6,1))
     &        +(mavgen(nset,0,2)+mavgen(nset,0,3)
     &        +mavgen(nset,0,6)+mavgen(nset,0,7)
     &        +mavgen(nset,1,1)+mavgen(nset,3,1)+mavgen(nset,4,1)
     &        +mavgen(nset,7,1)+mavgen(nset,8,1)
     &        +mavgen(nset,9,1)+mavgen(nset,10,1)),
     &        sqrt((msigen(nset,0,1))**2
     &         +(msigen(nset,0,5))**2+(msigen(nset,0,4))**2
     &         +(msigen(nset,0,8))**2+(msigen(nset,0,9))**2
     &         +(msigen(nset,2,1))**2
     &         +(msigen(nset,5,1))**2+(msigen(nset,6,1))**2
     &         +(msigen(nset,0,2))**2+(msigen(nset,0,3))**2
     &         +(msigen(nset,0,6))**2+(msigen(nset,0,7))**2
     &         +(msigen(nset,1,1))**2
     &         +(msigen(nset,3,1))**2+(msigen(nset,4,1))**2
     &         +(msigen(nset,7,1))**2+(msigen(nset,8,1))**2
     &         +(msigen(nset,9,1))**2+(msigen(nset,10,1))**2)
      else
         write(nout,500)' complete cross section:',
     &        (mavgen(nset,0,1)+mavgen(nset,0,5)+mavgen(nset,0,4)
     &        +mavgen(nset,0,8)+mavgen(nset,0,9)
     &        +mavgen(nset,2,1)+mavgen(nset,5,1)+mavgen(nset,6,1))
     &        +(mavgen(nset,0,2)+mavgen(nset,0,3)
     &        +mavgen(nset,0,10)
     &        +mavgen(nset,0,6)+mavgen(nset,0,7)
     &        +mavgen(nset,1,1)+mavgen(nset,3,1)+mavgen(nset,4,1)
     &        +mavgen(nset,7,1)+mavgen(nset,8,1)
     &        +mavgen(nset,9,1)+mavgen(nset,10,1)),
     &        sqrt((msigen(nset,0,1))**2
     &        +(msigen(nset,0,5))**2+(msigen(nset,0,4))**2
     &        +(msigen(nset,0,8))**2+(msigen(nset,0,9))**2
     &        +(msigen(nset,2,1))**2
     &        +(msigen(nset,5,1))**2+(msigen(nset,6,1))**2
     &        +(msigen(nset,0,2))**2+(msigen(nset,0,3))**2
     &        +(msigen(nset,0,10))**2
     &        +(msigen(nset,0,6))**2+(msigen(nset,0,7))**2
     &        +(msigen(nset,1,1))**2
     &        +(msigen(nset,3,1))**2+(msigen(nset,4,1))**2
     &        +(msigen(nset,7,1))**2+(msigen(nset,8,1))**2
     &        +(msigen(nset,9,1))**2+(msigen(nset,10,1))**2)
      endif
C      endif
      write(nout,500)' QCD corr. to cross section:',
     &     (mavgen(nset,0,5)+mavgen(nset,0,4)
     &     +mavgen(nset,0,8)+mavgen(nset,0,9)
     &     +mavgen(nset,2,1)+mavgen(nset,5,1)+mavgen(nset,6,1)),
     &     sqrt((msigen(nset,0,5))**2+(msigen(nset,0,4))**2
     &     +(msigen(nset,0,8))**2+(msigen(nset,0,9))**2
     &     +(msigen(nset,2,1))**2
     &     +(msigen(nset,5,1))**2+(msigen(nset,6,1))**2)
      write(nout,500)' incoming gluon cross section:',
     &     (mavgen(nset,0,8)+mavgen(nset,0,9)
     &     +mavgen(nset,5,1)+mavgen(nset,6,1)),
     &     sqrt((msigen(nset,0,8))**2+(msigen(nset,0,9))**2
     &     +(msigen(nset,5,1))**2+(msigen(nset,6,1))**2)
      write(nout,500)' virt+real QCD cross section:',
     &     (mavgen(nset,0,5)+mavgen(nset,0,4)+mavgen(nset,2,1)),
     &     sqrt((msigen(nset,0,5))**2+(msigen(nset,0,4))**2
     &     +(msigen(nset,2,1))**2)
      write(nout,500)' virt. QCD cross section:',
     &     mavgen(nset,0,5),msigen(nset,0,5)
      write(nout,500)' outgoing gluon cross section:',
     &     (mavgen(nset,0,4)+mavgen(nset,2,1)),
     &     sqrt((msigen(nset,0,4))**2+(msigen(nset,2,1))**2)
      if (sbarelep.eq.0) then
         write(nout,500) 'ew corr. to cross section:',
     &        (mavgen(nset,0,2)+mavgen(nset,0,3)
     &        +mavgen(nset,0,6)+mavgen(nset,0,7)
     &        +mavgen(nset,1,1)+mavgen(nset,3,1)+mavgen(nset,4,1)
     &        +mavgen(nset,7,1)+mavgen(nset,8,1)
     &        +mavgen(nset,9,1)+mavgen(nset,10,1)),
     &        sqrt((msigen(nset,0,2))**2+(msigen(nset,0,3))**2
     &        +(msigen(nset,0,6))**2+(msigen(nset,0,7))**2
     &        +(msigen(nset,1,1))**2
     &        +(msigen(nset,3,1))**2+(msigen(nset,4,1))**2
     &        +(msigen(nset,7,1))**2+(msigen(nset,8,1))**2
     &        +(msigen(nset,9,1))**2+(msigen(nset,10,1))**2)
      else
         write(nout,500) 'ew corr. to cross section:',
     &        (mavgen(nset,0,2)+mavgen(nset,0,3)
     &        +mavgen(nset,0,10)
     &        +mavgen(nset,0,6)+mavgen(nset,0,7)
     &        +mavgen(nset,1,1)+mavgen(nset,3,1)+mavgen(nset,4,1)
     &        +mavgen(nset,7,1)+mavgen(nset,8,1)
     &        +mavgen(nset,9,1)+mavgen(nset,10,1)),
     &        sqrt((msigen(nset,0,2))**2+(msigen(nset,0,3))**2
     &        +(msigen(nset,0,10))**2
     &        +(msigen(nset,0,6))**2+(msigen(nset,0,7))**2
     &        +(msigen(nset,1,1))**2
     &        +(msigen(nset,3,1))**2+(msigen(nset,4,1))**2
     &        +(msigen(nset,7,1))**2+(msigen(nset,8,1))**2
     &        +(msigen(nset,9,1))**2+(msigen(nset,10,1))**2)
      endif
      write(nout,500) 'incoming photon cross section:',
     &    (mavgen(nset,0,6)+mavgen(nset,0,7)
     &     +mavgen(nset,3,1)+mavgen(nset,4,1)
     &     +mavgen(nset,7,1)+mavgen(nset,8,1)
     &     +mavgen(nset,9,1)+mavgen(nset,10,1)),
     &    sqrt( (msigen(nset,0,6))**2+(msigen(nset,0,7))**2
     &         +(msigen(nset,3,1))**2+(msigen(nset,4,1))**2
     &         +(msigen(nset,7,1))**2+(msigen(nset,8,1))**2
     &         +(msigen(nset,9,1))**2+(msigen(nset,10,1))**2)
      if (sbarelep.eq.0) then
         write(nout,500) ' virt+real ew cross section:',
     &        (mavgen(nset,0,2)+mavgen(nset,0,3)
     &        +mavgen(nset,1,1)),
     &        sqrt((msigen(nset,0,2))**2+(msigen(nset,0,3))**2
     &        +(msigen(nset,1,1))**2)
      else
         write(nout,500) ' virt+real ew cross section:',
     &        (mavgen(nset,0,2)+mavgen(nset,0,3)
     &        +mavgen(nset,0,10)+mavgen(nset,1,1)),
     &        sqrt((msigen(nset,0,2))**2+(msigen(nset,0,3))**2
     &        +(msigen(nset,0,10))**2+(msigen(nset,1,1))**2)
      endif
      write(nout,500) ' virt. ew cross section:',
     &     mavgen(nset,0,2),msigen(nset,0,2)
      if (sbarelep.eq.0) then
         write(nout,500) 'outgoing photon cross section:',
     &        (mavgen(nset,0,3)+mavgen(nset,1,1)),
     &        sqrt((msigen(nset,0,3))**2+(msigen(nset,1,1))**2)
      else
       write(nout,500) 'outgoing photon cross section:',
     &        (mavgen(nset,0,3)+mavgen(nset,1,1)+mavgen(nset,0,10)),
     &        sqrt((msigen(nset,0,3))**2+(msigen(nset,1,1))**2
     &        +(msigen(nset,0,10))**2)
      endif
      enddo
#endif

#ifdef EXTRAOUT
      if(outputfile2(1:1).eq.' ')then
        nout2=6
      else
        nout2=23
        open(unit=nout2,file=outputfile2,position='append',
     &        status='unknown')
      endif

 2222 format(f12.5,4g18.10,6f14.7)
 2223 format(g12.5,4g18.10,6f14.7)

      if(outputvar.eq.'mh') then
         write(nout2,2222) mcm(25),
c      ' Born cross section:',
     &    avborn,siborn,
c      ' complete cross section:',
     &    average,sigma,
c      ' total correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &    +avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt( (sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2
     &         +(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &         +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &         +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &         +avgen(0,2)+avgen(0,3)
     &         +avgen(0,10)
     &         +avgen(0,6)+avgen(0,7)
     &         +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &         +avgen(7,1)+avgen(8,1)
     &         +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)
     &       /avborn,
c      'total ew correction:',
     &    (avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &        +((avgen(0,2)+avgen(0,3)
     &        +avgen(0,10)
     &        +avgen(0,6)+avgen(0,7)
     &        +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &        +avgen(7,1)+avgen(8,1)
     &        +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn,
c      'total QCD correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))/avborn,
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &        +(sigen(0,8))**2+(sigen(0,9))**2
     &        +(sigen(2,1))**2
     &        +(sigen(5,1))**2+(sigen(6,1))**2
     &        +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &        +avgen(2,1)+avgen(5,1)+avgen(6,1))*siborn/avborn)**2)
     &        /avborn
      else if(outputvar.eq.'qcdfacscalefac2') then
         write(nout2,2222) qcdfacscalefac*qcdfacscalefac,
c      ' Born cross section:',
     &    avborn,siborn,
c      ' complete cross section:',
     &    average,sigma,
c      ' total correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &    +avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt( (sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &         +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &         +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &         +avgen(0,2)+avgen(0,3)
     &         +avgen(0,10)
     &         +avgen(0,6)+avgen(0,7)
     &         +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &         +avgen(7,1)+avgen(8,1)
     &         +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)
     &       /avborn,
c      'total ew correction:',
     &    (avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &        +(sigen(0,10))**2
     &        +(sigen(0,6))**2+(sigen(0,7))**2
     &        +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &        +(sigen(7,1))**2+(sigen(8,1))**2
     &        +(sigen(9,1))**2+(sigen(10,1))**2
     &        +((avgen(0,2)+avgen(0,3)
     &        +avgen(0,10)
     &        +avgen(0,6)+avgen(0,7)
     &        +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &        +avgen(7,1)+avgen(8,1)
     &        +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn,
c      'total QCD correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))/avborn,
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &        +(sigen(0,8))**2+(sigen(0,9))**2
     &        +(sigen(2,1))**2
     &        +(sigen(5,1))**2+(sigen(6,1))**2
     &        +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &        +avgen(2,1)+avgen(5,1)+avgen(6,1))*siborn/avborn)**2)
     &        /avborn
      elseif(outputvar.eq.'deltac') then
         write(nout2,2223) mcslideltac,
c      ' Born cross section:',
     &    avborn,siborn,
c      ' complete cross section:',
     &    average,sigma,
c      ' total correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &    +avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt( (sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &         +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &         +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &         +avgen(0,2)+avgen(0,3)
     &         +avgen(0,10)
     &         +avgen(0,6)+avgen(0,7)
     &         +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &         +avgen(7,1)+avgen(8,1)
     &         +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)
     &       /avborn,
c      'total ew correction:',
     &    (avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &        +(sigen(0,10))**2
     &        +(sigen(0,6))**2+(sigen(0,7))**2
     &        +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &        +(sigen(7,1))**2+(sigen(8,1))**2
     &        +(sigen(9,1))**2+(sigen(10,1))**2
     &        +((avgen(0,2)+avgen(0,3)
     &        +avgen(0,10)
     &        +avgen(0,6)+avgen(0,7)
     &        +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &        +avgen(7,1)+avgen(8,1)
     &        +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn,
c      'total QCD correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))/avborn,
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &        +(sigen(0,8))**2+(sigen(0,9))**2
     &        +(sigen(2,1))**2
     &        +(sigen(5,1))**2+(sigen(6,1))**2
     &        +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &        +avgen(2,1)+avgen(5,1)+avgen(6,1))*siborn/avborn)**2)
     &        /avborn
      elseif(outputvar.eq.'deltas2') then
         write(nout2,2223) mcslideltas*mcslideltas,
c      ' Born cross section:',
     &    avborn,siborn,
c      ' complete cross section:',
     &    average,sigma,
c      ' total correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &    +avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt( (sigen(0,5))**2+(sigen(0,4))**2
     &         +(sigen(0,8))**2+(sigen(0,9))**2
     &         +(sigen(2,1))**2+(sigen(5,1))**2+(sigen(6,1))**2
     &         +(sigen(0,2))**2+(sigen(0,3))**2
     &         +(sigen(0,10))**2
     &         +(sigen(0,6))**2+(sigen(0,7))**2
     &         +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &         +(sigen(7,1))**2+(sigen(8,1))**2
     &         +(sigen(9,1))**2+(sigen(10,1))**2
     &         +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &         +avgen(2,1)+avgen(5,1)+avgen(6,1)
     &         +avgen(0,2)+avgen(0,3)
     &         +avgen(0,10)
     &         +avgen(0,6)+avgen(0,7)
     &         +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &         +avgen(7,1)+avgen(8,1)
     &         +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)
     &       /avborn,
c      'total ew correction:',
     &    (avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
     &    +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))/avborn,
     &    sqrt((sigen(0,2))**2+(sigen(0,3))**2
     &        +(sigen(0,10))**2
     &        +(sigen(0,6))**2+(sigen(0,7))**2
     &        +(sigen(1,1))**2+(sigen(3,1))**2+(sigen(4,1))**2
     &        +(sigen(7,1))**2+(sigen(8,1))**2
     &        +(sigen(9,1))**2+(sigen(10,1))**2
     &    +((avgen(0,2)+avgen(0,3)
     &    +avgen(0,10)
          +avgen(0,6)+avgen(0,7)
     &    +avgen(1,1)+avgen(3,1)+avgen(4,1)
     &    +avgen(7,1)+avgen(8,1)
     &    +avgen(9,1)+avgen(10,1))*siborn/avborn)**2)/avborn,
c      'total QCD correction:',
     &    (avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &    +avgen(2,1)+avgen(5,1)+avgen(6,1))/avborn,
     &    sqrt((sigen(0,5))**2+(sigen(0,4))**2
     &        +(sigen(0,8))**2+(sigen(0,9))**2
     &        +(sigen(2,1))**2
     &        +(sigen(5,1))**2+(sigen(6,1))**2
     &        +((avgen(0,5)+avgen(0,4)+avgen(0,8)+avgen(0,9)
     &        +avgen(2,1)+avgen(5,1)+avgen(6,1))*siborn/avborn)**2)
     &        /avborn
       endif	

      if(nout2.ne.6) close(nout2)
#endif

#ifdef mpiuse
       endif
       If(iMyRank.nE.0) gen=1
#endif

c output histograms
      call settings(pbeam,pevent,weight(1,1),nevt,gen,4)

#ifdef mpiuse
      if(iMyRank.EQ.0) Then
#endif

      if(lnoutmc.ge.3)then
c output exit statistics
 600  format(/' Rejected events')
 700  format(a10,i2,a8,i3,a1,i12,a3,f10.5,a3)
 701  format(a10,i2,a12,i12,a3,f10.5,a3)
 702  format(a10,i2,a12,i12,a3,f10.5,a3/)
      write(nout,600)
      endif
      do gen=0,ngenerator
        if(sgen(gen).ne.0)then
          count=0
          do iexit=0,maxexit
            if(mcexitevt(iexit,gen).ne.0)then
              if(lnoutmc.ge.3)then
              write(nout,700)' generator',gen,
     &            ',   exit',iexit,':',mcexitevt(iexit,gen)
     &            ,' = ',
     &            mcexitevt(iexit,gen)
     &            /dble(nevents)*freqgencon(gen,1)*1d2
     &            ,' % '
              endif
              count=count+mcexitevt(iexit,gen)
            endif
          enddo

          if( (count-mcexitevt(0,gen))
     &        /dble(nevents)*freqgencon(gen,1).gt.1d-2) then
            write(*,*) ' HAWK: run unreliable'
            write(*,*) '       please contact the authors'//
     &          ' providing input and output'
            stop
          elseif( (count-mcexitevt(0,gen))
     &        /dble(nevents)*freqgencon(gen,1).gt.1d-3) then
            write(*,*) ' HAWK: run potentially unreliable'
            write(*,*) '       use results with caution'
          endif

          if(lnoutmc.ge.3)then
          write(nout,701)' generator',gen,
     &        ', all exits:',count-mcexitevt(0,gen)
     &            ,' = ',
     &            (count-mcexitevt(0,gen))
     &        /dble(nevents)*freqgencon(gen,1)*1d2
     &            ,' % '
          write(nout,702)' generator',gen,
     &        ', sum exits:',count
     &            ,' = ',
     &            count/dble(nevents)*freqgencon(gen,1)*1d2
     &            ,' % '
          endif
        endif
      enddo

      if(lnoutmc.ge.3)then
c output cut statistics
 800  format(/' Cutted events')
 900  format(a10,i2,a7,i2,a7,i3,a1,i12,a3,f10.6,a3)
 901  format(a10,i2,a7,i2,a11,i12,a3,f10.6,a3)
 902  format(a21,a7,i3,a1,i12,a3,f10.6,a3)
 903  format(a21,a11,i12,a3,f10.6,a3)
 904  format(a10,i2,a7,i2,a11,i12,a3,f10.6,a3/)
 905  format(a21,a11,i12,a3,f10.6,a3/)
      write(nout,800)
      count=0
      if(swborn(0).ne.0)then
        do icut=1,maxcut
          if(mccutevt_0(icut).ne.0)then
            write(nout,902)' Born generator',
     &          ',   cut',icut,':',mccutevt_0(icut)
     &      ,' = '
     &      ,mccutevt_0(icut)/dble(nevents)*freqgencon(0,1)*1d2
     &      ,' % '
            count=count+mccutevt_0(icut)
          endif
        enddo
        write(nout,903)' Born generator',
     &      ', all cuts:',count
     &      ,' = ',
     &      count/dble(nevents)*freqgencon(0,1)*1d2
     &      ,' % '
        write(nout,903)' Born generator',
     &      ',  no cut :',mccutevt_0(0)
     &      ,' = ',
     &      mccutevt_0(0)
     &      /dble(nevents)*freqgencon(0,1)*1d2
     &      ,' % '
        write(nout,905)' Born generator',
     &      ', sum cuts:',count+mccutevt_0(0)
     &      ,' = ',
     &      (count+mccutevt_0(0))/dble(nevents)*freqgencon(0,1)*1d2
     &      ,' % '
      endif
      do gen=0,ngenerator
        if(sgen(gen).ne.0)then
          do icon=1,ncontr(gen)
            count=0
            if(sgencon(gen,icon).ne.0)then
              if(gen.eq.0)then
                do icut=1,maxcut
                  if(mccutevt(icut,gen,icon).ne.0)then
                    write(nout,900)' generator',gen,', con',icon,
     &                  ',   cut',icut,':',mccutevt(icut,gen,icon)
     &            ,' = ',
     &            mccutevt(icut,gen,icon)
     &                  /dble(nevents)*freqgencon(gen,icon)*1d2
     &            ,' % '
                    count=count+mccutevt(icut,gen,icon)
                  endif
                enddo
              else
                do icut=-2*nexternal(gen),maxcut
                  if(mccutevt(icut,gen,icon).ne.0.and.icut.ne.0)then
                    write(nout,900)' generator',gen,', con',icon,
     &                  ',   cut',icut,':',mccutevt(icut,gen,icon)
     &            ,' = ',
     &            mccutevt(icut,gen,icon)
     &                  /dble(nevents)*freqgencon(gen,1)*1d2
     &            ,' % '
                    count=count+mccutevt(icut,gen,icon)
                  endif
                enddo
              endif
              write(nout,901)' generator',gen,', con',icon,
     &            ', all cuts:',count
     &            ,' = ',
     &            (count)/dble(nevents)*freqgencon(gen,icon)*1d2
     &            ,' % '
               write(nout,901)' generator',gen,', con',icon,
     &            ',  no cut :',mccutevt(0,gen,icon)
     &            ,' = ',
     &            mccutevt(0,gen,icon)
     &            /dble(nevents)*freqgencon(gen,icon)*1d2
     &            ,' % '
             write(nout,904)' generator',gen,', con',icon,
     &            ', sum cuts:',count+mccutevt(0,gen,icon)
     &            ,' = ',
     &            (count+mccutevt(0,gen,icon))
     &             /dble(nevents)*freqgencon(gen,icon)*1d2
     &            ,' % '
            endif
          enddo
        endif
      enddo
      endif

 97   format(i2.2,'.',i2.2,'.20',i2.2)
 98   format(i2.2,'.',i2.2,'.',i2.2)
      call date_and_time(datum)
      da1=datum(7:8)
      da2=datum(5:6)
      da3=datum(1:4)
      datum=da1//'.'//da2//'.'//da3
      call itime(iarray)
      write(chtime,98)iarray(1),iarray(2),iarray(3)
      call cpu_time(rcpu)

      icpuh=int(rcpu/3600.)
      icpum=int(rcpu/60.-icpuh*60)
      icpus=int(rcpu-icpuh*3600-icpum*60)
      icpu =icpuh*60+icpum
      write(nout,100) datum,chtime,icpuh,icpum,icpus
100   format(/'  date: ',a10,' ,    time: ',a8,' hour ,',
     &       '    cpu: ',i3,' h',i3,' m',i3,' s')

#ifdef mpiuse
       endif
#endif

      if(outputfile(1:2).ne.'  ')close(nout)	

      if(lprintweightmax) then
        close(20)
      endif

      end

***********************************************************************
      subroutine integration()
***********************************************************************
*     phase space integration                                         *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     21.02.06  Ansgar Denner     last changed  14.11.12              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

#ifdef mpiuse
      include "mpif.h"
      integer status(MPI_STATUS_SIZE)
      double precision g(maxch,maxps)
      real*8 wi_ext(maxch,0:maxg),wran_ext(maxgr,maxr,0:maxg),dummy
      double precision keep_av,keep_si
      integer do_optimize,proccount
      common/do_optimize_common/do_optimize
#endif

#ifdef multpdf
#include "multpdf.h"
      integer nset
#endif
      integer iNumProcs, iErr, iMyRank
      common/mpistuff/iNumProcs,iMyRank

c local variables
      real*8 srandom(3*maxe+1)
      real*8 random(3*maxe+1),weighttot,sigmatot,weightgen(maxcon)
#ifdef SUBPROCESSES
      real*8 weighttotii(-4:4,-4:4),weightgenii(maxcon,-4:4,-4:4)
#endif
      integer switch,channel,phasespace
      integer i1,i2

      real*8 sevent(maxe,maxe),kevent(maxe,0:3)

      logical misnan



c variables for testing
      real*8 weighttest(maxps,maxcon),weightgentest(maxps)
      character*5 part(-4:4)

      integer selprocess
      common/selprocess/selprocess




 1050           format(1x,a28,": phasespace,weight = ",i4,g20.10)
 1051           format(1x,a28,":            weight = ",4x,g20.10)
 1029           format(/1x,a25,13x,3a24,a29)
 1030           format(/1x,a50/1x,a25,13x,3a24,a29)
 1031           format(1x,a5,a5,a21,' = ',5g24.14)
 1032           format(1x,a31,' = ',5g24.14)

* initialize hadrons


* initialize parton distributions

#ifdef mpiuse
      do_optimize=0
      proccount=0
#endif

c numerical integration
      do evt=1,nevents
      nevt=evt

c generate random numbers

c random1 choice of channel
c random3 choice of proton direction
c random4 optimization of x1
c random5 optimization of x2
c random6 x1  parton energy fraction
c random7 x2  parton energy fraction
c random8 z   energy fraction after collinear photon emission
c (random9 (1-z) for exclusive final states)
c random10 -- (3*maxe+1)  (3*maxe-8) phase space variable for maxe-2
c                    final state particles [(3maxe - 10) needed]?


c random1 choice of channel
c random2 choice of proton direction
c random5 optimization of x1
c random6 optimization of x2
c random7 x1
c random8 x2
c random9 z
c (random10 (1-z) for exclusive final states)
c random 10- (3*maxe+1)  phase space variable for maxe-1
c                       final state particles


#ifdef RANLUX
        call ranlux(random(1),3*maxe+1)
#else
        do rnr=1,3*maxe+1
          call ransrcarry(random(rnr))
        enddo
#endif

#ifdef mpiuse
        proccount=proccount+1
        if(proccount.EQ.iNumProcs) proccount=0
        if(iNumProcs.EQ.1) proccount=0

        if(proccount.EQ.iMyRank) Then
#endif


        weighttot=0d0
        sigmatot=0d0
#ifdef multpdf
	do nset=1,pdfmemberto-pdfmemberfrom+1
          mweighttot(nset)=0D0
          msigmatot(nset)=0D0
	enddo
#endif
#ifdef SUBPROCESSES
        do i1=-4,4
          do i2=-4,4
            weighttotii(i1,i2)=0d0
          enddo
        enddo
#endif

* partonic energy fractions
* put here???
ccc        print *, '6. ngenerator =  ',ngenerator

        do gen=0,ngenerator
ccc           print *, 'ngen = ',ngenerator, '   gen = ', gen,
ccc     &          '   sgen(',gen,') = ', sgen(gen),
ccc     &          '   scuts(',gen,') = ', scuts(gen),
ccc     &          '   ptcutl = ', ptcutl
          cogen=gen
ccc         print *, 'sgen(', gen, ') =  ', sgen(gen)
         if(sgen(gen).ne.0)then

* if switch gets > 0 event does not contribute owing to numerical problems.
            switch=0

* choice of channel


            do cha=1,nchannel(gen)
              if(random(1).ge.betaopt(cha-1,gen).and.
     &          random(1).le.betaopt(cha,gen))then
                channel=cha
              endif
            enddo


            cocha=channel




      if (selprocess.eq.0)then
         call vbfh_integrand(random(2),evt,channel,gen,switch)
      elseif (selprocess.ge.1) then
ccc         print *,'gen',gen
         call whzh_integrand(random(2),evt,channel,gen,switch)
      endif


c discard unreasonable events
            do phasespace=1,nphsp(gen)
              do icon=1,ncontr(gen)
                if(switch.eq.0.and.
     &              rejectedcon(icon,phasespace).eq.0.and.
     &              misnan(weight(icon,phasespace)))then
                  switch=50
                  if(numout.lt.maxout.and.lnoutmc.ge.5)then
 105                format(/" event = ",i9,", gen = ",i4,
     &                  ", channel = ",i4,
     &                  ", exit = ",i4,(a))
                    write(nout,105) nevt,cogen,cocha,switch,
     &                  ', integration: NaN weight'
                    numout=numout+1
                  endif
                endif
              enddo
            enddo
c--<

            if(switch.eq.0)then
              if(swborn(gen).eq.1)then
                if(rejected_0.eq.0)then
                  avborn=avborn+weight_0
                  siborn=siborn+weight_0**2


#ifdef multpdf
                  do nset=1,pdfmemberto-pdfmemberfrom+1
                    mavborn(nset)=mavborn(nset)+mweight_0(nset)
                    msiborn(nset)=msiborn(nset)+mweight_0(nset)**2
		  enddo
#endif
#ifdef SUBPROCESSES
                  do i1=-4,4
                    do i2=-4,4
                      avbornii(i1,i2)=avbornii(i1,i2)
     &                    +weightii_0(i1,i2)
                      sibornii(i1,i2)=sibornii(i1,i2)
     &                    +weightii_0(i1,i2)**2
                    enddo
                  enddo
#endif

                else
                  weight_0=0d0
#ifdef multpdf
      do nset=1,pdfmemberto-pdfmemberfrom+1
                  mweight_0(nset)=0d0
      enddo
#endif
                endif
                mccutevt_0(rejected_0)=mccutevt_0(rejected_0)+1
              endif
              do icon=1,ncontr(gen)
                if(mod(nevt,iNumProcs*freqgencon(gen,icon)).lt.
     &             iNumProcs)
     &              mccutevt(rejected(icon),gen,icon)=
     &                  mccutevt(rejected(icon),gen,icon)+1
              enddo


              do icon=1,ncontr(gen)
                weightgen(icon)=0d0
#ifdef multpdf
		do nset=1,pdfmemberto-pdfmemberfrom+1
                  mweightgen(nset,icon)=0d0
		enddo
#endif
#ifdef SUBPROCESSES
                do i1=-4,4
                  do i2=-4,4
                    weightgenii(icon,i1,i2)=0d0
                  enddo
                enddo
#endif
              enddo

              do phasespace=1,nphsp(gen)
                do icon=1,ncontr(gen)

                  if(rejectedcon(icon,phasespace).eq.0)then
                    weightgen(icon)=weightgen(icon)
     &                  +weight(icon,phasespace)
#ifdef multpdf
		    do nset=1,pdfmemberto-pdfmemberfrom+1
                      mweightgen(nset,icon)=mweightgen(nset,icon)
     &                    +mweight(nset,icon,phasespace)
		    enddo
#endif

#ifdef SUBPROCESSES
                    do i1=-4,4
                    do i2=-4,4
                      weightgenii(icon,i1,i2)=weightgenii(icon,i1,i2)
     &                    +weightii(icon,i1,i2,phasespace)
                    enddo
                    enddo
#endif
                  endif
                  mccutsubevt(rejectedcon(icon,phasespace),
     &                                          gen,icon,phasespace)=
     &                mccutsubevt(rejectedcon(icon,phasespace),
     &                                         gen,icon,phasespace)+1
               enddo            ! icon



*     fill histograms with subevents
                call settings(pbeam,pevent(1,0,phasespace),
     &              weight(1,phasespace),evt,gen,2)


              enddo             ! phasespace






c calculation of total weight

              do icon=1,ncontr(gen)
                weighttot=weighttot+weightgen(icon)
                sigmatot=sigmatot+weightgen(icon)**2
#ifdef multpdf
		do nset=1,pdfmemberto-pdfmemberfrom+1
                mweighttot(nset)=mweighttot(nset)+mweightgen(nset,icon)
                msigmatot(nset)=msigmatot(nset)+mweightgen(nset,icon)**2
		enddo
#endif
#ifdef SUBPROCESSES
                do i1=-4,4
                  do i2=-4,4
                    weighttotii(i1,i2)=weighttotii(i1,i2)
     &                  +weightgenii(icon,i1,i2)
                    avgenii(gen,icon,i1,i2)=avgenii(gen,icon,i1,i2)
     &                  +weightgenii(icon,i1,i2)
                    sigenii(gen,icon,i1,i2)=sigenii(gen,icon,i1,i2)
     &                  +weightgenii(icon,i1,i2)**2
                  enddo
                enddo
#endif
                avgencha(gen,icon,channel)=avgencha(gen,icon,channel)
     &              +weightgen(icon)
                sigencha(gen,icon,channel)=sigencha(gen,icon,channel)
     &              +weightgen(icon)**2

#ifdef multpdf
                do nset=1,pdfmemberto-pdfmemberfrom+1
                mavgencha(nset,gen,icon,channel)=
     &               mavgencha(nset,gen,icon,channel)
     &              +mweightgen(nset,icon)
                msigencha(nset,gen,icon,channel)=
     &               msigencha(nset,gen,icon,channel)
     &              +mweightgen(nset,icon)**2
	        enddo
#endif
                nevtgencha(gen,icon,channel)=nevtgencha(gen,icon,channel
     &              )+ 1
c rejected weight
                if(weightgen(icon).eq.0d0.and.
     &              mod(nevt,iNumProcs*freqgencon(gen,icon)).lt.
     &                  iNumProcs)then
                  rejcha(gen,icon,channel)=rejcha(gen,icon,channel)+1
                endif


c     for testing print event
c end testing print event

c maximal weight of generator
                if(dabs(weightgen(icon)).gt.
     &              weightmaxcha(gen,icon,channel))then
                  weightmaxcha(gen,icon,channel)=dabs(weightgen(icon))
                  nmaxcha(gen,icon,channel)=evt
#ifdef PRINTWEIGHTMAX
                  if(lprintweightmax) then
                    if(weightmaxcha(gen,icon,channel).gt.
     &                  weightmax(gen,icon)) then
                      weightmax(gen,icon)=weightmaxcha(gen,icon,channel)
                      channelmax(gen,icon) = channel
                      nmax(gen,icon) = evt

                      write(20,*) ' new maximal weight'
                      write(20,*) '      generator = ',gen
                      write(20,*) '   contribution = ',icon
                      write(20,*) '      channel   = ',channel
                      write(20,*) '      event     = ',evt
                      write(20,*) '      weight    = ',
     &                    weightmaxcha(gen,icon,channel)
                      write(20,*) '      gsum      = ',gsumm(1)
                      write(20,*) '      factor    = ',factor(gen)
                      write(20,*) '      wconv     = ',wconv
                      write(20,*) '      m2pdfpdf  = ',
     &                  weightmaxcha(gen,icon,channel)/wconv/factor(gen)
     &                    *gsumm(1)
     &                    , weightmaxcha(gen,icon,channel)/wconv
     &                    /factor(gen)
                      write(20,*) '            x1  = ',x1print
                      write(20,*) '            x2  = ',x2print
                      if (gen.eq.0)then
                        write(20,*) '            x   = ',xprint,omxprint
                      endif
                      do ext=1,nextgen(gen)
                        do vcomp=0,3
                          write(20,100) ext,vcomp,pprint(ext,vcomp)
                        enddo
                      enddo
                      call dets(pprint,sevent,nextgen(gen),0)
                      do ext=1,nextgen(gen)
                        do ext2=ext,nextgen(gen)
                          write(20,110) ext,ext2,sevent(ext,ext2)
                        enddo
                      enddo
 101                  format(1x,a18,d24.18)
 100                  format(6x,' p(',i1,',',i1,')    = ',d24.18)
 110                  format(6x,' s(',i1,',',i1,')    = ',d24.18)
                      write(20,*)
                    endif
                  endif
#endif
                endif


c     end loop over icons
              enddo



c test of convolution contribution to subtraction

c     end test of convolution contribution to subtraction

c     test of convolution contribution to slicing

c     end test of convolution contribution to slicing

c test of incoming photon convolution contribution to subtraction

c test of incoming photon convolution contribution to slicing

c test of incoming photon convolution contribution to eff. coll. fact
c     end test of incoming photon convolution contributions

c test of incoming gluon convolution contribution to subtraction

c     end test of incoming gluon convolution contributions

c test of incoming gluon convolution contribution with eff. coll. factors
c     end test of incoming gluon convolution contributions

c     end output of collinear events for test purposes

              if(mod(nevt,iNumProcs*freqgencon(gen,1)).lt.iNumProcs)
     &            mcexitevt(switch,gen)=mcexitevt(switch,gen)+1
            else
              mcexitevt(switch,gen)=mcexitevt(switch,gen)+1
            endif               ! switch.eq.0
          else
            weight_0=0d0
#ifdef multpdf
      do nset=1,pdfmemberto-pdfmemberfrom+1
            mweight_0(nset)=0d0
      enddo
#endif
          endif                 ! sgen(gen).ne.0
c     maximal weight
          if(dabs(weighttot).gt.weighttotmax)then
            weighttotmax=dabs(weighttot)
            nmaxtot=evt
          endif

c filling histograms with full event
          weight(1,1)=weight_0
#ifdef multpdf
      do nset=1,pdfmemberto-pdfmemberfrom+1
      enddo
#endif
          call settings(pbeam,pborn,weight(1,1),evt,gen,3)

        enddo
c end loop over generators

c calculating results
c total results
        average=average+weighttot
c weighttot = sum_icon weightgen(icon)
        sigma=sigma+weighttot**2
c sigmatot = sum_icon weightgen(icon)**2
        sigmaa=sigmaa+sigmatot
#ifdef multpdf
        do nset=1,pdfmemberto-pdfmemberfrom+1
	  maverage(nset)=maverage(nset)+mweighttot(nset)
          msigma(nset)=msigma(nset)+mweighttot(nset)**2
          msigmaa(nset)=msigmaa(nset)+msigmatot(nset)
	enddo
#endif
#ifdef SUBPROCESSES
        do i1=-4,4
        do i2=-4,4
          averageii(i1,i2)=averageii(i1,i2)+weighttotii(i1,i2)
          sigmaii(i1,i2)=sigmaii(i1,i2)+weighttotii(i1,i2)**2
        enddo
        enddo
#endif
c partial results
        av=av+weighttot
        si=si+weighttot**2

#ifdef mpiuse
        endif

        do gen=0,ngenerator	
	  if(nevt.eq.nopt(opt(gen),gen).AND.iNumProcs.NE.1)then
	     if(iMyRank.NE.0) then
	       call MPI_Send(wi, maxch*(maxg+1),
     &              MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, iErr)
               do i1=1,nchannel(gen)
		  wi(i1,gen)=0D0
	       enddo	
	     endif

	     if (iMyRank.EQ.0) then
	       do i2=1,iNumProcs-1
		 call MPI_Recv(wi_ext, maxch*(maxg+1),
     &                MPI_DOUBLE_PRECISION, i2, 1,
     &                MPI_COMM_WORLD, status, iErr)
        	 do i1=1,nchannel(gen)
		   wi(i1,gen)=wi(i1,gen)+wi_ext(i1,gen)
		 enddo
	       enddo

               do_optimize=1
	       call optimization(0D0,g,0D0,nevt,gen,nchannel(gen))
               do_optimize=0

	       do i2=1,iNumProcs-1
		 call MPI_Send(betaopt, (maxch+1)*(maxg+1),
     &                MPI_DOUBLE_PRECISION, i2, 2, MPI_COMM_WORLD, iErr)
		 call MPI_Send(alphaopt, maxch*maxo*(maxg+1),
     &                MPI_DOUBLE_PRECISION, i2, 3, MPI_COMM_WORLD, iErr)
	       enddo	
	     endif

             if(iMyRank.NE.0) then
	       call MPI_Recv(betaopt, (maxch+1)*(maxg+1),
     &              MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD,
     &              status, iErr)
	       call MPI_Recv(alphaopt, maxch*maxo*(maxg+1),
     &              MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD,
     &              status, iErr)
               if(nevents.gt.2*nevt) opt(gen)=opt(gen)+1
	     endif	 	 	
	  endif		

	  if(noptran(gen).NE.0) then
            if(mod(nevt,noptran(gen)).eq.0.AND.iNumProcs.NE.1)then
	      if(iMyRank.NE.0) then
	        call MPI_Send(wran, maxgr*maxr*(maxg+1),
     &              MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, iErr)
                do igrv=1,maxr
        	  do igrp=1,maxgr
        	    wran(igrp,igrv,gen)=0d0
        	  enddo
                enddo
	      endif

	      if (iMyRank.EQ.0) then
	        do i2=1,iNumProcs-1
	 	  call MPI_Recv(wran_ext, maxgr*maxr*(maxg+1),
     &                MPI_DOUBLE_PRECISION, i2, 1, MPI_COMM_WORLD,
     &                status, iErr)
        	  do igrv=1,maxr
        	    do igrp=1,maxgr
        	      wran(igrp,igrv,gen)=wran(igrp,igrv,gen)+
     &                                   wran_ext(igrp,igrv,gen)
        	    enddo
        	  enddo
	        enddo

                do_optimize=1
                call ranopt(0.5D0,0.5D0,dummy,dummy,nevt,1,gen,switch)
                call ranopt(0.5D0,0.5D0,dummy,dummy,nevt,2,gen,switch)
                do_optimize=0

	        do i2=1,iNumProcs-1
		  call MPI_Send(beran, (maxgr+1)*maxr*(maxg+1),
     &                MPI_DOUBLE_PRECISION, i2, 2, MPI_COMM_WORLD, iErr)
		  call MPI_Send(alran, maxgr*maxr*(maxg+1),
     &                MPI_DOUBLE_PRECISION, i2, 3, MPI_COMM_WORLD, iErr)
	        enddo	
	      endif

              if(iMyRank.NE.0) then
	        call MPI_Recv(beran, (maxgr+1)*maxr*(maxg+1),
     &              MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD,
     &              status, iErr)
	        call MPI_Recv(alran, maxgr*maxr*(maxg+1),
     &              MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD,
     &              status, iErr)
	      endif	 	 	
	    endif	
	  endif	
	enddo
#endif

c preliminary results
        if(mod(evt,nsteps).eq.0)then
#ifdef mpiuse
	  if(iMyRank.EQ.0.AND.iNumProcs.NE.1) Then
	    do i2=1,iNumProcs-1
	      call MPI_Recv(keep_av, 1, MPI_DOUBLE_PRECISION,
     &                      i2, 10, MPI_COMM_WORLD, status, iErr)
	      call MPI_Recv(keep_si, 1, MPI_DOUBLE_PRECISION,
     &                      i2, 11, MPI_COMM_WORLD, status, iErr)
              av=av+keep_av
              si=si+keep_si
	      average=average+keep_av
              sigma=sigma+keep_si
            enddo
	  else if(iMyRank.NE.0.AND.iNumProcs.NE.1) Then
	    call MPI_Send(av, 1, MPI_DOUBLE_PRECISION,
     &                    0, 10, MPI_COMM_WORLD, iErr)
	    call MPI_Send(si, 1, MPI_DOUBLE_PRECISION,
     &                    0, 11, MPI_COMM_WORLD, iErr)
            av=0D0
            si=0D0
	  endif
        if(iMyRank.EQ.0) Then
#endif
          if(evt.eq.nsteps) then
            if (lnoutmc.ge.3) write(nout,*)' '
            write(nout,'(a)')'    events  :     intermediate results'
     &          //'     :     preliminary results '
          endif
          write(nout,'(i10,2("  : ",1f12.5," +-",1f12.5))')evt,
     &        av/dble(nsteps),
     &        sqrt(dabs(si-av**2/dble(nsteps)))/dble(nsteps),
     &        average/dble(evt),
     &        sqrt(dabs(sigma-average**2/dble(evt)))/dble(evt)
          call flush(nout)
          av=0d0
          si=0d0
#ifdef mpiuse
        endif
#endif
        endif


 999    nevt = evt
      enddo
c end loop over events

      do gen=0,ngenerator
         do icon=1,ncontr(gen)
            do cha=1,nchannel(gen)
               avgen(gen,icon)=avgen(gen,icon)+avgencha(gen,icon,cha)
               sigen(gen,icon)=sigen(gen,icon)+sigencha(gen,icon,cha)
#ifdef multpdf
               do nset=1,pdfmemberto-pdfmemberfrom+1
                  mavgen(nset,gen,icon) = mavgen(nset,gen,icon) +
     &                 mavgencha(nset,gen,icon,cha)
                  msigen(nset,gen,icon) = msigen(nset,gen,icon) +
     &                 msigencha(nset,gen,icon,cha)
               enddo
#endif
               nevtgen(gen,icon) = nevtgen(gen,icon)
     &              + nevtgencha(gen,icon,cha)
               rejgen(gen,icon) = rejgen(gen,icon)
     &              + rejcha(gen,icon,cha)
               if (weightmax(gen,icon).lt.weightmaxcha(gen,icon,cha))
     &              then
                  weightmax(gen,icon) = weightmaxcha(gen,icon,cha)
                  channelmax(gen,icon) = cha
                  nmax(gen,icon) = nmaxcha(gen,icon,cha)
               endif
            enddo
         enddo
      enddo


      end
***********************************************************************


***********************************************************************
      subroutine dets(k,sk,lnext,switch)
***********************************************************************
*     calculation of invariants  sij = 2pi*pj!!!                      *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     21.02.06  Ansgar Denner     last changed  16.01.07              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 sk(maxe,maxe),k(maxe,0:3)
      integer lnext,switch
      if(switch.gt.0)return
      do ext=1,lnext
      do ext2=ext+1,lnext
        sk(ext,ext2)=2d0*(k(ext,0)*k(ext2,0)-k(ext,1)*k(ext2,1)
     &    -k(ext,2)*k(ext2,2)-k(ext,3)*k(ext2,3))
        sk(ext2,ext)=sk(ext,ext2)
      enddo
      enddo
      end

***********************************************************************
      subroutine ranopt(rand,random,ranout,
     &    lweight,ln,nrand,generator,switch)
***********************************************************************
*     adaptive optimization                                           *
*     adapted from COFFERGG by Markus Roth                            *
*     divides interval from 0 to 1 into maxgr equal intervals         *
*     the interval is chosen according to its contribution using      *
*     the random number rand and weights alpha that are determined    *
*     after each noptran events                                       *
*     the output random number ranout is determined from random       *
*     that runs throug the interval determined by rand                *
*---------------------------------------------------------------------*
*     by Markus Roth         last changed   29.03.07   Ansgar Denner  *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

#ifdef mpiuse
      integer do_optimize
      common/do_optimize_common/do_optimize
#endif
      integer iNumProcs,iMyRank
      common/mpistuff/iNumProcs,iMyRank

c local variables
      real*8 lweight,rand,random,ranout
      integer ngrid,nrand,generator,ln,switch

c to avoid compiler warning
      data ngrid /0/

      if(switch.gt.0)return
c importance sampling
        if(rand.eq.beran(0,nrand,generator))then
          write(*,*) 'ranopt rand=beran',rand,
     &        beran(0,nrand,generator),0,nrand,generator
        endif
      do igrp=1,maxgr
        if(rand.ge.beran(igrp-1,nrand,generator).and.
     &    rand.le.beran(igrp,nrand,generator))then
          ngrid=igrp
        endif
      enddo
      ranout=(dble(ngrid-1)+random)/dble(maxgr)
c local density
      lweight=lweight/(dble(maxgr)*alran(ngrid,nrand,generator))
c optimization
#ifdef mpiuse
        if((iMyRank.EQ.0 .AND. do_optimize.EQ.1) .OR.
     &      iNumProcs.EQ.1) Then
#endif
      if(mod(ln,noptran(generator)).eq.0)then
        do igrp=1,maxgr
          alran(igrp,nrand,generator)=dabs(wran(igrp,nrand,generator))
          beran(igrp,nrand,generator)=beran(igrp-1,nrand,generator)
     &      +alran(igrp,nrand,generator)
        enddo
        do igrp=1,maxgr
          if(beran(maxgr,nrand,generator).ne.0d0)then
            alran(igrp,nrand,generator)=alminran/dble(maxgr)
     &        +(1d0-alminran)*alran(igrp,nrand,generator)
     &          /beran(maxgr,nrand,generator)
            beran(igrp,nrand,generator)=beran(igrp-1,nrand,generator)
     &        +alran(igrp,nrand,generator)
          else
            alran(igrp,nrand,generator)=1d0/dble(maxgr)
            beran(igrp,nrand,generator)=dble(igrp)/dble(maxgr)
          endif
        enddo
      endif
#ifdef mpiuse
        endif
#endif
      end

***********************************************************************
      subroutine rancum(random,lweight,nrand,generator,switch)
***********************************************************************
*     cumulate weights for adaptive optimization                      *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     30.01.06  Ansgar Denner     last changed  08.12.06              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 lweight,random
      integer nrand,ngrid,generator,switch

      if(switch.gt.0)return

c cumulate weight for optimization
      if(random.lt.1d0)then
        ngrid=int(random*dble(maxgr))+1
      else
        ngrid=maxgr
        write(*,*) 'rancum random = 1',random,
     &   lweight,nrand,generator,switch
      endif

      wran(ngrid,nrand,generator)=wran(ngrid,nrand,generator)
     &    +lweight
      end

***********************************************************************
      function dli2(x,switch)
***********************************************************************
*     dilogarithmus with real argument less than 1                    *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     27.09.06  Ansgar Denner     last changed  27.09.06              *
***********************************************************************
      implicit none
c local variables
      real*8 x,dli2,dli2gen,pi2_6
      integer switch

      data pi2_6 /1.6449340668482264364724151666d0/

c mcoutput
      integer nout,numout,maxout,noutgen,nsteps,lnoutmc
      common/mcoutput/nout,numout,maxout,noutgen,nsteps,lnoutmc

      if(x.ge.-1d0.and.x.le.0.5d0)then
        dli2=dli2gen(x,switch)
      elseif(x.gt.0.5d0.and.x.lt.1d0)then
        dli2=pi2_6-dli2gen(1d0-x,switch)-dlog(x)*dlog(1d0-x)
      elseif(x.lt.-1d0)then
        dli2=-pi2_6-dli2gen(1d0/x,switch)-0.5d0*log(-x)**2
      elseif(x.eq.1d0)then
        dli2=pi2_6
      else
        write(*,'(a24,d16.10)')' dli2: no case found: x=',x
        stop
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     dilogarithm                                                 c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function dli2gen(x,switch)
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      integer n,nb
      parameter(nb=22)
      real*8 b(0:nb),error,x,dli2gen,x1,x2,xn,new,dli2new
      integer switch
      save b

c error is lower limit of the uncertainty, no message with error=0d0
      data b/1d0,-0.25d0,
     &       2.7777777777777777777777777777777776d-2,0d0,
     &       -2.7777777777777777777777777777777778d-4,0d0,
     &       4.7241118669690098261526832955404383d-6,0d0,
     &       -9.1857730746619635508524397413286296d-8,0d0,
     &       1.8978869988970999072009173019274030d-9,0d0,
     &       -4.0647616451442255268059093862919667d-11,0d0,
     &       8.9216910204564525552179873167527489d-13,0d0,
     &       -1.9939295860721075687236443477937898d-14,0d0,
     &       4.5189800296199181916504765528555929d-16,0d0,
     &       -1.0356517612181247014483411542218656d-17,0d0,
     &       2.3952186210261867457402837430009803d-19/
      x1=-dlog(1d0-x)
      x2=x1*x1
      xn=x2*x1
      dli2gen=b(0)*x1+b(1)*x2+b(2)*xn
      do n=4,nb,2
        xn=xn*x2
        new=b(n)*xn
        dli2new=dli2gen+new
        if(dli2new.eq.dli2gen)return
        dli2gen=dli2new
      enddo
      if(dabs(new/dli2gen).gt.error.and.error.ne.0d0)then
        error=dabs(new/dli2gen)
        switch=33
        if(numout.lt.maxout.and.lnoutmc.ge.5)then
 100        format(/" event = ",i9,", gen = ",i4,", channel = ",i4,
     &        ", exit = ",i4,(a))
            write(nout,100) nevt,cocha,switch
          write(nout,'(a25,d16.10)')' dli2gen: converges badly at ',x
          write(nout,'(a)')'         with error ',error
          numout=numout+1
        endif
        return
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                 c
c     gamma function                                              c
c                                                                 c
c     written by Markus Roth                                      c
c                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function gamma(x)
      implicit none
c local variable
      integer maxa
      parameter(maxa=10)
      real*8 gamma,x,omx,a(maxa),prod,pi
      integer i1

      pi=4d0*datan(1d0)
      data a/0.5772156649015328606d0,0d0,
     &       0.4006856343865314285d0,0d0,
     &       0.2073855510286739853d0,0d0,
     &       0.1440498967688461181d0,0d0,
     &       0.1113342658695646905d0,0d0/
c good convergence for small values of 1-x
      omx=x-1d0
      gamma=dexp(0.5d0*dlog(omx*pi/dsin(pi*omx)))
      do i1=1,maxa,2
        prod=dexp(-(a(i1)*omx**i1))
        gamma=gamma*prod
        if(dabs(prod-1d0).lt.1d-10)return
      enddo
      write(*,'(a)')' gamma: converges badly'
      stop
      end

***********************************************************************
      subroutine histogram(filename,y,x,xminh,xmaxh,nbin,histo,levt,
     *  addpoint,lgen,step)
***********************************************************************
*     filling of generic histogram                                    *
*     filename: name of output file                                   *
*     y       : y-value to be filled in                               *
*     x       : x-value to be filled in                               *
*     xminh   : minimal x of histogram                                *
*     xmaxh   : minimal x of histogram                                *
*     nbin   : number of bins                                         *
*     histo   : number of histogram                                   *
*     levt    : number of event                                       *
*     addpoint: if=1 extra point at xmaxh is added                    *
*     gen     : generator                                             *
*     steps: flag for subroutine histogram                            *
*           steps=1: initialization                                   *
*           steps=2: start filling histograms                         *
*           steps=3: finish filling of histograms                     *
*           steps=4: print data for histograms                        *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     31.10.06  Ansgar Denner     last changed  03.05.07              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

#ifdef mpiuse
      include "mpif.h"
      integer status(MPI_STATUS_SIZE)
      real*8 keep_haverage(0:maxh,maxb,0:maxg,maxcon),
     &       keep_hsigma(0:maxh,maxb,0:maxg,maxcon),
     &       keep_havborn(0:maxh,maxb),keep_hsiborn(0:maxh,maxb)
#endif

      integer iNumProcs,iMyRank,iErr
      common/mpistuff/iNumProcs,iMyRank

c local variable
      real*8 haverage(0:maxh,maxb,0:maxg,maxcon),
     &    hsigma(0:maxh,maxb,0:maxg,maxcon),
     &    ysum(0:maxh,maxb,maxcon)
      real*8 havborn(0:maxh,maxb),hsiborn(0:maxh,maxb)
      real*8 xminh,xmaxh,x,y(maxcon)
      real*8 hav(0:maxg,maxcon),hsi(0:maxg,maxcon)
      real*8 havtot,hsitot,havb,hsib,xshift
      real*8 havqcd,hsiqcd,havew,hsiew,havinp,hsiinp,having,hsiing
      real*8 havsum(0:maxg,maxcon),hsisum(0:maxg,maxcon)
      real*8 havtotsum,hsitotsum,havbsum,hsibsum
      integer i,i2,isum(0:maxh,maxps),imax(0:maxh)
      integer histo,nbin,levt,histout,step,addpoint,lgen,column
      integer histoutv,histoutr
      character*30 filename
      character*40 filenamev,filenamer
      save havtot,hsitot,haverage,hsigma,havborn,hsiborn,ysum,isum,imax
      save hav,hsi

      
      
      


      
      if(shisto.eq.0.or.lhisto(histo).ne.1)return
c initialization
      if(step.eq.1)then
        if(nbin.gt.maxb)then
          write(*,'(a20,i6)')' error: reset maxb >',nbin
          stop
        endif
        if(histo.gt.maxh)then
          write(*,'(a20,i6)')' error: reset maxh >',histo
          stop
        endif
        havtot=0d0
        hsitot=0d0
        do gen=0,maxg
          do icon=1,maxcon
            hav(gen,icon)=0d0
            hsi(gen,icon)=0d0
          enddo
        enddo
        do ihisto=0,maxh
          imax(ihisto)=0
          do ibin=1,maxb
            do gen=0,maxg
              do icon=1,ncontr(gen)
                ysum(ihisto,ibin,icon)=0d0
                haverage(ihisto,ibin,gen,icon)=0d0
                hsigma(ihisto,ibin,gen,icon)=0d0
              enddo
            enddo
            havborn(ihisto,ibin)=0d0
            hsiborn(ihisto,ibin)=0d0
          enddo
          do iphsp=1,maxps
            isum(ihisto,iphsp)=0
          enddo
        enddo
c input of single subevents
      elseif(step.eq.2)then
        if(x.ge.xminh.and.x.lt.xmaxh)then
c bin
          i=int((x-xminh)/(xmaxh-xminh)*nbin)+1
          do icon=1,ncontr(lgen)
            ysum(histo,i,icon)=y(icon)+ysum(histo,i,icon)
          enddo
c subevent #
          imax(histo)=imax(histo)+1
c bin for subevent
          isum(histo,imax(histo))=i
        endif
c input born event and final calculation of whole event
      elseif(step.eq.3)then
        if(x.ge.xminh.and.x.lt.xmaxh)then
          i=int((x-xminh)/(xmaxh-xminh)*nbin)+1
          havborn(histo,i)=y(1)+havborn(histo,i)
          hsiborn(histo,i)=y(1)*y(1)+hsiborn(histo,i)
        endif
        if(imax(histo).ge.1)then
          do iphsp=1,imax(histo)
            i=isum(histo,iphsp)
            do icon=1,ncontr(lgen)
              haverage(histo,i,lgen,icon)=haverage(histo,i,lgen,icon)
     &            +ysum(histo,i,icon)
              hsigma(histo,i,lgen,icon)=hsigma(histo,i,lgen,icon)
     &            +ysum(histo,i,icon)*ysum(histo,i,icon)
              ysum(histo,i,icon)=0d0
            enddo
          enddo
        endif
        imax(histo)=0
c output
      elseif(step.eq.4)then
        xshift=1d0

#ifdef mpiuse
	if(iNumProcs.NE.1) Then
	  if(iMyRank.EQ.0) Then
	    do i2=1,iNumProcs-1
	      call MPI_Recv(keep_haverage, maxh*maxb*(maxg+1)*maxcon,
     &             MPI_DOUBLE_PRECISION, i2, 101, MPI_COMM_WORLD,
     &             status, iErr)
	      call MPI_Recv(keep_hsigma, maxh*maxb*(maxg+1)*maxcon,
     &             MPI_DOUBLE_PRECISION, i2, 102, MPI_COMM_WORLD,
     &             status, iErr)
	      call MPI_Recv(keep_havborn, maxh*maxb,
     &             MPI_DOUBLE_PRECISION, i2, 102, MPI_COMM_WORLD,
     &             status, iErr)
	      call MPI_Recv(keep_hsiborn, maxh*maxb,
     &             MPI_DOUBLE_PRECISION, i2, 102, MPI_COMM_WORLD,
     &             status, iErr)
              do ibin=1,maxb
        	havborn(histo,ibin)=
     &                havborn(histo,ibin)+
     &                keep_havborn(histo,ibin)
        	hsiborn(histo,ibin)=
     &                hsiborn(histo,ibin)+
     &                keep_hsiborn(histo,ibin)
        	do gen=0,maxg
        	  do icon=1,ncontr(gen)
        	    haverage(histo,ibin,gen,icon)=
     &                   haverage(histo,ibin,gen,icon)+
     &                   keep_haverage(histo,ibin,gen,icon)
        	    hsigma(histo,ibin,gen,icon)=
     &                   hsigma(histo,ibin,gen,icon)+
     &                   keep_hsigma(histo,ibin,gen,icon)
		  enddo
		enddo
	      enddo
            enddo
	  else
	    call MPI_Send(haverage, maxh*maxb*(maxg+1)*maxcon,
     &           MPI_DOUBLE_PRECISION, 0, 101, MPI_COMM_WORLD, iErr)
	    call MPI_Send(hsigma, maxh*maxb*(maxg+1)*maxcon,
     &           MPI_DOUBLE_PRECISION, 0, 102, MPI_COMM_WORLD, iErr)
	    call MPI_Send(havborn, maxh*maxb,
     &           MPI_DOUBLE_PRECISION, 0, 102, MPI_COMM_WORLD, iErr)
	    call MPI_Send(hsiborn, maxh*maxb,
     &           MPI_DOUBLE_PRECISION, 0, 102, MPI_COMM_WORLD, iErr)
            do ibin=1,maxb
              havborn(histo,ibin)=0D0
              hsiborn(histo,ibin)=0D0
              do gen=0,maxg
        	do icon=1,ncontr(gen)
        	  haverage(histo,ibin,gen,icon)=0D0
        	  hsigma(histo,ibin,gen,icon)=0D0
		enddo
	      enddo
	    enddo
	  endif
	endif

	if(iMyRank.EQ.0) Then
#endif


        do i=35,2,-1
          if(filenamev(i:i).eq.' ')
     &        filenamev=filenamev(1:i-1)//filenamev(i+1:35)
          if(filenamer(i:i).eq.' ')
     &        filenamer=filenamer(1:i-1)//filenamer(i+1:35)
        enddo

        

        histout=25
        open(unit=histout,file=filename,status='unknown')
        write(histout,'("# histogram:",i3,", bins:",i4)')histo,nbin
        write(histout,'(a)')'# '


 1002   format("# column ",I2,": ",a53)
 1003   format("# column ",I2,": ",a50,a20)
 1004   format("# column ",I2,": ",a57,"\n",13x,a27,I1)
 1005   format("# column ",I2,": kinematical quantity on x-axis")
 1006   format("# column ",I2,": ",a63)

        column=0
        write(histout,1005)column+1
        column=column+2
        write(histout,1002)column,
     &      "Monte Carlo average for the Born cross-section       "
        write(histout,1002)column+1,
     &      "Monte Carlo error   for the Born cross-section       "
        column=column+2
        write(histout,1002)column,
     &      "Monte Carlo average for the complete cross-section   "
        write(histout,1002)column+1,
     &      "Monte Carlo error   for the complete cross-section   "
        column=column+2
        write(histout,1002)column,
     &      "Monte Carlo average for virt+real ew cross-section   "
        write(histout,1002)column+1,
     &      "Monte Carlo error   for virt+real ew cross-section   "
        column=column+2
        write(histout,1002)column,
     &      "Monte Carlo average for virt+real QCD cross-section  "
        write(histout,1002)column+1,
     &      "Monte Carlo error   for virt+real QCD cross-section  "
        column=column+2
        write(histout,1002)column,
     &      "Monte Carlo average for incoming photon cross-section"
        write(histout,1002)column+1,
     &      "Monte Carlo error   for incoming photon cross-section"
        column=column+2
        write(histout,1002)column,
     &      "Monte Carlo average for incoming gluon cross-section "
        write(histout,1002)column+1,
     &      "Monte Carlo error   for incoming gluon cross-section "


        write(histout,'(a)') '# '

        havbsum=0d0
        hsibsum=0d0
        havtotsum=0d0
        hsitotsum=0d0
        do gen=0,maxg
          do icon=1,ncontr(gen)
            havsum(gen,icon)=0d0
            hsisum(gen,icon)=0d0
          enddo
        enddo
        do i=1,nbin
          havtot=0d0
          hsitot=0d0
          do gen=0,maxg
            do icon=1,ncontr(gen)
              hav(gen,icon)=haverage(histo,i,gen,icon)
     &            /dble(levt)/(xmaxh-xminh)*dble(nbin)
              hsi(gen,icon)=hsigma(histo,i,gen,icon)
     &            /(dble(levt)*(xmaxh-xminh)/dble(nbin))**2
              havtot=havtot+hav(gen,icon)
              hsitot=hsitot+hsi(gen,icon)
              hsi(gen,icon)=sqrt(dabs(hsi(gen,icon)
     &            -hav(gen,icon)*hav(gen,icon)/dble(levt)))
            enddo
          enddo
          hsitot=sqrt(dabs(hsitot-havtot*havtot/dble(levt)))
          havb=havborn(histo,i)/dble(levt)/(xmaxh-xminh)*dble(nbin)
          hsib=hsiborn(histo,i)/(dble(levt)*(xmaxh-xminh)/dble(nbin))**2
          hsib=sqrt(dabs(hsib-havb*havb/dble(levt)))
          if (sbarelep.eq.0) then
             havew=hav(0,2)+hav(0,3)+hav(1,1)
             hsiew=sqrt(hsi(0,2)**2+hsi(0,3)**2+hsi(1,1)**2)
          else
             havew=hav(0,2)+hav(0,3)+hav(0,10)+hav(1,1)
             hsiew=sqrt(hsi(0,2)**2+hsi(0,3)**2+hsi(0,10)**2
     &            +hsi(1,1)**2)
          endif
          havqcd=hav(0,4)+hav(0,5)+hav(2,1)
          hsiqcd=sqrt(hsi(0,4)**2+hsi(0,5)**2+hsi(2,1)**2)
          havinp=hav(0,6)+hav(0,7)+hav(3,1)+hav(4,1)
     &         +hav(7,1)+hav(8,1)
     &         +hav(9,1)+hav(10,1)
          hsiinp=sqrt(hsi(0,6)**2+hsi(0,7)**2+hsi(3,1)**2+hsi(4,1)**2
     &         +hsi(7,1)**2+hsi(8,1)**2
     &         +hsi(9,1)**2+hsi(10,1)**2)
          having=hav(0,8)+hav(0,9)+hav(5,1)+hav(6,1)
          hsiing=sqrt(hsi(0,8)**2+hsi(0,9)**2+hsi(5,1)**2+hsi(6,1)**2)

          write(histout,'(e14.8,2x,20(e12.6,2x))')
     &        xminh+(xmaxh-xminh)/dble(nbin)*(dble(i)-xshift),
     &        havb,hsib,havtot,hsitot,havew,hsiew,havqcd,hsiqcd,
     &        havinp,hsiinp,having,hsiing



          havbsum=havbsum+havb
          hsibsum=hsibsum+hsib*hsib
          havtotsum=havtotsum+havtot
          hsitotsum=hsitotsum+hsitot*hsitot
        enddo
        havbsum=havbsum*(xmaxh-xminh)/dble(nbin)
        hsibsum=sqrt(hsibsum)*(xmaxh-xminh)/dble(nbin)
        havtotsum=havtotsum*(xmaxh-xminh)/dble(nbin)
        hsitotsum=sqrt(hsitotsum)*(xmaxh-xminh)/dble(nbin)
        do gen=0,maxg
          do icon=1,ncontr(gen)
            havsum(gen,icon)=havsum(gen,icon)*(xmaxh-xminh)/dble(nbin)
            hsisum(gen,icon)=sqrt(hsisum(gen,icon))
     &          *(xmaxh-xminh)/dble(nbin)
          enddo
        enddo

        if(addpoint.eq.1)then
          write(histout,'(e14.8,2x,20(e12.6,2x))')
     &        xmaxh,
     &        havb,hsib,havtot,hsitot,havew,hsiew,havqcd,hsiqcd,
     &        havinp,hsiinp,having,hsiing
          write(histout,'(a)') '# '
          write(histout,'(a31,1e12.6)')
     &        '# additional data point at x = ',xmaxh


      endif

      write(histout ,'(a)') '#'
      write(histout ,'(a)') '# useful for plotting with gnuplot'


      if(xshift.eq.0d0)then
        write(histout ,'(a)')'# x = right corner'
      elseif(xshift.eq.0.5d0)then
        write(histout ,'(a)')'# x = mid'
      else if(xshift.eq.1d0)then
        write(histout ,'(a)')'# x = left corner'
      endif
      close(histout)
      
 300  format(a11,i3,2x,a21,':  ',f8.2,' < x < ',f8.2)
 400  format(a12,i2,a14,i2,a1,1g18.10,' +/- ',1g18.10)
 500  format(a31,1g18.10,' +/- ',1g18.10)
      if(lhisto(0).eq.1.and.histo.eq.0.or.lhisto(0).eq.0
     &  .and.histo.EQ.1) then
        write(nout,*) ""
        write(nout,*) "The following histograms have been produced:"
        write(nout,*) ""
      endif
      write(nout,300)' histogram',histo,filename,xminh,xmaxh
      if(lnoutmc.GE.2) then
        write(nout,500)' Born cross section:',
     &    havbsum,hsibsum
        write(nout,500)' complete cross section:',
     &    havtotsum,hsitotsum
      endif
      do gen=0,ngenerator
        do icon=1,ncontr(gen)
          if(havsum(gen,icon).ne.0d0)then
            write(nout,400)' generator',gen,
     &    ', contribution',icon,':',havsum(gen,icon),hsisum(gen,icon)
          endif
        enddo
      enddo

      endif

#ifdef mpiuse
      endif
#endif

      end


