***********************************************************************
      subroutine whzh_recombination(kbeam,k,krec,generator,lnext,cut,
     *  switch, z, phasespace, em)
***********************************************************************
*---------------------------------------------------------------------*
*     renamed and shifted from 'subroutine recombination'             *
*     09.03.11  Stefan Kallweit   last changed  21.02.13              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8  kbeam(2,0:3),k(maxe,0:3),krec(maxe,0:3),p(maxe,0:3)
      real*8  cost(maxe),s(maxe),srecmin,cmax,kvec(maxe)
      real*8  y(maxe),kt2(maxe),dR2(maxe),dR2min,ktiktj,dphi
      integer nrecomb,generator,lnext,switch,cut
      integer recombine
      integer typepart(maxe),typejet(maxe),npartin,njet
      integer ii
      integer em, phasespace, pnumber
      real*8 z
      integer selprocess
      common/selprocess/selprocess

c to avoid compiler warning
      data nrecomb /0/

      if(switch.ne.0)return
      recombine=0

      do ext=1,lnext
         do vcomp=0,3
            krec(ext,vcomp)=k(ext,vcomp)
         enddo
      enddo
      do ext=lnext+1,maxe
         do vcomp=0,3
            krec(ext,vcomp)=0d0
         enddo
      enddo
      if(srecomb(generator).eq.1.or.srecomb(generator).eq.2) then

         if (sbarelep.ne.0) then
            if (generator.eq.0.and.phasespace.ge.6) then
               do vcomp=0,3
                  krec(em,vcomp)=z*k(em,vcomp)
                  krec(lnext,vcomp)=
     &                 (1d0-z)*k(em,vcomp)
               enddo
            elseif (generator.eq.1.and.phasespace.ge.2.and.
     &           (shtr.eq.0.and.lnext.eq.6.or.
     &            shtr.gt.0.and.lnext.eq.7)

     &           ) then
               do vcomp=0,3
                  krec(em,vcomp)=z*k(em,vcomp)
                  krec(lnext,vcomp)=
     &                 (1d0-z)*k(em,vcomp)
               enddo
            endif
         elseif (sbarelep.eq.0.and.generator.eq.1.and.
     &           phasespace.eq.1) then
            call lepgammakluster(kbeam,k,krec,lnext,
     &           njet,switch)
         endif

         call applylepcut(kbeam,krec,generator,njet,cut,switch)
      else
         print *, 'Not yet modified: applycut'
      endif
      end



***********************************************************************
      subroutine lepgammakluster(kbeam,kpart,krec,npartin,
     &     npartout,switch)
***********************************************************************
#include "mcparams.h"
#include "mccommon.h"

      real*8  kbeam(2,0:3),kpart(maxe,0:3),krec(maxe,0:3)
      integer typelep,npartin,npartout,switch
      real*8  ygamma, ylp, ylm
      real*8  kt2gamma, kt2lp, kt2lm
      real*8  ktgammalp, ktgammalm
      real*8  dR2lp, dR2lm
      real*8  dRlp, dRlm

      integer selprocess
      common/selprocess/selprocess

      npartout = npartin

      ygamma=.5d0*log((kpart(npartin,0)+kpart(npartin,3))
     &     /(kpart(npartin,0)-kpart(npartin,3)))
      kt2gamma=kpart(npartin,1)**2+kpart(npartin,2)**2

      if (selprocess.eq.1)then
         ylp=.5d0*log((kpart(3,0)+kpart(3,3))/(kpart(3,0)-kpart(3,3)))
         kt2lp=kpart(3,1)**2+kpart(3,2)**2
         ktgammalp=(kpart(3,1)*kpart(npartin,1)
     &        +kpart(3,2)*kpart(npartin,2))/sqrt(kt2lp*kt2gamma)
         dphilp=acos(ktgammalp)
         dR2lp= ((ylp-ygamma)**2+dphilp**2)
         dRlp = sqrt(dR2lp)
         if (dRlp.le.dgammaparameter)then
            do vcomp=0,3
               krec(3,vcomp)=kpart(3,vcomp)+kpart(npartin,vcomp)
            enddo
            npartout=npartout-1
         endif
      elseif (selprocess.eq.2)then
         ylm=.5d0*log((kpart(4,0)+kpart(4,3))/(kpart(4,0)-kpart(4,3)))
         kt2lm=kpart(4,1)**2+kpart(4,2)**2
         ktgammalm=(kpart(4,1)*kpart(npartin,1)
     &        +kpart(4,2)*kpart(npartin,2))/sqrt(kt2lm*kt2gamma)
         dphilm=acos(ktgammalm)
         dR2lm= ((ylm-ygamma)**2+dphilm**2)
         dRlm = sqrt(dR2lm)
         if (dRlm.le.dgammaparameter)then
            do vcomp=0,3
               krec(4,vcomp)=kpart(4,vcomp)+kpart(npartin,vcomp)
            enddo
            npartout=npartout-1
         endif
      elseif (selprocess.eq.3)then
         ylp=.5d0*log((kpart(3,0)+kpart(3,3))/(kpart(3,0)-kpart(3,3)))
         ylm=.5d0*log((kpart(4,0)+kpart(4,3))/(kpart(4,0)-kpart(4,3)))
         kt2lp=kpart(3,1)**2+kpart(3,2)**2
         ktgammalp=(kpart(3,1)*kpart(npartin,1)
     &        +kpart(3,2)*kpart(npartin,2))/sqrt(kt2lp*kt2gamma)
         kt2lm=kpart(4,1)**2+kpart(4,2)**2
         ktgammalm=(kpart(4,1)*kpart(npartin,1)
     &        +kpart(4,2)*kpart(npartin,2))/sqrt(kt2lm*kt2gamma)
         dphilp=acos(ktgammalp)
         dR2lp= ((ylp-ygamma)**2+dphilp**2)
         dphilm=acos(ktgammalm)
         dR2lm= ((ylm-ygamma)**2+dphilm**2)
         dRlp = sqrt(dR2lp)
         dRlm = sqrt(dR2lm)
         if (dRlp.le.dgammaparameter.and.dRlp.le.dRlm)then
            do vcomp=0,3
               krec(3,vcomp)=kpart(3,vcomp)+kpart(npartin,vcomp)
            enddo
            npartout=npartout-1
         elseif (dRlm.le.dgammaparameter.and.dRlm.le.dRlp)then
            do vcomp=0,3
               krec(4,vcomp)=kpart(4,vcomp)+kpart(npartin,vcomp)
            enddo
            npartout=npartout-1
         endif
      elseif (selprocess.eq.4)then
      endif

      end



***********************************************************************
      subroutine applylepcut(kbeam,krec,generator,njet,cut,switch)
***********************************************************************
*     application of separation cuts on leptons                       *
*                                                                     *
*     selprocess = 1     Hln     p p -> l+ nu  H                      *
*     selprocess = 2     Hnl     p p -> nu~ l- H                      *
*     selprocess = 3     Hll     p p -> l+ l-  H                      *
*     selprocess = 4     Hnn     p p -> nu~ nu H                      *
*---------------------------------------------------------------------*
*     08.06.10  Stefan Kallweit   last changed  21.02.13              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),p(0:3),kvec(maxe),kbeamvec(2)
      real*8 kh(0:3)
      real*8 krec(maxe,0:3),kll(0:3),kllh(0:3)
      real*8 kt2l(maxe),ktl(maxe),yk(maxe)
      real*8 kt2v,ktv,kth
      real*8 s,prod,ktiktj,dR,dy,dphi,sh,mll,mllh
      integer njet,cut,generator,switch,prevcuts

      integer selprocess
      common/selprocess/selprocess

      if(switch.gt.0)return

      if(scuts(generator).gt.0.or.shtr.gt.0)then

ccccccccccccccccccccccccccccccccccccccccccc
c  Higgs cuts for process p p -> x  x  H  c
ccccccccccccccccccccccccccccccccccccccccccc

         call applylephiggscut(kbeam,krec,generator,njet,cut,switch)
         if(cut.gt.0) return

         if (selprocess.eq.1) then

cccccccccccccccccccccccccccccccccccccccccccccc
c  leptonic cuts for process p p -> l+ nu H  c
cccccccccccccccccccccccccccccccccccccccccccccc
            prevcuts=20
            cut=prevcuts

cccccccccccccccc
c  cuts on l+  c
cccccccccccccccc
            ext=3
ccccccccccccccccccc
c  apply pT cuts  c
ccccccccccccccccccc
            kt2l(ext)=krec(ext,1)**2+krec(ext,2)**2
            ktl(ext)=sqrt(kt2l(ext))
            cut=cut+1
            if(ktl(ext).lt.ptcutl)return
c     cut = 21

ccccccccccccccccccccccc
c  apply energy cuts  c
ccccccccccccccccccccccc
            cut=cut+1
            if(krec(ext,0).lt.ecutl)return
c     cut = 22

ccccccccccccccccccccccccc
c  apply rapidity cuts  c
ccccccccccccccccccccccccc
            cut=cut+1
            if(krec(ext,1).eq.0d0.and.krec(ext,2).eq.0d0) then
               if (krec(ext,3).ge.0) then
                  yk(ext)=.1d0*ycutl
               elseif (krec(ext,3).lt.0) then
                  yk(ext)=-.1d0*ycutl
               endif
            else if(krec(ext,0)-krec(ext,3).ne.0d0)then
               yk(ext)=.5d0*log((krec(ext,0)+krec(ext,3))
     &              /(krec(ext,0)-krec(ext,3)))
            else
               return
            endif
            if(abs(yk(ext)).gt.ycutl)return
c     cut = 23



cccccccccccccccc
c  cuts on nu  c
cccccccccccccccc
            ext=4
cccccccccccccccccccccccccccccccccccc
c  apply transverse momentum cuts  c
cccccccccccccccccccccccccccccccccccc
            kt2l(ext)=krec(ext,1)**2+krec(ext,2)**2
            ktl(ext)=sqrt(kt2l(ext))
            cut=cut+1
            if(ktl(ext).lt.ptcutm)return
c     cut = 24



cccccccccccccccccccccccccccccccccccccccccccccc
c  cuts involving lepton and a possible jet  c
cccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccc
c  apply dR cuts  c
ccccccccccccccccccc
            if(generator.ge.2.and.generator.le.6) then
               if (shtr.eq.0.and.njet.eq.6)then
                  ext=3
                  dR=sqrt((yk(ext)-yk(6))**2+
     &                 acos((krec(ext,1)*krec(6,1)
     &                 +krec(ext,2)*krec(6,2))
     &                 /sqrt(kt2l(ext)*kt2l(6)))**2)
                  cut=cut+1
                  if(dR.lt.dRcutjl)return
c     cut = 25
               elseif (shtr.gt.0.and.njet.eq.7)then
                  ext=3
                  dR=sqrt((yk(ext)-yk(7))**2+
     &                 acos((krec(ext,1)*krec(7,1)
     &                 +krec(ext,2)*krec(7,2))
     &                 /sqrt(kt2l(ext)*kt2l(7)))**2)
                  cut=cut+1
                  if(dR.lt.dRcutjl)return
c     cut = 25
               endif
            else
               cut=cut+1
c     cut = 25
            endif



cccccccccccccccccccccccc
c  cuts on W+ = l+ nu  c
cccccccccccccccccccccccc
            kt2v=(krec(3,1)+krec(4,1))**2
     &           +(krec(3,2)+krec(4,2))**2
            ktv=sqrt(kt2v)
            cut=cut+1
            if(ktv.lt.ptcutv)return
c     cut = 26



cccccccccccccccccccc
c  veto cuts on l+ c
cccccccccccccccccccc
            ext=3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  apply maximum-pT cuts .and. minimum-rapidity cuts  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            cut=cut+1
            if((ktl(ext).gt.ptmaxl).and.
     &           (abs(yk(ext)).lt.yminl))return
c     cut = 27



         elseif (selprocess.eq.2) then

ccccccccccccccccccccccccccccccccccccccccccccccc
c  leptonic cuts for process p p -> nu~ l- H  c
ccccccccccccccccccccccccccccccccccccccccccccccc
            prevcuts=20
            cut=prevcuts

cccccccccccccccc
c  cuts on l-  c
cccccccccccccccc
            ext=4
ccccccccccccccccccc
c  apply pT cuts  c
ccccccccccccccccccc
            kt2l(ext)=krec(ext,1)**2+krec(ext,2)**2
            ktl(ext)=sqrt(kt2l(ext))
            cut=cut+1
            if(ktl(ext).lt.ptcutl)return
c     cut = 21

ccccccccccccccccccccccc
c  apply energy cuts  c
ccccccccccccccccccccccc
            cut=cut+1
            if(krec(ext,0).lt.ecutl)return
c     cut = 22

ccccccccccccccccccccccccc
c  apply rapidity cuts  c
ccccccccccccccccccccccccc
            cut=cut+1
            if(krec(ext,1).eq.0d0.and.krec(ext,2).eq.0d0) then
               if (krec(ext,3).ge.0) then
                  yk(ext)=.1d0*ycutl
               elseif (krec(ext,3).lt.0) then
                  yk(ext)=-.1d0*ycutl
               endif
            else if(krec(ext,0)-krec(ext,3).ne.0d0)then
               yk(ext)=.5d0*log((krec(ext,0)+krec(ext,3))
     &              /(krec(ext,0)-krec(ext,3)))
            else
               return
            endif
            if(abs(yk(ext)).gt.ycutl)return
c     cut = 23



ccccccccccccccccc
c  cuts on nu~  c
ccccccccccccccccc
            ext=3
cccccccccccccccccccccccccccccccccccc
c  apply transverse momentum cuts  c
cccccccccccccccccccccccccccccccccccc
            kt2l(ext)=krec(ext,1)**2+krec(ext,2)**2
            ktl(ext)=sqrt(kt2l(ext))
            cut=cut+1
            if(ktl(ext).lt.ptcutm)return
c     cut = 24



cccccccccccccccccccccccccccccccccccccccccccccc
c  cuts involving lepton and a possible jet  c
cccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccc
c  apply dR cuts  c
ccccccccccccccccccc
            if(generator.ge.2.and.generator.le.6) then
               if (shtr.eq.0.and.njet.eq.6)then
                  ext=4
                  dR=sqrt((yk(ext)-yk(6))**2+
     &                 acos((krec(ext,1)*krec(6,1)
     &                 +krec(ext,2)*krec(6,2))
     &                 /sqrt(kt2l(ext)*kt2l(6)))**2)
                  cut=cut+1
                  if(dR.lt.dRcutjl)return
c     cut = 25
               elseif (shtr.gt.0.and.njet.eq.7)then
                  ext=4
                  dR=sqrt((yk(ext)-yk(7))**2+
     &                 acos((krec(ext,1)*krec(7,1)
     &                 +krec(ext,2)*krec(7,2))
     &                 /sqrt(kt2l(ext)*kt2l(7)))**2)
                  cut=cut+1
                  if(dR.lt.dRcutjl)return
c     cut = 25
               endif
            else
               cut=cut+1
            endif



ccccccccccccccccccccccccc
c  cuts on W- = nu~ l-  c
ccccccccccccccccccccccccc
            kt2v=(krec(3,1)+krec(4,1))**2
     &           +(krec(3,2)+krec(4,2))**2
            ktv=sqrt(kt2v)
            cut=cut+1
            if(ktv.lt.ptcutv)return
c     cut = 26



cccccccccccccccccccc
c  veto cuts on l- c
cccccccccccccccccccc
            ext=4
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  apply maximum-pT cuts .and. minimum-rapidity cuts  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            cut=cut+1
            if((ktl(ext).gt.ptmaxl).and.
     &           (abs(yk(ext)).lt.yminl))return
c     cut = 27



         elseif (selprocess.eq.3) then

cccccccccccccccccccccccccccccccccccccccccccccc
c  leptonic cuts for process p p -> l+ l- H  c
cccccccccccccccccccccccccccccccccccccccccccccc
            prevcuts=20
            cut=prevcuts
            do ext=3,4

ccccccccccccccccccccc
c  cuts on l+ / l-  c
ccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccc
c  apply transverse momentum cuts  c
cccccccccccccccccccccccccccccccccccc
               kt2l(ext)=krec(ext,1)**2+krec(ext,2)**2
               ktl(ext)=sqrt(kt2l(ext))
               cut=cut+1
               if(ktl(ext).lt.ptcutl)return
c     cut = 21 / 24

ccccccccccccccccccccccc
c  apply energy cuts  c
ccccccccccccccccccccccc
               cut=cut+1
               if(krec(ext,0).lt.ecutl)return
c     cut = 22 / 25

ccccccccccccccccccccccccc
c  apply rapidity cuts  c
ccccccccccccccccccccccccc
               cut=cut+1
               if(krec(ext,1).eq.0d0.and.krec(ext,2).eq.0d0) then
                  if (krec(ext,3).ge.0) then
                     yk(ext)=.1d0*ycutl
                  elseif (krec(ext,3).lt.0) then
                     yk(ext)=-.1d0*ycutl
                  endif
               else if(krec(ext,0)-krec(ext,3).ne.0d0)then
                  yk(ext)=.5d0*log((krec(ext,0)+krec(ext,3))
     &                 /(krec(ext,0)-krec(ext,3)))
               else
                  return
               endif
               if(abs(yk(ext)).gt.ycutl)return
c     cut = 23 / 26
            enddo

ccccccccccccccccccccccccccccccccc
c  cuts involving both leptons  c
ccccccccccccccccccccccccccccccccc

ccccccccccccccccccc
c  apply dR cuts  c
ccccccccccccccccccc
c            prevcuts=30
c            cut=prevcuts

            dR=sqrt((yk(3)-yk(4))**2+
     &           acos((krec(3,1)*krec(4,1)+krec(3,2)*krec(4,2))/
     &           sqrt(kt2l(3)*kt2l(4)))**2)
            cut=cut+1
            if(dR.lt.dRcutll)return
c     cut = 27

ccccccccccccccccccccccccccccccc
c  apply invariant mass cuts  c
ccccccccccccccccccccccccccccccc
            do vcomp=0,3
               kll(vcomp)=krec(3,vcomp)+krec(4,vcomp)
            enddo
            mll=sqrt(kll(0)*kll(0)-kll(1)*kll(1)-kll(2)*kll(2)
     &           -kll(3)*kll(3))
            cut=cut+1
            if(mll.lt.mlcutll.or.mll.gt.mucutll) return
c     cut = 28

cccccccccccccccccccccccccccccccccccccccccccccc
c  cuts involving lepton and a possible jet  c
cccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccc
c  apply dR cuts  c
ccccccccccccccccccc
            if(generator.ge.2.and.generator.le.6) then
               if (shtr.eq.0.and.njet.eq.6)then
                  do ext=3,4
                     dR=sqrt((yk(ext)-yk(6))**2+
     &                    acos((krec(ext,1)*krec(6,1)
     &                    +krec(ext,2)*krec(6,2))
     &                    /sqrt(kt2l(ext)*kt2l(6)))**2)
                     cut=cut+1
                     if(dR.lt.dRcutjl)return
c     cut = 29 / 30
                  enddo
               elseif (shtr.gt.0.and.njet.eq.7)then
                  do ext=3,4
                     dR=sqrt((yk(ext)-yk(7))**2+
     &                    acos((krec(ext,1)*krec(7,1)
     &                    +krec(ext,2)*krec(7,2))
     &                    /sqrt(kt2l(ext)*kt2l(7)))**2)
                     cut=cut+1
                     if(dR.lt.dRcutjl)return
c     cut = 29 / 30
                  enddo
               endif
            endif

ccccccccccccccccccccccc
c  cuts on Z = l+ l-  c
ccccccccccccccccccccccc
            kt2v=(krec(3,1)+krec(4,1))**2
     &           +(krec(3,2)+krec(4,2))**2
            ktv=sqrt(kt2v)
            cut=cut+1
            if(ktv.lt.ptcutv)return
c     cut = 31

cccccccccccccccccccccccccccccccccccccc
c  apply invariant mass cuts  on llh c
cccccccccccccccccccccccccccccccccccccc
            if(shtr.eq.0)then
              do vcomp=0,3
                kllh(vcomp)=krec(5,vcomp)+krec(3,vcomp)+krec(4,vcomp)
              enddo
            else
              do vcomp=0,3
                kllh(vcomp)=krec(5,vcomp)+krec(6,vcomp)
     &              +krec(3,vcomp)+krec(4,vcomp)
              enddo
            endif

            mllh=sqrt(kllh(0)*kllh(0)-kllh(1)*kllh(1)-kllh(2)*kllh(2)
     &           -kllh(3)*kllh(3))
            cut=cut+1
            if(mllh.lt.mlcutllh.or.mllh.gt.mucutllh) return
c     cut = 32 / 33



         elseif (selprocess.eq.4) then

cccccccccccccccccccccccccccccccccccccccccccccccc
c  leptonic cuts for process p p -> nu~ nu  H  c
cccccccccccccccccccccccccccccccccccccccccccccccc
            prevcuts=20
            cut=prevcuts
            ext = 3
            ext2 = 4

cccccccccccccccccccccccc
c  cuts on missing pT  c
cccccccccccccccccccccccc
            do vcomp=0,3
               kll(vcomp)=krec(3,vcomp)+krec(4,vcomp)
            enddo
            kt2l(ext)=kll(1)**2+kll(2)**2
            ktl(ext)=sqrt(kt2l(ext))
            cut=cut+1
            if(ktl(ext).lt.ptcutm)return
c     cut = 21
         endif

         if(cut.gt.maxcut)then
            write(*,*) 'reset maxcut = ',cut
            stop
         endif

cccccccccccccccccccccccc
c  cuts on Z = nu~ nu  c
cccccccccccccccccccccccc
            kt2v=(krec(3,1)+krec(4,1))**2
     &           +(krec(3,2)+krec(4,2))**2
            ktv=sqrt(kt2v)
            cut=cut+1
            if(ktv.lt.ptcutv)return
c     cut = 22

      endif

      cut=0

      end



***********************************************************************
      subroutine applylephiggscut(kbeam,k,generator,lnext,cut,switch)
***********************************************************************
*     application of separation cuts on Higgs boson                   *
*     adapted from applylephiggscut in public.F                       *
*---------------------------------------------------------------------*
*     08.04.13  Stefan Kallweit   last changed  08.02.10              *
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

ccccccccccccccc
c  cuts on H  c
ccccccccccccccc
        ext=5                   ! higgs

c apply transverse momentum cuts

ccccccccccccccccccc
c  apply pT cuts  c
ccccccccccccccccccc
        kth=sqrt(kh(1)**2+kh(2)**2)
        cut=cut+1
        if(kth.lt.ptcuth)return

ccccccccccccccccccccccc
c  apply energy cuts  c
ccccccccccccccccccccccc
        cut=cut+1
        if(kh(0).lt.ecuth)return

ccccccccccccccccccccccccc
c  apply rapidity cuts  c
ccccccccccccccccccccccccc
        cut=cut+1
        if(kh(0)-kh(3).ne.0d0)then
          yh=.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))
        else
          return
        endif
        if(abs(yh).gt.ycuth)return

ccccccccccccccccccccccccccccccc
c  apply invariant mass cuts  c
ccccccccccccccccccccccccccccccc
        sh=kh(0)*kh(0)-kh(1)*kh(1)-kh(2)*kh(2)
     &      -kh(3)*kh(3)

        cut=cut+1
        if(sh.gt.sucuth.or.sh.lt.slcuth) return

cccccccccccccccccccccccccccccccccc
c  cuts on Higgs decay singlets  c
cccccccccccccccccccccccccccccccccc
        if(shtr.eq.2)then

cccccccccccccccccccccccccccccccccccc
c  apply transverse momentum cuts  c
cccccccccccccccccccccccccccccccccccc
          do ext=5,6
            kt(ext)=sqrt(k(ext,1)**2+k(ext,2)**2)

            cut=cut+1
            if(kt(ext).lt.ptcut(ext,generator))return
          enddo

ccccccccccccccccccccccccc
c  apply rapidity cuts  c
ccccccccccccccccccccccccc
          do ext=5,6
            if(k(ext,0)-k(ext,3).ne.0d0)then
              y(ext)=.5d0*log((k(ext,0)+k(ext,3))/(k(ext,0)-k(ext,3)))
            else
              return
            endif
            cut=cut+1
            if(abs(y(ext)).gt.ycut(ext,generator))return
          enddo
        endif
        if(cut.gt.maxcut)then
          write(*,*) 'reset maxcut = ',cut
          stop
        endif

      endif

      cut=0

      end


***********************************************************************
      subroutine init_cuts_Hll(gen)
***********************************************************************
* initialization of cuts for process p p -> l+ l- H (selprocess = 3)  *
*---------------------------------------------------------------------*
*     04.06.10  Stefan Kallweit   last changed  21.02.13              *
*---------------------------------------------------------------------*
*     !!! definitions are -not- used in actual cut routines !!!       *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none
      real*8 mtotcut

#include "mcparams.h"
#include "mccommon.h"

cccccccccccccccccccccccccccccccccccccc
c  standard cuts for p p -> l+ l- H  c
cccccccccccccccccccccccccccccccccccccc
        mtotcut = 0d0
        if(scuts(gen).eq.0.and.shtr.eq.1)then
           sucuth=mucuth*mucuth
           slcuth=mlcuth*mlcuth
        elseif(scuts(gen).eq.1.or.scuts(gen).eq.2)then
           ptcut(3,gen)=ptcutl
           ptcut(4,gen)=ptcutl
           ycut(3,gen)=ycutl
           ycut(4,gen)=ycutl
           dycut(3,4,gen)=dycutll
           dycut(4,3,gen)=dycutll
           scut(3,4,gen)=mlcutll
           scut(4,3,gen)=mlcutll
           dRcut(3,4,gen)=dRcutll
           dRcut(4,3,gen)=dRcutll
           if(shtr.ge.2)then
              ptcut(5,gen)=ptcuts
              ptcut(6,gen)=ptcuts
              ycut(5,gen)=ycuts
              ycut(6,gen)=ycuts
           endif
           sucuth=mucuth*mucuth
           slcuth=mlcuth*mlcuth
        elseif(scuts(0).ne.0) then
          write(*,'(a30,i3)')' Unknown input: scuts = ',scuts(gen)
          stop
        endif

c filling cut parameters
c calculating cut parameters
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
           ccut(ext,ext2,gen)=dcos(ccut(ext,ext2,gen)*pi/180d0)
           scut(ext,ext2,gen)=scut(ext,ext2,gen)**2
        enddo
        enddo
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
           ccut(ext,ext2,gen)=min(ccut(ext,ext2,gen),ccut(ext2,ext,gen))
           scut(ext,ext2,gen)=max(scut(ext,ext2,gen),scut(ext2,ext,gen))
        enddo
        enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  calculating cut parameters to be given to the generator  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        mtotcutgen=0d0
        do ext=3,nextgen(gen)
           mtotcutgen=mtotcutgen+sqrt(ptcut(ext,gen)**2+mass2(ext,gen))
           do ext2=1,nextgen(gen)
c ccutgen must be set to 1 if there are no cuts!
              ccutgen(ext2,ext,gen)=ccut(ext2,ext,gen)
              scutgen(ext2,ext,gen)=scut(ext2,ext,gen)
           enddo
        enddo
        mtotcutgen=max(mtotcutgen,mtotcut)
      end

***********************************************************************
      subroutine init_cuts_Hln(gen)
***********************************************************************
*---------------------------------------------------------------------*
*     04.06.10  Stefan Kallweit   last changed  21.02.13              *
*---------------------------------------------------------------------*
*     !!! definitions are -not- used in actual cut routines !!!       *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none
      real*8 mtotcut

#include "mcparams.h"
#include "mccommon.h"
ccccccccccccccccccccccccccccccccccccccc
c  standard cuts for p p -> l+ nu  H  c
ccccccccccccccccccccccccccccccccccccccc
        mtotcut = 0d0
        if(scuts(gen).eq.0.and.shtr.eq.1)then
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(gen).eq.1.or.scuts(gen).eq.2)then
          ptcut(3,gen)=ptcutl
          ptcut(4,gen)=ptcutm
          ycut(3,gen)=ycutl
          if(shtr.ge.2)then
            ptcut(5,gen)=ptcuts
            ptcut(6,gen)=ptcuts
            ycut(5,gen)=ycuts
            ycut(6,gen)=ycuts
          endif
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(0).ne.0) then
          write(*,'(a30,i3)')' Unknown input: scuts = ',scuts(gen)
          stop
        endif

c filling cut parameters
c calculating cut parameters
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
          ccut(ext,ext2,gen)=dcos(ccut(ext,ext2,gen)*pi/180d0)
          scut(ext,ext2,gen)=scut(ext,ext2,gen)**2
        enddo
        enddo
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
          ccut(ext,ext2,gen)=min(ccut(ext,ext2,gen),ccut(ext2,ext,gen))
          scut(ext,ext2,gen)=max(scut(ext,ext2,gen),scut(ext2,ext,gen))
        enddo
        enddo

c calculating cut parameters to be given to the generator
        mtotcutgen=0d0
        do ext=3,nextgen(gen)
          mtotcutgen=mtotcutgen+sqrt(ptcut(ext,gen)**2+mass2(ext,gen))
          do ext2=1,nextgen(gen)
c ccutgen must be set to 1 if there are no cuts!
            ccutgen(ext2,ext,gen)=ccut(ext2,ext,gen)
            scutgen(ext2,ext,gen)=scut(ext2,ext,gen)
          enddo
        enddo
        mtotcutgen=max(mtotcutgen,mtotcut)

        end

***********************************************************************
      subroutine init_cuts_Hnl(gen)
***********************************************************************
*---------------------------------------------------------------------*
*     23.02.11  Stefan Kallweit   last changed  21.02.13              *
*---------------------------------------------------------------------*
*     !!! definitions are -not- used in actual cut routines !!!       *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none
      real*8 mtotcut

#include "mcparams.h"
#include "mccommon.h"

ccccccccccccccccccccccccccccccccccccccc
c  standard cuts for p p -> nu~ l- H  c
ccccccccccccccccccccccccccccccccccccccc
        mtotcut = 0d0
        if(scuts(gen).eq.0.and.shtr.eq.1)then
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(gen).eq.1.or.scuts(gen).eq.2)then
          ptcut(4,gen)=ptcutl
          ptcut(3,gen)=ptcutm
          ycut(4,gen)=ycutl
          if(shtr.ge.2)then
            ptcut(5,gen)=ptcuts
            ptcut(6,gen)=ptcuts
            ycut(5,gen)=ycuts
            ycut(6,gen)=ycuts
          endif
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(0).ne.0) then
          write(*,'(a30,i3)')' Unknown input: scuts = ',scuts(gen)
          stop
        endif

c filling cut parameters
c calculating cut parameters
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
          ccut(ext,ext2,gen)=dcos(ccut(ext,ext2,gen)*pi/180d0)
          scut(ext,ext2,gen)=scut(ext,ext2,gen)**2
        enddo
        enddo
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
          ccut(ext,ext2,gen)=min(ccut(ext,ext2,gen),ccut(ext2,ext,gen))
          scut(ext,ext2,gen)=max(scut(ext,ext2,gen),scut(ext2,ext,gen))
        enddo
        enddo

c calculating cut parameters to be given to the generator
        mtotcutgen=0d0
        do ext=3,nextgen(gen)
          mtotcutgen=mtotcutgen+sqrt(ptcut(ext,gen)**2+mass2(ext,gen))
          do ext2=1,nextgen(gen)
c ccutgen must be set to 1 if there are no cuts!
            ccutgen(ext2,ext,gen)=ccut(ext2,ext,gen)
            scutgen(ext2,ext,gen)=scut(ext2,ext,gen)
          enddo
        enddo
        mtotcutgen=max(mtotcutgen,mtotcut)
      end

***********************************************************************
      subroutine init_cuts_Hnn(gen)
***********************************************************************
*---------------------------------------------------------------------*
*     23.02.11  Stefan Kallweit   last changed  21.02.13              *
*---------------------------------------------------------------------*
*     !!! definitions are -not- used in actual cut routines !!!       *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none
      real*8 mtotcut

#include "mcparams.h"
#include "mccommon.h"

c standard cuts for pp -> lvH
        mtotcut = 0d0
        if(scuts(gen).eq.0.and.shtr.eq.1)then
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(gen).eq.1.or.scuts(gen).eq.2)then
c          ptcut(3,gen)=ptcutm
c          ptcut(4,gen)=ptcutm
          if(shtr.ge.2)then
            ptcut(5,gen)=ptcuts
            ptcut(6,gen)=ptcuts
            ycut(5,gen)=ycuts
            ycut(6,gen)=ycuts
          endif
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(0).ne.0) then
          write(*,'(a30,i3)')' Unknown input: scuts = ',scuts(gen)
          stop
        endif

c filling cut parameters
c calculating cut parameters
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
          ccut(ext,ext2,gen)=dcos(ccut(ext,ext2,gen)*pi/180d0)
          scut(ext,ext2,gen)=scut(ext,ext2,gen)**2
        enddo
        enddo
        do ext=1,nextgen(gen)
        do ext2=1,nextgen(gen)
          ccut(ext,ext2,gen)=min(ccut(ext,ext2,gen),ccut(ext2,ext,gen))
          scut(ext,ext2,gen)=max(scut(ext,ext2,gen),scut(ext2,ext,gen))
        enddo
        enddo

c calculating cut parameters to be given to the generator
        mtotcutgen=0d0
        do ext=3,nextgen(gen)
          mtotcutgen=mtotcutgen+sqrt(ptcut(ext,gen)**2+mass2(ext,gen))
          do ext2=1,nextgen(gen)
c ccutgen must be set to 1 if there are no cuts!
            ccutgen(ext2,ext,gen)=ccut(ext2,ext,gen)
            scutgen(ext2,ext,gen)=scut(ext2,ext,gen)
          enddo
        enddo
        mtotcutgen=max(mtotcutgen,mtotcut)
      end



***********************************************************************
      subroutine settings_Hll(kbeam,k,lweight,n,gen,step)
***********************************************************************
*     histogram routine for process p p -> l+ l- H (selprocess = 3)   *
*                                                                     *
*     filling of individual histograms                                *
*     definition of histograms done here                              *
*     kbeam: beam momenta                                             *
*     k    : momenta of event                                         *
*     lweight(maxon): weights of (sub)event                           *
*     n    : event number                                             *
*     gen  : generator                                                *
*     steps: flag for subroutine histogram                            *
*           steps=1: initialization                                   *
*           steps=2: start filling histograms                         *
*           steps=3: finish filling of histograms                     *
*           steps=4: print data for histograms                        *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     04.06.10  Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),lweight(maxcon),lenergy,x
      real*8 kt(maxe),y(maxe),ktH,ktl1,ktl2,yH
      real*8 denum,kh(0:3),kv(0:3),khvec
      complex*16 eiphl1,eiphl2,eiphH
      integer n,step,histo,switch
      integer lep1,lep2,jet1,jet2,nbin,count
      integer lastp
      character*30 hname
      character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,                                           
     &I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,                                             
     &X,1PE18.11,X,1F3.0)"
      character(len=*),parameter :: Fmt2 = "(I2,2X,I2,2X,1PE14.7,                                                  
     &2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)"

      save count
      data nbin /100/

      if(step.eq.2)then
         do icon=1,maxcon
            if (lweight(icon).ne.0d0) goto 100
         enddo
         return
      endif

 100  continue
      switch=0
      lenergy=kbeam(1,0)+kbeam(2,0)

ccccccccccccccccccccccccccccc
c  calculate momentum of H  c
ccccccccccccccccccccccccccccc
      if(shtr.eq.0)then
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)
         enddo
         lastp=6
      else
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)+k(6,vcomp)
         enddo
         lastp=7
      endif

ccccccccccccccccccccccccccccccccc
c  calculate momentum of V = Z  c
ccccccccccccccccccccccccccccccccc
      do vcomp=0,3
         kv(vcomp)=k(3,vcomp)+k(4,vcomp)
      enddo

c declaration: histogram(hname,y,x,min,max,parts,histo,n,addpoint,step)

      histo=0

      if(lnoutmc.ge.3)then
ccccccccccccccccccccccccc
c  total cross section  c
ccccccccccccccccccccccccc
      hname='dat.tot'
      call histogram(hname,lweight,lenergy+0.5d0,lenergy,lenergy+1d0,1,
     &     histo,n,0,gen,step)
      endif

      if(shtr.gt.0) then
ccccccccccccccccccccccccc
c invariant mass distribution of Higgs boson
ccccccccccccccccccccccccc
        histo=histo+1
        hname='dat.mh'
        x=kh(0)**2-kh(1)**2-kh(2)**2-kh(3)**2
        if(x.ge.0d0) then
          x=sqrt(x)
        else
          x=-1d0
        endif
        call histogram(hname,lweight,x,mcm(25)-2d0*mcw(25),
     &            mcm(25)+2d0*mcw(25),nbin,histo,n,1,gen,step)

        histo=histo+1
        hname='dat.mhwide'
        call histogram(hname,lweight,x,0d0,
     &            1d3,nbin,histo,n,1,gen,step)
      endif

cccccccccccccccccccccccccccccccccccccccc
c  transverse momentum of Higgs boson  c
cccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.pth'
      x=sqrt(kh(1)**2+kh(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccc
c  rapidity of Higgs boson  c
ccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.yh'
      if(kh(0)-kh(3).ne.0d0.and.kh(0)+kh(3).ne.0d0) then
        x=.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))
      else
        x=1d10
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)

cccccccccccccccccccccccccccccccccc
c pseudorapidity of Higgs boson  c
cccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etah'
      if(sqrt(kh(1)**2+kh(2)**2+kh(3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,kh(3)/sqrt(kh(1)**2+kh(2)**2+kh(3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccccccccccc
c  transverse momentum of V = Z  c
cccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptv'
      x=sqrt(kv(1)**2+kv(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccccc
c  transverse momentum of l+  c
ccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptlp'
      x=sqrt(k(3,1)**2+k(3,2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccc
c  rapidity of l+  c
cccccccccccccccccccc
      histo=histo+1
      hname='dat.ylp'
      if(k(3,0)-k(3,3).ne.0d0.and.k(3,0)+k(3,3).ne.0d0) then
         x=.5d0*log((k(3,0)+k(3,3))/(k(3,0)-k(3,3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccc
c  pseudorapidity of l+  c
cccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etalp'
      if(sqrt(k(3,1)**2+k(3,2)**2+k(3,3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,k(3,3)/sqrt(k(3,1)**2+k(3,2)**2+k(3,3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

ccccccccccccccccccccccccccccccc
c  transverse momentum of l-  c
ccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptlm'
      x=sqrt(k(4,1)**2+k(4,2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccc
c  rapidity of l-  c
cccccccccccccccccccc
      histo=histo+1
      hname='dat.ylm'
      if(k(4,0)-k(4,3).ne.0d0.and.k(4,0)+k(4,3).ne.0d0) then
         x=.5d0*log((k(4,0)+k(4,3))/(k(4,0)-k(4,3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccc
c  pseudorapidity of l-  c
cccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etalm'
      if(sqrt(k(4,1)**2+k(4,2)**2+k(4,3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,k(4,3)/sqrt(k(4,1)**2+k(4,2)**2+k(4,3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccccccccccccccccccccc
c  select tagging leptons with highest pt  c
cccccccccccccccccccccccccccccccccccccccccccc
      kt(3)=sqrt(k(3,1)**2+k(3,2)**2)
      kt(4)=sqrt(k(4,1)**2+k(4,2)**2)
      kt(lastp)=sqrt(k(lastp,1)**2+k(lastp,2)**2)
      if(kt(3).gt.kt(4).and.kt(3).gt.kt(lastp))then
         lep1=3
         if(kt(4).gt.kt(lastp))then
            lep2=4
         else
            lep2=lastp
         endif
      elseif(kt(4).gt.kt(3).and.kt(4).gt.kt(lastp))then
         lep1=4
         if(kt(3).gt.kt(lastp))then
            lep2=3
         else
            lep2=lastp
         endif
      else
         lep1=lastp
         if(kt(3).gt.kt(4))then
            lep2=3
         else
            lep2=4
         endif
      endif

      if(kt(3).gt.kt(4))then
         lep1=3
         lep2=4
      else
         lep1=4
         lep2=3
      endif



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in azimuthal-angle difference between H and l+   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.dphilpH'
      ktl1=sqrt(k(3,1)**2+k(3,2)**2)
      ktH =sqrt(kh(1)**2+kh(2)**2)
      yH =.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))

      if(ktl1*ktH.ne.0d0)then
c         x=(k(3,1)*kh(1)+k(3,2)*kh(2))/(ktl1*ktH)

	 eiphl1 = dcmplx(k(3,1),k(3,2))/kt(3)
	 eiphH  = dcmplx(kh(1),kh(2))/ktH
	 x = dimag(log(eiphl1/eiphH))*180d0/pi
         if (y(3).lt.yH) x=-x

      else
         x=1000d0
      endif
      call histogram(hname,lweight,x,-180d0,180d0,
     &	nbin,histo,n,1,gen,step)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in azimuthal-angle difference between H and l-   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.dphilmH'
      ktl2=sqrt(k(4,1)**2+k(4,2)**2)
      ktH =sqrt(kh(1)**2+kh(2)**2)
      yH =.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))

      if(ktl2*ktH.ne.0d0)then
c         x=(k(4,1)*kh(1)+k(4,2)*kh(2))/(ktl2*ktH)

	 eiphl2 = dcmplx(k(4,1),k(4,2))/kt(4)
	 eiphH  = dcmplx(kh(1),kh(2))/ktH
	 x = dimag(log(eiphl2/eiphH))*180d0/pi
         if (y(4).lt.yH) x=-x

      else
         x=1000d0
      endif

      call histogram(hname,lweight,x,-180d0,180d0,
     &	nbin,histo,n,1,gen,step)
c     Write event out                                                                                 
c      print *, lweight
c     print *, "dpne"
      
      IF ( ( lweight(1) .NE. 0) .and. (step .EQ. 2 ) ) THEN
      write(1988,"(A)") "<event>"
      write(1988,Fmt2) 6, 60, lweight(1) ,1.0 ,1.0 ,1.0

      write(1988,Fmt1)-1,-1,0,0,0,501,kbeam(1,1),kbeam(1,2),kbeam(1,3),
     &kbeam(1,0),0.0,0.0,0.0
      write(1988,Fmt1)3,-1,0,0,502,0,kbeam(2,1),kbeam(2,2),kbeam(2,3),
     &kbeam(2,0),0.0,0.0,0.0
      write(1988,Fmt1)25,1,1,2,0,0,kh(1),kh(2),kh(3),kh(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)23,1,1,2,0,0,kv(1),kv(2),kv(3),kv(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)11,1,4,4,0,0,k(3,1),k(3,2),k(3,3),
     &k(3,0),0.0,0.0,0.0
      write(1988,Fmt1)-11,1,4,4,0,0,k(4,1),k(4,2),k(4,3),
     &k(4,0),0.0,0.0,0.0

      write(1988,"(A)") "</event>"

      ELSE
      PRINT  *,"0Weight"

      ENDIF




      return

      end



***********************************************************************
      subroutine settings_Hln(kbeam,k,lweight,n,gen,step)
***********************************************************************
*     histogram routine for process p p -> l+ nu H (selprocess = 1)   *
*                                                                     *
*     filling of individual histograms                                *
*     definition of histograms done here                              *
*     kbeam: beam momenta                                             *
*     k    : momenta of event                                         *
*     lweight(maxon): weights of (sub)event                           *
*     n    : event number                                             *
*     gen  : generator                                                *
*     steps: flag for subroutine histogram                            *
*           steps=1: initialization                                   *
*           steps=2: start filling histograms                         *
*           steps=3: finish filling of histograms                     *
*           steps=4: print data for histograms                        *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     04.06.10  Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),lweight(maxcon),lenergy,x
      real*8 kt(maxe),y(maxe),ktH,ktl1,yH
      real*8 denum,kh(0:3),kv(0:3),khvec
      complex*16 eiphl1,eiphl2,eiphH
      integer n,step,histo,switch
      integer lep1,lep2,jet1,jet2,nbin,count
      integer lastp
      character*30 hname
      character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,                                           
     &I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,                                             
     &X,1PE18.11,X,1F3.0)"
      character(len=*),parameter :: Fmt2 = "(I2,2X,I2,2X,1PE14.7,                                                  
     &2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)"

      
      save count
      data nbin /100/

      if(step.eq.2)then
         do icon=1,maxcon
            if (lweight(icon).ne.0d0) goto 100
         enddo
         return
      endif

 100  continue
      switch=0
      lenergy=kbeam(1,0)+kbeam(2,0)

ccccccccccccccccccccccccccccc
c  calculate momentum of H  c
ccccccccccccccccccccccccccccc
      if(shtr.eq.0)then
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)
         enddo
         lastp=6
      else
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)+k(6,vcomp)
         enddo
         lastp=7
      endif

cccccccccccccccccccccccccccccccccc
c  calculate momentum of V = W+  c
cccccccccccccccccccccccccccccccccc
      do vcomp=0,3
         kv(vcomp)=k(3,vcomp)+k(4,vcomp)
      enddo

c declaration: histogram(hname,y,x,min,max,parts,histo,n,addpoint,step)

      histo=0

      if(lnoutmc.ge.3)then
ccccccccccccccccccccccccc
c  total cross section  c
ccccccccccccccccccccccccc
      hname='dat.tot'
      call histogram(hname,lweight,lenergy+0.5d0,lenergy,lenergy+1d0,1,
     &     histo,n,0,gen,step)
      endif

      if(shtr.gt.0) then
ccccccccccccccccccccccccc
c invariant mass distribution of Higgs boson
ccccccccccccccccccccccccc
        histo=histo+1
        hname='dat.mh'
        x=kh(0)**2-kh(1)**2-kh(2)**2-kh(3)**2
        if(x.ge.0d0) then
          x=sqrt(x)
        else
          x=-1d0
        endif
        call histogram(hname,lweight,x,mcm(25)-2d0*mcw(25),
     &            mcm(25)+2d0*mcw(25),nbin,histo,n,1,gen,step)

        histo=histo+1
        hname='dat.mhwide'
        call histogram(hname,lweight,x,0d0,
     &            1d3,nbin,histo,n,1,gen,step)
      endif

cccccccccccccccccccccccccccccccccccccccc
c  transverse momentum of Higgs boson  c
cccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.pth'
      x=sqrt(kh(1)**2+kh(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccc
c  rapidity of Higgs boson  c
ccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.yh'
      if(kh(0)-kh(3).ne.0d0.and.kh(0)+kh(3).ne.0d0) then
         x=.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccccccccc
c  pseudorapidity of Higgs boson  c
ccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etah'
      if(sqrt(kh(1)**2+kh(2)**2+kh(3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,kh(3)/sqrt(kh(1)**2+kh(2)**2+kh(3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,
     &     step)

ccccccccccccccccccccccccccccccccccc
c  transverse momentum of V = W+  c
ccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptv'
      x=sqrt(kv(1)**2+kv(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccccc
c  transverse momentum of l+  c
ccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptlp'
      x=sqrt(k(3,1)**2+k(3,2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccc
c  rapidity of l+  c
cccccccccccccccccccc
      histo=histo+1
      hname='dat.ylp'
      if(k(3,0)-k(3,3).ne.0d0.and.k(3,0)+k(3,3).ne.0d0) then
         x=.5d0*log((k(3,0)+k(3,3))/(k(3,0)-k(3,3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccc
c  pseudorapidity of l+  c
cccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etalp'
      if(sqrt(k(3,1)**2+k(3,2)**2+k(3,3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,k(3,3)/sqrt(k(3,1)**2+k(3,2)**2+k(3,3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

ccccccccccccccccccccccccccccccc
c  transverse momentum of nu  c
ccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptmiss'
      x=sqrt(k(4,1)**2+k(4,2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccc
c  rapidity of nu  c
cccccccccccccccccccc
      histo=histo+1
      hname='dat.ymiss'
      if(k(4,0)-k(4,3).ne.0d0.and.k(4,0)+k(4,3).ne.0d0) then
         x=.5d0*log((k(4,0)+k(4,3))/(k(4,0)-k(4,3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccc
c  pseudorapidity of nu  c
cccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etamiss'
      if(sqrt(k(4,1)**2+k(4,2)**2+k(4,3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,k(4,3)/sqrt(k(4,1)**2+k(4,2)**2+k(4,3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in azimuthal-angle difference between H and l+   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.dphilpH'
      ktl1=sqrt(k(3,1)**2+k(3,2)**2)
      kt(3)=sqrt(k(3,1)**2+k(3,2)**2)
      y(3)=.5d0*log((k(3,0)+k(3,3))/(k(3,0)-k(3,3)))
      ktH =sqrt(kh(1)**2+kh(2)**2)
      yH =.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))

      if(ktl1*ktH.ne.0d0)then
c         x=(k(3,1)*kh(1)+k(3,2)*kh(2))/(ktl1*ktH)
	 eiphl1 = dcmplx(k(3,1),k(3,2))/kt(3)
	 eiphH  = dcmplx(kh(1),kh(2))/ktH
	 x = dimag(log(eiphl1/eiphH))*180d0/pi
         if (y(3).lt.yH) x=-x

      else
         x=1000d0
      endif

      call histogram(hname,lweight,x,-180d0,180d0,
     &	nbin,histo,n,1,gen,step)




c     Write event out                                                                                 
c      print *, lweight
c      print *, "dpne"

      IF ( ( lweight(1) .NE. 0) .and. (step .EQ. 2 ) ) THEN    
      write(1988,"(A)") "<event>"
      write(1988,Fmt2) 6, 60, lweight(1) ,1.0 ,1.0 ,1.0

      write(1988,Fmt1)-1,-1,0,0,0,501,kbeam(1,1),kbeam(1,2),kbeam(1,3),
     &kbeam(1,0),0.0,0.0,0.0
      write(1988,Fmt1)3,-1,0,0,502,0,kbeam(2,1),kbeam(2,2),kbeam(2,3),
     &kbeam(2,0),0.0,0.0,0.0
      write(1988,Fmt1)25,1,1,2,0,0,kh(1),kh(2),kh(3),kh(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)24,1,1,2,0,0,kv(1),kv(2),kv(3),kv(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)11,1,4,4,0,0,k(3,1),k(3,2),k(3,3),
     &k(3,0),0.0,0.0,0.0
      write(1988,Fmt1)-12,1,4,4,0,0,k(4,1),k(4,2),k(4,3),
     &k(4,0),0.0,0.0,0.0

      write(1988,"(A)") "</event>"

      ELSE
      PRINT  *,"0Weight"

      ENDIF






      
      end



***********************************************************************
      subroutine settings_Hnl(kbeam,k,lweight,n,gen,step)
***********************************************************************
*     histogram routine for process p p -> nu~ l- H (selprocess = 2)  *
*                                                                     *
*     filling of individual histograms                                *
*     definition of histograms done here                              *
*     kbeam: beam momenta                                             *
*     k    : momenta of event                                         *
*     lweight(maxon): weights of (sub)event                           *
*     n    : event number                                             *
*     gen  : generator                                                *
*     steps: flag for subroutine histogram                            *
*           steps=1: initialization                                   *
*           steps=2: start filling histograms                         *
*           steps=3: finish filling of histograms                     *
*           steps=4: print data for histograms                        *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     04.06.10  Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),lweight(maxcon),lenergy,x
      real*8 kt(maxe),y(maxe),ktH,ktl2,yH
      real*8 denum,kh(0:3),kv(0:3),khvec
      complex*16 eiphl1,eiphl2,eiphH
      integer n,step,histo,switch
      integer jet1,jet2,nbin,count
      integer lastp
      character*30 hname
      character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,                                           
     &I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,                                             
     &X,1PE18.11,X,1F3.0)"
      character(len=*),parameter :: Fmt2 = "(I2,2X,I2,2X,1PE14.7,                                                  
     &2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)"

      save count
      data nbin /100/

      if(step.eq.2)then
         do icon=1,maxcon
            if (lweight(icon).ne.0d0) goto 100
         enddo
         return
      endif

 100  continue
      switch=0
      lenergy=kbeam(1,0)+kbeam(2,0)

ccccccccccccccccccccccccccccc
c  calculate momentum of H  c
ccccccccccccccccccccccccccccc
      if(shtr.eq.0)then
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)
         enddo
         lastp=6
      else
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)+k(6,vcomp)
         enddo
         lastp=7
      endif

cccccccccccccccccccccccccccccccccc
c  calculate momentum of V = W-  c
cccccccccccccccccccccccccccccccccc
      do vcomp=0,3
         kv(vcomp)=k(3,vcomp)+k(4,vcomp)
      enddo

c declaration: histogram(hname,y,x,min,max,parts,histo,n,addpoint,step)

      histo=0

      if(lnoutmc.ge.3)then
ccccccccccccccccccccccccc
c  total cross section  c
ccccccccccccccccccccccccc
      hname='dat.tot'
      call histogram(hname,lweight,lenergy+0.5d0,lenergy,lenergy+1d0,1,
     &     histo,n,0,gen,step)
      endif

      if(shtr.gt.0) then
ccccccccccccccccccccccccc
c invariant mass distribution of Higgs boson
ccccccccccccccccccccccccc
        histo=histo+1
        hname='dat.mh'
        x=kh(0)**2-kh(1)**2-kh(2)**2-kh(3)**2
        if(x.ge.0d0) then
          x=sqrt(x)
        else
          x=-1d0
        endif
        call histogram(hname,lweight,x,mcm(25)-2d0*mcw(25),
     &            mcm(25)+2d0*mcw(25),nbin,histo,n,1,gen,step)

        histo=histo+1
        hname='dat.mhwide'
        call histogram(hname,lweight,x,0d0,
     &            1d3,nbin,histo,n,1,gen,step)
      endif

cccccccccccccccccccccccccccccccccccccccc
c  transverse momentum of Higgs boson  c
cccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.pth'
      x=sqrt(kh(1)**2+kh(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccc
c  rapidity of Higgs boson  c
ccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.yh'
      if(kh(0)-kh(3).ne.0d0.and.kh(0)+kh(3).ne.0d0) then
         x=.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccccccccc
c  pseudorapidity of Higgs boson  c
ccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etah'
      if(sqrt(kh(1)**2+kh(2)**2+kh(3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,kh(3)/sqrt(kh(1)**2+kh(2)**2+kh(3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,
     &     step)

ccccccccccccccccccccccccccccccccccc
c  transverse momentum of V = W-  c
ccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptv'
      x=sqrt(kv(1)**2+kv(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccccccccccccccc
c  transverse momentum of nu~  c
cccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptmiss'
      x=sqrt(k(3,1)**2+k(3,2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccc
c  rapidity of nu~  c
cccccccccccccccccccc
      histo=histo+1
      hname='dat.ymiss'
      if(k(3,0)-k(3,3).ne.0d0.and.k(3,0)+k(3,3).ne.0d0) then
         x=.5d0*log((k(3,0)+k(3,3))/(k(3,0)-k(3,3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

ccccccccccccccccccccccccccc
c  pseudorapidity of nu~  c
ccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etamiss'
      if(sqrt(k(3,1)**2+k(3,2)**2+k(3,3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,k(3,3)/sqrt(k(3,1)**2+k(3,2)**2+k(3,3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

ccccccccccccccccccccccccccccccc
c  transverse momentum of l-  c
ccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptlm'
      x=sqrt(k(4,1)**2+k(4,2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

cccccccccccccccccccc
c  rapidity of l-  c
cccccccccccccccccccc
      histo=histo+1
      hname='dat.ylm'
      if(k(4,0)-k(4,3).ne.0d0.and.k(4,0)+k(4,3).ne.0d0) then
         x=.5d0*log((k(4,0)+k(4,3))/(k(4,0)-k(4,3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccc
c  pseudorapidity of l-  c
cccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etalm'
      if(sqrt(k(4,1)**2+k(4,2)**2+k(4,3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,k(4,3)/sqrt(k(4,1)**2+k(4,2)**2+k(4,3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in azimuthal-angle difference between H and l-   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.dphilmH'
      ktl2=sqrt(k(4,1)**2+k(4,2)**2)
      kt(4)=sqrt(k(4,1)**2+k(4,2)**2)
      y(4)=.5d0*log((k(4,0)+k(4,3))/(k(4,0)-k(4,3)))
      ktH =sqrt(kh(1)**2+kh(2)**2)
      yH =.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))

      if(ktl2*ktH.ne.0d0)then
c         x=(k(4,1)*kh(1)+k(4,2)*kh(2))/(ktl2*ktH)
	 eiphl2 = dcmplx(k(4,1),k(4,2))/kt(4)
	 eiphH  = dcmplx(kh(1),kh(2))/ktH
	 x = dimag(log(eiphl2/eiphH))*180d0/pi
         if (y(4).lt.yH) x=-x

      else
         x=1000d0
      endif

      call histogram(hname,lweight,x,-180d0,180d0,
     &	nbin,histo,n,1,gen,step)

c     Write event out                                                                                 
c      print *, lweight
c      print *, "dpne"

      
      IF ( ( lweight(1) .NE. 0) .and. (step .EQ. 2 ) ) THEN    

      write(1988,"(A)") "<event>"
      write(1988,Fmt2) 6, 60, lweight(1) ,1.0 ,1.0 ,1.0

      write(1988,Fmt1)-1,-1,0,0,0,501,kbeam(1,1),kbeam(1,2),kbeam(1,3),
     &kbeam(1,0),0.0,0.0,0.0
      write(1988,Fmt1)3,-1,0,0,502,0,kbeam(2,1),kbeam(2,2),kbeam(2,3),
     &kbeam(2,0),0.0,0.0,0.0
      write(1988,Fmt1)25,1,1,2,0,0,kh(1),kh(2),kh(3),kh(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)-24,1,1,2,0,0,kv(1),kv(2),kv(3),kv(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)12,1,4,4,0,0,k(3,1),k(3,2),k(3,3),
     &k(3,0),0.0,0.0,0.0
      write(1988,Fmt1)-11,1,4,4,0,0,k(4,1),k(4,2),k(4,3),
     &k(4,0),0.0,0.0,0.0

      write(1988,"(A)") "</event>"

      ELSE
      PRINT  *,"0Weight"

      ENDIF

      




      
      end

***********************************************************************
      subroutine settings_Hnn(kbeam,k,lweight,n,gen,step)
***********************************************************************
*     histogram routine for process p p -> nu~ nu H (selprocess = 4)  *
*                                                                     *
*     filling of individual histograms                                *
*     definition of histograms done here                              *
*     kbeam: beam momenta                                             *
*     k    : momenta of event                                         *
*     lweight(maxon): weights of (sub)event                           *
*     n    : event number                                             *
*     gen  : generator                                                *
*     steps: flag for subroutine histogram                            *
*           steps=1: initialization                                   *
*           steps=2: start filling histograms                         *
*           steps=3: finish filling of histograms                     *
*           steps=4: print data for histograms                        *
*     adapted from COFFERGG by Markus Roth                            *
*---------------------------------------------------------------------*
*     04.06.10  Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),lweight(maxcon),lenergy,x
      real*8 kmiss(0:3)
      real*8 kt(maxe),y(maxe)
      real*8 denum,kh(0:3),kv(0:3),khvec
      real*8 ktmiss
      integer n,step,histo,switch
      integer jet1,jet2,nbin,count
      integer lastp
      character*30 hname

      character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,                                           
     &I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,                                             
     &X,1PE18.11,X,1F3.0)"
      character(len=*),parameter :: Fmt2 = "(I2,2X,I2,2X,1PE14.7,                                                  
     &2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)"


      
      save count
      data nbin /100/

      if(step.eq.2)then
         do icon=1,maxcon
            if (lweight(icon).ne.0d0) goto 100
         enddo
         return
      endif

 100  continue
      switch=0
      lenergy=kbeam(1,0)+kbeam(2,0)

ccccccccccccccccccccccccccccc
c  calculate momentum of H  c
ccccccccccccccccccccccccccccc
      if(shtr.eq.0)then
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)
         enddo
         lastp=6
      else
         do vcomp=0,3
            kh(vcomp)=k(5,vcomp)+k(6,vcomp)
         enddo
         lastp=7
      endif

ccccccccccccccccccccccccccccccccc
c  calculate momentum of V = Z  c
ccccccccccccccccccccccccccccccccc
      do vcomp=0,3
         kv(vcomp)=k(3,vcomp)+k(4,vcomp)
      enddo

c declaration: histogram(hname,y,x,min,max,parts,histo,n,addpoint,step)

      histo=0

      if(lnoutmc.ge.3)then
ccccccccccccccccccccccccc
c  total cross section  c
ccccccccccccccccccccccccc
      hname='dat.tot'
      call histogram(hname,lweight,lenergy+0.5d0,lenergy,lenergy+1d0,1,
     &     histo,n,0,gen,step)
      endif

      if(shtr.gt.0) then
ccccccccccccccccccccccccc
c invariant mass distribution of Higgs boson
ccccccccccccccccccccccccc
        histo=histo+1
        hname='dat.mh'
        x=kh(0)**2-kh(1)**2-kh(2)**2-kh(3)**2
        if(x.ge.0d0) then
          x=sqrt(x)
        else
          x=-1d0
        endif
        call histogram(hname,lweight,x,mcm(25)-2d0*mcw(25),
     &            mcm(25)+2d0*mcw(25),nbin,histo,n,1,gen,step)

        histo=histo+1
        hname='dat.mhwide'
        call histogram(hname,lweight,x,0d0,
     &            1d3,nbin,histo,n,1,gen,step)
      endif

cccccccccccccccccccccccccccccccccccccccc
c  transverse momentum of Higgs boson  c
cccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.pth'
      x=sqrt(kh(1)**2+kh(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccc
c  rapidity of Higgs boson  c
ccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.yh'
      if(kh(0)-kh(3).ne.0d0.and.kh(0)+kh(3).ne.0d0) then
         x=.5d0*log((kh(0)+kh(3))/(kh(0)-kh(3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccccccccccccc
c  pseudorapidity of Higgs boson  c
ccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etah'
      if(sqrt(kh(1)**2+kh(2)**2+kh(3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,kh(3)/sqrt(kh(1)**2+kh(2)**2+kh(3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccccccccccc
c  transverse momentum of V = Z  c
cccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptv'
      x=sqrt(kv(1)**2+kv(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccc
c  calculate kmiss  c
ccccccccccccccccccccc
      do vcomp=0,3
         kmiss(vcomp)=k(3,vcomp)+k(4,vcomp)
      enddo

ccccccccccccccccccccccccccccccccc
c  missing transverse momentum  c
ccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptmiss'
      x=sqrt(kmiss(1)**2+kmiss(2)**2)
      call histogram(hname,lweight,x,0d0,5d2,nbin,histo,n,1,gen,step)

ccccccccccccccccccccccccc
c  missing rapidity of  c
ccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ymiss'
      if(kmiss(0)-kmiss(3).ne.0d0.and.kmiss(0)+kmiss(3).ne.0d0) then
         x=.5d0*log((kmiss(0)+kmiss(3))/(kmiss(0)-kmiss(3)))
      else
         x=1d10
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)

cccccccccccccccccccccccccccccc
c  pseudorapidity of nu nu~  c
cccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.etamiss'

c      print *,kh(1),k(1,1)," ",kh(2),k(1,2)," ",kh(3),k(1,3),"         "
c     &,kh(0),k(1,0)

      
      if(sqrt(kmiss(1)**2+kmiss(2)**2+kmiss(3)**2).ne.0d0) then
         x=min(1d0,max(-1d0,kmiss(3)/sqrt(kmiss(1)**2
     &        +kmiss(2)**2+kmiss(3)**2)))
         x=acos(x)
         if(x.ne.0d0)then
            x=-log(tan(x/2))
         else
            x=1d10
         endif
      endif
      call histogram(hname,lweight,x,-2.5d0,2.5d0,nbin,histo,n,1,gen,
     &     step)


c     Write event out                                                                                 
c      print *, lweight
c      print *, "dpne"

      
      IF ( ( lweight(1) .NE. 0) .and. (step .EQ. 2 ) ) THEN    
      write(1988,"(A)") "<event>"
      write(1988,Fmt2) 6, 60, lweight(1) ,1.0 ,1.0 ,1.0

      write(1988,Fmt1)-1,-1,0,0,0,501,kbeam(1,1),kbeam(1,2),kbeam(1,3),
     &kbeam(1,0),0.0,0.0,0.0
      write(1988,Fmt1)3,-1,0,0,502,0,kbeam(2,1),kbeam(2,2),kbeam(2,3),
     &kbeam(2,0),0.0,0.0,0.0
      write(1988,Fmt1)25,1,1,2,0,0,kh(1),kh(2),kh(3),kh(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)23,1,1,2,0,0,kv(1),kv(2),kv(3),kv(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)12,1,4,4,0,0,k(3,1),k(3,2),k(3,3),
     &k(3,0),0.0,0.0,0.0
      write(1988,Fmt1)-12,1,4,4,0,0,k(4,1),k(4,2),k(4,3),
     &k(4,0),0.0,0.0,0.0

      write(1988,"(A)") "</event>"

      ELSE
      PRINT  *,"0Weight"

      ENDIF








      
      end

