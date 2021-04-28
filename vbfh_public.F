#define PRIVATE
#define NEWCUTS
#define TEST
***********************************************************************
      subroutine vbfh_recombination(kbeam,k,krec,generator,lnext,cut,
     *  switch)
***********************************************************************
*     input:     including incoming partons!!!                        *
*     kpart:    momenta of partons                                    *
*     typepart: types of partons                                      *
*     npartin: number of outgoing partons                             *
*     typepart: types of particles                                    *
*     particle types:  1 are jets and merge to jets (quarks,gluons)   *
*     particle types: -1 are NO jets but merge to jets (photons)      *
*     particle types:  0 are NO jets and do not merge to jets         *
*                                                  (Higgs,leptons)    *
*     output:                                                         *
*     kjet:     momenta of jets                                       *
*     typejet:  types of jets          result from particle types     *
*       jet types:  1 are jets                                        *
*       jet types: -1 are photons                                     *
*       jet types:  0 are Higgs,leptons                               *
*     njet:     number of jets                                        *
*     Higgs has ijet=3                                                *
*---------------------------------------------------------------------*
*     21.02.07  Ansgar Denner     last changed  11.05.07              *
*---------------------------------------------------------------------*
*     renamed and shifted from 'public.F'                             *
*     09.03.11  Stefan Kallweit   last changed  09.03.11              *
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

      real*8  kjet(maxe,0:3)
      integer typepart(maxe),typejet(maxe),npartin,njet

      integer ii

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

c jet recombination

      if(srecomb(generator).eq.1.or.srecomb(generator).eq.2) then

c  900      format(6x,' p(',i1,',',i1,')    = ',d20.14)

         if(generator.eq.0.or.lnext.eq.5.or.shtr.gt.0.and.lnext.eq.6)
     &        then
            typepart(1)=0
            typepart(2)=0
            typepart(3)=1
            typepart(4)=1
            if(shtr.eq.0)then
               typepart(5)=0
               npartin=5
            else
               typepart(5)=0
               typepart(6)=0
               npartin=6
            endif

            call ktcluster(kbeam,k,typepart,kjet,typejet,
     &           npartin,njet,switch)

            call applyjetcutn(kbeam,kjet,typejet,krec,
     &           generator,njet,cut,switch)

         elseif(generator.eq.1) then
            typepart(1)=0
            typepart(2)=0
            typepart(3)=1
            typepart(4)=1
            if(shtr.eq.0)then
               typepart(5)=0
               typepart(6)=-1
               npartin=6
            else
               typepart(5)=0
               typepart(6)=0
               typepart(7)=-1
               npartin=7
            endif

            call ktcluster(kbeam,k,typepart,kjet,typejet,
     &           npartin,njet,switch)

            call applyjetcutn(kbeam,kjet,typejet,krec,
     &           generator,njet,cut,switch)

         elseif(generator.ge.2.and.generator.le.6) then
            typepart(1)=0
            typepart(2)=0
            typepart(3)=1
            typepart(4)=1
            if(shtr.eq.0)then
               typepart(5)=0
               typepart(6)=1
               npartin=6
            else
               typepart(5)=0
               typepart(6)=0
               typepart(7)=1
               npartin=7
            endif
            call ktcluster(kbeam,k,typepart,kjet,typejet,
     &           npartin,njet,switch)

            call applyjetcutn(kbeam,kjet,typejet,krec,
     &           generator,njet,cut,switch)

         else
            write(*,*) 'recombination: recombination not applicable'
            write(*,*) 'recombination',generator,lnext
            stop
         endif
      endif

      end


***********************************************************************
      subroutine ktcluster(kbeam,kpart,typepart,kjet,typejet,
     &    npartin,njet,switch)
***********************************************************************
*     run II kT clustering algorithm  (hep-ex/0005012)                *
*     input:     including incoming partons!!!                        *
*     kpart:    momenta of partons                                    *
*     typepart: types of partons                                      *
*     npartin: number of outgoing partons                             *
*     typepart: types of particles                                    *
*     particle types:  1 are jets and merge to jets (quarks,gluons)   *
*     particle types: -1 are NO jets but merge to jets (photons)      *
*     particle types:  0 are NO jets and do not merge to jets         *
*                                                  (Higgs,leptons)    *
*     output:                                                         *
*     kjet:     momenta of jets                                       *
*     typejet:  types of jets          result from particle types     *
*       jet types:  1 are jets                                        *
*       jet types: -1 are photons                                     *
*       jet types:  0 are Higgs,leptons                               *
*     njet:     number of jets                                        *
*     Higgs has ijet=3                                                *
*---------------------------------------------------------------------*
*     21.02.07  Ansgar Denner     last changed  11.05.07              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

      real*8  kbeam(2,0:3),kpart(maxe,0:3),kjet(maxe,0:3)
      integer typepart(maxe),typejet(maxe)
      integer npartin,njet,switch

c local variables

      real*8  kpartl(maxe,0:3),kt2(maxe),yk(maxe),dphi
      real*8  ktiktj,dmin,d2(maxe,maxe),d1(maxe)
      integer n1min,n2min,npart,part,part2,mergestep,steps
      integer pointer(maxe)

 900  format(6x,' p(',i1,',',i1,')    = ',d20.14)

      if(switch.ne.0)return

c do not touch incoming partons
      do vcomp=0,3
        kjet(1,vcomp)=kpart(1,vcomp)
        kjet(2,vcomp)=kpart(2,vcomp)
      enddo

      njet=2
      part=0
      do ext=3,npartin
c extract nonclustering partons
        if(typepart(ext).eq.0)then
          njet=njet+1
          typejet(njet)=0
          do vcomp=0,3
            kjet(njet,vcomp)=kpart(ext,vcomp)
          enddo
        else
          if(kpart(ext,0)+kpart(ext,3).ne.0d0.and.
     &        kpart(ext,0)-kpart(ext,3).ne.0d0)then
            yk(ext)=.5d0*log((kpart(ext,0)+kpart(ext,3))/
     &          (kpart(ext,0)-kpart(ext,3)))
            if(yk(ext).lt.-10d0)then
              yk(ext)=-10d0
            elseif(yk(ext).gt.10d0)then
              yk(ext)=10d0
            endif
          else
            yk(ext)=1d10
          endif
c take only clustering partons within acceptance
          if(dabs(yk(ext)).lt.partetacut) then
            part=part+1
            pointer(part)=ext
            do vcomp=0,3
              kpartl(ext,vcomp)=kpart(ext,vcomp)
            enddo
          endif
        endif
      enddo
      npart=part



c calculate cluster parameters
      do part=1,npart
        ext=pointer(part)
        kt2(ext)=kpartl(ext,1)**2+kpartl(ext,2)**2
        d1(ext)=kt2(ext)**ktpower


        yk(ext)=.5d0*log((kpartl(ext,0)+kpartl(ext,3))/
     &      (kpartl(ext,0)-kpartl(ext,3)))
        if(yk(ext).lt.-10d0)then
          yk(ext)=-10d0
        elseif(yk(ext).gt.10d0)then
          yk(ext)=10d0
        endif
        do part2=1,part-1
          ext2=pointer(part2)
          ktiktj=(kpartl(ext,1)*kpartl(ext2,1)
     &        +kpartl(ext,2)*kpartl(ext2,2))/
     &        sqrt(kt2(ext)*kt2(ext2))


          if(abs(ktiktj).gt.1d0)then
            switch=45
            if(numout.lt.maxout.and.lnoutmc.ge.5)then
 100          format(/" event = ",i9,", gen = ",i4,
     &            ", channel = ",i4,", exit = ",i4,(a))
              write(nout,100) nevt,cogen,cocha,switch,
     &            ', ktcluster:  |cos(deltaphi)| > 1'
              numout=numout+1
            endif
            return
          endif
          dphi=acos(ktiktj)
c          d2(ext2,ext)=min(kt2(ext),kt2(ext2))*
c     &        ((yk(ext)-yk(ext2))**2+dphi**2)/ktd**2
          d2(ext2,ext)=min(d1(ext),d1(ext2))*
     &        ((yk(ext)-yk(ext2))**2+dphi**2)/ktd**2
          d2(ext,ext2)=d2(ext2,ext)

        enddo
      enddo



c cluster remaining partons
      steps=npart
      do mergestep=1,steps


c determine dmin
c        dmin=kt2(pointer(1))
        dmin=d1(pointer(1))
        n1min=1
        n2min=0
        do part=1,npart
          ext=pointer(part)
          if(d1(ext).lt.dmin)then
            dmin=d1(ext)
c          if(kt2(ext).lt.dmin)then
c            dmin=kt2(ext)
            n1min=part
            n2min=0
          endif
          do part2=1,part-1
            ext2=pointer(part2)
            if(d2(ext2,ext).lt.dmin.and.
     &          (typepart(ext).gt.0.or.typepart(ext2).gt.0))then
              dmin=d2(ext2,ext)
              n1min=part
              n2min=part2
            endif
          enddo
        enddo

        if(npart.gt.1)then
          if (n2min.eq.0) then
c move parton to jets
            njet=njet+1
            ext=pointer(n1min)
            typejet(njet)=typepart(ext)
            do vcomp=0,3
              kjet(njet,vcomp)=kpartl(ext,vcomp)
            enddo
c shift maximal value of pointer to empty slot
            if(n1min.ne.npart) pointer(n1min)=pointer(npart)
            npart=npart-1
          else
c join two partons
            ext=pointer(n1min)
            ext2=pointer(n2min)
            if(typepart(ext).ne.1.and.typepart(ext2).ne.1) then
              write(*,*) 'ktcluster: error2',ext,ext2,
     &            typepart(ext),typepart(ext2)
              stop
            endif
            typepart(ext)=1
            do vcomp=0,3
              kpartl(ext,vcomp)=kpartl(ext,vcomp)+kpartl(ext2,vcomp)
            enddo


c shift maximal value of pointer to empty slot
            if(n2min.ne.npart) pointer(n2min)=pointer(npart)
            if(n1min.eq.npart) n1min=n2min
            npart=npart-1
c recalculate parameters of new parton
            kt2(ext)=kpartl(ext,1)**2+kpartl(ext,2)**2
            d1(ext)=kt2(ext)**ktpower
            yk(ext)=.5d0*log((kpartl(ext,0)+kpartl(ext,3))/
     &          (kpartl(ext,0)-kpartl(ext,3)))

            do part2=1,npart
              if(part2.ne.n1min)then
                ext2=pointer(part2)
                ktiktj=(kpartl(ext,1)*kpartl(ext2,1)
     &              +kpartl(ext,2)*kpartl(ext2,2))/
     &              sqrt(kt2(ext)*kt2(ext2))
                if(abs(ktiktj).gt.1d0)then
                  switch=45
                  if(numout.lt.maxout.and.lnoutmc.ge.5)then
                    write(nout,100) nevt,cogen,cocha,switch,
     &                  ', ktcluster:  |cos(deltaphi)| > 1'
                    numout=numout+1
                  endif
                  return
                endif
                dphi=acos(ktiktj)
c                d2(ext,ext2)=min(kt2(ext),kt2(ext2))*
c     &              ((yk(ext)-yk(ext2))**2+dphi**2)/ktd**2
                d2(ext,ext2)=min(d1(ext),d1(ext2))*
     &              ((yk(ext)-yk(ext2))**2+dphi**2)/ktd**2
                 d2(ext2,ext)=d2(ext,ext2)

              endif

            enddo

          endif
        else
          njet=njet+1
          ext=pointer(1)
          typejet(njet)=typepart(ext)
          do vcomp=0,3
            kjet(njet,vcomp)=kpartl(ext,vcomp)
          enddo
        endif
      enddo

      do ext=njet+1,maxe
        do vcomp=0,3
          kjet(ext,vcomp)=0d0
        enddo
      enddo

      end


***********************************************************************
      subroutine applyjetcutn(kbeam,kjet,typejet,kevent,
     &    generator,njet,cut,switch)
***********************************************************************
*     applies cuts to jets and Higgs after jet algorithm              *
*     input:                                                          *
*     kjet:     momenta of jets  (including incoming partons)         *
*     typejet:  types of jets                                         *
*     njet:     number of jets                                        *
*     output:                                                         *
*     kevent:   momenta of event                                      *
*     cut:      number of relevant cut                                *
*---------------------------------------------------------------------*
*     11.10.13  Ansgar Denner     last changed  13.10.13              *
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

      real*8  kbeam(2,0:3),kjet(maxe,0:3),kevent(maxe,0:3)
      integer typejet(maxe)
      integer generator,njet,cut,switch

c local variables
      real*8  kt2(maxe),y(maxe),swap,sjj,cj1h,cj2h,dphij1h,dphij2h
      real*8  ptvis2,sjjh
      integer nqjet,nsjet,pointer(maxe)

      if(switch.ne.0)return

      if(typejet(3).ne.0.or.shtr.gt.0.and.typejet(4).ne.0)then
        write(*,*) 'applyjetcut: Higgs is not first jet'
        stop
      endif

      if(shtr.eq.0)then

      if(scuts(generator).eq.1.or.scuts(generator).eq.2)then

c calculate kt2 and y of jets
        do ext=3,njet
          kt2(ext)=kjet(ext,1)**2+kjet(ext,2)**2
          if(kjet(ext,0)-kjet(ext,3).eq.0) then
            y(ext) = 1d10
          elseif(kjet(ext,0)+kjet(ext,3).eq.0) then
            y(ext) = -1d10
          else
            y(ext)=.5d0*log((kjet(ext,0)+kjet(ext,3))
     &          /(kjet(ext,0)-kjet(ext,3)))
          endif
        enddo

        do ext=njet+1,6
          kt2(ext)=0d0
          y(ext) = 1d10
        enddo

c initialize pointers
        pointer(3)=0
        pointer(4)=0
        pointer(5)=0
        pointer(6)=0

c set pointers to specific jet types
        nqjet=0
        do ext=3,njet
          if(typejet(ext).eq.0) pointer(5)=ext ! Higgs
          if(typejet(ext).eq.-1) pointer(6)=ext ! Photon
          if(typejet(ext).eq.1) then
            nqjet=nqjet+1
            if(nqjet.eq.1)then
              pointer(3)=ext    ! Jet1
            elseif(nqjet.eq.2)then
              pointer(4)=ext    ! Jet2
            elseif(nqjet.eq.3)then
              pointer(6)=ext    ! Jet3
            endif
          endif
        enddo
        if(nqjet.eq.1) then
          pointer(4)=njet+1
        elseif(nqjet.eq.0) then
          pointer(3)=njet+1
          pointer(4)=njet+2
        endif


c          write(*,*) 'nqjet=',nqjet
c          write(*,*) 'pointer a=',pointer


c order quark/gluon jets
        if(scuts(generator).eq.1) then
          if(nqjet.eq.2)then
            if(kt2(pointer(3)).lt.kt2(pointer(4)))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          elseif(nqjet.eq.3)then
            if(kt2(pointer(3)).lt.kt2(pointer(4)))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
            if(kt2(pointer(4)).lt.kt2(pointer(6)))then
              swap=pointer(4)
              pointer(4)=pointer(6)
              pointer(6)=swap
            endif
            if(kt2(pointer(3)).lt.kt2(pointer(4)))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          endif
        elseif(scuts(generator).eq.2) then
          if(nqjet.eq.2)then
            if(kjet(pointer(3),0).lt.kjet(pointer(4),0))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          elseif(nqjet.eq.3)then
            if(kjet(pointer(3),0).lt.kjet(pointer(4),0))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
            if(kjet(pointer(4),0).lt.kjet(pointer(6),0))then
              swap=pointer(4)
              pointer(4)=pointer(6)
              pointer(6)=swap
            endif
            if(kjet(pointer(3),0).lt.kjet(pointer(4),0))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          endif
        endif

c reorder jets if leading jets outside cuts
        if(abs(y(pointer(3))).gt.ycutj1.or.
     &      kt2(pointer(3)).lt.ptcutj1**2)then
          swap=pointer(3)
          pointer(3)=pointer(4)
          pointer(4)=swap
          if(nqjet.eq.3.and.(abs(y(pointer(3))).gt.ycutj1.or.
     &        kt2(pointer(3)).lt.ptcutj1**2))then
            swap=pointer(3)
            pointer(3)=pointer(6)
            pointer(6)=swap
          endif
        endif

        if(nqjet.eq.3.and.(abs(y(pointer(4))).gt.ycutj2.or.
     &      kt2(pointer(4)).lt.ptcutj2**2))then
          swap=pointer(4)
          pointer(4)=pointer(6)
          pointer(6)=swap
        endif

c determine momenta of event
        if(nqjet.eq.0)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(3,vcomp)
            kevent(3,vcomp)=0d0
            kevent(4,vcomp)=0d0
            kevent(6,vcomp)=0d0
          enddo
        elseif(nqjet.eq.1)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(3,vcomp)
            kevent(3,vcomp)=kjet(pointer(3),vcomp)
            kevent(4,vcomp)=0d0
            kevent(6,vcomp)=0d0
          enddo
        elseif(nqjet.eq.2)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(3,vcomp)
            kevent(3,vcomp)=kjet(pointer(3),vcomp)
            kevent(4,vcomp)=kjet(pointer(4),vcomp)
            kevent(6,vcomp)=0d0
          enddo
        elseif(nqjet.eq.3)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(3,vcomp)
            kevent(3,vcomp)=kjet(pointer(3),vcomp)
            kevent(4,vcomp)=kjet(pointer(4),vcomp)
            kevent(6,vcomp)=kjet(pointer(6),vcomp)
          enddo
        else
          write(*,*) 'applyjetcutn: not defined nqjet=',nqjet
        endif

c apply cuts

        cut=11
        if(abs(y(pointer(3))).gt.ycutj1) return
        cut=12
        if(abs(y(pointer(4))).gt.ycutj2) return

        cut=13
        if(kt2(pointer(3)).lt.ptcutj1**2) return
        cut=14
        if(kt2(pointer(4)).lt.ptcutj2**2) return

        cut=15
        if(nqjet.eq.1) then
          if (kt2(pointer(3)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(3))).lt.yminj1) return
        elseif(nqjet.eq.2) then
          if((kt2(pointer(3)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(3))).lt.yminj1).or.
     &      (kt2(pointer(4)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(4))).lt.yminj1)) return
        elseif(nqjet.eq.3) then
          if( (kt2(pointer(3)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(3))).lt.yminj1).or.
     &      (kt2(pointer(4)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(4))).lt.yminj1).or.
     &      (kt2(pointer(6)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(6))).lt.yminj1)) return
        endif

        cut=16
        if(nqjet.eq.2) then
          if (kt2(pointer(4)).gt.ptmaxj2**2
     &      .and.abs(y(pointer(4))).lt.yminj2) return
        elseif(nqjet.eq.3) then
          if((kt2(pointer(4)).gt.ptmaxj2**2
     &      .and.abs(y(pointer(4))).lt.yminj2).or.
     &      (kt2(pointer(6)).gt.ptmaxj2**2
     &      .and.abs(y(pointer(6))).lt.yminj2)) return
        endif

        cut=17
        if(abs(y(pointer(3))-y(pointer(4))).lt.dycutjj) return

        cut=18
        if(lhemicut.and.y(pointer(3))*y(pointer(4)).gt.0d0) return

        sjj=(kjet(pointer(3),0)+kjet(pointer(4),0))**2
     &      - (kjet(pointer(3),1)+kjet(pointer(4),1))**2
     &      - (kjet(pointer(3),2)+kjet(pointer(4),2))**2
     &      - (kjet(pointer(3),3)+kjet(pointer(4),3))**2

        cut=19
        if(mlcutjj.gt.0d0.and.sjj.lt.mlcutjj*mlcutjj)return

        cut=20
        if(sjj.gt.mucutjj*mucutjj)return

        if(kt2(pointer(3))*kt2(pointer(5)).ne.0d0)then
          cj1h= (kjet(pointer(3),1)*kjet(pointer(5),1)
     &          +kjet(pointer(3),2)*kjet(pointer(5),2))/
     &        sqrt(kt2(pointer(3))*kt2(pointer(5)))

          dphij1h=acos(cj1h)   !   *180d0/pi

          cut=21
          if(dphij1h.lt.dphicutj1h) return
        endif

        if(kt2(pointer(4))*kt2(pointer(5)).ne.0d0)then
          cj2h= (kjet(pointer(4),1)*kjet(pointer(5),1)
     &          +kjet(pointer(4),2)*kjet(pointer(5),2))/
     &        sqrt(kt2(pointer(4))*kt2(pointer(5)))

          cut=22
          dphij2h=acos(cj2h)   !   *180d0/pi

          if(dphij2h.lt.dphicutj2h) return
        endif

        ptvis2=(kjet(pointer(5),1)
     &      +kjet(pointer(3),1)+kjet(pointer(4),1))**2
     &      + (kjet(pointer(5),2)
     &      +kjet(pointer(3),2)+kjet(pointer(4),2))**2

        cut=23
        if(ptvis2.lt.ptcutvis**2) return

        cut=24
        if(ptvis2.gt.ptmaxvis**2) return

        cut=25
        if(nqjet.eq.3) then
          if  (kt2(pointer(6)).gt.ptmaxj3**2
     &      .and.abs(y(pointer(6))).lt.yminj3) return
        endif


        sjjh=(kjet(pointer(3),0)+kjet(pointer(4),0)
     &      +kjet(pointer(5),0))**2
     &      - (kjet(pointer(3),1)+kjet(pointer(4),1)
     &      +kjet(pointer(5),1))**2
     &      - (kjet(pointer(3),2)+kjet(pointer(4),2)
     &      +kjet(pointer(5),2))**2
     &      - (kjet(pointer(3),3)+kjet(pointer(4),3)
     &      +kjet(pointer(5),3))**2

        cut=26
        if(mlcutjjh.gt.0d0.and.sjjh.lt.mlcutjjh*mlcutjjh)return

        cut=27
        if(sjjh.gt.mucutjjh*mucutjjh)return


        cut=0

        call applyhiggscut(kbeam,kevent,generator,njet,cut,switch)

        endif

      else   ! if(shtr.gt.0)

      if(scuts(generator).eq.1.or.scuts(generator).eq.2)then


c calculate kt2 and y of jets
        do ext=3,njet
          kt2(ext)=kjet(ext,1)**2+kjet(ext,2)**2
          if(kjet(ext,0)-kjet(ext,3).eq.0) then
            y(ext) = 1d10
          elseif(kjet(ext,0)+kjet(ext,3).eq.0) then
            y(ext) = -1d10
          else
            y(ext)=.5d0*log((kjet(ext,0)+kjet(ext,3))
     &          /(kjet(ext,0)-kjet(ext,3)))
          endif
        enddo

        do ext=njet+1,7
          kt2(ext)=0d0
          y(ext) = 1d10
        enddo

c initialize pointers
        pointer(3)=0
        pointer(4)=0
        pointer(5)=0
        pointer(6)=0
        pointer(7)=0

c set pointers to specific jet types
        nqjet=0
        nsjet=0
        do ext=3,njet
          if(typejet(ext).eq.0) then
            nsjet=nsjet+1
            if(nsjet.eq.1)then
              pointer(5)=ext    ! singlet1
            elseif(nsjet.eq.2)then
              pointer(6)=ext    ! singlet2
            endif
          endif
          if(typejet(ext).eq.-1) pointer(7)=ext ! Photon
          if(typejet(ext).eq.1) then
            nqjet=nqjet+1
            if(nqjet.eq.1)then
              pointer(3)=ext    ! Jet1
            elseif(nqjet.eq.2)then
              pointer(4)=ext    ! Jet2
            elseif(nqjet.eq.3)then
              pointer(7)=ext    ! Jet3
            endif
          endif
        enddo
        if(nqjet.eq.1) then
          pointer(4)=njet+1
        elseif(nqjet.eq.0) then
          pointer(3)=njet+1
          pointer(4)=njet+2
        endif


c order quark/gluon jets
        if(scuts(generator).eq.1) then
          if(nqjet.eq.2)then
            if(kt2(pointer(3)).lt.kt2(pointer(4)))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          elseif(nqjet.eq.3)then
            if(kt2(pointer(3)).lt.kt2(pointer(4)))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
            if(kt2(pointer(4)).lt.kt2(pointer(7)))then
              swap=pointer(4)
              pointer(4)=pointer(7)
              pointer(7)=swap
            endif
            if(kt2(pointer(3)).lt.kt2(pointer(4)))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          endif
        elseif(scuts(generator).eq.2) then
          if(nqjet.eq.2)then
            if(kjet(pointer(3),0).lt.kjet(pointer(4),0))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          elseif(nqjet.eq.3)then
            if(kjet(pointer(3),0).lt.kjet(pointer(4),0))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
            if(kjet(pointer(4),0).lt.kjet(pointer(7),0))then
              swap=pointer(4)
              pointer(4)=pointer(7)
              pointer(7)=swap
            endif
            if(kjet(pointer(3),0).lt.kjet(pointer(4),0))then
              swap=pointer(3)
              pointer(3)=pointer(4)
              pointer(4)=swap
            endif
          endif
        endif

c reorder jets if leading jets outside cuts
        if(abs(y(pointer(3))).gt.ycutj1.or.
     &      kt2(pointer(3)).lt.ptcutj1**2)then
          swap=pointer(3)
          pointer(3)=pointer(4)
          pointer(4)=swap
          if(nqjet.eq.3.and.(abs(y(pointer(3))).gt.ycutj1.or.
     &        kt2(pointer(3)).lt.ptcutj1**2))then
            swap=pointer(3)
            pointer(3)=pointer(7)
            pointer(7)=swap
          endif
        endif

        if(nqjet.eq.3)then
        if(abs(y(pointer(4))).gt.ycutj2.or.
     &      kt2(pointer(4)).lt.ptcutj2**2)then
          swap=pointer(4)
          pointer(4)=pointer(7)
          pointer(7)=swap
        endif
        endif

c apply cuts


        cut=11
        if(abs(y(pointer(3))).gt.ycutj1) return
        cut=12
        if(abs(y(pointer(4))).gt.ycutj2) return

        cut=13
        if(kt2(pointer(3)).lt.ptcutj1**2) return
        cut=14
        if(kt2(pointer(4)).lt.ptcutj2**2) return

        cut=15
        if(nqjet.eq.1.and.
     &      (kt2(pointer(3)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(3))).lt.yminj1)
     &      .or.nqjet.eq.2.and.(
     &      (kt2(pointer(3)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(3))).lt.yminj1).or.
     &      (kt2(pointer(4)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(4))).lt.yminj1))
     &      .or.nqjet.eq.3.and.(
     &      (kt2(pointer(3)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(3))).lt.yminj1).or.
     &      (kt2(pointer(4)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(4))).lt.yminj1).or.
     &      (kt2(pointer(7)).gt.ptmaxj1**2
     &      .and.abs(y(pointer(7))).lt.yminj1))
     &      ) return

        cut=16
        if(nqjet.eq.2.and.
     &      (kt2(pointer(4)).gt.ptmaxj2**2
     &      .and.abs(y(pointer(4))).lt.yminj2)
     &      .or.nqjet.eq.3.and.
     &      ((kt2(pointer(4)).gt.ptmaxj2**2
     &      .and.abs(y(pointer(4))).lt.yminj2).or.
     &      (kt2(pointer(7)).gt.ptmaxj2**2
     &      .and.abs(y(pointer(7))).lt.yminj2))
     &      ) return

        cut=17
        if(abs(y(pointer(3))-y(pointer(4))).lt.dycutjj) return

        cut=18
        if(lhemicut.and.y(pointer(3))*y(pointer(4)).gt.0d0) return

        sjj=(kjet(pointer(3),0)+kjet(pointer(4),0))**2
     &      - (kjet(pointer(3),1)+kjet(pointer(4),1))**2
     &      - (kjet(pointer(3),2)+kjet(pointer(4),2))**2
     &      - (kjet(pointer(3),3)+kjet(pointer(4),3))**2

        cut=19
        if(mlcutjj.gt.0d0.and.sjj.lt.mlcutjj*mlcutjj)return

        cut=20
        if(sjj.gt.mucutjj*mucutjj)return

        if(kt2(pointer(3))*((kjet(pointer(5),1)+kjet(pointer(6),1))**2
     &       +(kjet(pointer(5),2)+kjet(pointer(6),2))**2).ne.0d0)then
          cj1h=
     &      (kjet(pointer(3),1)*(kjet(pointer(5),1)+kjet(pointer(6),1))
     &      +kjet(pointer(3),2)*(kjet(pointer(5),2)+kjet(pointer(6),2)))
     &      /sqrt(kt2(pointer(3))
     &      *((kjet(pointer(5),1)+kjet(pointer(6),1))**2
     &       +(kjet(pointer(5),2)+kjet(pointer(6),2))**2))

          dphij1h=acos(cj1h)  !    *180d0/pi

          cut=21
          if(dphij1h.lt.dphicutj1h) return
        endif

        if(kt2(pointer(4))*((kjet(pointer(5),1)+kjet(pointer(6),1))**2
     &       +(kjet(pointer(5),2)+kjet(pointer(6),2))**2).ne.0d0)then
          cj2h=
     &       (kjet(pointer(4),1)*(kjet(pointer(5),1)+kjet(pointer(6),1))
     &      +kjet(pointer(4),2)*(kjet(pointer(5),2)+kjet(pointer(6),2)))
     &        /sqrt(kt2(pointer(4))
     &        *((kjet(pointer(5),1)+kjet(pointer(6),1))**2
     &         +(kjet(pointer(5),2)+kjet(pointer(6),2))**2))

          dphij2h=acos(cj2h)   !    *180d0/pi

          cut=22
          if(dphij2h.lt.dphicutj2h) return
        endif	

        ptvis2=(kjet(pointer(5),1)+kjet(pointer(6),1)
     &      +kjet(pointer(3),1)+kjet(pointer(4),1))**2
     &      + (kjet(pointer(5),2)+kjet(pointer(6),2)
     &      +kjet(pointer(3),2)+kjet(pointer(4),2))**2

        cut=23
        if(ptvis2.lt.ptcutvis**2) return

        cut=24
        if(ptvis2.gt.ptmaxvis**2) return

        cut=25
        if(nqjet.eq.3.and.
     &      (kt2(pointer(7)).gt.ptmaxj3**2
     &      .and.abs(y(pointer(7))).lt.yminj3)
     &      ) return

        sjjh=(kjet(pointer(3),0)+kjet(pointer(4),0)
     &      +kjet(pointer(5),0)+kjet(pointer(6),0))**2
     &      - (kjet(pointer(3),1)+kjet(pointer(4),1)
     &      +kjet(pointer(5),1)+kjet(pointer(6),1))**2
     &      - (kjet(pointer(3),2)+kjet(pointer(4),2)
     &      +kjet(pointer(5),2)+kjet(pointer(6),2))**2
     &      - (kjet(pointer(3),3)+kjet(pointer(4),3)
     &      +kjet(pointer(5),3)+kjet(pointer(6),3))**2

        cut=26
        if(mlcutjjh.gt.0d0.and.sjjh.lt.mlcutjjh*mlcutjjh)return

        cut=27
        if(sjjh.gt.mucutjjh*mucutjjh)return


        cut=0


c determine momenta of event
        if(nqjet.eq.0)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(pointer(5),vcomp)
            kevent(6,vcomp)=kjet(pointer(6),vcomp)
            kevent(3,vcomp)=0d0
            kevent(4,vcomp)=0d0
            kevent(7,vcomp)=0d0
          enddo
        elseif(nqjet.eq.1)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(pointer(5),vcomp)
            kevent(6,vcomp)=kjet(pointer(6),vcomp)
            kevent(3,vcomp)=kjet(pointer(3),vcomp)
            kevent(4,vcomp)=0d0
            kevent(7,vcomp)=0d0
          enddo
        elseif(nqjet.eq.2)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(pointer(5),vcomp)
            kevent(6,vcomp)=kjet(pointer(6),vcomp)
            kevent(5,vcomp)=kjet(3,vcomp)
            kevent(3,vcomp)=kjet(pointer(3),vcomp)
            kevent(4,vcomp)=kjet(pointer(4),vcomp)
            kevent(7,vcomp)=0d0
          enddo
        elseif(nqjet.eq.3)then
          do vcomp=0,3
            kevent(1,vcomp)=kjet(1,vcomp)
            kevent(2,vcomp)=kjet(2,vcomp)
            kevent(5,vcomp)=kjet(pointer(5),vcomp)
            kevent(6,vcomp)=kjet(pointer(6),vcomp)
            kevent(3,vcomp)=kjet(pointer(3),vcomp)
            kevent(4,vcomp)=kjet(pointer(4),vcomp)
            kevent(7,vcomp)=kjet(pointer(7),vcomp)
          enddo
        else
          write(*,*) 'applyjetcutn: not defined nqjet=',nqjet
        endif

        call applyhiggscut(kbeam,kevent,generator,njet,cut,switch)

      endif
      endif                     ! if shtr.eq.0
      end

***********************************************************************
      subroutine init_cuts_Hjj(gen)
***********************************************************************
*---------------------------------------------------------------------*
*     renamed and shifted from 'public.F'                             *
*     04.06.10  Stefan Kallweit   last changed  09.06.10              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none
      real*8 mtotcut

#include "mcparams.h"
#include "mccommon.h"

c      print *, 'standard cuts for pp -> qqH'

c standard cuts for pp -> qqH
        if(scuts(gen).eq.0.and.shtr.eq.1)then
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth
        elseif(scuts(gen).eq.1.or.scuts(gen).eq.2)then
          ptcut(3,gen)=ptcutj1
          ptcut(4,gen)=ptcutj2
          ycut(3,gen)=ycutj1
          ycut(4,gen)=ycutj2
          dycut(3,4,gen)=dycutjj
          dycut(4,3,gen)=dycutjj
          scut(3,4,gen)=mlcutjj
          scut(4,3,gen)=mlcutjj
          if(shtr.ge.2)then
            ptcut(5,gen)=ptcuts
            ptcut(6,gen)=ptcuts
            ycut(5,gen)=ycuts
            ycut(6,gen)=ycuts
          endif
          sucuth=mucuth*mucuth
          slcuth=mlcuth*mlcuth

        elseif(scuts(gen).ne.0)then
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
      subroutine settings_Hjj(kbeam,k,lweight,n,gen,step)
***********************************************************************
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
*     31.10.06  Ansgar Denner     last changed  15.02.12              *
*---------------------------------------------------------------------*
*     edited by Stefan Kallweit   last changed  11.10.11              *
*---------------------------------------------------------------------*
***********************************************************************
      implicit none

#include "mcparams.h"
#include "mccommon.h"

c local variables
      real*8 kbeam(2,0:3),k(maxe,0:3),lweight(maxcon),lenergy,x
      real*8 dot1,dot2
      real*8 kt(maxe),y(maxe),ktH
      real*8 denum,kh(0:3),khvec
      complex*16 eiphj1,eiphj2
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

c      print *, 'settings_Hjj for pp -> qqH'

      if(step.eq.2)then
        do icon=1,maxcon
          if (lweight(icon).ne.0d0) goto 100
        enddo
        return
      endif

 100  continue
      switch=0
      lenergy=kbeam(1,0)+kbeam(2,0)

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



c declaration: histogram(hname,y,x,min,max,parts,histo,n,addpoint,step)

      histo=0

      if(lnoutmc.ge.3)then
ccccccccccccccccccccccccc
c  total cross section  c
ccccccccccccccccccccccccc
      histo=0
      hname='dat.tot'
      call histogram(hname,lweight,lenergy+0.5d0,lenergy,lenergy+1d0,1,
     *  histo,n,0,gen,step)
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
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)


c jets are pt-ordered by cut routine

      kt(3)=sqrt(k(3,1)**2+k(3,2)**2)
      kt(4)=sqrt(k(4,1)**2+k(4,2)**2)
      kt(lastp)=sqrt(k(lastp,1)**2+k(lastp,2)**2)


      jet1=3
      jet2=4


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distributions in pt of leading and subleading jet  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.ptjmax1pt'
      x=sqrt(k(jet1,1)**2 + k(jet1,2)**2)
c     if(x.eq.0d0) x=-1d0
      call histogram(hname,lweight,x,0d0,400d0,nbin,histo,n,1,gen,step)
      histo=histo+1
      hname='dat.ptjmax2pt'
c     if(x.eq.0d0) x=-1d0
      x=sqrt(k(jet2,1)**2 + k(jet2,2)**2)
      call histogram(hname,lweight,x,0d0,150d0,nbin,histo,n,1,gen,step)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distributions in y of leading and subleading jet  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.yjmax1pt'
      if(k(jet1,0)-k(jet1,3).ne.0d0.and.k(jet1,0)+k(jet1,3).ne.0d0) then
         y(jet1)=.5d0*log((k(jet1,0)+k(jet1,3))/(k(jet1,0)-k(jet1,3)))
      else
         y(jet1)=1d10
      endif
      x=y(jet1)
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)
      histo=histo+1
      hname='dat.yjmax2pt'
      if(k(jet2,0)-k(jet2,3).ne.0d0.and.k(jet2,0)+k(jet2,3).ne.0d0) then
         y(jet2)=.5d0*log((k(jet2,0)+k(jet2,3))/(k(jet2,0)-k(jet2,3)))
      else
         y(jet2)=-1d10
      endif
      x=y(jet2)
      call histogram(hname,lweight,x,-5d0,5d0,nbin,histo,n,1,gen,step)




cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in dy of leading and subleading jet  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.dyjjmaxpt'
      if(y(jet1).lt.1d10.and.y(jet2).gt.-1d10)then
         x=y(jet1)-y(jet2)
      else
         x=2d10
      endif


      call histogram(hname,lweight,x,-10d0,10d0,nbin,histo,n,1,gen,step)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in dphi of leading and subleading jet  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.dphijjmaxpt'
      if(kt(jet1)*kt(jet2).ne.0d0)then
c         x=(k(jet1,1)*k(jet2,1)+k(jet1,2)*k(jet2,2))/(kt(jet1)*kt(jet2))

	 eiphj1 = dcmplx(k(jet1,1),k(jet1,2))/kt(jet1)
	 eiphj2 = dcmplx(k(jet2,1),k(jet2,2))/kt(jet2)
	 x = dimag(log(eiphj1/eiphj2))*180d0/pi
         if (y(jet1).lt.y(jet2)) x=-x

      else
         x=1000d0
      endif

      call histogram(hname,lweight,x,-180d0,180d0,
     &     nbin,histo,n,1,gen,step)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in invariant mass of leading and subleading jet  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.mjjmaxpt'
      x=((k(jet1,0)+k(jet2,0))**2  -(k(jet1,1)+k(jet2,1))**2
     &     -(k(jet1,2)+k(jet2,2))**2  -(k(jet1,3)+k(jet2,3))**2)
      if(x.le.0d0) then
         x=0d0
      else
         x=sqrt(x)
      endif


      call histogram(hname,lweight,x,0d0,3000d0,nbin,
     &     histo,n,1,gen,step)





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in invariant mass of leading and subleading jet  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.q2v1'
      dot1=k(1,1)*k(jet1,1)+k(1,2)*k(jet1,2)
     & +k(1,3)*k(jet1,3)  
      dot2=k(1,1)*k(jet2,1)+k(1,2)*k(jet2,2)
     & +k(1,3)*k(jet2,3)  


      if(dot1.ge.0) then
         x=((k(1,0) - k(jet1,0))**2  -(k(1,1)-k(jet1,1))**2
     & -(k(1,2)+k(jet1,2))**2  -(k(1,3)-k(jet1,3))**2)
      else
         x=((k(1,0) - k(jet2,0))**2  -(k(1,1)-k(jet2,1))**2
     & -(k(1,2)+k(jet2,2))**2  -(k(1,3)-k(jet2,3))**2)
      endif
      
      if(x.le.0d0) then
            
         x=-x
      else
         x=x
      endif

      x=sqrt(x)
      call histogram(hname,lweight,x,0d0,1000d0,nbin,
     &     histo,n,1,gen,step)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  distribution in invariant mass of leading and subleading jet  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      histo=histo+1
      hname='dat.q2v2'
      dot1=k(2,1)*k(jet1,1)+k(2,2)*k(jet1,2)
     & +k(2,3)*k(jet1,3)  
      dot2=k(2,1)*k(jet2,1)+k(2,2)*k(jet2,2)
     & +k(2,3)*k(jet2,3)  


      if(dot1.ge.0) then
         x=((k(2,0) - k(jet1,0))**2  -(k(2,1)-k(jet1,1))**2
     & -(k(2,2)+k(jet1,2))**2  -(k(2,3)-k(jet1,3))**2)
      else
         x=((k(2,0) - k(jet2,0))**2  -(k(2,1)-k(jet2,1))**2
     & -(k(2,2)+k(jet2,2))**2  -(k(2,3)-k(jet2,3))**2)
      endif

      
      if(x.le.0d0) then
         x=-x
      else
         x=x
      endif

      x=sqrt(x)

      
      
      
c      print *,"----------------------------------------------"
c      print *, "pt comp",kt(jet1),kt(jet2),sqrt(kh(1)**2+kh(2)**2)
c      print *, "Sum pt x: ", k(jet1,1) + k(jet2,1) + kh(1)
c      print *, "Sum pt y: ", k(jet1,2) + k(jet2,2) + kh(2)
      call histogram(hname,lweight,x,0d0,1000d0,nbin,
     &     histo,n,1,gen,step)

c     Write event out
c      print *, lweight
c      print *, "dpne"

      IF ( ( lweight(1) .NE. 0) .and. (step .EQ. 2 ) )   THEN
      print *,"writing event with weight :",lweight(1)   
      write(1988,"(A)") "<event>"
      write(1988,Fmt2) 5, 60, lweight(1) ,1.0 ,1.0 ,1.0
   
      write(1988,Fmt1)-1,-1,0,0,0,501,k(1,1),k(1,2),k(1,3),k(1,0),0.0,
     &0.0,0.0
      write(1988,Fmt1)3,-1,0,0,502,0,k(2,1),k(2,2),k(2,3),k(2,0),0.0,
     &0.0,0.0
      write(1988,Fmt1)25,1,1,2,0,0,kh(1),kh(2),kh(3),kh(0),125.0,
     &0.0,0.0
      write(1988,Fmt1)-2,1,1,2,0,501,k(jet1,1),k(jet1,2),k(jet1,3),
     &k(jet1,0),0.0,0.0,0.0
      write(1988,Fmt1)4,1,1,2,502,0,k(jet2,1),k(jet2,2),k(jet2,3),
     &k(jet2,0),0.0,0.0,0.0
     
      write(1988,"(A)") "</event>"

      ELSE
      PRINT  *,"0Weight"   

      ENDIF
      

      return
   
      





      end

      


