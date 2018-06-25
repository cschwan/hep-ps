!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialization generator                                    c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lusifer_initphasespace(name,generator,lightfermions, &
        includecuts,sout)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      real*8 power(maxv),m2,m3,e2,e3,scutinv
      integer idhep(maxe,maxe),binary(maxe,maxe),idhep2,idhep3
      integer i1,i2,i3,ns,nt,maxns,maxnt,binary1,binary2,binary3
      integer in1(maxe),in2(maxe),out1(maxe),out2(maxe),lusifer_id
      integer lightfermions,virt(maxe),channel,generator,includecuts
      integer prop2,prop3,sout
      character*3 gname(-maxv:maxv),name(maxe)
      logical lusifer_vertex,lusifer_included,exist,same, &
        lusifer_compareinv
      logical lusifer_comparedecay,lusifer_compareprocess
! general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
! cinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
! cprocess
      real*8 powerprocess(maxe,maxch,maxg),ccutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
! cdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
! cdensity
      integer nsinv(maxch,maxg),chinv(maxch,maxg),maxinv(maxg)
      integer ntprocess(maxch,maxg),chprocess(maxch,maxg)
      integer maxprocess(maxg),nsdecay(maxch,maxg),chdecay(maxch,maxg)
      integer maxdecay(maxg),numinv(maxe,maxch,maxg)
      integer numprocess(maxe,maxch,maxg),numdecay(maxe,maxch,maxg)
! output
      integer nout,numout,maxout
! techparam
      real*8 a,techcut
! cuts
      real*8 ecut(maxe),scut(maxe,maxe),ccut(maxe,maxe)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width, &
        nchannel,nexternal,allbinary
      common/lusifer_cinv/powerinv,mcutinv,ininv,idhepinv,ninv,lmin,lmax
      common/lusifer_cdecay/indecay,out1decay,out2decay,ndecay
      common/lusifer_cprocess/powerprocess,ccutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      common/lusifer_cdensity/nsinv,chinv,maxinv,ntprocess,chprocess, &
        maxprocess,nsdecay,chdecay,maxdecay,numinv,numprocess,numdecay
      common/lusifer_output/nout,numout,maxout
      common/lusifer_techparam/a,techcut
      common/lusifer_cuts/ecut,scut,ccut
! technical parameter in h function and subroutine process
      a=1d-6
      techcut=1d-7
! generic particle names, see also subroutine vertexg
      do i1=-maxv,maxv
        gname(i1)='   '
      enddo
      gname(1)='dq '
      gname(2)='uq '
      gname(3)='sq '
      gname(4)='cq '
      gname(5)='bq '
      gname(6)='tq '
      gname(11)='el '
      gname(12)='ne '
      gname(13)='mu '
      gname(14)='nm '
      gname(15)='ta '
      gname(16)='nt '
      gname(22)='ga '
      gname(23)='Z0 '
      gname(24)='W+ '
      gname(25)='H0 '
      gname(26)='gl '
      gname(30)='v1 ' ! auxiliary particle, 4-particle vertex
      gname(31)='v2 ' ! auxiliary particle, 4-particle vertex
      gname(32)='v3 ' ! auxiliary particle, 4-particle vertex
      do i1=1,maxv
        gname(-i1)=gname(i1)
        if(gname(i1).ne.'   '.and.gname(i1).ne.'ga '.and. &
          gname(i1).ne.'Z0 '.and.gname(i1).ne.'H0 '.and. &
          gname(i1).ne.'gl ')then
          gname(-i1)(3:3)='~'
        endif
      enddo
! particle identity
      do i1=1,nexternal(generator)
        idhep(i1,1)=0
        do i2=-maxv,maxv
          if(i1.le.2.and.name(i1).eq.gname(i2))idhep(i1,1)=i2
          if(i1.ge.3.and.name(i1).eq.gname(-i2))idhep(i1,1)=i2
        enddo
        if(idhep(i1,1).eq.0)then
          write(nout,'(a)')' Wrong particle name!'
          stop
        endif
      enddo
! initializing invariants
      do i1=1,nexternal(generator)
        s(2**(i1-1))=mass(abs(idhep(i1,1)))**2
      enddo
! mapping
      do i1=1,maxv
        power(i1)=0.9d0
        if(width(i1).gt.0d0)power(i1)=0d0
        if(gname(i1)(1:1).eq.'v')then
          mass(i1)=0d0
          width(i1)=0d0
          power(i1)=0d0
        endif
      enddo
! check with Madgraph
!      do i1=1,25
!        width(i1)=0d0
!        power(i1)=a+1d-5*i1
!      enddo
! initialize binary counting of particle number
      do i1=1,maxe
      do i2=1,maxe
        binary(i1,i2)=0
      enddo
      enddo
      allbinary(generator)=0
      do i1=1,nexternal(generator)
        binary(i1,1)=2**(i1-1)
        allbinary(generator)=allbinary(generator)+2**(i1-1)
      enddo
! initializing external masses
      do i1=1,nexternal(generator)
        s(allbinary(generator)-2**(i1-1))=s(2**(i1-1))
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize cuts                                             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i1=0,allbinary(generator)
        mcutinv(i1,generator)=0d0
        if(i1.ge.1)ccutprocess(i1,generator)=1d0
      enddo
! angular cut
      if(includecuts.eq.1)then
        do i1=1,2
        do i2=3,nexternal(generator)
          if(ccut(i1,i2).gt.0d0.and.ccut(i1,i2).lt.1d0)then
            ccutprocess(2**(i1-1)+2**(i2-1),generator)=ccut(i1,i2)
          endif
        enddo
        enddo
      endif
! invariant-mass cut
      do i1=1,allbinary(generator)
        scutinv=0d0
        do i2=1,nexternal(generator)
          if(lusifer_included(i1,2**(i2-1),nexternal(generator)))then
             scutinv=scutinv+mass(abs(idhep(i2,1)))
           endif
        enddo
        scutinv=scutinv**2
        if(includecuts.eq.1)then
          do i2=1,nexternal(generator)
          do i3=i2+1,nexternal(generator)
            if(lusifer_included(i1,2**(i2-1),nexternal(generator)).and. &
              lusifer_included(i1,2**(i3-1),nexternal(generator)))then
              m2=mass(abs(idhep(i2,1)))
              m3=mass(abs(idhep(i3,1)))
              e2=max(ecut(i2),m2)
              e3=max(ecut(i3),m3)
              scutinv=scutinv-(m2+m3)**2 &
                +max(scut(i2,i3),m2**2+m3**2+2d0*e2*e3 &
                  -2d0*dsqrt((e2**2-m2**2)*(e3**2-m3**2))*ccut(i2,i3))
            endif
          enddo
          enddo
        endif
        mcutinv(i1,generator)=dsqrt(scutinv)
        do i2=2,i1-1
          if(lusifer_included(i1,i2,nexternal(generator)))then
            mcutinv(i1,generator)=max(mcutinv(i1,generator), &
              mcutinv(i2,generator)+mcutinv(i1-i2,generator))
          endif
        enddo
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize channels                                         c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      channel=nchannel(generator)+1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     s-channel propagators                                       c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! number of s-channel propagators
      do maxns=0,nexternal(generator)-3
        ns=0
        if(maxns.eq.0)goto 500
! ns'th s-channel propagator
 100    ns=ns+1
        out1(ns)=2
! decay product 1 for ns'th decay
 200    out1(ns)=out1(ns)+1
        out2(ns)=out1(ns)
! decay product 2 for ns'th decay
 300    out2(ns)=out2(ns)+1
        virt(ns)=-maxv-1
! virtual particle for ns'th decay
 400    virt(ns)=virt(ns)+1
        if(gname(virt(ns)).eq.'   '.and.virt(ns).lt.maxv)goto 400
        if(gname(virt(ns)).eq.'   '.and.virt(ns).eq.maxv)goto 1500
        if(virt(ns).lt.0.and. &
          gname(virt(ns)).eq.gname(-virt(ns)))goto 400
! checking whether 3-particle vertex exists
        if(.not.lusifer_vertex(gname(idhep(out1(ns),ns)), &
          gname(idhep(out2(ns),ns)),gname(virt(ns)), &
          lightfermions))goto 1500
! checking last 3-particle vertex
        if(ns.eq.nexternal(generator)-3.and. &
          .not.lusifer_vertex(gname(idhep(in1(ns),ns)), &
            gname(idhep(in2(ns),ns)),gname(-virt(ns)), &
            lightfermions))goto 1500
! initializing next step
        do i2=1,nexternal(generator)
          binary(i2,ns+1)=binary(i2,ns)
          idhep(i2,ns+1)=idhep(i2,ns)
        enddo
! combining particle out1 and out2 into new external particle out2
        binary(out1(ns),ns+1)=0
        idhep(out1(ns),ns+1)=0
        binary(out2(ns),ns+1)= &
          binary(out1(ns),ns)+binary(out2(ns),ns)
        idhep(out2(ns),ns+1)=-virt(ns)
! find next s-channel propagator for ns < maxns
        if(ns.lt.maxns)goto 100
 500    continue
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     t-channel propagators                                       c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! number of t-channel propagators
        maxnt=nexternal(generator)-3
        if(maxns.eq.maxnt)goto 1000
        nt=maxns
! nt'th t-channel propagator
 600    nt=nt+1
        in1(nt)=0
! incoming particle 1 for nt's process
 700    in1(nt)=in1(nt)+1
! incoming particle 2 for nt's process
        in2(nt)=3-in1(nt)
        out1(nt)=2
! outgoing particle 1 for nt's process
 800    out1(nt)=out1(nt)+1
! virtual particle for nt's process
        virt(nt)=-maxv-1
 900    virt(nt)=virt(nt)+1
        if(gname(virt(nt)).eq.'   '.and.virt(nt).lt.maxv)goto 900
        if(gname(virt(nt)).eq.'   '.and.virt(nt).eq.maxv)goto 1300
        if(virt(nt).lt.0.and. &
          gname(virt(nt)).eq.gname(-virt(nt)))goto 900
! avoid doube counting of diagrams
        if(nt.ge.2.and.in1(nt-1).eq.in1(nt))goto 1300
! checking whether 3-particle vertex exists
        if(.not.lusifer_vertex(gname(idhep(in1(nt),nt)), &
          gname(idhep(out1(nt),nt)),gname(virt(nt)), &
          lightfermions))goto 1300
! initializing next step
        do i2=1,nexternal(generator)
          binary(i2,nt+1)=binary(i2,nt)
          idhep(i2,nt+1)=idhep(i2,nt)
        enddo
! combining particle in1 and out1 into new external particle in1
        binary(out1(nt),nt+1)=0
        idhep(out1(nt),nt+1)=0
        binary(in1(nt),nt+1)=binary(in1(nt),nt)+binary(out1(nt),nt)
        idhep(in1(nt),nt+1)=-virt(nt)
! find next t-channel propagator for nt < maxnt
        if(nt.lt.maxnt)goto 600
 1000   continue
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     checking last vertex                                        c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        i1=nexternal(generator)-2
        do i2=3,nexternal(generator)
          if(idhep(i2,i1).ne.0)i3=i2
        enddo
        if(.not.lusifer_vertex(gname(idhep(1,i1)),gname(idhep(2,i1)), &
          gname(idhep(i3,i1)),lightfermions))goto 1300
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initializing channels                                       c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ninv(channel,generator)=nexternal(generator)-4
        nprocess(channel,generator)=maxnt-maxns
        ndecay(channel,generator)=maxns
! initializing inv
        do i1=1,maxns
          ininv(i1,channel,generator)= &
            binary(out1(i1),i1)+binary(out2(i1),i1)
          idhepinv(i1,channel,generator)=abs(virt(i1))
          powerinv(i1,channel,generator)=power(abs(virt(i1)))
        enddo
        do i1=maxns+1,maxnt-1
          i2=maxns+maxnt-i1
          ininv(i2,channel,generator)=allbinary(generator)  &
            -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          idhepinv(i2,channel,generator)=0
          powerinv(i2,channel,generator)=0d0
        enddo
! minimal invariant-mass cuts for subroutine inv
        do i1=1,ninv(channel,generator)
        do i2=1,i1-1
          binary1=ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmin(i2,i1,channel,generator)=.false.
          if(lusifer_included(binary1,binary2,nexternal(generator)))then
            lmin(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(lusifer_included(binary3,binary2,nexternal(generator)) &
                .and.i2.ne.i3)then
                lmin(i2,i1,channel,generator)=.false.
              endif
            enddo
          endif
        enddo
        enddo
! maximal invariant-mass cuts for subroutine inv
        do i1=1,ninv(channel,generator)
        do i2=1,i1-1
          binary1=allbinary(generator)-3-ininv(i1,channel,generator)
          binary2=ininv(i2,channel,generator)
          lmax(i2,i1,channel,generator)=.false.
          if(lusifer_included(binary1,binary2,nexternal(generator)))then
            lmax(i2,i1,channel,generator)=.true.
            do i3=1,i1-1
              binary3=ininv(i3,channel,generator)
              if(lusifer_included(binary3,binary2,nexternal(generator)) &
                .and.i2.ne.i3)then
                lmax(i2,i1,channel,generator)=.false.
              endif
            enddo
          endif
        enddo
        enddo
! initializing process
        do i1=maxns+1,maxnt
          i2=i1-maxns
          in1process(i2,channel,generator)=allbinary(generator) &
            -binary(in1(i1),i1)
          in2process(i2,channel,generator)=allbinary(generator) &
            -binary(in2(i1),i1)
          out1process(i2,channel,generator)=binary(out1(i1),i1)
          out2process(i2,channel,generator)=allbinary(generator)  &
            -binary(in1(i1),i1)-binary(in2(i1),i1)-binary(out1(i1),i1)
          inprocess(i2,channel,generator)=allbinary(generator) &
            -binary(in1(i1),i1)-binary(in2(i1),i1)
          virtprocess(i2,channel,generator)=allbinary(generator) &
            -binary(in1(i1),i1)-binary(out1(i1),i1)
          idhepprocess(i2,channel,generator)=abs(virt(i1))
          powerprocess(i2,channel,generator)=power(abs(virt(i1)))
        enddo
! initializing decay
        do i1=1,maxns
          i2=maxns-i1+1
          indecay(i2,channel,generator)= &
            binary(out1(i1),i1)+binary(out2(i1),i1)
          out1decay(i2,channel,generator)=binary(out1(i1),i1)
          out2decay(i2,channel,generator)=binary(out2(i1),i1)
        enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether generator already exists                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i1=1,channel-1
          same=.true.
! checking number of s-channel propagators
          prop2=ndecay(channel,generator)
          do i2=1,ndecay(channel,generator)
            idhep2=idhepinv(i2,channel,generator)
            if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and. &
              power(idhep2).eq.0d0)prop2=prop2-1
          enddo
          prop3=ndecay(i1,generator)
          do i3=1,ndecay(i1,generator)
            idhep3=idhepinv(i3,i1,generator)
            if(mass(idhep3).eq.0d0.and.width(idhep3).eq.0d0.and. &
              power(idhep3).eq.0d0)prop3=prop3-1
          enddo
! checking s-channel propagators
          if(prop2.eq.prop3)then
            do i2=1,ndecay(i1,generator)
              binary2=ininv(i2,i1,generator)
              idhep2=idhepinv(i2,i1,generator)
              if(mass(idhep2).ne.0d0.or.width(idhep2).ne.0d0.or. &
                power(idhep2).ne.0d0)then
                exist=.false.
                do i3=1,ndecay(channel,generator)
                  binary3=ininv(i3,channel,generator)
                  idhep3=idhepinv(i3,channel,generator)
                  if(binary2.eq.binary3.and. &
                    mass(idhep2).eq.mass(idhep3).and. &
                    width(idhep2).eq.width(idhep3).and. &
                    power(idhep2).eq.power(idhep3))then
                      exist=.true.
                  endif
                enddo
                if(.not.exist)same=.false.
              endif
            enddo
          else
            same=.false.
          endif
! checking number of t-channel propagators
          prop2=nprocess(channel,generator)
          do i2=1,nprocess(channel,generator)
            idhep2=idhepprocess(i2,channel,generator)
            if(mass(idhep2).eq.0d0.and.width(idhep2).eq.0d0.and. &
              power(idhep2).eq.0d0)prop2=prop2-1
          enddo
          prop3=nprocess(i1,generator)
          do i3=1,nprocess(i1,generator)
            idhep3=idhepprocess(i3,i1,generator)
            if(mass(idhep3).eq.0d0.and.width(idhep3).eq.0d0.and. &
              power(idhep3).eq.0d0)prop3=prop3-1
          enddo
! checking t-channel propagators
          if(prop2.eq.prop3)then
            do i2=1,nprocess(i1,generator)
              binary2=min(virtprocess(i2,i1,generator), &
                allbinary(generator)-virtprocess(i2,i1,generator))
              idhep2=idhepprocess(i2,i1,generator)
              if(mass(idhep2).ne.0d0.or.width(idhep2).ne.0d0.or. &
                power(idhep2).ne.0d0)then
                exist=.false.
                do i3=1,nprocess(channel,generator)
                  binary3=min(virtprocess(i3,channel,generator), &
                    allbinary(generator) &
                    -virtprocess(i3,channel,generator))
                  idhep3=idhepprocess(i3,channel,generator)
                  if(binary2.eq.binary3.and. &
                    mass(idhep2).eq.mass(idhep3).and. &
                    width(idhep2).eq.width(idhep3).and. &
                    power(idhep2).eq.power(idhep3))then
                      exist=.true.
                  endif
                enddo
                if(.not.exist)same=.false.
              endif
            enddo
          else
            same=.false.
          endif
          if(same)goto 1300
        enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     output                                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(sout.eq.2)then
          write(nout,*)'channel=',channel
          do i2=1,maxns
            write(nout,*)'s channel: ', &
              lusifer_id(binary(out1(i2),i2)),lusifer_id( &
              binary(out2(i2),i2)),' ->', &
              lusifer_id(binary(out2(i2),i2+1)),'   ', &
              gname(-idhep(out1(i2),i2)), &
              gname(-idhep(out2(i2),i2)),' -> ', &
              gname(-idhep(out2(i2),i2+1))
          enddo
          do i2=maxns+1,maxnt
            write(nout,*)'t channel: ', &
              lusifer_id(binary(in1(i2),i2)),lusifer_id( &
              binary(out1(i2),i2)),' ->', &
              lusifer_id(binary(in1(i2),i2+1)),'   ', &
              gname(idhep(in1(i2),i2)), &
              gname(-idhep(out1(i2),i2)),' -> ', &
              gname(idhep(in1(i2),i2+1))
          enddo
          write(nout,'(a)')'   '
        endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     end of loops                                                c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        channel=channel+1
 1300   if(maxns.lt.maxnt)then
 1400     if(virt(nt).lt.maxv)goto 900
          if(out1(nt).lt.nexternal(generator))goto 800
          if(in1(nt).lt.2)goto 700
          if(nt.gt.maxns+1)then
            nt=nt-1
            goto 1400
          endif
        endif
! end of t-channel propagators
 1500   if(maxns.gt.0)then
 1600     if(virt(ns).lt.maxv)goto 400
          if(out2(ns).lt.nexternal(generator))goto 300
          if(out1(ns).lt.nexternal(generator)-1)goto 200
          if(ns.gt.1)then
            ns=ns-1
            goto 1600
          endif
        endif
! end of s-channel propagators
      enddo
! number of channels
      nchannel(generator)=channel-1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     initialize densities                                        c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initializing inv
      maxinv(generator)=0
      do channel=1,nchannel(generator)
      do ns=1,ninv(channel,generator)
        numinv(ns,channel,generator)=0
        do i1=1,channel-1
        do i2=1,ninv(i1,generator)
          if(lusifer_compareinv(ns,channel,i2,i1,generator)) &
            numinv(ns,channel,generator)=numinv(i2,i1,generator)
        enddo
        enddo
        if(numinv(ns,channel,generator).eq.0)then
          maxinv(generator)=maxinv(generator)+1
          numinv(ns,channel,generator)=maxinv(generator)
          nsinv(maxinv(generator),generator)=ns
          chinv(maxinv(generator),generator)=channel
        endif
      enddo
      enddo
! process
      maxprocess(generator)=0
      do channel=1,nchannel(generator)
      do nt=1,nprocess(channel,generator)
        numprocess(nt,channel,generator)=0
        do i1=1,channel-1
        do i2=1,nprocess(i1,generator)
          if(lusifer_compareprocess(nt,channel,i2,i1,generator)) &
            numprocess(nt,channel,generator)=numprocess(i2,i1,generator)
        enddo
        enddo
        if(numprocess(nt,channel,generator).eq.0)then
          maxprocess(generator)=maxprocess(generator)+1
          numprocess(nt,channel,generator)=maxprocess(generator)
          ntprocess(maxprocess(generator),generator)=nt
          chprocess(maxprocess(generator),generator)=channel
        endif
      enddo
      enddo
! decay
      maxdecay(generator)=0
      do channel=1,nchannel(generator)
      do ns=1,ndecay(channel,generator)
        numdecay(ns,channel,generator)=0
        do i1=1,channel-1
        do i2=1,ndecay(i1,generator)
          if(lusifer_comparedecay(ns,channel,i2,i1,generator)) &
            numdecay(ns,channel,generator)=numdecay(i2,i1,generator)
        enddo
        enddo
        if(numdecay(ns,channel,generator).eq.0)then
          maxdecay(generator)=maxdecay(generator)+1
          numdecay(ns,channel,generator)=maxdecay(generator)
          nsdecay(maxdecay(generator),generator)=ns
          chdecay(maxdecay(generator),generator)=channel
        endif
      enddo
      enddo
      if(sout.ne.0)then
        write(nout,'(a)')' '
        write(nout,'(a)')' Phase-space generator:'
        write(nout,'(" number of channels            =",i6)') &
          nchannel(generator)
        write(nout,'(" calculation of invariants     =",i6)') &
          maxinv(generator)
        write(nout,'(" calculation of 2->2 processes =",i6)') &
          maxprocess(generator)
        write(nout,'(" calculation of 1->2 decays    =",i6)') &
          maxdecay(generator)
      endif
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether two 1->2 decays are equal                     c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_comparedecay(ns1,ch1,ns2,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg
      parameter(maxe=9,maxch=20000,maxg=1)
      integer ns1,ns2,ch1,ch2,generator
      logical lusifer_comparedecay
! cdecay
      integer indecay(maxe,maxch,maxg),out1decay(maxe,maxch,maxg)
      integer out2decay(maxe,maxch,maxg),ndecay(maxch,maxg)
      common/lusifer_cdecay/indecay,out1decay,out2decay,ndecay
      lusifer_comparedecay=.false.
      if(indecay(ns1,ch1,generator).ne.indecay(ns2,ch2,generator)) &
        return
      if(out1decay(ns1,ch1,generator).ne.out1decay(ns2,ch2,generator) &
        .and.out1decay(ns1,ch1,generator).ne. &
        out2decay(ns2,ch2,generator))return
      lusifer_comparedecay=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether two 2->2 processes are equal                  c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_compareprocess(ns1,ch1,ns2,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      integer ns1,ns2,ch1,ch2,generator,idhep1,idhep2
      logical lusifer_compareprocess
! general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
! cprocess
      real*8 powerprocess(maxe,maxch,maxg),ccutprocess(2**maxe,maxg)
      integer in1process(maxe,maxch,maxg),in2process(maxe,maxch,maxg)
      integer out1process(maxe,maxch,maxg),out2process(maxe,maxch,maxg)
      integer inprocess(maxe,maxch,maxg),virtprocess(maxe,maxch,maxg)
      integer idhepprocess(maxe,maxch,maxg),nprocess(maxch,maxg)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width, &
        nchannel,nexternal,allbinary
      common/lusifer_cprocess/powerprocess,ccutprocess,in1process, &
        in2process,out1process,out2process,inprocess,virtprocess, &
        idhepprocess,nprocess
      lusifer_compareprocess=.false.
      if(inprocess(ns1,ch1,generator).ne.inprocess(ns2,ch2,generator)) &
        return
      if(virtprocess(ns1,ch1,generator).ne. &
        virtprocess(ns2,ch2,generator).and. &
        virtprocess(ns1,ch1,generator).ne. &
        allbinary(generator)-virtprocess(ns2,ch2,generator))return
      if(powerprocess(ns1,ch1,generator).ne. &
        powerprocess(ns2,ch2,generator))return
      idhep1=idhepprocess(ns1,ch1,generator)
      idhep2=idhepprocess(ns2,ch2,generator)
      if(mass(idhep1).ne.mass(idhep2).or.width(idhep1).ne.width(idhep2)) &
        return
      lusifer_compareprocess=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check whether caculation of invariants are equal            c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_compareinv(ns1,ch1,ns2,ch2,generator)
      implicit none
! local variables
      integer maxe,maxch,maxg,maxv
      parameter(maxe=9,maxch=20000,maxg=1,maxv=40)
      integer i1,i2,ns1,ns2,ch1,ch2,generator,idhep1,idhep2
      logical lusifer_compareinv,lusifer_included
! general
      real*8 alphaisr,scale,meisr,s(2**maxe),p(0:3,2**maxe)
      real*8 mass(0:maxv),width(0:maxv)
      integer nchannel(maxg),nexternal(maxg),allbinary(maxg)
! cinv
      real*8 powerinv(maxe,maxch,maxg),mcutinv(0:2**maxe,maxg)
      integer ininv(maxe,maxch,maxg),idhepinv(maxv,maxch,maxg)
      integer ninv(maxch,maxg)
      logical lmin(maxe,maxe,maxch,maxg),lmax(maxe,maxe,maxch,maxg)
      common/lusifer_general/alphaisr,scale,meisr,s,p,mass,width, &
        nchannel,nexternal,allbinary
      common/lusifer_cinv/powerinv,mcutinv,ininv,idhepinv,ninv,lmin,lmax
      lusifer_compareinv=.false.
      if(ininv(ns1,ch1,generator).ne.ininv(ns2,ch2,generator))return
      if(powerinv(ns1,ch1,generator).ne.powerinv(ns2,ch2,generator)) &
        return
      idhep1=idhepinv(ns1,ch1,generator)
      idhep2=idhepinv(ns2,ch2,generator)
      if(mass(idhep1).ne.mass(idhep2).or.width(idhep1).ne.width(idhep2)) &
        return
      do i1=1,ns1-1
        if(lmin(i1,ns1,ch1,generator))then
          lusifer_included=.false.
          do i2=1,ns2-1
            if(lmin(i2,ns2,ch2,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              lusifer_included=.true.
          enddo
          if(.not.lusifer_included)return
        endif
        if(lmax(i1,ns1,ch1,generator))then
          lusifer_included=.false.
          do i2=1,ns2-1
            if(lmax(i2,ns2,ch2,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              lusifer_included=.true.
          enddo
          if(.not.lusifer_included)return
        endif
      enddo
      do i2=1,ns2-1
        if(lmin(i2,ns2,ch2,generator))then
          lusifer_included=.false.
          do i1=1,ns1-1
            if(lmin(i1,ns1,ch1,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              lusifer_included=.true.
          enddo
          if(.not.lusifer_included)return
        endif
        if(lmax(i2,ns2,ch2,generator))then
          lusifer_included=.false.
          do i1=1,ns1-1
            if(lmax(i1,ns1,ch1,generator).and. &
              ininv(i1,ch1,generator).eq.ininv(i2,ch2,generator)) &
              lusifer_included=.true.
          enddo
          if(.not.lusifer_included)return
        endif
      enddo
      lusifer_compareinv=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     check if binary2 is included in binary1                     c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_included(binary1,binary2,nexternal)
      implicit none
! local variables
      integer binary1,binary2,i1,b1,b2,nexternal
      logical lusifer_included
      lusifer_included=.false.
      do i1=1,nexternal
        b1=binary1/2**(i1-1)
        b2=binary2/2**(i1-1)
        if(2*(b2/2).ne.b2.and.2*(b1/2).eq.b1)return
      enddo
      lusifer_included=.true.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     3-particle vertex                                           c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_vertex(p1,p2,p3,lightfermions)
      implicit none
! local variables
      character*3 p1,p2,p3
      integer lightfermions
      logical lusifer_vertex,lusifer_vertexg
      lusifer_vertex=.true.
      if(lusifer_vertexg(p1,p2,p3,lightfermions))return
      if(lusifer_vertexg(p1,p3,p2,lightfermions))return
      if(lusifer_vertexg(p2,p1,p3,lightfermions))return
      if(lusifer_vertexg(p2,p3,p1,lightfermions))return
      if(lusifer_vertexg(p3,p1,p2,lightfermions))return
      if(lusifer_vertexg(p3,p2,p1,lightfermions))return
      lusifer_vertex=.false.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     generic 3-particle vertex                                   c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_vertexg(p1,p2,p3,lightfermions)
      implicit none
! local variables
      character*3 p1,p2,p3
      integer lightfermions
      logical lusifer_vertexg
      lusifer_vertexg=.true.
! ew three-particle vertices
      if(p1.eq.'el '.and.p2.eq.'ne~'.and.p3.eq.'W+ ')return
      if(p1.eq.'el~'.and.p2.eq.'ne '.and.p3.eq.'W+~')return
      if(p1.eq.'mu '.and.p2.eq.'nm~'.and.p3.eq.'W+ ')return
      if(p1.eq.'mu~'.and.p2.eq.'nm '.and.p3.eq.'W+~')return
      if(p1.eq.'ta '.and.p2.eq.'nt~'.and.p3.eq.'W+ ')return
      if(p1.eq.'ta~'.and.p2.eq.'nt '.and.p3.eq.'W+~')return
      if(p1.eq.'dq '.and.p2.eq.'uq~'.and.p3.eq.'W+ ')return
      if(p1.eq.'dq~'.and.p2.eq.'uq '.and.p3.eq.'W+~')return
      if(p1.eq.'sq '.and.p2.eq.'cq~'.and.p3.eq.'W+ ')return
      if(p1.eq.'sq~'.and.p2.eq.'cq '.and.p3.eq.'W+~')return
      if(p1.eq.'bq '.and.p2.eq.'tq~'.and.p3.eq.'W+ ')return
      if(p1.eq.'bq~'.and.p2.eq.'tq '.and.p3.eq.'W+~')return
      if(p1.eq.'ne '.and.p2.eq.'ne~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'el '.and.p2.eq.'el~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'el '.and.p2.eq.'el~'.and.p3.eq.'ga ')return
      if(p1.eq.'nm '.and.p2.eq.'nm~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'mu '.and.p2.eq.'mu~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'mu '.and.p2.eq.'mu~'.and.p3.eq.'ga ')return
      if(p1.eq.'nt '.and.p2.eq.'nt~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'ta '.and.p2.eq.'ta~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'ta '.and.p2.eq.'ta~'.and.p3.eq.'ga ')return
      if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'ga ')return
      if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'ga ')return
      if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'ga ')return
      if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'ga ')return
      if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'ga ')return
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'Z0 ')return
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'ga ')return
      if(p1.eq.'ga '.and.p2.eq.'W+ '.and.p3.eq.'W+~')return
      if(p1.eq.'Z0 '.and.p2.eq.'W+ '.and.p3.eq.'W+~')return
      if(p1.eq.'H0 '.and.p2.eq.'H0 '.and.p3.eq.'H0 ')return
      if(p1.eq.'H0 '.and.p2.eq.'W+ '.and.p3.eq.'W+~')return
      if(p1.eq.'H0 '.and.p2.eq.'Z0 '.and.p3.eq.'Z0 ')return
      if(lightfermions.eq.1)then
        if(p1.eq.'el '.and.p2.eq.'el~'.and.p3.eq.'H0 ')return
        if(p1.eq.'mu '.and.p2.eq.'mu~'.and.p3.eq.'H0 ')return
        if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'H0 ')return
        if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'H0 ')return
        if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'H0 ')return
        if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'H0 ')return
      endif
      if(lightfermions.eq.1.or.lightfermions.eq.2)then
        if(p1.eq.'ta '.and.p2.eq.'ta~'.and.p3.eq.'H0 ')return
        if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'H0 ')return
      endif
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'H0 ')return
! ew four-particle vertices (auxiliary particles v1,v2)
      if(p1.eq.'ga '.and.p2.eq.'ga '.and.p3.eq.'v1 ')return
      if(p1.eq.'ga '.and.p2.eq.'Z0 '.and.p3.eq.'v1 ')return
      if(p1.eq.'Z0 '.and.p2.eq.'Z0 '.and.p3.eq.'v1 ')return
      if(p1.eq.'W+ '.and.p2.eq.'W+~'.and.p3.eq.'v1 ')return
      if(p1.eq.'W+ '.and.p2.eq.'W+~'.and.p3.eq.'v1~')return
      if(p1.eq.'Z0 '.and.p2.eq.'Z0 '.and.p3.eq.'v2 ')return
      if(p1.eq.'W+ '.and.p2.eq.'W+~'.and.p3.eq.'v2 ')return
      if(p1.eq.'H0 '.and.p2.eq.'H0 '.and.p3.eq.'v2 ')return
      if(p1.eq.'H0 '.and.p2.eq.'H0 '.and.p3.eq.'v2~')return
! QCD vertices
      if(p1.eq.'dq '.and.p2.eq.'dq~'.and.p3.eq.'gl ')return
      if(p1.eq.'uq '.and.p2.eq.'uq~'.and.p3.eq.'gl ')return
      if(p1.eq.'sq '.and.p2.eq.'sq~'.and.p3.eq.'gl ')return
      if(p1.eq.'cq '.and.p2.eq.'cq~'.and.p3.eq.'gl ')return
      if(p1.eq.'bq '.and.p2.eq.'bq~'.and.p3.eq.'gl ')return
      if(p1.eq.'tq '.and.p2.eq.'tq~'.and.p3.eq.'gl ')return
!     if(p1.eq.'gl '.and.p2.eq.'gl '.and.p3.eq.'gl ')return
! QCD four-gluon vertex (auxiliary particle v3)
      if(p1.eq.'gl '.and.p2.eq.'gl '.and.p3.eq.'v3 ')return
      if(p1.eq.'gl '.and.p2.eq.'gl '.and.p3.eq.'v3~')return
      lusifer_vertexg=.false.
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                 c
!     identifying particles                                       c
!                                                                 c
!     written by Markus Roth                                      c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lusifer_id(binary)
      implicit none
! local variables
      integer lusifer_id,binary,i1,i2,i3,nexternal
      i3=1
      nexternal=30
      lusifer_id=0
      do i1=nexternal,1,-1
        i2=binary/2**(i1-1)
        if((i2/2)*2.ne.i2)then
          lusifer_id=lusifer_id+i1*i3
          i3=10*i3
        endif
      enddo
      end
