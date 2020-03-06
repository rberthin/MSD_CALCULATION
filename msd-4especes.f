PROGRAM MSD

! Last altered 14.3.01

   IMPLICIT DOUBLE PRECISION(a-h,o-z)
   PARAMETER(nmsdcorrmax=10000)
   PARAMETER(nummax=4000)
   PARAMETER(nspmax=4)

   COMMON/msdarray/xdisp(nummax),xdispstore(nummax,nmsdcorrmax),
                 & ydisp(nummax),ydispstore(nummax,nmsdcorrmax),
                 & zdisp(nummax),zdispstore(nummax,nmsdcorrmax)
   COMMON/msdcorr/xmsd(nummax,0:nmsdcorrmax),
                & ymsd(nummax,0:nmsdcorrmax),
                & zmsd(nummax,0:nmsdcorrmax),
                & norm(nummax,0:nmsdcorrmax)
   COMMON/conduct/xds(0:nmsdcorrmax,nspmax,nspmax),
                & yds(0:nmsdcorrmax,nspmax,nspmax),
                & zds(0:nmsdcorrmax,nspmax,nspmax),
                & normtot(0:nmsdcorrmax )
   COMMON/msdint/nmsdlength,nmsdcalltime,mcorrtime
   COMMON/typsp/ntype(nummax),numspc(nspmax)
   COMMON/sizeof/num
   COMMON/step/nstep, nrun
   COMMON/over/overflow 
   COMMON/timest/dtime
   COMMON/ntype/ncation1, ncation2, ncation3, nanion
   COMMON/ntype2/z(nspmax), nspecies
   COMMON/files/rstfile, filename

   CHARACTER(LEN=200) :: filein
   CHARACTER(LEN=20) :: rstfile, filename
   rstfile='msdrst.dat'
 
   LOGICAL overflow
   LOGICAL restart
   LOGICAL endrun
   LOGICAL readfrominpt_log
   restart=.false.
   overflow=.false.
   endrun=.false.
   readfrominpt_log=.false.
   REAL bin
   INTEGER nbin

!----------------
! GET PARAMETERS
!----------------

      OPEN(10, FILE='msd.inpt')
      READ(10,'(a)') filein
      OPEN(11, FILE=filein, STATUS='old', FORM='formatted')
      READ(10,*) num
      READ(10,*) nanion
      READ(10,*) ncation1
      READ(10,*) ncation2
      READ(10,*) ncation3
      READ(10,*) nmsdlength
      READ(10,*) nspecies

      IF(nspecies.gt.nspmax.or.num.gt.nummax.or.nmsdlength.gt.nmsdcorrmax) THEN
          WRITE(6,*) '***Problems with array dimensions***'
          STOP
      ENDIF

      DO 10 ns = 1, nspecies
         READ(10,*) z(ns)
10    CONTINUE
 
      DO 12 ns = 1, nspecies
          filename = 'msd'//char(ns+48)//'.dat'
          OPEN (52 + ns, FILE=filename)
          DO 11 ns2 = ns, nspecies
              filename = 'msdcollect'//char(ns+48)//char(ns2+48)//'.dat'
              OPEN (70+(ns*nspecies)+ns2,FILE=filename)
11        CONTINUE
12    CONTINUE
      OPEN (98, FILE='work.dat')
      OPEN (99, FILE='nernst.dat')

!****************
! Zero arrays
!****************

      DO 20 n = 1, num
         xdisp(n) = 0.0d0
         ydisp(n) = 0.0d0
         xdisp(n) = 0.0d0

         DO 15 i = 1, nmsdcorrmax
            xdispstore(n,i) = 0.0d0
            ydispstore(n,i) = 0.0d0
            zdispstore(n,i) = 0.0d0
15       CONTINUE
20    CONTINUE

      DO 40 n=1,num
         DO 30 i = 0, nmsdcorrmax
            xmsd(n,i) = 0.0d0
            ymsd(n,i) = 0.0d0
            zmsd(n,i) = 0.0d0
            norm(n,i) = 0
30       CONTINUE
40    CONTINUE

      DO 50 i = 0, nmsdcorrmax
         DO 42 m = 1, nspecies
            DO 41 l = 1, nspecies
               xds(i,m,l) = 0.0d0
               yds(i,m,l) = 0.0d0
               zds(i,m,l) = 0.0d0
41          CONTINUE
42       CONTINUE
          normtot(i) = 0
50    CONTINUE
 
      read(10,*) restart

      if (restart) then

          open(12,file=rstfile,status='OLD',form='formatted')

          do 70 i=0,nmsdlength
              do 60 n=1,num

                  read(12,*) xmsd(n,i)
                  read(12,*) ymsd(n,i)
                  read(12,*) zmsd(n,i)
                  read(12,*) norm(n,i)

60            continue

              do 65 i1=1,nspecies
                  do 64 i2=1,nspecies

                      read(12,*) xds(i,i1,i2)
                      read(12,*) yds(i,i1,i2)
                      read(12,*) zds(i,i1,i2)

64                continue
65            continue

              read(12,*) normtot(i)

70        CONTINUE
          CLOSE(12)
      ENDIF

      read(10,*) readfrominpt_log

      if (readfrominpt_log) then

          read (10,*) nmsdcalltime
          read (10,*) dtime
          read (10,*) nrun

          do 80 n=1,nanion
              ntype(n)=1
80        continue

          if (ncation1.ne.0) then
              do 90 n=nanion+1,nanion+ncation1
                  ntype(n)=2
90            continue
          endif

          if (ncation2.ne.0) then
              do 100 n=nanion+ncation1+1,num
                  ntype(n)=3
100           continue
          endif

          if (ncation3.ne.0) then
              do 101 n=nanion+ncation1+ncation2+1,num
                  ntype(n)=4
101           continue
          endif

!          read (11,*) nbin,bin,nbin
!          do 105 i=1,num
!              read(11,*) nbin
!105       continue

!      else
!
!          read (11,*) nmsdcalltime,dtime,nrun
!
!          do 110 i=1,num
!              read(11,*) ntype(i)
!110       continue
!
      endif

      do 120 ispec=1,nspecies
          numspc(ispec)=0
120   continue
      
      do 130 n=1,num
          ispec = ntype(n)
          numspc(ispec) = numspc(ispec)+1
130   continue
      
      close(10)


      nstep=0
      mcorrtime=1

200   nstep=nstep+nmsdcalltime
c     write(6,*) nstep

      if (nstep.ge.nrun) endrun=.true.

      do 232 i=1,num

          read(11,*) xdisp(i),ydisp(i),zdisp(i)
 
232   continue

      print *,'entering msdcalc',nstep
      call msdcalc

      if (endrun) goto 500

      goto 200

500   close(11)

c=========================================================
c Write out restart file
c=========================================================

       call msdoutput

       do 600 ns=1,nspecies
           close(102+ns)
600    continue

       stop
       end


      subroutine msdcalc

      implicit double precision (a-h,o-z)

      parameter (nmsdcorrmax=10000)
      parameter (nummax=4000)
      parameter (nspmax=4)

      common/msdarray/xdisp(nummax),xdispstore(nummax,nmsdcorrmax),
     x                ydisp(nummax),ydispstore(nummax,nmsdcorrmax),
     x                zdisp(nummax),zdispstore(nummax,nmsdcorrmax)
      common/msdcorr/xmsd(nummax,0:nmsdcorrmax),
     x               ymsd(nummax,0:nmsdcorrmax),
     x               zmsd(nummax,0:nmsdcorrmax),
     x               norm(nummax,0:nmsdcorrmax)
      common/conduct/ xds(0:nmsdcorrmax,nspmax,nspmax),
     x                yds(0:nmsdcorrmax,nspmax,nspmax),
     x                zds(0:nmsdcorrmax,nspmax,nspmax),
     x                normtot(0:nmsdcorrmax )
      common/msdint/nmsdlength,nmsdcalltime,mcorrtime
      common/typsp/ntype(nummax),numspc(nspmax)
      common/sizeof/num
      common/over/overflow
      common/ntype/ncation1,ncation2,ncation3,nanion
      common/ntype2/z(nspmax),nspecies

      dimension xtot(0:nmsdcorrmax,nspecies), 
     x          ytot(0:nmsdcorrmax,nspecies),
     x          ztot(0:nmsdcorrmax,nspecies)
      logical overflow


      do 10 i=1,num

          xdispstore(i,mcorrtime)=0.0d0
          ydispstore(i,mcorrtime)=0.0d0
          zdispstore(i,mcorrtime)=0.0d0

10    continue

c accumulate displacements

      do 30 j=1,mcorrtime
          do 20 i=1,num

              xdispstore(i,j)=xdispstore(i,j)+xdisp(i)
              ydispstore(i,j)=ydispstore(i,j)+ydisp(i)
              zdispstore(i,j)=zdispstore(i,j)+zdisp(i)


20        continue

30    continue

      if (overflow) then

          do 80 j=mcorrtime+1,nmsdlength

              do 70 i=1,num

                  xdispstore(i,j)=xdispstore(i,j)+xdisp(i)
                  ydispstore(i,j)=ydispstore(i,j)+ydisp(i)
                  zdispstore(i,j)=zdispstore(i,j)+zdisp(i)

70            continue

80        continue

      endif

c zero displacement arrays

      do 90 i=1,num
          xdisp(i)=0.0d0
          ydisp(i)=0.0d0
          zdisp(i)=0.0d0
90    continue

 
      do 200 j=1,mcorrtime

          nt=mcorrtime-j
          do 110 ipt=1,nspecies
               xtot(j,ipt)=0.0d0
               ytot(j,ipt)=0.0d0
               ztot(j,ipt)=0.0d0
110       continue

          do 120 i=1,num
              ipoint=ntype(i)
                
              xtot(j,ipoint)= xtot(j,ipoint) + xdispstore(i,j)
              ytot(j,ipoint)= ytot(j,ipoint) + ydispstore(i,j)
              ztot(j,ipoint)= ztot(j,ipoint) + zdispstore(i,j)


              norm(i,nt)=norm(i,nt)+1
              xmsd(i,nt)=xmsd(i,nt)+xdispstore(i,j)**2
              ymsd(i,nt)=ymsd(i,nt)+ydispstore(i,j)**2
              zmsd(i,nt)=zmsd(i,nt)+zdispstore(i,j)**2
120       continue
          normtot(nt)=normtot(nt)+1
          do 190 i1=1,nspecies
              do 180 i2=i1,nspecies
  
                  xds(nt,i1,i2)=xds(nt,i1,i2)+xtot(j,i1)*xtot(j,i2)
                  yds(nt,i1,i2)=yds(nt,i1,i2)+ytot(j,i1)*ytot(j,i2)
                  zds(nt,i1,i2)=zds(nt,i1,i2)+ztot(j,i1)*ztot(j,i2)
180           continue
190       continue

200   continue

 
      if (overflow) then
 
          do 60 j=mcorrtime+1,nmsdlength

              nt=mcorrtime-j+nmsdlength

              do 27 ipt=1,nspecies
                  xtot(j,ipt)=0.0d0
                  ytot(j,ipt)=0.0d0
                  ztot(j,ipt)=0.0d0
27            continue

              do 50 i=1,num

                  ipoint=ntype(i)

                  xtot(j,ipoint)= xtot(j,ipoint) + xdispstore(i,j)
                  ytot(j,ipoint)= ytot(j,ipoint) + ydispstore(i,j)
                  ztot(j,ipoint)= ztot(j,ipoint) + zdispstore(i,j)

                  norm(i,nt)=norm(i,nt)+1
                  xmsd(i,nt)=xmsd(i,nt)+xdispstore(i,j)**2
                  ymsd(i,nt)=ymsd(i,nt)+ydispstore(i,j)**2
                  zmsd(i,nt)=zmsd(i,nt)+zdispstore(i,j)**2

50            continue
 
              normtot(nt)=normtot(nt)+1
 
              do 57 i1=1,nspecies
                  do 55 i2=i1,nspecies
  
                      xds(nt,i1,i2)=xds(nt,i1,i2)+ xtot(j,i1)*xtot(j,i2)
                      yds(nt,i1,i2)=yds(nt,i1,i2)+ ytot(j,i1)*ytot(j,i2)
                      zds(nt,i1,i2)=zds(nt,i1,i2)+ ztot(j,i1)*ztot(j,i2)
55                continue
57            continue

60        continue


      endif

c========================================================================
c Update array counters
c========================================================================
      if (mod(float(mcorrtime),float(nmsdlength)).eq.0) then
          overflow=.true.
      endif

      mcorrtime=int(mod(float(mcorrtime),float(nmsdlength)))
      mcorrtime=mcorrtime+1

      return
      end


      subroutine msdoutput

      implicit double precision(a-h,o-z)

      parameter (nmsdcorrmax=10000)
      parameter (nummax=4000)
      parameter (nspmax=4)
  

  
      common/conduct/ xds(0:nmsdcorrmax,nspmax,nspmax),
     x                yds(0:nmsdcorrmax,nspmax,nspmax),
     x                zds(0:nmsdcorrmax,nspmax,nspmax),
     x                normtot(0:nmsdcorrmax )

      common/msdcorr/xmsd(nummax,0:nmsdcorrmax),
     x               ymsd(nummax,0:nmsdcorrmax),
     x               zmsd(nummax,0:nmsdcorrmax),
     x               norm(nummax,0:nmsdcorrmax)
      common/msdint/nmsdlength,nmsdcalltime,mcorrtime
      common/typsp/ntype(nummax),numspc(nspmax)
      common/sizeof/num
      common/step/nstep,nrun
      common/timest/dtime
      common/ntype/ncation1,ncation2,ncation3,nanion
      common/ntype2/z(nspmax),nspecies
      common/files/rstfile,filename

      double precision msd(0:nmsdcorrmax,nspecies),
     x          msdcoll(0:nmsdcorrmax,nspecies,nspecies)
      double precision sum(nspecies)

      character*20 rstfile,filename

      do 30 i=0,nmsdlength
          do 20 ip=1,nspecies
              msd(i,ip)=0.0d0
              do 10 ip2=1,nspecies
                  msdcoll(i,ip,ip2)=0.0d0
10             continue
20        continue
30    continue
  
c   Writing out values for restart file

      open (14,file=rstfile,form='formatted')

      do 60 i=0,nmsdlength

          do 40 n=1,num

              write(14,*) xmsd(n,i)
              write(14,*) ymsd(n,i)
              write(14,*) zmsd(n,i)
              write(14,*) norm(n,i)

40        continue

          do 55 i1=1,nspecies
              do 50 i2=1,nspecies

                  write(14,*) xds(i,i1,i2)
                  write(14,*) yds(i,i1,i2)
                  write(14,*) zds(i,i1,i2)

50            continue
55        continue

          write(14,*) normtot(i)

60    continue

      close(14)

c======================================================================
c Average over components and number of molecules
c======================================================================

      do 80 i=1,num
          do 70 j=0,nmsdlength
              if (norm(i,j ).gt.0) then
                  xmsd(i,j)=xmsd(i,j)/float(norm(i,j))
                  ymsd(i,j)=ymsd(i,j)/float(norm(i,j))
                  zmsd(i,j)=zmsd(i,j)/float(norm(i,j))
              endif
70        continue

80    continue

      do 120 j=0,nmsdlength

          do 110 n1=1,nspecies

              do 100 n2=1,nspecies

                  if (normtot(j).gt.0) then
                      rn=sqrt(float(numspc(n1)*numspc(n2)))*
     &                      float(normtot(j))

                      xds(j,n1,n2)=xds(j,n1,n2)/rn
                      yds(j,n1,n2)=yds(j,n1,n2)/rn
                      zds(j,n1,n2)=zds(j,n1,n2)/rn
                      msdcoll(j,n1,n2)=xds(j,n1,n2)+
     &                    yds(j,n1,n2)+zds(j,n1,n2)

                  endif
100           continue

110       continue
             
120   continue


      do 150 i=0,nmsdlength

          do 130 ip=1,nspecies
              sum(ip)=0.0d0
130       continue
 

          do 140 j=1,num

              ip=ntype(j)
              sum(ip)=sum(ip) +xmsd(j,i)+ymsd(j,i)+zmsd(j,i)
 
140       continue

          do 145 ip=1,nspecies
              if (numspc(ip).gt.0) then
                  msd(i,ip)=sum(ip)/float(numspc(ip))
              endif
145       continue

150   continue

          dum1=0.0
          dum2=0.0

          do 161 ip=1,nspecies
              write (52+ip,*) dum1,dum2
161       continue

          write(98,*)dum1,dum2
          write(99,*)dum1,dum2


          do 181 ip1=1,nspecies

              if (numspc(ip1).gt.0) then
                  do 171 ip2=ip1,nspecies
                      write(70+(ip1*nspecies)+ip2,*) dum1,dum2
171               continue
              endif
181       continue


      do 200 i=0,nmsdlength-1
          time=(dble(i)+1)*dble(nmsdcalltime)*dtime*2.418d-5

          do 160 ip=1,nspecies
              write (52+ip,*) time,msd(i,ip)
160       continue

          work=0.0d0
          wnernst=0.0d0

          do 180 ip1=1,nspecies

              if (numspc(ip1).gt.0) then
                  do 170 ip2=ip1,nspecies
                      write(70+(ip1*nspecies)+ip2,*) time,
     &                               msdcoll(i,ip1,ip2)
                      msdcoll(i,ip1,ip2)=msdcoll(i,ip1,ip2)*
     &                                   z(ip1)*z(ip2)*
     &                sqrt(float(numspc(ip1)*numspc(ip2)))/float(num)
	              if (ip1.ne.ip2) then 
	                  work=work+msdcoll(i,ip1,ip2)
	              endif
	              work=work+msdcoll(i,ip1,ip2)
170               continue
              endif

          wnernst=wnernst+msd(i,ip1)*z(ip1)*z(ip1)
     &           *float(numspc(ip1))/float(num)

180       continue
  
          WRITE(98,*)time,work
          WRITE(99,*)time,wnernst

200   CONTINUE
  
      WRITE(6,*)
      WRITE(6,*) '*** Mean squared displacements written out. ***'
      WRITE(6,*)

      return
END
