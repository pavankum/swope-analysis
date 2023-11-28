c
c     program to compute rmsd between pairs of conformers
c     coordinates for the conformers come from two sdf files
c
c     Swope, September 2022
c
c     This program was produced by modification of getrmsd.f
c     which is in the DirtyDozen/RMSD directory
c     It is modified here to work with the OpenFF Public dataset.
c
c     Program is driven from a compound list that should have
c     a superset of all compounds being considered
c     Generate it from inside one of the directories with the sdf with 
c     ls *.sdf | cut -c 1-9 | sort | uniq > compound.list
c     should produce an alphabetized list
c
      implicit real*8 (a-h,o-z)
      character*1024 dir1,dir2,dir3
      character*1024 line,filename
      character*1024 sdffile1(10),sdffile2(10)
      dimension energy1(10),energy2(10)
      dimension rmsdarray(10,10)
      logical ifexist

      character*1024 confname1,confname2
      character*2 elname1(200),elname2(200)
      dimension coord1(3,200),coord2(3,200)
      dimension cent1(3),cent2(3)
      dimension w(3,3),scratch(3,3)
      logical ifexist2
c     needed for lapack
      dimension a(3,3),s(3,3),u(3,3),vt(3,3),work(100)
      dimension mates(0:9,2)

      maxatom=200  !size of coord arrays

c     sdffile1='qm-filtered.sdf'
c     sdffile2='openff-2.0.0.sdf'
c     confname1='GNT-00023-02'
c     confname2='GNT-00023-00'


      write(*,'(/,''Program to compute RMSD values '',
     x  ''between pairs of sdf files,'',/,
     x  ''one from each of two different directories.'',/)')
c     call getarg(1,dir1)
c     call getarg(2,dir2)
      dir1='./b3lyp-d3bj_dzvp/'

      dir2='./MM-geometries/'
      dir3='./RMSD-OpenFF/'

c      dir2='./ffld_opls4_def/'
c      dir3='./RMSD-OPLS4-def/'
c
c      dir2='./ffld_opls4_cst/'
c      dir3='./RMSD-OPLS4-cst/'

c     redo the proprietary set
c      dir1='../../4-compute-qm-filtered/b3lyp-d3bj/dzvp/'
c      dir2='../../4-compute-mm-filtered/openff-2.0.0-rc.2/'
c      dir3='./RMSD-OpenFF/'
c
c      dir2='../../4-compute-mm-filtered/opls4_custom/'
c      dir3='./RMSD-OPLS4-cst/'


      is1=lstchr(dir1)
      is2=lstchr(dir2)
      is3=lstchr(dir3)
      write(*,'(''Directory 1: '',a)') dir1(1:is1)
      write(*,'(''Directory 2: '',a)') dir2(1:is2)
      write(*,'(''Directory 3: '',a)') dir3(1:is3)

      open(10,file='compound.list')

5     continue
      read(10,'(a1024)',end=500) line
      is4=index(line,' ')-1
      write(*,'(''Compound '',a)') line(1:is4)
      niconf=0
      njconf=0
      do iconf=1,10 
        filename=dir1(1:is1)//line(1:is4)//'-  .sdf'
        write(filename(is1+is4+2:is1+is4+3),'(i2.2)') iconf-1
        write(*,'(''Looking for file *'',a,''*'')') 
     x    filename(1:is1+is4+7)
        inquire(file=filename(1:is1+is4+7),exist=ifexist)
        if(ifexist) then
          niconf=niconf+1
          sdffile1(niconf)=filename(1:is1+is4+7)
        endif
        filename=dir2(1:is2)//line(1:is4)//'-  '//'.sdf'
        write(filename(is2+is4+2:is2+is4+3),'(i2.2)') iconf-1
c       write(*,'(''Looking for file *'',a,''*'')') 
c    x    filename(1:is2+is4+7)
        inquire(file=filename(1:is2+is4+7),exist=ifexist)
        if(ifexist) then
          njconf=njconf+1
          sdffile2(njconf)=filename(1:is2+is4+7)
        endif
      enddo

      do iconf=1,niconf
        do jconf=1,njconf
          call getrmsd(sdffile1(iconf),sdffile2(jconf),rmsd,en1,en2)
          rmsdarray(iconf,jconf)=rmsd
          energy1(iconf)=en1
          energy2(jconf)=en2
        enddo
      enddo


      open(12,file=dir3(1:is3)//line(1:is4)//'.rmsd')
      write(12,'(i2)') niconf
      do i=1,niconf
        write(12,'(a)') sdffile1(i)(1:is1+is4+7)
        write(12,'(f25.10)') energy1(i)
      enddo
      write(12,'(i2)') njconf
      do i=1,njconf
        write(12,'(a)') sdffile2(i)(1:is2+is4+7)
        write(12,'(f25.10)') energy2(i)
      enddo
      do i=1,niconf
        write(12,'(10f8.4)') (rmsdarray(i,j),j=1,njconf)
      enddo
      close(12)

      go to 5

500   continue







c     do icomp=0,2000
c       do im=0,9
c         mates(im,1)=0
c         mates(im,2)=0
c       enddo
c       do iconf=0,9
c         write(confname1,'(a,''/GNT-'',i5.5,''-'',i2.2,
c    x    ''.sdf'')')
c    x    dir1(1:is),icomp,iconf
c         write(*,'(''Looking for file '',a)') 
c    x    confname1(1:lstchr(confname1))
c       inquire(file=confname1(1:lstchr(confname1)),exist=ifexist)
c       if(ifexist) then
c         mates(iconf,1)=1
c         write(*,'(''File found: '',a)') 
c    x      confname1(1:lstchr(confname1))    
c         look for any files in directory 2 that have the same
c         compound number
c         do jconf=0,9
c           write(confname2,'(a,''/GNT-'',i5.5,''-'',i2.2,
c    x        ''.sdf'')')
c    x        dir2(1:is2),icomp,jconf
c           inquire(file=confname2(1:lstchr(confname2)),
c    x        exist=ifexist2)
c           if(ifexist2) then
c             mates(jconf,2)=1
c             write(*,'(''File2 found:   '',a)') 
c    x        confname2(1:lstchr(confname2))    

c             call getrmsd(confname1,confname2,rmsd)
c             isx=index(confname1,'.sdf')
c             isy=index(confname2,'.sdf')
c             write(*,'(a,5x,a,5x,f12.8)')
c    x          confname1(isx-12:isx-1), 
c    x          confname2(isy-12:isy-1), 
c    x          rmsd

c           endif
c         enddo !end of loop over conformers in the second directory
c       endif
c       enddo  !end of loop over conformers in the first directory
c       do im=0,9
c         if(mates(im,1).ne.mates(im,2)) then
c           if(mates(im,1).eq.1) write(*,'(''Warning: '',
c    x        i5.5,''-'',i2.2,
c    x        '' exists in the first directory, '',
c    x        ''but not in the second.'')')
c    x        icomp,im
c           if(mates(im,2).eq.1) write(*,'(''Warning: '',
c    x        i5.5,''-'',i2.2,
c    x        '' exists in second directory, '',
c    x        ''but not in the first.'')')
c    x        icomp,im
c         endif
c       enddo
c     enddo   !end of loop over compounds 




      end



      subroutine getrmsd(sdffile1,sdffile2,rmsd,energy1,energy2)
      implicit real*8 (a-h,o-z)
      character*(*) sdffile1,sdffile2
      character*1024 line
      character*1024 confname1,confname2
      character*1024 dir1,dir2
      character*2 elname1(500),elname2(500)
      dimension coord1(3,500),coord2(3,500)
      dimension cent1(3),cent2(3)
      dimension w(3,3),scratch(3,3)
c     needed for lapack
      dimension a(3,3),s(3,3),u(3,3),vt(3,3),work(100)

      maxatom=500  !size of coord arrays


      open(20,file=sdffile1)
      open(21,file=sdffile2)

      is1=index(sdffile1,'.sdf')
      confname1=sdffile1(is1-12:is1-1)
      is2=index(sdffile2,'.sdf')
      confname2=sdffile2(is2-12:is2-1)



c     get coord set 1
10    read(20,'(a1024)',end=100) line
      if(index(line,confname1(1:lstchr(confname1)))
     x   .eq.0) go to 10
c     write(*,'(''First conformer found'')')
      read(20,'(a1024)') line
      read(20,'(a1024)') line
      read(20,'(a1024)') line
      read(line,'(i3)') nat1
c     write(*,'(''Number atoms in first conformer '',
c    x  i5)') nat1
      if(nat1.gt.maxatom) then
        write(*,'(''Increase maxatom (now '',i3,
     x    '') to at least '',i3)') maxatom,nat1
        stop 'Increase storage, too many atoms'
      endif
      do i=1,nat1
        read(20,'(a1024)') line
c       write(*,'(a)') line(1:lstchr(line))
        read(line,'(3F10.4,1x,a2)')  
     x    (coord1(j,i),j=1,3),elname1(i)
c       write(*,'(''Site '',i2,'' Elname '',a2,5x,3f10.4)')
c    x    i,elname1(i),(coord1(j,i),j=1,3)
      enddo
15    continue
      read(20,'(a1024)') line
      if(index(line,'final_energy').eq.0) go to 15
      read(20,*) energy1

100   continue


c     get coord set 2
20    read(21,'(a1024)',end=200) line
      if(index(line,confname2(1:lstchr(confname2)))
     x   .eq.0) go to 20
c     write(*,'(''Second conformer found'')')
      read(21,'(a1024)') line
      read(21,'(a1024)') line
      read(21,'(a1024)') line
      read(line,'(i3)') nat2
c     write(*,'(''Number atoms in second conformer '',
c    x  i5)') nat2
      do i=1,nat2
        read(21,'(a1024)') line
c       write(*,'(a)') line(1:lstchr(line))
        read(line,'(3F10.4,1x,a2)')  
     x    (coord2(j,i),j=1,3),elname2(i)
c       write(*,'(''Site '',i2,'' Elname '',a2,5x,3f10.4)')
c    x    i,elname2(i),(coord2(j,i),j=1,3)
      enddo
25    continue
      read(21,'(a1024)') line
      if(index(line,'final_energy').eq.0) go to 25
      read(21,*) energy2


200   continue
      if(nat1.ne.nat2) then
        write(*,'(''Error:  conformers have '',
     x     ''different numbers of atoms.'')')
        stop 'error - conformers have different number atoms'
      endif

      close(20)
      close(21)

c     compute centroids
      do i=1,3
        cent1(i)=0.d0
        cent2(i)=0.d0
      enddo

      do i=1,nat1
        do j=1,3
          cent1(j)=cent1(j)+coord1(j,i)/nat1
          cent2(j)=cent2(j)+coord2(j,i)/nat1
        enddo
      enddo

c     compute w
c     (if we want to ignore hydrogens in the rmsd computation,
c     here is a place to do it)

      do i=1,3
        do j=1,3
          dot=0.d0
          do k=1,nat1
            dot=dot+(coord1(i,k)-cent1(i))
     x             *(coord2(j,k)-cent2(j))
          enddo
c         write(*,'(''w('',i1,'','',i1,'')='',f10.4)') 
c    x      i,j,dot
          w(i,j)=dot  
        enddo
      enddo

c     write(*,'(''A=W:'')') 
      do i=1,3
        do j=1,3
          a(i,j)=w(i,j)
          s(i,j)=0.d0
          u(i,j)=0.d0
          vt(i,j)=0.d0
        enddo
c       write(*,'(3f10.4)') (w(i,j),j=1,3)
      enddo

      n=3
      m=3
      lda=3
      ldu=3
      ldvt=3
      lwork=100
      call dgesvd('A',   ! get m-by-m matrix U
     x            'A',   ! get n-by-n matrix transpose(V)
     x            m,n,a, ! m and n and m-by-n matrix a
     x            lda,
     x            s,
     x            u,ldu,vt,ldvt, !m-by-m matrix u; n-by-n matrix vT
     x            work,lwork,info)
c     write(*,'(''On output from dgesvd, info= '',i5,
c    x  '' singular values '',3f10.4)') 
c    x  info,s(1,1),s(2,1),s(3,1)

c     write(*,'(''On output from dgesvd, matrix u'')')
c     do i=1,3
c       write(*,'(3f10.4)') (u(i,j),j=1,3)
c     enddo

c     check column orthogonality of u
c     write(*,'(''Matrix U column orthog check'')')
c     do i=1,3
c       do j=1,3
c         dot=0.d0
c         do k=1,3
c           dot=dot+u(k,i)*u(k,j) 
c         enddo
c         scratch(i,j)=dot
c       enddo
c       write(*,'(3f10.4)') (scratch(i,j),j=1,3)
c     enddo

c     check row orthogonality of u
c     write(*,'(''Matrix U row orthog check'')')
c     do i=1,3
c       do j=1,3
c         dot=0.d0
c         do k=1,3
c           dot=dot+u(i,k)*u(j,k) 
c         enddo
c         scratch(i,j)=dot
c       enddo
c       write(*,'(3f10.4)') (scratch(i,j),j=1,3)
c     enddo





c     write(*,'(''On output from dgesvd, matrix vT'')')
c     do i=1,3
c       write(*,'(3f10.4)') (vt(i,j),j=1,3)
c     enddo

c     check column orthogonality of vT
c     write(*,'(''Matrix vT column orthog check'')')
c     do i=1,3
c       do j=1,3
c         dot=0.d0
c         do k=1,3
c           dot=dot+vt(k,i)*vt(k,j) 
c         enddo
c         scratch(i,j)=dot
c       enddo
c       write(*,'(3f10.4)') (scratch(i,j),j=1,3)
c     enddo

c     check row orthogonality of vt
c     write(*,'(''Matrix VT row orthog check'')')
c     do i=1,3
c       do j=1,3
c         dot=0.d0
c         do k=1,3
c           dot=dot+vt(i,k)*vt(j,k) 
c         enddo
c         scratch(i,j)=dot
c       enddo
c       write(*,'(3f10.4)') (scratch(i,j),j=1,3)
c     enddo



c     write(*,'(''On output from dgesvd, matrix s'')')
c     do i=1,3
c       write(*,'(3f10.4)') (s(i,j),j=1,3)
c     enddo


c     check that A = U S Vt
c     get scratch = S * Vt
c     do i=1,3
c       do j=1,3
c         scratch(i,j)=s(i,1)*vt(i,j)
c       enddo
c     enddo
c     get a = U * scratch
c     do i=1,3
c       do j=1,3
c         dot=0.d0
c         do k=1,3
c          dot=dot+u(i,k)*scratch(k,j)
c         enddo
c         a(i,j)=dot
c       enddo
c     enddo
c     compare a and w
c     do i=1,3
c       do j=1,3
c         write(*,'(2i2,'' W: '',f10.4,
c    x    '' A:'',f10.4,'' Dif: '',f10.4)')
c    x    i,j,w(i,j),a(i,j),w(i,j)-a(i,j)
c       enddo
c     enddo


c
c     compute det (V * U^T)
c
c     write(*,'(''Matrix V*Ut:'')')
      do i=1,3
        do j=1,3
          dot=0.d0
          do k=1,3
           dot=dot+vt(k,i)*u(j,k)
          enddo
          scratch(i,j)=dot
        enddo
c       write(*,'(3f10.4)') (scratch(i,j),j=1,3)
      enddo
      det=0.d0
      det = det + scratch(1,1)*
     x    (scratch(2,2)*scratch(3,3)
     x    -scratch(3,2)*scratch(2,3))
      det = det - scratch(1,2)*
     x    (scratch(2,1)*scratch(3,3)
     x    -scratch(3,1)*scratch(2,3))
      det = det + scratch(1,3)*
     x    (scratch(2,1)*scratch(3,2)
     x    -scratch(3,1)*scratch(2,2))
c     write(*,'(''Determinant = '',f10.4)') det

      do i=1,3
        do j=1,3
          s(i,j)=0.d0
        enddo
        s(i,i)=1.d0
      enddo
      if(det.le.0.d0) s(3,3)=-1.d0


c     construct rotation matrix R = V S UT
c     first scratch = S UT
      do i=1,3
        do j=1,3
          scratch(i,j)=s(i,i)*u(j,i)
        enddo 
      enddo
c     compute rotation matrix (in A) R = V * scratch
      do i=1,3
        do j=1,3
          dot=0.d0
          do k=1,3
            dot=dot + vt(k,i)*scratch(k,j)
          enddo
          a(i,j)=dot
        enddo 
      enddo
c     write(*,'(''Rotation matrix '')') 
c     do i=1,3
c       write(*,'(3f10.4)') (a(i,j),j=1,3)
c     enddo

c     write(*,'(''R * RT should be I '')') 
c     do i=1,3
c       do j=1,3
c         dot=0.d0
c         do k=1,3
c           dot=dot + a(i,k)*a(j,k)
c         enddo
c         scratch(i,j)=dot
c       enddo 
c     enddo
c     do i=1,3
c       write(*,'(3f10.4)') (scratch(i,j),j=1,3)
c     enddo
      
c
c     given the rotation matrix, compute rmsd
c     msd = sum over atoms n of R*(P-<P>)_k * (Q-<Q>)_k
      rmsd=0.d0
      do n=1,nat1
c       x comp:
        dx= a(1,1)*(coord1(1,n)-cent1(1)) 
     x     +a(1,2)*(coord1(2,n)-cent1(2))
     x     +a(1,3)*(coord1(3,n)-cent1(3)) 
     x     -      ( coord2(1,n)-cent2(1)) 
c       y comp:
        dy= a(2,1)*(coord1(1,n)-cent1(1)) 
     x     +a(2,2)*(coord1(2,n)-cent1(2))
     x     +a(2,3)*(coord1(3,n)-cent1(3)) 
     x     -      ( coord2(2,n)-cent2(2)) 

        dz= a(3,1)*(coord1(1,n)-cent1(1)) 
     x     +a(3,2)*(coord1(2,n)-cent1(2))
     x     +a(3,3)*(coord1(3,n)-cent1(3)) 
     x     -      ( coord2(3,n)-cent2(3)) 

        dist=dsqrt(dx**2+dy**2+dz**2)
c       write(*,'(''Atom '',i2,'' displacement '',
c    x    3f10.4,'' dist '',f10.4)')
c    x    n,dx,dy,dz,dist
        rmsd=rmsd+dist**2

      enddo
     
      rmsd=dsqrt(rmsd/nat1)
c     write(*,'(''RMSD '',f14.8)') rmsd             



      end

c
c     function to return the index of last nonblank character
c
      integer function lstchr(string)
      character*(*) string
      long=len(string)
      lstchr=0
      do i=long,1,-1
        if(string(i:i).ne.' ') then
          lstchr=i
          return
        endif
      enddo
      return
      end
