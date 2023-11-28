c
c  This program uses the *.rmsd files in ./RMSD to identify the
c  matching conformer pairs between the QM and FF sets,
c  then picks a refeceence pair from among the matches,
c  then computes the relative energies on the QM and FF
c  surfaces.  
c
c  Then it identifies and reports problems with each compound.
c  These problems include 
c  --orphans (especially bad if all conformers are orphans)
c  --orphans (when the lowest energy conformer on either surface
c            is an orphan)
c  --ddE too big for conformers that are low in energy on either
c            surface
c
c  This program was derived from getmatch.f in NEWBENCH/DirtyDozen/RMSD/
c
c
c  Swope, Sept 2022
c
c  Program is driven from a file "compound.list", with the list 
c  of compounds:
c  Program reads files produced by a prior run of computermsd.f
c  This deposits into the directory ./RMSD  one file for each 
c  compound with the RMSD matrix, the names of the conformer files
c  and their energies.
c
c
c  1) identify matches and orphans for each compound
c  based on the rmsd matrix 
c  2) identify reference match for establishing relative 
c  energy scales
c  3) compute the relative energies on each surface
c  4) collect statistics on various kinds of failure over the 
c  compound set 

c  loop over each compound, get the array of rmsd 
c  values computed
c  find the closest of each quantum conformer
c  find the parallel matches for each quantum conformer
c
      implicit real*8 (a-h,o-z)
      character*1024 line,compoundfile,rmsdfile
      character*1024 confname(10,2),cmpdname,rmsddir
      dimension energies(10,2)
      dimension rmsd2d(10,10)
      dimension qmconf(0:9,10),ffconf(0:9,10)
      dimension iqde(0:10),qde(0:10)  !for sorting the relative energies
      dimension ifde(0:10),fde(0:10)  !for sorting the relative energies
 


      character*1024 dir1,dir2
      dimension rmsd(200,3)


      dimension closest(0:20,2),closestref(0:20,2)
      dimension matches(0:9,2)
      dimension iscore(100)


      threshstruc=0.4d0 !threshold for the RMSD where we consider pair matched or not
      threshde=4.0d0   !threshold for |dE(QM or FF)| for what counts are low energy
      threshdde=2.0d0  !threshold for |ddE| for what we consider a good/bad ddE


c     200 is max number of QM-FF pairs; should be .ge. 100 here
c     since we have at most 10 qm conformers and 10 ff conformers
      do i=1,200
        do j=1,3
          rmsd(i,j)=0.d0
        enddo
      enddo

      do i=0,20
        closest(i,1)=-1.d0
        closest(i,2)=1000.d0
      enddo
 
      do i=0,9   !index qm conformers
        do j=0,9 !index ff conformers
          rmsd2d(i,j)=9.d0
        enddo
      enddo

 
 

c-----------------new code------------------------

c     program is driven from a file of unique compound names

      call getarg(1,rmsddir     )  !get the directory name with the rmsd files
      write(*,'(''Using RMSD files in directory '',a)') 
     x  rmsddir(1:lstchr(rmsddir))

      write(*,'(''Threshold for RMSD to be considered a match '',
     x  f10.2)') threshstruc
      write(*,'(''Threshold for relative energy for a '',
     x  ''conformer to be consideded low energy '',f10.2)') threshde
      write(*,'(''Threshold for ddE to be considered to be a '',
     x  ''successful dE value '',f10.2)') threshdde


      rmsddir='./RMSD-OpenFF/'
c     rmsddir='./RMSD-OPLS4-def/'
c     rmsddir='./RMSD-OPLS4-cst/'

      compoundfile='./compound.list'
      is=lstchr(compoundfile)
      open(10,file=compoundfile(1:is))

      ncomp=0
10    continue
      read(10,'(a1024)',end=200) line
      ncomp=ncomp+1
      is=lstchr(line)
      cmpdname=line(1:is)
      write(*,'(/,''Compound '',i4,'': '',a)') ncomp,line(1:is)
c     rmsdfile='./RMSD-OpenFF/'//line(1:is)//'.rmsd'
c     rmsdfile='./RMSD-OPLS4-def/'//line(1:is)//'.rmsd'
      rmsdfile=rmsddir(1:lstchr(rmsddir))//line(1:is)//'.rmsd'
      is=lstchr(rmsdfile)
c     write(*,'(''Opening rmsd file '',a)') rmsdfile(1:is)
      open(12,file=rmsdfile)
      read(12,*) nqmconf
      do i=1,nqmconf
        read(12,'(a1024)') confname(i,1)
        read(12,*) energies(i,1)
      enddo
      read(12,*) nffconf
      do i=1,nffconf
        read(12,'(a1024)') confname(i,2)
        read(12,*) energies(i,2)
      enddo
      do i=1,10
         do j=1,10
           rmsd2d(i,j)=0.d0
         enddo
      enddo

      do i=1,nqmconf
        read(12,*)  (rmsd2d(i,j),j=1,nffconf)
      enddo
      close(12)

      do i=0,9              ! xxconf(i,1) = 0/1 if conformer i is present
        do j=1,10           !          2  = index (1-10) into confname
          qmconf(i,j)=0.d0  !          3  = index (0-9) of closest conf
          ffconf(i,j)=0.d0  !          4  = rmsd to closest conf
        enddo               !          5  = index of match (-1 or 0-9)
        qmconf(i,1)=-1.d0   !          6  = final energy
        ffconf(i,1)=-1.d0   !          7  = relative energy
        qmconf(i,2)=-1.d0   !                             
        ffconf(i,2)=-1.d0   !                                
        qmconf(i,5)=-1.d0   ! set to -1 (not matched)        
        ffconf(i,5)=-1.d0   ! set to -1 (not matched)        
      enddo

      do i=1,100
        iscore(i)=0
      enddo

      do i=1,nqmconf
        is=index(confname(i,1),'.sdf')-2
        read(confname(i,1)(is:is+1),'(i2)') icindex
        write(*,'(''QM conf '',i2,'' icindex '',i2,
     x    '' File '',a,'' Energy '',f15.4)')
     x    i,icindex,confname(i,1)(1:lstchr(confname(i,1))),
     x    energies(i,1)
        qmconf(icindex,1)=1 !this conformer is present
        qmconf(icindex,2)=i !pointer into the confname  
      enddo

      do i=1,nffconf
        is=index(confname(i,2),'.sdf')-2
        read(confname(i,2)(is:is+1),'(i2)') icindex
        write(*,'(''FF conf '',i2,'' icindex '',i2,
     x    '' File '',a,'' Energy '',f15.4)')
     x    i,icindex,confname(i,2)(1:lstchr(confname(i,2))),
     x    energies(i,2)
        ffconf(icindex,1)=1 !this conformer is present
        ffconf(icindex,2)=i !pointer into the confname  
      enddo

      write(*,'(/,''RMSD Matrix '')')
      do i=1,nqmconf
        write(*,'(10f8.4)')(rmsd2d(i,j),j=1,nffconf)
      enddo

c     if(nqmconf.lt.nffconf) then
c       write(*,'(''Warning:  '',i2,'' missing QM conformers'')')
c    x  nffconf-nqmconf
c     endif
c     if(nffconf.lt.nqmconf) then
c       write(*,'(''Warning:  '',i2,'' missing FF conformers'')')
c    x  nqmconf-nffconf
c     endif

c     check if each parent qm conformer has a child present
c     xxconf(i,1) should be -1 (not present) or +1 (present)
      nmissing=0
      do i=0,9
        if(int(qmconf(i,1)).eq.1.and.int(ffconf(i,1)).ne.1) then
          write(*,'(''Missing child: QM conformer '',
     x      i2,'' has no child.'')') i
          nmissing=nmissing+1
        endif
        if(qmconf(i,1).ne.1.and.ffconf(i,1).eq.1) then
          write(*,'(''Missing parent: FF conformer '',
     x      i2,'' has no parent.'')') i
        endif
      enddo

      iscore(1)=nmissing !save number of missing FF conformers



c     find the index and rmsd to the closest ffmatch
      nqmmatch=0            
      do i=1,nqmconf
        ismall=-1
        rsmall=100.d0
        do j=1,nffconf
          if(rmsd2d(i,j).lt.rsmall) then
            ismall=j
            rsmall=rmsd2d(i,j)
          endif
        enddo
        is=index(confname(i,1),'.sdf')-2
        read(confname(i,1)(is:is+1),'(i2)') iqindex
        if(ismall.ne.-1) then
          is=index(confname(ismall,2),'.sdf')-2
          read(confname(ismall,2)(is:is+1),'(i2)') ifindex
          qmconf(iqindex,3)=ifindex
        endif
        qmconf(iqindex,4)=rsmall
        qmconf(iqindex,5)=-1  ! set to -1 => not matched
        if(int(qmconf(iqindex,1)).eq.1.and.
     x     qmconf(iqindex,4).le.threshstruc) then
          qmconf(iqindex,5)=qmconf(iqindex,3)
          nqmmatch=nqmmatch+1
        endif
        qmconf(iqindex,6)=energies(i,1)
      enddo

c     find the index and rmsd to the closest qmmatch
      nffmatch=0
      do i=1,nffconf
        ismall=1
        rsmall=rmsd2d(1,i)
        do j=2,nqmconf
          if(rmsd2d(j,i).lt.rsmall) then
            ismall=j
            rsmall=rmsd2d(j,i)
          endif
        enddo
        is=index(confname(i,2),'.sdf')-2
        read(confname(i,2)(is:is+1),'(i2)') ifindex
        is=index(confname(ismall,1),'.sdf')-2
        read(confname(ismall,1)(is:is+1),'(i2)') iqindex
        ffconf(ifindex,3)=iqindex
        ffconf(ifindex,4)=rsmall
        ffconf(ifindex,5)=-1  ! set to -1 => not matched
        if(ffconf(ifindex,1).eq.1.and.
     x     ffconf(ifindex,4).le.threshstruc) then
          ffconf(ifindex,5)=ffconf(ifindex,3)
          nffmatch=nffmatch+1
        endif
        ffconf(ifindex,6)=energies(i,2)
      enddo

c     write(*,'(/,''Number of matches to qm conformers: '',i2,
c    x  '' of '',i2)') nqmmatch,nqmconf
c     write(*,'(''Number of matches to ff conformers: '',i2,
c    x  '' of '',i2)') nffmatch,nffconf

c     write(*,'(/,''qmconf data structure'')')
c     do i=0,9
c       write(*,'(i2,3x,i3,i3,i3,f10.3,i3,f12.3)')
c    x  i,int(qmconf(i,1)),int(qmconf(i,2)),int(qmconf(i,3)),
c    x  qmconf(i,4),int(qmconf(i,5)),qmconf(i,6)
c     enddo

c     write(*,'(/,''ffconf data structure'')')
c     do i=0,9
c       write(*,'(i2,3x,i3,i3,i3,f10.3,i3,f12.3)')
c    x  i,int(ffconf(i,1)),int(ffconf(i,2)),int(ffconf(i,3)),
c    x  ffconf(i,4),int(ffconf(i,5)),ffconf(i,6)
c     enddo


c     report the best matches                            
c     write(*,'(/,''Best matches for each quantum conformation'')')
c     do i=0,9
c       if(qmconf(i,1).eq.1) then
c         write(*,'(''QM Conformer '',i2,
c    x      '' best match '',i2,5x,'' rmsd '',f12.6)')
c    x      i,int(qmconf(i,3)),qmconf(i,4)
c       endif
c     enddo 
c     write(*,'(/,''Best matches for each forcefield conformer'')')
c     do i=0,9
c       if(ffconf(i,1).eq.1) then
c         write(*,'(''FF Conformer '',i2,
c    x      '' best match '',i2,5x,'' rmsd '',f12.6)')
c    x      i,int(ffconf(i,3)),ffconf(i,4)
c       endif
c     enddo

c     write out matches and orphans
      write(*,'(/,''Matches/orphans of the QM conformers'')')
      nqmorph=0
      nqmmatch=0
      do i=0,9
        if(int(qmconf(i,1)).eq.1) then
          if(int(qmconf(i,5)).eq.-1) then
            write(*,'(''QM conf '',i1,'' is an orphan'')') i 
            nqmorph=nqmorph+1
          else
            write(*,'(''QM conf '',i1,'' matches FF conf '',i1)')
     x      i,int(qmconf(i,5))
            nqmmatch=nqmmatch+1
          endif
        endif
      enddo

      write(*,'(/,''Matches/orphans of the FF conformers'')')
      nfforph=0
      nffmatch=0
      do i=0,9
        if(int(ffconf(i,1)).eq.1) then
          if(int(ffconf(i,5)).eq.-1) then
            write(*,'(''FF conf '',i1,'' is an orphan'')') i
            nfforph=nfforph+1
          else
            write(*,'(''FF conf '',i1,'' matches QM conf '',i1)')
     x      i,int(ffconf(i,5))
            nffmatch=nffmatch+1
          endif
        endif
      enddo

      write(*,'(/,''QM matches/orphans/conformers  '',3i3,
     x       '' ; FF matches/orphans/conformers  '',3i3)')
     x    nqmmatch,nqmorph,nqmconf,
     x    nffmatch,nfforph,nffconf

      iscore(2)=nqmconf
      iscore(3)=nqmmatch
      iscore(4)=nqmorph
      iscore(5)=10d0*nqmorph/nqmconf

c     report on strange case when nqmmatch .ne. nffmatch
      if(nqmmatch.ne.nffmatch) then
        write(*,'(/,''Warning: number QM matches '',
     x    ''ne number FF matches'')')
      endif

c     report on cases where matches do not correspond to parent-child pairs
      nmis=0
      do i=0,9
        if(int(qmconf(i,1)).ne.0.and.
     x     int(qmconf(i,5)).ne.i.and.
     x     int(qmconf(i,5)).ne.-1) then
          write(*,'(''Match for QM conformer '',i2,
     x      '' is to FF conformer '',i2,'' (not a child)'')')
     x     i,int(qmconf(i,5))
          nmis=nmis+1
        endif
      enddo

      iscore(6)=nmis  !count of the number of cases where matches are not p-c

      do i=0,9
        if(int(ffconf(i,1)).ne.0.and.
     x     int(ffconf(i,5)).ne.i.and.
     x     int(ffconf(i,5)).ne.-1) then
          write(*,'(''Match for FF conformer '',i2,
     x      '' is to QM conformer '',i2,'' (not a parent)'')')
     x     i,int(ffconf(i,5))
        endif
      enddo

c     iscore(12)   !cases of >1 QM conformer matched with one FF conformer QS 
      iscore(12)=0
c     count cases where a FF structure is matched with more than one QM structure
      do i=0,9  !loop over the FF conformers
        if(int(ffconf(i,1)).ne.0.and.int(ffconf(i,5)).ne.-1) then
          nm=0
          do j=0,9
            if(int(qmconf(j,1)).ne.0.and.int(qmconf(j,5)).ne.-1) then
              if(int(qmconf(j,5)).eq.i) nm=nm+1
            endif
          enddo
          if(nm.gt.1) iscore(12)=iscore(12)+1
        endif
      enddo

c     iscore(13)   !cases of >1 FF conformer matched with one QM conformer FS
      iscore(13)=0
c     count cases where a QM structure is matched with more than one FF structure
      do i=0,9  !loop over the QM conformers
        if(int(qmconf(i,1)).ne.0.and.int(qmconf(i,5)).ne.-1) then
          nm=0
          do j=0,9
            if(int(ffconf(j,1)).ne.0.and.int(ffconf(j,5)).ne.-1) then
              if(int(ffconf(j,5)).eq.i) nm=nm+1
            endif
          enddo
          if(nm.gt.1) iscore(13)=iscore(13)+1
        endif
      enddo






c     find conformer with lowest energy for all qm and ff conformers
c     find the lowest energy conformers that have a qm-ff match
c     flag cases where the lowest energy qm/ff conformer is an orphans
c     the reference case will be the one with the lowest energy 
c     that also has a match.

      iqmlow=-1  !index of the lowest energy quantum conformer
      ifflow=-1  !index of the lowest energy force field conformer
      do i=0,9
        if(qmconf(i,1).eq.1) then
          if(iqmlow.eq.-1) iqmlow=i
          eqmlow=qmconf(iqmlow,6)
          if(qmconf(i,6).lt.eqmlow) iqmlow=i
        endif
        if(ffconf(i,1).eq.1) then
          if(ifflow.eq.-1) ifflow=i
          efflow=ffconf(ifflow,6)
          if(ffconf(i,6).lt.efflow) ifflow=i
        endif
      enddo

      iqmmatchlow=-1    !index of lowest energy quantum conformer among matched
      eqmmatchlow=1.d8  
      iffmatchlow=-1    !index of lowest energy forcefi conformer among matched
      effmatchlow=1.d8
      do i=0,9
        if(qmconf(i,5).ne.-1) then
          if(qmconf(i,6).lt.eqmmatchlow) then
            iqmmatchlow=i
            eqmmatchlow=qmconf(i,6)
          endif
        endif
        if(ffconf(i,5).ne.-1) then
          if(ffconf(i,6).lt.effmatchlow) then
            iffmatchlow=i
            effmatchlow=ffconf(i,6)
          endif
        endif
      enddo

      write(*,'(/,''Identifying reference pair and energy'')')
      write(*,'(''Lowest energy QM conformer is         '',i2,
     x          '' ; E(QM) '',f18.6)') iqmlow,qmconf(iqmlow,6)
      if(ifflow.ne.-1) then
        write(*,'(''Lowest energy FF conformer is         '',i2,
     x            '' ; E(FF) '',f18.6)') ifflow,ffconf(ifflow,6)
      endif

      nmatch=nqmmatch

c     see if the lowest energy QM or FF conformer is an orphan
      iscore(7)=0
      iscore(8)=0
      if(int(qmconf(iqmlow,5)).eq.-1) iscore(7)=1
      if(ifflow.ne.-1) then
        if(int(ffconf(ifflow,5)).eq.-1) iscore(8)=1
      endif

      if(nmatch.gt.0) then
        write(*,'(''Lowest energy matched QM confomers is '',i2,
     x    '' ; E(QM) '',f18.6)') 
     x    iqmmatchlow,qmconf(iqmmatchlow,6)
        write(*,'(''Lowest energy matched FF conformer is '',i2,
     x    '' ; E(FF) '',f18.6)') 
     x    iffmatchlow,ffconf(iffmatchlow,6)
c       if(qmconf(iqmmatchlow,5).ne.iffmatchlow) iscore(11)=1
      endif


c     use the lowest energy quantum conformation that is matched
c     (iqmmatchlow) and its matched ff conformation 
c     qmconf(iqmmatchlow,5) as the reference pair for
c     the energy spectrum

      if(nmatch.gt.0) then
        qmrefengy=qmconf(iqmmatchlow,6)
        ffrefengy=ffconf(int(qmconf(iqmmatchlow,5)),6)
        write(*,'(''QM reference energy from conformer    '',i2,
     x    '' ; E(QM) '',f18.6)') 
     x    iqmmatchlow,qmrefengy
        write(*,'(''FF reference energy from conformer    '',i2,
     x    '' ; E(FF) '',f18.6)') 
     x    int(qmconf(iqmmatchlow,5)),ffrefengy
        do i=0,9
          iqde(i)=-1
          ifde(i)=-1
          qde(i)=-100.d0
          fde(i)=-100.d0
          if(qmconf(i,5).ne.-1) then
            iqde(i)=i
            ifde(i)=i
            qmconf(i,7)=qmconf(i,6)-qmrefengy              
            qde(i)=qmconf(i,7)
            m=qmconf(i,5)
            ffrel      =ffconf(m,6)-ffrefengy
            fde(i)=ffrel
            dde=ffrel-qmconf(i,7)
c           if(dabs(dde).le.2.d0) iscore(12)=iscore(12)+1
          endif
        enddo

c       sort the dE(QM) and dE(FF) energies
        do i=0,8
          do j=i+1,9
            if(qde(j).ge.qde(i)) then
              xswap=qde(i)
              qde(i)=qde(j)
              qde(j)=xswap
              iswap=iqde(i)
              iqde(i)=iqde(j)
              iqde(j)=iswap
            endif
            if(fde(j).ge.fde(i)) then
              xswap=fde(i)
              fde(i)=fde(j)
              fde(j)=xswap
              iswap=ifde(i)
              ifde(i)=ifde(j)
              ifde(j)=iswap
            endif
          enddo
        enddo

        write(*,'(/,''Relative energy ordered by dE(QM)'')')
        nmis=0   !count cases where dE(FF)<0
        nflag=0  !count cases where |dE(QM or FF)| .le. 4 but |ddE|>2
        nlowengy=0 !count cases where |dE(QM or FF)| .le. 4
        nlowpass=0 !count cases where low energy AND passing |ddE|<2
        do i=0,9
          if(iqde(i).ne.-1) then    !this qm conf has a match
            m=qmconf(iqde(i),5)     !get the FF conformer that matches
            ffrel=ffconf(m,6)-ffrefengy !rel FF energy of the match
            if(ffrel.lt.0.d0) nmis=nmis+1
            dde=ffrel-qmconf(iqde(i),7)  !ddE = dE(FF) - dE(FF)
c           following line give the results reported in the Powerpoint presentation
c           if(dabs(dde).ge.2.d0) nflag=nflag+1
c           following is preferred since it only counts errors for low energy conformers
            iflag=0
            if(dabs(ffrel).le.threshde  .or.
     x         dabs(qmconf(iqde(i),7)).le.threshde) then
               nlowengy=nlowengy+1
              if(dabs(dde).gt.threshdde) then
                iflag=1
                nflag=nflag+1
              else
                nlowpass=nlowpass+1
              endif
            endif
            if(iflag.eq.0) then
              write(*,'(i5.5,2x,a9,'' QM Conf '',i1,'' dE(QM) '',f10.3,
     x        '' matches FF Conf '',i1,'' dE(FF) '',f10.3,
     x        '' ddE '',f10.3)')
     x        ncomp,cmpdname(1:lstchr(cmpdname)),
     x        iqde(i),qmconf(iqde(i),7),m,ffrel,dde
            else
              write(*,'(i5.5,2x,a9,'' QM Conf '',i1,'' dE(QM) '',f10.3,
     x        '' matches FF Conf '',i1,'' dE(FF) '',f10.3,
     x        '' ddE '',f10.3,'' X '')')
     x        ncomp,cmpdname(1:lstchr(cmpdname)),
     x        iqde(i),qmconf(iqde(i),7),m,ffrel,dde
            endif
          endif
        enddo

        write(*,'(/,''Number low energy conformers '',i2,
     x    '' (|dE(QM or FF)|.le.'',f5.2,'')'')') 
     x    nlowengy,threshde
        write(*,'(''Number passing low engy conf '',i2,
     x    '' (|ddE|.le.'',f5.2,'')'')')
     x    nlowpass,threshdde
        flowpassscore=0.d0
        if(nlowpass.ge.1.and.nlowengy.ge.2)then 
          flowpassscore=min(99.d0,100.d0*(nlowpass-1)/(nlowengy-1))
          flowpassscore=max(1.d0,flowpassscore)    
        endif
        write(*,'(''Low energy pass score '',i3,i3,f6.2)')
     x    nlowpass,nlowengy,flowpassscore

        iscore(9)=nmis  !save number of FF structures with lower energy than the FF ref
        iscore(11)=nflag !save number of cases where ddE is large but dE is small
c       SC score
        iscore(14)=min(99,int(100.d0*(nqmmatch-nflag)/nqmconf))
c       NL score
        iscore(15)=flowpassscore

        write(*,*)
        nmisord=0
        nmisordlow=0
        do i=0,9
          if(iqde(i).ne.ifde(i)) then 
            write(*,'(''FF conformer '',i2,
     x        '' is misordered '')') ifde(i)
            nmisord=nmisord+1
c           count only low energy states for misordering
c           the qm conformer is iqde(i);     
            m=qmconf(iqde(i),5)     !get the FF conformer that matches
            ffrel=ffconf(m,6)-ffrefengy !rel FF energy of the match
            if(ffrel.le.threshde .or.
     x         dabs(qmconf(iqde(i),7)).le.threshde) then
              nmisordlow=nmisordlow+1
            endif
          endif
        enddo
        write(*,'(''Number of misordered FF conformers:          '',
     x    i2)') nmisord
        write(*,'(''Number of misordered low-engy FF conformers: '',
     x    i2)') nmisordlow
      endif

      iscore(10)=nmisord  !save the number of FF structures out of order rel to QM order
      iscore(10)=nmisordlow  !save the number of FF structures out of order rel to QM order

      write(*,'(''MD Num missing FF conformers '',i2)') iscore(1)
      write(*,'(''NC Num QM conformers         '',i2)') iscore(2)
      write(*,'(''NM Num QM matches            '',i2)') iscore(3)
      write(*,'(''NO Num QM orphans            '',i2)') iscore(4)
      write(*,'(''EO 10*(num orph)/(num conf)  '',i2)') iscore(5)
      write(*,'(''PC Num cases match .ne. p-c  '',i2)') iscore(6)
      write(*,'(''LQ 1 if lowest QM is orphan  '',i2)') iscore(7)
      write(*,'(''LF 1 of lowest FF is orphan  '',i2)') iscore(8)
      write(*,'(''BF Num FF strucs below ref   '',i2)') iscore(9)
      write(*,'(''MO FF strucs our of order    '',i2)') iscore(10)
      write(*,'(''BE Num small dE large ddE    '',i2)') iscore(11)
      write(*,'(''QS Num 2Q:1F matches         '',i2)') iscore(12)
      write(*,'(''FS Num 1Q:2F matches         '',i2)') iscore(13)
      write(*,'(''SC 0-99 score                '',i2)') iscore(14)
      write(*,'(''NL 0-99 score lowE passes    '',i2)') iscore(15)
      write(*,'(i5.5,2x,a9,'' MD '',i2,'' NC '',i2,'' NM '',i2,
     x                     '' NO '',i2,'' EO '',i2,'' PC '',i2,
     x                     '' LQ '',i2,'' LF '',i2,'' BF '',i2,
     x                     '' MO '',i2,'' BE '',i2,'' QS '',i2,
     x                     '' FS '',i2,'' SC '',i2,'' NL '',i2)')
     x  ncomp,cmpdname(1:lstchr(cmpdname)),
     x  (iscore(j),j=1,15)

      go to 10

200   continue
      write(*,'(''Number of compounds processed: '',i5)') ncomp



      stop 'debug'

c------------end-new code------------------------



11    continue




        iscore(1)=nqmconf
        iscore(2)=nqmconf-nqmorph
        iscore(3)=nqmorph
        iscore(4)=nffconf
        iscore(5)=nffconf-nfforph
        iscore(6)=nfforph
        if(nqmconf-nqmorph.eq.0.or.nffconf-nfforph.eq.0)
     x    iscore(7)=1 

        call getenergies(dir1,dir2,icomp,qmconf,ffconf)

c       find conformer with lowest energy for qm and ff conformers
c       find the lowest energy conformer that has a qm-ff match
c       flag cases where the lowest energy qm/ff conformer is an orphans
c       the reference case will be the one with the lowest energy 
c       that also has a match.
        iqmlow=-1  
        ifflow=-1
        do i=0,9
          if(qmconf(i,1).eq.1) then
            if(iqmlow.eq.-1) iqmlow=i
            eqmlow=qmconf(iqmlow,6)
            if(qmconf(i,6).lt.eqmlow) iqmlow=i
          endif
          if(ffconf(i,1).eq.1) then
            if(ifflow.eq.-1) ifflow=i
            efflow=ffconf(ifflow,6)
            if(ffconf(i,6).lt.efflow) ifflow=i
          endif
        enddo

c       check if lowest energy structures are orphans
        if(qmconf(iqmlow,4).eq.0) iscore(8)=1
        if(ffconf(ifflow,4).eq.0) iscore(9)=1

c       check if lowest energy qm and ff conformers are matched
        if(qmconf(iqmlow,4).eq.0) iscore(10)=1
        if(ffconf(ifflow,4).eq.0) iscore(10)=1
        if(qmconf(iqmlow,4).eq.1.and.
     x     qmconf(iqmlow,5).ne.ifflow) iscore(10)=1
        if(ffconf(ifflow,4).eq.1.and.
     x     ffconf(ifflow,5).ne.iqmlow) iscore(10)=1
   


        do i=0,9
          if(qmconf(i,1).eq.1.or.ffconf(i,1).eq.1) then
            write(*,'(''Conformer '',i1,'' E(QM) '',e18.12,
     x        '' E(FF) '',e18.12)') i,qmconf(i,6),ffconf(i,6)
          endif
        enddo

         

        iqmmatchlow=-1
        eqmmatchlow=1.d8
        iffmatchlow=-1
        effmatchlow=1.d8
        do i=0,9
          if(qmconf(i,4).eq.1) then
            if(qmconf(i,6).lt.eqmmatchlow) then
              iqmmatchlow=i
              eqmmatchlow=qmconf(i,6)
            endif
          endif
          if(ffconf(i,4).eq.1) then
            if(ffconf(i,6).lt.effmatchlow) then
              iffmatchlow=i
              effmatchlow=ffconf(i,6)
            endif
          endif
        enddo

        write(*,'(''Conformer '',i1,
     x    '' has lowest E(QM) '',e18.12)') iqmlow,qmconf(iqmlow,6)
        write(*,'(''Conformer '',i1,
     x    '' has lowest E(FF) '',e18.12)') ifflow,ffconf(ifflow,6)
        if(nmatch.gt.0) then
          write(*,'(''Of matched conformers '',i1,
     x      '' has lowest E(QM) '',e18.12)') 
     x      iqmmatchlow,qmconf(iqmmatchlow,6)
          write(*,'(''Of matched conformers '',i1,
     x      '' has lowest E(FF) '',e18.12)') 
     x      iffmatchlow,ffconf(iffmatchlow,6)
          if(qmconf(iqmmatchlow,5).ne.iffmatchlow) iscore(11)=1
        endif

c       print warning if lowest energy conformer is an orphan
        if(qmconf(iqmlow,1).eq.1.and.qmconf(iqmlow,4).eq.0) 
     x    then
          write(*,'(''Warning: Lowest energy QM conformation '',
     x    i1,'' is an orphan'')')  iqmlow
        endif
        if(ffconf(ifflow,1).eq.1.and.ffconf(ifflow,4).eq.0) 
     x    then
          write(*,'(''Warning: Lowest energy FF conformation '',
     x    i1,'' is an orphan'')')  ifflow
        endif

c       print warning if the lowest enery QM structure is not matched
c       with the lowest energy FF structure
        if(qmconf(iqmmatchlow,2).ne.iffmatchlow
     x    .or. ffconf(iffmatchlow,2).ne.iqmmatchlow) then
          write(*,'(''Warning: Lowest energy matched '',
     x      ''QM conf '',i1,
     x      '' not matched with lowest energy matched FF conf '',
     x      i1)') iqmmatchlow,iffmatchlow
        endif

c       use the lowest energy quantum conformation that is matched
c       (iqmmatchlow) and its matched ff conformation 
c       qmconf(iqmmatchlow,2) as the reference pair for
c       the energy spectrum
        qmrefengy=qmconf(iqmmatchlow,6)
        ffrefengy=ffconf(int(qmconf(iqmmatchlow,2)),6)
        write(*,'(''QM Ref energy from conf '',i1,5x,e18.12)')
     x    iqmmatchlow,qmrefengy
        write(*,'(''FF Ref energy from conf '',i1,5x,e18.12)')
     x    int(qmconf(iqmmatchlow,2)),ffrefengy
        do i=0,9
          if(qmconf(i,4).eq.1) then
            qmconf(i,7)=qmconf(i,6)-qmrefengy              
            m=qmconf(i,2)
            ffrel      =ffconf(m,6)-ffrefengy
            dde=ffrel-qmconf(i,7)
            if(dabs(dde).le.2.d0) iscore(12)=iscore(12)+1
            write(*,'(i5.5,'' QM Conf '',i1,'' dE(QM) '',f10.3,
     x      '' matches FF Conf '',i1,'' dE(FF) '',f10.3,
     x      '' ddE '',f10.3)')
     x      icomp,i,qmconf(i,7),m,ffrel,dde
          endif
        enddo

c-----------------------done checking to here--------------

c       write(*,'(''Report on compound '',i5.5)') icomp
        write(*,'(i5.5,'' Num QM conf/match/orphan '',3i3,
     x      '' Num FF conf/match/orphan '',3i3)') icomp,
     x      (iscore(k),k=1,6)
        write(*,'(i5.5,'' No matches '',i1)') icomp,
     x      iscore(7)
        write(*,'(i5.5,'' Lowest QM/FF is orphan '',2i2)') 
     x      icomp,iscore(8),iscore(9)
        write(*,'(i5.5,'' Lowest QM not matched to lowest FF '',
     x    i1)') icomp,iscore(10)
        write(*,'(i5.5,'' Lowest match QM not matched to '',
     x    ''lowest matched FF '',i1)') icomp,iscore(11)
        write(*,'(i5.5,'' Number of good ddE matches '',i2,
     x    '' of '',i2)')
     x    icomp,iscore(12),nmatch-1
        write(*,'(i5.5,'' Structure score '',f6.4,
     x                 '' ; Energy score '',f6.4,
     x                 '' ; Overall score '',f6.4)')
     x    icomp,
     x    dfloat(iscore(2))/dfloat(iscore(1)),
     x    dfloat(iscore(12))/dfloat(iscore(2)),
     x    dfloat(iscore(12))/dfloat(iscore(1)) 
        


c
c
c        if(iqmrefconf.ne.-1) then
c          call getenergies(dir1,dir2,icomp,energies)
c          write(*,'('' QM energy QM relative  '',
c     x              '' FF energy FF relative  '')')
c          do i=0,9
c            if(energies(i,1).ne.0.d0) then
c            write(*,'(''Conf '',i2,
cc    x                             '' E(QM) '',e20.8,
c     x                             '' E(QM)rel '',f8.4,
cc    x                             '' E(FF) '',e20.8,
c     x                             '' E(FF)rel '',f8.4,
c     x                             '' DDE(F-Q)'',f8.4)')
c     x        i,
cc    x          energies(i,1),
c     x          energies(i,1)-energies(iqmrefconf,1),
cc    x          energies(i,2),
c     x          energies(i,2)-energies(iffrefconf,2),
c     x          energies(i,2)-energies(iffrefconf,2) 
c     x       - (energies(i,1)-energies(iqmrefconf,1))
c     x          
c            endif
c          enddo
c          write(*,'(/,''Energy spectrum of matched pairs'')')
c          do i=0,9
c            if(energies(i,1).ne.0.d0
c     x        .and.matches(i,1).ne.-1.and.matches(i,2).ne.-1) then
c            write(*,'(''Conf '',i2,
cc    x                             '' E(QM) '',e20.8,
c     x                             '' E(QM)rel '',f8.4,
cc    x                             '' E(FF) '',e20.8,
c     x                             '' E(FF)rel '',f8.4,
c     x                             '' DDE(F-Q)'',f8.4)')
c     x          i,
cc    x            energies(i,1),
c     x            energies(i,1)-energies(iqmrefconf,1),
cc    x              energies(i,2),
c     x            energies(i,2)-energies(iffrefconf,2),
c     x           (energies(i,2)-energies(iffrefconf,2))
c     x       -   (energies(i,1)-energies(iqmrefconf,1))
c            endif
c          enddo
c        endif
c
cc     ---------------initialize storage structures----

        ncomp=ncomp+1
        icomp=icompnext
        write(*,'(/,''RMSD Matrix for compound '',i4,
     x    '' with ID '',i5.5)') 
     x    ncomp,icomp

        nline=0

        do i=1,200
          do j=1,3
            rmsd(i,j)=0.d0
          enddo
        enddo

        do i=0,9
          do j=0,9
            rmsd2d(i,j)=9.d0
          enddo
        enddo
 
        do i=0,9
          do j=1,10
            qmconf(i,j)=0.d0
            ffconf(i,j)=0.d0
          enddo
          qmconf(i,2)=-1.d0
          ffconf(i,2)=-1.d0
        enddo

        do i=0,20
          closest(i,1)=-1.d0
          closest(i,2)=1000.d0
        enddo

        do i=1,100
          iscore(i)=0
        enddo

c     endif

c     ---------------process the next conformer-----

      nline=nline+1

c     write(*,'(''Line '',i2,'' for compound '',i4,
c    x   '' with ID '',i4,5x,a)')
c    x   nline,ncomp,icomp,line(1:40)
      read(line(11:12),*) in1
      read(line(29:29),*) in2
      read(line(32:),*) r
c     write(*,'(a,2x,i2,5x,a,i2,5x,f12.6)') 
c    x  line(11:12),in1,line(28:29),in2,r
      rmsd(nline,1)=in1
      rmsd(nline,2)=in2
      rmsd(nline,3)=r
      rmsd2d(in1,in2)=r

      qmconf(in1,1)=1
      ffconf(in2,1)=1
      if(int(qmconf(in1,2)).eq.-1
     x  .or.r.le.qmconf(in1,3))then
        qmconf(in1,2)=in2
        qmconf(in1,3)=r
      endif
      if(int(ffconf(in2,2)).eq.-1
     x  .or.r.le.ffconf(in2,3))then
        ffconf(in2,2)=in1
        ffconf(in2,3)=r
      endif

      if(r.le.closest(in1,2)) then
        closest(in1,1)=in2
        closest(in1,2)=r
      endif

      
      go to 11

100   continue

c-------------------------process the last compound------------
      
c     write(*,'(''All RMSDs for compound '',i5)') icomp
      mxconf=0
      do i=1,nline
        iconf=int(rmsd(i,1))
        jconf=int(rmsd(i,2))
        mxconf=max(mxconf,iconf,jconf)
c       write(*,'(2i3,f12.6)') iconf,jconf,rmsd(i,3)
      enddo

c     write out the low precison 2D array of rmsds
      write(*,'(''2D array of rmsd values; nconf '',i2)') 
     x  mxconf+1 
      write(*,'(4x,10i5)') ((j),j=0,mxconf)
      do i=0,mxconf
        write(*,'(i4,10f5.2)') i,(rmsd2d(i,j),j=0,mxconf)
      enddo


c     Parallel (parent-child) compares
c     loop over the rmsd file and find "parallel matches"
c     write(*,'(''Parallel RMSDs compared  with reference'')')
      do i=1,nline
        iconf=int(rmsd(i,1))
        jconf=int(rmsd(i,2))
        qmconf(iconf,1)=1
        ffconf(jconf,1)=1
c       if(iconf.eq.jconf) then
c         write(*,'(2i4,5x,f14.8,'' Ref '',f14.8,
c    x     '' Diff '',f14.8)') 
c    x     iconf,jconf,rmsd(i,3),rmsdref(iconf),
c    x     rmsd(i,3)-rmsdref(iconf)
c       endif
      enddo



c     report the best matches                            
      write(*,'(''Best matches for each quantum conf'')')
      do i=0,9
          if(qmconf(i,1).eq.1) then
            write(*,'(''QM Conformer '',i2,
     x        '' best match '',i2,5x,'' rmsd '',f12.6)')
     x        i,int(qmconf(i,2)),qmconf(i,3)
          endif
      enddo

      write(*,'(''Best matches for each forcefield conf'')')
      do i=0,9
        if(ffconf(i,1).eq.1) then
          write(*,'(''FF Conformer '',i2,
     x      '' best match '',i2,5x,'' rmsd '',f12.6)')
     x      i,int(ffconf(i,2)),ffconf(i,3)
        endif
      enddo


      call matchup(qmconf,ffconf,nmatch,
     x    iqmrefconf,iffrefconf)
c     call matchup(rmsd2d,matches,mxconf,iqmrefconf,iffrefconf)

c     reference pair based on rmsd alone
c     write(*,'(''Reference pair:  QM conf '',i2,
c    x  '' with FF conf '',i2)')
c    x  iqmrefconf,iffrefconf

c     count qm orphans, ff orphans, qm matches, ff matches
      nqmorph=0
      nfforph=0
      nqmconf=0
      nffconf=0
      do i=0,9
        if(qmconf(i,1).ne.0) then
          nqmconf=nqmconf+1
          if(qmconf(i,4).eq.0) nqmorph=nqmorph+1
        endif
        if(ffconf(i,1).ne.0) then
          nffconf=nffconf+1
          if(ffconf(i,4).eq.0) nfforph=nfforph+1
        endif
      enddo
      write(*,'(''QM matches/orphans/conformers  '',3i3,
     x       '' ; FF matches/orphans/conformers  '',3i2)')
     x  nqmconf-nqmorph,nqmorph,nqmconf,
     x  nffconf-nfforph,nfforph,nffconf


      iscore(1)=nqmconf
      iscore(2)=nqmconf-nqmorph
      iscore(3)=nqmorph
      iscore(4)=nffconf
      iscore(5)=nffconf-nfforph
      iscore(6)=nfforph
      if(nqmconf-nqmorph.eq.0.or.nffconf-nfforph.eq.0)
     x  iscore(7)=1 


      call getenergies(dir1,dir2,icomp,qmconf,ffconf)

c     find conformer with lowest energy for qm and ff conformers
c     find the lowest energy conformer that has a qm-ff match
c     flag cases where the lowest energy qm/ff conformer is an orphans
c     the reference case will be the one with the lowest energy 
c     that also has a match.
      iqmlow=-1  
      ifflow=-1
      do i=0,9
        if(qmconf(i,1).eq.1) then
          if(iqmlow.eq.-1) iqmlow=i
          eqmlow=qmconf(iqmlow,6)
          if(qmconf(i,6).lt.eqmlow) iqmlow=i
        endif
        if(ffconf(i,1).eq.1) then
          if(ifflow.eq.-1) ifflow=i
          efflow=ffconf(ifflow,6)
          if(ffconf(i,6).lt.efflow) ifflow=i
        endif
      enddo

c     check if lowest energy structures are orphans
      if(qmconf(iqmlow,4).eq.0) iscore(8)=1
      if(ffconf(ifflow,4).eq.0) iscore(9)=1

c     check if lowest energy qm and ff conformers are matched
      if(qmconf(iqmlow,4).eq.0) iscore(10)=1
      if(ffconf(ifflow,4).eq.0) iscore(10)=1
      if(qmconf(iqmlow,4).eq.1.and.
     x   qmconf(iqmlow,5).ne.ifflow) iscore(10)=1
      if(ffconf(ifflow,4).eq.1.and.
     x   ffconf(ifflow,5).ne.iqmlow) iscore(10)=1


      do i=0,9
        if(qmconf(i,1).eq.1.or.ffconf(i,1).eq.1) then
          write(*,'(''Conformer '',i1,'' E(QM) '',e18.12,
     x      '' E(FF) '',e18.12)') i,qmconf(i,6),ffconf(i,6)
        endif
      enddo




      iqmmatchlow=-1
      eqmmatchlow=1.d8
      iffmatchlow=-1
      effmatchlow=1.d8
      do i=0,9
        if(qmconf(i,4).eq.1) then
          if(qmconf(i,6).lt.eqmmatchlow) then
            iqmmatchlow=i
            eqmmatchlow=qmconf(i,6)
          endif
        endif
        if(ffconf(i,4).eq.1) then
          if(ffconf(i,6).lt.effmatchlow) then
            iffmatchlow=i
            effmatchlow=ffconf(i,6)
          endif
        endif
      enddo

      write(*,'(''Conformer '',i1,
     x  '' has lowest E(QM) '',e18.12)') iqmlow,qmconf(iqmlow,6)
      write(*,'(''Conformer '',i1,
     x  '' has lowest E(FF) '',e18.12)') ifflow,ffconf(ifflow,6)
      if(nmatch.gt.0) then
        write(*,'(''Of matched conformers '',i1,
     x    '' has lowest E(QM) '',e18.12)') 
     x    iqmmatchlow,qmconf(iqmmatchlow,6)
        write(*,'(''Of matched conformers '',i1,
     x    '' has lowest E(FF) '',e18.12)') 
     x    iffmatchlow,ffconf(iffmatchlow,6)
        if(qmconf(iqmmatchlow,5).ne.iffmatchlow) iscore(11)=1
      endif

c     print warning if lowest energy conformer is an orphan
      if(qmconf(iqmlow,1).eq.1.and.qmconf(iqmlow,4).eq.0) 
     x  then
        write(*,'(''Warning:  lowest energy QM conformation '',
     x  i1,'' is an orphan'')')  iqmlow
      endif
      if(ffconf(ifflow,1).eq.1.and.ffconf(ifflow,4).eq.0) 
     x  then
        write(*,'(''Warning:  lowest energy FF conformation '',
     x  i1,'' is an orphan'')')  ifflow
      endif


c     print warning if the lowest enery QM structure is not matched
c     with the lowest energy FF structure
      if(qmconf(iqmmatchlow,2).ne.iffmatchlow
     x  .or. ffconf(iffmatchlow,2).ne.iqmmatchlow) then
        write(*,'(''Warning: Lowest energy matched '',
     x    ''QM conf '',i1,
     x    '' not matched to lowest energy matched FF conf'',
     x    i1)') iqmmatchlow,iffmatchlow
      endif


c     use the lowest energy quantum conformation that is matched
c     (iqmmatchlow) and its matched ff conformation 
c     qmconf(iqmmatchlow,2) as the reference pair for
c     the energy spectrum
      qmrefengy=qmconf(iqmmatchlow,6)
      ffrefengy=ffconf(int(qmconf(iqmmatchlow,2)),6)
      write(*,'(''QM Ref energy from conf '',i1,5x,e18.12)')
     x  iqmmatchlow,qmrefengy
      write(*,'(''FF Ref energy from conf '',i1,5x,e18.12)')
     x  int(qmconf(iqmmatchlow,2)),ffrefengy
      do i=0,9
        if(qmconf(i,4).eq.1) then
          qmconf(i,7)=qmconf(i,6)-qmrefengy              
          m=qmconf(i,2)
          ffrel      =ffconf(m,6)-ffrefengy
          dde=ffrel-qmconf(i,7)
          if(dabs(dde).le.2.d0) iscore(12)=iscore(12)+1
          write(*,'(i5.5'' QM Conf '',i1,'' dE(QM) '',f10.3,
     x    '' matches FF Conf '',i1,'' dE(FF) '',f10.3,
     x    '' ddE '',f10.3)')
     x    icomp,i,qmconf(i,7),m,ffrel,dde
        endif
      enddo

c     write(*,'(''Report on compound '',i5.5)') icomp
      write(*,'(i5.5,'' Num QM conf/match/orphan '',3i3,
     x    '' Num FF conf/match/orphan '',3i3)') icomp,
     x    (iscore(k),k=1,6)
      write(*,'(i5.5,'' No matches '',i1)') icomp,
     x    iscore(7)
      write(*,'(i5.5,'' Lowest QM/FF is orphan '',2i2)') 
     x    icomp,iscore(8),iscore(9)
      write(*,'(i5.5,'' Lowest QM not matched to lowest FF '',
     x  i1)') icomp,iscore(10)
      write(*,'(i5.5,'' Lowest match QM not matched to '',
     x  ''lowest matched FF '',i1)') icomp,iscore(11)
      write(*,'(i5.5,'' Number of good ddE matches '',i2,
     x  '' of '',i2)')
     x  icomp,iscore(12),nmatch-1
      write(*,'(i5.5,'' Structure score '',f6.4,
     x               '' ; Energy score '',f6.4,
     x               '' ; Overall score '',f6.4)')
     x  icomp,
     x  dfloat(iscore(2))/dfloat(iscore(1)),
     x  dfloat(iscore(12))/dfloat(iscore(2)),
     x  dfloat(iscore(12))/dfloat(iscore(1)) 
        


c
c
c      if(iqmrefconf.ne.-1) then
c        call getenergies(dir1,dir2,icomp,energies)
c        write(*,'('' QM energy QM relative  '',
c     x            '' FF energy FF relative  '')')
c        do i=0,9
c          if(energies(i,1).ne.0.d0) then
c          write(*,'(''Conf '',i2,
cc    x                             '' E(QM) '',e20.8,
c     x                             '' E(QM)rel '',f8.4,
cc    x                             '' E(FF) '',e20.8,
c     x                             '' E(FF)rel '',f8.4,
c     x                             '' DDE(F-Q)'',f8.4)')
c     x        i,
cc    x          energies(i,1),
c     x          energies(i,1)-energies(iqmrefconf,1),
cc    x          energies(i,2),
c     x          energies(i,2)-energies(iffrefconf,2),
c     x         (energies(i,2)-energies(iffrefconf,2))
c     x       - (energies(i,1)-energies(iqmrefconf,1))
c          endif
c        enddo
c        write(*,'(/,''Energy spectrum of matched pairs'')')
c        do i=0,9
c          if(energies(i,1).ne.0.d0
c     x      .and.matches(i,1).ne.-1.and.matches(i,2).ne.-1) then
c          write(*,'(''Conf '',i2,
cc    x                             '' E(QM) '',e20.8,
c     x                             '' E(QM)rel '',f8.4,
cc    x                             '' E(FF) '',e20.8,
c     x                             '' E(FF)rel '',f8.4,
c     x                             '' DDE(F-Q)'',f8.4)')
c     x        i,
cc    x          energies(i,1),
c     x          energies(i,1)-energies(iqmrefconf,1),
cc    x          energies(i,2),
c     x          energies(i,2)-energies(iffrefconf,2),
c     x         (energies(i,2)-energies(iffrefconf,2))
c     x       - (energies(i,1)-energies(iqmrefconf,1))
c          endif
c        enddo
c      endif
       
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
c
c     subroutine to find matches
c     the rmsd with the closest conformer must be within
c     threshold
c
c     subroutine matchup(rmsd2d,matches,mxconf,iqmref,iffref)
c     dimension matches(0:9,2),rmsd2d(0:9,0:9)

      subroutine matchup(qmconf,ffconf,
     x  nmatch,iqmref,iffref)
      implicit real*8 (a-h,o-z)
      dimension qmconf(0:9,10),ffconf(0:9,10)

      thresh=0.4d0   !rmsd must be .le. to this to be a match



      qmclose=100.d0
      ffclose=100.d0
      iqmclose=-1
      iffclose=-1
      do i=0,9
        if(qmconf(i,4).eq.1
     x    .and.qmconf(i,3).lt.qmclose) then
          qmclose=qmconf(i,3) 
          iqmclose=i 
        endif
        if(ffconf(i,4).eq.1
     x    .and.ffconf(i,3).lt.ffclose) then
          ffclose=ffconf(i,3) 
          iffclose=i 
        endif
      enddo

      if(nqmmatch.gt.0) 
     x  write(*,'(
     x    '' Closest QM match is between QM conf '',i1,
     x    '' and FF conf '',i1,'' with rmsd '',f10.4)')
     x    iqmclose,int(qmconf(iqmclose,5)),
     x    qmconf(iqmclose,3)
      if(nffmatch.gt.0) 
     x  write(*,'(
     x  '' Closest FF match is between FF conf '',i1,
     x  '' and QM conf '',i1,'' with rmsd '',f10.4)')
     x  iffclose,int(ffconf(iffclose,5)),
     x  ffconf(iffclose,3)


c     flag strange conditions

      if(nffmatch.ne.nqmmatch) then
        write(*,'(''Warning: Strange condition, num qm/ff '',
     x    ''matches not equal'')') 
      endif
      nmatch=nqmmatch
      if(nffmatch.eq.0.or.nqmmatch.eq.0) then
        write(*,'(''Warning: All conformers are orphans   '',
     x    ''nqmmatch '',i2,'' nffmatch '',i2)') 
     x    nqmmatch,nffmatch
      endif


c     iqmclose  !qm conf closest to any ff conf
c     qmconf(iqmclose,2) !the ff conformer it is closest to
c     iffclose  !ff conf closest to any qm conf
c     ffconf(iffclose,2) !the qm conformer it is closest to

      if(nmatch.gt.0) then
        if(iqmclose.ne.int(ffconf(iffclose,2))) then
          write(*,'(''Warning: closest qm and ff pairs '',
     x      ''mismatched'')') 
          write(*,'(''QM:  closest match is  '',2i4)')
     x    iqmclose,int(qmconf(iqmclose,2))
          write(*,'(''FF:  closest match is  '',2i4)')
     x    iffclose,int(ffconf(iffclose,2))
        endif
        if(iffclose.ne.int(qmconf(iqmclose,2))) then
          write(*,'(''Warning: closest qm and ff pairs '',
     x      ''mismatched'')') 
          write(*,'(''QM:  closest match is  '',2i4)')
     x    iqmclose,int(qmconf(iqmclose,2))
          write(*,'(''FF:  closest match is  '',2i4)')
     x    iffclose,int(ffconf(iffclose,2))
        endif
      endif

c
c     note that if all these conditions pass, the candidate
c     for the reference pair, based on rmsd, should be between
c     these unique pairings of the QM and FF conformers that
c     are closest to each other
      iqmref=-1
      iffref=-1
      if(nmatch.gt.0) then
        iqmref=iqmclose
        iffref=iffclose
      endif
c



      return
      end

