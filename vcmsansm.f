c programme pour calculer le vitesse et la position du centre de masse
c d'une galaxie a partir d'un fichier ci.txt ou cf.txt
c
       double precision  vx1(1000),vy1(1000),vz1(1000),x1(1000),
     +y1(1000),z1(1000),vxcm1,vycm1,vzcm1,xcm1,ycm1,zcm1,m,mcm,v1,
     +bidon,mcin1,vmax
       double precision ecm1,epot1,r1,G,ecin1,l1,lx1,ly1,lz1,normx1,
     +normy1,normz1
       double precision norme1,erot1,miner1,rcm1,vcm1,nomoy1

       double precision  vx2(2000),vy2(2000),vz2(2000),x2(2000),
     +y2(2000),z2(2000),vxcm2,vycm2,vzcm2,xcm2,ycm2,zcm2,v2,
     +mcin2
       double precision ecm2,epot2,r2,ecin2,l2,lx2,ly2,lz2,normx2,
     +normy2,normz2
       double precision norme2,erot2,miner2,rcm2,vcm2,nomoy2,ecin

       double precision r,epot,etot,ecint,ek

       integer i,netoile,j,nmax
       character*2 nom
       character*40 nom1,nom2
       G=6.6725985e-11
       print*,'Numero du fichier (vecteur##.txt)'
       read*,nom
       nom1='gal1vecteur'//nom//'.txt'
       nom2='gal2vecteur'//nom//'.txt'

       xcm1=0.
       ycm1=0.
       zcm1=0.
       vxcm1=0.
       vycm1=0.
       vzcm1=0.
       xcm2=0.
       ycm2=0.
       zcm2=0.
       vxcm2=0.
       vycm2=0.
       vzcm2=0.


c       print*,'Entrez la masse d une etoile:'
c       read*,m
       open(unit=2,file='ci.txt',status='old')
          read(2,*) bidon
          read(2,*) bidon,bidon,bidon,bidon,bidon,bidon,m
       close(unit=2)
       print*,'Masse d une etoile=',m
       mcm=0
       open(unit=1,file=nom1,status='old')
          read(1,*) netoile
       ecin1=0.
       miner1=0.
       lx1=0.
       ly1=0.
       lz1=0.
       nomoy1=0.

       nmax=0
          do i=1,netoile
             read(1,*) x1(i),y1(i),z1(i),vx1(i),vy1(i),vz1(i)

             xcm1=xcm1+x1(i)*m
             ycm1=ycm1+y1(i)*m
             zcm1=zcm1+z1(i)*m
             vx1(i)=vx1(i)/1.e15
             vy1(i)=vy1(i)/1.e15
             vz1(i)=vz1(i)/1.e15


             vxcm1=vxcm1+vx1(i)*m
             vycm1=vycm1+vy1(i)*m
             vzcm1=vzcm1+vz1(i)*m

         ecin1=ecin1+0.5*m*(vx1(i)**2.+vy1(i)**2.+vz1(i)**2.)
c calcul du moment cinetique 
         normx1=(y1(i)-ycm1)*(vz1(i)-vzcm1)-(z1(i)-zcm1)*(vy1(i)-vycm1)
         normy1=(z1(i)-zcm1)*(vx1(i)-vxcm1)-(x1(i)-xcm1)*(vz1(i)-vzcm1)
         normz1=(x1(i)-xcm1)*(vy1(i)-vycm1)-(y1(i)-ycm1)*(vx1(i)-vxcm1)

         norme1=sqrt(normx1**2.+normy1**2.+normz1**2.)
         nomoy1=nomoy1+norme1
c         normx=normx/norme
c         normy=normy/norme
c         normz=normz/norme
         rcm1=sqrt((x1(i)-xcm1)**2.+(y1(i)-ycm1)**2.+(z1(i)-zcm1)
     +**2.)
         vcm1=sqrt(vx1(i)**2.+vy1(i)**2.+vz1(i)**2.)
         lx1=lx1+m*normx1
         ly1=ly1+m*normy1
         lz1=lz1+m*normz1
         miner1=miner1+m*rcm1**2.
             mcm=mcm+m
          enddo
          xcm1=xcm1/mcm
          ycm1=ycm1/mcm
          zcm1=zcm1/mcm
          vxcm1=vxcm1/mcm
          vycm1=vycm1/mcm
          vzcm1=vzcm1/mcm
            v1=sqrt(vxcm1**2.+vycm1**2.+vzcm1**2.)
            ecm1=0.5*dble(netoile)*m*v1**2.
                l1=sqrt(lx1**2.+ly1**2.+lz1**2.)
                erot1=l1**2./(2.*miner1)
                nomoy1=m*nomoy1/dble(netoile)
          print*,'Position CM gal1=',xcm1,ycm1,zcm1
          print*,'Vitesse CM gal1=',vxcm1,vycm1,vzcm1
          print*,'grandeur vitesse CM gal1=',v1
          print*,'energie cinetique du CM gal1=',ecm1
       close(unit=1)

          ek=0.
          do i=1,netoile
             vmax=sqrt((vx1(i)-vxcm1)**2.+(vy1(i)-vycm1)**2.+
     +(vz1(i)-vzcm1)**2.)
             if (vmax.gt.2.e5) then
                ek=ek+0.5*m*vmax**2.
                print*,'On depasse c!',ek, vmax
                nmax=nmax+1
                print*,nmax
             endif
         enddo
         ecin1=ecin1-ek

c calcul de l'energie potentielle
c
       epot1=0.
       do i=1,netoile
         do j=1,netoile
            if (i.ne.j) then
               r1=sqrt((x1(i)-x1(j))**2.+(y1(i)-y1(j))**2.+(z1(i)-
     +z1(j))**2.)
               epot1=epot1-G*m*m/r1
            endif
          enddo
        enddo
        epot1=epot1/2.
        print*,'Energie potentielle propre galaxie 1',epot1



       mcm=0
       open(unit=1,file=nom2,status='old')
          read(1,*) netoile
       ecin2=0.
       miner2=0.
       lx2=0.
       ly2=0.
       lz2=0.
       nomoy2=0.
          do i=1,netoile
             read(1,*) x2(i),y2(i),z2(i),vx2(i),vy2(i),vz2(i)
             xcm2=xcm2+x2(i)*m
             ycm2=ycm2+y2(i)*m
             zcm2=zcm2+z2(i)*m
             vx2(i)=vx2(i)/1.e15
             vy2(i)=vy2(i)/1.e15
             vz2(i)=vz2(i)/1.e15
             vxcm2=vxcm2+vx2(i)*m
             vycm2=vycm2+vy2(i)*m
             vzcm2=vzcm2+vz2(i)*m
         ecin2=ecin2+0.5*m*(vx2(i)**2.+vy2(i)**2.+vz2(i)**2.)
c calcul du moment cinetique 
         normx2=(y2(i)-ycm2)*(vz2(i)-vzcm2)-(z2(i)-zcm2)*(vy2(i)-vycm2)
         normy2=(z2(i)-zcm2)*(vx2(i)-vxcm2)-(x2(i)-xcm2)*(vz2(i)-vzcm2)
         normz2=(x2(i)-xcm2)*(vy2(i)-vycm2)-(y2(i)-ycm2)*(vx2(i)-vxcm2)

         norme2=sqrt(normx2**2.+normy2**2.+normz2**2.)
         nomoy2=nomoy2+norme2
c         normx=normx/norme
c         normy=normy/norme
c         normz=normz/norme
         rcm2=sqrt((x2(i)-xcm2)**2.+(y2(i)-ycm2)**2.+(z2(i)-zcm2)**2.)
         vcm2=sqrt(vx2(i)**2.+vy2(i)**2.+vz2(i)**2.)
         lx2=lx2+m*normx2
         ly2=ly2+m*normy2
         lz2=lz2+m*normz2
         miner2=miner2+m*rcm2**2.
             mcm=mcm+m
          enddo
          xcm2=xcm2/mcm
          ycm2=ycm2/mcm
          zcm2=zcm2/mcm
          vxcm2=vxcm2/mcm
          vycm2=vycm2/mcm
          vzcm2=vzcm2/mcm
            v2=sqrt(vxcm2**2.+vycm2**2.+vzcm2**2.)
            ecm2=0.5*dble(netoile)*m*v2**2.
                l2=sqrt(lx2**2.+ly2**2.+lz2**2.)
                erot2=l2**2./(2.*miner2)
                nomoy2=m*nomoy2/dble(netoile)
          print*,'Position CM gal2=',xcm2,ycm2,zcm2
          print*,'Vitesse CM gal2=',vxcm2,vycm2,vzcm2
          print*,'grandeur vitesse CM gal2=',v2
          print*,'energie cinetique du CM gal2=',ecm2
       close(unit=1)
          ek=0.
          do i=1,netoile
             vmax=sqrt((vx2(i)-vxcm2)**2.+(vy2(i)-vycm2)**2.+
     + (vz2(i)-vzcm2)**2.)
             if (vmax.gt.5.e5) then
                ek=ek+0.5*m*vmax**2.
                print*,'On depasse c!',ek, vmax
                nmax=nmax+1
                print*,nmax
             endif
         enddo
         ecin2=ecin2-ek



c calcul de l'energie potentielle
c
       epot2=0.
       do i=1,netoile
         do j=1,netoile
            if (i.ne.j) then
               r2=sqrt((x2(i)-x2(j))**2.+(y2(i)-y2(j))**2.+(z2(i)-
     +z2(j))**2.)
               epot2=epot2-G*m*m/r2
            endif
          enddo
        enddo
        epot2=epot2/2
        print*,'Energie potentielle propre galaxie 2',epot2

c calcul de l energie potentielle totale
       do i=1,netoile
         do j=1,netoile
            if (i.ne.j) then
               r=sqrt((x1(i)-x1(j))**2.+(y1(i)-y1(j))**2.+(z1(i)-
     +z1(j))**2.)
               epot=epot-G*m*m/r
            endif
          enddo
           do j=1,netoile
               r=sqrt((x1(i)-x2(j))**2.+(y1(i)-y2(j))**2.+(z1(i)-
     +z2(j))**2.)
               epot=epot-G*m*m/r
          enddo        
        enddo

       do i=1,netoile
         do j=1,netoile
            if (i.ne.j) then
               r=sqrt((x2(i)-x2(j))**2.+(y2(i)-y2(j))**2.+(z2(i)-
     +z2(j))**2.)
               epot=epot-G*m*m/r
            endif
          enddo
           do j=1,netoile
               r=sqrt((x2(i)-x1(j))**2.+(y2(i)-y1(j))**2.+(z2(i)-
     +z1(j))**2.)
               epot=epot-G*m*m/r
          enddo        
        enddo
        epot=epot/2.


          ecin=ecin1+ecin2
          etot=epot+ecin1+ecin2
          ecint=ecin1-ecm1+ecin2-ecm2

          print*,'Bilan total'
          print*,'Energie potentielle totale (les deux galaxies) =',
     +epot
         print*,'Energie cinetique totale (les deux galaxies) =',ecin
     +,ecin1,ecin2
          print*,'Energie totale (les deux galaxies) =',etot
          print*,'Energie cinetique interne (rot et disp.)=',
     +ecint
c          print*,'Energie de rotation interne',erot

       print*,vxcm2,vycm2,vzcm2
       stop
       end
