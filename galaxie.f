c
c  programme pour simuler le probleme du mouvement a N corps sous
c  l action de la force de gravite
c  programme par Martin Aube
c  utilisation de la methode de verlet
c
c =====
c Declaration des variables
c
       integer size
       parameter (size=2048)
       double precision xa(size),xn(size),ya(size),yn(size),za(size),
     + zin(size)
       double precision dt,dti,G,rroche1,rroche2,rroche
       double precision m(size),pi,ray
       double precision vxa(size),vya(size),vxn(size),vyn(size),
     + vzia(size),vzin(size)
       double precision axa(size),aya(size),axn(size),
     + ayn(size),aza(size),azin(size)
       integer step,nloop,i,j,mout,n,ne,nm,nn,ii
       character*30 nomvect
c
c =====
c Identification des variables
c
c m(n)=masse de la masse n
c nm=nombre de masses
c xa(n) ya(n) za(n) les position de la particile n a t_i
c xn(n) et yn(n) zn(n) sont les position de la particile n a t_n
c vxa(n) et vya(n) vzia(n) sont les vitesses de la particile n a t_i
c vxn(n) et vyn(n) vzn(n) sont les vitesses de la particile n a t_n
c axa(n) et aya(n) aza(n) sont les accelerations de la particile n a t_i
c axn(n) et ayn(n) azn(n) sont les accelerations de la particile n a t_n
c ray=rayon entre deux masses
c dt=pas de calcul
c dti=pas de calcul initial
c Tm,To,Ts periodes orbitalesza
c nloop=nombre de boucles de calcul
c mout=saut pour l impression de sorties
c
c
c Commentaires:
c Le programme attend en entree un fichier nomme ci.txt (conditions initiales)
c et retourne les conditions finales cf.txt dans le meme format
c le format est le nombre de masses sur la premiere ligne suivie sur les lignes
c suivantes (une ligne par masse) des x y z vx vy vz m
c
c Le programme fait 2E9 boucles de calculs temporels. La variable dti (a ajuster
c avant de compiler) est le pas de calcul. Pour connaitre le temps apres les
c 2E9 boucles il faut faire tf = 2E9 x dti - 0.99dti
c Dans sa version originale dti=1jour ce qui implique un tf>5E6 annees
c
c =====
c Valeurs initiales
c
c dti=86400*7 => une semaine
c
       dti=86400.*365.25*50.
       nloop=4000000
       mout=1000000
       n=0
       pi=3.14159
       G=6.6725985E-11
       dt=dti
       nn=0
c
c =====
c Lecture des positions et vitesses initiales
c
       open(unit=1,file='ci.txt',status='unknown')
       print*,'Reading initial conditions: file ci.txt'
         read(1,*) nm
         if (nm.gt.size) then
            print*,'Stars number exceed max value! Abort.'
            stop
         endif
         do i=1,nm
            read(1,*) xa(i),ya(i),za(i),vxa(i),vya(i),vzia(i),m(i)
         enddo
       close(unit=1)
c
c
c

c
c =====
c
c
c ouverture du fichier de sortie
       open(unit=2,file="cf.txt",status="unknown")
       do n=1,nloop
         do i=1,nm
           axn(i)=0.
           ayn(i)=0.
           azin(i)=0.
           do j=1,nm
             if (i.ne.j) then
               ray=((xa(j)-xa(i))**2.+(ya(j)-ya(i))**2.+
     +         (za(j)-za(i))**2.)**0.5

c              calcul du rayon de roche
c              (16^.333*rayon soleil/Masse du soleil ^0.333)*Masse ^0.333=0.14
               rroche1=m(i)**0.33333333333*0.14
               rroche2=m(j)**0.33333333333*0.14
               rroche=rroche1+rroche2
               if (ray.lt.rroche) then
                 print*,'limite de roche atteinte'
                 xa(i)=(m(i)*xa(i)+m(j)*xa(j))/(m(i)+m(j))
                 ya(i)=(m(i)*ya(i)+m(j)*ya(j))/(m(i)+m(j))
                 za(i)=(m(i)*za(i)+m(j)*za(j))/(m(i)+m(j))
                 vxa(i)=(m(i)*vxa(i)+m(j)*vxa(j))/(m(i)+m(j))
                 vya(i)=(m(i)*vya(i)+m(j)*vya(j))/(m(i)+m(j))
                 vzia(i)=(m(i)*vzia(i)+m(j)*vzia(j))/(m(i)+m(j))
                 m(i)=m(i)+m(j)
                 m(j)=0.
                 xa(j)=100.
                 ya(j)=0.
                 za(j)=0.
                 vxa(j)=0.
                 vya(j)=0.
                 vzia(j)=0.
               endif

               axn(i)=axn(i)+G*m(j)*(xa(j)-xa(i))/ray**3.
               ayn(i)=ayn(i)+G*m(j)*(ya(j)-ya(i))/ray**3.
               azin(i)=azin(i)+G*m(j)*(za(j)-za(i))/ray**3.
             endif
           enddo
         if (n.eq.1) then
c
c =====
c Ajuste la premiere valeur prededente a la valeur actuelle
c
           dt=dti/100.
           axa(i)=axn(i)
           aya(i)=ayn(i)
           aza(i)=azin(i)
         else
           dt=dti
         endif
         vxn(i)=vxa(i)+0.5*(axa(i)+axn(i))*dt
         vyn(i)=vya(i)+0.5*(aya(i)+ayn(i))*dt
         vzin(i)=vzia(i)+0.5*(aza(i)+azin(i))*dt
         xn(i)=xa(i)+vxa(i)*dt+0.5*axa(i)*dt**2.
         yn(i)=ya(i)+vya(i)*dt+0.5*aya(i)*dt**2.
         zin(i)=za(i)+vzia(i)*dt+0.5*aza(i)*dt**2.
         nn=nn+1
c
c =====
c Impression des resultats 
c
         if (nn.eq.mout) then
            nn=0
            print*," temps(jours)=",(dble(n-1)*dti+dti/100.)/3600./24.
     +      ,'/',dble(nloop)*dti/3600./24.,
     +      nint(100.*(dble(n-1)+1/100.)/dble(nloop)),'%'
         endif
          nomvect='toto'
         if (nloop/50.eq.n) nomvect='vecteur01.txt'
         if (2*nloop/50.eq.n) nomvect='vecteur02.txt'
         if (3*nloop/50.eq.n) nomvect='vecteur03.txt'
         if (4*nloop/50.eq.n) nomvect='vecteur04.txt'
         if (5*nloop/50.eq.n) nomvect='vecteur05.txt'
         if (6*nloop/50.eq.n) nomvect='vecteur06.txt'
         if (7*nloop/50.eq.n) nomvect='vecteur07.txt'
         if (8*nloop/50.eq.n) nomvect='vecteur08.txt'
         if (9*nloop/50.eq.n) nomvect='vecteur09.txt'
         if (10*nloop/50.eq.n) nomvect='vecteur10.txt'
         if (11*nloop/50.eq.n) nomvect='vecteur11.txt'
         if (12*nloop/50.eq.n) nomvect='vecteur12.txt'
         if (13*nloop/50.eq.n) nomvect='vecteur13.txt'
         if (14*nloop/50.eq.n) nomvect='vecteur14.txt'
         if (15*nloop/50.eq.n) nomvect='vecteur15.txt'
         if (16*nloop/50.eq.n) nomvect='vecteur16.txt'
         if (17*nloop/50.eq.n) nomvect='vecteur17.txt'
         if (18*nloop/50.eq.n) nomvect='vecteur18.txt'
         if (19*nloop/50.eq.n) nomvect='vecteur19.txt'
         if (20*nloop/50.eq.n) nomvect='vecteur20.txt'
         if (21*nloop/50.eq.n) nomvect='vecteur21.txt'
         if (22*nloop/50.eq.n) nomvect='vecteur22.txt'
         if (23*nloop/50.eq.n) nomvect='vecteur23.txt'
         if (24*nloop/50.eq.n) nomvect='vecteur24.txt'
         if (25*nloop/50.eq.n) nomvect='vecteur25.txt'
         if (26*nloop/50.eq.n) nomvect='vecteur26.txt'
         if (27*nloop/50.eq.n) nomvect='vecteur27.txt'
         if (28*nloop/50.eq.n) nomvect='vecteur28.txt'
         if (29*nloop/50.eq.n) nomvect='vecteur29.txt'
         if (30*nloop/50.eq.n) nomvect='vecteur30.txt'
         if (31*nloop/50.eq.n) nomvect='vecteur31.txt'
         if (32*nloop/50.eq.n) nomvect='vecteur32.txt'
         if (33*nloop/50.eq.n) nomvect='vecteur33.txt'
         if (34*nloop/50.eq.n) nomvect='vecteur34.txt'
         if (35*nloop/50.eq.n) nomvect='vecteur35.txt'
         if (36*nloop/50.eq.n) nomvect='vecteur36.txt'
         if (37*nloop/50.eq.n) nomvect='vecteur37.txt'
         if (38*nloop/50.eq.n) nomvect='vecteur38.txt'
         if (39*nloop/50.eq.n) nomvect='vecteur39.txt'
         if (40*nloop/50.eq.n) nomvect='vecteur40.txt'
         if (41*nloop/50.eq.n) nomvect='vecteur41.txt'
         if (42*nloop/50.eq.n) nomvect='vecteur42.txt'
         if (43*nloop/50.eq.n) nomvect='vecteur43.txt'
         if (44*nloop/50.eq.n) nomvect='vecteur44.txt'
         if (45*nloop/50.eq.n) nomvect='vecteur45.txt'
         if (46*nloop/50.eq.n) nomvect='vecteur46.txt'
         if (47*nloop/50.eq.n) nomvect='vecteur47.txt'
         if (48*nloop/50.eq.n) nomvect='vecteur48.txt'
         if (49*nloop/50.eq.n) nomvect='vecteur49.txt'
         if (50*nloop/50.eq.n) nomvect='vecteur50.txt'
c
c =====
c Ecrire les donnees en format vector pour gnuplot
c x,y,z,deltax,deltay,deltaz
         if (nomvect.ne.'toto') then
           open(unit=3,file=nomvect,status="unknown")
             write(3,*) nm
             do ii=1,nm
               write(3,*) xa(ii),ya(ii),za(ii),vxa(ii)*1.e15,v
     +    ya(ii)*1.e15,vzia(ii)*1.e15
             enddo
           close(unit=3)
         endif
c
c =====
c store la valeur precedentes
c
           axa(i)=axn(i)
           aya(i)=ayn(i)
           aza(i)=azin(i)
           vxa(i)=vxn(i)
           vya(i)=vyn(i)
           vzia(i)=vzin(i)
           xa(i)=xn(i)
           ya(i)=yn(i)
           za(i)=zin(i)
         enddo
       enddo
c
c =====
c Ecrire les donnees de sortie
c
           write(2,*) nm
           do i=1,nm
             write(2,*) xa(i),ya(i),za(i),vxa(i),vya(i),vzia(i),m(i)
           enddo
       close(unit=2)
       stop
       end
