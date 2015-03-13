c programme pour creer une galaxie spherique aleatoire
c declaration des variables
       double precision alpha(1000000),phi(1000000),alph,x,y,deltar,
     + a,f,g
       double precision alphaf(1000000),phif(1000000),posx(1000000),
     + posy(1000000)
       double precision posz(1000000),vx(1000000),vy(1000000),x0,z0,
     + vz(1000000)

       double precision xtmp,ytmp,ztmp,v0
c martin <
       double precision vmax,v,facteur,masse,rayonact,rayon,massei
c > martin
       integer k,i,j,n,jmax,nbre,nn,nmax,m,pos,ii,nbmasse,nbgal,ng,ne
c      a est un parametre fixe arbitrairement qui influe sur la concentration des etoiles au centre de la galaxie.
c      nn est le nombre total d'etoiles dans la galaxie
c      deltar est la largeur d'une couche de la galaxie
c      k permet de faire varier le nombre de couches de la galaxie
c      posx, posy et posz sont les positions en coordonnees cartesiennes des etoiles
c      m est un compteur qui compte le nombre d'etoiles total
c      vmax est la vitesse maximale d'une etoile
       masse=4.e+37
       a=19.46e+19

       deltar=4.73e+17
       facteur=1000000000000000000.
          m=0
        ne=0
         open(unit=1,file='3d',status='unknown')
         open(unit=2,file='vector',status='unknown')
         print*,'Nombre de galaxies?'
         read*,nbgal
         do ng=1,nbgal
         print*,'Position initiale de la galaxie ',ng,' en x et z'
         read*,x0,z0
         print*,'Vitesse initiale de la galaxie ',ng,' en z'
         read*,v0
       nn=1000
         do k=1,1000
             n=0
             nbre=dnint((deltar*dble(nn)*a)/((a+(dble(k)*deltar))**2.))
             nmax=dnint(sqrt(dble(nbre)/0.5914))
             do i=1,100
c avec une boucle de 100 et 100 cos alpha on obtiens 6227 valeurs     
c On tire ici la banque des angles pour la position   
             call RANDOM_NUMBER(x)
             alph=(x-0.5)*3.1415926
             jmax=dnint(100.*cos(alph))
             do j=1,jmax
                n=n+1
                alpha(n)=alph
                call RANDOM_NUMBER(y)
                phi(n)=(y*6.2831853)
c                print*,alpha(n),phi(n)       
            enddo        
            enddo
c       print*,'n=',n
       do i=1,nbre
          call RANDOM_NUMBER(f)
          pos=dnint(f*dble(n))
          m=m+1
          alphaf(m)=alpha(pos)
          phif(m)=phi(pos)
c          print*,k,m, alphaf(m),phif(m)
          posx(m)=dble(k)*deltar*cos(alphaf(m))*cos(phif(m))
          posy(m)=dble(k)*deltar*cos(alphaf(m))*sin(phif(m))
          posz(m)=dble(k)*deltar*sin(alphaf(m))
c          write(1,*) posx(m),posy(m),posz(m)
c martin <     enddo
c             do i=1,m
c   > martin


c martin <
          pos=dnint(g*dble(n))
          alphaf(m)=alpha(pos)
          phif(m)=phi(pos)
c > martin

c           alphaf(m)=alpha(v)
c           phif(m)=phi(v)
c           vx(m)=dble(k)*deltar*cos(alphaf(m))*cos(phif(m))
c           vy(m)=dble(k)*deltar*cos(alphaf(m))*sin(phif(m))
c           vz(m)=dble(k)*deltar*sin(alphaf(m))
c           write(1,*) posx(m),posy(m),posz(m),vx(m),vy(m),vz(m)
c           print*,k,m, alphaf(m),phif(m),vx(m),vy(m),vz(m)
c martin <

           print*,m,'/',nn

c > martin
           enddo
       enddo
       nn=m
       m=0
       n=0
             do i=1,100
c avec une boucle de 100 et 100 cos alpha on obtiens 6227 valeurs        
             call RANDOM_NUMBER(x)
             alph=(x-0.5)*3.1415926
             jmax=dnint(100.*cos(alph))
             do j=1,jmax
                n=n+1
                alpha(n)=alph
                call RANDOM_NUMBER(y)
                phi(n)=(y*6.2831853)
c                print*,alpha(n),phi(n)       
             enddo        
             enddo
       do i=1,nn
           call RANDOM_NUMBER(f)
           pos=dnint(f*dble(n))
           m=m+1
          alphaf(m)=alpha(pos)
          phif(m)=phi(pos)
           call RANDOM_NUMBER(g)
          rayonact=sqrt(posx(m)**2.+posy(m)**2.+posz(m)**2.)
          nbmasse=0
          do ii=1,nn
             rayon=sqrt(posx(ii)**2.+posy(ii)**2.+posz(ii)**2.)
             if (rayon.lt.rayonact) then
                nbmasse=nbmasse+1
             endif
          enddo
          massei=dble(nbmasse)*masse
          vmax=sqrt((2.*6.67428e-11*massei)/rayonact)
          v=g*vmax
           vx(m)=v*cos(alphaf(m))*cos(phif(m))
           vy(m)=v*cos(alphaf(m))*sin(phif(m))
           vz(m)=v*sin(alphaf(m))+v0
           xtmp=posx(m)+x0
           ytmp=posy(m)
           ztmp=posz(m)+z0
           write(2,*) xtmp,ytmp,ztmp,vx(m)*facteur,
     +     vy(m)*facteur,vz(m)*facteur
           write(1,*) xtmp,ytmp,ztmp,vx(m),vy(m),vz(m),masse
           ne=ne+1
           print*,'Etoile no.',ne
           print*,i,'/',nn



       enddo
           m=0
       enddo  
       close(unit=1)
       close(unit=2)
       stop
       end
