module modalgo

  !contient touts les algorithmes necessaires a la resolution d'un systeme lineaire Ax=b

  implicit none

  integer,parameter::pr=8

contains

!=======================================================================================================================================================

    !fonction pour la methode de gradient a  pas optimal

    function grad_pas_optimal(A,b,x0,kmax,e) result(x)

      !variables d'entree
      real(kind=Pr), intent(in):: e                     !precision
      integer, intent(in)::kmax                         !test d'arret
      real(kind=Pr), dimension(:,:), intent(in)::A      !matrice du systeme lineaire
      real(kind=Pr),dimension(:), intent(in):: b, x0    !second membre + vecteur initial

      !sortie
      real(kind=Pr), dimension(:),allocatable::x        !solution approchee

      !variables locales
      real(kind=Pr), dimension(:), allocatable :: r   !residu reel
      real(kind=Pr), dimension(:), allocatable :: z
      real(kind=Pr) :: alpha, beta
      integer :: k, n

      !initialisation
      n=size(A(:,1))
      allocate(z(n))
      allocate(r(n))
      allocate(x(n))

      x=x0
      r=b-matmul(A,x0)
      beta=NORM2(r)
      k=0


      !algorithme de la methode 
      do while (k<=kmax .and.beta>e )
         z=matmul(A,r)
         alpha=DOT_PRODUCT(r,r)/DOT_PRODUCT(z,r)
         x=x+alpha*r
         r=r-alpha*z
         beta=NORM2(r)
         k=k+1

      end do

      if (k>kmax) then
         print*, "tolerance non atteinte", NORM2(r)
      else
         print*, "convergence en :", k, "iterations"
      end if

      deallocate(r,z)
      
    end function grad_pas_optimal
      
      
!==============================================================================================================================================    
    !fonction pour la methode du residu minimum

    function res_min(A,b,x0,kmax,e)result(x)

      !Declaration des variables
      real(kind=Pr), intent(in):: e                      !precision
      integer, intent(in)::kmax                          !test d'arret
      real(kind=Pr), dimension(:,:), intent(in)::A       !matrice du syteme lineaire
      real(kind=Pr),dimension(:), intent(in):: b, x0     !second membre + vecteur initial

      !sortie
      real(kind=Pr), dimension(:),allocatable::x         !solution approchee 

      !variables locales
      real(kind=Pr), dimension(:), allocatable :: r,z
      real(kind=Pr) :: alpha, beta
      integer :: k, n


      !initialisation
      n=size(A(:,1))
      allocate(z(n),r(n),x(n))
      x=x0
      r=b-matmul(A,x0)
      beta=NORM2(r)
      z=0._pr
      k=0

      !algorithme de la methode
      do while (k<=kmax .and. beta> e)    !Fonction norme Ã  dÃ©finir
         z=matmul(A,r)
         alpha=DOT_PRODUCT(r,z)/DOT_PRODUCT(z,z)
         x=x+alpha*r
         r=r-alpha*z
         beta=NORM2(r)
         k=k+1

      end do

      

      if (k>kmax) then
         print*, "tolerance non atteinte", NORM2(r)
      else
         print*, "convergence en :", k, "iterations"
      end if

      deallocate(z,r)
  end function res_min



!=====================================================================================================================================================

  
!fonction pour la methode du gradient conjugue
  
function gradient_conjugue(A,b,x0,kmax,e) result(x)
    
    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice de résolution du probleme
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision

    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchée du systeme 

    !variables locales
    real(kind=pr),dimension(:),allocatable:: r, rplus, p, z
    real(kind=pr):: beta, alpha, gamma, w
    integer:: k, n

    !initialisation
  
    n=size(A(:,1))
    !taille du systeme
    allocate(x(n),r(n),p(n),z(n),rplus(n))
    x=x0
    r=b-MATMUL(A,x0)
    beta=NORM2(r)
    p=r
    
    k=0

    !algorithme de la methode 
    do while (beta>e.and.k <= kmax)
       z=MATMUL(A,p)
       w=DOT_PRODUCT(r,r)
       alpha=w/DOT_PRODUCT(z,p)
       x=x+alpha*p
       rplus=r-alpha*z
       gamma=DOT_PRODUCT(rplus,rplus)/w
       p=rplus+gamma*p
       r=rplus
       beta=NORM2(r)
       k=k+1
    end do

    

    if(k>kmax) then
       print*, "tolerance non atteinte:" ,  beta
    else
       print*, "convergence en :", k, "iterations"
    end if
    
    deallocate(r,rplus,p,z)

    
  end function gradient_conjugue

!===============================================================================================================================================

!fonction qui contient la methode de resolution FOM

  
  function FOM(A,b,x0,kmax,e,m) result(x)

    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice du systeme lineaire
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    integer,intent(in):: m                            !caracterise espace de krylov dans lequel on se place

    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchée du systeme

    !variables locales
    real(kind=pr),dimension(:),allocatable:: r        !residu reel
    real(kind=pr),dimension(:),allocatable:: y, u, betae1
    real(kind=pr):: beta, somme
    integer:: k, n, i, p
    real(kind=pr),dimension(:,:),allocatable:: Hm, Vmplus    !matrice de Hessenberg et matrice de la nouvelle base orthogonal extrait de A à partir du résidu
    real(kind=pr),dimension(:,:),allocatable:: Qm, Rm    ! matrice de decomposition QR de Hmbarre 

    !initialisation
    n=size(A(:,1))
    allocate(r(n),x(n))
    x=x0
    r=b-MATMUL(A,x0)     !residu initial
    beta=NORM2(r)
    
    allocate(y(m))
    y=0._pr
    
    allocate(Hm(m+1,m))
    allocate(Vmplus(n,m+1))
    allocate(Qm(m,m),Rm(m,m)) 

    allocate(u(m),betae1(m))
    u=0._pr
    betae1=0._pr
    betae1(1)=beta

    k=0
    
    !algorithme de la methode
    
    do while (beta>e.and.k<= kmax)

       Hm=0._pr
       Vmplus=0._pr
       Qm=0._pr
       Rm=0._pr

       betae1=0._pr
       betae1(1)=beta
       
       !obtention de Hm et Vm par la methode d'Arnoldi en partant de r et A
       call Arnoldi(r,A,Hm,Vmplus)
       
       !resolution de Hmbarre*y=beta*e1

       !obtention de la decomposition QR de Hmbarre pour pouvoir réaliser la résolution
       call QR(Hm(1:m,1:m),Qm,Rm)

       !resolution de Qm*u=betae1
       !Q appartient au groupe orthogonal donc TQmQm=I donc u=TQm*betae1
       u=MATMUL(transpose(Qm),betae1)  
       !resolution de Rm*y=u
       !Rm est triangulaire superieure a diagonale non nulle
       !methode de descente
       y(m)=1._pr/Rm(m,m)*u(m)
       do i=m-1,1,-1
          somme=0._pr
          do p=i+1,m
             somme=somme+Rm(i,p)*y(p)
          end do
          y(i)=(u(i)-somme)/Rm(i,i)
       end do
        
       x=x+MATMUL(Vmplus(1:n,1:m),y)   
       !r=-Hm(m+1,m)*y(m)*Vmplus(1:n,m+1)
       r=r-MATMUL(A,MATMUL(Vmplus(1:n,1:m),y))
       beta=NORM2(r)
       k=k+1    
    end do
    
    if (k>kmax) then
       print*, 'tolerance non atteinte' , beta
    else
       print*, "convergence en :", k, "iterations"
    end if

    deallocate(r,y,Hm,Vmplus,Qm,Rm,u,betae1)
    
  end function FOM

!=============================================================================================================================================

!fonction pour la methode GMRes

  
  function GMRes(A,b,x0,kmax,e,m)  result(x)

    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice de résolution du probleme
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    integer,intent(in):: m

    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchée du systeme

    !variables locales
    real(kind=pr),dimension(:),allocatable:: r, y, betae1, u
    real(kind=pr):: beta
    integer:: n,k, p, i
    real(kind=pr):: somme
    real(kind=pr),dimension(:,:),allocatable:: Hm, Vmplus, Qm, Rm
    
    
    !initialisation
    n=size(A(:,1))
    
    allocate(x(n),r(n))
    x=x0
    r=b-MATMUL(A,x0)
    beta=NORM2(r)
   
    allocate(betae1(m),u(m))
    u=0._pr
     !definition de betae1
    betae1=0._pr
    betae1(1)=beta

    allocate(y(m))
    y=0._pr
    
    allocate(Hm(m+1,m))
    allocate(Vmplus(n,m+1))
    allocate(Rm(m,m),Qm(m,m))
    
    k=0

    !algorithme de la  methode

    do while (beta>e.and.k<= kmax)

       
       Hm=0._pr
       Vmplus=0._pr
       Qm=0._pr
       Rm=0._pr

       betae1=0._pr
       betae1(1)=beta
       !obtention de Hm et Vm+1
       call Arnoldi(r,A,Hm,Vmplus)

       !calcul de y=argmin(beta e1 - Hmy)
       !resolution de Hmbarre*y=beta*e1
       !obtention de la decomposition QR de Hmbarre pour pouvoir réaliser la résolution
       call QR(Hm(1:m,1:m),Qm,Rm)

       !resolution de Qm*u=betae1
       !Q appartient au groupe orthogonal donc TQmQm=I donc u=TQm*betae1
       u=MATMUL(transpose(Qm),betae1)
       
       !resolution de Rm*y=u
       !Rm est triangulaire superieure a diagonale non nulle
       !methode de descente
       y(m)=1._pr/Rm(m,m)*u(m)
       do i=m-1,1,-1
          somme=0._pr
          do p=i+1,m
             somme=somme+Rm(i,p)*y(p)
          end do
          y(i)=(u(i)-somme)/Rm(i,i)
       end do
       
       x=x+MATMUL(Vmplus(1:n,1:m),y)  
       r=r-MATMUL(A,MATMUL(Vmplus(1:n,1:m),y))
       beta=NORM2(betae1-MATMUL(Hm(1:m,1:m),y))
       k=k+1
       
    end do

    
    if (k>kmax) then
       print*, 'tolerance non atteinte', beta
    else
       print*, "convergence en :", k, "iterations"
    end if

    deallocate(r,betae1,y,Hm,Vmplus,Rm,Qm,u)
    
  end function GMRes

!==================================================================================================================================================
  !fonction qui contient la methode de resolution FOM dans le cas particulier où la matrice du ssytéme linéaire est sdp

  
  function FOM_sdp(A,b,x0,kmax,e,m) result(x)

    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice du systeme lineaire
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    integer,intent(in):: m                            !caracterise espace de krylov dans lequel on se place

    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchée du systeme

    !variables locales
    real(kind=pr),dimension(:),allocatable:: r        !residu reel
    real(kind=pr),dimension(:),allocatable:: y, u, betae1
    real(kind=pr):: beta, somme
    integer:: k, n, i, p
    real(kind=pr),dimension(:,:),allocatable:: Hm, Vmplus    !matrice de Hessenberg et matrice de la nouvelle base orthogonal extrait de A à partir du résidu
    real(kind=pr),dimension(:,:),allocatable:: L, tL    ! matrice de decomposition de cholesky de Hmbarre 

    !initialisation
    n=size(A(:,1))
    allocate(r(n),x(n))
    x=x0
    r=b-MATMUL(A,x0)     !residu initial
    beta=NORM2(r)
    
    allocate(y(m))
    y=0._pr
    
    allocate(Hm(m+1,m))
    allocate(Vmplus(n,m+1))

    
    allocate(L(m,m),tL(m,m))
    
    allocate(u(m),betae1(m))
    u=0._pr
    betae1=0._pr
    betae1(1)=beta

    k=0
    
    !algorithme de la methode
    
    do while (beta>e.and.k<= kmax)

       Hm=0._pr
       Vmplus=0._pr
       L=0._pr
       tL=0._pr

       betae1=0._pr
       betae1(1)=beta
       !obtention de Hm et Vm par la methode d'Arnoldi ameliore dans le cas ou A est sdp en partant de r et A
       call Arnoldi_tridiag(r,A,Hm,Vmplus)

       !obtention de la decomposition de Cholesky amelioree de Hm
       ! L=cholesky_tridiag(Hm(1:m,1:m))
       L=cholesky_tridiag(Hm(1:m,1:m))
       
       !resolution de Hmbarre*y=beta*e1
       !resolution de Lu=betae1
       !L est triangulaiure inferieure
       !methode de montee
       u(1)=1._pr/L(1,1)*betae1(1)
       do i=2,m
          somme=0._pr
          do p=1,i-1
             somme=somme+L(i,p)*u(p)
          end do
          u(i)=-somme/L(i,i)
       end do
       
       !resolution de tLy=u
       !tL est triangulaire superieure a diagonale non nulle
       !methode de descente
       tL=transpose(L)
       y(m)=1._pr/tL(m,m)*u(m)
       do i=m-1,1,-1
          somme=0._pr
          do p=i+1,m
             somme=somme+tL(i,p)*y(p)
          end do
          y(i)=(u(i)-somme)/tL(i,i)
       end do
        
       x=x+MATMUL(Vmplus(1:n,1:m),y)   
       ! r=-Hm(m+1,m)*y(m)*Vmplus(1:n,m+1)
       r=r-MATMUL(Vmplus,MATMUL(Hm,y))
       beta=NORM2(r)
       k=k+1    
    end do
    
    if (k>kmax) then
       print*, 'tolerance non atteinte' , beta
    else
       print*, "convergence en :", k, "iterations"
    end if

    deallocate(r,y,Hm,Vmplus,L,tL,u,betae1)
    
  end function FOM_SDP
!==================================================================================================================================================

!fonction pour la methode GMRes dans le cas ou la matrice du ssyteme lineaire est sdp

  
  function GMRes_sdp(A,b,x0,kmax,e,m)  result(x)

    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice de résolution du probleme
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    integer,intent(in):: m

    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchee du systeme

    !variables locales
    real(kind=pr),dimension(:),allocatable:: r, y, betae1, u
    real(kind=pr):: beta
    integer:: n, k, p, i
    real(kind=pr):: somme
    real(kind=pr),dimension(:,:),allocatable:: Hm, Vmplus, L, tL     
    
    
    !initialisation
    n=size(A(:,1))
    
    allocate(x(n),r(n))
    x=x0
    r=b-MATMUL(A,x0)
    beta=NORM2(r)
    !definition de betae1
    allocate(betae1(m),u(m))
    u=0._pr
    betae1=0._pr
    betae1(1)=beta
    
    allocate(y(m))
    y=0._pr
    
    allocate(Hm(m+1,m))
    allocate(Vmplus(n,m+1))
    allocate(L(m,m),tL(m,m))
    
    k=0

    !algorithme de la  methode

    do while (beta>e.and.k<= kmax)
       
       Hm=0._pr
       Vmplus=0._pr
       L=0._pr
       tL=0._pr

       betae1=0._pr
       betae1(1)=beta

       !obtention de Hm et Vm+1 par la methode d'Arnoldi amelioree dans le cas ou Hmbarre est une matrice tridiagonale
     
       call Arnoldi_tridiag(r,A,Hm,Vmplus)

       !obtention de la decomposition de Cholesky amelioree de Hm
       L=cholesky_tridiag(Hm(1:m,1:m))
     
       !calcul de y=argmin(beta e1 - Hmy)
       !resolution de Hmbarre*y=beta*e1

       !resolution de Lu=betae1
       !L est triangulaiure inferieure
       !methode de montee
       u(1)=1._pr/L(1,1)*betae1(1)
       do i=2,m
          somme=0._pr
          do p=1,i-1
             somme=somme+L(i,p)*u(p)
          end do
          u(i)=-somme/L(i,i)
       end do
  
       !resolution de tLy=u
       !tL est triangulaire superieure a diagonale non nulle
       !methode de descente
       tL=transpose(L)
       y(m)=1._pr/tL(m,m)*u(m)
       do i=m-1,1,-1
          somme=0._pr
          do p=i+1,m
             somme=somme+tL(i,p)*y(p)
          end do
          y(i)=(u(i)-somme)/tL(i,i)
       end do
       
       x=x+MATMUL(Vmplus(1:n,1:m),y)  
       r=r-MATMUL(Vmplus,MATMUL(Hm,y))
       beta=NORM2(r)
       k=k+1
       
    end do

    
    if (k>kmax) then
       print*, 'tolerance non atteinte', beta
    else
       print*, "convergence en :", k, "iterations"
    end if

    deallocate(r,betae1,y,Hm,Vmplus,u,L,tL)
    
  end function GMRes_Sdp

!===================================================================================================================================================

  !algorithme d'Arnoldi: methode d'orthogonalisation de Gram-Schmidt modifiee
  !obtention de la matrice Vmplus qui contient tout les nouveaux vecteurs de la nouvelle base
  
  subroutine Arnoldi(r,A,Hm,Vmplus)

    !variables d'entrees
    real(kind=pr),dimension(:,:),intent(in):: A
    real(kind=pr),dimension(:),intent(in)::r                !vecteur qui definit l'espace de Krylov {v,Av,.....,A(m-1)v}
   
    !variables de sortie
    real(kind=pr),dimension(:,:),intent(inout):: Vmplus    !matrice qui contient les vecteurs de la nouvelle base orthonormale
    real(kind=pr),dimension(:,:),intent(inout):: Hm     !matrice de Hessenberg qui contient les coefficients de la methode de Gram-Schmidt
    
    !variables locales
    integer:: n, m, i, j
    real(kind=pr),dimension(:),allocatable:: wj, vj
   
    !initialisation
    n=size(A(:,1))
    m=size(Hm(1,:))
    
    allocate(wj(n),vj(n))
    Hm=0._pr
    Vmplus=0._pr

    !algorithme d'Arnoldi

    Vmplus(1:n,1)=r/NORM2(r)
    do j=1,m
       vj=Vmplus(:,j)
       wj=MATMUL(A,vj)
       do i=1,j
          Hm(i,j)=DOT_PRODUCT(wj,Vmplus(:,i))
          wj=wj-Hm(i,j)*Vmplus(:,i)
       end do
      
       Hm(j+1,j)=NORM2(wj)

       if( Hm(j+1,j)==0) then
          stop
       end if
       Vmplus(:,j+1)=wj/Hm(j+1,j)
       
    end do

    deallocate(wj,vj)
    
  end subroutine Arnoldi

!==========================================================================================================================================
  !contient l'algorithme d'Arnoldi ameliore dans les cas ou la matrice obtenue est tridiagonale+derniere ligne
   subroutine Arnoldi_tridiag(r,A,Hm,Vmplus)

    !variables d'entrees
    real(kind=pr),dimension(:,:),intent(in):: A
    real(kind=pr),dimension(:),intent(in)::r                !vecteur qui definit l'espace de Krylov {v,Av,.....,A(m-1)v}
   
    !variables de sortie
    real(kind=pr),dimension(:,:),intent(inout):: Vmplus    !matrice qui contient les vecteurs de la nouvelle base orthonormale
    real(kind=pr),dimension(:,:),intent(inout):: Hm     !matrice de Hessenberg qui contient les coefficients de la methode de Gram-Schmidt
    
    !variables locales
    integer:: n, m, i, j
    real(kind=pr),dimension(:),allocatable:: wj, vj
   
    !initialisation
    n=size(A(:,1))
    m=size(Hm(1,:))
    
    allocate(wj(n),vj(n))
    Hm=0._pr
    Vmplus=0._pr

    !algorithme d'Arnoldi

    Vmplus(1:n,1)=r/NORM2(r)
    do j=1,m
       vj=Vmplus(:,j)
       wj=MATMUL(A,vj)
       if (j==1) then
          Hm(j,j)=DOT_PRODUCT(wj,Vmplus(:,j))
          wj=wj-Hm(j,j)*Vmplus(:,j)
       else
          do i=j-1,j
             Hm(i,j)=DOT_PRODUCT(wj,Vmplus(:,i))
             wj=wj-Hm(i,j)*Vmplus(:,i)
          end do
       end if
      
       Hm(j+1,j)=NORM2(wj)

       if( Hm(j+1,j)==0) then
          stop
       end if
       Vmplus(:,j+1)=wj/Hm(j+1,j)
       
    end do
    
    deallocate(wj,vj)
    
  end subroutine Arnoldi_Tridiag

  

!===========================================================================================================================================

  !fonction qui contient la decomposition polaire matrice A de taille m*m dans le cas particulier d'une matrice de Hessenberg 
  !la methode de Givens est utilisee 

  subroutine QR(A,Q,R)

    !variables d'entree et de sortie
    real(kind=pr),dimension(:,:),intent(in):: A           !matrice dont on effectue la decomposition A=QR
    real(kind=pr),dimension(:,:),intent(inout):: Q        !matrice du groupe orthogonal
    real(kind=pr),dimension(:,:),intent(inout):: R        !matrice triangulaire superieure

    !variables internes
    integer:: m, i, j, k
    real(kind=pr),dimension(:,:),allocatable:: P
    real(kind=pr),dimension(:,:),allocatable:: TQ
    real(kind=pr):: c, s


    !initialisation
    m=size(A(1,:))
    allocate(P(m,m),TQ(m,m))
    TQ=0._pr
    do i=1,m
       TQ(i,i)=1._pr
    end do
    
    P=0._pr
    R=0._pr

    c=0._pr
    s=0._pr

    !m-1 rotations a effectuer
    do k=1,m-1

       !remise à 0 de P
       P=0._pr

       !initialisation de la diagonale de 1
       do i=1, m
          P(i,i)=1._pr
       end do
       
       !remplisassage des matrices de rotation
       c=A(k,k)/sqrt(A(k,k)**2+A(k+1,k)**2)
       s=A(k+1,k)/sqrt(A(k,k)**2+A(k+1,k)**2)
       
       P(k+1,k+1)=c
       P(k,k)=c
       P(k+1,k)=-s
       P(k,k+1)=s

       TQ=MATMUL(P,TQ)

    end do

    
    !obtention de Q
    Q=transpose(TQ)

    !obtention de R
    R=MATMUL(TQ,A)

    deallocate(P,TQ)

  end subroutine QR

!===========================================================================================================================================
  !contient la decomposition de Cholesky d'une matrice tridiagonale sdp

  function chol(A )result(L)

     Implicit None

     !les variables

     integer                              :: i,j,k,n
     real(Pr),dimension(:,:),intent(in)   :: A
     real(Pr),dimension(:,:),allocatable  :: L
     real(Pr)                             :: s,p

     !initialisations
     n=size(A,1)
     Allocate(L(1:n,1:n))
     L=0._pr

     !calcul de L
     !boucle sur toutes les colonnes
     do i=1,n
        s=0._PR
        !element diagonal
        do k=1,i-1
           s= s + (L(i,k))**2
        end do

        !ici on vérifie que la matrice est sdp

        if (A(i,i)-s<=0) then
           print*, "attention la matrice n'est pas définie positive"
           print*, A(i,i)-s
        end if

        L(i,i)=sqrt(A(i,i)-s)

        !boucle sur tous les éléments extra-diagonaux
        do j=i+1,n
           p=0._PR
           do k=1,i-1
              p=p+ L(i,k)*L(j,k)
           end do
           L(j,i)=(A(j,i)-p)/L(i,i)
        end do
     end do

   end function chol


  function cholesky_tridiag(A) result(L)

    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A
    !variables de sortie
    real(kind=pr),dimension(:,:),allocatable:: L
    !variables locales
    integer:: m, j, k
    real(kind=pr):: somme1, somme2

    !initialisation
    m=size(A(1,:))
    allocate(L(m,m))
    L=0._pr
   
    !algorithme de decomposition remplissage de L
    
    do j=1,m-1
        somme1=0._pr
       do k=1,j-1
          somme1=somme1+L(j,k)**2
          somme2=somme2+L(j+1,j)*L(j,k)
       end do
       L(j,j)=sqrt(A(j,j)-somme1)
       somme2=0._pr
       do k=1,j-1
          somme2=somme2+ L(j,k)*L(j+1,k)
       end do
       L(j+1,j)=(A(j+1,j)-somme2)/L(j,j)
    end do
    !cas j=m
    somme1=0._pr
    do k=1,m-1
       somme1=somme1+L(m,k)**2
    end do
    L(m,m)=sqrt(A(m,m)-somme1)

  end function cholesky_tridiag

  
end module modalgo
