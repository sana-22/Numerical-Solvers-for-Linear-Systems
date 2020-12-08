module modalgo

  !contient touts les algorithmes necessaires a la resolution d'un systeme lineaire Ax=b

  implicit none

  integer,parameter::pr=8

contains

!=========================================================================================================================================

    !fonction pour la methode de gradient a  pas optimal

    function grad_pas_optimal(A,b,x0,kmax,e) result(x)

      !declaration des arguments
      real(kind=Pr), intent(in):: e                     !precision
      integer, intent(in)::kmax                         !test d'arret
      real(kind=Pr), dimension(:,:), intent(in)::A      !matrice du systeme lineaire
      real(kind=Pr),dimension(:), intent(in):: b, x0    !second membre + vecteur initial
      real(kind=Pr), dimension(:),allocatable::x        !solution approchee

      !varaibles locales
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
      
      
!===================================================================================================================================    
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


  !===================================================================================================================================================

!algorithme d'Arnoldi: construction de Vm et Hm
  
  subroutine Arnoldi(r,A,m,Hm,Vm)

    !variables d'entrees
    real(kind=pr),dimension(:,:),intent(in):: A
    real(kind=pr),dimension(:),intent(in)::r                !vecteur qui definit l'espace de Krylov {v,Av,.....,A(m-1)v}
    integer,intent(in)::m
   
    !variables de sortie
    real(kind=pr),dimension(:,:),allocatable,intent(out):: Vm     !matrice qui contient les vecteurs de la nouvelle base orthonormale
    real(kind=pr),dimension(:,:),allocatable,intent(out):: Hm     !matrice de Hessenberg qui contient les coefficients de la methode de Gram-Schmidt
    
    !variables locales
    integer:: n, i, j
    real(kind=pr),dimension(:),allocatable:: wj, vj
   
    !initialisation
    n=size(A(:,1))
    allocate(Hm(m+1,m))
    allocate(Vm(n,m+1))
    allocate(wj(n),vj(n))
    Hm=0._pr
    Vm=0._pr

    !algorithme d'Arnoldi

    Vm(:,1)=1._pr/NORM2(r)*r
    do j=1,m
       vj=Vm(1:n,j)
       wj=MATMUL(A,vj)   
       do i=1,j
          Hm(i,j)=DOT_PRODUCT(wj,Vm(1:n,i))
          wj=wj-Hm(i,j)*Vm(1:n,i)
       end do
       Hm(j+1,j)=NORM2(wj)

       if( Hm(j+1,j)==0) then
          stop
       end if
       Vm(1:n,j+1)=1._pr/Hm(j+1,j)*wj
       
    end do
    deallocate(wj)
    
  end subroutine Arnoldi


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
    real(kind=pr),dimension(:,:),allocatable:: Hm, Vm    !matrice de Hessenberg et matrice de la nouvelle base orthogonal extrait de A à partir du résidu
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
    allocate(Vm(n,m+1))
    Hm=0._pr
    Vm=0._pr

    
    allocate(Qm(m,m),Rm(m,m))
    Qm=0._pr
    Rm=0._pr
    allocate(u(m),betae1(m))
    u=0._pr
    betae1=0._pr
    betae1(1)=beta

    k=0

    call Arnoldi(r,A,m,Hm,Vm)
    call QR(Hm(1:m,1:m),Qm,Rm)
    print*, Qm(:,1)
    print*, Qm(:,2)
    print*, Qm(:,3)
    
    !algorithme de la methode
    
    do while (beta>e.and.k<= kmax)
       
       !obtention de Hm et Vm par la methode d'Arnoldi en partant de r et A
       call Arnoldi(r,A,m,Hm,Vm)
       
       !resolution de Hmbarre*y=beta*e1
       
       !obtention de la decomposition QR de Hmbarre pour pouvoir réaliser la résolution
       call QR(Hm(1:m,1:m),Qm,Rm)

       !resolution de Qm*u=betae1
       !Q appartient au groupe orthogonal donc TQmQm=I donc u=TQm*betae1
       u=MATMUL(transpose(Qm),betae1)
       
       !resolution de Rm*y=u
       !Rm est triangulaire superieure a diagonale non nulle
       !methode de remontee
       y(m)=1._pr/Rm(m,m)*u(m)
       do i=m-1,1,-1
          somme=0._pr
          do p=0,i-1
             somme=somme+Rm(i,m-p)*y(m-p)
          end do
          y(i)=1._pr/Rm(i,i)*(u(i)-somme)
       end do
       
       x=x+MATMUL(Vm(1:n,1:m),y)
       r=-Hm(m+1,m)*y(m)*Vm(1:n,m+1)
       beta=NORM2(r)
       k=k+1    
    end do
    
    if (k>kmax) then
       print*, 'tolerance non atteinte' , beta
    else
       print*, "convergence en :", k, "iterations"
    end if

    deallocate(r,y,Hm,Vm,Qm,Rm,u,betae1)
    
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
    real(kind=pr),dimension(:),allocatable:: r, y, betae1
    real(kind=pr):: beta
    integer:: n,k
    real(kind=pr),dimension(:,:),allocatable:: Hm, Vm
    
    !initialisation
    n=size(A(:,1))
    
    allocate(x(n),r(n),betae1(n))
    x=x0
    r=b-MATMUL(A,x0)
    beta=NORM2(r)

    betae1=0._pr
    
    allocate(y(m))
    y=0._pr
    
    allocate(Hm(m+1,m))
    allocate(Vm(n,m+1))
    Hm=0._pr
    Vm=0._pr
    
    k=0

    !algorithme de la  methode

    do while (beta>e.and.k<= kmax)
       call Arnoldi(r,A,m,Hm,Vm)
       !calcul de y=argmin(beta e1 - Hmy)

       x=x+MATMUL(Vm(1:n,1:m),y)
       !definition de betae1
       betae1(1)=beta
       r=betae1-MATMUL(Hm(1:m,1:m),y)
       beta=NORM2(r)
       k=k+1
       
    end do

    
    if (k>kmax) then
       print*, 'tolerance non atteinte', beta
    else
       print*, "convergence en :", k, "iterations"
    end if

    deallocate(r,betae1,y,Hm,Vm)
    
  end function GMRes

  !===========================================================================================================================================

  !fonction qui contient la decomposition polaire matrice A de taille m*m dans le cas particulier d'une matrice de Hessenberg A=Hmbarre
  !methode de Givens

  subroutine QR(A,Q,R)

    !variables d'entree de sortie
    real(kind=pr),dimension(:,:),intent(in):: A
    real(kind=pr),dimension(:,:),intent(out):: Q
    real(kind=pr),dimension(:,:),intent(out):: R
    !variables internes
    integer:: m, i, j
    integer(kind=pr),dimension(:,:),allocatable:: P
    integer(kind=pr),dimension(:,:),allocatable:: TQ
    real(kind=pr):: c, s


    !initialisation
    m=size(A(1,:))

    allocate(P(m,m),TQ(m,m))
    TQ=0._pr
    P=0._pr

    c=0._pr
    s=0._pr

    do j=1,m
       P=0._pr
       do i=1,m
       P(i,i)=1._pr
          if (i<=j) then
             !définition des coefficients de la rotation
             c=A(i,i)/sqrt(A(i,i)**2+A(i,j)**2)
             s=-A(j,i)/sqrt(A(i,i)**2+A(i,j)**2)
             !remplissage de la matrice de rotation i, j qui annule le terme A(i,j)
             P(i,i)=c
             P(i,j)=-s
             P(j,i)=s
             P(j,j)=c
             TQ=MATMUL(P,TQ)
          end if
          
       end do

       Q=transpose(TQ)
       R=MATMUL(TQ,A)
    end do
    
    deallocate(P,TQ)
    
  end subroutine QR
  
!===========================================================================================================================================
     !Fonction pour le calcul de la normeInf
     
   ! function norme_inf (a) result(norme)
        !implicit none
        ! --- arguments
        !type (element), dimension(:), intent(in) :: a
        !real :: norme
        ! --- variables locales
       ! integer :: i, taille
        !real, dimension (:) , allocatable :: y
        ! --- calcul taille de la matrice pleine associee
        !taille=0
            !do i=1, size (a)
              !  taille=max (taille, a(i)%indl, a (i)%indc )
           ! end do
        ! --- calcul de la norme
       ! norme=0.
       ! allocate (y(1:taille))
       ! y=0.
            !do i=1, size (a)
               ! y(a(i)%indl)=y (a(i)%indl)+abs (a(i)%coef)
          !  end do
       ! norme=maxval(y)
        !deallocate(y)
   ! end function norme_inf
    

end module modalgo
