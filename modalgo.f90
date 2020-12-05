module modalgo

  !contient touts les algorithmes necessaire a la resolution d'un systeme lineaire Ax=b

  implicit none

  integer,parameter::pr=8

contains

  !produit matrice matrice
  !produit matrice vecteur
  !norme p
  !norme infinie
  !transposée d'une matrice
  !trace

    function gradient_conjugue(A,b,x0,kmax,e) result(x)

    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice de r�solution du probleme
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donn�e intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approch�e du systeme
    !variables locales
    real(kind=pr),dimension(:),allocatable:: r, p, z
    real(kind=pr):: beta, alpha
    integer:: k, n

    !initialisation

    n=size(A)                                          !taille du systeme
    allocate(x(n),r(n),p(n),z(n))
    x=0._pr
    beta=3._pr
    r=MATMUL(A,x0)
    r=b-r

    do while (beta>e.and.k <= kmax)
       z=MATMUL(A,p)
       !alpha=

    end do



  end function gradient_conjugue



    !fonction pour la méthode de gradient à pas optimal

    function grad_pas_optimal(A,b,x0,kmax,e)result(x)

      !déclaration des arguments
      !précision
      real(kind=Pr), intent(in):: e
      !test d'arrêt
      integer, intent(in)::kmax
      real(kind=Pr), dimension(:,:), intent(in)::A
      real(kind=Pr),dimension(:), intent(in):: b, x0
      real(kind=Pr), dimension(:),allocatable, intent(out)::x
      real(kind=Pr), dimension(:), allocatable :: r,z
      real(kind=Pr) :: alpha
      integer :: k, n

      n=size(A)
      Allocate(A(n))
      Allocate(z(n))
      Allocate(x(n))

      r=b-matmul(A,x0)
      k=0

      do while (k<=kmax .and. )
         z=matmul(A,r)
         alph a= dotproduct(r,r)/dotproduct(z,r)
         x=x+alpha*r
         r=r-alpha*z
         k=k+1

      end do
      return x

      if (k>kmax) then
         return "tolérance atteinte"
      end if
      end function

!====================================================================================    
    !fonction pour la mÃ©thode du rÃ©sidu minimum

    function res_min(A,b,x0,kmax,e) result(x)

      !DÃ©claration des variables
      real(kind=Pr), intent(in):: e
      !ItÃ©ration maximale
      integer, intent(in)::kmax
      real(kind=Pr), dimension(:,:), intent(in)::A
      real(kind=Pr),dimension(:), intent(in):: b, x0
      real(kind=Pr), dimension(:),allocatable, intent(out)::x
      real(kind=Pr), dimension(:), allocatable :: r,z
      real(kind=Pr) :: alpha
      integer :: k, n

      n=size(A)
      Allocate(A(n))
      Allocate(z(n))
      Allocate(x(n))

      r=b-matmul(A,x0)
      k=0

      do while (k<=kmax .and. NORM2(r) > e)    !Fonction norme Ã  dÃ©finir
         z=matmul(A,r)
         alpha= dotproduct(r,z)/dotproduct(z,z)
         x=x+alpha*r
         r=r-alpha*z
         k=k+1

      end do
      return x

      if (k>kmax) then
         return "tolÃ©rance atteinte", NORM2(r)
      end if
  end function res_min






  !===================================================================================
function gradient_conjugue(A,b,x0,kmax,e) result(x)
    
    !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice de résolution du probleme
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchée du systeme 
    !variables locales
    real(kind=pr),dimension(:),allocatable:: r, rplus,  p, z
    real(kind=pr):: beta, alpha, gamma, w
    integer:: k, n

    !initialisation
  
    n=size(A)                                          !taille du systeme
    allocate(x(n),r(n),p(n),z(n),rplus(n))
    x=x0
    r=b-MATMUL(A,x0)
    beta=NORM2(r)
    p=r
    k=0

    !algorithme 
    do while (beta>e.and.k <= kmax)
       z=MATMUL(A,p)
       w=DOT_PRODUCT(r,r)
       alpha=w/DOT_PRODUCT(z,p)
       x=x-alpha*p
       rplus=r-alpha*z
       gamma=DOT_PRODUCT(rplus,rplus)/w
       p=rplus+gamma*p
       beta=NORM2(r)
       r=rplus
       k=k+1
    end do

    if(k>kmax) then
       print*, "tolérance non atteinte:" ,  beta
    end if
    
    deallocate(r,rplus,p,z)

    
  end function gradient_conjugue

  

  subroutine Arnoldi(v,A,Hm,Vm)

    !variables d'entrees 
    real(kind=pr),dimension(:,:),intent(in)::v
    real(kind=pr),dimension(:,:),intent(in):: A
    !variable de sortie
    real(kind=pr),dimension(:,:),allocatable,intent(out):: Vm
    real(kind=pr),dimension(:,:),allocatable,intent(out):: Hm
    !variables locales
    integer:: m, n, i, j
    real(kind=pr),dimension(:),allocatable:: wj

    !initialisation
    n=size(v(:,1))
    m=size(v(1,:))
    allocate(Hm(m+1,m))
    allocate(Vm(n,m))
    allocate(wj(n))
    Hm=0._pr

    !algorithme d'Arnoldi

    Vm(:,1)=v(:,1)
    do j=1,m
       wj=MATMUL(A,v(:,j))
       do i=1,j
          Hm(i,j)=DOT_PRODUCT(wj,v(:,i))
          wj=wj-Hm(i,j)*v(:,i)
       end do
       Hm(j+1,j)=NORM2(wj)

       if( Hm(j+1,j)==0) then
          stop
       end if
       Vm(:,j+1)=1._pr/Hm(j+1,j)*wj
       
    end do

    
  end subroutine Arnoldi


  function FOM(A,b,x0,kmax,e) result(x)
     !variables d'entree
    real(kind=pr),dimension(:,:),intent(in):: A       !matrice de résolution du probleme
    real(kind=pr),dimension(:),intent(in):: b, x0     !vecteur second membre et donnée intiale
    integer,intent(in):: kmax                         !test d'arret
    real(kind=pr),intent(in):: e                      !precision
    !variables de sortie
    real(kind=pr),dimension(:),allocatable:: x        !solution approchée du systeme
    !variables locales
    real(kind=pr),dimension(:),allocatable:: r
    real(kind=pr):: beta
    integer:: k

    !initialisation
    r=b-MATMUL(A,x0)
    beta=NORM2(r)
    k=0

    do while (beta>e.and.k <= kmax)
    end do
  end function FOM

  
 !====================================================================================
     !Fonction pour le calcul de la normeInf
     
    function norme_inf (a) result(norme)
        implicit none
        ! --- arguments
        type (element), dimension(:), intent(in) :: a
        real :: norme
        ! --- variables locales
        integer :: i, taille
        real, dimension (:) , allocatable :: y
        ! --- calcul taille de la matrice pleine associee
        taille=0
            do i=1, size (a)
                taille=max (taille, a(i)%indl, a (i)%indc )
            end do
        ! --- calcul de la norme
        norme=0.
        allocate (y(1:taille))
        y=0.
            do i=1, size (a)
                y(a(i)%indl)=y (a(i)%indl)+abs (a(i)%coef)
            end do
        norme=maxval(y)
        deallocate(y)
    end function norme_inf


end module modalgo
