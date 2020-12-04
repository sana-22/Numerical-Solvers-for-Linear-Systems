module modalgo

  !contient touts les algorithmes necessaire a la resolution d'un systeme lineaire Ax=b

  implicit none

  integer,parameter::pr=8

contains


!calcul le produit scalaire de deux vecteurs
  function produit_scalaire(u,v) result(w)
    implicit none
    !variables d'entrees
    real(kind=pr),dimension(:),intent(in):: u, v
    real(kind=pr),dimension(:),allocatable:: w
  end function produit_scalaire

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

end module modalgo
