program main

  use modalgo

  implicit none

  !données d'entree du probléme 
  real(kind=pr),dimension(:,:),allocatable:: A           !matrice du sytème lineaire
  real(kind=pr),dimension(:),allocatable:: b             !vecteur second membre
  real(kind=pr),dimension(:),allocatable:: x0            !vecteur initial de depart
  real(kind=pr):: e                                      !precision voulue
  integer:: kmax, m, n                                   !kmax test d'arret des methodes et m caractérisant espace de Krylov n dimension du syteme à resoudre

  !variables
  real(kind=pr),dimension(:),allocatable:: x
  real(kind=pr),dimension(:,:),allocatable:: Bn, An, Tn
  real(kind=pr),dimension(:),allocatable:: tab
  real(kind=pr):: alpha 
  integer:: i, j

!==================================================================================================================================
  !premiers tests
  !initialisation

  e=0.0001
  kmax=1000000000
  m=2
  n=3
  allocate(A(n,n),b(n),x0(n),x(n))

  A(1,1)=2._pr
  A(2,2)=2._pr
  A(3,3)=2._pr
  
  A(1,2)=1._pr
  A(2,1)=1._pr
  A(2,3)=1._pr
  A(3,2)=1._pr
  A(3,1)=1._pr
  A(1,3)=1._pr


  b(1)=1._pr
  b(2)=2._pr
  b(3)=3._pr


  x0(1)=1._pr
  x0(2)=1._pr
  x0(3)=1._pr
  !test des methodes
  print*, "test des methodes de resolution Ax=b  pour une matrice simple"
  print*, "A:", A
  print*, "b:", b
  print*, "x0:", x0
  print*,"precision:", e
  print*, "m:", m
  print*, "methode du gradient a pas optimal"
  x=grad_pas_optimal(A,b,x0,kmax,e)
  print*, x

  print*, "methode du residu minimum"
  x=res_min(A,b,x0,kmax,e)
  print*, x

  print*, "methode du gradient conjugue"
  x=gradient_conjugue(A,b,x0,kmax,e)
  print*, x
  
  print*, "methode FOM"
  x=FOM(A,b,x0,kmax,e,m)
  print*, x

  print*, "methode GMRes"
 ! x=GMRes(A,b,x0,kmax,e,m)
 ! print*, x 
  
  deallocate(A,b,x0,x)


!====================================================================================================================================================
  !tests sur An=In-alpha*tBn* Bn

  print*, "saisir n la taille du systeme:"
  read*, n
  print*, "saisir m:"
  read*, m

  
  e=0.0001
  kmax=1000000000


  allocate(An(n,n),Bn(n,n),Tn(n,n),x0(n),x(n),b(n),tab(n*n))

  An=0._pr
  Bn=0._pr
  !remplissage de Bn
  do i=1,n
     do j=1,n
        Bn(i,j)=rand(0)
     end do
  end do


  !Tn=tBn Bn
  Tn=MATMUL(transpose(Bn),Bn)
  
  !alpha=max(tBnBN)(i,j)
  !remplissage de tab
  do i=1,n
     do j=1,n
        tab(i+(j-1)*n)=Tn(i,j)
     end do
  end do

  alpha=MAXVAL(tab)
    
  !remplissage de An à partir de Bn
  do i=1,n
     An(i,i)=1._pr
  end do

  An=alpha*Tn+An
 
  b=1._pr
  x0=0._pr

  print*, "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"
  print*, "test des methodes pour An=In+alpha*tBn Bn "
  print*,"precision:", e
  print*, "m:", m
  print*, "methode du gradient a pas optimal"
  x=grad_pas_optimal(An,b,x0,kmax,e)
  print*, x

  print*, "methode du residu minimum"
  x=res_min(An,b,x0,kmax,e)
  print*, x

  print*, "methode du gradient conjugue"
  x=gradient_conjugue(An,b,x0,kmax,e)
  print*, x
  
  print*, "methode FOM"
  x=FOM(An,b,x0,kmax,e,m)
  print*, x

  print*, "methode GMRes"
 ! x=GMRes(An,b,x0,kmax,e,m)
 ! print*, x
  

  deallocate(An,Bn,Tn,x0,x,b)
  
end program main
