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
  
  !initialisation
  n=3
  allocate(A(n,n),b(n),x0(n),x(n))
  A(1,1)=8._pr
  A(2,2)=7._pr
  A(3,3)=6._pr
  A(1,2)=1._pr
  A(2,1)=1._pr
  A(2,3)=2._pr
  A(3,2)=2._pr
  A(3,1)=3._pr
  A(1,3)=3._pr


  b(1)=1._pr
  b(2)=2._pr
  b(3)=3._pr


  x0(1)=6._pr
  x0(2)=5._pr
  x0(3)=4._pr

  e=0.001_pr
  kmax=10000
  m=2
  

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
  x=GMRes(A,b,x0,kmax,e,m)
  print*, x
  
  
  deallocate(A,b,x0,x)
end program main
