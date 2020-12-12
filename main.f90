program main

  use modalgo

  implicit none

  !données d'entree du probléme 
  real(kind=pr),dimension(:,:),allocatable:: A           !matrice du sytème lineaire
  real(kind=pr),dimension(:),allocatable:: b             !vecteur second membre
  real(kind=pr),dimension(:),allocatable:: x0            !vecteur initial de depart
  real(kind=pr):: e                                      !precision voulue
  integer:: n                                            !dimension du syteme à resoudre
  integer:: m                                            !permet de definir un espace de Krylov d'une certaine taille
  integer:: kmax                                         !test d'arret des methodes
  !variables de resolution
  real(kind=pr),dimension(:),allocatable:: x
  real(kind=pr),dimension(:,:),allocatable:: Bn, An, Tn
  real(kind=pr),dimension(:),allocatable:: tab
  real(kind=pr):: alpha, sommetBnBn
  integer:: i, j
!=================================================================================================================================================================


  print*, "ce programme test des methodes de resolution de systeme lineaire de la forme: Ax=b"
  print*, "ou A est une matrice symetrique definie positive de taille n"
  print*, "les methodes calculent une solution approchee du systeme avec une certaine precision "

  
!=====================================================================================================================================================================
  !premiers tests sur une matrice petite et simple A symetrique et definie positive 
  !initialisation
  e=0.0001
  kmax=100000
  m=4
  n=10
  allocate(A(n,n),b(n),x0(n),x(n))

  !remplissage de la matrice
  A=1._pr
  do i=1,n
     A(i,i)=2._pr
  end do

  !vecteur second membre
  do i=1,n
     b(i)=i
  end do

  !vecteur initial de depart
  x0=0._pr
  
  !affichage des resultats
  print*, "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"
  print*, "test des methodes de resolution Ax=b  pour une matrice simple"
  print*, "taille du systeme:", n
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
  !x=FOM(A,b,x0,kmax,e,m)
 ! print*, x
  print*, "methode GMRes"
  x=GMRes(A,b,x0,kmax,e,m)
  print*, x
  print*, "methode GMRes amelioree pour une matrice sdp"
  x=GMRes_sdp(A,b,x0,kmax,e,m)
  print*,x
  
  deallocate(A,b,x0,x)
!====================================================================================================================================================
  !tests de FOM et GMRes classiques sur des matrices non sdp de petite taille
  print*, "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"
  print*, "test de FOM et GMRes sur une matrice quelconque non sdp de petite taille"

!====================================================================================================================================================
  !tests sur An=In-alpha*tBn* Bn 
  !Bn est remplie de coefficents aleatoires compris entre 0 et 1 et alpha coefficient suffisament grand pour que la matrice soit sdp

  
  !entree de l'utilisateur
  print*, "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"
  print*, "tests sur une matrice de la forme An=In-alpha*tBn* Bn "
  print*, "saisir n la taille du systeme:"
  read*, n
  print*, "saisir m pour definir l'espace de Krylov :"
  read*, m
  print*, "saisir la precision voulue sur la solution:"
  read*, e

  kmax=1000000
  allocate(An(n,n),Bn(n,n),Tn(n,n),x0(n),x(n),b(n),tab(n))

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
  
  !alpha=max sur i somme sur j(tBnBN)(i,j)
  !remplissage de tab
  do i=1,n
     sommetBnBn=0._pr
     do j=1,n
        sommetBnBn=sommetBnBn+Tn(i,j)
     end do
     tab(i)=sommetBnBN
  end do

  alpha=MAXVAL(tab)
    
  !remplissage de An à partir de Bn
  do i=1,n
     An(i,i)=1._pr
  end do

  An=alpha*Tn+An

  !vecteur second membre
  b=1._pr

  !vecteur initial de depart
  x0=0._pr

  !affichage des resultats
  print*, "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"
  print*, "methode du gradient a pas optimal"
  x=grad_pas_optimal(An,b,x0,kmax,e)
  print*, "methode du residu minimum"
  x=res_min(An,b,x0,kmax,e)
  print*, "methode du gradient conjugue"
  x=gradient_conjugue(An,b,x0,kmax,e)
 ! print*, x 
  print*, "methode FOM"
 ! x=FOM(An,b,x0,kmax,e,m)
  ! print*, x
  print*, "methode FOM amelioree pour une matrice sdp"
  !x=FOM_sdp(An,b,x0,kmax,e,m)
  !print*, x
  print*, "methode GMRes"
  x=GMRes(An,b,x0,kmax,e,m)
  print*, x
  print*, "methode GMRes amelioree pour une matrice sdp"
  x=GMRes_sdp(An,b,x0,kmax,e,m)
  print*,x
 
  deallocate(An,Bn,Tn,x0,x,b)
  
end program main
