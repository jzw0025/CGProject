    subroutine distance(Pos, Tar, rc, rnum, Dim, Dlim, np, sumiden)
    implicit none
    integer,intent(in):: Dim, rnum
    real(8),intent(in):: rc
    real(8),intent(in), dimension(:,:) :: Pos
    real(8),intent(in), dimension(:,:) :: Tar
    real(8),intent(in), dimension(1:3,1:2):: Dlim

    real(8),intent(out), dimension(1:rnum,Dim):: np
    integer,intent(out):: sumiden
    real(8), dimension(Dim)::rij,ranp
    !defined interior variables
    real(8), dimension(:,:), allocatable :: Posa
    real(8), dimension(:,:), allocatable :: Tara
    real(8), dimension(:,:), allocatable :: TemS
    real(8), dimension(:,:), allocatable :: TemA
    real(8) :: rc2,d2,radius,angle1,angle2,u,r1,r2,r3
    real(8), parameter :: Pi = 3.1415927
    integer, dimension(1:rnum):: Identify
    integer :: j, ri

    rc2=rc*rc
    ! set up the size of array for the Samples and Active Points
    allocate(Posa(lbound(Pos,dim=1):ubound(Pos,dim=1),lbound(Pos,dim=2):ubound(Pos,dim=2)))
    allocate(Tara(lbound(Pos,dim=1):ubound(Pos,dim=1),lbound(Pos,dim=2):ubound(Pos,dim=2)))
    Posa = Pos
    Tara = Tar
    ! iterative arrays
    ! randomly choose point from active list
    !call random_number(u)
    !chotar = floor(NAcore*u)
    !Tari(1) = Tar(chotar,1)
    !Tari(2) = Tar(chotar,2)

    ! store Pos(i,:) in a temporary array for faster access in j loop
    ! initilize the random samples around the active point
    ! initilize the center for the chosen point
    ! initialize the array
    do ri = 1,rnum
        ! random_number function generating the uniform random number in (0,1)

        call random_number(r1)
        call random_number(r2)
	call random_number(r3)
        radius = rc*(r1+1)
        
        angle1 = 2*Pi*r2
	angle2 = 2*Pi*r3
        ! x coordinate and y coordinate
        ! old terms
        ranp(1) = Tar(1,1) +  rc*cos(angle1)*sin(angle2)
        ranp(2) = Tar(1,2) +  rc*sin(angle1)*sin(angle2)
	ranp(3) = Tar(1,3) +  rc*cos(angle2)
	
	! maybe the fixed terms(using radius instead of rc to get more tighten pack points)
	! ranp(1) = Tar(1,1) +  radius *cos(angle1)*sin(angle2)
        ! ranp(2) = Tar(1,2) +  radius *sin(angle1)*sin(angle2)
	! ranp(3) = Tar(1,3) +  radius *cos(angle2)

	! setup the simulation boundaries
	if (ranp(1) < Dlim(1,1)) then
            ranp(1) = Dlim(1,1)
        end if
	if (ranp(1) > Dlim(1,2)) then
            ranp(1) = Dlim(1,2)
        end if
	if (ranp(2) < Dlim(2,1)) then
            ranp(2) = Dlim(2,1)
        end if
	if (ranp(2) > Dlim(2,2)) then
            ranp(2) = Dlim(2,2)
        end if
	if (ranp(3) < Dlim(3,1)) then
            ranp(3) = Dlim(3,1)
        end if
	if (ranp(3) > Dlim(3,2)) then
            ranp(3) = Dlim(3,2)
        end if

        ! try to catch if any random point from existing sample points:
        ! fast starts from the newer sample point
        ! NAtom = size(Posa,1)
        do j = ubound(Posa,dim=1),lbound(Posa,dim=1),-1
            rij = Posa(j,:) - ranp
            !compute only the squared distance and compare to squared cut
            d2=sum(rij * rij)
            if (d2 < rc2) then
                identify(ri) = 1
                exit
            end if
            identify(ri) = 0
        enddo
        if (identify(ri) == 0) then
            ! adding to sample list
            allocate(TemS(lbound(Posa,dim=1):ubound(Posa,dim=1),lbound(Posa,dim=2):ubound(Posa,dim=2)))
            allocate(TemA(lbound(Tara,dim=1):ubound(Tara,dim=1),lbound(Tara,dim=2):ubound(Tara,dim=2)))
            TemS = Posa
            TemA = Tara
            Deallocate(Posa)
            Deallocate(Tara)
            Allocate(Posa(lbound(TemS,dim=1):(ubound(TemS,dim=1)+1),lbound(TemS,dim=2):ubound(TemS,dim=2)))
            Allocate(Tara(lbound(TemA,dim=1):(ubound(TemA,dim=1)+1),lbound(TemA,dim=2):ubound(TemA,dim=2)))
            Posa(lbound(TemS,dim=1):ubound(TemS,dim=1),lbound(TemS,dim=2):ubound(TemS,dim=2)) = TemS
            Tara(lbound(TemA,dim=1):ubound(TemA,dim=1),lbound(TemA,dim=2):ubound(TemA,dim=2)) = TemA
            ! adding the new point into the allocatable array
            Posa((ubound(TemS,dim=1)+1),lbound(TemS,dim=2))=ranp(1)
            Posa((ubound(TemS,dim=1)+1),lbound(TemS,dim=2)+1)=ranp(2)
			Posa((ubound(TemS,dim=1)+1),ubound(TemS,dim=2))=ranp(3)

            Tara((ubound(TemA,dim=1)+1),lbound(TemA,dim=2))=ranp(1)
            Tara((ubound(TemA,dim=1)+1),lbound(TemA,dim=2)+1)=ranp(2)
			Tara((ubound(TemA,dim=1)+1),ubound(TemA,dim=2))=ranp(3)

            Deallocate(TemS)
            Deallocate(TemA)
            ! update the sample array and activate array
        endif
    enddo
	! calculate the summation over the identify vector
	sumiden = rnum - sum(identify)
    !allocate(np(lbound(Posa,dim=1):ubound(Posa,dim=1),lbound(Posa,dim=2):ubound(Posa,dim=2)))
    np(lbound(Tara,dim=1):ubound(Tara,dim=1),lbound(Tara,dim=2):ubound(Tara,dim=2)) = Tara
    end subroutine

