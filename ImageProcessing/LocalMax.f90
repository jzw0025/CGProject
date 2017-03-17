    subroutine LocalMaximum(subset1, radius, sx, sy, sz, rij)
    implicit none
    ! define external variables
    ! subset1 --- the large region
    ! radius
    
    integer,intent(in):: radius, sx, sy, sz
    
    real(8),intent(in), dimension(:, :, :) :: subset1
    
    real(8),intent(out), dimension(1:sx-2*radius, 1:sy-2*radius, 1:sz-2*radius):: rij
    
    ! define internal variables
    integer :: x, y, z
    real(8) :: localVal, globalVal, meanVal
    
    real(8), dimension(1:2*radius+1,1:2*radius+1,1:2*radius+1) :: cube

    !loopping index
    !$OMP PARALLEL DO
    do z = radius+1, sz-radius, 1
        do y = radius+1, sy-radius,1
            do x = radius+1, sx-radius,1
                
                cube(:,:,:) = subset1(x+1-radius:x+1+radius, y+1-radius:y+1+radius, z+1-radius:z+1+radius)
                
                localVal = subset1(x, y, z)
                globalVal = maxval(cube)
                meanVal = sum(cube)/(2*radius+1)/(2*radius+1)/(2*radius+1)
                
                ! adding the if condition 
                if ((localVal==globalVal) .and. (localVal  /= meanVal)) then
                    rij(x-radius,y-radius,z-radius) = localVal
                else 
                    rij(x-radius,y-radius,z-radius) = 0.0          
                end if 
                
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO
    end subroutine


