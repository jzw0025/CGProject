    subroutine DCC(subset1, subset2, dsub1, dsub2, drij, rij )
    implicit none
    ! define external variables
    integer,intent(in):: drij
    ! dsub1 and dsub2 are the size of subset1 and subset2
    ! drij is the size of unbiased(unpadded) subset2 size
    real(8),intent(in), dimension(:, :, :) :: subset1
    real(8),intent(in), dimension(:, :, :) :: subset2
    real(8),intent(out), dimension(0:drij-1, 0:drij-1, 0:drij-1):: rij
    ! define internal variables
    integer :: x, y, z, padrad, dsub1, dsub2
    real(8) :: temp1_sqare_sum, temp2_sqare_sum, nomi,deno
    real(8), dimension(0:dsub1-1,0:dsub1-1,0:dsub1-1) :: temp1_mis
    real(8), dimension(0:dsub1-1,0:dsub1-1,0:dsub1-1) :: temp2
    real(8), dimension(0:dsub1-1,0:dsub1-1,0:dsub1-1) :: temp2_mis
    real(8), dimension(0:dsub1-1,0:dsub1-1,0:dsub1-1) :: A
    !dsub1 = size(subset1,1) ! get the size of subset1 in first dimension
    !dsub2 = size(subset2,1) ! get the size of subset2 in first dimension
    !loopping index
    padrad = (dsub1 - 1)/2
    temp1_mis = subset1 - sum(subset1)/dsub1/dsub1/dsub1
    !temp1_mis_pack = PACK(temp1_mis,.TRUE.)    
    temp1_sqare_sum = sum(temp1_mis**2)
    !$OMP PARALLEL DO
    do z = padrad,padrad+drij-1,1
        do y = padrad,padrad+drij-1,1
            do x = padrad,padrad+drij-1,1
                temp2(:,:,:) = subset2(x+1-padrad:x+1+padrad, y+1-padrad:y+1+padrad, z+1-padrad:z+1+padrad)
                temp2_mis = temp2-sum(temp2)/dsub1/dsub1/dsub1
                temp2_sqare_sum = sum(temp2_mis**2)
                A= temp1_mis*temp2_mis
                nomi = sum(A)
                ! adding the if condition to check the zeros of division
                if((temp1_sqare_sum==0.0).OR.(temp2_sqare_sum==0.0))then
                    rij(x-padrad,y-padrad,z-padrad) = 2.0  
                else
                    deno = sqrt(temp1_sqare_sum)*sqrt(temp2_sqare_sum)
                    !if(isnan(nomi/deno))then
                    !rij(x-padrad,y-padrad,z-padrad) = (x-padrad) 
                    !else
                    !1.0-nomi/deno
                    rij(x-padrad,y-padrad,z-padrad) = 1.0-nomi/deno
                    !end if               
                end if 
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO
    end subroutine


