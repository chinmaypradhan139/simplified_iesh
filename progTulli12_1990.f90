program tulli12
use mod_global
implicit none
integer :: i,j,total_trajectories,number_of_cores,ntraj,TT
real*8 :: start,finish
real*8, dimension(:,:),allocatable :: population_matrix
real*8, dimension(:),allocatable :: final_pop

open(1, file ='input.txt')
read(1,*) number_of_cores,total_trajectories,iseed
close(1)

call setup_initial_values1
ntraj=int(total_trajectories/number_of_cores)


TT=int(total_time/dtc)

allocate(population_matrix(TT,ntraj))
allocate(final_pop(TT))

do i=1,ntraj
call setup_initial_values2


call CPU_TIME(start)
call classical_evolution 


population_matrix(:,i)=pop_mat


call CPU_time(finish)
write(119,*)finish-start
enddo



final_pop=0

do j=1,ntraj
   final_pop(:)=final_pop(:)+population_matrix(:,j)
enddo


do j=1,int(total_time/dtc)
   write(21,*) j*dtc,final_pop(j)/real(ntraj)
enddo
end program
!..............................................................................


