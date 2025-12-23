program demo_advance
implicit none

!THIS IS A DEMO PROGRAM TO PRINT CHARECTERS IN THE SAME line
!CAUTION : print*, don't suppress new line

integer :: i

print*, "NOW WE PRINT CHARECTERS IN THE SAME LINE"

write(*,'(A)',advance='no') '|'  !advance ='no' means dont advance to next line/ no new line

do I =1,5
    write(*,'(A)',advance='no') '0'
end do 

write(*,'(A)',advance='no') '>'
write(*,*)

print*,"================="
print*, "NOW WE PRINT CHARECTERS NOT IN THE SAME LINE"

do I =1,5
    write(*,'(A)') '0'
end do 



end program





