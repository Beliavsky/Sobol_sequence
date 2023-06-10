program main
! compare Sobol and uniform variates in 1-D
use sobol_mod, only: i8_sobol
use stats_mod, only: mean, sd, bin_counts
implicit none
integer, parameter :: ikind=selected_int_kind(15), dp = kind(1.0d0), nbins=10
integer(kind=ikind), parameter :: dim_num = 1, nobs = 10**6, niter=5, &
                                  counts_exp = nobs/nbins
integer(kind=ikind) :: i, iter, seed, counts(nbins+1)
integer :: iran_type
logical, parameter :: print_x = .false.
real(kind=dp) :: xmat(dim_num, nobs), x(nobs), xthresh(nbins)
call random_seed()
print*,"#obs =",nobs
xthresh = [(i, i=1,nbins)]/real(nbins, kind=dp)
seed = 0
do iran_type=1,2
print "(/,a,a)", "random number type: ", merge("Sobol  ", "uniform", iran_type==1)
write (*,"(*(a10))") "mean", "sd", "dev_sq", "counts"
do iter=1, niter
   do i=1,nobs
      if (iran_type == 1) then
         call i8_sobol(dim_num, seed, xmat(:, i))
      else
         call random_number(xmat(:, i))
      end if
   end do
   x = xmat(1, :)
   if (print_x) print*,x
   counts = bin_counts(x, xthresh)
   print "(2f10.6, *(i10))", mean(x), sd(x), sum((counts(:nbins)-counts_exp)**2), counts(:nbins)
end do
end do
end program main
