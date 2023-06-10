module stats_mod
use kind_mod, only: dp
implicit none
public :: mean, sd, bin_counts
contains
pure function mean(x) result(y)
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: y
y = sum(x)/max(1,size(x))
end function mean
!
pure function sd(x) result(y)
! sample standard deviation
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: y
real(kind=dp)             :: xmean
integer                   :: n
n = size(x)
if (n < 2) then
   y = -1.0_dp
   return
end if
xmean = sum(x)/n
y = sqrt(sum((x-xmean)**2)/(n-1))
end function sd
!
function bin_counts(x, thresh) result(counts)
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(in) :: thresh(:)
integer                   :: counts(size(thresh) + 1)
integer                   :: i, j, n, nthresh
n = size(x)
nthresh = size(thresh)
counts = 0
do i=1,n
   do_j: do j=1,nthresh
      if (x(i) < thresh(j)) exit do_j
   end do do_j
   counts(j) = counts(j) + 1
end do
end function bin_counts
end module stats_mod
