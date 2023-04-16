subroutine grid_setup

  use grid
  implicit none

  real :: S_init,T_init,depth_in,lat_in,lon_in
  call grid_read(S_init,T_init,depth_in,lat_in,lon_in)
  call grid_allocate
  call grid_init(S_init,T_init,depth_in,lat_in,lon_in)

return
end subroutine

