subroutine mysub(someargument,file,line)
  character(*) file
  integer someargument, line
  print *, 'mysub(',someargument,')'
  print *, 'I was called in file "',file,'" at line ',line
end subroutine mysub

#define mysub_wrapped(arg) mysub(arg,__FILE__,__LINE__)

program test
  integer x
  x = 12
  call mysub_wrapped(1)
end program test
