module debug

#define mysub_wrapped(arg) mysub(arg,__FILE__,__LINE__)

contains

subroutine mysub(someargument,file,line)
  character(*) file
  integer someargument, line
  print *, 'mysub(',someargument,')'
  print *, 'I was called in file "',file,'" at line ',line
end subroutine mysub

end module
