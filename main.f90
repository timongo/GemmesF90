program main

  ! load module gemmes
  use gemmes_mod, only : init, read_namelist, initial_conditions, &
                         solve, output, endgemmes

  ! INITIALISATION
  ! Initialisation of variables to default
  call init
  ! Change default variables according to namelist
  call read_namelist
  ! Find all state initial conditions
  call initial_conditions

  ! BOUCLE TEMPORELLE
  ! Main program
  call solve
  ! Outputs
  call output

  ! END PROGRAM
  ! Free memory
  call endgemmes
  
end program main
