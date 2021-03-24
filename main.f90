program main

  ! Initialisation of variables to default
  call init
  ! Change default variables according to namelist
  call read_namelist
  ! Find all state initial conditions
  call initial_conditions
  ! Main program
  call solve
  ! Outputs
  call output
  ! Free memory
  call endgemmes
  
end program main

subroutine initial_conditions
  use state
  use other_vars
  use model_pars
  use ini_cond
  implicit none
  real(8) :: t2016, t2100
  real(8) :: tc
  real(8) :: transfers
  real(8) :: kappapi
  real(8) :: deltapik
  real(8) :: philam
  real(8) :: deltad
 

  ! time intermediate variables
  t2016 = 1.
  t2100 = 2100. - 2016.

  ! carbon price
  pcar = pbs_ini*n_ini**(theta - 1.)

  ! emissivity
  sigma = eind_ini/((1.-n_ini)*y_ini)

  ! reduction emission factor and abatement
  abat = 0.001*sigma*pbs_ini*n_ini**theta/theta

  select case(dam_type)
     case(1)
        dam = 0.
     case(2)
        dam = 1. - 1./(1 + pi1*temp_ini + pi2*temp_ini**2)
     case(3)
        dam = 1. - 1./(1 + pi1*temp_ini + pi2*temp_ini**2 + pi3*temp_ini**zeta3)
  end select
  
  ! Damage
  dam_k = fdamk*dam
  dam_y = 1. - (1.-dam)/(1.-dam_k)
  deltad = delta + dam_k

  ! Total cost of climate change
  tc = (1.-dam_y)*(1.-abat)
  ! GDP
  gdp0 = y_ini/tc
  capital = gdp0*nu
  ! workforce
  workforce = lambda_ini*npop_ini
  ! Wage share
  productivity = gdp0/workforce
  ! Initial price index
  price = 1.
  ! inflation
  inflation = eta*(mu*(omega_ini + omitted) - 1.)
  ! central bank interest rate
  call Taylor(inflation,rcb)

  ! Net transfers between public and private sectors
  transfers = sa*abat*gdp0 - pcar*conv10to15*eind_ini

  ! wages
  wage = omega_ini*price*y_ini/workforce
  ! Nominal debt
  debt = d_ini*price*y_ini

  ! profits
  pi = price*y_ini - wage*workforce - rcb*debt - deltad*price*capital + price*transfers
  smallpi_k = pi/(price*capital)
  smallpi = pi/(price*y_ini)

  ! emissions
  find = f2co2*log(co2at_ini/cat_pind)/log(2.)
  fexo = fexo0
  forcing = find + fexo
  emissions = eind_ini + eland_ini

  ! leverage
  leverage = d_ini*tc/nu
  call Tau(leverage,cr)

  ! growth rate of gpd0
  call kappa(smallpi,kappapi)
  call Dividends(smallpi_k,deltapik)
  g0 = ((1.-CR)*(1.25*kappapi*tc/nu - deltad) &
       & +CR*(smallpi_k - deltapik - tc*srep*d_ini/nu))

  ! initial carbon concentrations
  co2at = co2at_ini
  co2up = co2up_ini
  co2lo = co2lo_ini

  ! other variables that were not set yet
  npop = npop_ini
  eland = eland_ini
  gsigma = gsigma_ini
  temp = temp_ini
  temp0 = temp0_ini
  pbs = pbs_ini

end subroutine initial_conditions

subroutine solve
  use state
  use other_vars
  use model_pars
  use num_pars
  use ini_cond
  use solution
  implicit none
  integer :: i
  real(8) :: t
  real(8), dimension(16) :: u
  
  ! initialisation of solution
  time(1) = t_ini
  sol(1,1)    = capital
  sol(1,2)    = npop_ini
  sol(1,3)    = debt
  sol(1,4)    = wage
  sol(1,5)    = productivity
  sol(1,6)    = price
  sol(1,7)    = eland_ini
  sol(1,8)    = sigma
  sol(1,9)    = gsigma_ini
  sol(1,10)   = co2at_ini
  sol(1,11)   = co2up_ini
  sol(1,12)   = co2lo_ini
  sol(1,13)   = temp_ini
  sol(1,14)   = temp0_ini
  sol(1,15)   = pbs_ini
  sol(1,16)   = pcar
  sol(1,17)   = omega_ini
  sol(1,18)   = lambda_ini
  sol(1,19)   = d_ini
  sol(1,20)   = gdp0
  sol(1,21)   = y_ini
  sol(1,22)   = eind_ini
  sol(1,23)   = inflation
  sol(1,24)   = abat
  sol(1,25)   = n_ini
  sol(1,26)   = smallpi
  sol(1,27)   = smallpi_k
  sol(1,28)   = dam
  sol(1,29)   = dam_k
  sol(1,30)   = dam_y
  sol(1,31)   = fexo
  sol(1,32)   = find
  sol(1,33)   = rcb

  ! state init
  u(1)  = capital
  u(2)  = npop
  u(3)  = debt
  u(4)  = wage
  u(5)  = productivity
  u(6)  = price
  u(7)  = eland
  u(8)  = sigma
  u(9)  = gsigma
  u(10) = co2at
  u(11) = co2up
  u(12) = co2lo
  u(13) = temp
  u(14) = temp0
  u(15) = pbs
  u(16) = pcar

  do i=2,nts
     t = time(i-1)
     call rk4(t,u)
     time(i) = t_ini+dt*real(i-1,8)
     sol(i,1:16) = u
     sol(i,17)   = omega
     sol(i,18)   = lambda
     sol(i,19)   = debtratio
     sol(i,20)   = gdp0
     sol(i,21)   = gdp
     sol(i,22)   = eind
     sol(i,23)   = inflation
     sol(i,24)   = abat
     sol(i,25)   = n_red_fac
     sol(i,26)   = smallpi
     sol(i,27)   = smallpi_k
     sol(i,28)   = dam
     sol(i,29)   = dam_k
     sol(i,30)   = dam_y
     sol(i,31)   = fexo
     sol(i,32)   = find
     sol(i,33)   = rcb
  end do

end subroutine solve

subroutine system_function(t,u,k)
  use state
  use model_pars
  use other_vars
  implicit none

  real(8) :: t
  real(8), dimension(16), intent(in) :: u
  real(8), dimension(16), intent(out) :: k

  ! Intermediate variables
  real(8) :: t2016, t2100
  real(8) :: tc
  real(8) :: transfers
  real(8) :: kappapi
  real(8) :: deltapik
  real(8) :: philam
  real(8) :: deltad
  real(8) :: beta
  real(8) :: rho
  
  ! dstate/dt
  real(8) :: kdot
  real(8) :: ndot
  real(8) :: ddot
  real(8) :: wdot
  real(8) :: adot
  real(8) :: pdot
  real(8) :: elanddot
  real(8) :: sigmadot
  real(8) :: gsigmadot
  real(8) :: co2atdot
  real(8) :: co2updot
  real(8) :: co2lodot
  real(8) :: tdot
  real(8) :: t0dot
  real(8) :: pbsdot
  real(8) :: pcdot

  rho = f2co2/climate_sens

  ! expand state
  capital = u(1)
  npop = u(2)
  debt = u(3)
  wage = u(4)
  productivity = u(5)
  price = u(6)
  eland = u(7)
  sigma = u(8)
  gsigma = u(9)
  co2at = u(10)
  co2up = u(11)
  co2lo = u(12)
  temp = u(13)
  temp0 = u(14)
  pbs = u(15)
  pcar = u(16)

  ! time intermediate variables
  t2016 = 1.
  t2100 = 2100. - 2016.

  ! reduction emission factor and abatement
  n_red_fac = min((pcar/((1.-sa)*pbs))**(1./(theta-1.)),1.)
  abat = 0.001*sigma*pbs*n_red_fac**theta/theta

  select case(dam_type)
     case(1)
        dam = 0.
     case(2)
        dam = 1. - 1./(1 + pi1*temp + pi2*temp**2)
     case(3)
        dam = 1. - 1./(1 + pi1*temp + pi2*temp**2 + pi3*temp**zeta3)
  end select
  
  ! Damage
  dam_k = fdamk*dam
  dam_y = 1. - (1.-dam)/(1.-dam_k)
  deltad = delta + dam_k

  ! Total cost of climate change
  tc = (1.-dam_y)*(1.-abat)
  ! GDP
  gdp0 = capital/nu
  gdp = gdp0*tc
  ! workforce
  workforce = gdp0/productivity
  ! Wage share
  omega = wage*workforce/(price*gdp)
  ! employment parameter
  lambda = workforce/npop
  ! debt ratio
  debtratio = debt/(price*gdp)
  ! Industrial emission
  eind = gdp0*sigma*(1.-n_red_fac)

  ! inflation
  inflation = eta*(mu*(omega + omitted) - 1.)
  ! central bank interest rate
  call Taylor(inflation,rcb)

  ! Net transfers between public and private sectors
  transfers = sa*abat*gdp0 - pcar*conv10to15*eind
  ! profits
  pi = price*gdp - wage*workforce - rcb*debt - deltad*price*capital + price*transfers
  smallpi_k = pi/(price*capital)
  smallpi = pi/(price*gdp)
  ! leverage
  leverage = debtratio*tc/nu
  call Tau(leverage,cr)

  ! population growth
  beta = deltanpop*(1.-npop/npopbar)

  ! emissions
  find = f2co2*log(co2at/cat_pind)/log(2.)
  fexo = min(fexo0 + (fexo1 - fexo0)*t/84.,fexo1)
  forcing = find + fexo
  emissions = eind + eland
  
  ! investment
  call kappa(smallpi,kappapi)
  id = 1.25*kappapi*gdp
  call Dividends(smallpi_k,deltapik)
  pir = pi - deltapik*price*capital
  investment = cr*(pir/price + deltad*capital - srep*debt/price) &
       &     + (1. - cr)*id

  ! debt demand
  debtd = price*id - pir + srep*debt - deltad*price*capital

  ! computation of final results
  kdot = investment - deltad*capital
  ndot = beta*npop
  ddot = (1.-cr)*debtd - srep*debt
  call Phi(lambda,philam)
  wdot = philam*wage
  adot = productivity*alpha
  pdot = inflation*price
  elanddot = delta_eland*eland
  sigmadot = gsigma*sigma
  gsigmadot = deltagsigma*gsigma
  call co2dot(emissions,co2at,co2up,co2lo,co2atdot,co2updot,co2lodot)
  tdot = (forcing - rho*temp - gammastar*(temp-temp0))/heat_cap_at
  t0dot = gammastar*(temp-temp0)/heat_cap_lo
  pbsdot = pbs*deltapbs
  pcdot = pcar*(apc + bpc/(t+t2016))

  g0 = ((1.-CR)*(1.25*kappapi*tc/nu - deltad) &
       & +CR*(smallpi_k - deltapik - tc*srep*debtratio/nu))

  k(1)  = kdot
  k(2)  = ndot
  k(3)  = ddot
  k(4)  = wdot
  k(5)  = adot
  k(6)  = pdot
  k(7)  = elanddot
  k(8)  = sigmadot
  k(9)  = gsigmadot
  k(10) = co2atdot
  k(11) = co2updot
  k(12) = co2lodot
  k(13) = tdot
  k(14) = t0dot
  k(15) = pbsdot
  k(16) = pcdot

end subroutine system_function

subroutine Taylor(inflation,rcb)
  use model_pars
  implicit none
  
  real(8) :: inflation
  real(8) :: rcb
  
  rcb = min(max(rstar - phitaylor*istar + (1.+phitaylor)*inflation,0.),1.)

end subroutine Taylor

subroutine Tau(leverage,cr)
  use model_pars
  implicit none
  
  real(8) :: leverage
  real(8) :: cr

  cr = min(cr0 + crlev*leverage,1.)
  
end subroutine Tau

subroutine Kappa(pi,kappapi)
  use model_pars
  implicit none
  
  real(8) :: pi
  real(8) :: kappapi
  
  kappapi = min(max(kappa0 + kappa1*pi,0.),1.)

end subroutine Kappa

subroutine Dividends(pik,deltapik)
  use model_pars
  implicit none
  
  real(8) :: pik
  real(8) :: deltapik
  
  deltapik = min(max(div0 + div1*pik,divmin),divmax)

end subroutine Dividends

subroutine Phi(lambda,philam)
  use model_pars
  implicit none
  
  real(8) :: lambda
  real(8) :: philam
  
  philam = phi0 + phi1*lambda
  
end subroutine Phi

subroutine AX(n,A,X,Y)
  implicit none
  integer, intent(in) :: n
  real(8), dimension(n,n) :: A
  real(8), dimension(n) :: X,Y
  integer :: i,j

  Y = 0.

  do i=1,n
     do j=1,n
        Y(i) = Y(i) + A(i,j)*X(j)
     end do
  end do

end subroutine AX

subroutine co2dot(emissions,co2at,co2up,co2lo,co2atdot,co2updot,co2lodot)
  use model_pars
  implicit none
  
  real(8), intent(in) :: emissions
  real(8), intent(in) :: co2at,co2up,co2lo
  real(8), intent(out) :: co2atdot,co2updot,co2lodot
  
  real(8), dimension(3)   :: co2,co2dotvec

  co2(1) = co2at
  co2(2) = co2up
  co2(3) = co2lo

  call AX(3,phimat,co2,co2dotvec)

  co2atdot = emissions*convco2toc + co2dotvec(1)
  co2updot = co2dotvec(2)
  co2lodot = co2dotvec(3)

end subroutine co2dot

subroutine read_namelist
  use model_pars
  use ini_cond
  use num_pars
  use solution
  implicit none
  integer :: mp
  real(8) :: catup,cuplo

  namelist /model_parameters/ &
       alpha, &
       deltanpop, &
       npopbar, &
       nu, &
       pi1, &
       pi2, &
       pi3, &
       zeta3, &
       fdamk, &
       delta, &
       sa, &
       apc, &
       bpc, &
       conv10to15, &
       deltagsigma, &
       eta, &
       etar, &
       mu, &
       omitted, &
       rstar, &
       phitaylor, &
       istar, &
       srep, &
       climate_sens, &
       gammastar, &
       f2co2, &
       cat_pind, &
       cup_pind, &
       clo_pind, &
       fexo0, &
       fexo1, &
       phi12, &
       phi23, &
       deltapbs, &
       theta, &
       phi0, &
       phi1, &
       kappa0, &
       kappa1, &
       kappamin, &
       kappamax, &
       div0, &
       div1, &
       divmin, &
       divmax, &
       heat_cap_at, &
       heat_cap_lo, &
       delta_eland, &
       cr0, &
       crlev, &
       dam_type
  namelist /initial_conditions/ &
       co2at_ini, &
       co2up_ini, &
       co2lo_ini, &
       d_ini, &
       eind_ini, &
       eland_ini, &
       gsigma_ini, &
       pbs_ini, &
       n_ini, &
       npop_ini, &
       temp_ini, &
       temp0_ini, &
       y_ini, &
       lambda_ini, &
       omega_ini, &
       r_ini
  namelist /numerical_parameters/ &
       dt, &
       tmax
       
  mp=1
  open(mp, &
       & file='gemmes.dat', &
       & delim= 'apostrophe', &
       & form = 'formatted', &
       & action = 'read', &
       & status = 'old')
  read(mp,model_parameters)
  read(mp,initial_conditions)
  read(mp,numerical_parameters)
  close(mp)

  catup = cat_pind/cup_pind
  cuplo = cup_pind/clo_pind

  phimat(1,1) = -phi12
  phimat(1,2) = phi12*catup
  phimat(2,1) = phi12
  phimat(2,2) = -phi12*catup - phi23
  phimat(2,3) = phi23*cuplo
  phimat(3,2) = phi23
  phimat(3,3) = -phi23*cuplo  

  nts = floor(tmax/dt)+1
  allocate(sol(nts,33),time(nts))

end subroutine read_namelist

subroutine output
  use solution
  implicit none
  integer :: i,j
  integer :: mp

  mp = 1
  open(mp, file='gemmes.out', form = 'formatted')
  write(mp,'(34A)') 'time ', 'capital  ','npop  ','debt  ','wage  ','productivity  ','price  ','eland  ','sigma  ','gsigma  ','co2at  ','co2up  ','co2lo  ','temp  ','temp0  ','pbs  ','pcar  ','omega  ','lambda  ','debtratio  ','gdp0  ','gdp  ','eind  ','inflation  ','abat  ','n_red_fac  ','smallpi  ','smallpi_k  ','dam  ','dam_k  ','dam_y  ','fexo  ','find  ','rcb  '
  do i=1,nts
     write(mp,'(34E24.16)') time(i), (sol(i,j), j=1,33)
  end do
  close(mp)

end subroutine output

subroutine endgemmes
  use solution
  implicit none

  deallocate(sol,time)

end subroutine endgemmes

subroutine rk4(t,u)
  use num_pars, only : dt
  implicit none
  real(8), intent(in) :: t
  real(8) :: t_aux
  real(8), dimension(16), intent(inout) :: u
  real(8), dimension(16) :: u_aux,k1,k2,k3,k4

  call system_function(t,u,k1)
  u_aux = u + k1*0.5*dt
  t_aux = t + 0.5*dt
  
  call system_function(t_aux,u_aux,k2)
  u_aux = u + k2*0.5*dt
  
  call system_function(t_aux,u_aux,k3)
  u_aux = u + k3*dt
  t_aux = t + dt
  
  call system_function(t_aux,u_aux,k4)

  u = u + (k1+2.*(k2+k3)+k4)*dt/6.

end subroutine rk4

subroutine init
  use model_pars
  use ini_cond
  use num_pars
  implicit none
  ! Default variables
  alpha=0.02 ! growth rate of productivity
  deltanpop=0.0305 ! leading growth rate of workforce
  npopbar=7.055925493 ! maximum population in the logistic evolution
  nu=2.7 ! Constant capital-to-output ratio
  pi1=0. ! damage function parameter
  pi2=0.00236 ! damage function parameter
  pi3=0.0000819 ! damage function parameter
  zeta3=6.754 ! damage function parameter
  fdamk=0. ! (paper = 1./3) fraction of environmental damage allocated to the stock of capital
  delta=0.04 ! capital depreciation rate
  sa=0 ! Fraction of abatement costs that are subsidized
  apc=0.02 ! carbon price parameter
  bpc=0. ! carbon price parameter
  conv10to15=1.160723971/1000. ! conversion factor
  deltagsigma=-0.001 ! dynamics of emissivity
  eta=0.5 ! relaxation parameter of inflation 
  etar=10. ! relaxation parameter of the interest rate
  mu=1.3 ! markup of prices over the average cost
  omitted=0.3 ! offset for the production cost in the inflation
  rstar=0.01 ! Long-term interest rate target of the economy
  phitaylor=0.5 ! parameter characterizing the reactivity of the monetary policy
  istar=0.02 !interest rate targeted by the monetary policy
  srep=0.1 ! Fraction of the outstanding debt repaid yearly
  climate_sens=3.1 ! Climate sensitivity
  gammastar=0.0176 ! Heat exchange coefficient between temperature layers
  f2co2=3.6813 ! Change in radiative forcing resulting from doubling of CO2
  cat_pind=588. ! CO2 preindustrial concentration in atmosphere
  cup_pind=360. ! CO2 preindustrial concentration in upper layer of ocean and biosphere
  clo_pind=1720. ! CO2 preindustrial concentration in bottom layer of the ocean
  fexo0=0.5 ! Initial value of the exogenous radiative forcing
  fexo1=1. ! value of the exogenous radiative forcing in 2100
  phi12=0.0239069 ! Transfer coefficient for carbon from the atmosphere to the upper ocean
  phi23=0.0013409 ! Transfer coefficient for carbon from the upper ocean/biosphere to the lower ocean
  deltapbs=-.00505076337 ! Exogenous growth rate of the back-stop technology price
  theta=2.6 ! parameter of the abatement cost function
  phi0=-.291535421 ! Constant of the short-term Philips curve
  phi1=.468777035 ! Slope of the short-term Philips curve
  kappa0=.031774 ! Constant of the investment function
  kappa1=.575318414 ! Slope of the investment function
  kappamin=0. ! Minimum of the investment function
  kappamax=0.3 ! Maximum of the investment function
  div0=0.0512863 ! Constant of the dividend function
  div1=0.4729 ! Slope of the dividend function
  divmin=0. ! Minimum of the dividend function
  divmax=1.! Maximum of the dividend function
  heat_cap_at=49.75817526819656 ! Heat capacity of the atmosphere biosphere and upper ocean
  heat_cap_lo=3.52 ! Heat capacity of the deeper ocean
  delta_eland=-.0220096 ! Growth rate of land-use change CO2 of emissions
  cr0=0.17 ! Constant of the leverage function
  crlev=.1! slope of the leverage function
  dam_type = 1 ! by default, no damage
  
  ! Initial conditions
  co2at_ini=851.
  co2up_ini=460.
  co2lo_ini=1740.
  d_ini=1.53282171
  eind_ini=35.85
  eland_ini=2.6
  gsigma_ini=-0.0105
  pbs_ini=547.2220801465
  n_ini=0.03
  npop_ini=4.825484061
  temp_ini=0.85
  temp0_ini=0.0068
  y_ini=59.7387
  lambda_ini=0.674828
  omega_ini=0.5782323
  r_ini=0.10627649250000004

  ! initial time
  t_ini = 0.

  ! Numerical parameters
  dt = 0.05
  tmax = 84.
end subroutine init


