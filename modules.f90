module state
  real(8) :: capital
  real(8) :: npop
  real(8) :: debt
  real(8) :: wage
  real(8) :: productivity 
  real(8) :: price
  real(8) :: eland
  real(8) :: sigma
  real(8) :: gsigma
  real(8) :: co2at
  real(8) :: co2up
  real(8) :: co2lo
  real(8) :: temp
  real(8) :: temp0
  real(8) :: pbs
  real(8) :: pcar
end module state

module other_vars
  real(8) :: n_red_fac ! reduction factor (n)
  real(8) :: abat ! abatement (A)
  real(8) :: dam ! damage
  real(8) :: dam_k ! capital damage
  real(8) :: dam_y ! output damage
  real(8) :: gdp ! GDP (Y)
  real(8) :: gdp0 ! GDP before abatement and damage (Y0)
  real(8) :: workforce ! Workforce (L)
  real(8) :: omega ! wage share wL/pY
  real(8) :: lambda ! employent rate L/N
  real(8) :: debtratio ! debta ratio d = D/(pY)
  real(8) :: eind ! industrial emissions
  real(8) :: inflation ! inflation
  real(8) :: rcb ! central bank interest rate
  real(8) :: pi ! profits
  real(8) :: smallpi ! profit ratio pi/(pY)
  real(8) :: smallpi_K ! return on assets
  real(8) :: leverage ! leverage
  real(8) :: cr
  real(8) :: forcing ! radiative forcing
  real(8) :: find ! industrial forcing
  real(8) :: fexo ! exogenous forcing
  real(8) :: emissions ! total emissions
  real(8) :: id ! demand for investment
  real(8) :: debtd ! demand for debt
  real(8) :: pir ! remaining profits after payment of dividends
  real(8) :: investment ! investment
end module other_vars

module model_pars
  real(8) :: alpha ! growth rate of productivity
  real(8) :: deltanpop ! leading growth rate of workforce
  real(8) :: npopbar ! maximum population in the logistic evolution
  real(8) :: nu ! Constant capital-to-output ratio
  real(8) :: pi1 ! damage function parameter
  real(8) :: pi2 ! damage function parameter
  real(8) :: pi3 ! damage function parameter
  real(8) :: zeta3 ! damage function parameter
  real(8) :: fdamk ! (paper = 1./3) fraction of environmental damage allocated to the stock of capital
  real(8) :: delta ! capital depreciation rate
  real(8) :: sa ! Fraction of abatement costs that are subsidized
  real(8) :: apc ! carbon price parameter
  real(8) :: bpc ! carbon price parameter
  real(8) :: conv10to15 ! conversion factor
  real(8) :: convco2toc=1./3.666 ! conversion factor
  real(8) :: deltagsigma ! dynamics of emissivity
  real(8) :: eta ! relaxation parameter of inflation 
  real(8) :: etar ! relaxation parameter of the interest rate
  real(8) :: mu ! markup of prices over the average cost
  real(8) :: omitted ! offset for the production cost in the inflation
  real(8) :: rstar ! Long-term interest rate target of the economy
  real(8) :: phitaylor ! parameter characterizing the reactivity of the monetary policy
  real(8) :: istar !interest rate targeted by the monetary policy
  real(8) :: srep ! Fraction of the outstanding debt repaid yearly
  real(8) :: climate_sens ! Climate sensitivity
  real(8) :: gammastar ! Heat exchange coefficient between temperature layers
  real(8) :: f2co2 ! Change in radiative forcing resulting from doubling of CO2
  real(8) :: cat_pind ! CO2 preindustrial concentration in atmosphere
  real(8) :: cup_pind ! CO2 preindustrial concentration in upper layer of ocean and biosphere
  real(8) :: clo_pind ! CO2 preindustrial concentration in bottom layer of the ocean
  real(8) :: fexo0 ! Initial value of the exogenous radiative forcing
  real(8) :: fexo1 ! value of the exogenous radiative forcing in 2100
  real(8) :: phi12 ! Transfer coefficient for carbon from the atmosphere to the upper ocean
  real(8) :: phi23 ! Transfer coefficient for carbon from the upper ocean/biosphere to the lower ocean
  real(8) :: deltapbs ! Exogenous growth rate of the back-stop technology price
  real(8) :: theta ! parameter of the abatement cost function
  real(8) :: phi0 ! Constant of the short-term Philips curve
  real(8) :: phi1 ! Slope of the short-term Philips curve
  real(8) :: kappa0 ! Constant of the investment function
  real(8) :: kappa1 ! Slope of the investment function
  real(8) :: kappamin ! Minimum of the investment function
  real(8) :: kappamax ! Maximum of the investment function
  real(8) :: div0 ! Constant of the dividend function
  real(8) :: div1 ! Slope of the dividend function
  real(8) :: divmin ! Minimum of the dividend function
  real(8) :: divmax ! Maximum of the dividend function
  real(8) :: heat_cap_at ! Heat capacity of the atmosphere, biosphere and upper ocean (C)
  real(8) :: heat_cap_lo ! Heat capacity of the deeper ocean (C0)
  real(8) :: delta_eland ! Growth rate of land-use change CO2 of emissions
  real(8) :: cr0 ! Constant of the leverage function
  real(8) :: crlev ! slope of the leverage function
  real(8) :: g0 ! growth rate of gdp0
  integer :: dam_type ! Type of damage (1 = 'No', 2 = 'Q', 3 = 'Extreme')
  integer :: rate_type ! Type of rate (1 = constant, 2 = Taylor)
  real(8), dimension(3,3) :: phimat
end module model_pars
        
module ini_cond
  real(8) :: co2at_ini
  real(8) :: co2up_ini
  real(8) :: co2lo_ini
  real(8) :: d_ini
  real(8) :: eind_ini
  real(8) :: eland_ini
  real(8) :: gsigma_ini
  real(8) :: pbs_ini
  real(8) :: n_ini
  real(8) :: npop_ini
  real(8) :: temp_ini
  real(8) :: temp0_ini
  real(8) :: y_ini
  real(8) :: lambda_ini
  real(8) :: omega_ini
  real(8) :: r_ini
  real(8) :: t_ini ! time
end module ini_cond

module num_pars
  real(8) :: dt
  real(8) :: tmax
end module num_pars

module solution
  integer :: nts
  real(8), allocatable, dimension(:,:) :: sol
  real(8), allocatable, dimension(:)   :: time
end module solution
