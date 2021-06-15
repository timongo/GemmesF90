!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of the iLOVECLIM coupled climate model under the LUDUS framework.
!!      global_constants_mod is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
!!      License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
!!      version.
!!
!!      global_constants_mod is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
!!      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: gemmes_mod
!
!>     @author  Hugo Martin, Aurélien Quiquet, Didier M. Roche (dmr)
!
!
!>     @brief This module gemmes_mod is used to couple iLOVECLIM with the GEMMES Integrated Assessment Model
!
!>     @date Creation date: December, 4th, 2019
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : ham
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!     Now this module contains all the gemmes code.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! iLVC = 0 : Run Gemmes as an independant program.
!      = 1 : Run the iLOVECLIM/GEMMES model full fortran.
!      = 2 : Run the iLOVECLIM/GEMMES model fortran/R. 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#define iLVC 2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   stdin, stdout definition ... global variables and subroutines
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module gemmes_mod

#if ( iLVC == 0 )
        implicit none
        
        public :: init, read_namelist, initial_conditions, &
                  solve, output, endgemmes

        private
#else
        use global_constants_mod, only: dblp=>dp

        implicit none
        
        public :: gemmes_init, gemmes_step, gemmes_accum_tglob, &
                  gemmes_recup_emissions, gemmes_emissions, timelov, &
                  init, read_namelist, initial_conditions, &
                  solve, output, endgemmes

        private
        
        intrinsic :: get_environment_variable,execute_command_line
        
        real(kind=dblp),save :: timelov
        real(kind=dblp),save :: glob_t2m_init
        real(kind=dblp),save :: glob_t2m_accum
        real(kind=dblp),save :: glob_t2m_anom
        real(kind=dblp),save :: t2m_to_gemmes
        real(kind=dblp),save :: gemmes_emissions
        character(len=128), parameter :: gemmespath = "/home/giraud/&
             &Code_Hugo/iloveclim_gemmes/gemmes/"
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Declaration of variables and constants.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! module state
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
        !end module state
       
        !module other_vars
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
        !end module other_vars
        
        !module model_pars
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
        integer :: infla_type ! Type of inflation 1 constant, 2 with omitted
        real(8) :: infla ! conatnt value of inflation
        real(8), dimension(3,3) :: phimat
        !end module model_pars
               
        !module ini_cond
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
        !end module ini_cond
       
        !module num_pars
        real(8) :: dt
        real(8) :: tmax
        !end module num_pars
       
        !module solution
        integer :: nts
        real(8), allocatable, dimension(:,:) :: sol
        real(8), allocatable, dimension(:)   :: time
        !end module solution       

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Declaration of subroutines.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! subroutine init
          subroutine init
       
            implicit none
          
            ! Default variables
          
            ! Climate damages
            pi1=0. !* damage function parameter
            pi2=0.00236 !* damage function parameter
            pi3=0.0000819 !* damage function parameter
            zeta3=6.754 !* damage function parameter
            fdamk=1./3 !* (paper = 1./3) fraction of environmental damage allocated
                       !  to the stock of capital
            dam_type = 1 !  by default, 1: no damage
                         !              2: medium damage
                         !*             3: high damage
            ! Workforce
            deltanpop=0.0305 !* leading growth rate of workforce
            npopbar=7.055925493 !* maximum population in the logistic evolution
          
            ! Global population
            !H il manque ces variables
          
            ! Capital
            delta=0.04 !* capital depreciation rate
            nu=2.7 !* Constant capital-to-output ratio
          
            ! Dividends
            div0=0.0512863 ! Constant of the dividend function
            div1=0.4729 !* Slope of the dividend function
            divmin=0. !* Minimum of the dividend function
            divmax=.3 !* Maximum of the dividend function; last value 1.
            !H Je ne retrouve pas div0 de la meme maniere
          
            ! Investment
            kappa0=.031774 !* Constant of the investment function
            kappa1=.575318414 !* Slope of the investment function
            kappamax=0.3 !* Maximum of the investment function
            kappamin=0. !* Minimum of the investment function
            !H il manque I_Dam mais il est a 0 dans le code R
          
            ! Inflation
            infla_type = 1 !1 constant, 2 with omitted, 3 with WACC
            infla = 0.02 ! constant inflation value
            eta=0.5 !* relaxation parameter of inflation 
            mu=1.3 !* markup of prices over the average cost
            omitted=0.3 ! offset for the production cost in the inflation
          
            ! Productivity
            alpha=0.02 !* growth rate of productivity
            !H Il me manque a_Tp_1, a_Tp_2, min max
            !H On est dans le cas constant du code R
          
            ! Wages
            phi0=-.291535421 !* Constant of the short-term Philips curve
            phi1=.468777035 !* Slope of the short-term Philips curve
            !H  il me manque m mais j'ai l'impression qu'il ne sert à rien dans le code R
          
            ! Interest rate
            phitaylor=0.5 ! param characterizing the reactivity of the monetary policy
            etar=10. ! relaxation parameter of the interest rate
            rstar=0.01 ! Long-term interest rate target of the economy
                       ! This is also the rate value if constant
            istar=0.02 ! interest rate targeted by the monetary policy
            srep=0.1 ! Fraction of the outstanding debt repaid yearly
          
            ! CO2 emissions
            delta_eland=-.0220096 !* Growth rate of land-use change CO2 of emissions
            deltagsigma=-0.001 !* dynamics of emissivity
          
            ! Abatement and control prices
            theta=2.6 !* parameter of the abatement cost function
            deltapbs=-.00505076337 !* Exo growth rate of the back-stop technology price
          
            ! Carbon price
            conv10to15=1.160723971/1000. !* conversion factor
            apc=0.125 !* carbon price parameter
            bpc=0.625 !* carbon price parameter
          
            ! public - private sector
            sa=0 ! Fraction of abatement costs that are subsidized
          
            ! leverage
            !H pas de leverage dans le code R
            cr0=0. ! Constant of the leverage function
            crlev=0.! slope of the leverage function
          
            ! Climate constants for climate module
            climate_sens=3.1 ! Climate sensitivity
            gammastar=0.0176 ! Heat exchange coefficient between temperature layers
            f2co2=3.6813 ! Change in radiative forcing resulting from doubling of CO2
            cat_pind=588. ! CO2 preind conc in atmosphere
            cup_pind=360. ! CO2 preind conc in upper layer of ocean and biosphere
            clo_pind=1720. ! CO2 preindustrial concentration in bottom layer of the ocean
            fexo0=0.5 ! Initial value of the exogenous radiative forcing
            fexo1=1. ! value of the exogenous radiative forcing in 2100
            phi12=0.0239069 ! Transfer coef for carbon from the atmo to the upper ocean
            phi23=0.0013409 ! Transfer coef for carbon from the upper ocean/biosphere 
                            ! to the lower ocean
            heat_cap_at=49.75817526819656 ! Heat capacity of the atmo biosphere and 
                                          ! upper ocean
            heat_cap_lo=3.52 ! Heat capacity of the deeper ocean
          
            ! Initial conditions, values with * are equal as in the code R
            !H par rapport au code R, il me manque NG, F_exo ; r_ini dans R.
            y_ini=59.7387 !*
            npop_ini=4.825484061 !*
            lambda_ini=0.674828 !*
            omega_ini=0.5782323 !*
            d_ini=1.53282171 !*
            eland_ini=2.6 !*
            eind_ini=35.85 !*
            n_ini=0.03 !*
            gsigma_ini=-0.0152 !* replace last value =-0.0105
            pbs_ini=547.2220801465 !*
          
            co2at_ini=851. !*
            co2up_ini=460. !*
            co2lo_ini=1740. !*
            temp_ini=0.85 !*
            temp0_ini=0.0068 !*
          
            r_ini=0.10627649250000004 !H a quoi sert r_ini
          
            ! initial time
            t_ini = 0.
          
            ! Numerical parameters
            dt = 0.05
            tmax = 84.
          end subroutine init

! subroutine read_namelist
          subroutine read_namelist

            implicit none

            integer :: mp
            real(8) :: catup,cuplo
          
            !namelist /model_parameters/ &
            namelist /gemmes_mod/ &
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
                 dam_type, &
                 infla_type, &
                 infla, &

            !namelist /initial_conditions/ &
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
                 r_ini, &

            !namelist /numerical_parameters/ &
                 dt, &
                 tmax
                 
            mp=1
            open(mp, &
                 & file='gemmes.dat', &
                 & delim= 'apostrophe', &
                 & form = 'formatted', &
                 & action = 'read', &
                 & status = 'old')
            read(mp,gemmes_mod)
            !read(mp,initial_conditions)
            !read(mp,numerical_parameters)
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

! subroutine initial_conditions
          subroutine initial_conditions

            implicit none

            real(8) :: t2016, t2100
            real(8) :: tc
            real(8) :: transfers
            real(8) :: kappapi
            real(8) :: deltapik
            real(8) :: philam
            real(8) :: deltad
            real(8) :: costprod
            real(8) :: WACC
          
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
                  dam = 1.-1./(1 + pi1*temp_ini + pi2*temp_ini**2 & 
                        & + pi3*temp_ini**zeta3)
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
          
            ! wages
            wage = omega_ini*price*y_ini/workforce
          
            ! Net transfers between public and private sectors
            transfers = sa*abat*gdp0 - pcar*conv10to15*eind_ini
          
            ! Nominal debt
            debt = d_ini*price*y_ini
          
            ! inflation for infla_type=1, 2 and rate
            if (infla_type<3) then
                if (infla_type.eq.1) then
                    inflation = infla
                else
                    inflation = eta*(mu*(omega_ini + omitted) - 1.)
                endif
               ! central bank interest rate
               call Taylor(inflation,rcb)
            else
               rcb =  rstar
            endif
          
            ! profits
            pi = price*y_ini - wage*workforce - rcb*debt &
                 & - deltad*price*capital + price*transfers
            smallpi_k = pi/(price*capital)
            smallpi = pi/(price*y_ini)
          
            ! Dividends
            call Dividends(smallpi_k,deltapik)
          
            ! inflation for infla_type = 3
            if (infla_type.eq.3) then
               WACC = (rcb*debt + deltapik*price*capital)/price/capital
               costprod = 1./tc*(wage/price/productivity &
                              & + pcar*conv10to15*sigma*(1-n_ini) &
                              & + nu*(WACC + deltad))
               inflation = eta*(mu*costprod - 1.)
            endif
          
            ! leverage
            leverage = debtratio*tc/nu
            call Tau(leverage,cr)
          
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
          
            g0 = ((1.-cr)*(1.25*kappapi*tc/nu - deltad) &
                 & +cr*(smallpi_k - deltapik - tc*srep*d_ini/nu))
          
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

! subroutine solve
          subroutine solve

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

! subroutine rk4
          subroutine rk4(t,u)

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
          
            !H carbon price if like in R code
            if (u(16).gt.u(15)) then
                u(16) = u(15)
            end if
          
          end subroutine rk4

! subroutine system_function
          subroutine system_function(t,u,k)

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
            real(8) :: costprod
            real(8) :: WACC
          
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
          
            !H if like in R code
            if (pcar.gt.pbs) then
                pcar = pbs
            end if
          
            ! time intermediate variables
            t2016 = 1.
            t2100 = 2100. - 2016.
          
            ! Say's law
            call Sayslaw(npop,capital,productivity,gdp0,workforce,lambda)
          
            ! reduction emission factor and abatement
            n_red_fac = min((pcar/((1.-sa)*pbs))**(1./(theta-1.)),1.)
            abat = 0.001*sigma*pbs*n_red_fac**theta/theta
          
            select case(dam_type)
               case(1)
                  dam = 0.
               case(2)
                  dam = 1. - 1./(1 + pi1*temp + pi2*temp**2)
               case(3)
                  dam = 1. - 1./(1 + pi1*temp + pi2*temp**2 &
                        & + pi3*temp**zeta3)
            end select
            
            ! Damage
            dam_k = fdamk*dam
            dam_y = 1. - (1.-dam)/(1.-dam_k)
            deltad = delta + dam_k
          
            ! Total cost of climate change
            tc = (1.-dam_y)*(1.-abat)
          
            ! GDP
            gdp = gdp0*tc
          
            ! Wage share
            omega = wage*workforce/(price*gdp)
          
            ! debt ratio
            debtratio = debt/(price*gdp)
          
            ! Industrial emission
            eind = gdp0*sigma*(1.-n_red_fac)
          
            ! inflation and rate for infla_type=1, 2
            if (infla_type<3) then
                if (infla_type.eq.1) then
                    inflation = infla
                else
                    inflation = eta*(mu*(omega + omitted) - 1.)
                endif
               ! central bank interest rate
               call Taylor(inflation,rcb)
            else
               rcb =  rstar
            endif
          
            ! Net transfers between public and private sectors
            transfers = sa*abat*gdp0 - pcar*conv10to15*eind
          
            ! profits
            pi = price*gdp - wage*workforce - rcb*debt & 
                 & - deltad*price*capital + price*transfers
            smallpi_k = pi/(price*capital)
            smallpi = pi/(price*gdp)
          
            ! Dividends
            call Dividends(smallpi_k,deltapik)
          
            ! inflation for infla_type=3
            if (infla_type.eq.3) then
               WACC = (rcb*debt + deltapik*price*capital)/price/capital
               costprod = 1./tc*(wage/price/productivity & 
                              & + pcar*conv10to15*sigma*(1-n_red_fac) &
                              & + nu*(WACC + deltad))
               inflation = eta*(mu*costprod - 1.)
            endif
          
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
            id = kappapi*gdp
            pir = pi - deltapik*price*capital
            investment = cr*(pir/price + deltad*capital & 
                         & - srep*debt/price) + (1. - cr)*id
          
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
            call co2dot(emissions,co2at,co2up, &
                        & co2lo,co2atdot,co2updot,co2lodot)
            tdot = (forcing - rho*temp & 
                   & - gammastar*(temp-temp0))/heat_cap_at
            t0dot = gammastar*(temp-temp0)/heat_cap_lo
            pbsdot = pbs*deltapbs
            pcdot = pcar*(apc + bpc/(t+t2016))
          
            g0 = ((1.-cr)*(1.25*kappapi*tc/nu - deltad) &
                 & +cr*(smallpi_k - deltapik - tc*srep*debtratio/nu))
          
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

! subroutine Taylor
          subroutine Taylor(inflation,rcb)

            implicit none
            
            real(8) :: inflation
            real(8) :: rcb
            
            rcb = min(max(rstar - phitaylor*istar &
                  & + (1.+phitaylor)*inflation,0.),1.)
          
          end subroutine Taylor

! subroutine Sayslaw
          subroutine Sayslaw(npop,capital,productivity,gdp0, &
                             & workforce,lambda)

            implicit none
          
            real(8) :: npop
            real(8) :: capital
            real(8) :: productivity
            real(8) :: gdp0
            real(8) :: workforce
            real(8) :: lambda
          
            ! GDP
            gdp0 = capital/nu
          
            ! workforce
            workforce = gdp0/productivity
          
            ! check Say's law and set employment parameter
            if (workforce>npop) then
                lambda = 1
                workforce = npop
                gdp0 = productivity*workforce
                capital = nu*gdp0
            else
                lambda = workforce/npop
            endif
          
          end subroutine Sayslaw
          
! subroutine Tau
          subroutine Tau(leverage,cr)

            implicit none
            
            real(8) :: leverage
            real(8) :: cr
          
            cr = min(cr0 + crlev*leverage,1.)
            
          end subroutine Tau
          
! subroutine Kappa
          subroutine Kappa(pi,kappapi)

            implicit none
            
            real(8) :: pi
            real(8) :: kappapi
            
            kappapi = min(max(kappa0 + kappa1*pi,kappamin),kappamax)
          
          end subroutine Kappa
          
! subroutine Dividends
          subroutine Dividends(pik,deltapik)

            implicit none
            
            real(8) :: pik
            real(8) :: deltapik
            
            deltapik = min(max(div0 + div1*pik,divmin),divmax)
           
          end subroutine Dividends
          
! subroutine Phi
          subroutine Phi(lambda,philam)

            implicit none
            
            real(8) :: lambda
            real(8) :: philam
            
            philam = phi0 + phi1*lambda
            
          end subroutine Phi
          
! subroutine AX
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
          
! subroutine co2dot
          subroutine co2dot(emissions,co2at,co2up,co2lo, &
                            & co2atdot,co2updot,co2lodot)

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

! subroutine output
          subroutine output

            implicit none
            integer :: i,j
            integer :: mp
          
            mp = 1
            open(mp, file='gemmes.out', form = 'formatted')
            write(mp,'(34A)') 'time ', &
                 &            'capital  ', &
                 &            'npop  ', &
                 &            'debt  ', &
                 &            'wage  ', &
                 &            'productivity  ', &
                 &            'price  ', &
                 &            'eland  ', &
                 &            'sigma  ', &
                 &            'gsigma  ', &
                 &            'co2at  ', &
                 &            'co2up  ', &
                 &            'co2lo  ', &
                 &            'temp  ', &
                 &            'temp0  ', &
                 &            'pbs  ', &
                 &            'pcar  ', &
                 &            'omega  ', &
                 &            'lambda  ', &
                 &            'debtratio  ', &
                 &            'gdp0  ', &
                 &            'gdp  ', &
                 &            'eind  ', &
                 &            'inflation  ', &
                 &            'abat  ', &
                 &            'n_red_fac  ', &
                 &            'smallpi  ', &
                 &            'smallpi_k  ', &
                 &            'dam  ', &
                 &            'dam_k  ', &
                 &            'dam_y  ', &
                 &            'fexo  ', &
                 &            'find  ', &
                 &            'rcb  '
            do i=1,nts
               write(mp,'(34E24.16)') time(i), (sol(i,j), j=1,33)
            end do
            close(mp)
          
          end subroutine output
          
! subroutine endgemmes          
          subroutine endgemmes

            implicit none
          
            deallocate(sol,time)
          
          end subroutine endgemmes

#if ( iLVC == 0 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here if independant GEMMES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module gemmes_mod

#else

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! subroutines to interact with iLOVECLIM
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! subroutine gemmes_init
        subroutine gemmes_init

          character(len=256) :: runinfo ! to get the local directory
          
          timelov = 2016. ! start year of climate-economy simulations (GEMMES reference)
          
          call get_environment_variable('PWD',runinfo)
          write(*,*) runinfo
          open(520,file=trim(gemmespath)//"runinfo.txt")
          write(520,*) trim(runinfo)
          close(520)
          open(520,file='timelovinit.txt')
          write(520,*) timelov
          close(520)
          open(520,file='timelov.txt')
          write(520,*) timelov
          close(520)
          ! glob_t2m_accum is the mean global temperature after 1 year of simulated time in iLoveclim
          glob_t2m_init = glob_t2m_accum
          glob_t2m_anom = 0. 
          t2m_to_gemmes = 0.85 + glob_t2m_anom
          ! t2m_to_GEMMES: 2016 t2m anom relative to preind (0.85: GEMMES default)
          open(520,file='t2mlov.txt')
          write(520,*) t2m_to_gemmes
          close(520)
          call execute_command_line(                                   &
         &  trim(gemmespath)//"gemmes_init.sh",wait=.true.)
          write(*,*) "on a appele l'init", t2m_to_gemmes, glob_t2m_init

          gemmes_emissions = 0.
          
        end subroutine gemmes_init
        
! subroutine gemmes_step
        subroutine gemmes_step
           
          timelov = timelov+1 ! yearly coupling
          
          open(520,file='timelov.txt')
          write(520,*) timelov
          close(520)
          glob_t2m_anom = glob_t2m_accum - glob_t2m_init 
          t2m_to_gemmes = 0.85 + glob_t2m_anom
          ! t2m_to_GEMMES: 2016 t2m anom relative to preind (0.85: GEMMES default)
          open(520,file='t2mlov.txt')
          write(520,*) t2m_to_gemmes
          close(520)
          call execute_command_line(                                   &
         &  trim(gemmespath)//"gemmes_step.sh",wait=.true.)
          write(*,*) "on a appele le step", t2m_to_gemmes,             &
          glob_t2m_accum

        end subroutine gemmes_step
         
! subroutine gemmes_accum_tglob
        subroutine gemmes_accum_tglob(eoy)
           
          use ipcc_output_mod, only: tsurfmean

          implicit none

          logical,intent(in) :: eoy
          !real(8) :: tsurfmean

          if (eoy) then !end of year, we reset the global mean
             glob_t2m_accum = 0.
          else          !accumulation frequency: day
             glob_t2m_accum = glob_t2m_accum + tsurfmean/360.
          endif

        end subroutine gemmes_accum_tglob

! subroutine gemmes_recup_emissions
        subroutine gemmes_recup_emissions

          implicit none

          real(kind=dblp) :: gemmes_emissions_yearly
          !real(8) :: gemmes_emissions_yearly          

          open(520,file=trim(gemmespath)//"emissions.txt")
          read(520,*) gemmes_emissions_yearly
          close(520)
          gemmes_emissions = gemmes_emissions + gemmes_emissions_yearly
          write(*,*) "GEMMES emissions:", gemmes_emissions, & 
          gemmes_emissions_yearly
          
        end subroutine gemmes_recup_emissions

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here for iLOVECLIM/GEMMES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module gemmes_mod

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
