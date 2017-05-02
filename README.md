# HyperSIM
Nonlinear simulation of hypersonic flight dynamics, made with MATLAB and FlightGear

The MATLAB functions included are:

pvnrt.m: STP atmospheric properties, calculates air density
rayleigh_flow.m: Rayleigh flow inside the scramjet combustor, calculates heat added and exit Mach Number at the combustor
comb_inlet.m: Evaluates flow properties at the combustor inlet (Mach, Temperature, Pressure)
ctrlsurf.m: Control surface scramjet model (calculated aerodynamic lift/drag), Prandtl-Meyer Expansion and Oblique Shock computations
main_SINGLE_LOOP.m: ONE iteration of the main program
shockangle.m: Calculates oblique shock angle numerically (False Position method)
runfirst.m: Run this first to use in Simulink and FlightGear (SIMULINK AND FLIGHTGEAR ONLY)
oblique_down.m: Calculates oblique shock pressure and temperature.
exp_and_thrust.m: Calculates thrust of the scramjet
fnctns.m: EOM for a 3-DOF vehicle.
main.m: The main program (utilizes all subsystems above). Includes a 4th order RK solver for the 6 ODE'S

runfg.bat: FLIGHTGEAR ONLY, use to load FlightGear with MATLAB.

Run main.m first, it will output 6 graphs of the vehicle states (position, force etc.)

If you have SIMULINK and FlightGear, use the runfirst.m and the model.slx files to run your simulation. The CAD file is also included.
