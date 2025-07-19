# solveForDuctDiamterGivenPressureLoss

This function will solve for the hydraulic diameter of a duct such that a desired pressure loss across the duct. As inputs this function requires

 - deltaPTarget: the desired pressure loss across the duct (Pa)
 - mdot: mass of air moving per second (kg/s)
 - L : the length of the pipe (m)
 - rho : the density of the air (kg/m^3)
 - nu : kinematic viscosity (m^2/s)
 - k : surface roughness (m)

This function numerically solves for d (hydraulic diameter) of the duct using the following steps

Given a candidate d

1. Directly evaluate
    - A: area
    - v: velocity of the air
    - Re: the Reynolds number
2. Numerically solve the Colebrook equation for the friction factor f using a (modified) Newton's method
3. Directly solve for the pressure loss

The final value for d is solved for using the above procedure and the bisection method.


## Usage

1. Open Matlab/Octave
2. Open main_solveForHydraulicDiameter.m
3. Set usingOctave to 1 if you're using Octave and to 0 for Matlab.
4. Run the script
5. The file ductSolutionLog.txt will appear in the output directory
6. If flag_generatePlot = 1, then a pdf plot will also appear in the output folder

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)
