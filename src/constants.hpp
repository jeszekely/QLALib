/**
 * \file "constants.hpp"
 * \author  J. Szekely
 */

#ifndef QLALIB_CONSTANTS
#define QLALIB_CONSTANTS

using namespace CONSTANTS {
	const double HBAR       = 1.0;
	const double EMASS      = 1.0;
	const double ECHARGE    = 1.0;
	const double VACPERM    = (1.0/(4.0*M_PI));
	const double LEN        = 0.0529177; // nm
	const double VEL        = 2.18e8; // cm/s
	const double EN         = 27.21; // eV
	const double TIME       = 2.42e-17; // s
	const double FREQ       = 4.13e16; // Hz
	const double ELECFIELD  = 5.14e9; // V/cm^2
	const double LASERINTEN = 3.51e16; // W/cm^2 (0.5*vacuum_permitivity*speed_of_light*EField^2)
	const double C          = 2.998e10/VEL; //speed of light
}

#endif