#ifndef COOLING_H
#define COOLING_H

#define N_METALLICITIES 8
#define N_TEMPS 91
#define MIN_TEMP 4.0 // log10(T/Kelvin)
#define MAX_TEMP 8.5 // log10(T/Kelvin)

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  double gas_cooling(struct galaxy_t* gal);
  void cool_gas_onto_galaxy(struct galaxy_t* gal, double cooling_mass);
  void read_cooling_functions(void);
  double interpolate_cooling_rate(double logTemp, double logZ);

#ifdef __cplusplus
}
#endif

#endif
