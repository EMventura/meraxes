#ifndef STAR_FORMATION_H
#define STAR_FORMATION_H

#include "meraxes.h"

struct FR_parameters
{
  double a;
  double b;
  double c;
  double d;
};

typedef enum SFtype
{
  INSITU,
  MERGER
} SFtype;

#ifdef __cplusplus
extern "C"
{
#endif
#if USE_2DISK_MODEL
  void update_reservoirs_from_sf(struct galaxy_t* gal, double new_starsD1, double new_starsD2, int snapshot, SFtype type);
#else
  void update_reservoirs_from_sf(struct galaxy_t* gal, double new_stars, int snapshot, SFtype type);
#endif
  void insitu_star_formation(struct galaxy_t* gal, int snapshot);
  double pressure_dependent_star_formation(struct galaxy_t* gal, int snapshot);

#ifdef __cplusplus
}
#endif

#endif
