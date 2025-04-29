#ifndef SUPERNOVA_FEEDBACK_H
#define SUPERNOVA_FEEDBACK_H

#include "meraxes.h"
#ifdef __cplusplus
extern "C"
{
#endif
#if USE_2DISK_MODEL
  void update_reservoirs_from_sn_feedback(struct galaxy_t* gal,
                                          double m_reheat,
                                          double m_eject,
                                          double m_recycled_II,
                                          double m_recycled_III,
                                          double m_remnant,
                                          double new_metals);
#else
  void update_reservoirs_from_sn_feedback(struct galaxy_t* gal,
                                          double m_reheat,
                                          double m_eject,
                                          double m_recycled,
                                          double m_recycled_III,
                                          double m_recycled_II,
                                          double m_remnant,
                                          double new_metals);
#endif
  void delayed_supernova_feedback(struct galaxy_t* gal, int snapshot);
  void calc_metal_bubble(struct galaxy_t* gal, int snapshot);
#if USE_2DISK_MODEL
  void contemporaneous_supernova_feedback(struct galaxy_t* gal,
                                          double* m_stars,
                                          double* m_stars2,
                                          int snapshot,
                                          double* m_reheat,
                                          double* m_eject,
                                          double* m_recycled,
                                          double* m_recycled2,
                                          double* m_remnant,
                                          double* new_metals);
#else
  void contemporaneous_supernova_feedback(struct galaxy_t* gal,
                                          double* m_stars,
                                          int snapshot,
                                          double* m_reheat,
                                          double* m_eject,
                                          double* m_recycled,
                                          double* m_remnant,
                                          double* new_metals);
#endif

#ifdef __cplusplus
}
#endif

#endif
