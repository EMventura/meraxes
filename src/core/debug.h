#ifndef DEBUG_H
#define DEBUG_H

#include "meraxes.h"

#ifdef __cplusplus
extern "C"
{
#endif

  void mpi_debug_here(void);
  void check_mhysa_pointer(void);
  void check_counts(fof_group_t* fof_group, int NGal, int NFof);

#ifdef DEBUG
  void check_pointers(halo_t* halos, fof_group_t* fof_groups, trees_info_t* trees_info);
#endif

#ifdef __cplusplus
}
#endif

#endif
