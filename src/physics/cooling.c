#include <hdf5_hl.h>
#include <assert.h>
#include <math.h>

#include "blackhole_feedback.h"
#include "cooling.h"
#include "misc_tools.h"
#include "virial_properties.h"
#include "meraxes.h"
#include "reionization.h"

// *** This code is taken from the model of Croton+ 2006 with minimal modifications. ***
//
// Future todo list:
//  - Read in the below hard coded values from the HDF5 file itself and malloc appropriately
//  - Set up more flexible interpolation (e.g. without assumption of fixed temp step)
//  - Use GSL routines for interpolation

// metallicies with repect to solar.
// These will be converted to absolute metallicities by adding log10(Z_sun), Zsun=0.02
static double metallicities[N_METALLICITIES] = {
  -5.0, // This is actually primordial in the tables but that is log10(0) = -infinity
  -3.0, -2.0, -1.5, -1.0, -0.5, +0.0, +0.5
};

static char group_name[N_METALLICITIES][6] = { "mzero", "m-30", "m-20", "m-15", "m-10", "m-05", "m-00", "m+05" };

static double cooling_rate[N_METALLICITIES][N_TEMPS];

void read_cooling_functions()
{
  if (run_globals.mpi_rank == 0) {
    hid_t fd;
    char dset_name[30];
    char fname[STRLEN + 11];

    sprintf(fname, "%s/SD93.hdf5", run_globals.params.CoolingFuncsDir);
    fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    for (int i_m = 0; i_m < N_METALLICITIES; i_m++) {
      sprintf(dset_name, "%s/log(lambda_norm)", group_name[i_m]);
      H5LTread_dataset_double(fd, dset_name, cooling_rate[i_m]);
    }

    H5Fclose(fd);
  }

  // broadcast the values to all cores
  MPI_Bcast(&cooling_rate, sizeof(cooling_rate), MPI_BYTE, 0, run_globals.mpi_comm);

  // add solar metallicity to all metallicity values
  for (int i_m = 0; i_m < N_METALLICITIES; i_m++)
    metallicities[i_m] += log10(0.02);
}

static double interpolate_temp_dependant_cooling_rate(int i_m, double logTemp)
{
  int i_t;
  double temp_step;
  double rate, rate_below, rate_above;
  double logT_below;

  // First deal with boundary conditions
  if (logTemp < MIN_TEMP)
    return -27.0;

  // Now find the index of the tabulated temp immediately below our value
  temp_step = (MAX_TEMP - MIN_TEMP) / (double)(N_TEMPS - 1);
  i_t = (int)((logTemp - MIN_TEMP) / temp_step);
  if (i_t >= (N_TEMPS - 1))
    i_t = N_TEMPS - 2;

  // Now grab the cooling rates for the temp immediately above and below
  rate_below = cooling_rate[i_m][i_t];
  rate_above = cooling_rate[i_m][i_t + 1];

  // Calculate the tabulated temperature immediately below our value
  logT_below = MIN_TEMP + temp_step * i_t;

  // Now linearly interpolate the cooling rate
  rate = rate_below + (rate_above - rate_below) / temp_step * (logTemp - logT_below);

  return rate;
}

double interpolate_cooling_rate(double logTemp, double logZ)
{
  int i_m;
  double rate_below, rate_above, rate;

  // First deal with boundary conditions
  if (logZ < metallicities[0])
    logZ = metallicities[0];
  if (logZ > metallicities[N_METALLICITIES - 1])
    logZ = metallicities[N_METALLICITIES - 1];

  // Now find the indices of the metallicity values which bound our input
  i_m = 0;
  while (logZ > metallicities[i_m + 1])
    i_m++;

  // Get the cooling rates for this temperature value
  rate_below = interpolate_temp_dependant_cooling_rate(i_m, logTemp);
  rate_above = interpolate_temp_dependant_cooling_rate(i_m + 1, logTemp);

  // Finally, linearly interpolate the cooling rates
  rate = rate_below +
         (rate_above - rate_below) / (metallicities[i_m + 1] - metallicities[i_m]) * (logZ - metallicities[i_m]);

  return pow(10, rate);
}

double gas_cooling(galaxy_t* gal)
{
  double cooling_mass = 0.0;

  // we only need to do cooling if there is anything to cool!
  if (gal->HotGas > 1e-10) {
    fof_group_t* fof_group = gal->Halo->FOFGroup;
    int halo_type = 1;
    
    // calculate the halo virial temperature and log10 metallicity value
    double Tvir = Vvir_to_Tvir(fof_group->Vvir, halo_type);
    double log10Tvir = log10(Tvir);
    double logZ;
    double t_cool, max_cooling_mass;
    double lambda, x, rho_r_cool, r_cool, isothermal_norm;
    run_units_t* units = &(run_globals.units);
    double max_cooling_mass_factor = run_globals.params.physics.MaxCoolingMassFactor;

    if (gal->MetalsHotGas > 0)
      logZ = log10(calc_metallicity(gal->HotGas, gal->MetalsHotGas));
    else
      logZ = -10.0;

    if (Tvir >= 1e4) {

      t_cool = fof_group->Rvir / fof_group->Vvir; // internal units

      // interpolate the temperature and metallicity dependant cooling rate (lambda)
      lambda = interpolate_cooling_rate(log10Tvir, logZ);

      x = PROTONMASS * BOLTZMANN * Tvir / lambda;              // now this has units sec g/cm^3
      x /= (units->UnitDensity_in_cgs * units->UnitTime_in_s); // now in internal units
      rho_r_cool = x / t_cool * 0.885; // 0.885 = 3/2 * mu, mu=0.59 for a fully ionized gas

      assert(rho_r_cool > 0);
      isothermal_norm = gal->HotGas / (4. * M_PI * fof_group->Rvir);
      r_cool = sqrt(isothermal_norm / rho_r_cool);
      gal->Rcool = r_cool;

      max_cooling_mass = max_cooling_mass_factor * gal->HotGas / t_cool * gal->dt;

      if (r_cool > fof_group->Rvir)
        cooling_mass = max_cooling_mass;
      else {
        cooling_mass = max_cooling_mass / fof_group->Rvir * r_cool;
        if (cooling_mass > max_cooling_mass)
          cooling_mass = max_cooling_mass;
      }

      if (cooling_mass > gal->HotGas)
        cooling_mass = gal->HotGas;

      if (run_globals.params.physics.Flag_BHFeedback)
        cooling_mass -= radio_mode_BH_heating(gal, cooling_mass, x);

      if (cooling_mass < 0)
        cooling_mass = 0.0;
    }
  }
  return cooling_mass;
}

void cool_gas_onto_galaxy(galaxy_t* gal, double cooling_mass)
{
  double cooling_metals;

  if (cooling_mass > gal->HotGas)
    cooling_mass = gal->HotGas;

  // what mass of metals is coming along with this cooling gas?
  cooling_metals = cooling_mass * calc_metallicity(gal->HotGas, gal->MetalsHotGas);

  // save the cooling mass
  gal->Mcool = cooling_mass;

  // update the galaxy reservoirs
  gal->HotGas -= cooling_mass;
  gal->MetalsHotGas -= cooling_metals;
  gal->ColdGas += cooling_mass;
  gal->MetalsColdGas += cooling_metals;
}
