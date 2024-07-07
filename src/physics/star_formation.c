#include <assert.h>
#include <gsl/gsl_integration.h>

#include "misc_tools.h"
#include "meraxes.h"
#include "star_formation.h"
#include "supernova_feedback.h"

static void backfill_ghost_star_formation(galaxy_t* gal, double m_stars, double sfr, double metallicity, int snapshot)
{
  double* LTTime = run_globals.LTTime;
  double burst_time = LTTime[gal->LastIdentSnap] - gal->dt * 0.5;

  for (int snap = snapshot - 1, ii = 1; snap >= gal->LastIdentSnap; --snap, ++ii) {

    if (LTTime[snap] > burst_time) {
      if (ii < N_HISTORY_SNAPS) {
        gal->NewStars[ii] += m_stars;
        gal->NewMetals[0] += m_stars * metallicity;
      }
      break;
    }
  }
}

void update_reservoirs_from_sf(galaxy_t* gal, double new_stars, int snapshot, SFtype type)
{
  if (new_stars > 0) {
    double metallicity;
    bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);

    // instantaneous recycling approximation of stellar mass
    metallicity = calc_metallicity(gal->ColdGas, gal->MetalsColdGas);

    // update the galaxy's SFR value
    double sfr = new_stars / gal->dt;
    gal->Sfr += sfr;
    assert(gal->Sfr >= 0);

    gal->ColdGas -= new_stars;
    gal->MetalsColdGas -= new_stars * metallicity;
    gal->StellarMass += new_stars;
    gal->MetalsStellarMass += new_stars * metallicity;

    gal->GrossStellarMass += new_stars; 

    if ((type == INSITU) && !Flag_IRA && (gal->LastIdentSnap < (snapshot - 1))) {
      // If this is a reidentified ghost, then back fill NewStars and
      // escape fraction dependent properties to reflect this new insitu
      // SF burst.
      backfill_ghost_star_formation(gal, new_stars, sfr, metallicity, snapshot);
    } else {
      // update the stellar mass history assuming the burst is happening in this snapshot
      gal->NewStars[0] += new_stars;
      gal->NewMetals[0] += new_stars * metallicity;
    }

    // Check the validity of the modified reservoir values.
    // Note that the ColdGas reservers *can* be negative at this point.  This
    // is because some fraction of the stars in this burst will go nova and
    // return mass to the ISM.  This will be accounted for when we update the
    // reservoirs due to supernova feedback.
    if (gal->StellarMass < 0)
      gal->StellarMass = 0.0;
    if (gal->MetalsStellarMass < 0)
      gal->MetalsStellarMass = 0.0;
  }
}

void insitu_star_formation(galaxy_t* gal, int snapshot)
{
  // there is no point doing anything if there is no cold gas!
  if (gal->ColdGas > 1e-10) {
    double r_disk;
    double v_disk;
    double m_crit;
    double m_stars;
    double m_reheat;
    double m_eject;
    double m_recycled;
    double new_metals;
    double m_remnant;
    double zplus1;
    double zplus1_n;

    zplus1 = 1.0 + run_globals.ZZ[snapshot];
    zplus1_n = pow(zplus1, run_globals.params.physics.SfEfficiencyScaling);

    double SfEfficiency_II = run_globals.params.physics.SfEfficiency;
    double SfCriticalSDNorm = run_globals.params.physics.SfCriticalSDNorm;
    int SfDiskVelOpt = run_globals.params.physics.SfDiskVelOpt;
    int SfPrescription = run_globals.params.physics.SfPrescription;

    // What velocity are we going to use as a proxy for the disk rotation velocity?
    switch (SfDiskVelOpt) {
      case 1:
        v_disk = gal->Vmax;
        break;
      case 2:
        v_disk = gal->Vvir;
        break;
      default:
        mlog_error("Unrecognised value for SfVelocityOpt parameter! Defaulting to v_disk=Vmax.");
        v_disk = gal->Vmax;
        break;
    }

    // calculate disk scalelength using Mo, Mao & White (1998) eqn. 12 and
    // multiply it by 3 to approximate the star forming region size (ala
    // Croton+ 2006).
    r_disk = gal->DiskScaleLength * 3.0;

    switch (SfPrescription) {
      case 1:
        // what is the critical mass within r_crit?
        // from Kauffmann (1996) eq7 x piR^2, (Vvir in km/s, reff in Mpc/h) in units of 10^10Msun/h
        m_crit = SfCriticalSDNorm * v_disk * r_disk;
        if (gal->ColdGas > m_crit)
          m_stars = zplus1_n * SfEfficiency_II * (gal->ColdGas - m_crit) / r_disk * v_disk * gal->dt;
        else
          // no star formation
          return;
        break;

      case 2:
        // GALFORM
        m_stars = gal->ColdGas / (r_disk / v_disk / 0.029 * pow(200. / v_disk, 1.5)) * gal->dt;
        break;

      default:
        m_stars = 0;
        mlog_error("Unknown SfPrescription!");
        ABORT(EXIT_FAILURE);
        break;
    }
    if (m_stars > gal->ColdGas)
      m_stars = gal->ColdGas;
    // calculate the total supernova feedback which would occur if this star
    // formation happened continuously and evenly throughout the snapshot
    contemporaneous_supernova_feedback(
      gal, &m_stars, snapshot, &m_reheat, &m_eject, &m_recycled, &m_remnant, &new_metals);
    // update the baryonic reservoirs (note that the order we do this in will change the result!)
    update_reservoirs_from_sf(gal, m_stars, snapshot, INSITU);
    update_reservoirs_from_sn_feedback(gal, m_reheat, m_eject, m_recycled, 0, m_recycled, m_remnant, new_metals);
  }
}
