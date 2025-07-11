#ifndef _INIT_MERAXES
#define _INIT_MERAXES

#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <hdf5.h>
#include <mlog.h>
#include <stdio.h>

/*
 * Definitions
 */

#include "meraxes_conf.h"

/*
 * Units (cgs)
 */

#define GRAVITY 6.672e-8
#define SOLAR_MASS 1.989e33
#define SOLAR_LUM 3.826e33
#define RAD_CONST 7.565e-15
#define AVOGADRO 6.0222e23
#define BOLTZMANN 1.3806e-16
#define GAS_CONST 8.31425e7
#define SPEED_OF_LIGHT 2.9979e10 // [cm/s]
#define PLANCK 6.6262e-27        //! [erg/s]
#define PROTONMASS 1.6726e-24
#define HUBBLE 3.2407789e-18 //! [h/sec]
#define SEC_PER_MEGAYEAR 3.155e13
#define SEC_PER_YEAR 3.155e7
#define MPC 3.086e24
#define TCMB 2.728
#define NU_LL (double)(3.29e15)
#define NU_LW (double)(2.71e15)
#define PLANCK_EV (double)(4.1357e-15)

// ======================================================
// Don't change these unless you know what you are doing!
#define STRLEN 256 //!< Default string length
#define REL_TOL (float)1e-5
#define ABS_TOL (float)1e-8
// ======================================================

// Define things used for aborting exceptions
#ifdef __cplusplus
extern "C"
{
#endif
  void myexit(int signum);
#ifdef __cplusplus
}
#endif
#define ABORT(sigterm)                                                                                                 \
  do {                                                                                                                 \
    fprintf(stderr, "\nIn file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__);                          \
    myexit(sigterm);                                                                                                   \
  } while (0)

//! Physics parameter values: duplicate these for Pop. III / Pop. II SF
typedef struct physics_params_t
{
  double SfEfficiency;
  double SfEfficiencyScaling;
  double SfCriticalSDNorm;
  double SfRecycleFraction;
  int SnModel;
  double SnReheatRedshiftDep;
  double SnReheatEff;
  double SnReheatLimit;
  double SnReheatScaling;
  double SnReheatScaling2;
  double SnReheatNorm;
  double SnEjectionRedshiftDep;
  double SnEjectionEff;
  double SnEjectionScaling;
  double SnEjectionScaling2;
  double SnEjectionNorm;
  double MaxCoolingMassFactor;
  int ReincorporationModel;
  double ReincorporationEff;
  double Yield;
  double RadioModeEff;
  double QuasarModeEff;
  double BlackHoleGrowthRate;
  double EddingtonRatio;
  double quasar_mode_scaling;
  double quasar_open_angle;
  double quasar_fobs;

  double SfCriticalSDNorm_III;
  double SfEfficiency_III;
  double SfEfficiencyScaling_III;
  double SfRecycleFraction_III;
  double SnReheatRedshiftDep_III;
  double SnReheatEff_III;
  double SnReheatLimit_III;
  double SnReheatScaling_III;
  double SnReheatScaling2_III;
  double SnReheatNorm_III;
  double SnEjectionRedshiftDep_III;
  double SnEjectionEff_III;
  double SnEjectionScaling_III;
  double SnEjectionScaling2_III;
  double SnEjectionNorm_III;
  double Yield_III;

  double ThreshMajorMerger;
  double MinMergerStellarMass;
  double MinMergerRatioForBurst;
  double MergerBurstScaling;
  double MergerBurstFactor;
  double MergerTimeFactor;

  // TODO: These parameters should be used to set the TOCF HII_EFF_FACTOR value
  double ReionEfficiency;
  double ReionEfficiencyIII;
  double ReionNionPhotPerBary;
  double ReionNionPhotPerBaryIII;
  double BlackHoleSeed;
  double BlackHoleMassLimitReion;
  double ReionTcool;
  double Y_He;

  // Pop III IMF parameters
  double MminIMF;
  double MmaxIMF;
  double AlphaIMF;
  double McharIMF; // for LogNorm
  double SigmaIMF;
  int PopIII_IMF;
  int PopIIIAgePrescription;

  // Parameters for metallicity evolution. These include the parameters for the fitting function of clustering and the
  // critical metallicity
  double AlphaCluster;
  double BetaCluster;
  double GammaCluster;
  double NormCluster;

  double ZCrit;

  // Parameters to describe the X-ray properties of the sources
  double LXrayGal;
  double NuXrayGalThreshold;
  double SpecIndexXrayGal;
  double LXrayGalIII;
  double SpecIndexXrayIII;
  double NuXraySoftCut;
  double NuXrayMax;

  double ReionMaxHeatingRedshift;

  double ReionGammaHaloBias;
  double ReionAlphaUV;
  double ReionAlphaUVBH;
  double ReionRBubbleMin;
  double ReionRBubbleMax;
  double ReionRBubbleMaxRecomb;

  double EscapeFracNorm;
  double EscapeFracNormIII;
  double EscapeFracRedshiftOffset;
  double EscapeFracRedshiftScaling;
  double EscapeFracPropScaling;
  double EscapeFracBHNorm;
  double EscapeFracBHScaling;

  // global reionization prescription
  double ReionSobacchi_Zre;
  double ReionSobacchi_DeltaZre;
  double ReionSobacchi_DeltaZsc;
  double ReionSobacchi_T0;

  // global reionization prescription
  double ReionGnedin_z0;
  double ReionGnedin_zr;

  // filtering mass fit
  double ReionSMParam_m0;
  double ReionSMParam_a;
  double ReionSMParam_b;
  double ReionSMParam_c;
  double ReionSMParam_d;

  // options
  int EscapeFracDependency;
  int SfDiskVelOpt;
  int SfPrescription;
  bool InstantSfIII;

  // Flags
  int Flag_ReionizationModifier;
  int Flag_BHFeedback;
  int Flag_IRA;
  int Flag_FixDiskRadiusOnInfall;
  int Flag_FixVmaxOnInfall;
  int Flag_ReheatToFOFGroupTemp;
} physics_params_t;

enum tree_ids
{
  VELOCIRAPTOR_TREES,
  GBPTREES_TREES,
  VELOCIRAPTOR_TREES_AUG
};

//! Run params
//! Everything in this structure is supplied by the user...
typedef struct run_params_t
{
  char DefaultsFile[STRLEN];
  char SimParamsFile[STRLEN];
  char OutputDir[STRLEN];
  char FileNameGalaxies[STRLEN];
  char SimName[STRLEN];
  char SimulationDir[STRLEN];
  char CatalogFilePrefix[STRLEN];
  char OutputSnapsString[STRLEN];
  char PhotometricTablesDir[STRLEN];
  char TargetSnaps[STRLEN];
  char BetaBands[STRLEN];
  char RestBands[STRLEN];
  double BirthCloudLifetime;
  double DeltaT; // New Parameter added to consider different time of observation! Very important for Pop. III
  char CoolingFuncsDir[STRLEN];
  char StellarFeedbackDir[STRLEN];
  char TablesForXHeatingDir[STRLEN];
  char IMF[STRLEN];
  char MagSystem[STRLEN];
  char MagBands[STRLEN];
  char ForestIDFile[STRLEN];
  char MvirCritFile[STRLEN];
  char MvirCritMCFile[STRLEN];
  char MassRatioModifier[STRLEN];
  char BaryonFracModifier[STRLEN];
  char FFTW3WisdomDir[STRLEN];

  physics_params_t physics;

  double BoxSize;
  double VolumeFactor;
  double Hubble_h;
  double BaryonFrac;
  double OmegaM;
  double OmegaK;
  double OmegaR;
  double OmegaLambda;
  double Sigma8;
  double wLambda;
  double SpectralIndex;
  double PartMass;
  long long NPart;

  double* MvirCrit;
  double* MvirCrit_MC;

  double ReionDeltaRFactor;
  double ReionPowerSpecDeltaK;
  int ReionGridDim;
  int ReionFilterType;
  int TsHeatingFilterType;
  int ReionRtoMFilterType;
  int ReionUVBFlag;
  int MetalGridDim;

  enum tree_ids TreesID;
  int FirstFile;
  int LastFile;
  int NSteps;
  int SnaplistLength;
  int RandomSeed;
  int FlagSubhaloVirialProps;
  int FlagInteractive;
  int FlagMCMC;
  int Flag_PatchyReion;
  int Flag_IncludeSpinTemp;
  int Flag_IncludeLymanWerner;
  int Flag_IncludeStreamVel;
  int Flag_IncludeMetalEvo; // New for Metallicity
  int Flag_IncludeRecombinations;
  int Flag_Compute21cmBrightTemp;
  int Flag_ComputePS;
  int Flag_IncludePecVelsFor21cm;
  int Flag_ConstructLightcone;

  int TsVelocityComponent;
  int TsNumFilterSteps;

  double ReionSfrTimescale;

  double EndRedshiftLightcone;
  int EndSnapshotLightcone;
  long long LightconeLength;
  long long CurrentLCPos;
  int PS_Length;
  int Flag_OutputGrids;
  int Flag_OutputGridsPostReion;
  int FlagIgnoreProgIndex;
} run_params_t;

typedef struct run_units_t
{
  double UnitTime_in_s;
  double UnitLength_in_cm;
  double UnitVelocity_in_cm_per_s;
  double UnitTime_in_Megayears;
  double UnitMass_in_g;
  double UnitDensity_in_cgs;
  double UnitPressure_in_cgs;
  double UnitCoolingRate_in_cgs;
  double UnitEnergy_in_cgs;
  // TOTAL : 72  (must be multiple of 8)
} run_units_t;

typedef struct hdf5_output_t
{
  char** params_tag;
  void** params_addr;
  int* params_type;
  size_t* dst_offsets;
  size_t* dst_field_sizes;
  const char** field_names;
  const char** field_units;
  const char** field_h_conv;
  hid_t* field_types;
  size_t dst_size;
  hid_t array3f_tid; // sizeof(hid_t) = 4
  hid_t array_nmag_f_tid;
  hid_t array_nhist_f_tid;
  int n_props;
  int params_count;

  // TOTAL : 52 + 4 padding (must be multiple of 8)
} hdf5_output_t;

#if USE_MINI_HALOS
typedef struct metal_grids_t // New stuff for MetalEvo, probably you will need to add new stuff
{
  ptrdiff_t* slab_nix_metals;
  ptrdiff_t* slab_ix_start_metals;

  float* buffer_metals;

  struct gal_to_slab_t* galaxy_to_slab_map_metals;

  double volume_ave_ZIGM;        // Is it necessary? Maybe this one could be used as a log message
  double volume_ave_mass_metals; // Is it necessary?

  int buffer_size_metals;

  float* N_bubbles;
  float* mass_IGM; // It's the total! Computed from the overdensity and you sum the galaxy net contribution (mass_gas)
  float* mass_metals;
  float* mass_gas; // The one coming from galaxies (Ejected - Hot - Cold)
  float* Zigm_box;
  float* Probability_metals;
  float* R_ave; // Average radius in the bubble
  float* R_max; // Max radius in the bubble

} metal_grids_t;
#endif

typedef struct reion_grids_t
{
  ptrdiff_t* slab_nix;
  ptrdiff_t* slab_ix_start;
  ptrdiff_t* slab_n_complex;

  float* buffer;

  float* stars;
  fftwf_complex* stars_unfiltered;
  fftwf_complex* stars_filtered;
  fftwf_plan stars_forward_plan;
  fftwf_plan stars_filtered_reverse_plan;

  float* deltax;
  fftwf_complex* deltax_unfiltered;
  fftwf_complex* deltax_filtered;
  fftwf_plan deltax_forward_plan;
  fftwf_plan deltax_filtered_reverse_plan;

  float* sfr;
  float* sfr_histories;
  float* weighted_sfr;
  fftwf_complex* sfr_unfiltered;
  fftwf_complex* sfr_filtered;
  fftwf_complex* weighted_sfr_unfiltered;
  fftwf_complex* weighted_sfr_filtered;
  fftwf_plan sfr_forward_plan;
  fftwf_plan weighted_sfr_forward_plan;
  fftwf_plan sfr_filtered_reverse_plan;
  fftwf_plan weighted_sfr_filtered_reverse_plan;

#if USE_MINI_HALOS
  float* starsIII;
  fftwf_complex* starsIII_unfiltered;
  fftwf_complex* starsIII_filtered;
  fftwf_plan starsIII_forward_plan;
  fftwf_plan starsIII_filtered_reverse_plan;

  float* sfrIII;
  float* sfrIII_histories;
  float* weighted_sfrIII;
  fftwf_complex* sfrIII_unfiltered;
  fftwf_complex* sfrIII_filtered;
  fftwf_complex* weighted_sfrIII_unfiltered;
  fftwf_complex* weighted_sfrIII_filtered;
  fftwf_plan sfrIII_forward_plan;
  fftwf_plan weighted_sfrIII_forward_plan;
  fftwf_plan sfrIII_filtered_reverse_plan;
  fftwf_plan weighted_sfrIII_filtered_reverse_plan;
#endif

  float* xH;
  float* z_at_ionization;
  float* J_21_at_ionization;
  float* J_21;
  float* Mvir_crit;
#if USE_MINI_HALOS
  float* Mvir_crit_MC;
#endif
  float* r_bubble;

  // Grids necessary for the IGM spin temperature
  float* x_e_box;
  fftwf_complex* x_e_unfiltered;
  fftwf_complex* x_e_filtered;
  fftwf_plan x_e_box_forward_plan;
  fftwf_plan x_e_filtered_reverse_plan;

  float* x_e_box_prev;
  float* Tk_box;
  float* Tk_box_prev;
  float* TS_box;

#if USE_MINI_HALOS
  float* Tk_boxII;
  float* Tk_box_prevII;
  float* TS_boxII;
#endif

  double* SMOOTHED_SFR_GAL;
#if USE_MINI_HALOS
  double* SMOOTHED_SFR_III;
#endif

  // Grids necessary for LW background and future disentangling between MC/AC Pop3/Pop2 stuff

#if USE_MINI_HALOS
  float* JLW_box;
  float* JLW_boxII;
#endif

  // Grids necessary for inhomogeneous recombinations
  float* z_re;

  float* N_rec;
  fftwf_complex* N_rec_unfiltered;
  fftwf_complex* N_rec_filtered;
  fftwf_plan N_rec_forward_plan;
  fftwf_plan N_rec_filtered_reverse_plan;

  float* Gamma12;

  // Grids necessary for the 21cm brightness temperature
  float* delta_T;
  float* delta_T_prev;
#if USE_MINI_HALOS
  float* delta_TII;
  float* delta_TII_prev;
#endif
  float* vel;
  fftwf_complex* vel_gradient;
  fftwf_plan vel_forward_plan;
  fftwf_plan vel_gradient_reverse_plan;

  // Grid for the lightcone (cuboid) box
  float* LightconeBox;
  float* Lightcone_redshifts;

  // Data for the power spectrum
  float* PS_k;
  float* PS_data;
  float* PS_error;

#if USE_MINI_HALOS
  float* PSII_data;
  float* PSII_error;
#endif

  struct gal_to_slab_t* galaxy_to_slab_map;

  double volume_weighted_global_xH;
  double volume_weighted_global_J_21;
  double mass_weighted_global_xH;

  double volume_ave_J_alpha;
  double volume_ave_xalpha;
  double volume_ave_Xheat;
  double volume_ave_Xion;
  double volume_ave_TS;
  double volume_ave_TK;
  double volume_ave_xe;
  double volume_ave_Tb;
#if USE_MINI_HALOS
  double volume_ave_J_alphaII;
  double volume_ave_J_LW;
  double volume_ave_J_LWII;
  double volume_ave_XheatII;
  double volume_ave_TSII;
  double volume_ave_TKII;
  double volume_ave_TbII;
#endif

  int started;
  int finished;
  int buffer_size;
  bool flag_wisdom;
} reion_grids_t;

typedef struct galaxy_t
{
  double NewStars[N_HISTORY_SNAPS];
#if USE_MINI_HALOS
  double NewStars_II[N_HISTORY_SNAPS]; // New
  double NewStars_III[N_HISTORY_SNAPS];
#endif
  double NewMetals[N_HISTORY_SNAPS];

#ifdef CALC_MAGS
  double inBCFlux[MAGS_N];
  double outBCFlux[MAGS_N];
#if USE_MINI_HALOS
  double inBCFluxIII[MAGS_N];
  double outBCFluxIII[MAGS_N];
#endif
#endif

  // Unique ID for the galaxy
  unsigned long ID;

  // properties of subhalo at the last time this galaxy was a central galaxy
  float Pos[3];
  float Vel[3];
  double Mvir;
  double Rvir;
  double Vvir;
  double Vmax;
  double Spin;

  double dt; //!< Time between current snapshot and last identification

  struct halo_t* Halo;
  struct galaxy_t* FirstGalInHalo;
  struct galaxy_t* NextGalInHalo;
  struct galaxy_t* Next;
  struct galaxy_t* MergerTarget;

  // baryonic reservoirs
  double HotGas;
  double MetalsHotGas;
  double ColdGas;
  double MetalsColdGas;
  double H2Frac;
  double H2Mass;
  double HIMass;
  double Mcool;
  double StellarMass;
  double GrossStellarMass;
  double Fesc;
  double FescWeightedGSM;
  double MetalsStellarMass;
  double DiskScaleLength;
  double Sfr;
  double EjectedGas;
  double MetalsEjectedGas;
  double BlackHoleMass;
  double FescBH;
  double BHemissivity;
  double EffectiveBHM;
  double BlackHoleAccretedHotMass;
  double BlackHoleAccretedColdMass;
  double BlackHoleAccretingColdMass;

  int Galaxy_Population; // You need it also if you are not disentangling PopIII/PopII (when Mini_halos is off, this is
                         // = 2)
#if USE_MINI_HALOS
  // Differentiation Pop III / Pop II
  double StellarMass_II;
  double StellarMass_III;
  double GrossStellarMassIII;
  double FescIII;
  double FescIIIWeightedGSM;

  double Remnant_Mass; // Coming from Pop III with M between 40 and 140 and larger than 260 Msol and remnant of CCSN
                       // [8,40]Msun. Atm those are silent.

  // Metallicity stuff
  double RmetalBubble; // New for MetalEvo
  double PrefactorBubble;
  double TimeBubble;
  double Metal_Probability;    // Probability to be polluted by other metal bubbles
  double GalMetal_Probability; // Random number between 0 and 1 associated to the galaxy.
  double Metals_IGM;
  double Gas_IGM;
  double Metallicity_IGM; //
  double MaxBubble;       // Need this for Boost probability
  double AveBubble;   // Same (you will likely use only one of the two). You could actually save only the boost factor
  int Flag_ExtMetEnr; // 0 if not enriched from a bubble, 1 yes.

  double Prefactor[N_HISTORY_SNAPS]; // here you store the prefactors of the metal bubbles
  double Times[N_HISTORY_SNAPS];     // Time at which the SN explode!
  double Radii[N_HISTORY_SNAPS];
#endif

  // baryonic hostories
  double mwmsa_num;
  double mwmsa_denom;

  // misc
  double Rcool;
  double Cos_Inc;
  double MergTime;
  double MergerStartRadius;
  double BaryonFracModifier;
  double FOFMvirModifier;
  double MvirCrit;
  double MvirCrit_MC;
  double MergerBurstMass;

  int Type;
  int OldType;
  int Len;
  int MaxLen;
  int SnapSkipCounter;
  int HaloDescIndex;
  int TreeFlags;
  int LastIdentSnap; //!< Last snapshot at which the halo in which this galaxy resides was identified
  int output_index;  //!< write index

  bool ghost_flag;

  // N.B. There will be padding present at the end of this struct, but amount
  // is dependent on CALC_MAGS, MAX_PHOTO_NBANDS, NOUT and N_HISTORY_SNAPS.
} galaxy_t;

//! The meraxes halo structure
typedef struct halo_t
{
  struct fof_group_t* FOFGroup;
  struct halo_t* NextHaloInFOFGroup;
  galaxy_t* Galaxy;

  float Pos[3]; //!< Most bound particle position [Mpc/h]
  float Vel[3]; //!< Centre of mass velocity [Mpc/h]
  float AngMom; //!< Specific angular momentum length [Mpc/h *km/s]

  double Mvir; //!< virial mass [M_sol/h]
  double Rvir; //!< Virial radius [Mpc/h]
  double Vvir; //!< Virial velocity [km/s]

  float Vmax;       //!< Maximum circular velocity [km/s]
  unsigned long ID; //!< Halo ID
  int Type;         //!< Type (0 for central, 1 for satellite)
  int SnapOffset;   //!< Number of snapshots this halo skips before reappearing
  int DescIndex;    //!< Index of descendant in next relevant snapshot
  int ProgIndex;    //!< Index of progenitor in previous relevant snapshot
  int TreeFlags;    //!< Bitwise flag indicating the type of match in the trees
  int Len;          //!< Number of particles in the structure
} halo_t;

typedef struct fof_group_t
{
  halo_t* FirstHalo;
  halo_t* FirstOccupiedHalo;
  double Mvir;
  double Rvir;
  double Vvir;
  double FOFMvirModifier;
  int TotalSubhaloLen;
} fof_group_t;

//! Tree info struct
typedef struct trees_info_t
{
  int n_halos;
  int n_halos_max;
  int max_tree_id;
  int n_fof_groups;
  int n_fof_groups_max;
} trees_info_t;

// This structure carries the information about
//   the GPU allocated to this CPU's scope. It
//   needs to be declared by all compilers since
//   it is used by run_globals.
#ifdef USE_CUDA
#include <cuda_runtime.h>
typedef struct gpu_info
{
  int device;                       // the ordinal of the current context's device
  struct cudaDeviceProp properties; // Properties of this context's assigned device
  int n_threads;                    // No. of threads to use in kernal calls
  int n_contexts;                   // No. of ranks with successfully allocated GPU contexts
} gpu_info;
#else
typedef char gpu_info;
#endif

#ifdef CALC_MAGS
typedef struct mag_params_t
{
  int targetSnap[MAGS_N_SNAPS];
  int nBeta;
  int nRest;
  int minZ;
  int maxZ;
  int nMaxZ;
  double tBC;
  int iAgeBC[MAGS_N_SNAPS];
  size_t totalSize;
  double* working;
  double* inBC;
  double* outBC;
  double* centreWaves;
  double allcentreWaves[MAGS_N_SNAPS][MAGS_N_BANDS];
  double* logWaves;
#ifdef USE_MINI_HALOS
  size_t totalSizeIII;
  double* workingIII;
#endif
} mag_params_t;
#endif

//! Global variables which will will be passed around
typedef struct run_globals_t
{
  struct run_params_t params;
  char FNameOut[STRLEN];
  reion_grids_t reion_grids;
#if USE_MINI_HALOS
  metal_grids_t metal_grids;
#endif
  struct run_units_t units;
  hdf5_output_t hdf5props;

  MPI_Comm mpi_comm;
  int mpi_rank;
  int mpi_size;
  gpu_info* gpu;

  double* AA;
  double* ZZ;
  double* LTTime;
  double* rhocrit;
  long* RequestedForestId;
  int RequestedMassRatioModifier;
  int RequestedBaryonFracModifier;
  int* ListOutputSnaps;
  halo_t** SnapshotHalo;
  fof_group_t** SnapshotFOFGroup;
  int** SnapshotIndexLookup;
  float** SnapshotDeltax;
  float** SnapshotVel;
  trees_info_t* SnapshotTreesInfo;
  galaxy_t* FirstGal;
  galaxy_t* LastGal;
  gsl_rng* random_generator;
  void* mhysa_self;
  double Hubble;
  double RhoCrit;
  double G;
  double Csquare;
  // PopIII stuff

  double IMFnorm;
  double NumberPISN;
  double MassPISN;
  double NumberSNII;
  double MassSNII;
  double MassBHs;

  float* Mass_Values;
  float* Time_Values;

#ifdef CALC_MAGS
  struct mag_params_t mag_params;
#endif

  int NOutputSnaps;
  int LastOutputSnap;
  int NGhosts;
  int NHalosMax;
  int NFOFGroupsMax;
  int NRequestedForests;
  int NStoreSnapshots;
  int NstoreSnapshots_SFR;

  bool SelectForestsSwitch;
  struct Modifier* mass_ratio_modifier;
  struct Modifier* baryon_frac_modifier;
} run_globals_t;
#ifdef _MAIN
run_globals_t run_globals;
#else
extern run_globals_t run_globals;
#endif

/*
 * Functions
 */

#ifdef __cplusplus
extern "C"
{
#endif

  // core/dracarys.c
  void dracarys(void);

  // core/read_parameter_file.c
  void read_parameter_file(char* fname, int mode);

  // core/init.c
  void init_storage(void);
  void init_meraxes(void);

  int set_sfr_history(void);
  // core/cleanup.c
  void cleanup(void);

  // core/magnitudes.c
  void get_output_magnitudes(float* mags, float* dusty_mags, galaxy_t* gal, int snapshot);
  void get_output_magnitudesIII(float* mags, galaxy_t* gal, int snapshot);

// MCMC related
// meraxes_mhysa_hook must be implemented by the calling code (Mhysa)!
#ifdef _MAIN
  int (*meraxes_mhysa_hook)(void* self, int snapshot, int ngals);
#else
extern int (*meraxes_mhysa_hook)(void* self, int snapshot, int ngals);
#endif

#ifdef USE_CUDA
// This stuff is needed by the GPU routines.  This
//    needs to be included after mlog.h is included
//    and after ABORT(), myexit() & run_globals are
//    defined, since they are used within.
#include "meraxes_gpu.h"
#endif

#ifdef __cplusplus
}
#endif

#endif // _INIT_MERAXES
