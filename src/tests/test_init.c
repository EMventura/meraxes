#define _MAIN
#include <criterion/criterion.h>
#include <criterion/parameterized.h>
#include <meraxes.h>

// This gives us access to the static functions
#include "../core/init.c"

void setup(void)
{
  const int n_snaps = 201;
  char alist[] = "0.03999999989 0.04065227126 0.04131517907 "
                 "0.04198889605 0.04267359992 0.0433694691 0.04407668567 0.04479543467 "
                 "0.04552590415 0.04626828524 0.04702277218 0.04778956237 0.04856885645 "
                 "0.04936085831 0.05016577517 0.05098381763 0.05181519974 0.05266013902 "
                 "0.05351885653 0.05439157604 0.05527852773 0.05617994275 0.05709605694 "
                 "0.05802711002 0.05897334558 0.05993501119 0.06091235848 0.06190564316 "
                 "0.06291512512 0.06394106848 0.06498374168 0.06604341753 0.06712037328 "
                 "0.06821489072 0.06932725621 0.07045776081 0.07160669908 0.07277437406 "
                 "0.07396109007 0.0751671576 0.07639289221 0.07763861461 0.07890465074 "
                 "0.08019133184 0.08149899458 0.08282798109 0.0841786391 0.08555132199 "
                 "0.08694638893 0.08836420493 0.08980514094 0.09126957399 0.09275788563 "
                 "0.09427046845 0.09580771663 0.09737003238 0.09895782448 0.1005715084 "
                 "0.1022115062 0.1038782472 0.1055721674 0.1072937099 0.1090433253 "
                 "0.1108214712 0.112628613 0.1144652235 0.1163317832 0.1182287805 "
                 "0.1201567096 0.1221160791 0.1241073997 0.1261311922 0.1281879864 "
                 "0.1302783202 0.1324027406 0.1345618035 0.1367560738 0.1389861256 "
                 "0.1412525423 0.1435559171 0.1458968524 0.1482759609 0.150693865 "
                 "0.1531511973 0.1556486008 0.1581867262 0.160766243 0.1633878234 "
                 "0.1660521534 0.16875993 0.1715118618 0.1743086687 0.1771510825 "
                 "0.1800398469 0.1829757179 0.1859594634 0.1889918642 0.1920737137 "
                 "0.1952058183 0.1983889975 0.201624084 0.2049119209 0.2082533754 "
                 "0.2116493182 0.2151006379 0.2186082376 0.2221730349 0.2257959627 "
                 "0.2294779687 0.2332200164 0.2370230848 0.2408881691 0.2448162805 "
                 "0.2488084467 0.2528657123 0.2569891389 0.2611798053 0.2654388035 "
                 "0.2697672567 0.2741662932 0.2786370637 0.2831807382 0.2877985053 "
                 "0.2924915734 0.2972611704 0.3021085441 0.3070349629 0.3120417158 "
                 "0.3171301127 0.3223014849 0.3275571857 0.33289859 0.3383270954 "
                 "0.3438441223 0.3494511082 0.355149532 0.3609408787 0.3668266637 "
                 "0.372808427 0.3788877335 0.385066174 0.3913453649 0.3977269493 "
                 "0.4042125968 0.4108040044 0.4175028966 0.4243110262 0.4312301745 "
                 "0.4382621518 0.4454087981 0.4526719755 0.4600535997 0.4675555944 "
                 "0.4751799225 0.482928579 0.490803591 0.4988070192 0.5069409576 "
                 "0.5152075344 0.5236089124 0.5321472899 0.5408249009 0.5496440158 "
                 "0.5586069421 0.5677160249 0.5769736476 0.5863822222 0.5959442305 "
                 "0.6056621644 0.6155385666 0.6255760211 0.6357771543 0.6461446351 "
                 "0.6566811761 0.6673895343 0.6782725114 0.6893329548 0.7005737585 "
                 "0.7119978635 0.723608259 0.7354079826 0.7474001218 0.7595878143 "
                 "0.7719742355 0.7845626527 0.7973563465 0.8103586643 0.8235730081 "
                 "0.8370028354 0.8506516599 0.8645230529 0.8786206437 0.8929481209 "
                 "0.9075092331 0.9223077902 0.9373476643 0.9526327903 0.9681671675 "
                 "0.9839548605 0.9999999828 1.0";

  run_globals.params.SnaplistLength = n_snaps;
  run_globals.params.Hubble_h = 0.6571;
  run_globals.params.OmegaM = 0.3121;
  run_globals.params.OmegaK = 0.0;
  run_globals.params.OmegaLambda = 0.6879;
  run_globals.units.UnitLength_in_cm = 3.08568e+24;
  run_globals.units.UnitMass_in_g = 3.08568e+24;
  run_globals.units.UnitVelocity_in_cm_per_s = 3.08568e+24;

  set_units();

  run_globals.LTTime = malloc(sizeof(double) * n_snaps);
  run_globals.ZZ = malloc(sizeof(double) * n_snaps);
  char* head = alist;
  for (int ii = 0; ii < n_snaps + 1; ii++) {
    int seek = 0;
    double expfac = 0;
    if (sscanf(head, "%lf%n", &expfac, &seek)) {
      run_globals.ZZ[ii] = 1 / expfac - 1;
      run_globals.LTTime[ii] = time_to_present(run_globals.ZZ[ii]);
      head += seek + 1;
    }
  }
}

void teardown(void)
{
  free(run_globals.ZZ);
  free(run_globals.LTTime);
}

// THESE TESTS ARE CURRENTLY DEFUNCT AND NEED UPDATED
ParameterizedTestParameters(init, find_min_dt)
{
  int n_history_snaps[] = { 9, 10, 11, 12 };
  return cr_make_param_array(int, n_history_snaps, 4);
}

ParameterizedTest(int* n_history_snaps, init, find_min_dt, .init = setup, .fini = teardown)
{
  double min_dt = find_min_dt(*n_history_snaps);
  if (*n_history_snaps >= 11)
    cr_expect(min_dt >= 40.0);
  else
    cr_expect(min_dt < 40.0);
}