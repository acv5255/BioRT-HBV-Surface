#ifndef BIORT_HEADER
#define BIORT_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <stdarg.h>
#include <float.h>
#include <sys/stat.h>
#include <stdbool.h>

// SUNDIAL Header Files
#include "sundials/sundials_dense.h"        // Prototypes for small dense fcts.

#include "custom_io.h"

#define VERSION    "0.1.0-alpha"

#define BADVAL                  -999
#define NWS                     6           // number of water storages
#define NQ                      9           // number of water fluxes

#define STORAGE_MIN             1.0         // minimum water storage (mm)
#define ZERO_CONC               1.0E-20     // minimum concentration

#define TOLERANCE               1E-7        // tolerance for reaction calculation

static const int PRECIP = 0;
static const int RECHG = 1;
static const int PERC = 2;
static const int Q0 = 3;
static const int Q1 = 4;
static const int Q2 = 5;
static const int Prain = 6;
static const int Psnow = 7;
static const int snowmelt = 8;

#define SNOW                    0
#define SURFACE                 1          // add surface Q0 reaction, 2021-05-14
#define UZ                      2
#define LZ                      3
#define STREAM                  4
#define SM                      5
// static const int SNOW = 0;
// static const int SURFACE = 1;
// static const int UZ = 2;
// static const int LZ = 3;
// static const int STREAM = 4;
// static const int SM = 5;

#define MAXSPS 20
#define MAXDEP 4
#define SKIP_JACOB 1

// RT simulation mode
static const int KIN_REACTION = 0;
static const int TRANSPORT_ONLY = 1;

// RT primary species types
#define AQUEOUS                 1
#define ADSORPTION              2
#define CATION_ECHG             3
#define MINERAL                 4
#define SECONDARY               5

// RT mass action types
#define IMMOBILE_MA             0
#define MOBILE_MA               1
#define MIXED_MA                2

// RT kinetic reaction types
#define TST                     1
#define PRCP_ONLY               2
#define DISS_ONLY               3
#define MONOD                   4

extern int     verbose_mode;


typedef struct ctrl_struct
{
    int             recycle;                // number of times to recycle forcing
    int             read_restart;           // flag to read rt restart file
    int             use_activity;           // activity coefficient mode: 0 = unity coefficient, 1 = DH equation
    int             transport_only;         // transport only flag: 0 = simulate kinetic reaction, 1 = transport only
    int             variable_precipchem;    // precipitation chemistry mode: 0 = constant precipitation chemistry, 1 = time-series precipitation chemistry   2021-05-20
    int             precipchem_numexp;      // Numerical experiment mode: 0 = same precipitation chemistry during warm-up and simulation
                                            // 1 = different precipitation chemistry during warm-up and simulation useful for numerical experiment  2021-09-09
    double         *steps;                  // model steps
} ctrl_struct;

typedef struct rttbl_struct
{
    int             num_stc;                // number of total species
    int             num_spc;                // number of primary species
    int             num_ssc;                // number of secondary speices
    int             num_sdc;                // number of independent species
    int             num_min;                // number of minerals
    int             num_ads;                // number of adsorption species
    int             num_cex;                // number of cation exchange
    int             num_mkr;                // number of mineral kinetic reactions
    int             num_akr;                // number of aqueous kinetic reactions
    double          tmp;                    // temperature of the moment
    double          dep_mtx[MAXSPS][MAXSPS];// dependency of secondary species on primary species
    double          dep_kin[MAXSPS][MAXSPS];// dependency of kinetic species on primary species
    double          conc_contrib[MAXSPS][MAXSPS];   // contribution of each species to total concentration
    double          keq[MAXSPS];            // Keq's of secondary species
    double          keq_kin[MAXSPS];        // Keq's of kinetic species
    double          adh;                    // Debye Huckel parameter
    double          bdh;                    // Debye Huckel parameter
    double          bdt;                    // Debye Huckel parameter
} rttbl_struct;

typedef struct chemtbl_struct
{
    char            name[MAXSTRING];        // molecular formula or name
    double          molar_mass;             // (g mol-1)
    double          molar_vol;              // (cm3 mol-1)
    double          charge;                 // charge
    double          size_fac;               // size factor for DH equation
    int             itype;                  // type of primary species
                                            // 1 = primary aqueous, 2 = primary adsorption,
                                            // 3 = primary cation exchange, 4 = primary mineral
    int             mtype;                  // type of the mass action species
                                            // 0 = immobile mass action, 1 = mobile mass action,
                                            // 2 = mixed mobility mass action
} chemtbl_struct;

typedef struct kintbl_struct
{
    int             position;               // position of target mineral in the array of primary species
    char            label[MAXSTRING];       // label of kinetic reaction
    int             type;                   // type of the kinetic reaction
                                            // 1: tst, 2: precipitation-only, 3: dissolution-only, 4: monod
    double          rate;                   // rate of kinetic reaction
    double          actv;                   // activation energy, used to calculate kinetic rate of reaction under
                                            // different temperatures
    double          keq;                    // equilibrium constant
    int             ndep;                   // number of dependency
    int             dep_index[MAXDEP];      // position of species that kinetic reaction depends on
    double          dep_power[MAXDEP];      // power of dependency
    int             biomass_index;          // position of biomass species
    int             nmonod;                 // number of monod species
    int             monod_index[MAXDEP];    // position of monod species
    double          monod_para[MAXDEP];     // parameter for monod dependency
    int             ninhib;                 // number of inhibition species
    int             inhib_index[MAXDEP];    // position of inhibition species
    double          inhib_para[MAXDEP];     // parameters that controls this inhibition
} kintbl_struct;

typedef struct chmstate_struct
{
    double          tot_conc[MAXSPS];       // concentration (mol kgH2O-1)
    double          prim_conc[MAXSPS];      // primary concentration (mol kgH2O-1)
    double          sec_conc[MAXSPS];       // secondary concentration (mol kgH2O-1)
    double          prim_actv[MAXSPS];      // activity of primary species
    double          ssa[MAXSPS];            // specific surface area (m2 g-1)
    double          sw_thld[MAXSPS];        // threshold in soil moisture function (-)
    double          sw_exp[MAXSPS];         // exponent in soil moisture function (-)
    double          q10[MAXSPS];            // Q10 factor (-)
    double          n_alpha[MAXSPS];        // n*alpha in depth function (-)
    double          tot_mol[MAXSPS];        // total moles (mol m-2)
} chmstate_struct;

typedef struct calib_struct
{
    double          xsorption;
    double          rate;
    double          ssa;
} calib_struct;

typedef struct subcatch_struct
{
    double        **ws;                     // water storages (mm)
    double        **q;                      // water fluxes (mm day-1)
    double        **prcp_conc_time;         // time-series precipitation concentration (mol L-1)  2021-05-20
    double         *tmp;                    // air temperature (degree C)
    double          prcp_conc[MAXSPS];      // concentration in precipitation (mol kgH2O-1)
    //double          d_surface;              // surface zone maximum water (passive + dynamic) storage capacity (mm)
    double          d_uz;                   // upper zone maximum water (passive + dynamic) storage capacity (mm)
    double          d_lz;                   // lower zone maximum water (passive + dynamic) storage capacity (mm)
    double          k1;                     // recession coefficient for upper zone (day -1)
    double          k2;                     // recession coefficient for lower zone (day -1)
    double          maxbas;                 // routing parameter
    double          perc;                   // percolation rate (mm day-1)
    //double          porosity_surface;       // surface zone porosity (m3 m-3), 2021-05-14
    double          porosity_uz;            // upper zone porosity (m3 m-3)
    double          porosity_lz;            // lower zone porosity (m3 m-3)
    //double          res_surface;            // surface zone passive water storage (mm), 2021-05-14
    double          res_uz;                 // upper zone passive water storage (mm)
    double          res_lz;                 // lower zone passive water storage (mm)
    double          sfcf;                   // snow fall correction factor
    double          tt;                     // threshold temperature (degree C)
    double          react_rate[NWS][MAXSPS];// reaction rate (mol m-2 day-1)
    chmstate_struct chms[NWS];
    chmstate_struct river_chms;
} subcatch_struct;

#define MIN(x, y)               (((x) < (y)) ? (x) : (y))
#define MAX(x, y)               (((x) > (y)) ? (x) : (y))

#define fopen                   _custom_fopen
#define biort_printf(...)       _custom_printf(verbose_mode, __VA_ARGS__)
#if defined(_WIN32) || defined(_WIN64)
# define mkdir(path)            _mkdir((path))
#else
# define mkdir(path)            mkdir(path, 0755)
#endif

void            CopyConstSubcatchProp(const subcatch_struct* subcatch, subcatch_struct* subcatch_numexp);
void            CopyInitChemSubcatch(rttbl_struct *, const subcatch_struct* subcatch, subcatch_struct* subcatch_numexp );
int             CountLeapYears(int, int);
int             FindChem(const char [MAXSTRING], int, const chemtbl_struct[]);
void            FreeStruct(int nsteps, int *steps[], subcatch_struct subcatch[]);
int             GetDifference(int, int);
void            InitChem(const char [], const calib_struct *, const ctrl_struct *, chemtbl_struct [],
    kintbl_struct [], rttbl_struct *, subcatch_struct []);
void            InitChemState(double, double, const chemtbl_struct [], const rttbl_struct *, const ctrl_struct *,
    chmstate_struct *);
void            Lookup(FILE *, const calib_struct *, chemtbl_struct [], kintbl_struct [], rttbl_struct *);
int             MatchWrappedKey(const char [], const char []);
void            ParseCmdLineParam(int argc, char *argv[], char dir[]);
void            ParseLine(const char [], char [], double *);
void            PrintDailyResults(FILE *, int, int, const rttbl_struct *, const subcatch_struct* subcatch);
void            PrintHeader(FILE *, int, const rttbl_struct *, const chemtbl_struct chemtbl[]);
double          ReactControl(const chemtbl_struct [], const kintbl_struct [], const rttbl_struct *, double, double, double,
    double, double, double, double [], chmstate_struct *);
void            Reaction(int, double, const int [], const chemtbl_struct [], const kintbl_struct [],
    const rttbl_struct *, subcatch_struct* subcatch);
void            ReadAdsorption(const char [], int, int, chemtbl_struct [], rttbl_struct *);
void            ReadCationEchg(const char [], double, chemtbl_struct [], rttbl_struct *);
void            ReadChem(const char [], ctrl_struct *, rttbl_struct *, chemtbl_struct [], kintbl_struct []);
void            ReadCini(const char [], const chemtbl_struct *, rttbl_struct *, subcatch_struct* subcatch);
void            ReadConc(FILE *, int, const chemtbl_struct [], int *, double [], double [], double [], double [], double[], double[]);
void            ReadDHParam(const char [], int, double *);
void            ReadHbvParam(const char [], subcatch_struct* subcatch);
void            ReadHbvResults(const char [], int *, int **, subcatch_struct* subcatch, int);
void            ReadPrecipChem(const char [], int *, int **, subcatch_struct* subcatch, int, const chemtbl_struct [], int);
void            ReadMinerals(const char [], int, int, double [MAXSPS][MAXSPS], double [], chemtbl_struct [],
    rttbl_struct *);
void            ReadMinKin(FILE *, int, double, int *, char [], chemtbl_struct [], kintbl_struct *);
int             ReadParam(const char [], const char [], char, const char [], int, void *);
void            ReadPrimary(const char [], int, chemtbl_struct []);
void            ReadSecondary(const char [], int, int, chemtbl_struct [], rttbl_struct *);
void            ReadSoil(const char [], subcatch_struct* sub);
void            ReadTempPoints(const char [], double, int *, int *);
int             roundi(double);
//void            RunNumExp(int, int, const chemtbl_struct [], const kintbl_struct [], rttbl_struct *, const ctrl_struct *,
//    subcatch_struct [], FILE *);//2021-09-09
double          SoilTempFactor(double, double);
double          SoilMoistFactor(double, double, double);
int             SolveReact(double, const chemtbl_struct [], const kintbl_struct [], const rttbl_struct *, double, double,
    double, double, chmstate_struct *);
int             SolveSpeciation(const chemtbl_struct [], const ctrl_struct *, const rttbl_struct *, int, chmstate_struct *);
void            SortChem(char [MAXSPS][MAXSTRING], const int [MAXSPS], int, chemtbl_struct []);
void            Speciation(const chemtbl_struct [], const ctrl_struct *, const rttbl_struct *, subcatch_struct* subcatch);
int             SpeciesType(const char [], const char []);
void            StreamSpeciation(int, const chemtbl_struct [], const ctrl_struct *, const rttbl_struct *,
    subcatch_struct* subcatch);
void            Transpt(int step, const chemtbl_struct *chemtbl, rttbl_struct *rttbl, const ctrl_struct *ctrl, subcatch_struct* subcatch);   // 2021-05-21
void            WrapInParentheses(char *str);
double          WTDepthFactor(double ,double );
void            UnwrapParentheses(const char wrapped_str[], char str[]);
void            UpdatePrimConc(const rttbl_struct *rttbl, const ctrl_struct *ctrl, subcatch_struct* subcatch);

// Andrew's functions
void            ComputeDependence(double tmpconc[MAXSPS], const double dep_mtx[MAXSPS][MAXSPS], const double keq[MAXSPS], int num_rows, int num_cols, int offset);
void            GetLogActivity(double tmpconc[MAXSPS], double gamma[MAXSPS], const double dep_mtx[MAXSPS][MAXSPS], const double keq[MAXSPS], int num_rows, int num_cols, int offset);
void            Log10Arr(const double src[MAXSPS], double dst[MAXSPS], int num_species);
void            Pow10Arr(const double src[MAXSPS], double dst[MAXSPS], int num_species);
void            SetZero(double arr[MAXSPS]);
void            SetZeroRange(double arr[MAXSPS], int start, int end);
double          SumArr(const double arr[MAXSPS], int num_species);

void            SoilMoistFactorRange(double dst[MAXSPS], double satn, const double sw_threshold[MAXSPS], const double sw_exponent[MAXSPS], 
                    int start, int end, int offset);
void            GetSurfaceAreaRange(double area[MAXSPS], const double prim_conc[MAXSPS], const double ssa[MAXSPS], const chemtbl_struct chemtbl[], int start, int end, int offset);
void            GetTempFactorRange(double ftemp[MAXSPS], const double q10[MAXSPS], double temperature, int start, int end, int offset);
void            GetWTDepthFactorRange(double fzw[MAXSPS], double Zw, const double n_alpha[MAXSPS], int start, int end, int offset);
#endif
