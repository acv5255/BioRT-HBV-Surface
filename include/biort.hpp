#pragma once

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

// C++ includes
#include <cmath>
#include <array>

// Usings
using f64 = double;
using std::array;

// SUNDIAL Header Files
#include "sundials/sundials_dense.h"        // Prototypes for small dense fcts.

#include "custom_io.hpp"

#define VERSION    "0.1.0-alpha"

const double BADVAL = -999;
const int NWS = 6;
const int NQ = 9; 

const double STORAGE_MIN = 1.0;
const double ZERO_CONC = 1.0E-20;

const double TOLERANCE = 1E-8;       // Tolerance for iterative methods

const int PRECIP = 0;
const int RECHG = 1;
const int PERC = 2;
const int Q0 = 3;
const int Q1 = 4;
const int Q2 = 5;
const int Prain = 6;
const int Psnow = 7;
const int snowmelt = 8;

const int SNOW = 0;
const int SURFACE = 1;
const int UZ = 2;
const int LZ = 3;
const int STREAM = 4;
const int SM = 5;

const int MAXSPS = 20;
const int MAXDEP = 4;
const int SKIP_JACOB = 1;

// RT simulation mode
const int KIN_REACTION = 0;
const int TRANSPORT_ONLY = 1;

// RT primary species types
const int AQUEOUS = 1;
const int ADSORPTION = 2;
const int CATION_ECHG = 3;
const int MINERAL = 4;
const int SECONDARY = 5;

// RT mass action types
const int IMMOBILE_MA = 0;
const int MOBILE_MA = 1;
const int MIXED_MA = 2;

// RT kinetic reaction types
const int TST = 1;
const int PRCP_ONLY = 2;
const int DISS_ONLY = 3;
const int MONOD = 4;

extern int     verbose_mode;


class ControlData
{
    public:
        int             recycle;                // number of times to recycle forcing
        int             read_restart;           // flag to read rt restart file
        int             use_activity;           // activity coefficient mode: 0 = unity coefficient, 1 = DH equation
        int             transport_only;         // transport only flag: 0 = simulate kinetic reaction, 1 = transport only
        int             variable_precipchem;    // precipitation chemistry mode: 0 = constant precipitation chemistry, 1 = time-series precipitation chemistry   2021-05-20
        int             precipchem_numexp;      // Numerical experiment mode: 0 = same precipitation chemistry during warm-up and simulation
                                                // 1 = different precipitation chemistry during warm-up and simulation useful for numerical experiment  2021-09-09
};

class ReactionNetwork
{
    public:
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
        array<f64, MAXSPS> keq;
        array<f64, MAXSPS> keq_kin;
        double          adh;                    // Debye Huckel parameter
        double          bdh;                    // Debye Huckel parameter
        double          bdt;                    // Debye Huckel parameter
};

class ChemTableEntry
{
    public:
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
};

class KineticTableEntry
{
    public:
        int             position;               // position of target mineral in the array of primary species
        char            label[MAXSTRING];       // label of kinetic reaction
        int             type;                   // type of the kinetic reaction
                                                // 1: tst, 2: precipitation-only, 3: dissolution-only, 4: monod
        double          rate;                   // rate of kinetic reaction
        double          actv;                   // activation energy, used to calculate kinetic rate of reaction under
                                                // different temperatures
        double          keq;                    // equilibrium constant
        int             ndep;                   // number of dependency
        array<int, MAXDEP> dep_index;           // position of species that kinetic reaction depends on
        array<f64, MAXDEP> dep_power;           // power of dependency
        int             biomass_index;          // position of biomass species
        int             nmonod;                 // number of monod species
        array<int, MAXDEP> monod_index;         // position of monod species
        array<f64, MAXDEP> monod_para;          // parameter for monod dependency
        int             ninhib;                 // number of inhibition species
        array<int, MAXDEP> inhib_index;         // position of inhibition species
        array<f64, MAXDEP> inhib_para;          // parameters that controls this inhibition
};

class ChemicalState
{
    public:
        array<f64, MAXSPS> tot_conc;       // concentration (mol kgH2O-1)
        array<f64, MAXSPS> prim_conc;      // primary concentration (mol kgH2O-1)
        array<f64, MAXSPS> sec_conc;       // secondary concentration (mol kgH2O-1)
        array<f64, MAXSPS> prim_actv;      // activity of primary species
        array<f64, MAXSPS> ssa;            // specific surface area (m2 g-1)
        array<f64, MAXSPS> sw_thld;        // threshold in soil moisture function (-)
        array<f64, MAXSPS> sw_exp;         // exponent in soil moisture function (-)
        array<f64, MAXSPS> q10;            // Q10 factor (-)
        array<f64, MAXSPS> n_alpha;        // n*alpha in depth function (-)
        array<f64, MAXSPS> tot_mol;        // total moles (mol m-2)
};

class CalibrationStruct
{
    public:
        double          xsorption;
        double          rate;
        double          ssa;
};

class SoilConstants {
    public:
        double porosity;
        double depth;
        double ws_passive;
};

class Subcatchment
{
    public:
        double         (*ws)[NWS];
        double         (*q)[NQ];
        double        **prcp_conc_time;         // time-series precipitation concentration (mol L-1)  2021-05-20
        double         *tmp;                    // air temperature (degree C)
        array<f64, MAXSPS> prcp_conc;           // concentration in precipitation (mol kgH2O-1)
        SoilConstants   soil_surface;           // Surface soil parameters
        SoilConstants   soil_sz;                // Shallow zone soil parameters
        SoilConstants   soil_dz;                // Deep zone soil parameters
        double          k1;                     // recession coefficient for upper zone (day -1)
        double          k2;                     // recession coefficient for lower zone (day -1)
        double          maxbas;                 // routing parameter
        double          perc;                   // percolation rate (mm day-1)
        double          sfcf;                   // snow fall correction factor
        double          tt;                     // threshold temperature (degree C)
        // double          react_rate[NWS][MAXSPS];
        array<array<f64, MAXSPS>, NWS> react_rate;  // reaction rate (mol m-2 day-1)
        ChemicalState chms[NWS];
        ChemicalState river_chms;
};


#define fopen                   _custom_fopen
#define biort_printf(...)       _custom_printf(verbose_mode, __VA_ARGS__)
#if defined(_WIN32) || defined(_WIN64)
# define mkdir(path)            _mkdir((path))
#else
# define mkdir(path)            mkdir(path, 0755)
#endif

void            CopyConstSubcatchProp(const Subcatchment* subcatch, Subcatchment* subcatch_numexp);
void            CopyInitChemSubcatch(ReactionNetwork *, const Subcatchment* subcatch, Subcatchment* subcatch_numexp );
int             CountLeapYears(int, int);
int             FindChem(const char [MAXSTRING], int, const ChemTableEntry[]);
void            FreeStruct(int *steps[], Subcatchment subcatch[]);
int             GetDifference(int, int);
void            InitChem(const char [], const CalibrationStruct *, const ControlData ctrl, ChemTableEntry [],
    KineticTableEntry [], ReactionNetwork *, Subcatchment* subcatch);
void            InitChemState(double, double, const ChemTableEntry [], const ReactionNetwork *, const ControlData ctrl,
    ChemicalState *);
void            Lookup(FILE *, const CalibrationStruct *, ChemTableEntry [], KineticTableEntry [], ReactionNetwork *);
int             MatchWrappedKey(const char [], const char []);
void            ParseCmdLineParam(int argc, char *argv[], char dir[]);
void            ParseLine(const char [], char [], double *);
void            PrintDailyResults(FILE *, int, int, const ReactionNetwork *, const Subcatchment* subcatch);
void            PrintHeader(FILE *, int, const ReactionNetwork *, const ChemTableEntry chemtbl[]);
double          ReactControl(const ChemTableEntry [], const KineticTableEntry [], const ReactionNetwork *, double, double, double,
    double, double, double, array<f64, MAXSPS>&, ChemicalState *);
void            Reaction(int, double, const ChemTableEntry [], const KineticTableEntry [],
    const ReactionNetwork *, Subcatchment* subcatch);
void            ReadAdsorption(const char [], int, int, ChemTableEntry [], ReactionNetwork *);
void            ReadCationEchg(const char [], double, ChemTableEntry [], ReactionNetwork *);
void            ReadChem(const char [], ControlData *, ReactionNetwork *, ChemTableEntry [], KineticTableEntry []);
void            ReadCini(const char [], const ChemTableEntry *, ReactionNetwork *, Subcatchment* subcatch);
void            ReadConc(FILE *, int, const ChemTableEntry [], int *, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&);
void            ReadDHParam(const char [], int, double *);
void            ReadHbvParam(const char [], Subcatchment* subcatch);
void            ReadHbvResults(const char [], int *, int **, Subcatchment* subcatch, int);
void            ReadPrecipChem(const char [], int *, int **, Subcatchment* subcatch, int, const ChemTableEntry [], int);
void            ReadMinerals(const char [], int, int, double [MAXSPS][MAXSPS], double [], ChemTableEntry [],
    ReactionNetwork *);
void            ReadMinKin(FILE *, int, double, int *, char [], ChemTableEntry [], KineticTableEntry *);
void            ReadPrimary(const char [], int, ChemTableEntry []);
void            ReadSecondary(const char [], int, int, ChemTableEntry [], ReactionNetwork *);
void            ReadSoil(const char [], Subcatchment* sub);
void            ReadTempPoints(const char [], double, int *, int *);
int             roundi(double);
double          SoilTempFactor(double, double);
double          SoilMoistFactor(double, double, double);
int             SolveReact(double, const ChemTableEntry [], const KineticTableEntry [], const ReactionNetwork *, double, double,
    double, double, ChemicalState *);
int             SolveSpeciation(const ChemTableEntry [], const ControlData ctrl, const ReactionNetwork *, int, ChemicalState *);
void            SortChem(char [MAXSPS][MAXSTRING], const int [MAXSPS], int, ChemTableEntry []);
void            Speciation(const ChemTableEntry [], const ControlData ctrl, const ReactionNetwork *, Subcatchment* subcatch);
int             SpeciesType(const char [], const char []);
void            StreamSpeciation(int, const ChemTableEntry [], const ControlData ctrl, const ReactionNetwork *,
    Subcatchment* subcatch);
void            Transport(int step, const ChemTableEntry *chemtbl, ReactionNetwork *rttbl, const ControlData ctrl, Subcatchment* subcatch);   // 2021-05-21
void            WrapInParentheses(char *str);
double          WTDepthFactor(double ,double );
void            UnwrapParentheses(const char wrapped_str[], char str[]);
void            UpdatePrimConc(const ReactionNetwork *rttbl, const ControlData ctrl, Subcatchment* subcatch);

// Andrew's functions
bool            CheckArrayForNan(const array<f64, MAXSPS>& arr);
void            CheckChmsForNonFinite(const ChemicalState* chms, const char* filename, const int line_number);
bool            CheckNonzeroRanged(const array<f64, MAXSPS>& arr, const int num);
void            ErrOnZeroRanged(const char* filename, const char* arr_name, const int line_number, const array<f64, MAXSPS>& arr, const int num);
void            ComputeDependence(array<f64, MAXSPS>& tmpconc, const double dep_mtx[MAXSPS][MAXSPS], const array<f64, MAXSPS>& keq, int num_rows, int num_cols, int offset);
void            ErrorOnArrayNan(const array<f64, MAXSPS>& arr, const char* array_name, const char* filename, const int line_number);
void            ErrorOnArrayNanIter(const array<f64, MAXSPS>& arr, const char* array_name, const char* filename, const int line_number, const int iter);
void            GetLogActivity(array<f64, MAXSPS>& tmpconc, double gamma[MAXSPS], const double dep_mtx[MAXSPS][MAXSPS], const array<f64, MAXSPS>& keq, int num_rows, int num_cols, int offset);
void            GetSecondarySpecies(double conc[MAXSPS], const double gamma[MAXSPS], const ReactionNetwork* rttbl);
void            GetSurfaceAreaRange(double area[MAXSPS], const array<f64, MAXSPS>& prim_conc, const array<f64, MAXSPS>& ssa, const ChemTableEntry chemtbl[], int start, int end, int offset);
void            GetTempFactorRange(double ftemp[MAXSPS], const array<f64, MAXSPS>& q10, double temperature, int start, int end, int offset);
void            GetWTDepthFactorRange(double fzw[MAXSPS], double Zw, const array<f64, MAXSPS>& n_alpha, int start, int end, int offset);
void            Log10Arr(const array<f64, MAXSPS>& src, array<f64, MAXSPS>& dst, int num_species);
void            Pow10Arr(const double src[MAXSPS], double dst[MAXSPS], int num_species);
void            PrintArray(const array<f64, MAXSPS>& arr);
void            PrintMatrix(const realtype** mat, const int nrows, const int ncols);
void            SetZero(double arr[MAXSPS]);
void            SetZeroRange(double arr[MAXSPS], int start, int end);
void            SoilMoistFactorRange(double dst[MAXSPS], double satn, const array<f64, MAXSPS>& sw_threshold, const array<f64, MAXSPS>& sw_exponent, 
                    int start, int end, int offset);
double          SumArr(const double arr[MAXSPS], int num_species);

void            TransportDeepZone(const ReactionNetwork* rttbl, const ChemTableEntry chemtbl[], Subcatchment* subcatch, const int step);
void            TransportShallowZone(const ReactionNetwork* rttbl, const ChemTableEntry chemtbl[], Subcatchment* subcatch, const int step);
void            TransportSurfaceZone(const ReactionNetwork* rttbl, const ControlData ctrl, Subcatchment* subcatch, const int step);

void            PrintChemicalState(const ChemicalState* chms);

double          ReadParamToDouble(const char buffer[], const char keyword[], const char fn[], int line_number);
void            ResetReactionRates(Subcatchment * subcatch);
int             ReadParamToInt(const char buffer[], const char keyword[], const char fn[], int line_number);
void            ReactSurfaceZone(const double temp, const SoilConstants soil, const double tot_water, const ChemTableEntry chemtbl[],
                    const KineticTableEntry kintbl[], const ReactionNetwork* rttbl, double stepsize, Subcatchment* subcatch);
int             SolveSurfaceReact(double stepsize, const ChemTableEntry chemtbl[], const KineticTableEntry kintbl[], const ReactionNetwork *rttbl,
                    double tot_water, double temp, double porosity, ChemicalState *chms);
double          ReactSurfaceControl(const ChemTableEntry chemtbl[], const KineticTableEntry kintbl[], const ReactionNetwork *rttbl,
                    double stepsize, double porosity, double depth, double satn, double temp, array<f64, MAXSPS>& react_rate,
                    ChemicalState *chms);

//===== Constant definitions =====//
static const char CHEM_FILE_DIR[] = "input/%s/chem.txt";
static const char CHEM_RECYCLE_ID[] = "RECYCLE";
static const char CHEM_ACTIVITY_ID[] = "ACTIVITY";
static const char CHEM_TRANSPORT_ONLY_ID[] = "TRANSPORT_ONLY";
static const char CHEM_PRECIPCHEM_ID[] = "PRECIPCHEM";
static const char CHEM_NUMEXP_ID[] = "NUMEXP";
static const char CHEM_TEMPERATURE_ID[] = "TEMPERATURE";

static const char CHEM_PRIMARY_SPECIES_ID[] = "PRIMARY_SPECIES";
static const char CHEM_SECONDARY_SPECIES_ID[] = "SECONDARY_SPECIES";
static const char CHEM_MINERAL_KINETICS_ID[] = "MINERAL_KINETICS";

static const double SATN_MINIMUM = 1e-4;       // Minimum saturation to initiate kinetic reactions
static const int MAX_ITERATIONS = 25;
static const double MINIMUM_SUBSTEP = 1e-20;