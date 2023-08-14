#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>

// C++ includes
#include <cmath>
#include <array>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <tuple>
#include <optional>

// Not implemented error
class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};

// Usings
using f64 = double;
using std::array;
using std::vector;
using std::string;
using std::tuple;
using std::optional;

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
enum SimMode {
    KINETIC_REACTIONS = 0,
    TRANSPORT_ONLY = 1,
};

// RT primary species types
enum PrimarySpeciesType {
    AQUEOUS = 1,
    ADSORPTION = 2,
    CATION_ECHG = 3,
    MINERAL = 4,
    // SECONDARY = 5,
    SPECIES_TYPE_ERROR = -1
};

// RT mass action types
enum MassActionType {
    IMMOBILE = 0,
    MOBILE = 1,
    MIXED = 2,
    ERROR = -1
};

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
        SimMode use_activity;           // activity coefficient mode: 0 = unity coefficient, 1 = DH equation
        int             transport_only;         // transport only flag: 0 = simulate kinetic reaction, 1 = transport only
        int             variable_precipchem;    // precipitation chemistry mode: 0 = constant precipitation chemistry, 1 = time-series precipitation chemistry   2021-05-20
        int             precipchem_numexp;      // Numerical experiment mode: 0 = same precipitation chemistry during warm-up and simulation
                                                // 1 = different precipitation chemistry during warm-up and simulation useful for numerical experiment  2021-09-09
        
        // Methods
        ControlData();
        ControlData(const int recycle, const SimMode use_activity, const int transport_only, const int variable_precipchem, const int precipchem_numexp) :
            recycle(recycle), use_activity(use_activity), transport_only(transport_only), variable_precipchem(variable_precipchem), precipchem_numexp(precipchem_numexp) { };
        ControlData copy();
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
        array<array<f64, MAXSPS>, MAXSPS> dep_mtx; // dependency of secondary species on primary species
        array<array<f64, MAXSPS>, MAXSPS> dep_kin;// dependency of kinetic species on primary species
        array<array<f64, MAXSPS>, MAXSPS> conc_contrib;   // contribution of each species to total concentration
        array<f64, MAXSPS> keq;
        array<f64, MAXSPS> keq_kin;
        double          adh;                    // Debye Huckel parameter
        double          bdh;                    // Debye Huckel parameter
        double          bdt;                    // Debye Huckel parameter

        // Methods
        ReactionNetwork();
        // ReactionNetwork copy();
};

class ChemTableEntry
{
    public:
        string          name;
        double          molar_mass;             // (g mol-1)
        double          molar_vol;              // (cm3 mol-1)
        double          charge;                 // charge
        double          size_fac;               // size factor for DH equation
        PrimarySpeciesType     itype;                  // type of primary species
                                                // 1 = primary aqueous, 2 = primary adsorption,
                                                // 3 = primary cation exchange, 4 = primary mineral
        MassActionType  mtype;

        // Methods
        ChemTableEntry();
        ChemTableEntry(const string& name, const f64 molar_mass, const f64 molar_vol, const f64 charge, const f64 size_fac, const PrimarySpeciesType itype, const MassActionType mtype);
        ChemTableEntry copy();
};

class KineticTableEntry
{
    public:
        int             position;               // position of target mineral in the array of primary species
        string          label;       // label of kinetic reaction
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

        // Methods
        // KineticTableEntry();
        // KineticTableEntry copy();
};

class SoilParameters {
    public:
        array<f64, MAXSPS> ssa;
        array<f64, MAXSPS> sw_thld;
        array<f64, MAXSPS> sw_exp;
        array<f64, MAXSPS> q10;
        array<f64, MAXSPS> n_alpha;

        // Methods
        SoilParameters();
        SoilParameters(const array<f64, MAXSPS> ssa, const array<f64, MAXSPS>  sw_thld, const array<f64, MAXSPS> sw_exp, const array<f64, MAXSPS> q10, const array<f64, MAXSPS> n_alpha) :
            ssa(ssa), sw_thld(sw_thld), sw_exp(sw_exp), q10(q10), n_alpha(n_alpha) { };
        SoilParameters copy() const;
};

class ChemicalState
{
    public:
        array<f64, MAXSPS> tot_conc;       // concentration (mol kgH2O-1)
        array<f64, MAXSPS> prim_conc;      // primary concentration (mol kgH2O-1)
        array<f64, MAXSPS> sec_conc;       // secondary concentration (mol kgH2O-1)
        array<f64, MAXSPS> prim_actv;      // activity of primary species
        array<f64, MAXSPS> tot_mol;        // total moles (mol m-2)
        SoilParameters soil_parameters;

        // Methods
        ChemicalState();
        ChemicalState(const array<f64, MAXSPS> tot_conc, const array<f64, MAXSPS> prim_conc, const array<f64, MAXSPS> sec_conc, const array<f64, MAXSPS> prim_actv, const array<f64, MAXSPS> tot_mol, SoilParameters params) : 
            tot_conc(tot_conc), prim_conc(prim_conc), sec_conc(sec_conc), prim_actv(prim_actv), tot_mol(tot_mol), soil_parameters(params.copy()) { };
        ChemicalState copy() const;
};

class CalibrationStruct {
    public:
        double          xsorption;
        double          rate;
        double          ssa;

        // Methods
        CalibrationStruct();
        CalibrationStruct(const f64 xsorption, const f64 rate, const f64 ssa) :
            xsorption(xsorption), rate(rate), ssa(ssa) { };
        CalibrationStruct copy();
};

class SoilConstants {
    public:
        f64 porosity;
        f64 depth;
        f64 ws_passive;

        // Methods
        SoilConstants();
        SoilConstants(const f64 porosity, const f64 depth, const f64 ws_passive) : porosity(porosity), depth(depth), ws_passive(ws_passive) {};
        SoilConstants copy();
};

class HBVParameters {
    public:
        f64 k1;
        f64 k2;
        f64 maxbas;
        f64 perc;
        f64 sfcf;
        f64 tt;

        // Methods
        HBVParameters();
        HBVParameters(const f64 k1, const f64 k2, const f64 maxbas, const f64 perc, const f64 sfcf, const f64 tt) :
            k1(k1), k2(k2), maxbas(maxbas), perc(perc), sfcf(sfcf), tt(tt) { };
        HBVParameters copy();
};

class Subcatchment
{
    public:
        vector<array<f64, NWS>> ws;
        vector<array<f64, NQ>> q;
        vector<vector<f64>> prcp_conc_time;     // time-series precipitation concentration (mol L-1)  2021-05-20
        vector<f64> tmp;                        // air temperature (degree C)
        array<f64, MAXSPS> prcp_conc;           // concentration in precipitation (mol kgH2O-1)
        SoilConstants   soil_surface;           // Surface soil parameters
        SoilConstants   soil_sz;                // Shallow zone soil parameters
        SoilConstants   soil_dz;                // Deep zone soil parameters
        HBVParameters hbv_parameters;
        array<array<f64, MAXSPS>, NWS> react_rate;  // reaction rate (mol m-2 day-1)
        array<ChemicalState, NWS> chms;
        ChemicalState river_chms;

        // Methods
        // Subcatchment copy();
};


#define fopen                   _custom_fopen
#define biort_printf(...)       _custom_printf(verbose_mode, __VA_ARGS__)
#if defined(_WIN32) || defined(_WIN64)
# define mkdir(path)            _mkdir((path))
#else
# define mkdir(path)            mkdir(path, 0755)
#endif

void            CopyConstSubcatchProp(const Subcatchment& subcatch, Subcatchment& subcatch_numexp);
void            CopyInitChemSubcatch(ReactionNetwork& rttbl, const Subcatchment& subcatch, Subcatchment& subcatch_numexp);
int             CountLeapYears(int, int);
int             FindChem(const char [MAXSTRING], int, const array<ChemTableEntry, MAXSPS>& chemtbl);
int             GetDifference(int, int);
void            InitChem(const string& input_dir, const CalibrationStruct& calib, const ControlData ctrl, array<ChemTableEntry, MAXSPS>& chemtbl,
    array<KineticTableEntry, MAXSPS>& kintbl, ReactionNetwork& rttbl, Subcatchment& subcatch);
void            InitChemState(double, double, const array<ChemTableEntry, MAXSPS>& chemtbl, const ReactionNetwork& rttbl, const ControlData ctrl,
    ChemicalState& chms);
void            Lookup(FILE *, const CalibrationStruct& calib, array<ChemTableEntry, MAXSPS>& chemtbl, array<KineticTableEntry, MAXSPS>& kintbl, ReactionNetwork& rttbl);
int             MatchWrappedKey(const char [], const char []);
string          ParseCmdLineParam(int argc, char *argv[]);
void            ParseLine(const char [], char [], double *);
void            PrintDailyResults(FILE *, int, int, const ReactionNetwork& rttbl, const Subcatchment& subcatch);
void            PrintHeader(FILE *, int, const ReactionNetwork& rttbl, const array<ChemTableEntry, MAXSPS>& chemtbl);
double          ReactControl(const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl, double, double, double,
    double, double, double, array<f64, MAXSPS>&, ChemicalState&);
void            Reaction(int, double, const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl,
    const ReactionNetwork& rttbl, Subcatchment& subcatch);
void            ReadAdsorption(const char [], int, int, array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork& rttbl);
void            ReadCationEchg(const char [], double, array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork& rttbl);
void            ReadChem(const string& input_dir, ControlData& ctrl, ReactionNetwork& rttbl, array<ChemTableEntry, MAXSPS>& chemtbl, array<KineticTableEntry, MAXSPS>& kintbl);
void            ReadCini(const string& input_dir, const array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork& rttbl, Subcatchment& subcatch);
void            ReadConc(FILE *, int, const array<ChemTableEntry, MAXSPS>& chemtbl, int *, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&, array<f64, MAXSPS>&);
void            ReadDHParam(const char [], int, double *);
void            ReadHbvParam(const string& input_dir, Subcatchment& subcatch);
int             ReadHbvResults(const string& input_dir, vector<int>& steps, Subcatchment& subcatch, int);
vector<int>     ReadPrecipChem(const string& input_dir, Subcatchment& subcatch, int, const array<ChemTableEntry, MAXSPS>& chemtbl, int);
void            ReadMinerals(const char [], int, int, double [MAXSPS][MAXSPS], double [], array<ChemTableEntry, MAXSPS>& chemtbl,
    ReactionNetwork& rttbl);
void            ReadMinKin(FILE *, int, double, int *, char [], array<ChemTableEntry, MAXSPS>& chemtbl, KineticTableEntry& kintbl);
void            ReadPrimary(const char [], int, array<ChemTableEntry, MAXSPS>& chemtbl);
void            ReadSecondary(const char [], int, int, array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork& rttbl);
void            ReadSoil(const string& input_dir, Subcatchment& sub);
void            ReadTempPoints(const char [], double, int *, int *);
int             roundi(double);
double          SoilTempFactor(double, double);
double          SoilMoistFactor(double, double, double);
optional<ChemicalState>   SolveReact(const double, const array<ChemTableEntry, MAXSPS>& chemtbl, const KineticTableEntry [], const ReactionNetwork& rttbl, const double, const double,
    const double, const double, const ChemicalState& chms);
int             SolveSpeciation(const array<ChemTableEntry, MAXSPS>& chemtbl, const ControlData ctrl, const ReactionNetwork& rttbl, int, ChemicalState& chms);
void            SortChem(char [MAXSPS][MAXSTRING], const int [MAXSPS], int, array<ChemTableEntry, MAXSPS>& chemtbl);
void            Speciation(const array<ChemTableEntry, MAXSPS>& chemtbl, const ControlData ctrl, const ReactionNetwork& rttbl, Subcatchment& subcatch);
int             SpeciesType(const char [], const char []);
void            StreamSpeciation(int, const array<ChemTableEntry, MAXSPS>& chemtbl, const ControlData ctrl, const ReactionNetwork& rttbl,
    Subcatchment& subcatch);
void            Transport(int step, const array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork& rttbl, const ControlData ctrl, Subcatchment& subcatch);   // 2021-05-21
void            WrapInParentheses(char *str);
double          WTDepthFactor(double ,double );
void            UnwrapParentheses(const char wrapped_str[], char str[]);
void            UpdatePrimConc(const ReactionNetwork& rttbl, const ControlData ctrl, Subcatchment& subcatch);

// Andrew's functions
bool            CheckArrayForNan(const array<f64, MAXSPS>& arr);
void            CheckChmsForNonFinite(const ChemicalState& chms, const char* filename, const int line_number);
bool            CheckNonzeroRanged(const array<f64, MAXSPS>& arr, const int num);
void            ErrOnZeroRanged(const char* filename, const char* arr_name, const int line_number, const array<f64, MAXSPS>& arr, const int num);
void            ComputeDependence(array<f64, MAXSPS>& tmpconc, const array<array<f64, MAXSPS>, MAXSPS>& dep_mtx, const array<f64, MAXSPS>& keq, int num_rows, int num_cols, int offset);
void            ErrorOnArrayNan(const array<f64, MAXSPS>& arr, const char* array_name, const char* filename, const int line_number);
void            ErrorOnArrayNanIter(const array<f64, MAXSPS>& arr, const char* array_name, const char* filename, const int line_number, const int iter);
void            GetLogActivity(array<f64, MAXSPS>& tmpconc, double gamma[MAXSPS], const array<array<f64, MAXSPS>, MAXSPS>& dep_mtx, const array<f64, MAXSPS>& keq, int num_rows, int num_cols, int offset);
void            GetSurfaceAreaRange(array<f64, MAXSPS>& area, const array<f64, MAXSPS>& prim_conc, const array<f64, MAXSPS>& ssa, const array<ChemTableEntry, MAXSPS>& chemtbl, int start, int end, int offset);
void            GetTempFactorRange(array<f64, MAXSPS>& ftemp, const array<f64, MAXSPS>& q10, double temperature, int start, int end, int offset);
void            GetWTDepthFactorRange(array<f64, MAXSPS>& fzw, double Zw, const array<f64, MAXSPS>& n_alpha, int start, int end, int offset);
void            Log10Arr(const array<f64, MAXSPS>& src, array<f64, MAXSPS>& dst, int num_species);
void            Pow10Arr(const double src[MAXSPS], double dst[MAXSPS], int num_species);
void            PrintArray(const array<f64, MAXSPS>& arr);
void            PrintMatrix(const realtype** mat, const int nrows, const int ncols);
void            SetZero(double arr[MAXSPS]);
void            SetZeroRange(double arr[MAXSPS], int start, int end);
void            SoilMoistFactorRange(array<f64, MAXSPS>& dst, double satn, const array<f64, MAXSPS>& sw_threshold, const array<f64, MAXSPS>& sw_exponent, 
                    int start, int end, int offset);
double          SumArr(const double arr[MAXSPS], int num_species);

void            TransportDeepZone(const ReactionNetwork& rttbl, const array<ChemTableEntry, MAXSPS>& chemtbl, Subcatchment& subcatch, const int step);
void            TransportShallowZone(const ReactionNetwork& rttbl, const array<ChemTableEntry, MAXSPS>& chemtbl, Subcatchment& subcatch, const int step);
void            TransportSurfaceZone(const ReactionNetwork& rttbl, const ControlData ctrl, Subcatchment& subcatch, const int step);

void            PrintChemicalState(const ChemicalState& chms);

double          ReadParamToDouble(const char buffer[], const char keyword[], const char fn[], int line_number);
int             ReadParamToInt(const char buffer[], const char keyword[], const char fn[], int line_number);
void            ReactSurfaceZone(const double temp, const SoilConstants soil, const double tot_water, const array<ChemTableEntry, MAXSPS>& chemtbl,
                    const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl, double stepsize, Subcatchment& subcatch);
optional<ChemicalState> SolveSurfaceReact(const double stepsize, const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl,
                    const double tot_water, const double temp, const double porosity, const ChemicalState& chms);
double          ReactSurfaceControl(const array<ChemTableEntry, MAXSPS>& chemtbl, const array<KineticTableEntry, MAXSPS>& kintbl, const ReactionNetwork& rttbl,
                    double stepsize, double porosity, double depth, double satn, double temp, array<f64, MAXSPS>& react_rate,
                    ChemicalState&chms);

//===== Constant definitions =====//
const string CHEM_FILE_DIR = "input/%s/chem.txt";
const string CHEM_RECYCLE_ID = "RECYCLE";
const string CHEM_ACTIVITY_ID = "ACTIVITY";
const string CHEM_TRANSPORT_ONLY_ID = "TRANSPORT_ONLY";
const string CHEM_PRECIPCHEM_ID = "PRECIPCHEM";
const string CHEM_NUMEXP_ID = "NUMEXP";
const string CHEM_TEMPERATURE_ID = "TEMPERATURE";

const string CHEM_PRIMARY_SPECIES_ID = "PRIMARY_SPECIES";
const string CHEM_SECONDARY_SPECIES_ID = "SECONDARY_SPECIES";
const string CHEM_MINERAL_KINETICS_ID = "MINERAL_KINETICS";

const double SATN_MINIMUM = 1e-4;       // Minimum saturation to initiate kinetic reactions
const int MAX_ITERATIONS = 25;
const double MINIMUM_SUBSTEP = 1e-20;