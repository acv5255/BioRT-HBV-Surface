/*
    Input files for BioRT-HBV
*/
#include <string>
#include <map>
#include <vector>
#include <filesystem>

using f64 = double;
using std::string;
using std::map;
using std::vector;
using std::filesystem::path;

class CiniFile {
    public:
        // Class structs
        struct MineralCondition {
            f64 vol_frac;
            f64 ssa;
            f64 q10;
            f64 sw_thr;
            f64 sw_exp;
            f64 n_alpha;

            MineralCondition from_tokens(const string& line);
        };

        map<string, f64> aqueous;
        map<string, MineralCondition> mineral;

        // Methods
        static CiniFile FromFile(const path& file_path);
};

class ChemFile {
    struct  RunParameters {
        int recycle;
        bool use_activity;
        bool transport_only;
        bool variable_precipchem;
        bool do_numexp;
        f64 temperature;            // Temperature of the simluation (in celsius)

        RunParameters(int recycle, bool use_activity, bool transport_only, bool variable_precipchem, bool do_numexp, f64 temperature) :
            recycle(recycle), use_activity(use_activity), transport_only(transport_only), variable_precipchem(variable_precipchem), temperature(temperature) { };
    };


    public:
        vector<string> primary_species;
        vector<string> secondary_species;
        vector<string> cation_echg_species;
        map<string, string> mineral_kinetic_entries;

        static ChemFile FromFile(const path& chemfile_path);
};

class SoilFile {
    public:
        struct SoilParameters {
            f64 porosity;       // Porosity [0,1]
            f64 ws_passive;     // Passive water storage [mm*m^2/m^2]
            f64 depth;          // Depth of zone [mm]
        };

        SoilFile FromFile(const path& soilfile_path);

};

class PrecipChem {
    public:
        map<string, vector<f64>> data;

        static PrecipChem FromFile(const path& precipchem_path);
};

class HBVResults {

};

class HBVParameters {

};

class Database {

};