/*
    Input files for BioRT-HBV
*/
#include <string>
#include <map>

using f64 = double;
using std::string;
using std::map;

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
        static CiniFile FromFile(const string& file_path);
};

class ChemFile {

};

class SoilFile {
    public:
        struct SoilParameters {
            f64 porosity;       // Porosity [0,1]
            f64 ws_passive;     // Passive water storage [mm*m^2/m^2]
            f64 depth;          // Depth of zone [mm]
        };

};

class PrecipChem {

};

class HBVResults {

};

class HBVParameters {

};

class Database {

};