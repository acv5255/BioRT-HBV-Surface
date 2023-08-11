#include "biort.hpp"

void ReadCini(const string& input_dir, const array<ChemTableEntry, MAXSPS>& chemtbl, ReactionNetwork& rttbl,
    Subcatchment& subcatch)
{
    char            file_name[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *file_pointer;
    int             line_number = 0;
    array<f64, MAXSPS> dummy_arr;

    sprintf(file_name, "input/%s/cini.txt", input_dir.c_str());
    file_pointer = fopen(file_name, "r");

    // printf("INITIAL CONCENTRATIONS:\n");

    // Read precipitation concentration
    // printf("\tPRECIPITATION:\n");
    FindLine(file_pointer, "PRECIPITATION", &line_number, cmdstr);
    ReadConc(file_pointer, rttbl.num_stc, chemtbl, &line_number, subcatch.prcp_conc, dummy_arr, dummy_arr, dummy_arr, dummy_arr, dummy_arr);


    // Read surface concentration  2021-05-07
    // printf("\tSURFACE:\n");
    FindLine(file_pointer, "SURFACE", &line_number, cmdstr);
    ReadConc(
        file_pointer, 
        rttbl.num_stc, 
        chemtbl, 
        &line_number, 
        subcatch.chms[SURFACE].tot_conc, 
        subcatch.chms[SURFACE].soil_parameters.ssa, 
        subcatch.chms[SURFACE].soil_parameters.q10, 
        subcatch.chms[SURFACE].soil_parameters.sw_thld, 
        subcatch.chms[SURFACE].soil_parameters.sw_exp, 
        subcatch.chms[SURFACE].soil_parameters.n_alpha);


    // Read upper zone concentration
    // printf("\tSHALLOW ZONE:\n");
    FindLine(file_pointer, "UZ", &line_number, cmdstr);
    ReadConc(
        file_pointer, 
        rttbl.num_stc, 
        chemtbl, 
        &line_number, 
        subcatch.chms[UZ].tot_conc, 
        subcatch.chms[UZ].soil_parameters.ssa, 
        subcatch.chms[UZ].soil_parameters.q10, 
        subcatch.chms[UZ].soil_parameters.sw_thld, 
        subcatch.chms[UZ].soil_parameters.sw_exp, 
        subcatch.chms[UZ].soil_parameters.n_alpha);

    // Read lower zone concentration
    // printf("\tDEEP ZONE:\n");
    FindLine(file_pointer, "LZ", &line_number, cmdstr);
    ReadConc(
        file_pointer, 
        rttbl.num_stc, 
        chemtbl, 
        &line_number, 
        subcatch.chms[LZ].tot_conc, 
        subcatch.chms[LZ].soil_parameters.ssa, 
        subcatch.chms[LZ].soil_parameters.q10, 
        subcatch.chms[LZ].soil_parameters.sw_thld, 
        subcatch.chms[LZ].soil_parameters.sw_exp, 
        subcatch.chms[LZ].soil_parameters.n_alpha);

    printf("\n\n");

    fclose(file_pointer);
}

void ReadConc(FILE *fp, int num_stc, const array<ChemTableEntry, MAXSPS>& chemtbl, int *lno, array<f64, MAXSPS>& tot_conc, array<f64, MAXSPS>& ssa, array<f64, MAXSPS>& q10, array<f64, MAXSPS>& sw_thld, array<f64, MAXSPS>& sw_exp, array<f64, MAXSPS>& n_alpha)
{
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             ind;
    int             convert;

    for (int kspc = 0; kspc < num_stc; kspc++)
    {
        NextLine(fp, cmdstr, lno);
        sscanf(cmdstr, "%s", temp_str);
        if (strcmp("pH", temp_str) == 0)
        {
            convert = 1;
        }

        ind = FindChem(temp_str, num_stc, chemtbl);
        if (ind < 0)
        {
            biort_printf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
            exit(EXIT_FAILURE);
        }

        if (chemtbl[ind].itype == MINERAL)
        {
            if (sscanf(cmdstr, "%*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf ", &tot_conc[ind], &ssa[ind], &q10[ind], &sw_thld[ind], &sw_exp[ind], &n_alpha[ind]) !=6)
            {
                biort_printf(VL_ERROR, "Error reading initial condition in %s at Line %d.\n", "cini.txt", *lno);
                exit(EXIT_FAILURE);
            }
            // printf("\t\tMINERAL: %s\n", chemtbl[ind].name);
            // printf("\t\t\t%-20s: %g\n", "Volume Fraction", tot_conc[ind]);
            // printf("\t\t\t%-20s: %g\n", "SSA", ssa[ind]);
            // printf("\t\t\t%-20s: %g\n", "Q10", q10[ind]);
            // printf("\t\t\t%-20s: %g\n", "SW Threshold", sw_thld[ind]);
            // printf("\t\t\t%-20s: %g\n", "SW Exponent", sw_exp[ind]);
            // printf("\t\t\t%-20s: %g\n", "N_alpha", n_alpha[ind]);
        }
        else
        {
            if (sscanf(cmdstr, "%*s %lf", &tot_conc[ind]) != 1)
            {
                biort_printf(VL_ERROR, "Error reading initial condition in %s at Line %d.\n", "cini.txt", *lno);
            }
            // printf("\t\tAQUEOUS: %s\n", chemtbl[ind].name);
            // printf("\t\t\t%-20s: %g\n", "Concentration", tot_conc[ind]);
        }

        tot_conc[ind] = (strcmp(chemtbl[ind].name, "pH") == 0 && convert == 1) ?
            pow(10, -tot_conc[ind]) : tot_conc[ind];
    }
}
