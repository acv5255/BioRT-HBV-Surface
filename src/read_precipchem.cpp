#include "biort.hpp"

// Returns vector of date ints
vector<int> ReadPrecipChem(const char dir[], Subcatchment& subcatch, int num_stc, const array<ChemTableEntry, MAXSPS>& chemtbl,int mode)
{
    FILE           *file_pointer;
    char            file_name[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             pH_index = 0;
    int             pH_convert = 0;
    int             ind;
    

    if (mode == 0)
    {
        sprintf(file_name, "input/%s/precipchem.txt", dir);
        biort_printf(VL_NORMAL, "\nREADING TIME-SERIES PRECIPITATION CHEMISTRY from \"precipchem.txt\"\n");
    }

    else if (mode == 1)

    {
        sprintf(file_name, "input/%s/Numexp_precipchem.txt", dir);
        biort_printf(VL_NORMAL, "\nREADING TIME-SERIES PRECIPITATION CHEMISTRY from \"Numexp_precipchem.txt\"\n");
    }

    file_pointer = fopen(file_name, "r");

    const int ntime = CountLines(file_pointer, cmdstr, 0) - 1;
    rewind(file_pointer);

    // *steps = (int *)malloc(ntime * sizeof(int));
    vector<int> steps = vector<int>(ntime, BADVAL);

    subcatch.prcp_conc_time = vector(ntime, vector(num_stc, BADVAL));

    // read header to locate pH position
    for (int kspc = 0; kspc < num_stc + 1; kspc++)  // add one more column of date
    {
        fscanf(file_pointer, "%s", temp_str);
        printf("Looking for species: %s\n", temp_str);

        if (strcmp("pH", temp_str) == 0)
        {
            pH_convert = 1;
            pH_index = kspc - 1;
        }

        // also check chemical species, 0629
        if (kspc > 0)
        {
            ind = FindChem(temp_str, num_stc, chemtbl);
            if (ind < 0)
            {
                biort_printf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
                exit(EXIT_FAILURE);
            }
        }

    }

    for (int kstep = 0; kstep < ntime; kstep++)
    {

        fscanf(file_pointer, "%d", &steps[kstep]);    // Read model steps

        for (int kspc = 0; kspc < num_stc; kspc++)  // Read precipitation chemistry
        {
            if (kspc == pH_index && pH_convert == 1)
            {
                fscanf(file_pointer, "%lf", &subcatch.prcp_conc_time[kstep][kspc]);
                //printf("  step = %d, converting time-series precipitation pH (%lf) to ", kstep, subcatch.prcp_conc_time[kstep][kspc]);
                subcatch.prcp_conc_time[kstep][kspc] = pow(10, -subcatch.prcp_conc_time[kstep][kspc]);
                //printf("H+ concentration (%lf) \n", subcatch.prcp_conc_time[kstep][kspc]);
            }
            else
            {
                fscanf(file_pointer, "%lf", &subcatch.prcp_conc_time[kstep][kspc]);
            }
        }

    }

    fclose(file_pointer);

    // return ntime;
    return steps;
}
