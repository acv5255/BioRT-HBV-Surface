#include "biort.hpp"

void ReadPrecipChem(const char dir[], int *nsteps, int *steps[], Subcatchment subcatch[], int num_stc, const ChemTableEntry chemtbl[],int mode)
{
    FILE           *file_pointer;
    char            file_name[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             pH_index = 0;
    int             pH_convert = 0;
    int             ind;
    int             ntime;

    ntime = *nsteps;

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

    *nsteps = CountLines(file_pointer, cmdstr, 0) - 1;

    rewind(file_pointer);

    if ((ntime != *nsteps) & (mode == 1)){
        biort_printf(VL_ERROR,"\nNumber of time steps in \"Numexp_precipchem.txt\" should be same as in \"Numexp_Results.txt\" file.\n");
        exit(EXIT_FAILURE);
    }
    *steps = (int *)malloc(*nsteps * sizeof(int));


    subcatch->prcp_conc_time = (double **)malloc(*nsteps * sizeof(double *));

    for (int kstep = 0; kstep < *nsteps; kstep++)
    {
        subcatch->prcp_conc_time[kstep] = (double *)malloc(num_stc * sizeof(double));
    }


    // read header to locate pH position
    for (int kspc = 0; kspc < num_stc + 1; kspc++)  // add one more column of date
    {
        fscanf(file_pointer, "%s", temp_str);

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

    for (int kstep = 0; kstep < *nsteps; kstep++)
    {

        fscanf(file_pointer, "%d", &((*steps)[kstep]));    // Read model steps

        for (int kspc = 0; kspc < num_stc; kspc++)  // Read precipitation chemistry
        {
            if (kspc == pH_index && pH_convert == 1)
            {
                fscanf(file_pointer, "%lf", &subcatch->prcp_conc_time[kstep][kspc]);
                //printf("  step = %d, converting time-series precipitation pH (%lf) to ", kstep, subcatch->prcp_conc_time[kstep][kspc]);
                subcatch->prcp_conc_time[kstep][kspc] = pow(10, -subcatch->prcp_conc_time[kstep][kspc]);
                //printf("H+ concentration (%lf) \n", subcatch->prcp_conc_time[kstep][kspc]);
            }
            else
            {
                fscanf(file_pointer, "%lf", &subcatch->prcp_conc_time[kstep][kspc]);
            }
        }

    }

    fclose(file_pointer);
}