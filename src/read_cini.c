#include "biort.h"

void ReadCini(const char dir[], const ChemTableEntry *chemtbl, ReactionNetwork *rttbl,
    Subcatchment subcatch[])
{
    char            file_name[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *file_pointer;
    int             line_number = 0;
    double          dummy[MAXSPS];

    sprintf(file_name, "input/%s/cini.txt", dir);
    file_pointer = fopen(file_name, "r");

    // Read precipitation concentration
    FindLine(file_pointer, "PRECIPITATION", &line_number, cmdstr);
    ReadConc(file_pointer, rttbl->num_stc, chemtbl, &line_number, subcatch->prcp_conc, dummy, dummy, dummy, dummy, dummy);

    // Read surface concentration  2021-05-07
    //FindLine(fp, "SURFACE", &lno, cmdstr);
    //ReadConc(fp, rttbl->num_stc, chemtbl, &lno, subcatch->chms[SURFACE].tot_conc, subcatch->chms[SURFACE].ssa);

    // Read upper zone concentration
    FindLine(file_pointer, "UZ", &line_number, cmdstr);
    ReadConc(file_pointer, rttbl->num_stc, chemtbl, &line_number, subcatch->chms[UZ].tot_conc, subcatch->chms[UZ].ssa, subcatch->chms[UZ].q10, subcatch->chms[UZ].sw_thld, subcatch->chms[UZ].sw_exp, subcatch->chms[UZ].n_alpha);

    // Read lower zone concentration
    FindLine(file_pointer, "LZ", &line_number, cmdstr);
    ReadConc(file_pointer, rttbl->num_stc, chemtbl, &line_number, subcatch->chms[LZ].tot_conc, subcatch->chms[LZ].ssa, subcatch->chms[LZ].q10, subcatch->chms[LZ].sw_thld, subcatch->chms[LZ].sw_exp, subcatch->chms[LZ].n_alpha);

    fclose(file_pointer);
}

void ReadConc(FILE *fp, int num_stc, const ChemTableEntry chemtbl[], int *lno, double tot_conc[], double ssa[], double q10[], double sw_thld[], double sw_exp[], double n_alpha[])
{
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             ind;
    int             convert;
    int             kspc;

    for (kspc = 0; kspc < num_stc; kspc++)
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
        }
        else
        {
            if (sscanf(cmdstr, "%*s %lf", &tot_conc[ind]) != 1)
            {
                biort_printf(VL_ERROR, "Error reading initial condition in %s at Line %d.\n", "cini.txt", *lno);
            }
        }

        tot_conc[ind] = (strcmp(chemtbl[ind].name, "pH") == 0 && convert == 1) ?
            pow(10, -tot_conc[ind]) : tot_conc[ind];
    }
}
