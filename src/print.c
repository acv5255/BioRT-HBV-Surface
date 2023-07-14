#include "biort.h"

void PrintHeader(FILE *file_pointer, int transpt, const ReactionNetwork *rttbl, const ChemTableEntry chemtbl[])
{
    char            chemn[MAXSTRING];
    char            tempstr[MAXSTRING];

    // Soil concentration file header
    fprintf(file_pointer, "%-15s",  "TIME");

    if (transpt == KIN_REACTION)
    {
        // Snow chemistry
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            // sprintf(tempstr, "%s_SNOW", chemn);
            sprintf(tempstr, "%s_SNOW", substr);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }

        // Surface zone
        for (int kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)  // 2021-06-29
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_SURFACE", substr);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }

        // Shallow zone
        for (int kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_UZ", substr);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }

        // Deep zone
        for (int kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_LZ", substr);
            // sprintf(tempstr, "%s_LZ", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }

        // Stream chemistry
        for (int kspc = 0; kspc < rttbl->num_stc + rttbl->num_ssc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_riv", substr);
            // sprintf(tempstr, "%s_riv", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }

        //===== Reaction rates =====//
        for (int kspc = 0; kspc < rttbl->num_min; kspc++)  // 2021-05-14
        {
            UnwrapParentheses(chemtbl[kspc + rttbl->num_stc - rttbl->num_min].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_rate_SF", substr);
            fprintf(file_pointer, "\t%-23s", tempstr);
        }

        for (int kspc = 0; kspc < rttbl->num_min; kspc++)  // 2021-05-14
        {
            UnwrapParentheses(chemtbl[kspc + rttbl->num_stc - rttbl->num_min].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_rate_UZ", substr);
            fprintf(file_pointer, "\t%-23s", tempstr);
        }

        for (int kspc = 0; kspc < rttbl->num_min; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc + rttbl->num_stc - rttbl->num_min].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_rate_LZ", substr);
            fprintf(file_pointer, "\t%-23s", tempstr);
        }
        //==========================//
    }
    else    // In transport mode, only print primary species
    {

        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_SNOW", substr);
            // sprintf(tempstr, "%s_SNOW", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_SURFACE", substr);
            // sprintf(tempstr, "%s_SURFACE", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_UZ", substr);
            // sprintf(tempstr, "%s_UZ", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_LZ", substr);
            // sprintf(tempstr, "%s_LZ", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            UnwrapParentheses(chemtbl[kspc].name, chemn);

            char substr[SUBSTRING_SIZE];
            strncpy(substr, chemn, SUBSTRING_SIZE);

            sprintf(tempstr, "%s_riv", substr);
            // sprintf(tempstr, "%s_riv", chemn);
            fprintf(file_pointer, "\t%-15s", tempstr);

        }
    }
    fprintf(file_pointer, "\n");

    // UNITS
    fprintf(file_pointer, "%-15s",  "YYYYMMDD");

    if (transpt == KIN_REACTION)
    {
        // Snow zones
        for (int kspc = 0; kspc < 1 * (rttbl->num_spc); kspc++)   // 2021-09-27, SNOW +SURFACE
        {
            fprintf(file_pointer, "\t%-15s", "mol/L");
        }
        // Surface, upper, lower, and stream chemistry
        for (int kspc = 0; kspc < 4 * (rttbl->num_stc + rttbl->num_ssc); kspc++)   // UZ + LZ + STREAM
        {
            fprintf(file_pointer, "\t%-15s", "mol/L");
        }
        // Kinetic reaction rates for surface, upper, and lower zones
        for (int kspc = 0; kspc < 3 * rttbl->num_min; kspc++)  // Surface + UZ + LZ
        {
            fprintf(file_pointer, "\t%-23s", "mol/m2/day");
        }
    }
    else    // In transport mode, only print primary species
    {
        for (int kspc = 0; kspc < 5 * rttbl->num_spc; kspc++)  // 2021-06-29
        {
            fprintf(file_pointer, "\t%-15s", "mol/L");
        }
    }
    fprintf(file_pointer, "\n");

    fflush(file_pointer);
}

void PrintDailyResults(FILE *fp, int transpt, int step, const ReactionNetwork *rttbl,
    const Subcatchment* subcatch)
{
    // Soil concentration file header
    fprintf(fp, "%-15d",  step);

    if (transpt == KIN_REACTION)
    {

        // Snow concentration
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)   // 2021-06-29
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[SNOW].prim_conc[kspc]);
        }
        // Surface concentration
        for (int kspc = 0; kspc < rttbl->num_stc; kspc++)   // 2021-06-29
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[SURFACE].prim_conc[kspc]);
        }
        // Surface secondary concentration
        for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)   // 2021-06-29
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[SURFACE].sec_conc[kspc]);
        }

        // Upper zone primary concentration
        for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[UZ].prim_conc[kspc]);
        }
        // Upper zone secondary concentration
        for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[UZ].sec_conc[kspc]);
        }

        // Lower zone primary concentration
        for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[LZ].prim_conc[kspc]);
        }

        // Lower zone secondary concentration
        for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[LZ].sec_conc[kspc]);
        }

        // Stream primary concentration
        for (int kspc = 0; kspc < rttbl->num_stc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[STREAM].prim_conc[kspc]);
        }

        // Stream secondary concentration
        for (int kspc = 0; kspc < rttbl->num_ssc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[STREAM].sec_conc[kspc]);
        }

        //===== Reaction rates =====//
        for (int kspc = 0; kspc < rttbl->num_min; kspc++)  // 2021-05-14
        {
           fprintf(fp, "\t%-23lg", subcatch->react_rate[SURFACE][kspc]);
        }
        for (int kspc = 0; kspc < rttbl->num_min; kspc++)
        {
            fprintf(fp, "\t%-23lg", subcatch->react_rate[UZ][kspc]);
        }
        for (int kspc = 0; kspc < rttbl->num_min; kspc++)
        {
            fprintf(fp, "\t%-23lg", subcatch->react_rate[LZ][kspc]);
        }
    }
    else    // In transport mode, only print primary species
    {
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[SNOW].prim_conc[kspc]);
        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)  // 2021-06-29
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[SURFACE].prim_conc[kspc]);
        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[UZ].prim_conc[kspc]);
        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
            fprintf(fp, "\t%-15lg", subcatch->chms[LZ].prim_conc[kspc]);
        }
        for (int kspc = 0; kspc < rttbl->num_spc; kspc++)
        {
        fprintf(fp, "\t%-15lg", subcatch->chms[STREAM].prim_conc[kspc]);
        }
    }
    fprintf(fp, "\n");

    fflush(fp);
}
