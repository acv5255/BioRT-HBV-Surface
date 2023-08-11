#include "biort.hpp"

void ReadChem(const string& input_dir, ControlData& ctrl, ReactionNetwork& rttbl, array<ChemTableEntry, MAXSPS>& chemtbl,
    array<KineticTableEntry, MAXSPS>& kintbl)
{
    int             i;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    char            chemn[MAXSPS][MAXSTRING];
    char            file_name[MAXSTRING];
    int             p_type[MAXSPS];
    int             lno = 0;
    FILE           *file_pointer;

    // READ CHEM.TXT FILE
    sprintf(file_name, CHEM_FILE_DIR, input_dir.c_str());
    file_pointer = fopen(file_name, "r");

    biort_printf(VL_NORMAL, "\nBIORT CONTROL PARAMETERS\n");

    NextLine(file_pointer, cmdstr, &lno);
    ctrl.recycle = ReadParamToInt(cmdstr, CHEM_RECYCLE_ID, file_name, lno);
    biort_printf(VL_NORMAL, "  Forcing recycle %d time(s). \n", ctrl.recycle);

    NextLine(file_pointer, cmdstr, &lno);
    ctrl.use_activity = ReadParamToInt(cmdstr, CHEM_ACTIVITY_ID, file_name, lno);
    biort_printf(VL_NORMAL, "  Activity correction is set to %d. \n", ctrl.use_activity);

    NextLine(file_pointer, cmdstr, &lno);
    ctrl.transport_only = ReadParamToInt(cmdstr, CHEM_TRANSPORT_ONLY_ID, file_name, lno);

    if (ctrl.transport_only == KIN_REACTION) {
        biort_printf(VL_NORMAL, "  Transport only mode disabled.\n");
    }
    else if (ctrl.transport_only == TRANSPORT_ONLY) {
        biort_printf(VL_NORMAL, "  Transport only mode enabled. \n");
    }
    else {
        printf("Unkown simulation mode flag, must be 1 (transport only) or 0 (kinetic reactions)");
        exit(-1);
    }

    NextLine(file_pointer, cmdstr, &lno);  // 2021-05-20
    ctrl.variable_precipchem = ReadParamToInt(cmdstr, CHEM_PRECIPCHEM_ID, file_name, lno);
    switch (ctrl.variable_precipchem)
    {
        case 0:
            biort_printf(VL_NORMAL, "  Using constant precipitation chemistry in cini.txt. \n");
            break;
        case 1:
            biort_printf(VL_NORMAL, "  Using time-series precipitation chemistry in precipchem.txt. \n");
            break;
        default:
            break;
    }

    NextLine(file_pointer, cmdstr, &lno);  // 2021-09-09
    ctrl.precipchem_numexp = ReadParamToInt(cmdstr, CHEM_NUMEXP_ID, file_name, lno);
    switch (ctrl.precipchem_numexp)
    {
        case 0:
            biort_printf(VL_NORMAL, "  Using same precipitation chemistry during warmup and simulation run. \n");
            break;
        case 1:
            biort_printf(VL_NORMAL, "  Using different precipitation chemistry during warmup and simulation run. \n");
            break;
        default:
            break;
    }

    NextLine(file_pointer, cmdstr, &lno);
    rttbl.tmp = ReadParamToDouble(cmdstr, CHEM_TEMPERATURE_ID, file_name, lno);
    biort_printf(VL_NORMAL, "  Temperature = %3.1f \n", rttbl.tmp);

    // Count numbers of species and reactions
    FindLine(file_pointer, CHEM_PRIMARY_SPECIES_ID, &lno, file_name);
    rttbl.num_stc = CountLines(file_pointer, cmdstr, 1, "SECONDARY_SPECIES");
    rttbl.num_ssc = CountLines(file_pointer, cmdstr, 1, "MINERAL_KINETICS");
    rttbl.num_mkr = CountLines(file_pointer, cmdstr, 1, "PRECIPITATION_CONC");
    rttbl.num_akr = 0;     // Not implemented yet

    // Primary species block
    biort_printf(VL_NORMAL, "\nPRIMARY SPECIES\n");
    biort_printf(VL_NORMAL, "  %d chemical species specified. \n", rttbl.num_stc);
    FindLine(file_pointer, "BOF", &lno, file_name);
    FindLine(file_pointer, CHEM_PRIMARY_SPECIES_ID, &lno, file_name);

    rttbl.num_spc = 0;
    rttbl.num_ads = 0;
    rttbl.num_cex = 0;
    rttbl.num_min = 0;

    for (i = 0; i < rttbl.num_stc; i++)
    {
        NextLine(file_pointer, cmdstr, &lno);
        if (sscanf(cmdstr, "%s", chemn[i]) != 1)
        {
            biort_printf(VL_ERROR, "Error reading primary_species in %s near Line %d.\n", file_name, lno);
        }
        p_type[i] = SpeciesType(input_dir.c_str(), chemn[i]);

        switch (p_type[i])
        {
            case 0:     // Species type is 0 when it is not found in the database.
                biort_printf(VL_ERROR, "Error finding primary species %s in the database.\n", chemn[i]);
                exit(EXIT_FAILURE);
            case AQUEOUS:
                rttbl.num_spc++;
                break;
            case ADSORPTION:
                rttbl.num_ads++;
                break;
            case CATION_ECHG:
                rttbl.num_cex++;
                break;
            case MINERAL:
                rttbl.num_min++;
                break;
            case SECONDARY:
                biort_printf(VL_ERROR, "%s is a secondary species, but is listed as a primary species.\n"
                    "Error at Line %d in %s.\n", chemn[i], lno, file_name);
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
    }

    biort_printf(VL_NORMAL, "  %d aqueous species specified. \n", rttbl.num_spc);
    biort_printf(VL_NORMAL, "  %d surface complexation specified. \n", rttbl.num_ads);
    biort_printf(VL_NORMAL, "  %d cation exchange specified. \n", rttbl.num_cex);
    biort_printf(VL_NORMAL, "  %d minerals specified. \n", rttbl.num_min);

    SortChem(chemn, p_type, rttbl.num_stc, chemtbl);

    // Number of species that others depend on
    rttbl.num_sdc = rttbl.num_stc - rttbl.num_min;

    // Secondary_species block
    biort_printf(VL_NORMAL, "\nSECONDARY SPECIES\n");
    biort_printf(VL_NORMAL, "  %d secondary species specified. \n", rttbl.num_ssc);
    FindLine(file_pointer, CHEM_SECONDARY_SPECIES_ID, &lno, file_name);
    for (i = 0; i < rttbl.num_ssc; i++)
    {
        NextLine(file_pointer, cmdstr, &lno);
        chemtbl[rttbl.num_stc + i].name.resize(MAXSTRING);
        if (sscanf(cmdstr, "%s", chemtbl[rttbl.num_stc + i].name.data()) != 1)
        {
            biort_printf(VL_ERROR, "Error reading secondary_species in %s near Line %d.\n", file_name, lno);
        }

        if (SpeciesType(input_dir.c_str(), chemtbl[rttbl.num_stc + i].name.c_str()) == 0)
        {
            biort_printf(VL_ERROR, "Error finding secondary species %s in the database.\n",
                chemtbl[rttbl.num_stc + i].name);
            exit(EXIT_FAILURE);
        }
    }

    // Minerals block
    biort_printf(VL_NORMAL, "\nMINERAL KINETIC REACTIONS\n");
    biort_printf(VL_NORMAL, "  %d mineral kinetic reaction(s) specified. \n", rttbl.num_mkr);
    FindLine(file_pointer, CHEM_MINERAL_KINETICS_ID, &lno, file_name);

    for (i = 0; i < rttbl.num_mkr; i++)
    {
        NextLine(file_pointer, cmdstr, &lno);
        if (sscanf(cmdstr, "%s %*s %s", temp_str, kintbl[i].label) != 2)
        {
            biort_printf(VL_ERROR, "Error reading mineral information in %s near Line %d.\n", file_name, lno);
            exit(EXIT_FAILURE);
        }

        biort_printf(VL_NORMAL, "  Kinetic reaction on '%s' is specified, label '%s'.\n", temp_str, kintbl[i].label);

        kintbl[i].position = FindChem(temp_str, rttbl.num_stc, chemtbl);

        if (kintbl[i].position < 0)
        {
            biort_printf(VL_ERROR, "Error finding mineral %s in species table.\n", temp_str);
            exit(EXIT_FAILURE);
        }
        else
        {
            biort_printf(VL_NORMAL, "  Position_check (num_mkr[i] vs num_stc[j]) (%d, %d)\n", i, kintbl[i].position);
        }
    }

    fclose(file_pointer);
}

// This subroutine is used to find out what the input species is.
//   0) not found within database
//   1) aqueous
//   2) adsorption
//   3) cation exchange
//   4) mineral
int SpeciesType(const char dir[], const char chemn[])
{
    char            fn[MAXSTRING];
    char            tempn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    FILE           *fp;

    sprintf(fn, "input/%s/cdbs.txt", dir);
    fp = fopen(fn, "r");

    if (strcmp(chemn, "pH") == 0)
    {
        fclose(fp);
        return AQUEOUS;
    }

    sprintf(tempn, "'%s'", chemn);

    FindLine(fp, "BOF", &lno, "cdbs.txt");

    NextLine(fp, cmdstr, &lno);
    while (MatchWrappedKey(cmdstr, "'End of primary'") != 0)
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            fclose(fp);
            return AQUEOUS;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (MatchWrappedKey(cmdstr, "'End of secondary'") != 0)
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            fclose(fp);
            return 5;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (MatchWrappedKey(cmdstr, "'End of minerals'") != 0)
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            fclose(fp);
            return MINERAL;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (strcmp(cmdstr, "End of surface complexation\r\n") != 0 &&
        strcmp(cmdstr, "End of surface complexation\n") != 0)
    {
        // Notice that in CrunchFlow database, starting from surface complexation, there is not apostrophe marks around
        // blocking keywords
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            fclose(fp);
            return ADSORPTION;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (!feof(fp))
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            fclose(fp);
            return CATION_ECHG;
        }
        NextLine(fp, cmdstr, &lno);
    }

    fclose(fp);

    return 0;
}

int MatchWrappedKey(const char cmdstr[], const char key[])
{
    char            optstr[MAXSTRING];

    if (sscanf(cmdstr, "'%[^']'", optstr) != 1)
    {
        return 1;
    }
    else
    {
        WrapInParentheses(optstr);
        return (strcmp(optstr, key) == 0) ? 0 : 1;
    }
}

void SortChem(char chemn[MAXSPS][MAXSTRING], const int p_type[MAXSPS], int nsps, array<ChemTableEntry, MAXSPS>& chemtbl)
{
    int             i, j;
    int             temp;
    int             rank[MAXSPS];
    int             ranked_type[MAXSPS];

    for (i = 0; i < nsps; i++)
    {
        rank[i] = i;
        ranked_type[i] = p_type[i];
    }

    for (i = 0; i < nsps - 1; i++)
    {
        for (j = 0; j < nsps - i - 1; j++)
        {
            if (ranked_type[j] > ranked_type[j + 1])
            {
                temp = rank[j];
                rank[j] = rank[j + 1];
                rank[j + 1] = temp;

                temp = ranked_type[j];
                ranked_type[j] = ranked_type[j + 1];
                ranked_type[j + 1] = temp;
            }
        }
    }

    for (i = 0; i < nsps; i++)
    {
        strcpy(chemtbl[i].name.data(), chemn[rank[i]]);
        chemtbl[i].itype = p_type[rank[i]];
    }
}

int FindChem(const char chemn[MAXSTRING], int nsps, const array<ChemTableEntry, MAXSPS>& chemtbl)
{
    int             i;
    int             ind = BADVAL;

    for (i = 0; i < nsps; i++)
    {
        if (strcmp(chemn, chemtbl[i].name.c_str()) == 0)
        {
            ind = i;
            break;
        }
    }

    return ind;
}
