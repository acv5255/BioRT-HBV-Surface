#include "biort.hpp"
#include "optparse.hpp"

int roundi(double x)
{
    return (int)((x < 0.0) ? x - 0.5 : x + 0.5);
}

void SetZero(double arr[MAXSPS]) {
    for (size_t i = 0; i < MAXSPS; i++) {
        arr[i] = 0.0;
    }
}

void SetZeroRange(double arr[MAXSPS], int start, int end) {
    /* Set some subset of the data to equal zero */
    for (int i = start; i < end; i++) {
        arr[i] = 0.0;
    }
}

void Log10Arr(const array<f64, MAXSPS>& src, array<f64, MAXSPS>& dst, int num_species) {
    for (int i = 0; i < num_species; i++) {
        dst[i] = log10(src[i]);
    }
}

void Pow10Arr(const double src[MAXSPS], double dst[MAXSPS], int num_species) {
    for (int i = 0; i < num_species; i++) {
        dst[i] = pow(10.0, src[i]);
    }
}

double SumArr(const double arr[MAXSPS], int num_species) {
    double res;

    for (int i = 0; i < num_species; i++) {
        res += arr[i];
    }

    return res;
}

void WrapInParentheses(char *str)
{
    char            word[MAXSTRING];

    sprintf(word, "'%s'", str);
    strcpy(str, word);
}

void UnwrapParentheses(const char wrapped_str[], char str[])
{
    int             i, j = 0;

    for (i = 0; i < (int)strlen(wrapped_str); i++)
    {
        if (wrapped_str[i] != '\'')
        {
            str[j] = wrapped_str[i];
            j++;
        }
    }

    str[j] = '\0';
}

void FreeStruct(int *steps[])
{
    free(*steps);
}

void ParseCmdLineParam(int argc, char *argv[], char dir[])
{
    int             option;
    struct optparse options;
    struct optparse_long longopts[] = {
        {"brief",      'b', OPTPARSE_NONE},
        {"silent",     's', OPTPARSE_NONE},
        {"version",    'V', OPTPARSE_NONE},
        {"verbose",    'v', OPTPARSE_NONE},
        {0, 0, OPTPARSE_NONE}
    };

    optparse_init(&options, argv);

    while ((option = optparse_long(&options, longopts, NULL)) != -1)
    {
        switch (option)
        {
            case 'v':
                // Verbose mode
                verbose_mode = VL_VERBOSE;
                break;
            case 'b':
                // Brief mode
                verbose_mode = VL_BRIEF;
                break;
            case 's':
                // Silent mode
                verbose_mode = VL_SILENT;
                break;
            case 'V':
                // Print version number
                printf("HBV-BioRT Version %s\n", VERSION);
                exit(EXIT_SUCCESS);
                break;
            case '?':
                biort_printf(VL_ERROR, "Option not recognizable %s\n", options.errmsg);
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }

        fflush(stdout);
    }

    if (options.optind >= argc)
    {
        biort_printf(VL_ERROR, "Error:You must specify the name of input directory!\n"
            "Usage: ./biort [-b] [-v] [-V]"
            " <project directory>\n"
            "    -b Brief mode\n"
            "    -V Version number\n"
            "    -v Verbose mode\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        // Parse remaining arguments
        strcpy(dir, optparse_arg(&options));
    }
}

void ComputeDependence(array<f64, MAXSPS>& tmpconc, const array<array<f64, MAXSPS>, MAXSPS>& dep_mtx, const array<f64, MAXSPS>& keq, int num_rows, int num_cols, int offset) {
    for (int i = 0; i < num_rows; i++)
    {
        tmpconc[i + offset] = 0.0;
        for (int j = 0; j < num_cols; j++)
        {
            tmpconc[i + offset] += tmpconc[j] * dep_mtx[i][j];
        }
        tmpconc[i + offset] -= keq[i];
    }
}

void GetLogActivity(array<f64, MAXSPS>& tmpconc, double gamma[MAXSPS], const array<array<f64, MAXSPS>, MAXSPS>& dep_mtx, const array<f64, MAXSPS>& keq, int num_rows, int num_cols, int offset) {
    for (int i = 0; i < num_rows; i++)
    {
        tmpconc[i + offset] = 0.0;
        for (int j = 0; j < num_cols; j++)
        {
            tmpconc[i + offset] += (tmpconc[j] + gamma[j]) * dep_mtx[i][j];
        }
        tmpconc[i + offset] -= keq[i] + gamma[i + offset];
    }
}

bool CheckArrayForNan(const array<f64, MAXSPS>& arr) {
    for (int i = 0; i < MAXSPS; i++) {
        if (std::isinf(arr[i])) {
            return false;
        }
    }

    return true;
}

void ErrorOnArrayNan(const array<f64, MAXSPS>& arr, const char* array_name, const char* filename, const int line_number) {
    /* If an array contains a nan value, exit, or else do nothing */
    if (!CheckArrayForNan(arr)) {
        printf("Array '%s' in file '%s' near line %d contains nonfinite value, exiting...\n", array_name, filename, line_number);
        PrintArray(arr);
        exit(-1);
    }

    return;
}

void PrintArray(const array<f64, MAXSPS>& arr) {
    // Print a constant sized array
    for (int i = 0; i < MAXSPS; i++) {
        printf("%.4g  ", arr[i]);
    }
    printf("\n");
}

void CheckChmsForNonFinite(const ChemicalState* chms, const char* filename, const int line_number) {
    /* Check a chemical state for any nan values */
    const char* msg = "Chemical state entry ChemicalState->%s contained nonfinite value in file '%s' near line %d: \n";

    if (!CheckArrayForNan(chms->tot_conc)) {
        printf(msg, "tot_conc", filename, line_number);
        PrintArray(chms->tot_conc);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->prim_conc)) {
        printf(msg, "prim_conc", filename, line_number);
        PrintArray(chms->prim_conc);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->sec_conc)) {
        printf(msg, "sec_conc", filename, line_number);
        PrintArray(chms->sec_conc);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->prim_actv)) {
        printf(msg, "prim_actv", filename, line_number);
        PrintArray(chms->prim_actv);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->soil_parameters.ssa)) {
        printf(msg, "ssa", filename, line_number);
        PrintArray(chms->soil_parameters.ssa);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->soil_parameters.sw_thld)) {
        printf(msg, "sw_thld", filename, line_number);
        PrintArray(chms->soil_parameters.sw_thld);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->soil_parameters.sw_exp)) {
        printf(msg, "sw_exp", filename, line_number);
        PrintArray(chms->soil_parameters.sw_exp);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->soil_parameters.q10)) {
        printf(msg, "q10", filename, line_number);
        PrintArray(chms->soil_parameters.q10);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->soil_parameters.n_alpha)) {
        printf(msg, "n_alpha", filename, line_number);
        PrintArray(chms->soil_parameters.n_alpha);
        exit(-1);
    }
    if (!CheckArrayForNan(chms->tot_mol)) {
        printf(msg, "tot_mol", filename, line_number);
        PrintArray(chms->tot_mol);
        exit(-1);
    }
}

void ErrorOnArrayNanIter(const array<f64, MAXSPS>& arr, const char* array_name, const char* filename, const int line_number, const int iter) {
    /* If an array contains a nan value, exit, or else do nothing */
    if (!CheckArrayForNan(arr)) {
        printf("Array '%s' in file '%s' near line %d on iter %d contains NaN value, exiting...\n", array_name, filename, line_number, iter);
        PrintArray(arr);
        exit(-1);
    }

    return;   
}

bool CheckNonzeroRanged(const array<f64, MAXSPS>& arr, const int num) {
    const double TOL = 1e-15;
    for (int i = 0; i < num; i++) {
        if (fabs(arr[i]) < TOL) return false;
    }
    return true;
}

void ErrOnZeroRanged(const char* filename, const char* arr_name, const int line_number, const array<f64, MAXSPS>& arr, const int num) {
    if (!CheckNonzeroRanged(arr, num)) {
        printf("Array '%s' in file '%s' near line '%d' contains zero value (ranged) when it shouldn't, exiting...\n", arr_name, filename, line_number);
        PrintArray(arr);
        exit(-1);
    }
    return;
}

bool CompareFloats(const double lhs, const double rhs) {
    const double TOL = 1e-12;

    return fabs(lhs - rhs) < TOL;
}

bool CompareArray(const array<f64, MAXSPS>& lhs, const array<f64, MAXSPS>& rhs) {
    for (unsigned int i = 0; i < lhs.size(); i++) {
        if (CompareFloats(lhs[i], rhs[i])) return false;
    }

    return true;
}

bool CompareChemicalState(const ChemicalState* lhs, const ChemicalState* rhs) {
    return CompareArray(lhs->tot_conc, rhs->tot_conc)
        && CompareArray(lhs->prim_conc, rhs->prim_conc)
        && CompareArray(lhs->prim_actv, rhs->prim_actv)
        && CompareArray(lhs->sec_conc, rhs->sec_conc)
        && CompareArray(lhs->soil_parameters.ssa, rhs->soil_parameters.ssa)
        && CompareArray(lhs->soil_parameters.sw_thld, rhs->soil_parameters.sw_thld)
        && CompareArray(lhs->soil_parameters.sw_exp, rhs->soil_parameters.sw_exp)
        && CompareArray(lhs->soil_parameters.q10, rhs->soil_parameters.q10)
        && CompareArray(lhs->soil_parameters.n_alpha, rhs->soil_parameters.n_alpha)
        && CompareArray(lhs->tot_mol, rhs->tot_mol);
}

// Compare two subcatchments and determine whether the data in each is the same

void PrintMatrix(const realtype** mat, const int nrows, const int ncols) {
    // Print a matrix
    for (int i =0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            printf("%8g\t", mat[i][j]);
        }
        printf("\n");
    }
}

// Print all the values in a chemical state
void PrintChemicalState(const ChemicalState* chms) {
    printf("Total concentration: ");
    PrintArray(chms->tot_conc);

    printf("Total moles: ");
    PrintArray(chms->tot_mol);
    
    printf("Primary concentration: ");
    PrintArray(chms->prim_conc);

    printf("Secondary concentration: ");
    PrintArray(chms->sec_conc);

    // printf("SSA: ");
    // PrintArray(chms->ssa);

    // printf("SW Threshold: ");
    // PrintArray(chms->sw_thld);

    // printf("SW Exponent: ");
    // PrintArray(chms->sw_exp);

    // printf("Q10: ");
    // PrintArray(chms->q10);

    // printf("N_alpha: ");
    // PrintArray(chms->n_alpha);

    return;
}
