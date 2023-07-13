#include "biort.h"
#include "optparse.h"

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

void Log10Arr(const double src[MAXSPS], double dst[MAXSPS], int num_species) {
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

void FreeStruct(int *steps[], Subcatchment subcatch[])
{
    const int ksub = 0;
    free(*steps);

    free(subcatch[ksub].ws);
    free(subcatch[ksub].q);
    free(subcatch[ksub].tmp);
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
        {0, 0, 0}
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

void ComputeDependence(double tmpconc[MAXSPS], const double dep_mtx[MAXSPS][MAXSPS], const double keq[MAXSPS], int num_rows, int num_cols, int offset) {
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

void GetLogActivity(double tmpconc[MAXSPS], double gamma[MAXSPS], const double dep_mtx[MAXSPS][MAXSPS], const double keq[MAXSPS], int num_rows, int num_cols, int offset) {
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