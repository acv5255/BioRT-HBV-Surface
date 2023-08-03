#include "biort.hpp"

double ReadParamToDouble(const char buffer[], const char keyword[], const char fn[], int line_number)
{
    char            optstr[MAXSTRING];
    int             match;
    double          value;

    match = sscanf(buffer, "%s %lf", optstr, &value);
    if (strcasecmp(keyword, optstr) != 0)
    {
        printf("Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
        exit(-1);
    }
    else if (match != 2)
    {
        printf("Failed to parse double value from line %s:%d with keyword '%s', exiting...\n", fn, line_number, keyword);
        exit(-1);
    }

    return value;
}

int ReadParamToInt(const char buffer[], const char keyword[], const char fn[], int line_number)
{
    char            optstr[MAXSTRING];
    int             match;
    int          value;

    match = sscanf(buffer, "%s %d", optstr, &value);
    if (strcasecmp(keyword, optstr) != 0)
    {
        printf("Expected keyword \"%s\", detected keyword \"%s\".\n", keyword, optstr);
        exit(-1);
    }
    else if (match != 2)
    {
        printf("Failed to parse int value from line %s:%d with keyword '%s', exiting...\n", fn, line_number, keyword);
        exit(-1);
    }

    return value;
}