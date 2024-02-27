#include <stdbool.h>
#include <stdlib.h>
#include "biort.h"
#include "unity.h"

// Compare two floating point numbers
bool cmp(double a, double b) {
    return fabs(a - b) < 1e-14;
}

void test_GetIAP() {
    // Test 1
    {
        double iap[MAXSPS];
        double actv[MAXSPS];
        double dep_kin[MAXSPS][MAXSPS];
    }

    TEST_FAIL();
}

void test_GetMonodTerm() {
    TEST_FAIL();
}

void test_GetInhibTerm() {
    TEST_FAIL();
}

void test_GetDependenceTerm() {
    TEST_FAIL();
}

void test_GetSecondarySpecies() {
    TEST_FAIL();
}

void test_SoilTempFactor() {
    TEST_FAIL();
}

void test_SoilMoistFactor() {
    TEST_FAIL();
}

void test_WTDepthFactor() {
    TEST_FAIL();
}