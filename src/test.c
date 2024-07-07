#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "cvectors.h"

#define EPSILON 1e-6
#define TEST_RUNS 1000

void test_add_vectors() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Vector v1 = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector v2 = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector* result = add_vectors(&v1, &v2);

        assert(result != NULL);
        assert(fabs(result->elements[0] - (v1.elements[0] + v2.elements[0])) < EPSILON);
        assert(fabs(result->elements[1] - (v1.elements[1] + v2.elements[1])) < EPSILON);

        free(result->elements);
        free(result);
    }
}

void test_subtract_vectors() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Vector v1 = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector v2 = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector* result = subtract_vectors(&v1, &v2);

        assert(result != NULL);
        assert(fabs(result->elements[0] - (v1.elements[0] - v2.elements[0])) < EPSILON);
        assert(fabs(result->elements[1] - (v1.elements[1] - v2.elements[1])) < EPSILON);

        free(result->elements);
        free(result);
    }
}

void test_dot_product() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Vector v1 = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector v2 = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        float result = dot_product(&v1, &v2);

        float expected = v1.elements[0] * v2.elements[0] + v1.elements[1] * v2.elements[1];
        assert(fabs(result - expected) < EPSILON);
    }
}

void test_cross_product() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Vector v1 = {3, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector v2 = {3, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector* result = cross_product(&v1, &v2);

        assert(result != NULL);
        assert(fabs(result->elements[0] - ((v1.elements[1] * v2.elements[2]) - (v1.elements[2] * v2.elements[1]))) < EPSILON);
        assert(fabs(result->elements[1] - (-((v1.elements[0] * v2.elements[2]) - (v1.elements[2] * v2.elements[0])))) < EPSILON);
        assert(fabs(result->elements[2] - ((v1.elements[0] * v2.elements[1]) - (v1.elements[1] * v2.elements[0]))) < EPSILON);

        free(result->elements);
        free(result);
    }
}

void test_transform_vector() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Matrix m = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Vector v = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        Vector* result = transform_vector(&m, &v);

        assert(result != NULL);
        assert(fabs(result->elements[0] - (m.elements[0][0] * v.elements[0] + m.elements[0][1] * v.elements[1])) < EPSILON);
        assert(fabs(result->elements[1] - (m.elements[1][0] * v.elements[0] + m.elements[1][1] * v.elements[1])) < EPSILON);

        free(result->elements);
        free(result);
    }
}

void test_scale_vector() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Vector v = {2, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}};
        float c = (float)rand() / RAND_MAX;
        Vector* result = scale_vector(c, &v);

        assert(result != NULL);
        assert(fabs(result->elements[0] - (v.elements[0] * c)) < EPSILON);
        assert(fabs(result->elements[1] - (v.elements[1] * c)) < EPSILON);

        free(result->elements);
        free(result);
    }
}

void test_normalize_vector() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Vector v = {2, (float[]){(float)rand() / RAND_MAX + 0.1, (float)rand() / RAND_MAX + 0.1}}; // Adding 0.1 to avoid zero vector
        Vector* result = normalize_vector(&v);

        float magnitude = sqrt(v.elements[0] * v.elements[0] + v.elements[1] * v.elements[1]);

        assert(result != NULL);
        assert(fabs(result->elements[0] - (v.elements[0] / magnitude)) < EPSILON);
        assert(fabs(result->elements[1] - (v.elements[1] / magnitude)) < EPSILON);

        free(result->elements);
        free(result);
    }
}

void test_add_matrices() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Matrix m1 = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Matrix m2 = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Matrix* result = add_matrices(&m1, &m2);

        assert(result != NULL);
        assert(fabs(result->elements[0][0] - (m1.elements[0][0] + m2.elements[0][0])) < EPSILON);
        assert(fabs(result->elements[0][1] - (m1.elements[0][1] + m2.elements[0][1])) < EPSILON);
        assert(fabs(result->elements[1][0] - (m1.elements[1][0] + m2.elements[1][0])) < EPSILON);
        assert(fabs(result->elements[1][1] - (m1.elements[1][1] + m2.elements[1][1])) < EPSILON);

        free_matrix(result);
    }
}

void test_subtract_matrices() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Matrix m1 = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Matrix m2 = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Matrix* result = subtract_matrices(&m1, &m2);

        assert(result != NULL);
        assert(fabs(result->elements[0][0] - (m1.elements[0][0] - m2.elements[0][0])) < EPSILON);
        assert(fabs(result->elements[0][1] - (m1.elements[0][1] - m2.elements[0][1])) < EPSILON);
        assert(fabs(result->elements[1][0] - (m1.elements[1][0] - m2.elements[1][0])) < EPSILON);
        assert(fabs(result->elements[1][1] - (m1.elements[1][1] - m2.elements[1][1])) < EPSILON);

        free_matrix(result);
    }
}

void test_multiply_matrices() {
    for (int i = 0; i < TEST_RUNS; i++) {
        Matrix m1 = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Matrix m2 = {2, (float*[]){(float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}, (float[]){(float)rand() / RAND_MAX, (float)rand() / RAND_MAX}}};
        Matrix* result = multiply_matrices(&m1, &m2);

        assert(result != NULL);
        assert(fabs(result->elements[0][0] - (m1.elements[0][0] * m2.elements[0][0] + m1.elements[0][1] * m2.elements[1][0])) < EPSILON);
        assert(fabs(result->elements[0][1] - (m1.elements[0][0] * m2.elements[0][1] + m1.elements[0][1] * m2.elements[1][1])) < EPSILON);
        assert(fabs(result->elements[1][0] - (m1.elements[1][0] * m2.elements[0][0] + m1.elements[1][1] * m2.elements[1][0])) < EPSILON);
        assert(fabs(result->elements[1][1] - (m1.elements[1][0] * m2.elements[0][1] + m1.elements[1][1] * m2.elements[1][1])) < EPSILON);

        free_matrix(result);
    }
}

int main() {
    srand(time(NULL));
    #ifdef VECTOR_TEST
    test_add_vectors();
    test_subtract_vectors();
    test_dot_product();
    test_cross_product();
    test_transform_vector();
    test_scale_vector();
    test_normalize_vector();
    #endif

    #ifdef MATRIX_TEST
    test_add_matrices();
    test_subtract_matrices();
    test_multiply_matrices();
    #endif

    printf("All tests passed successfully!\n");
    return 0;
}
