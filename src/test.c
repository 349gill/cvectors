#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "cvectors.h"

#define EPSILON 1e-6

void test_vector_operations() {
    Vector* v1 = create_3D_vector(1, 2, 3);
    Vector* v2 = create_3D_vector(4, 5, 6);
    Vector* sum = add_vectors(v1, v2);

    assert(fabs(sum->elements[0] - 5) < EPSILON);
    assert(fabs(sum->elements[1] - 7) < EPSILON);
    assert(fabs(sum->elements[2] - 9) < EPSILON);

    Vector* diff = subtract_vectors(v2, v1);
    assert(fabs(diff->elements[0] - 3) < EPSILON);
    assert(fabs(diff->elements[1] - 3) < EPSILON);
    assert(fabs(diff->elements[2] - 3) < EPSILON);
    float dot = dot_product(v1, v2);
    assert(fabs(dot - 32) < EPSILON);

    Vector* cross = cross_product(v1, v2);
    assert(fabs(cross->elements[0] - (-3)) < EPSILON);
    assert(fabs(cross->elements[1] - 6) < EPSILON);
    assert(fabs(cross->elements[2] - (-3)) < EPSILON);
    float mag = magnitude_vector(v1);
    assert(fabs(mag - sqrt(14)) < EPSILON);

    Vector* norm = normalize_vector(v1);
    assert(fabs(magnitude_vector(norm) - 1) < EPSILON);

    free_vector(v1);
    free_vector(v2);
    free_vector(sum);
    free_vector(diff);
    free_vector(cross);
    free_vector(norm);
}

void test_matrix_operations() {

    Matrix* m1 = create_matrix(4);
    Matrix* m2 = create_matrix(4);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            m1->elements[i][j] = i * 4 + j + 1;
            m2->elements[i][j] = (i * 4 + j + 1) * 2;
        }
    }
    Matrix* sum = add_matrices(m1, m2);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            assert(fabs(sum->elements[i][j] - (m1->elements[i][j] + m2->elements[i][j])) < EPSILON);
        }
    }

    Matrix* prod = multiply_matrices(m1, m2);
    assert(fabs(prod->elements[0][0] - 90) < EPSILON);
    assert(fabs(prod->elements[1][1] - 290) < EPSILON);
    assert(fabs(prod->elements[2][2] - 490) < EPSILON);
    assert(fabs(prod->elements[3][3] - 690) < EPSILON);

    Matrix* trans = transpose_matrix(m1);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            assert(fabs(trans->elements[i][j] - m1->elements[j][i]) < EPSILON);
        }
    }
    float det = determinant_matrix(m1);
    assert(fabs(det) < EPSILON);

    free_matrix(m1);
    free_matrix(m2);
    free_matrix(sum);
    free_matrix(prod);
    free_matrix(trans);
}

void test_matrix_from_quaternion() {
    float angle = M_PI / 2;
    Quaternion* q = create_quaternion(0, sin(angle/2), 0, cos(angle/2));

    Matrix* m = matrix_from_quaternion(q);
    float expected[4][4] = {
        { 0, 0, 1, 0},
        { 0, 1, 0, 0},
        {-1, 0, 0, 0},
        { 0, 0, 0, 1}
    };

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            assert(fabs(m->elements[i][j] - expected[i][j]) < EPSILON);
        }
    }

    free_quaternion(q);
    free_matrix(m);
}

void test_quaternion_operations() {

    Quaternion* q1 = create_quaternion(1, 2, 3, 4);
    Quaternion* q2 = create_quaternion(5, 6, 7, 8);
    Quaternion* sum = add_quaternions(q1, q2);
    assert(fabs(sum->x - 6) < EPSILON);
    assert(fabs(sum->y - 8) < EPSILON);
    assert(fabs(sum->z - 10) < EPSILON);
    assert(fabs(sum->w - 12) < EPSILON);

    float dot = dot_product_quaternion(q1, q2);
    assert(fabs(dot - 70) < EPSILON);

    Quaternion* norm = normalize_quaternion(q1);
    float mag = sqrt(norm->x*norm->x + norm->y*norm->y + norm->z*norm->z + norm->w*norm->w);
    assert(fabs(mag - 1) < EPSILON);

    test_matrix_from_quaternion();
    
    free_quaternion(q1);
    free_quaternion(q2);
    free_quaternion(sum);
    free_quaternion(norm);
}

int main() {
    #ifdef VECTOR
    test_vector_operations();
    #endif
    #ifdef MATRIX
    test_matrix_operations();
    #endif
    #ifdef QUAT
    test_quaternion_operations();
    #endif
    printf("\nAll tests passed.\n");
    return 0;
}