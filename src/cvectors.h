/*
 *  cvectors.h
 *  
 *  First Commit: May 27, 2024
 *      Author: Harsh Gill
 *
 */

#ifndef CVECTORS_H
#define CVECTORS_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int dimension;
    float *elements;
} Vector;

typedef struct {
    int dimension;
    float **elements;
} Matrix;

typedef struct {
    float x, y, z, w;
} Quaternion;

Vector* create_vector(int dimension);
Vector* create_3D_vector(float x, float y, float z);
Vector* create_4D_vector(float x, float y, float z, float w);
void free_vector(Vector* v);
Vector* duplicate_vector(Vector* v);
Vector* add_vectors(Vector* v1, Vector* v2);
Vector* subtract_vectors(Vector* v1, Vector* v2);
float angle_between_vectors(Vector* v1, Vector* v2);
float distance_between_vectors(Vector* v1, Vector* v2);
float dot_product(Vector* v1, Vector* v2);
Vector* cross_product(Vector* v1, Vector* v2);
Vector* transform_vector(Matrix* m, Vector* v);
Vector* scale_vector(float c, Vector* v);
float magnitude_vector(Vector* v);
Vector* normalize_vector(Vector* v);
Vector* reflection_vector(Vector* i, Vector* n);
Vector* max_vector(Vector* v1, Vector* v2);
Vector* min_vector(Vector* v1, Vector* v2);
Vector* projection_vector(Vector* v1, Vector* v2);

Matrix* create_matrix(int dimension);
void free_matrix(Matrix* m);
Matrix* add_matrices(Matrix* m1, Matrix* m2);
Matrix* subtract_matrices(Matrix* m1, Matrix* m2);
Matrix* multiply_matrices(Matrix* m1, Matrix* m2);
Matrix* transpose_matrix(Matrix* m);
Matrix* scale_matrix(float c, Matrix* m);
Matrix* scale_matrix_anisotropic(float x, float y, float z, Matrix* m);
float determinant_matrix(Matrix* m);
Matrix* translate_matrix(float x, float y, float z);
Matrix* translate_matrix_in_place(Matrix* m, float x, float y, float z);
Matrix* outer_product(Vector* v1, Vector* v2);
Matrix* duplicate_matrix(Matrix* m);
Matrix* identity_matrix(int dimension);
Matrix* rotate_matrix(Matrix* m, float x, float y, float z, float angle);
Matrix* rotate_matrix_X(Matrix* m, float angle);
Matrix* rotate_matrix_Y(Matrix* m, float angle);
Matrix* rotate_matrix_Z(Matrix* m, float angle);
Matrix* invert_matrix(Matrix* m);
Matrix* look_at_matrix(Vector* eye, Vector* center, Vector* up);
Matrix* orthographic_projection_matrix(float left, float right, float bottom, float top, float near, float far);
Matrix* perspective_projection_matrix(float fov, float aspect, float near, float far);

Quaternion* create_quaternion(float x, float y, float z, float w);
void free_quaternion(Quaternion* q);
Quaternion* add_quaternions(Quaternion* q1, Quaternion* q2);
Quaternion* subtract_quaternions(Quaternion* q1, Quaternion* q2);
float dot_product_quaternion(Quaternion* q1, Quaternion* q2);
Quaternion* scale_quaternion(float c, Quaternion* q);
Quaternion* normalize_quaternion(Quaternion* q);
Quaternion* identity_quaternion();
Quaternion* quaternion_conjugate(Quaternion* q);

Vector* projection_vector(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return NULL;
    float n = dot_product(v1, v2) / dot_product(v1, v1);
    return scale_vector(n, v1);
}

Vector* max_vector(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return NULL;
    Vector* max = create_vector(v1->dimension);
    for (int i = 0; i < v1->dimension; i++)
        max->elements[i] = v1->elements[i] > v2->elements[i] ? v1->elements[i] : v2->elements[i];
    return max;
}

Vector* create_3D_vector(float x, float y, float z) {
    Vector* v = create_vector(3);
    v->elements[0] = x;
    v->elements[1] = y;
    v->elements[2] = z;
    return v;
}

Vector* create_4D_vector(float x, float y, float z, float w) {
    Vector* v = create_vector(4);
    v->elements[0] = x;
    v->elements[1] = y;
    v->elements[2] = z;
    v->elements[3] = w;
    return v;
}

float angle_between_vectors(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return 0.f;

    float magnitude = magnitude_vector(v1) * magnitude_vector(v2);
    if (magnitude == 0)
        return 0.f;

    return acos(dot_product(v1, v2) / magnitude);
}

float distance_between_vectors(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return -1;

    float sum = 0, diff;
    for (int i = 0; i < v1->dimension; i++) {
        diff = v1->elements[i] - v2->elements[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

Vector* min_vector(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return NULL;
    Vector* min = create_vector(v1->dimension);
    for (int i = 0; i < v1->dimension; i++)
        min->elements[i] = v1->elements[i] < v2->elements[i] ? v1->elements[i] : v2->elements[i];
    return min;
}

Vector* create_vector(int dimension) {
    Vector* v = (Vector*) malloc(sizeof(Vector));
    v->dimension = dimension;
    v->elements = (float *) malloc(dimension * sizeof(float));
    return v;
}

void free_vector(Vector* v) {
    free(v->elements);
    free(v);
}

Vector* duplicate_vector(Vector* v) {
    Vector* d = create_vector(v->dimension);
    for (int i = 0; i < d->dimension; i++)
        d->elements[i] = v->elements[i];
    return d;
}

Vector* add_vectors(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return NULL;

    Vector* v = (Vector*) malloc(sizeof(Vector));
    v->dimension = v1->dimension;
    v->elements = (float *) malloc(v->dimension * sizeof(float));

    for (int i = 0; i < v->dimension; i++)
        v->elements[i] = v1->elements[i] + v2->elements[i];

    return v;
}

Vector* subtract_vectors(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension)
        return NULL;
    
    Vector* v = (Vector*) malloc(sizeof(Vector));
    v->dimension = v1->dimension;
    v->elements = (float *) malloc(v->dimension * sizeof(float));

    for (int i = 0; i < v->dimension; i++)
        v->elements[i] = v1->elements[i] - v2->elements[i];

    return v;
}

float dot_product(Vector* v1, Vector* v2) {
    if (v1->dimension != v2->dimension) 
        return 0;

    float x = 0;
    for (int i = 0; i < v1->dimension; i++)
        x += v1->elements[i] * v2->elements[i];

    return x;
}

Vector* cross_product(Vector* v1, Vector* v2) {
    if (v1->dimension == 3 && v2->dimension == 3) {
        Vector* v = (Vector*) malloc(sizeof(Vector));
        v->dimension = v1->dimension;
        v->elements = (float *) malloc(v->dimension * sizeof(float));
        v->elements[0] = (v1->elements[1] * v2->elements[2]) - (v1->elements[2] * v2->elements[1]);
        v->elements[1] =  - ((v1->elements[0] * v2->elements[2]) - (v1->elements[2] * v2->elements[0]));
        v->elements[2] =  (v1->elements[0] * v2->elements[1]) - (v1->elements[1] * v2->elements[0]);
        return v;
    }
    if ((v1->dimension == 4 && v2->dimension == 4)) {
        Vector* v = (Vector*) malloc(sizeof(Vector));
        v->dimension = v1->dimension;
        v->elements = (float *) malloc(v->dimension * sizeof(float));
        v->elements[0] = (v1->elements[1] * v2->elements[2]) - (v1->elements[2] * v2->elements[1]);
        v->elements[1] =  - ((v1->elements[0] * v2->elements[2]) - (v1->elements[2] * v2->elements[0]));
        v->elements[2] =  (v1->elements[0] * v2->elements[1]) - (v1->elements[1] * v2->elements[0]);
        v->elements[3] = 1.f;
        return v;
    }

    return NULL;
}


Vector* transform_vector(Matrix* m, Vector* v) {
    if (m->dimension != v->dimension)
        return NULL;
    Vector* t = (Vector*) malloc(sizeof(Vector));
    t->dimension = v->dimension;
    t->elements = (float *) malloc(v->dimension * sizeof(float));
    for (int i = 0; i < t->dimension; i++) {
        t->elements[i] = 0;
        for (int j = 0; j < t->dimension; j++) {
            t->elements[i] += (m->elements[i][j]) * (v->elements[j]);
        }
    }
    return t;
}

Vector* scale_vector(float c, Vector* v) {
    Vector* s = (Vector*) malloc(sizeof(Vector));
    s->dimension = v->dimension;
    s->elements = (float *) malloc(v->dimension * sizeof(float));
    for (int i = 0; i < v->dimension; i++){
        s->elements[i] = v->elements[i] * c;
    }
    return s;
}

float magnitude_vector(Vector* v) {
    float x = 0;
    for (int i = 0; i < v->dimension; i++){
        x += pow(v->elements[i], 2);
    }

    return sqrt(x);
}

Vector* normalize_vector(Vector* v) {
    float magnitude = magnitude_vector(v);
    
    if (magnitude == 0)
        return 0;

    Vector* n = (Vector*) malloc(sizeof(Vector));
    n->dimension = v->dimension;
    n->elements = (float *) malloc(v->dimension * sizeof(float));
    for (int i = 0; i < v->dimension; i++){
        n->elements[i] = v->elements[i] / magnitude;
    }

    return n;
}

Vector* reflection_vector(Vector* i, Vector* n) {
    Vector* r = (Vector*) malloc(sizeof(Vector));
    r->dimension = i->dimension;
    r->elements = (float *) malloc(r->dimension * sizeof(float));

    Vector* temp = scale_vector(2.f * dot_product(i, n), n);
    r = subtract_vectors(i, temp);

    free(temp->elements);
    free(temp);

    return r;
}

Matrix* create_matrix(int dimension) {
    Matrix* m = (Matrix*) malloc(sizeof(Matrix));
    m->dimension = dimension;
    m->elements = (float **) malloc(dimension * sizeof(float*));
    for (int i = 0; i < dimension; i++)
        m->elements[i] = (float *) malloc(dimension * sizeof(float));
    return m;
}

void free_matrix(Matrix* m) {
    for (int i = 0; i < m->dimension; i++) {
        free(m->elements[i]);
    }
    free(m->elements);
    free(m);
}

Matrix* add_matrices(Matrix* m1, Matrix* m2) {
    if (m1->dimension != m2->dimension)
        return NULL;
    Matrix* m = create_matrix(m1->dimension);
    for (int i = 0; i < m->dimension; i++)
        for (int j = 0; j < m->dimension; j++)
            m->elements[i][j] = m1->elements[i][j] + m2->elements[i][j];
    return m;
}


Matrix* subtract_matrices(Matrix* m1, Matrix* m2) {
    if (m1->dimension != m2->dimension)
        return NULL;
    Matrix* m = create_matrix(m1->dimension);
    for (int i = 0; i < m->dimension; i++)
        for (int j = 0; j < m->dimension; j++)
            m->elements[i][j] = m1->elements[i][j] - m2->elements[i][j];
    return m;
}

Matrix* multiply_matrices(Matrix* m1, Matrix* m2) {
    if (m1->dimension != m2->dimension)
        return NULL;

    Matrix* result = create_matrix(m1->dimension);
    for (int i = 0; i < m1->dimension; i++) {
        for (int j = 0; j < m1->dimension; j++) {
            result->elements[i][j] = 0;
            for (int k = 0; k < m1->dimension; k++) {
                result->elements[i][j] += m1->elements[i][k] * m2->elements[k][j];
            }
        }
    }
    return result;
}

Matrix* transpose_matrix(Matrix* m) {
    Matrix* t = create_matrix(m->dimension);
    for (int i = 0; i < t->dimension; i++) {
        for (int j = 0; j < t->dimension; j++) {
            t->elements[i][j] = m->elements[j][i];
        }
    }
    return t;
}

Matrix* scale_matrix(float c, Matrix* m) {
    Matrix* s = create_matrix(m->dimension);
    for (int i = 0; i < m->dimension; i++) {
        for (int j = 0; j < m->dimension; j++) {
            s->elements[i][j] = m->elements[i][j] * c;
        }
    }
    return s;
}

Matrix* scale_matrix_anisotropic(float x, float y, float z, Matrix* m) {
    Matrix* s = create_matrix(m->dimension);
    for (int i = 0; i < m->dimension; i++) {
        for (int j = 0; j < m->dimension; j++) {
            if (j == 0) s->elements[i][j] = m->elements[i][j] * x;
            else if (j == 1) s->elements[i][j] = m->elements[i][j] * y;
            else if (j == 2) s->elements[i][j] = m->elements[i][j] * z;
            else if (j == 3) s->elements[i][j] = m->elements[i][j];
        }
    }
    return s;
}

float determinant_matrix(Matrix* m) {
    if (m->dimension == 2)
        return (m->elements[0][0] * m->elements[1][1]) - (m->elements[0][1] * m->elements[1][0]);

    if (m->dimension == 1)
        return m->elements[0][0];
    
    if (m->dimension == 3) {
        return m->elements[0][0] * (m->elements[1][1] * m->elements[2][2] - m->elements[1][2] * m->elements[2][1])
             - m->elements[0][1] * (m->elements[1][0] * m->elements[2][2] - m->elements[1][2] * m->elements[2][0])
             + m->elements[0][2] * (m->elements[1][0] * m->elements[2][1] - m->elements[1][1] * m->elements[2][0]);
    }

    if (m->dimension == 4) {
        return m->elements[0][0] * (
                m->elements[1][1] * (m->elements[2][2] * m->elements[3][3] - m->elements[2][3] * m->elements[3][2])
              - m->elements[1][2] * (m->elements[2][1] * m->elements[3][3] - m->elements[2][3] * m->elements[3][1])
              + m->elements[1][3] * (m->elements[2][1] * m->elements[3][2] - m->elements[2][2] * m->elements[3][1])
            )
          - m->elements[0][1] * (
                m->elements[1][0] * (m->elements[2][2] * m->elements[3][3] - m->elements[2][3] * m->elements[3][2])
              - m->elements[1][2] * (m->elements[2][0] * m->elements[3][3] - m->elements[2][3] * m->elements[3][0])
              + m->elements[1][3] * (m->elements[2][0] * m->elements[3][2] - m->elements[2][2] * m->elements[3][0])
            )
          + m->elements[0][2] * (
                m->elements[1][0] * (m->elements[2][1] * m->elements[3][3] - m->elements[2][3] * m->elements[3][1])
              - m->elements[1][1] * (m->elements[2][0] * m->elements[3][3] - m->elements[2][3] * m->elements[3][0])
              + m->elements[1][3] * (m->elements[2][0] * m->elements[3][1] - m->elements[2][1] * m->elements[3][0])
            )
          - m->elements[0][3] * (
                m->elements[1][0] * (m->elements[2][1] * m->elements[3][2] - m->elements[2][2] * m->elements[3][1])
              - m->elements[1][1] * (m->elements[2][0] * m->elements[3][2] - m->elements[2][2] * m->elements[3][0])
              + m->elements[1][2] * (m->elements[2][0] * m->elements[3][1] - m->elements[2][1] * m->elements[3][0])
            );
    }
    return 0.0f;
}

Matrix* translate_matrix(float x, float y, float z) {
    Matrix* t = create_matrix(4);
    for (int i = 0; i < 4; i++)
        t->elements[i][i] = 1;

    t->elements[3][0] = x;
    t->elements[3][1] = y;
    t->elements[3][2] = z;
    return t;
}

Matrix* translate_matrix_in_place(Matrix* m, float x, float y, float z) {
    if (m->dimension != 4)
        return NULL;

    Vector* t = create_vector(4);
    t->elements[0] = x;
    t->elements[1] = y;
    t->elements[2] = z;
    t->elements[3] = 0;

    for (int i = 0; i < 4; i++) {
        float dot_product = 0;
        for (int j = 0; j < 4; j++) {
            dot_product += m->elements[j][i] * t->elements[j];
        }
        m->elements[3][i] += dot_product;
    }
    free_vector(t);
    return m;
}

Matrix* outer_product(Vector* v1, Vector* v2) {
    Matrix* m = create_matrix(4);

	for (int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
		    m->elements[i][j] = (i < 3 && j < 3) ? v1->elements[i] * v2->elements[j] : 0.f;

    return m; 
}

Matrix* duplicate_matrix(Matrix* m) {
    Matrix* n = create_matrix(m->dimension);
    for (int i = 0; i < m->dimension; i++)
        for (int j = 0; j < m->dimension; j++)
            n->elements[i][j] = m->elements[i][j];
    return n;
}

Matrix* identity_matrix(int dimension) {
    Matrix* m = create_matrix(dimension);
    for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
            m->elements[i][j] = i == j ? 1 : 0;
    return m;
}

Matrix* rotate_matrix(Matrix* m, float x, float y, float z, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    float t = 1 - c;
    
    float mag = sqrt(x*x + y*y + z*z);
    x /= mag;
    y /= mag;
    z /= mag;

    Matrix* r = create_matrix(4);
    
    r->elements[0][0] = t*x*x + c;
    r->elements[0][1] = t*x*y - s*z;
    r->elements[0][2] = t*x*z + s*y;
    r->elements[0][3] = 0;

    r->elements[1][0] = t*x*y + s*z;
    r->elements[1][1] = t*y*y + c;
    r->elements[1][2] = t*y*z - s*x;
    r->elements[1][3] = 0;

    r->elements[2][0] = t*x*z - s*y;
    r->elements[2][1] = t*y*z + s*x;
    r->elements[2][2] = t*z*z + c;
    r->elements[2][3] = 0;

    r->elements[3][0] = 0;
    r->elements[3][1] = 0;
    r->elements[3][2] = 0;
    r->elements[3][3] = 1;

    Matrix* result = multiply_matrices(m, r);
    free_matrix(r);
    
    return result;
}

Matrix* rotate_matrix_X(Matrix* m, float angle) {
    float c = cos(angle);
    float s = sin(angle);

    Matrix* r = create_matrix(4);
    
    r->elements[0][0] = 1;
    r->elements[0][1] = 0;
    r->elements[0][2] = 0;
    r->elements[0][3] = 0;

    r->elements[1][0] = 0;
    r->elements[1][1] = c;
    r->elements[1][2] = -s;
    r->elements[1][3] = 0;

    r->elements[2][0] = 0;
    r->elements[2][1] = s;
    r->elements[2][2] = c;
    r->elements[2][3] = 0;

    r->elements[3][0] = 0;
    r->elements[3][1] = 0;
    r->elements[3][2] = 0;
    r->elements[3][3] = 1;

    Matrix* result = multiply_matrices(m, r);
    free_matrix(r);
    
    return result;
}

Matrix* rotate_matrix_Y(Matrix* m, float angle) {
    float c = cos(angle);
    float s = sin(angle);

    Matrix* r = create_matrix(4);
    
    r->elements[0][0] = c;
    r->elements[0][1] = 0;
    r->elements[0][2] = s;
    r->elements[0][3] = 0;

    r->elements[1][0] = 0;
    r->elements[1][1] = 1;
    r->elements[1][2] = 0;
    r->elements[1][3] = 0;

    r->elements[2][0] = -s;
    r->elements[2][1] = 0;
    r->elements[2][2] = c;
    r->elements[2][3] = 0;

    r->elements[3][0] = 0;
    r->elements[3][1] = 0;
    r->elements[3][2] = 0;
    r->elements[3][3] = 1;

    Matrix* result = multiply_matrices(m, r);
    free_matrix(r);
    
    return result;
}

Matrix* rotate_matrix_Z(Matrix* m, float angle) {
    float c = cos(angle);
    float s = sin(angle);

    Matrix* r = create_matrix(4);
    
    r->elements[0][0] = c;
    r->elements[0][1] = -s;
    r->elements[0][2] = 0;
    r->elements[0][3] = 0;

    r->elements[1][0] = s;
    r->elements[1][1] = c;
    r->elements[1][2] = 0;
    r->elements[1][3] = 0;

    r->elements[2][0] = 0;
    r->elements[2][1] = 0;
    r->elements[2][2] = 1;
    r->elements[2][3] = 0;

    r->elements[3][0] = 0;
    r->elements[3][1] = 0;
    r->elements[3][2] = 0;
    r->elements[3][3] = 1;

    Matrix* result = multiply_matrices(m, r);
    free_matrix(r);
    
    return result;
}

Matrix* invert_matrix(Matrix* m) {
    if (m->dimension != 4) {
        return NULL;  // This function only works for 4x4 matrices
    }

    float s[6];
    float c[6];

    s[0] = m->elements[0][0] * m->elements[1][1] - m->elements[1][0] * m->elements[0][1];
    s[1] = m->elements[0][0] * m->elements[1][2] - m->elements[1][0] * m->elements[0][2];
    s[2] = m->elements[0][0] * m->elements[1][3] - m->elements[1][0] * m->elements[0][3];
    s[3] = m->elements[0][1] * m->elements[1][2] - m->elements[1][1] * m->elements[0][2];
    s[4] = m->elements[0][1] * m->elements[1][3] - m->elements[1][1] * m->elements[0][3];
    s[5] = m->elements[0][2] * m->elements[1][3] - m->elements[1][2] * m->elements[0][3];

    c[0] = m->elements[2][0] * m->elements[3][1] - m->elements[3][0] * m->elements[2][1];
    c[1] = m->elements[2][0] * m->elements[3][2] - m->elements[3][0] * m->elements[2][2];
    c[2] = m->elements[2][0] * m->elements[3][3] - m->elements[3][0] * m->elements[2][3];
    c[3] = m->elements[2][1] * m->elements[3][2] - m->elements[3][1] * m->elements[2][2];
    c[4] = m->elements[2][1] * m->elements[3][3] - m->elements[3][1] * m->elements[2][3];
    c[5] = m->elements[2][2] * m->elements[3][3] - m->elements[3][2] * m->elements[2][3];

    float det = s[0]*c[5] - s[1]*c[4] + s[2]*c[3] + s[3]*c[2] - s[4]*c[1] + s[5]*c[0];

    if (fabs(det) < 1e-6) {
        return NULL;  // Matrix is not invertible
    }

    float invdet = 1.0f / det;

    Matrix* inv = create_matrix(4);

    inv->elements[0][0] = ( m->elements[1][1] * c[5] - m->elements[1][2] * c[4] + m->elements[1][3] * c[3]) * invdet;
    inv->elements[0][1] = (-m->elements[0][1] * c[5] + m->elements[0][2] * c[4] - m->elements[0][3] * c[3]) * invdet;
    inv->elements[0][2] = ( m->elements[3][1] * s[5] - m->elements[3][2] * s[4] + m->elements[3][3] * s[3]) * invdet;
    inv->elements[0][3] = (-m->elements[2][1] * s[5] + m->elements[2][2] * s[4] - m->elements[2][3] * s[3]) * invdet;

    inv->elements[1][0] = (-m->elements[1][0] * c[5] + m->elements[1][2] * c[2] - m->elements[1][3] * c[1]) * invdet;
    inv->elements[1][1] = ( m->elements[0][0] * c[5] - m->elements[0][2] * c[2] + m->elements[0][3] * c[1]) * invdet;
    inv->elements[1][2] = (-m->elements[3][0] * s[5] + m->elements[3][2] * s[2] - m->elements[3][3] * s[1]) * invdet;
    inv->elements[1][3] = ( m->elements[2][0] * s[5] - m->elements[2][2] * s[2] + m->elements[2][3] * s[1]) * invdet;

    inv->elements[2][0] = ( m->elements[1][0] * c[4] - m->elements[1][1] * c[2] + m->elements[1][3] * c[0]) * invdet;
    inv->elements[2][1] = (-m->elements[0][0] * c[4] + m->elements[0][1] * c[2] - m->elements[0][3] * c[0]) * invdet;
    inv->elements[2][2] = ( m->elements[3][0] * s[4] - m->elements[3][1] * s[2] + m->elements[3][3] * s[0]) * invdet;
    inv->elements[2][3] = (-m->elements[2][0] * s[4] + m->elements[2][1] * s[2] - m->elements[2][3] * s[0]) * invdet;

    inv->elements[3][0] = (-m->elements[1][0] * c[3] + m->elements[1][1] * c[1] - m->elements[1][2] * c[0]) * invdet;
    inv->elements[3][1] = ( m->elements[0][0] * c[3] - m->elements[0][1] * c[1] + m->elements[0][2] * c[0]) * invdet;
    inv->elements[3][2] = (-m->elements[3][0] * s[3] + m->elements[3][1] * s[1] - m->elements[3][2] * s[0]) * invdet;
    inv->elements[3][3] = ( m->elements[2][0] * s[3] - m->elements[2][1] * s[1] + m->elements[2][2] * s[0]) * invdet;

    return inv;
}

Matrix* look_at_matrix(Vector* eye, Vector* center, Vector* up) {
    Vector* temp = subtract_vectors(center, eye);

    Vector* f = normalize_vector(temp);
    Vector* u = normalize_vector(up);
    Vector* s = normalize_vector(cross_product(f, u));

    free_vector(u);
    free_vector(temp);
    u = cross_product(s, f);

    Matrix* result = create_matrix(4);
    result->elements[0][0] = s->elements[0];
    result->elements[0][1] = s->elements[1];
    result->elements[0][2] = s->elements[2];
    result->elements[0][3] = 0;

    result->elements[1][0] = u->elements[0];
    result->elements[1][1] = u->elements[1];
    result->elements[1][2] = u->elements[2];
    result->elements[1][3] = 0;

    result->elements[2][0] = -f->elements[0];
    result->elements[2][1] = -f->elements[1];
    result->elements[2][2] = -f->elements[2];
    result->elements[2][3] = 0;

    result->elements[3][0] = 0;
    result->elements[3][1] = 0;
    result->elements[3][2] = 0;
    result->elements[3][3] = 1;
    Matrix* translation = translate_matrix(-eye->elements[0], -eye->elements[1], -eye->elements[2]);
    Matrix* final_result = multiply_matrices(result, translation);

    free_matrix(result); free_matrix(translation);
    free_vector(f); free_vector(u); free_vector(s);

    return final_result;
}

Matrix* perspective_projection_matrix(float fov, float aspect, float near, float far) {
    Matrix* result = create_matrix(4);
    float t = tan(fov / 2.0f);

    result->elements[0][0] = 1.0f / (aspect * t);
    result->elements[1][1] = 1.0f / t;
    result->elements[2][2] = -(far + near) / (far - near);
    result->elements[2][3] = -1.0f;
    result->elements[3][2] = -(2.0f * far * near) / (far - near);
    result->elements[3][3] = 0.0f;

    return result;
}

Matrix* orthographic_projection_matrix(float left, float right, float bottom, float top, float near, float far) {
    Matrix* result = create_matrix(4);

    result->elements[0][0] = 2.0f / (right - left);
    result->elements[1][1] = 2.0f / (top - bottom);
    result->elements[2][2] = -2.0f / (far - near);
    result->elements[3][0] = -(right + left) / (right - left);
    result->elements[3][1] = -(top + bottom) / (top - bottom);
    result->elements[3][2] = -(far + near) / (far - near);
    result->elements[3][3] = 1.0f;

    return result;
}


Quaternion* multiply_quaternions(Quaternion* q1, Quaternion* q2) {
    Quaternion* result = create_quaternion(
        q1->w * q2->x + q1->x * q2->w + q1->y * q2->z - q1->z * q2->y,
        q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x,
        q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w,
        q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z
    );
    return result;
}

Quaternion* rotate_quaternion(Quaternion* q, Vector* axis, float angle) {
    float half_angle = angle / 2.0f;
    float sin_half = sin(half_angle);
    
    Quaternion* rotation = create_quaternion(
        axis->elements[0] * sin_half,
        axis->elements[1] * sin_half,
        axis->elements[2] * sin_half,
        cos(half_angle)
    );
    
    Quaternion* rotated = multiply_quaternions(rotation, q);
    Quaternion* conjugate = quaternion_conjugate(rotation);
    Quaternion* result = multiply_quaternions(rotated, conjugate);
    
    free_quaternion(rotation);
    free_quaternion(rotated);
    free_quaternion(conjugate);
    
    return result;
}

Quaternion* quaternion_conjugate(Quaternion* q) {
    return create_quaternion(-q->x, -q->y, -q->z, q->w);
}

Matrix* matrix_from_quaternion(Quaternion* q) {
    Matrix* m = create_matrix(4);
    
    float xx = q->x * q->x;
    float xy = q->x * q->y;
    float xz = q->x * q->z;
    float xw = q->x * q->w;
    
    float yy = q->y * q->y;
    float yz = q->y * q->z;
    float yw = q->y * q->w;
    
    float zz = q->z * q->z;
    float zw = q->z * q->w;
    
    m->elements[0][0] = 1 - 2 * (yy + zz);
    m->elements[0][1] = 2 * (xy - zw);
    m->elements[0][2] = 2 * (xz + yw);
    m->elements[0][3] = 0;
    
    m->elements[1][0] = 2 * (xy + zw);
    m->elements[1][1] = 1 - 2 * (xx + zz);
    m->elements[1][2] = 2 * (yz - xw);
    m->elements[1][3] = 0;
    
    m->elements[2][0] = 2 * (xz - yw);
    m->elements[2][1] = 2 * (yz + xw);
    m->elements[2][2] = 1 - 2 * (xx + yy);
    m->elements[2][3] = 0;
    
    m->elements[3][0] = 0;
    m->elements[3][1] = 0;
    m->elements[3][2] = 0;
    m->elements[3][3] = 1;
    
    return m;
}

Quaternion* create_quaternion(float x, float y, float z, float w) {
    Quaternion* q = (Quaternion*) (malloc(sizeof(Quaternion)));
    q->x = x;
    q->y = y;
    q->z = z;
    q->w = w;
    return q;
}

Quaternion* identity_quaternion() {
    return create_quaternion(0, 0, 0, 1);
}

void free_quaternion(Quaternion* q) {
    free(q);
    return;
}

Quaternion* add_quaternions(Quaternion* q1, Quaternion* q2) {
    Quaternion* q = (Quaternion*) (malloc(sizeof(Quaternion)));
    q->x = q1->x + q2->x;
    q->y = q1->y + q2->y;
    q->z = q1->z + q2->z;
    q->w = q1->w + q2->w;
    return q;
}

Quaternion* subtract_quaternions(Quaternion* q1, Quaternion* q2) {
    Quaternion* q = (Quaternion*) (malloc(sizeof(Quaternion)));
    q->x = q1->x - q2->x;
    q->y = q1->y - q2->y;
    q->z = q1->z - q2->z;
    q->w = q1->w - q2->w;
    return q;
}

float dot_product_quaternion(Quaternion* q1, Quaternion* q2) {
    float sum = 0;
    sum += q1->x * q2->x;
    sum += q1->y * q2->y;
    sum += q1->z * q2->z;
    sum += q1->w * q2->w;
    return sum;
}

Quaternion* scale_quaternion(float c, Quaternion* q) {
    Quaternion* s = (Quaternion*) (malloc(sizeof(Quaternion)));
    s->x = q->x * c;
    s->y = q->y * c;
    s->z = q->z * c;
    s->w = q->w * c;
    return s;
}

Quaternion* normalize_quaternion(Quaternion* q) {
    float magnitude = sqrt(dot_product_quaternion(q, q));
    Quaternion* n = (Quaternion*) (malloc(sizeof(Quaternion)));
    n->x = q->x / magnitude;
    n->y = q->y / magnitude;
    n->z = q->z / magnitude;
    n->w = q->w / magnitude;
    return n;
}


#endif