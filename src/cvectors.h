/*
 *  cvectors.h
 *  
 *  First Commit: May 27, 2024
 *      Author: Harsh Gill
 * 
 *  Linear Algebra Library Supporting upto 4-dimensional
 *  Vectors and Matrices. Intended for Computer Graphics purposes.
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

Vector* add_vectors(Vector* v1, Vector* v2);
Vector* subtract_vectors(Vector* v1, Vector* v2);
float dot_product(Vector* v1, Vector* v2);
Vector* cross_product(Vector* v1, Vector* v2);
Vector* transform_vector(Matrix* m, Vector* v);
Vector* scale_vector(float c, Vector* v);
float magnitude_vector(Vector* v);
Vector* normalize_vector(Vector* v);
Vector* reflection_vector(Vector* i, Vector* n);

Matrix* create_matrix(int dimension);
void free_matrix(Matrix* m);
Matrix* add_matrices(Matrix* m1, Matrix* m2);
Matrix* subtract_matrices(Matrix* m1, Matrix* m2);
Matrix* multiply_matrices(Matrix* m1, Matrix* m2);
Matrix* transpose_matrix(Matrix* m);
Matrix* scale_matrix(float c, Matrix* m);
Matrix* scale_matrix_anisotropic(float x, float y, float z, Matrix* m);
float determinant_matrix(Matrix* m);
Matrix* translate_matrix(Matrix* m, float x, float y, float z);
Matrix* translate_matrix_in_place(Matrix* m, float x, float y, float z); // TO DO
Matrix* outer_product(Vector* v1, Vector* v2);

// TO DO
Matrix* rotate_matrix(Matrix* m, float x, float y, float z, float angle);
Matrix* rotate_matrix_X(Matrix* m, float angle);
Matrix* rotate_matrix_Y(Matrix* m, float angle);
Matrix* rotate_matrix_Z(Matrix* m, float angle);

// Rotating matrices, Translating matrices, Inverting matrices, Orthonormalization, Eigen stuff
// Gram–Schmidt process for orthogonal matrices, Vector projections, Frustum projections, Quaternions
// Matrix raised to some power. Sine, Cosine, and Tangent of a Matrix
// Maximum and Minimum vectors, Identity Matrix, Duplicate Vector/Matrix;

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

Matrix* translate_matrix(Matrix* m, float x, float y, float z) {
    Matrix* t = create_matrix(4);
    for (int i = 0; i < 4; i++)
        t->elements[i][i] = 1;

    t->elements[3][0] = x;
    t->elements[3][1] = y;
    t->elements[3][2] = z;
    return t;
}

Matrix* translate_matrix_in_place(Matrix* m, float x, float y, float z) {
    return NULL;
}

Matrix* outer_product(Vector* v1, Vector* v2) {
    Matrix* m = create_matrix(4);

	for (int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
		    m->elements[i][j] = (i < 3 && j < 3) ? v1->elements[i] * v2->elements[j] : 0.f;

    return m; 
}

Matrix* invert_matrix(Matrix* m) {
	return NULL;
}

#endif