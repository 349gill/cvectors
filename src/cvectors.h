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

Matrix* create_matrix(int dimension);
void free_matrix(Matrix* m);
Matrix* add_matrices(Matrix* m1, Matrix* m2);
Matrix* subtract_matrices(Matrix* m1, Matrix* m2);
Matrix* multiply_matrices(Matrix* m1, Matrix* m2);
Matrix* transpose_matrix(Matrix* m);
Matrix* scale_matrix(float c, Matrix* m);
float determinant_matrix(Matrix* m);

// Matrix Helper Functions
Matrix* strassen(Matrix* A, Matrix* B, int dim);
Matrix* pad_matrix(Matrix* m, int new_dim);
Matrix* unpad_matrix(Matrix* m, int original_dim);
Matrix* determinant_helper(Matrix* m, int exclude_row, int exclude_col);


// Rotating matrices, Translating matrices, Inverting matrices, Orthonormalization, Eigen stuff
// Gram–Schmidt process for orthogonal matrices, Vector projections, Frustum projections, Quaternions
// Matrix raised to some power. Sine, Cosine, and Tangent of a Matrix
// Maximum and Minimum vectors, reflected vectors

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
    if (v1->dimension != 3 || v2->dimension != 3)
        return NULL;

    Vector* v = (Vector*) malloc(sizeof(Vector));
    v->dimension = v1->dimension;
    v->elements = (float *) malloc(v->dimension * sizeof(float));
    v->elements[0] = (v1->elements[1] * v2->elements[2]) - (v1->elements[2] * v2->elements[1]);
    v->elements[1] =  - ((v1->elements[0] * v2->elements[2]) - (v1->elements[2] * v2->elements[0]));
    v->elements[2] =  (v1->elements[0] * v2->elements[1]) - (v1->elements[1] * v2->elements[0]);
    return v;
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

// Strassen’s Matrix Multiplication
Matrix* pad_matrix(Matrix* m, int new_dim) {
    Matrix* padded = create_matrix(new_dim);
    for (int i = 0; i < new_dim; i++) {
        for (int j = 0; j < new_dim; j++) {
            if (i < m->dimension && j < m->dimension) {
                padded->elements[i][j] = m->elements[i][j];
            } else {
                padded->elements[i][j] = 0;
            }
        }
    }
    return padded;
}

Matrix* unpad_matrix(Matrix* m, int original_dim) {
    Matrix* unpadded = create_matrix(original_dim);
    for (int i = 0; i < original_dim; i++) {
        for (int j = 0; j < original_dim; j++) {
            unpadded->elements[i][j] = m->elements[i][j];
        }
    }
    return unpadded;
}

Matrix* strassen(Matrix* A, Matrix* B, int dim) {
    if (dim == 1) {
        Matrix* C = create_matrix(1);
        C->elements[0][0] = A->elements[0][0] * B->elements[0][0];
        return C;
    }

    int new_dim = dim / 2;
    Matrix *a = create_matrix(new_dim);
    Matrix *b = create_matrix(new_dim);
    Matrix *c = create_matrix(new_dim);
    Matrix *d = create_matrix(new_dim);
    Matrix *e = create_matrix(new_dim);
    Matrix *f = create_matrix(new_dim);
    Matrix *g = create_matrix(new_dim);
    Matrix *h = create_matrix(new_dim);

    for (int i = 0; i < new_dim; i++) {
        for (int j = 0; j < new_dim; j++) {
            a->elements[i][j] = A->elements[i][j];
            b->elements[i][j] = A->elements[i][j + new_dim];
            c->elements[i][j] = A->elements[i + new_dim][j];
            d->elements[i][j] = A->elements[i + new_dim][j + new_dim];
            e->elements[i][j] = B->elements[i][j];
            f->elements[i][j] = B->elements[i][j + new_dim];
            g->elements[i][j] = B->elements[i + new_dim][j];
            h->elements[i][j] = B->elements[i + new_dim][j + new_dim];
        }
    }

    Matrix *p1 = strassen(a, subtract_matrices(f, h), new_dim);
    Matrix *p2 = strassen(add_matrices(a, b), h, new_dim);
    Matrix *p3 = strassen(add_matrices(c, d), e, new_dim);
    Matrix *p4 = strassen(d, subtract_matrices(g, e), new_dim);
    Matrix *p5 = strassen(add_matrices(a, d), add_matrices(e, h), new_dim);
    Matrix *p6 = strassen(subtract_matrices(b, d), add_matrices(g, h), new_dim);
    Matrix *p7 = strassen(subtract_matrices(a, c), add_matrices(e, f), new_dim);

    Matrix *C11 = add_matrices(subtract_matrices(add_matrices(p5, p4), p2), p6);
    Matrix *C12 = add_matrices(p1, p2);
    Matrix *C21 = add_matrices(p3, p4);
    Matrix *C22 = add_matrices(subtract_matrices(add_matrices(p5, p1), p3), p7);

    Matrix *C = create_matrix(dim);
    for (int i = 0; i < new_dim; i++) {
        for (int j = 0; j < new_dim; j++) {
            C->elements[i][j] = C11->elements[i][j];
            C->elements[i][j + new_dim] = C12->elements[i][j];
            C->elements[i + new_dim][j] = C21->elements[i][j];
            C->elements[i + new_dim][j + new_dim] = C22->elements[i][j];
        }
    }

    free_matrix(a); free_matrix(b); free_matrix(c); free_matrix(d);
    free_matrix(e); free_matrix(f); free_matrix(g); free_matrix(h);
    free_matrix(p1); free_matrix(p2); free_matrix(p3); free_matrix(p4);
    free_matrix(p5); free_matrix(p6); free_matrix(p7);
    free_matrix(C11); free_matrix(C12); free_matrix(C21); free_matrix(C22);

    return C;
}

Matrix* multiply_matrices(Matrix* m1, Matrix* m2) {
    int original_dim = m1->dimension;
    int new_dim = pow(2, ceil(log2(original_dim)));

    Matrix* padded_m1 = pad_matrix(m1, new_dim);
    Matrix* padded_m2 = pad_matrix(m2, new_dim);

    Matrix* padded_result = strassen(padded_m1, padded_m2, new_dim);

    Matrix* result = unpad_matrix(padded_result, original_dim);

    free_matrix(padded_m1);
    free_matrix(padded_m2);
    free_matrix(padded_result);

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

float determinant_matrix(Matrix* m) {

    if (m->dimension == 2)
        return (m->elements[0][0] * m->elements[1][1]) - (m->elements[0][1] * m->elements[1][0]);

    if (m->dimension == 1)
        return m->elements[0][0];

    int factor = 1;
    float determinant = 0;

    for (int i = 0; i < m->dimension; i++) {
        Matrix* submatrix = determinant_helper(m, 0, i);
        determinant += factor * m->elements[0][i] * determinant_matrix(submatrix);
        factor *= -1;
        free_matrix(submatrix);
    }
    return determinant;
}

Matrix* determinant_helper(Matrix* m, int exclude_row, int exclude_col) {
    int new_dim = m->dimension - 1;
    Matrix* sub = create_matrix(new_dim);
    for (int i = 0, subi = 0; i < m->dimension; i++) {
        if (i == exclude_row) continue;
        for (int j = 0, subj = 0; j < m->dimension; j++) {
            if (j == exclude_col) continue;
            sub->elements[subi][subj] = m->elements[i][j];
            subj++;
        }
        subi++;
    }
    return sub;
}



#endif