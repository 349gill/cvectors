#include <math.h>
#include <stdlib.h>

typedef struct {
    int dimension;
    float *elements;
} Vector;

// Only 1x1, 2x2, 3x3, and 4x4 Matrices are supported. In the case of 
// Square matrices with dimension > 4, all defined operations will work as 
// expected except Matrix Multiplication and computation of Determinants.
typedef struct {
    int dimension;
    float **elements;
} Matrix;

Vector* add_vectors(Vector v1, Vector v2);
Vector* subtract_vectors(Vector v1, Vector v2);
float dot_product(Vector v1, Vector v2);
Vector* cross_product(Vector v1, Vector v2);
Vector* transform_vector(Matrix m, Vector v);
Vector* scale_vector(float c, Vector v);
float magnitude_vector(Vector v);
Vector* normalize_vector(Vector v);

Matrix* create_matrix(int dimension);
void free_matrix(Matrix* m);
Matrix* add_matrices(Matrix m1, Matrix m2);
Matrix* subtract_matrices(Matrix m1, Matrix m2);
Matrix* multiply_matrices(Matrix m1, Matrix m2);
Matrix* transpose_matrix(Matrix m);
Matrix* scale_matrix(float c, Matrix m);
float determinant_matrix(Matrix m);

// To implement: Rotating matrices, Translating matrices, Inverting matrices, Orthonormalization, Eigen stuff
// Gram–Schmidt process for orthogonal matrices, Vector projections, Frustum projections, Quaternions
// Matrix raised to some power. Sine, Cosine, and Tangent of a Matrix
// Maximum and Minimum vectors, reflected vectors

Vector* add_vectors(Vector v1, Vector v2) {
    if (v1.dimension != v2.dimension) {
        return NULL;
    }
    Vector v;
    v.dimension = v1.dimension;
    v.elements = (float *) malloc(v.dimension * sizeof(float));
    for (int i = 0; i < v.dimension; i++) {
        v.elements[i] = v1.elements[i];
        v.elements[i] += v2.elements[i];
    }
    return &v1;
}

Vector* subtract_vectors(Vector v1, Vector v2) {
    if (v1.dimension != v2.dimension) {
        return NULL;
    }
    Vector v;
    v.dimension = v1.dimension;
    v.elements = (float *) malloc(v.dimension * sizeof(float));
    for (int i = 0; i < v.dimension; i++) {
        v.elements[i] = v1.elements[i];
        v.elements[i] -= v2.elements[i];
    }
    return &v1;
}

float dot_product(Vector v1, Vector v2) {
    if (v1.dimension != v2.dimension) {
        return NULL;
    }
    int x = 0;
    for (int i = 0; i < v1.dimension; i++) {
        x += v1.elements[i] * v2.elements[i];
    }
    return x;
}

Vector* cross_product(Vector v1, Vector v2) {
    if ((v1.dimension != v2.dimension) && (v1.dimension != 3)) {
        return NULL;
    }

    Vector v;
    v.dimension = v1.dimension;
    v.elements = (float *) malloc(v.dimension * sizeof(float));
    v.elements[0] = (v1.elements[1] * v2.elements[2]) - (v1.elements[2] * v2.elements[1]);
    v.elements[1] =  - ((v1.elements[0] * v2.elements[2]) - (v1.elements[2] * v2.elements[0]));
    v.elements[2] =  (v1.elements[0] * v2.elements[1]) - (v1.elements[1] * v2.elements[0]);
    return &v;
}


Vector* transform_vector(Matrix m, Vector v) {
    if (m.dimension != v.dimension) {
        return NULL;
    }
    Vector t;
    t.dimension = v.dimension;
    t.elements = (float *) malloc(v.dimension * sizeof(float));
    for (int i = 0; i < t.dimension; i++) {
        t.elements[i] = 0;
        for (int j = 0; j < t.dimension; j++) {
            t.elements[i] += (m.elements[i][j]) * (v.elements[j]);
        }
    }
    return &t;
}

Vector* scale_vector(float c, Vector v) {
    Vector s;
    s.dimension = v.dimension;
    s.elements = (float *) malloc(v.dimension * sizeof(float));
    for (int i = 0; i < v.dimension; i++){
        s.elements[i] = v.elements[i] * c;
    }
    return &s;
}

float magnitude_vector(Vector v) {
    float x = 0;
    for (int i = 0; i < v.dimension; i++){
        x += pow(v.elements[i], 2);
    }

    return pow(x, 0.5);
}

Vector* normalize_vector(Vector v) {
    float magnitude = magnitude_vector(v);
    
    Vector n;
    n.dimension = v.dimension;
    n.elements = (float *) malloc(v.dimension * sizeof(float));
    for (int i = 0; i < v.dimension; i++){
        n.elements[i] = v.elements[i] / magnitude;
    }

    return &n;
}

Matrix* create_matrix(int dimension) {
    Matrix m;
    m.dimension = dimension;
    m.elements = (float **) malloc(dimension * dimension * sizeof(float));
    return &m;
}

void free_matrix(Matrix* m) {
    free(m->elements);
    return;
}

Matrix* add_matrices(Matrix m1, Matrix m2) {
    if (m1.dimension != m2.dimension) {
        return NULL;
    }
    Matrix m;
    m.dimension = m1.dimension;
    m.elements = (float **) malloc(m.dimension * m.dimension * sizeof(float));
    for (int i = 0; i < m.dimension; i++) {
        for (int j = 0; j < m.dimension; j++) {
            m.elements[i][j] = m1.elements[i][j] + m2.elements[i][j];
        }
    }

    return &m;
}


Matrix* subtract_matrices(Matrix m1, Matrix m2) {
    if (m1.dimension != m2.dimension) {
        return NULL;
    }
    Matrix m;
    m.dimension = m1.dimension;
    m.elements = (float **) malloc(m.dimension * m.dimension * sizeof(float));
    for (int i = 0; i < m.dimension; i++) {
        for (int j = 0; j < m.dimension; j++) {
            m.elements[i][j] = m1.elements[i][j] - m2.elements[i][j];
        }
    }

    return &m;
}

// Strassen’s Matrix Multiplication
Matrix* multiply_matrices(Matrix m1, Matrix m2) {
    if (m1.dimension == m2.dimension) {
        if (m1.dimension == 3) {
            Matrix m1_4d, m2_4d;
            m1_4d.dimension = 4; m2_4d.dimension = 4;
            m1_4d.elements = (float **) malloc(16 * sizeof(float)); m2_4d.elements = (float **) malloc(16 * sizeof(float));
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    if ((i == 4) || (j == 4)) {
                        m1_4d.elements[i][j] = 0;
                        m2_4d.elements[i][j] = 0;
                    }
                    else {
                        m1_4d.elements[i][j] = m1.elements[i][j];
                        m2_4d.elements[i][j] = m2.elements[i][j];
                    }
                }
            }
            Matrix m_4d = *multiply_matrices(m1_4d, m2_4d);

            Matrix m;
            m.dimension = 3;
            m.elements = (float **) malloc(9 * sizeof(float));
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m.elements[i][j] = m_4d.elements[i][j];
                }
            }
            free(m1_4d.elements); free(m2_4d.elements); free(m_4d.elements);
            return &m;

        }
        Matrix m;
        m.dimension = m1.dimension;
        m.elements =  (float **) malloc(m.dimension * m.dimension * sizeof(float));
        if (m1.dimension == 1) {
            m.elements[0][0] = m1.elements[0][0] * m2.elements[0][0];
        }
        else if (m1.dimension == 2) {
            int mult1, mult2, mult3, mult4 , mult5, mult6, mult7;
            mult1 = (m1.elements[0][0] + m1.elements[1][1]) * (m2.elements[0][0] + m2.elements[1][1]);
            mult2 = (m1.elements[1][0] + m1.elements[1][1]) * m2.elements[0][0];
            mult3 = m1.elements[0][0] * (m2.elements[0][1] - m2.elements[1][1]);
            mult4 = m1.elements[1][1] * (m2.elements[1][0] - m2.elements[0][0]);
            mult5 = (m1.elements[0][0] + m1.elements[0][1]) * m2.elements[1][1];
            mult6 = (m1.elements[1][0] - m1.elements[0][0]) * (m2.elements[0][0] + m2.elements[0][1]);
            mult7 = (m1.elements[0][1] - m1.elements[1][1]) * (m2.elements[1][0] + m2.elements[1][1]);
            m.elements[0][0] = mult1 + mult4 - mult5 + mult7;
            m.elements[0][1] = mult3 + mult5;
            m.elements[1][0] = mult2 + mult4;
            m.elements[1][1] = mult1 - mult2 + mult3 + mult6;
        }
        else if (m1.dimension == 4) {
            Matrix mult1, mult2, mult3, mult4 , mult5, mult6, mult7;
            Matrix A1, A2, A3, A4, B1, B2, B3, B4;

            A1.dimension = 2, A2.dimension = 2, A3.dimension = 2, A4.dimension = 2;
            B1.dimension = 2, B2.dimension = 2, B3.dimension = 2, B4.dimension = 2;

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    A1.elements[i][j] = m1.elements[i][j];
                    A2.elements[i][j] = m1.elements[i][j + 2];
                    A3.elements[i][j] = m1.elements[i + 2][j];
                    A4.elements[i][j] = m1.elements[i + 2][j + 2];
                    B1.elements[i][j] = m2.elements[i][j];
                    B2.elements[i][j] = m2.elements[i][j + 2];
                    B3.elements[i][j] = m2.elements[i + 2][j];
                    B4.elements[i][j] = m2.elements[i + 2][j + 2];
                }
            }

            Matrix *M1 = multiply_matrices(*add_matrices(A1, A4), *add_matrices(B1, B4));
            Matrix *M2 = multiply_matrices(*add_matrices(A3, A4), B1);
            Matrix *M3 = multiply_matrices(A1, *subtract_matrices(B2, B4));
            Matrix *M4 = multiply_matrices(A4, *subtract_matrices(B3, B1));
            Matrix *M5 = multiply_matrices(*add_matrices(A1, A2), B4);
            Matrix *M6 = multiply_matrices(*subtract_matrices(A3, A1), *add_matrices(B1, B2));
            Matrix *M7 = multiply_matrices(*subtract_matrices(A2, A4), *add_matrices(B3, B4));

            Matrix *C11 = add_matrices(*subtract_matrices(*add_matrices(*M1, *M4), *M5), *M7);
            Matrix *C12 = add_matrices(*M3, *M5);
            Matrix *C21 = add_matrices(*M2, *M4);
            Matrix *C22 = add_matrices(*subtract_matrices(*add_matrices(*M1, *M3), *M2), *M6);

            Matrix *C = create_matrix(m1.dimension);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    C->elements[i][j] = C11->elements[i][j];
                    C->elements[i][j + 2] = C12->elements[i][j];
                    C->elements[i + 2][j] = C21->elements[i][j];
                    C->elements[i + 2][j + 2] = C22->elements[i][j];
                }
            }

            free_matrix(&A1); free_matrix(&A2); free_matrix(&A3); free_matrix(&A4);
            free_matrix(&B1); free_matrix(&B2); free_matrix(&B3); free_matrix(&B4);
            free_matrix(M1); free_matrix(M2); free_matrix(M3); free_matrix(M4); free_matrix(M5); free_matrix(M6); free_matrix(M7);
            free_matrix(C11); free_matrix(C12); free_matrix(C21); free_matrix(C22);

            return C;
        }
        else {
            free(m.elements);
            return NULL;
        }
        return &m;
    }
    return NULL;
}

Matrix* transpose_matrix(Matrix m) {
    Matrix t;
    t.dimension = m.dimension;
    t.elements = (float **) malloc(m.dimension * m.dimension * sizeof(float));
    for (int i = 0; i < m.dimension; i++) {
        for (int j = 0; j < m.dimension; j++) {
            t.elements[i][j] = m.elements[j][i];
        }
    }
    return &t;
}

Matrix* scale_matrix(float c, Matrix m) {
    Matrix s;
    s.dimension = m.dimension;
    s.elements = (float **) malloc(m.dimension * m.dimension * sizeof(float));
    for (int i = 0; i < m.dimension; i++) {
        for (int j = 0; j < m.dimension; j++) {
            s.elements[i][j] = m.elements[i][j] * c;
        }
    }
    return &m;
}

float determinant_matrix(Matrix m) {

    if (m.dimension == 2) {
        return (m.elements[0][0] * m.elements[2][2]) - (m.elements[1][2] * m.elements[2][1]);
    }
    else if (m.dimension == 1) {
        return m.elements[0][0];
    }
    else if (m.dimension == 3 || m.dimension == 4) {
        float determinant; int factor = 1;
        Matrix *new_matrix;
        for(int i = 0; i < m.dimension; i++) 
        {
            new_matrix->elements = (float **) malloc((m.dimension-1) * (m.dimension - 1) * sizeof(float));
            // Skip first row
            for(int j = 1; j < m.dimension; j++) 
            {
                for (int k = 0; k < m.dimension; k++)
                {
                    if (k == i) {
                        continue;
                    };
                    new_matrix->elements[j-1][k < i ? k :(k - 1)] = m.elements[j][k];
                }
            }
            determinant += factor * m.elements[0][i] * determinant_matrix(*new_matrix);
            factor *= -1;
            free(new_matrix);
        }
        return determinant;
    }

    return 0;
}
