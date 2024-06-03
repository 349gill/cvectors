#include <math.h>
#include <string.h>

typedef struct {
    int dimension;
    float *elements;
} Vector;

typedef struct {
    int dimension;
    float **elements;
} Matrix;

Vector add_vectors(Vector v1, Vector v2);
Vector subtract_vectors(Vector v1, Vector v2);
Vector dot_product(Vector v1, Vector v2);
Vector cross_product(Vector v1, Vector v2);
Vector transform_vector(Matrix m, Vector v);
Vector scale_vector(float c, Vector v);
Vector magnitude_vector(Vector v);
Vector normalize_vector(Vector v);
Vector duplicate_vector(Vector v1, Vector v2);
Vector minimum_vector(Vector v1, Vector v2);
Vector maximum_vector(Vector v1, Vector v2);
Vector reflected_vector(Vector v1, Vector v2);

Matrix add_matrices(Matrix m1, Matrix m2);
Matrix subtract_matrices(Matrix m1, Matrix m2);
Matrix multiply_matrices(Matrix m1, Matrix m2);
Matrix transpose_matrix(Matrix m);
Matrix scale_matrix(float c, Matrix m);
Matrix determinant_matrix(Matrix m);

// To implement: Rotating matrices, Translating matrices, Inverting matrices, Orthonormalization, Eigen stuff
// Gram–Schmidt process for orthogonal matrices, Vector projections, Frustum projections, Quaternions
// Matrix raised to some power. Sine, Cosine, and Tangent of a Matrix
