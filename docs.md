# CVector Library Documentation

## Vector Functions

### create_vector(int dimension)
Creates a vector with the specified dimension.

### create_3D_vector(float x, float y, float z)
Creates a 3D vector with the given x, y, and z components.

### create_4D_vector(float x, float y, float z, float w)
Creates a 4D vector with the given x, y, z, and w components.

### free_vector(Vector* v)
Frees the memory allocated for the vector.

### duplicate_vector(Vector* v)
Creates a copy of the given vector.

### add_vectors(Vector* v1, Vector* v2)
Adds two vectors element-wise.
result[i] = v1[i] + v2[i]

### subtract_vectors(Vector* v1, Vector* v2)
Subtracts two vectors element-wise.
result[i] = v1[i] - v2[i]

### angle_between_vectors(Vector* v1, Vector* v2)
Calculates the angle between two vectors in radians.
θ = arccos((v1 · v2) / (|v1| * |v2|))

### distance_between_vectors(Vector* v1, Vector* v2)
Calculates the Euclidean distance between two vectors.
d = √((v1[0] - v2[0])² + (v1[1] - v2[1])² + ... + (v1[n] - v2[n])²)

### dot_product(Vector* v1, Vector* v2)
Calculates the dot product of two vectors.
v1 · v2 = v1[0]*v2[0] + v1[1]*v2[1] + ... + v1[n]*v2[n]

### cross_product(Vector* v1, Vector* v2)
Calculates the cross product of two 3D vectors.
v1 × v2 = (v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0])

### transform_vector(Matrix* m, Vector* v)
Transforms a vector by a matrix.
result = M * v

### scale_vector(float c, Vector* v)
Scales a vector by a scalar.
result[i] = c * v[i]

### magnitude_vector(Vector* v)
Calculates the magnitude (length) of a vector.
|v| = √(v[0]² + v[1]² + ... + v[n]²)

### normalize_vector(Vector* v)
Normalizes a vector to unit length.
v_norm = v / |v|

### reflection_vector(Vector* i, Vector* n)
Calculates the reflection of a vector i about a normal vector n.
r = i - 2(i · n)n

### max_vector(Vector* v1, Vector* v2)
Returns a vector with each component being the maximum of the corresponding components of v1 and v2.

### min_vector(Vector* v1, Vector* v2)
Returns a vector with each component being the minimum of the corresponding components of v1 and v2.

### projection_vector(Vector* v1, Vector* v2)
Calculates the projection of v1 onto v2.
proj_v2(v1) = ((v1 · v2) / (v2 · v2)) * v2

## Matrix Functions

### create_matrix(int dimension)
Creates a square matrix with the specified dimension.

### free_matrix(Matrix* m)
Frees the memory allocated for the matrix.

### add_matrices(Matrix* m1, Matrix* m2)
Adds two matrices element-wise.
result[i][j] = m1[i][j] + m2[i][j]

### subtract_matrices(Matrix* m1, Matrix* m2)
Subtracts two matrices element-wise.
result[i][j] = m1[i][j] - m2[i][j]

### multiply_matrices(Matrix* m1, Matrix* m2)
Multiplies two matrices.
result[i][j] = Σ(m1[i][k] * m2[k][j]) for k = 0 to n-1

### transpose_matrix(Matrix* m)
Transposes a matrix.
result[i][j] = m[j][i]

### scale_matrix(float c, Matrix* m)
Scales a matrix by a scalar.
result[i][j] = c * m[i][j]

### scale_matrix_anisotropic(float x, float y, float z, Matrix* m)
Scales a 4x4 matrix anisotropically along x, y, and z axes.

### determinant_matrix(Matrix* m)
Calculates the determinant of a matrix (up to 4x4).

### translate_matrix(float x, float y, float z)
Creates a 4x4 translation matrix.

### translate_matrix_in_place(Matrix* m, float x, float y, float z)
Applies a translation to an existing 4x4 matrix.

### outer_product(Vector* v1, Vector* v2)
Calculates the outer product of two vectors.
result[i][j] = v1[i] * v2[j]

### duplicate_matrix(Matrix* m)
Creates a copy of the given matrix.

### identity_matrix(int dimension)
Creates an identity matrix of the specified dimension.

### rotate_matrix(Matrix* m, float x, float y, float z, float angle)
Rotates a matrix around an arbitrary axis (x, y, z) by the given angle.

### rotate_matrix_X(Matrix* m, float angle)
Rotates a matrix around the X-axis by the given angle.

### rotate_matrix_Y(Matrix* m, float angle)
Rotates a matrix around the Y-axis by the given angle.

### rotate_matrix_Z(Matrix* m, float angle)
Rotates a matrix around the Z-axis by the given angle.

### invert_matrix(Matrix* m)
Calculates the inverse of a 4x4 matrix.

### look_at_matrix(Vector* eye, Vector* center, Vector* up)
Creates a view matrix for a camera.

### orthographic_projection_matrix(float left, float right, float bottom, float top, float near, float far)
Creates an orthographic projection matrix.

### perspective_projection_matrix(float fov, float aspect, float near, float far)
Creates a perspective projection matrix.

## Quaternion Functions

### create_quaternion(float x, float y, float z, float w)
Creates a quaternion with the given components.

### free_quaternion(Quaternion* q)
Frees the memory allocated for the quaternion.

### add_quaternions(Quaternion* q1, Quaternion* q2)
Adds two quaternions.
result = (q1.x + q2.x, q1.y + q2.y, q1.z + q2.z, q1.w + q2.w)

### subtract_quaternions(Quaternion* q1, Quaternion* q2)
Subtracts two quaternions.
result = (q1.x - q2.x, q1.y - q2.y, q1.z - q2.z, q1.w - q2.w)

### dot_product_quaternion(Quaternion* q1, Quaternion* q2)
Calculates the dot product of two quaternions.
q1 · q2 = q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w

### scale_quaternion(float c, Quaternion* q)
Scales a quaternion by a scalar.
result = (c*q.x, c*q.y, c*q.z, c*q.w)

### normalize_quaternion(Quaternion* q)
Normalizes a quaternion to unit length.
q_norm = q / |q|

### identity_quaternion()
Creates an identity quaternion (0, 0, 0, 1).

### quaternion_conjugate(Quaternion* q)
Calculates the conjugate of a quaternion.
q* = (-q.x, -q.y, -q.z, q.w)

### multiply_quaternions(Quaternion* q1, Quaternion* q2)
Multiplies two quaternions.
result = 
(q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w,
q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z)

### rotate_quaternion(Quaternion* q, Vector* axis, float angle)
Rotates a quaternion around an axis by the given angle.

### matrix_from_quaternion(Quaternion* q)
Converts a quaternion to a 4x4 rotation matrix.