test: test.c cvectors.h
	gcc -std=c99 -Wall -Werror -Wextra test.c -o test -D VECTOR_TEST -D MATRIX_TEST -D QUAT

vector: test.c cvectors.h
	gcc -std=c99 -Wall -Werror -Wextra test.c -o test -D VECTOR_TEST

matrix: test.c cvectors.h
	gcc -std=c99 -Wall -Werror -Wextra test.c -o test -D MATRIX_TEST

quat: test.c cvectors.h
	gcc -std=c99 -Wall -Werror -Wextra test.c -o test -D QUAT
