#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int** init_2d_array(int num_rows, int num_cols);
void* Malloc(size_t size);


/* 
 * Function: matrix_multiply
 * ---------------------------- 
 *   Computes X * Y with an interative algorithm.
 * 
 *   X: 2D matrix to left-multiply
 *   Y: 2D matrix to right-multiply
 *   x_rows: the number of rows in X
 *   x_cols: the number of columns in X
 *   y_rows: the number of rows in Y
 *   y_cols: the number of columns in Y
 * 
 *   returns: the matrix X * Y as a 2D array
 */
int** matrix_multiply(int** X, int** Y, int x_rows, int x_cols, int y_rows, int y_cols) {
	// Check whether X and Y are compatible
	if (x_cols != y_rows) {
		fprintf(stderr, "Matrix sizes are incompatible for multiplication.\n");
		exit(EXIT_FAILURE);
	}

	int z_rows = x_rows;
	int z_cols = y_cols;

	// Allocate the product array Z
	int **Z = init_2d_array(z_rows, z_cols);

	// Assume X and Y are initialized and filled such that for X[i][j] or
	// Y[i][j], i refers to the row number and j refers to the column number
	// This also means that X and Y are 2D arrays
	for (int cur_z_row = 0; cur_z_row < z_rows; cur_z_row++) {
		for (int cur_z_col = 0; cur_z_col < z_cols; cur_z_col++) {
			int dot_product = 0;

			for (int dot_index = 0; dot_index < x_cols; dot_index++)
				// Assume that these are valid indices, that is, X contains
				// exactly x_rows and x_cols, and Y contains exactly y_rows
				// and y_cols
				dot_product += X[cur_z_row][dot_index] * Y[dot_index][cur_z_col];

			Z[cur_z_row][cur_z_col] = dot_product;
		}
	}

	return Z;
}

/* 
 * Function: init_2d_array
 * ---------------------------- 
 *   Allocates a 2D array of size num_rows x num_cols.
 * 
 *   num_rows: the number of rows
 *   num_cols: the number of columns
 * 
 *   returns: the allocated 2D array
 */
int** init_2d_array(int num_rows, int num_cols) {
	int **R = (int **) Malloc(num_rows * sizeof(int *));
	for (int i = 0; i < num_rows; i++)
		R[i] = (int *) Malloc(num_cols * sizeof(int));

	return R;
}

/* 
 * Function: Malloc
 * ---------------------------- 
 *   Calls malloc with error checking. Exits with error code on failure.
 * 
 *   size: the size of memory to allocate
 * 
 *   returns: the allocated memory
 */
void* Malloc(size_t size) {
	void *to_ret;
	if ((to_ret = malloc(size)) == NULL) {
		perror("Malloc");
		exit(EXIT_FAILURE);
	}
	return to_ret;
}


/* --------------------------------------------------------- */
/* Below are additional functions that were used for testing */
/* --------------------------------------------------------- */

void fill_rand_2d_array(int** R, int num_rows, int num_cols, int upper) {
	for (int i = 0; i < num_rows; i++)
		for (int j = 0; j < num_cols; j++)
			R[i][j] = rand() % upper;
}

void free_2d_array(int** R, int num_rows, int num_cols) {
	for (int i = 0; i < num_rows; i++)
		free(R[i]);
	free(R);
}

void print_matrix(int **R, int r_rows, int r_cols) {
	for (int i = 0; i < r_rows; i++) {
		for (int j = 0; j < r_cols; j++) {
			printf("%d ", R[i][j]);
		}
		printf("\n");
	}
}

int main() {
	int x_rows = 4, x_cols = 5;
	int y_rows = 5, y_cols = 3;

	srand(time(NULL)); 

	int **X = init_2d_array(x_rows, x_cols);
	fill_rand_2d_array(X, x_rows, x_cols, 10);

	int **Y = init_2d_array(y_rows, y_cols);
	fill_rand_2d_array(Y, y_rows, y_cols, 10);

	int **Z = matrix_multiply(X, Y, x_rows, x_cols, y_rows, y_cols);

	printf("---X---\n");
	print_matrix(X, x_rows, x_cols);
	printf("---Y---\n");
	print_matrix(Y, y_rows, y_cols);
	printf("---Z---\n");
	print_matrix(Z, x_rows, y_cols);

	free_2d_array(X, x_rows, x_cols);
	free_2d_array(Y, y_rows, y_cols);
	free_2d_array(Z, x_rows, y_cols);
}
