#include <stdio.h>
#include <stdlib.h>

struct CSR_Matrix *init_CSR_matrix(int num_val, int num_rows, int num_cols);
void* Malloc(size_t size);


/* 
 * A matrix in Compressed Row Storage format, which stores just the non-zero
 * values of a matrix in row-major order.
 * Saves on a lot of memory over a 2D array representation when the number of
 * non-zero values are less than (m * (n - 1) - 1) / 2 for an m x n matrix,
 * which would normally require m x n space.
 * Saves time on matrix computations as not all comparisons have to be made,
 * just the comparisons of non-zero values.
 */
struct CSR_Matrix {
	int *val;  /* The non-zero values in the matrix */
	int *col_ind;  /* The column indices of the corresponding values in val */

	/* Points to the non-zero values at the start of each row */
	/* Defined recursively as
		i) row_ptr[0] = 0
		ii) row_ptr[i] = row_ptr[i - 1] + the number of non-zero values in the
			ith row */
	int *row_ptr;
	int num_rows;
	int num_cols;
};

/* 
 * A matrix in Compressed Column Storage format. Similar to CSR except in
 * column-major order.
 */
struct CCS_Matrix {
	int *val;
	int *row_ind;
	int *col_ptr;
	int num_rows;
	int num_cols;
};


/* 
 * Function: sparse_matrix_multiply
 * ---------------------------- 
 *   Computes X * Y sparse matrices. The benefits of multiplying a CSR matrix
 *   with a CCS matrix are optimal; CSR matrices are fast at directly traversing
 *   rows, and CCS matrices are fast at directly traversing columns. Both of
 *   these corresponding tasks are required in the compuation X * Y.
 * 
 *   X: 2D sparse matrix to left-multiply
 *   Y: 2D sparse matrix to right-multiply
 * 
 *   returns: the matrix X * Y as a CSR matrix
 */
struct CSR_Matrix *sparse_matrix_multiply(struct CSR_Matrix *X, struct CCS_Matrix *Y) {
	// Check whether X and Y are compatible
	if (X->num_cols != Y->num_rows) {
		fprintf(stderr, "Matrix sizes are incompatible for multiplication.\n");
		exit(EXIT_FAILURE);
	}

	int z_rows = X->num_rows;
	int z_cols = Y->num_cols;

	/* Create these on the stack and fill them in to a struct later */
	/* Although more time-consuming, ensures the structs are space efficient */
	/* Space efficiency is important in compressed formats for future
	*  computations */
	int z_val[z_rows * z_cols];
	int z_col_ind[z_rows * z_cols];
	int z_row_ptr[z_rows + 1];
	int z_val_count = 0;

	z_row_ptr[0] = 0;

	for (int cur_z_row = 0; cur_z_row < z_rows; cur_z_row++) {
		int cur_row_val_count = 0;

		for (int cur_z_col = 0; cur_z_col < z_cols; cur_z_col++) {
			int dot_product = 0;

			/* For traversing the row in Y */
			int y_row = Y->col_ptr[cur_z_col];

			/* Traverse the column in X */
			for (int x_col = X->row_ptr[cur_z_row];
				x_col < X->row_ptr[cur_z_row + 1]; x_col++) {
				while (Y->row_ind[y_row] < X->col_ind[x_col] && \
					y_row < Y->col_ptr[cur_z_col + 1]) y_row++;

				/* There are no more non-zero values in the row in Y, no need to
				*  continue */
				if (y_row >= Y->col_ptr[cur_z_col + 1]) break;

				if (Y->row_ind[y_row] == X->col_ind[x_col])
					dot_product += X->val[x_col] * Y->val[y_row];
			}

			if (dot_product != 0) {
				z_val[z_val_count] = dot_product;
				z_col_ind[z_val_count] = cur_z_col;
				z_val_count++;
				cur_row_val_count++;
			}
		}

		z_row_ptr[cur_z_row + 1] = z_row_ptr[cur_z_row] + cur_row_val_count;
	}

	struct CSR_Matrix *Z = init_CSR_matrix(z_val_count, z_rows, z_cols);

	/* Fill in the CSR matrix */
	for (int i = 0; i < z_val_count; i++) {
		Z->val[i] = z_val[i];
		Z->col_ind[i] = z_col_ind[i];
	}

	for (int i = 0; i < z_rows + 1; i++) Z->row_ptr[i] = z_row_ptr[i];

	return Z;
}

/* 
 * Function: init_CSR_matrix
 * ---------------------------- 
 *   Allocates a CSR matrix of size num_rows x num_cols and num_val non-zero
 *   values.
 * 
 *   num_val: the number of non-zero values
 *   num_rows: the number of rows
 *   num_cols: the number of columns
 * 
 *   returns: the allocated CSR matrix
 */
struct CSR_Matrix *init_CSR_matrix(int num_val, int num_rows, int num_cols) {
	struct CSR_Matrix *R = (struct CSR_Matrix *) Malloc(sizeof(struct CSR_Matrix));
	R->val = (int *) Malloc(num_val * sizeof(int));
	R->col_ind = (int *) Malloc(num_val * sizeof(int));
	R->row_ptr = (int *) Malloc((num_rows + 1) * sizeof(int));
	R->num_rows = num_rows;
	R->num_cols = num_cols;

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

struct CCS_Matrix *init_CCS_matrix(int num_val, int num_rows, int num_cols) {
	struct CCS_Matrix *R = (struct CCS_Matrix *) Malloc(sizeof(struct CCS_Matrix));
	R->val = (int *) Malloc(num_val * sizeof(int));
	R->row_ind = (int *) Malloc(num_val * sizeof(int));
	R->col_ptr = (int *) Malloc((num_rows + 1) * sizeof(int));
	R->num_rows = num_rows;
	R->num_cols = num_cols;

	return R;
}

void free_CSR_matrix(struct CSR_Matrix *R) {
	free(R->val);
	free(R->col_ind);
	free(R->row_ptr);
	free(R);
}

void free_CCS_matrix(struct CCS_Matrix *R) {
	free(R->val);
	free(R->row_ind);
	free(R->col_ptr);
	free(R);
}

void print_CSR_matrix(struct CSR_Matrix *R) {
	for (int cur_row = 0; cur_row < R->num_rows; cur_row++) {
		int cur_ptr = R->row_ptr[cur_row];

		for (int cur_col = 0; cur_col < R->num_cols; cur_col++) {
			if (cur_col == R->col_ind[cur_ptr] && \
				cur_ptr < R->row_ptr[cur_row + 1]) {
				printf("%d ", R->val[cur_ptr]);
				cur_ptr++;
			} else printf("0 ");
		}

		printf("\n");
	}
}

/* Not very time efficient, but sufficient for the purposes of testing small
*  matrices */
void print_CCS_matrix(struct CCS_Matrix *R) {
	for (int cur_row = 0; cur_row < R->num_rows; cur_row++) {
		for (int cur_col = 0; cur_col < R->num_cols; cur_col++) {
			int val = 0;

			for (int cur_ptr = R->col_ptr[cur_col];
				cur_ptr < R->col_ptr[cur_col + 1] && \
				R->row_ind[cur_ptr] <= cur_row;
				cur_ptr++) {
				if (R->row_ind[cur_ptr] == cur_row) val = R->val[cur_ptr];
			}

			printf("%d ", val);
		}

		printf("\n");
	}
}

int main() {
	int x_rows = 7, x_cols = 5;
	int y_rows = 5, y_cols = 6;

	/* Manually create X and Y matrices for one test case */
	struct CSR_Matrix *X = init_CSR_matrix(6, x_rows, x_cols);
	X->val[0] = 2;
	X->val[1] = 4;
	X->val[2] = 3;
	X->val[3] = 1;
	X->val[4] = 6;
	X->val[5] = 2;
	X->col_ind[0] = 0;
	X->col_ind[1] = 3;
	X->col_ind[2] = 2;
	X->col_ind[3] = 0;
	X->col_ind[4] = 1;
	X->col_ind[5] = 4;
	X->row_ptr[0] = 0;
	X->row_ptr[1] = 2;
	X->row_ptr[2] = 2;
	X->row_ptr[3] = 3;
	X->row_ptr[4] = 4;
	X->row_ptr[5] = 4;
	X->row_ptr[6] = 5;
	X->row_ptr[7] = 6;

	struct CCS_Matrix *Y = init_CCS_matrix(9, y_rows, y_cols);
	Y->val[0] = 3;
	Y->val[1] = 11;
	Y->val[2] = 2;
	Y->val[3] = 3;
	Y->val[4] = 5;
	Y->val[5] = 4;
	Y->val[6] = 2;
	Y->val[7] = 6;
	Y->val[8] = 5;
	Y->row_ind[0] = 0;
	Y->row_ind[1] = 4;
	Y->row_ind[2] = 1;
	Y->row_ind[3] = 1;
	Y->row_ind[4] = 3;
	Y->row_ind[5] = 0;
	Y->row_ind[6] = 1;
	Y->row_ind[7] = 2;
	Y->row_ind[8] = 4;
	Y->col_ptr[0] = 0;
	Y->col_ptr[1] = 2;
	Y->col_ptr[2] = 3;
	Y->col_ptr[3] = 5;
	Y->col_ptr[4] = 6;
	Y->col_ptr[5] = 8;
	Y->col_ptr[6] = 9;

	struct CSR_Matrix *Z = sparse_matrix_multiply(X, Y);

	printf("---X---\n");
	print_CSR_matrix(X);
	printf("---Y---\n");
	print_CCS_matrix(Y);
	printf("---Z---\n");
	print_CSR_matrix(Z);

	free_CSR_matrix(X);
	free_CCS_matrix(Y);
	free_CSR_matrix(Z);
}
