// Gaussian elimination with finding the largest element in a submatrix.

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

#define EPS 0.1E-18

using namespace std;

int p = 0;
int index_in_answer[1000000];

// A function that reads data from a file and returns a coefficient matrix and a vector of free terms.
pair<vector<vector<long double>>, vector<long double>> open_file (const string& PATH) {

    ifstream my_file(PATH);

    if (my_file.is_open()) {

        int n;

        my_file >> n;

        vector<vector<long double>> matrix(n, vector<long double>(n));
        vector<long double> b(n);

        for (int i = 0; i < n; i++) {

            for (int j = 0; j < n; j++) {

                if (j == n - 1) {

                    my_file >> matrix[i][j];
                    my_file >> b[i];
                } else my_file >> matrix[i][j];
            }
        }

        return make_pair(matrix, b);
    } else {

        cout << "Error!";

        vector<long double> empty_free_terms;
        vector<vector<long double >> empty_matrix;

        return make_pair(empty_matrix, empty_free_terms);
    }
}

// A function that performs the Gauss method forward.
pair<vector<vector<long double>>, vector<long double>> gauss (vector<vector<long double>> matrix, vector<long double> b) {

    int n = matrix.size();

    // We remember the order of the roots because in this method
    // we will change the columns and the order of the roots from that changes.
    for(int i  = 0; i < n; i++){
        index_in_answer[i] = i;
    }

    for(int i = 0; i < n - 1; i++){

        long double max = abs(matrix[i][i]);
        int i_index = i;
        int j_index = i;
        bool flag = false;

        // Search for the largest element in the submatrix.
        for(int ii = i; ii < n; ii++){

            for(int j = i; j < n; j++){

                if(max < abs(matrix[ii][j])){

                    max = abs(matrix[ii][j]);
                    i_index = ii;
                    j_index = j;
                    flag = true;
                }
            }
        }

        if(flag){

            // If we found max element in another row and column we must swap row and column and determinant will change twice.
            if(i_index != i && j_index != i)
                p += 2;
            // Else we swap row or column and determinant will change once.
            else
                p++;
        }

        // Swap row.
        for(int j = 0; j < n; j++){

            swap(matrix[i_index][j], matrix[i][j]);
        }

        // Swap column.
        for(int j = 0; j < n; j++){

            swap(matrix[j][j_index], matrix[j][i]);
        }

        // Swap order of the root if we swapped column.
        swap(index_in_answer[i], index_in_answer[j_index]);

        swap(b[i_index], b[i]);

        // We look for coefficients to zero the elements under the main diagonal and do it.
        for (int j = i + 1; j < n; j++) {

            long double tik = -1 * (matrix[j][i] / matrix[i][i]);
            // Change a free terms.
            b[j] = b[j] + b[i] * tik;

            for (int k = i; k < n; k++) {

                // Change every element in j row.
                matrix[j][k] = matrix[j][k] + matrix[i][k] * tik;
            }
        }
    }

    // Return a pair of triangle matrix and free terms.
    return make_pair(matrix, b);
}

// A function which, having a diagonal matrix and a vector of free terms,
// looks for the solution of the equation system and the determinant of the matrix
// and also writes it all to a file.
void write_answer(vector<vector<long double>> matrix, vector<long double> b, const string& PATH) {

    int n = matrix.size();
    long double determinant = 1.0;

    ofstream my_file(PATH);

    // Look for values close to zero on the main diagonal.
    for (int i = 0; i < n; i++) {

        if (abs(matrix[i][i]) <= EPS) {

            my_file << "The matrix has no unambiguous solution";
            my_file << "\nThe determinant of matrix is: 0";
            return;
        }
    }

    // If main diagonal has not value close to zero, calculating solution
    vector<long double> result(n);

    result[index_in_answer[n - 1]] = b[n - 1] / matrix[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--) {

        long double adder_to_free_elements = 0.0;

        for (int k = i + 1; k < n; k++) {
            adder_to_free_elements += matrix[i][k] * result[index_in_answer[k]];
        }

        result[index_in_answer[i]] = (b[i] - adder_to_free_elements) / matrix[i][i];
    }

    // Write in file.
    for (int i = 0; i < n; i++){

        my_file << fixed << setprecision(6) << "X" << i+1 << " = " << result[i] << endl;
        determinant *= matrix[i][i];
    }

    my_file << "The determinant of matrix is: " << determinant*pow(-1, p);
}

int main() {

    const string PATH = R"(D:\Projects\Numerical Methods\AdvancedGauss\Matrix4.txt)", PATH_ANSWER = R"(D:\Projects\Numerical Methods\AdvancedGauss\Answer4.txt)";

    pair<vector<vector<long double>>, vector<long double>> matrix_and_b = open_file(PATH);

    pair<vector<vector<long double>>, vector<long double >> triangle_matrix = gauss(matrix_and_b.first,
                                                                                    matrix_and_b.second);

    write_answer(triangle_matrix.first, triangle_matrix.second, PATH_ANSWER);

    return 0;
}
