// Gaussian elimination with finding the largest element in a column.

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

#define EPS 1E-18

using namespace std;

int p = 0;

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

        my_file.close();
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

    for (int i = 0; i < n - 1; i++) {

        // Find max element in the column.
        long double max = abs(matrix[i][i]);
        long double max_index = i;

        // Flag to see if we found a larger item than the diagonal item.
        bool flag = false;

        for (int j = i + 1; j < n; j++) {

            if (max < abs(matrix[j][i])) {

                max = abs(matrix[j][i]);
                max_index = j;
                flag = true;
            }
        }

        // If we found larger item than the diagonal item we will swap and we will change value of determinant on opposite.
        // And we have counter how much we did this.
        if(flag){
            p++;
        }

        // Swap row.
        for (int j = 0; j < n; j++) {

            swap(matrix[max_index][j], matrix[i][j]);
        }

        // Swap free terms
        swap(b[max_index], b[i]);

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
            my_file.close();

            return;
        }
    }

    // If main diagonal has not value close to zero, calculating solution
    vector<long double> result(n);

    result[n - 1] = b[n - 1] / matrix[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--) {

        long double adder_to_free_elements = 0.0;

        for (int k = i + 1; k < n; k++) {
            adder_to_free_elements += matrix[i][k] * result[k];
        }

        result[i] = (b[i] - adder_to_free_elements) / matrix[i][i];
    }

    // Write in file.
    for (int i = 0; i < n; i++){

        my_file << fixed << setprecision(6) << "X" << i+1 << " = " << result[i] << endl;
        determinant *= matrix[i][i];
    }

    my_file << "The determinant of matrix is: " << determinant*pow(-1, p);
    my_file.close();
}

int main() {

    const string PATH = R"(D:\Projects\Numerical Methods\Gauss\Matrix4.txt)", PATH_ANSWER = R"(D:\Projects\Numerical Methods\Gauss\Answer4.txt)";

    pair<vector<vector<long double>>, vector<long double>> matrix_and_b = open_file(PATH);

    pair<vector<vector<long double>>, vector<long double >> triangle_matrix = gauss(matrix_and_b.first,
                                                                                        matrix_and_b.second);

    write_answer(triangle_matrix.first, triangle_matrix.second, PATH_ANSWER);

    return 0;
}
