#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>


using namespace std;    

const int n =10;

void fill(
    std::vector<double>& A,
    std::vector<double>& B)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = cos(i + j);
            B[i * n + j] = sin(i - j);
        }
}


void fill(
    std::vector<double>& A)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            A[i * n + j] = cos(i + j);
        }
}


void zero(std::vector<double>& A)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i * n + j] = 0.0;
}


void DisplayFlattenMatrix(const vector<double>& A) {
    int n = sqrt(A.size());

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[n * i + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


void FlattenMatrixDiff(const vector<double>& A, const vector<double>& B, vector<double>& res) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        res[i] = A[i] - B[i];
    }
}


void FlattenMatrixSum(const vector<double>& A, const vector<double>& B, vector<double>& res) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        res[i] = A[i] + B[i];
    }
}

void LU(vector<double>& A) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i];
        }

        for (int j = i + 1; j < n; j++) {
            for (int k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
}

void ProductLU(const vector<double>& A, vector<double>& res) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) {
                for (int k = 0; k < i; k++) {
                    //res[i][j] += A[i][k] * A[k][j];
                    res[i * n + j] += A[i * n + k] * A[k * n + j];
                }
                res[i * n + j] += A[i * n + j];
            }
            else {
                for (int k = 0; k <= j; k++) {
                    //res[i][j] += A[i][k] * A[k][j];
                    res[i * n + j] += A[i * n + k] * A[k * n + j];
                }
            }
        }
    }

}


void Residual(const vector<double>& A, vector<double>& LU) {
    vector<double> productLU(n * n);
    ProductLU(LU, productLU);
    double maxRes = -1;
    double minRes = 1e17;

    for (int i = 0; i < n * n; i++) {
        maxRes = max(maxRes, abs(productLU[i] - A[i]));
        minRes = min(minRes, abs(productLU[i] - A[i]));
    }
    cout << "maxRes = " << maxRes << ", minRes = " << minRes << endl;
}

int main()
{
    vector<double> A(n * n), copymatrixA(n*n), productLU(n * n), raznost(n*n);
    fill(A);
    fill(copymatrixA);

    DisplayFlattenMatrix(A);

    LU(A);

    DisplayFlattenMatrix(A);
    Residual(copymatrixA, A);
}

