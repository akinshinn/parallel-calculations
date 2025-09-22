#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>


using namespace std;    

const int n = 5;

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

int main()
{
    vector<double> A(n * n);
    fill(A);

    DisplayFlattenMatrix(A);

    LU(A);

    DisplayFlattenMatrix(A);
}

