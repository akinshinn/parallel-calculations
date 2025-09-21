#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

const int n = 4096;
//const int b = 32;


using namespace std;



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

void zero(std::vector<double>& A)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i * n + j] = 0.0;
}


double mulIKJ(
    const std::vector<double>& A,
    const std::vector<double>& B,
    std::vector<double>& C,
    int n)
{
    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                C[i * n + j] += A[i * n + k] * B[k * n + j];

    /*if (display) DisplayFlattenMatrix(C);*/
    auto finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    double time = elapsed.count();
    return time;
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
    //vector<double> res(n);
    for (int i = 0; i < n; i++) {
        res[i] = A[i] - B[i];
    }
}


void FlattenMatrixSum(const vector<double>& A, const vector<double>& B, vector<double>& res) {
    int n = A.size();
    //vector<double> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = A[i] + B[i];
    }
}


//void Strassen(const vector<double>& A, const vector<double>& B, int n, vector<double>& res) {
//    //vector<double> res(n*n);
//    //cout << n << endl;
//    if (n == 1) {
//        res[0] = A[0] * B[0];
//        //return res;
//    }
//
//    int d = n / 2; // размер блока 
//    vector<double> A11(d*d), A12(d*d), A21(d*d), A22(d*d),
//        B11(d*d), B12(d*d), B21(d*d), B22(d*d);
//    for (int i = 0; i < d; i++) {
//        for (int j = 0; j < d; j++) {
//            //cout << "i = " << i << " j = " << j << endl;
//            A11[i*d + j] = A[0 * d + 2*i*d + j];
//            A12[i*d + j] = A[1 * d + 2 * i * d + j];
//            A21[i * d + j] = A[2*d * d + 2 * i * d + j];
//            A22[i * d + j] = A[(2 *d + 1)* d + 2 * i * d + j];
//
//            B11[i * d + j] = B[0 * d + 2 * i * d + j];
//            B12[i * d + j] = B[1 * d + 2 * i * d + j];
//            B21[i * d + j] = B[2 * d * d + 2 * i * d + j];
//            B22[i * d + j] = B[(2 * d + 1) * d + 2 * i * d + j];
//        }
//    }
//    vector<double> P1(d*d), P2(d * d), P3(d * d), P4(d * d), P5(d * d), P6(d * d), P7(d * d), C11(d * d), C12(d * d), C21(d * d), C22(d * d), S1(d*d), S2(d * d), S3(d * d), S4(d * d), S5(d * d), S6(d * d), S7(d * d), S8(d * d), S9(d * d), S10(d * d) ;
//    
//    FlattenMatrixDiff(A12, A22, S1);
//    FlattenMatrixSum(B21, B22, S2);
//    Strassen(S1, S2, d, P1);
//
//    FlattenMatrixSum(A11, A22, S3);
//    FlattenMatrixSum(B11, B22, P4);
//    Strassen(S3, S4, d, P2);
//
//    FlattenMatrixDiff(A11, A21, S5);
//    FlattenMatrixSum(B11, B12, S6);
//    Strassen(S5, S6, d, P3);
//
//    FlattenMatrixSum(A11, A12, S7);
//    Strassen(S7, B22, d, P4);
//
//    FlattenMatrixDiff(B12, B22, S8);
//    Strassen(A11, S8, d, P5);
//
//    FlattenMatrixDiff(B21, B11, S9);
//    Strassen(A22, S9, d, P6);
//
//    FlattenMatrixSum(A21, A22, S10);
//    Strassen(S10, B11, d, P7);
//
//    FlattenMatrixSum(P1, P2, C11);
//    FlattenMatrixDiff(C11, P4, C11);
//    FlattenMatrixSum(C11, P6, C11);
//
//    FlattenMatrixSum(P4, P5, C12);
//    FlattenMatrixSum(P6, P7, C21);
//
//    FlattenMatrixDiff(P2, P3, C22);
//    FlattenMatrixSum(C22, P5, C22);
//    FlattenMatrixDiff(C22, P7,C22);
//
//    for (int i = 0; i < n / 2; i++) {
//        for (int j = 0; j < n / 2; j++) {
//            res[i * n + j] = C11[i*n/2 + j];
//            res[i * n + j + n / 2] = C12[i * n / 2 + j];
//            res[(i + n / 2) * n + j] = C21[i * n / 2 + j];
//            res[(i + n / 2) * n + j + n / 2] = C22[i * n / 2 + j];
//        }
//    }
//    //return res;
//}

void StrassenThreshold(const vector<double>& A, const vector<double>& B, int n, vector<double>& res) {
    //vector<double> res(n * n);
    //cout << n << endl;
    if (n <= 128) {
        //res[0] = A[0] * B[0];
        mulIKJ(A, B, res, n);
        return;
    }
    else {
        int d = n / 2; // размер блока 
        vector<double> A11(d * d), A12(d * d), A21(d * d), A22(d * d),
            B11(d * d), B12(d * d), B21(d * d), B22(d * d);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                //cout << "i = " << i << " j = " << j << endl;
                A11[i * d + j] = A[0 * d + 2 * i * d + j];
                A12[i * d + j] = A[1 * d + 2 * i * d + j];
                A21[i * d + j] = A[2 * d * d + 2 * i * d + j];
                A22[i * d + j] = A[(2 * d + 1) * d + 2 * i * d + j];
                B11[i * d + j] = B[0 * d + 2 * i * d + j];
                B12[i * d + j] = B[1 * d + 2 * i * d + j];
                B21[i * d + j] = B[2 * d * d + 2 * i * d + j];
                B22[i * d + j] = B[(2 * d + 1) * d + 2 * i * d + j];
            }
        }


        vector<double> P1(d * d), P2(d * d), P3(d * d), P4(d * d), P5(d * d), P6(d * d), P7(d * d), C11(d * d), C12(d * d), C21(d * d), C22(d * d), S1(d * d), S2(d * d);

        FlattenMatrixDiff(A12, A22, S1);
        FlattenMatrixSum(B21, B22, S2);


        StrassenThreshold(S1, S2, d, P1);

        FlattenMatrixSum(A11, A22, S1);
        FlattenMatrixSum(B11, B22, S2);
        StrassenThreshold(S1, S2, d, P2);

        FlattenMatrixDiff(A11, A21, S1);
        FlattenMatrixSum(B11, B12, S2);
        StrassenThreshold(S1, S2, d, P3);

        FlattenMatrixSum(A11, A12, S1);
        StrassenThreshold(S1, B22, d, P4);

        FlattenMatrixDiff(B12, B22, S1);
        StrassenThreshold(A11, S1, d, P5);

        FlattenMatrixDiff(B21, B11, S1);
        StrassenThreshold(A22, S1, d, P6);

        FlattenMatrixSum(A21, A22, S1);
        StrassenThreshold(S1, B11, d, P7);

        FlattenMatrixSum(P1, P2, C11);
        FlattenMatrixDiff(C11, P4, C11);
        FlattenMatrixSum(C11, P6, C11);

        FlattenMatrixSum(P4, P5, C12);
        FlattenMatrixSum(P6, P7, C21);

        FlattenMatrixDiff(P2, P3, C22);
        FlattenMatrixSum(C22, P5, C22);
        FlattenMatrixDiff(C22, P7, C22);

        for (int i = 0; i < n / 2; i++) {
            for (int j = 0; j < n / 2; j++) {
                res[i * n + j] = C11[i * n / 2 + j];
                res[i * n + j + n / 2] = C12[i * n / 2 + j];
                res[(i + n / 2) * n + j] = C21[i * n / 2 + j];
                res[(i + n / 2) * n + j + n / 2] = C22[i * n / 2 + j];
            }
        }
        return;

    }

}





int main()
{
    std::vector<double> A(n * n), B(n * n), C(n * n);

    fill(A, B);
    double time;

    time = mulIKJ(A, B, C, n);

    cout << "time for IKJ = " << time;
    cout << endl;
    auto start = std::chrono::steady_clock::now();
    //vector<double> res(n*n);
    StrassenThreshold(A, B, n, C);
    auto finish = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed = finish - start;
     time = elapsed.count();


    cout << "time for Strassen = " << time;

}

