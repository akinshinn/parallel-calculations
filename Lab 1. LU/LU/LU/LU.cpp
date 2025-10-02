#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>


using namespace std;

const int n = 4;
const int b = 2;
const int max_block_size = b;

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


void fill_random(
    vector<double>& A
) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = rand();

        }
    }
}


void zero(std::vector<double>& A)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i * n + j] = 0.0;
}


void DisplayVector(const vector<double>& A) {
    int n = A.size();
    cout << "{ ";
    for (int i = 0; i < n; ++i) {
        cout << A[i] << ", ";
    }
    cout << "} ";
    cout << endl;
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


// A: m x n
void DisplayFlattenMatrix(const vector<double>& A, int m, int n) {

    for (int i = 0; i < m; i++) {
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
    int n = sqrt(A.size());
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


void FindLowerTriangleMatrix(vector<double>& L, const vector<double>& U, const vector<double>& A, int local_n) {
    // U: d x d, L: n-d x d, A: n-d x d,
    double sum;
    DisplayFlattenMatrix(U);
    for (int i = 0; i < local_n - b; i++) {
        for (int j = 0; j < b; j++) {
            // L[i][j] = A[i][j] - sum( L[i][k] * U[k][j],k=0, j-1)
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i * b + k] * U[k * b + j];
            }
            cout << sum << endl;
            cout << A[i * b + j] << endl;
            cout << U[j * b + j] << endl;
            L[i * b + j] = (A[i * b + j] - sum) / U[j * b + j];
        }
    }
}


void FindUpperTriangleMatrix(const vector<double>& L, vector<double>& U, const vector<double>& A) {
    // L: b x b (lower triangular), U: b x (n-b), A: b x (n-b)
    // U[i][j] = (A[i][j] - sum(L[i][k] * U[k][j], k=0, i-1)) / L[i][i]
    int rows = b;
    int cols = U.size() / b;
    double sum = 0.0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[i * b + k] * U[k * cols + j];
            }
            U[i * cols + j] = (A[i * cols + j] - sum);
        }
    }
}


// A: m x l, B: l x n, res: m x n
void MultiplyAB(const vector<double>& A, const vector<double>& B, int m, int l, int n, vector<double>& res) {
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < l; ++k)
            for (int j = 0; j < n; ++j)
                res[i * n + j] += A[i * l + k] * B[k * n + j];
}


void display_minmax(const vector<double>& A) {
    int n = A.size();
    double minn = 1e16, maxx = -1e16;


    for (int i = 0; i < n; i++) {
        minn = min(minn, A[i]);
        maxx = max(maxx, A[i]);
    }
    cout << "maximum = " << maxx << endl;
    cout << "minimun = " << minn << endl;
}


void LU_bl(vector<double>& A, int n) {

    if (n <= max_block_size) {
        cout << "in if" << endl;
        DisplayFlattenMatrix(A);
        LU(A);
        return;
    }
    vector<double> A11(b * b), A12(b * (n - b)), A21((n - b) * b), A22((n - b) * (n - b));
    // A11: b * b, A12: b * (n-b), A21: (n-b) * b, A22: (n-b)*(n-b)

    // заполняем блоки
    for (int i = 0; i < b; i++) {
        for (int j = 0; j < b; j++) {
            A11[i * b + j] = A[i * n + j];
        }
        for (int j = 0; j < n - b; j++) {
            A12[i * (n - b) + j] = A[i * n + j + b];
        }
    }
    for (int i = 0; i < n - b; i++) {
        for (int j = 0; j < b; j++) {
            A21[i * b + j] = A[(i + b) * n + j];
        }
        for (int j = 0; j < n - b; j++) {
            A22[i * (n - b) + j] = A[(i + b) * n + j + b];
        }
    }

    // Вычисляем LU для A11
    LU(A11);
    // Получаем A12 = L11 * U12
    std::vector<double> U12(b * (n - b));


    FindUpperTriangleMatrix(A11, U12, A12);


    // Получаем A21 = L21 * U11

    cout << "A11" << endl;
    DisplayFlattenMatrix(A11);
    vector<double> L21(b * (n - b));
    FindLowerTriangleMatrix(L21, A11, A21, n);
    cout << "L21" << endl;
    DisplayFlattenMatrix(L21, n - b, b);


    cout << "A21" << endl;
    DisplayFlattenMatrix(A21, n - b, b);
    // Преобразовываем матрицу A22
    vector<double> A22_mul((n - b) * (n - b)), A22_diff((n - b) * (n - b));
    MultiplyAB(L21, U12, n - b, b, n - b, A22_mul);
    cout << "A22" << endl;
    DisplayFlattenMatrix(A22);
    cout << "mult" << endl;
    DisplayFlattenMatrix(A22_mul);
    FlattenMatrixDiff(A22, A22_mul, A22_diff);
    cout << "diff" << endl;
    DisplayFlattenMatrix(A22_diff);
    A22 = A22_diff;
    LU_bl(A22, n - b);
    //LU(A22);


    // все собрать в одну матрицу 



    for (int i = 0; i < b; i++) {

        //A11
        for (int j = 0; j < b; j++) {
            A[i * n + j] = A11[i * b + j];
        }
        //A12
        for (int j = 0; j < n - b; j++) {
            A[i * n + j + b] = U12[i * (n - b) + j];
        }
    }

    for (int i = 0; i < n - b; i++) {

        //A21
        for (int j = 0; j < b; j++) {
            A[(i + b) * n + j] = L21[i * b + j];
        }

        //A22
        for (int j = 0; j < n - b; j++) {
            A[(i + b) * n + j + b] = A22[i * (n - b) + j];
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
    cout << "product LU:" << endl;
    DisplayFlattenMatrix(productLU);
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
    vector<double> A(n * n), A_1(n * n), A_2(n * n);


    //cout << "Simple LU:" << endl;

    fill(A);
    DisplayFlattenMatrix(A);

    //fill_random(A);
    A_1 = A;

    LU(A_1);
    DisplayFlattenMatrix(A_1);
    A_2 = A;
    //fill(A_2);
    LU_bl(A_2, n);
    DisplayFlattenMatrix(A_2);
    //cout << "Block LU" << endl;


    //vector<double> A(n*n), A_2;
    ////fill_random(A);
    //fill(A);
    //// 
    ////A = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    //A_2 = A;
    ////DisplayFlattenMatrix(A);


    ////cout << "Simple LU" << endl;
    //auto start = std::chrono::steady_clock::now();
    //LU(A);
    //auto finish = std::chrono::steady_clock::now();
    //std::chrono::duration<double> elapsed = finish - start;
    //auto time = elapsed.count();
    ////cout << "time for Simple LU = " << time << endl;
    ////DisplayFlattenMatrix(A);

    ////DisplayFlattenMatrix(A);
    //start = std::chrono::steady_clock::now();
    //LU_bl(A_2, n);
    //finish = std::chrono::steady_clock::now();
    //elapsed = finish - start;
    //time = elapsed.count();
    //cout << "Block LU" << endl;
    ////cout << "time for Block LU = " << time << endl;
    //DisplayFlattenMatrix(A_2);

    //vector<double> diff(n * n);
    //FlattenMatrixDiff(A_1, A_2, diff);
    //display_minmax(diff);
    //DisplayFlattenMatrix(diff);

    cout << endl;
    cout << "Simple LU residual:" << endl;
    Residual(A, A_1);
    cout << endl;
    cout << "Block LU residual" << endl;
    Residual(A, A_2);


}
