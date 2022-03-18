#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <bits/stdc++.h>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;

// Command to Compile the c++ file: c++ -std=gnu++2a -Wall -g -O3 strassen.cc -o strassen
/**
1) Analytically determine crossover point
2) Generate matrices of a given size for testing matrix multiplication
3) Implement naive matrix multiplication
2) Implement Strassen's algorithm
3) Benchmark multiplication time for different sizes of matrix for Strassen's vs naive implementation
6) Implement Strassen's with cutoff to use naive multiplication below cutoff
*/

int getNextPowerOf2(int n) {
    return pow(2, int(ceil(log2(n))));
}

void populate_matrix_values(vector< vector<int> > &arr, int r, int c) {
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            arr[i][j] = rand()%2;
        }
    }
}

void print_matrix(vector< vector<int> > &matrix, int r, int c){
    cout << "-------Matrix--------" << endl;
    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            if (j != 0){
                cout << "\t";
            }
            cout << matrix[i][j];
        }
        cout << endl;
    }
}

void naive_matrix_multiplication(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int r1, int r2, int c1, int c2) {
    //vector< vector<int> > &res = get_zero_matrix(r1,c2);
    
    // Multiply arr1 and arr2 and store result
    for(int i = 0; i < r1; ++i)
        for(int j = 0; j < c2; ++j)
            for(int k = 0; k < c1; ++k)
            {
                res[i][j] += arr1[i][k] * arr2[k][j];
            }
    //return res;
}

void add_matrices(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int s) {
    // Add arr1 and arr2 and store result
    for(int i = 0; i < s; ++i)
        for(int j = 0; j < s; ++j)
            res[i][j] = arr1[i][j] + arr2[i][j];
}

void subtract_matrices(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int s) {
    // Add arr1 and arr2 and store result
    for(int i = 0; i < s; ++i)
        for(int j = 0; j < s; ++j)
            res[i][j] = arr1[i][j] - arr2[i][j];
}

//vector< vector<int> > &get_zero_matrix(int r, int c) {
//    vector< vector<int> > m(r, vector<int> (c, 0));
//    return &m;
    /**vector< vector<int> > &arr = new int*[r];
    for (int i = 0; i < r; ++i) {
        arr[i] = new int[c];
        for (int j = 0; j < c; ++j) {
            arr[i][j] = 0;
        }
    }
    return arr;*/
//}

/**

// Return an nxn matrix
vector< vector<int> > &get_matrix(int r, int c) {
    cout << "get_matrix()" << endl;
    vector< vector<int> > &arr = new int*[r];
    cout << "initialize matrix" << endl;
    for (int i = 0; i < r; ++i) {
        arr[i] = new int[c];
        for (int j = 0; j < c; ++j) {
            arr[i][j] = rand()%2;
        }
    }
    //cout << "finish initialize matrix" << endl;
    return arr;
}*/

// r is the number of rows in arr1, l is the number of cols in arr1 and rows in arr2, c is the number of cols in arr2
// Assumes input matrics are square and have dimensions that are powers of 2
void strassen(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int size) {
    // Initialize empty result matrix
    //vector< vector<int> > &res = get_zero_matrix(size,size);

    //base case
    if (size == 1)
    {
        res[0][0] = arr1[0][0] * arr2[0][0];
        return;// res;
    }

    // Create 8 submatrices of size n/2
    int ns = size/2; // new size of matrices split in half
    // NOTE: there's a more efficient way to initialize these matrices
    vector< vector<int> > a(ns, vector<int> (ns, 0));
    vector< vector<int> > b(ns, vector<int> (ns, 0));
    vector< vector<int> > c(ns, vector<int> (ns, 0));
    vector< vector<int> > d(ns, vector<int> (ns, 0));
    vector< vector<int> > e(ns, vector<int> (ns, 0));
    vector< vector<int> > f(ns, vector<int> (ns, 0));
    vector< vector<int> > g(ns, vector<int> (ns, 0));
    vector< vector<int> > h(ns, vector<int> (ns, 0));
    vector< vector<int> > p1(ns, vector<int> (ns, 0));
    vector< vector<int> > p2(ns, vector<int> (ns, 0));
    vector< vector<int> > p3(ns, vector<int> (ns, 0));
    vector< vector<int> > p4(ns, vector<int> (ns, 0));
    vector< vector<int> > p5(ns, vector<int> (ns, 0));
    vector< vector<int> > p6(ns, vector<int> (ns, 0));
    vector< vector<int> > p7(ns, vector<int> (ns, 0));
    vector< vector<int> > ae_bg(ns, vector<int> (ns, 0));
    vector< vector<int> > af_bh(ns, vector<int> (ns, 0));
    vector< vector<int> > ce_dg(ns, vector<int> (ns, 0));
    vector< vector<int> > cf_dh(ns, vector<int> (ns, 0));
    vector< vector<int> > temp1(ns, vector<int> (ns, 0));
    vector< vector<int> > temp2(ns, vector<int> (ns, 0));
    // Divide the matrices into sub matrices of size/2 by size/2
    int i,j;
    for (int i = 0; i < ns; i++)
            {
                for (int j = 0; j < ns; j++)
                {
                    a[i][j] = arr1[i][j];
                    b[i][j] = arr1[i][j + ns];
                    c[i][j] = arr1[i + ns][j];
                    d[i][j] = arr1[i + ns][j + ns];

                    e[i][j] = arr2[i][j];
                    f[i][j] = arr2[i][j + ns];
                    g[i][j] = arr2[i + ns][j];
                    h[i][j] = arr2[i + ns][j + ns];
                }
    }
    // Calculate 7 Subproblems
    subtract_matrices(f,h,temp1,ns);
    strassen(a,temp1,p1,ns); // P1
    add_matrices(a,b,temp1,ns);
    strassen(temp1,h,p2,ns); // P2 = (A+B)*H
    add_matrices(c,d,temp1,ns);
    strassen(temp1,e,p3,ns); // P3 = (c+d)*e
    subtract_matrices(g,e,temp1,ns);
    strassen(d,temp1,p4,ns); // P4
    add_matrices(a,d,temp1,ns);
    add_matrices(e,h,temp2,ns);
    strassen(temp1,temp2,p5,ns); // P5 = a*d+e*h
    subtract_matrices(b,d,temp1,ns);
    add_matrices(g,h,temp2,ns);
    strassen(temp1,temp2,p6,ns); // P6
    subtract_matrices(c,a,temp1,ns);
    add_matrices(e,f,temp2,ns);
    strassen(temp1,temp2,p7,ns); // P7

    /**
    • AE +BG = P4 +P5 +P6 - P2
    • AF +BH = P1 +P2
    • CE +DG = P3 +P4
    • CF +DH = P1 −P3 +P5 +P7
    */
    // Use 7 Subproblems to calculate the values for the four quadrants
    subtract_matrices(p6,p2,temp1,ns);
    add_matrices(temp1,p5,temp2,ns);
    add_matrices(temp2,p4,ae_bg,ns);
    add_matrices(p1,p2,af_bh,ns);
    add_matrices(p3,p4,ce_dg,ns);
    subtract_matrices(p1,p3,temp1,ns);
    add_matrices(temp1,p5,temp2,ns);
    add_matrices(temp2,p7,cf_dh,ns);

    // Combine the 4 quadrant subproblems into the results matrix
    for (i = 0; i < ns; i++)
    {
        for (j = 0; j < ns; j++)
        {
            res[i][j] = ae_bg[i][j];
            res[i][j + ns] = af_bh[i][j];
            res[i + ns][j] = ce_dg[i][j];
            res[i + ns][j + ns] = cf_dh[i][j];
        }
    }
}

// Run strassen on input matrices of any size
void runStrassen(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int r1, int c1, int r2, int c2)
{  
// Check to see if these matrices are already square and have dimensions of a power of 2. If not,
// the matrices must be resized and padded with zeroes to meet this criteria.
    int k = max({r1,c1,r2,c2});
    int s = getNextPowerOf2(k);
    if (s==k && r1==k && r2==k && c1==k && c2==k) {
        // Don't need to resize matrices
        strassen(arr1,arr2,res, s);
        return;
    }
    vector< vector<int> > arr1Resized(s, vector<int> (s, 0));
    vector< vector<int> > arr2Resized(s, vector<int> (s, 0));
    vector< vector<int> > resResized(s, vector<int> (s, 0));

    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c1; j++)
        {
            arr1Resized[i][j] = arr1[i][j];
        }
    }
    for (int i = 0; i < r2; i++)
    {
        for (int j = 0; j < c2; j++)
        {
            arr2Resized[i][j] = arr2[i][j];
        }
    }
    strassen(arr1Resized, arr2Resized, resResized, s);
    //vector< vector<int> > &res = get_zero_matrix(r1,c2);
    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c2; j++)
        {
            res[i][j] = resResized[i][j];
        }
    }
    //print_matrix(res, r1, c2);
    //return res;
}

bool verify_strassen_equals_naive_multiplication(vector< vector<int> > &naive_result, vector< vector<int> > &strassen_result, int r1, int c1, int c2) {
    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c2; j++)
        {
            if (naive_result[i][j] != strassen_result[i][j]) {
                return false;
            }
        }
    }
    return true;
}

void measure_multiplication_time() {
    int matrix_sizes_to_test = 7;
    int n_values[12] = {16,32,64,128, 256, 512, 1024};
    
    for (int i = 0; i < matrix_sizes_to_test; ++i) {
        int n = n_values[i];
        chrono::steady_clock::time_point beginM = chrono::steady_clock::now();
        vector< vector<int> > arr1(n, vector<int> (n,0));
        vector< vector<int> > arr2(n, vector<int> (n,0));
        vector< vector<int> > naive_result(n, vector<int> (n, 0));
        vector< vector<int> > strassen_result(n, vector<int> (n, 0));
        chrono::steady_clock::time_point endM = chrono::steady_clock::now();
        cout << "Time to initialize four  " << n << " by " << n << " matrices: " << chrono::duration_cast<chrono::milliseconds>(endM - beginM).count() << "[ms]" << endl;

        chrono::steady_clock::time_point beginN = chrono::steady_clock::now();
        naive_matrix_multiplication(arr1,arr2,naive_result,n,n,n,n);
        chrono::steady_clock::time_point endN = chrono::steady_clock::now();
        cout << "Naive multiplication -  " << n << " by " << n << " matrix: " << chrono::duration_cast<chrono::milliseconds>(endN - beginN).count() << "[ms]" << endl;
        
        chrono::steady_clock::time_point beginS = chrono::steady_clock::now();
        runStrassen(arr1,arr2,strassen_result, n,n,n,n);
        chrono::steady_clock::time_point endS = chrono::steady_clock::now();
        cout << "Strassen -  " << n << " by " << n << " matrix: " << chrono::duration_cast<chrono::milliseconds>(endS - beginS).count() << "[ms]" << endl;
        // Assertion fails if result of strassen multiplication doesn't equal result of naive multiplication
        assert(verify_strassen_equals_naive_multiplication(strassen_result,naive_result,n,n,n));
    }
    


}

// ./strassen 0 numpoints numtrials dimension
int main(int argc, char **argv) {
    // Get matrices
    int r1=2; // Define dimensions of input matrices
    int c1=3;
    int r2=c1;
    int c2=1;

    vector< vector<int> > arr1(r1, vector<int> (c1, 0));
    vector< vector<int> > arr2(r2, vector<int> (c2, 0));
    vector< vector<int> > naive_result(r1, vector<int> (c2, 0));
    vector< vector<int> > strassen_result(r1, vector<int> (c2, 0));

    populate_matrix_values(arr1, r1, c1); // Random 0's and 1's
    populate_matrix_values(arr2, r2, c2);
    
    // Run matrix multiplication 
    naive_matrix_multiplication(arr1,arr2,naive_result, r1,c1,c1,c2);
    runStrassen(arr1,arr2,strassen_result,r1,c1,c1,c2);
    
    
    // Display result of calculations
    print_matrix(arr1,r1,c1);
    print_matrix(arr2,c1,c2);
    print_matrix(naive_result,r1,c2);
    print_matrix(strassen_result,r1,c2);

    // Compare runtimes of strassen and ours
    measure_multiplication_time();
    return 0;

    
}
