#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <bits/stdc++.h>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;

// Command to Compile the c++ file: g++ -std=c++17 -O2 -Wall -Wextra strassen.cpp -o strassen -lm -lpthread ./strassen <args>
//          Deprecated: c++ -std=gnu++2a -Wall -g -O3 strassen.cc -o strassen
/**
1) Analytically determine crossover point
2) Generate matrices of a given size for testing matrix multiplication
3) Implement naive matrix multiplication
2) Implement Strassen's algorithm
3) Benchmark multiplication time for different sizes of matrix for Strassen's vs naive implementation
6) Implement Strassen's with cutoff to use naive multiplication below cutoff
*/

/**
 * Optimization roadmap
 * 1) Reduce the number of n/2 by n/2 allocations needed for the strassen() algorithm
 * 2) Avoid padding with extra zeros needlessly
 * 3) Modify wrapper code to take the input/output formats specified in the assignment
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

void output_values_along_diagonal(vector< vector<int> > &matrix, int s) {
    //cout << "-------Output Matrix Values Along Diagonal--------" << endl;
    for (int i = 0; i < s; ++i) {
        cout << matrix[i][i] << endl;
    }
    //cout << "-------Finish Outputting Matrix Values Along Diagonal--------" << endl;
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

/**

 */

//---------------------------Strassens----------------------------------
void naive_matrix_multiplication(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int n, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Multiply arr1 and arr2 and store result
    // Note: We use i k j ordering to speedup multiplication
    for(int i = 0; i < n; ++i)
        for(int k = 0; k < n; ++k)
            for(int j = 0; j < n; ++j)
            {
                res[i+r3][j+c3] += arr1[i+r1][k+c1] * arr2[k+r2][j+c2];
            }
}

void add_matrices(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int s, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Add arr1 and arr2 and store result
    for(int i = 0; i < s; ++i)
        for(int j = 0; j < s; ++j)
            res[i+r3][j+c3] = arr1[i+r1][j+c1] + arr2[i+r2][j+c2];
}

void subtract_matrices(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int s, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Add arr1 and arr2 and store result
    for(int i = 0; i < s; ++i)
        for(int j = 0; j < s; ++j)
            res[i+r3][j+c3] = arr1[i+r1][j+c1] - arr2[i+r2][j+c2];
}

// r is the number of rows in arr1, l is the number of cols in arr1 and rows in arr2, c is the number of cols in arr2
// Assumes input matrics are square and have dimensions that are powers of 2
void strassen(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int size, int r1, int c1, int r2, int c2, int r3, int c3) {
    //base case
    //cout << "start strassen " << size << " " << r1 << " " << c1 << " " << r2 << " " << c2 << " " << r3 << " " << c3 << endl;
    
    if (size == 1)
    {
        res[r3][c3] = arr1[r1][c1] * arr2[r2][c2];
        //res[0][0] = arr1[0][0] * arr2[0][0];
        //cout << "Return val 0" << endl;
        return;
    }
    //cout << "allocate submatrices" << endl;

    // Create 8 submatrices of size n/2
    int ns = size/2; // new size of matrices split in half
    // NOTE: there's a more efficient way to initialize these matrices
    // TODO: we are just reading the data in a-e, so we don't need to allocate/copy a-e submatrices data into 8 new matrices
    // TODO: there's a way to reallocate using fewer than 7 submatrices
    // We can calculate P1, add P1 to the submatrices that use it in the final result matrix, and then forget it, then
    //      start on the rest of the calculation
    // We don't need to allocate 4 submatrices ae_bg-cf_dh, we can just input values directly into the final matrix
    /**
     * @brief 
     * We can
     */
    //NOTE: 10^24 choose 3 number of triangles for part 3 of this assignment according to michael zhao
    // TODO: switch the strassen() function to take in r1,c1,r2,c2,r3,c3 input parameters and then
    //      remove all uses of the a-h submatrices in function calls
    /**
    vector< vector<int> > a(ns, vector<int> (ns, 0));
    vector< vector<int> > b(ns, vector<int> (ns, 0));
    vector< vector<int> > c(ns, vector<int> (ns, 0));
    vector< vector<int> > d(ns, vector<int> (ns, 0));
    vector< vector<int> > e(ns, vector<int> (ns, 0));
    vector< vector<int> > f(ns, vector<int> (ns, 0));
    vector< vector<int> > g(ns, vector<int> (ns, 0));
    vector< vector<int> > h(ns, vector<int> (ns, 0));*/
    vector< vector<int> > p1(ns, vector<int> (ns, 0));
    vector< vector<int> > p2(ns, vector<int> (ns, 0));
    vector< vector<int> > p3(ns, vector<int> (ns, 0));
    vector< vector<int> > p4(ns, vector<int> (ns, 0));
    vector< vector<int> > p5(ns, vector<int> (ns, 0));
    vector< vector<int> > p6(ns, vector<int> (ns, 0));
    vector< vector<int> > p7(ns, vector<int> (ns, 0));
    //vector< vector<int> > ae_bg(ns, vector<int> (ns, 0));
    //vector< vector<int> > af_bh(ns, vector<int> (ns, 0));
    //vector< vector<int> > ce_dg(ns, vector<int> (ns, 0));
    //vector< vector<int> > cf_dh(ns, vector<int> (ns, 0));
    vector< vector<int> > temp1(ns, vector<int> (ns, 0));
    vector< vector<int> > temp2(ns, vector<int> (ns, 0));
    // Divide the matrices into sub matrices of size/2 by size/2
    /**
    int i,j;
    for (i = 0; i < ns; i++)
            {
                for (j = 0; j < ns; j++)
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
    }*/
    //cout << "calc 7 subproblems" << endl;
    // Calculate 7 Subproblems
    //subtract_matrices(f,h,temp1,ns);
    //TODO: change subtract_matrices, everything else to not allocate new submatrices
    //subtract_matrices(&arr2[0][ns],&arr2[ns][ns]);
    subtract_matrices(arr2,arr2,temp1,ns,0+r2,ns+c2,ns+r2,ns+c2,0,0);
    // TODO: finish changing indexing for arr1, arr2 when calling add(), subtract(), multiply()
    //strassen(a,temp1,p1,ns); // P1
    strassen(arr1,temp1,p1,ns,r1,c1,0,0,0,0);
    add_matrices(arr1,arr1,temp1,ns,0+r1,0+c1,0+r1,ns+c1,0,0); // add(a,b)
    strassen(temp1,arr2,p2,ns,0,0,r2+ns,c2+ns,0,0); // P2 = (A+B)*H - temp1 h
    add_matrices(arr1,arr1,temp1,ns,ns+r1,0+c1,ns+r1,ns+c1,0,0); // add(c,d)
    strassen(temp1,arr2,p3,ns,0,0,r2,c2,0,0); // P3 = (c+d)*e - temp1 e
    //subtract_matrices(g,e,temp1,ns); // sub(g,e)
    subtract_matrices(arr2,arr2,temp1,ns,ns+r2,0+c2,0+r2,0+c2,0,0);
    strassen(arr1,temp1,p4,ns,r1+ns,c1+ns,0,0,0,0); // P4 - d temp1
    //add_matrices(a,d,temp1,ns);
    //add_matrices(e,h,temp2,ns);
    add_matrices(arr1,arr1,temp1,ns,0+r1,0+c1,ns+r1,ns+c1,0,0);
    add_matrices(arr2,arr2,temp2,ns,0+r2,0+c2,ns+r2,ns+c2,0,0);
    strassen(temp1,temp2,p5,ns,0,0,0,0,0,0); // P5 = a*d+e*h
    //subtract_matrices(b,d,temp1,ns);
    //add_matrices(g,h,temp2,ns);
    subtract_matrices(arr1,arr1,temp1,ns,0+r1,ns+c1,ns+r1,ns+c1,0,0);
    add_matrices(arr2,arr2,temp2,ns,ns+r2,0+c2,ns+r2,ns+c2,0,0);
    strassen(temp1,temp2,p6,ns,0,0,0,0,0,0); // P6
    //subtract_matrices(c,a,temp1,ns);
    //add_matrices(e,f,temp2,ns);
    subtract_matrices(arr1,arr1,temp1,ns,ns+r1,0+c1,0+r1,0+c1,0,0);
    add_matrices(arr2,arr2,temp2,ns,0+r2,0+c2,0+r2,ns+c2,0,0);
    strassen(temp1,temp2,p7,ns,0,0,0,0,0,0); // P7

    //cout << "finish calc 7 subproblems" << endl;

    /**
    • AE +BG = P4 +P5 +P6 - P2
    • AF +BH = P1 +P2
    • CE +DG = P3 +P4
    • CF +DH = P1 −P3 +P5 +P7
    */
    // Use 7 Subproblems to calculate the values for the four quadrants
    /**
    subtract_matrices(p6,p2,temp1,ns);
    add_matrices(temp1,p5,temp2,ns);
    add_matrices(temp2,p4,ae_bg,ns);
    add_matrices(p1,p2,af_bh,ns);
    add_matrices(p3,p4,ce_dg,ns);
    subtract_matrices(p1,p3,temp1,ns);
    add_matrices(temp1,p5,temp2,ns);
    add_matrices(temp2,p7,cf_dh,ns);*/
    subtract_matrices(p6,p2,temp1,ns,0,0,0,0,0,0);
    add_matrices(temp1,p5,temp2,ns,0,0,0,0,0,0);
    add_matrices(temp2,p4,res,ns,0,0,0,0,r3,c3);
    add_matrices(p1,p2,res,ns,0,0,0,0,r3,c3+ns);
    add_matrices(p3,p4,res,ns,0,0,0,0,r3+ns,c3+0);
    subtract_matrices(p1,p3,temp1,ns,0,0,0,0,0,0);
    add_matrices(temp1,p5,temp2,ns,0,0,0,0,0,0);
    add_matrices(temp2,p7,res,ns,0,0,0,0,r3+ns,c3+ns);

    // Combine the 4 quadrant subproblems into the results matrix
    /**
    for (i = 0; i < ns; i++)
    {
        for (j = 0; j < ns; j++)
        {
            res[i][j] = ae_bg[i][j];
            res[i][j + ns] = af_bh[i][j];
            res[i + ns][j] = ce_dg[i][j];
            res[i + ns][j + ns] = cf_dh[i][j];
        }
    }*/
}

// Run strassen on input matrices of any size
void runStrassen(vector< vector<int> > &arr1, vector< vector<int> > &arr2, vector< vector<int> > &res, int n)
{  
// Check to see if these matrices have dimensions of a power of 2. If not,
// the matrices must be resized and padded with zeroes to meet this criteria.
/**
 * @brief 
 * TODO: we don't need to pad by zeros
 * We could pass two indices tracking the last non-bogus row and the last non-bogus column
 * E.g. if we want to multiply two 17-17 matrices, then strassen() will split the problem into 4 16-16 matrices M1-M4,
 *      but call add, subtract, naive_matrix_multiplication functions with parameters
 *      add(M1, size=16, valid_row=16, valid_col=16)
 *      add(M2, 16, valid_row=16, valid_col=1)
 *      add(M3, 16, valid_row=1, valid_col=16)
 *      add(M4, 16, valid_row=1, valid_col=1)
 *      Then add, subtract, naive_matrix_multiplication() functions should just "generate" magic 
 *          0 values for the matrix values not within the valid_row, valid_col bounds
 *      This avoids allocating extra 0 values that we don't need --> speedup
 * 
 */
    int s = getNextPowerOf2(n);
    if (s==n) {
        // Don't need to resize matrices
        strassen(arr1,arr2,res, s,0,0,0,0,0,0);
        return;
    }
    vector< vector<int> > arr1Resized(s, vector<int> (s, 0));
    vector< vector<int> > arr2Resized(s, vector<int> (s, 0));
    vector< vector<int> > resResized(s, vector<int> (s, 0));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            arr1Resized[i][j] = arr1[i][j];
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            arr2Resized[i][j] = arr2[i][j];
        }
    }
    strassen(arr1Resized, arr2Resized, resResized, s,0,0,0,0,0,0);
    // NOTE: is there a way to avoid copying data over to the original matrix?
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = resResized[i][j];
        }
    }
}

// ---------------------------Extra Verification/Utility Code Below-------------------------------

bool verify_strassen_equals_naive_multiplication(vector< vector<int> > &naive_result, vector< vector<int> > &strassen_result, int n) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
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
    int n_values[12] = {16,32,64,128,256,512};
    
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
        naive_matrix_multiplication(arr1,arr2,naive_result,n,0,0,0,0,0,0);
        chrono::steady_clock::time_point endN = chrono::steady_clock::now();
        cout << "Naive multiplication -  " << n << " by " << n << " matrix: " << chrono::duration_cast<chrono::milliseconds>(endN - beginN).count() << "[ms]" << endl;
        
        chrono::steady_clock::time_point beginS = chrono::steady_clock::now();
        runStrassen(arr1,arr2,strassen_result, n);
        chrono::steady_clock::time_point endS = chrono::steady_clock::now();
        cout << "Strassen -  " << n << " by " << n << " matrix: " << chrono::duration_cast<chrono::milliseconds>(endS - beginS).count() << "[ms]" << endl;
        // Assertion fails if result of strassen multiplication doesn't equal result of naive multiplication
        assert(verify_strassen_equals_naive_multiplication(strassen_result,naive_result,n));
    }
    


}

void read_in_matrix_values(vector< vector<int> > &arr1, vector< vector<int> > &arr2,  int n, std::string input_file) {
    //cout << input_file << endl;
    std::ifstream pFile(input_file.c_str());
    string line;
    int x;
    int r = 0;
    int c = 0;
    while (r < n && pFile >> x) {
        arr1[r][c]=x;
        c++;
        if (c==n) {
            c=0;
            r++;
        }
    }
    r = 0;
    c = 0;
    while (r < n && pFile >> x) {
        arr2[r][c]=x;
        c++;
        if (c==n) {
            c=0;
            r++;
        }
    }
    pFile.close();
}

// ./strassen 0 dimension input_file
int main(int argc, char **argv) {
    assert(argc==4);
    int dimension = atoi(argv[2]);
    uint flag=atoi(argv[1]);
    std::string input_file = argv[3];

    // Get matrices
    (void) flag;
    int n = dimension;

    vector< vector<int> > arr1(n, vector<int> (n, 0));
    vector< vector<int> > arr2(n, vector<int> (n, 0));
    vector< vector<int> > naive_result(n, vector<int> (n, 0));
    vector< vector<int> > strassen_result(n, vector<int> (n, 0));

    read_in_matrix_values(arr1,arr2,n,input_file);

    //populate_matrix_values(arr1, n, n); // Random 0's and 1's
    //populate_matrix_values(arr2, n, n);
    
    // Run matrix multiplication 
    naive_matrix_multiplication(arr1,arr2,naive_result, n,0,0,0,0,0,0);
    runStrassen(arr1,arr2,strassen_result,n);
    
    
    // Display result of calculations
    
    //print_matrix(naive_result,n,n);
    //print_matrix(strassen_result,n,n);
    //output_values_along_diagonal(naive_result,n);
    if (flag == 0) {
        output_values_along_diagonal(strassen_result,n);
    } else if (flag == 1) {
        // Compare runtimes of strassen and ours
        measure_multiplication_time();
    } else {
        print_matrix(arr1,n,n);
        print_matrix(arr2,n,n);
        print_matrix(strassen_result,n,n);
        print_matrix(naive_result,n,n);
    }    

    
    return 0;

    
}

/**
Input/Output Description
Your code should take three arguments: a flag, a dimension, and an input file:
$ ./strassen 0 dimension inputfile
The flag 0 is meant to provide you some flexibility; you may use other values for your own testing,
debugging, or extensions. The dimension, which we refer to henceforth as d, is the dimension of the matrix
you are multiplying, so that 32 means you are multiplying two 32 by 32 matrices together. 
The inputfile is
an ASCII file with 2d2 integer numbers, one per line, representing two matrices A and B; you are to find
the product AB = C. The first integer number is matrix entry a0,0, followed by a0,1, a0,2, . . . , a0,d−1; next
comes a1,0, a1,1, and so on, for the first d2 numbers. The next d2 numbers are similar for matrix B.
Your program should put on standard output (in C: printf, cout, System.out, etc.) a list of the values
of the diagonal entries c0,0, n,1, . . . , cd−1,d−1, one per line, including a trailing newline. The output will
be checked by a script – add no clutter. (You should not output the whole matrix, although of course all
entries should be computed.)
 */
