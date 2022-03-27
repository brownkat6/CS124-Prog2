#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <bits/stdc++.h>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;

// TODO: run ./strassen 4 3 input_file.txt
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

// TODO: calculate number of triangles in a matrix

/**
 * Optimization roadmap
 * 1) Reduce the number of n/2 by n/2 allocations needed for the strassen() algorithm. Done.
 * 2) Avoid padding with extra zeros needlessly
 * 3) Modify wrapper code to take the input/output formats specified in the assignment. Done.
 */

int cross_over_point = 128;

void padMatrix(vector< vector<long long int> > &mat, int size, int new_size) {
    if (size % 2 == 1) {
        for (int i = 0; i < size; ++i) {
            vector<long long int> v_prime(new_size-size,0);
            mat[i].reserve(mat[i].size() + distance(v_prime.begin(),v_prime.end()));
            mat[i].insert(mat[i].end(),v_prime.begin(),v_prime.end());
        }
    }
    for (int i = 0; i < new_size-size; ++i) {
        mat.push_back(vector<long long int> (new_size,0));
    }
}

int getNextPowerOf2(int n) {
    return pow(2, int(ceil(log2(n))));
}

void populate_matrix_values(vector< vector<long long int> > &arr, int r, int c) {
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            arr[i][j] = rand()%2;
        }
    }
}

void populate_matrix_values_p(vector< vector<long long int> > &arr, int s, float p) {
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            arr[i][j] = (rand()%100 < p*100) ? 1 : 0;
        }
    }
}

void populate_matrix_zeros(vector< vector<long long int> > &arr, int s) {
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            arr[i][j] = 0;
        }
    }
}

void output_values_along_diagonal(vector< vector<long long int> > &matrix, int s) {
    //cout << "-------Output Matrix Values Along Diagonal--------" << endl;
    for (int i = 0; i < s; ++i) {
        cout << matrix[i][i] << endl;
    }
    //cout << "-------Finish Outputting Matrix Values Along Diagonal--------" << endl;
}

void print_matrix(vector< vector<long long int> > &matrix, int r, int c){
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
void naive_matrix_multiplication(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int n, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Multiply arr1 and arr2 and store result
    // Note: We use i k j ordering to speedup multiplication

    /**
    for(int i = 0; i < n; ++i) 
        for (int k = 0; k < n; ++k) 
            for(int j = 0; j < n; ++j) 
                res[i+r3][j+c3] += arr1[i+r1][k+c1] * arr2[k+r2][j+c2];
                //res[i+r3][j+c3] += ((i+r1<size1 && k+c1<size1) ? arr1[i+r1][k+c1] : 0) * ((k+r2<size2 && j+c2<size2) ? arr2[k+r2][j+c2] : 0);
            //*/
    ///**
    int size1 = arr1.size();
    int size2 = arr2.size();
    for(int i = 0; i < min(size1-r1,n); ++i) 
        for (int k = 0; k < min(size1-c1,min(size2-r2,n)); ++k) 
            for(int j = 0; j < min(size2-c2,n); ++j) 
                res[i+r3][j+c3] += arr1[i+r1][k+c1] * arr2[k+r2][j+c2];//*/
}

void add_matrices(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int s, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Add arr1 and arr2 and store result
    /**
    for(int i = 0; i < s; ++i)
        for(int j = 0; j < s; ++j)
            res[i+r3][j+c3] = arr1[i+r1][j+c1] + arr2[i+r2][j+c2];*/
    ///**
    int size1 = arr1.size();
    int size2 = arr2.size();
    for(int i = 0; i < s; ++i) {
        for(int j = 0; j < s; ++j) {
            //res[i+r3][j+c3] = arr1[i+r1][j+c1] + arr2[i+r2][j+c2];
            res[i+r3][j+c3] = ((i+r1<size1 && j+c1<size1) ? arr1[i+r1][j+c1] : 0) + ((i+r2<size2 && j+c2<size2) ? arr2[i+r2][j+c2] : 0);
        }
    }
    //*/
}

void subtract_matrices(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int s, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Add arr1 and arr2 and store result
    int size1 = arr1.size();
    int size2 = arr2.size();
    for(int i = 0; i < s; ++i) {
        for(int j = 0; j < s; ++j) {
            //res[i+r3][j+c3] = arr1[i+r1][j+c1] - arr2[i+r2][j+c2];
            res[i+r3][j+c3] = ((i+r1<size1 && j+c1<size1) ? arr1[i+r1][j+c1] : 0) - ((i+r2<size2 && j+c2<size2) ? arr2[i+r2][j+c2] : 0);
        }
    }
    /**
    for(int i = 0; i < s; ++i)
        for(int j = 0; j < s; ++j)
            res[i+r3][j+c3] = arr1[i+r1][j+c1] - arr2[i+r2][j+c2];*/
}

// r is the number of rows in arr1, l is the number of cols in arr1 and rows in arr2, c is the number of cols in arr2
// Assumes input matrics are square and have dimensions that are powers of 2
void strassen(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int size, int r1, int c1, int r2, int c2, int r3, int c3, int cross_over_point) {
    //base case
    if (size < cross_over_point)
    {
        naive_matrix_multiplication(arr1,arr2,res,size,r1,c1,r2,c2,r3,c3);
        return;
    }
    // new size of matrices split in half
    int ns = size/2;
    
    vector< vector<long long int> > p(ns, vector<long long int> (ns, 0));
    vector< vector<long long int> > temp1(ns, vector<long long int> (ns, 0));
    vector< vector<long long int> > temp2(ns, vector<long long int> (ns, 0));
    
    // Calculate 7 Subproblems
    // Calculate P1
    subtract_matrices(arr2,arr2,temp1,ns,0+r2,ns+c2,ns+r2,ns+c2,0,0);
    strassen(arr1,temp1,p,ns,r1,c1,0,0,0,0,cross_over_point);
    // Use P1 for result matrix
    add_matrices(p,res,res,ns,0,0,r3,c3+ns,r3,c3+ns);
    add_matrices(p,res,res,ns,0,0,r3+ns,c3+ns,r3+ns,c3+ns);
    populate_matrix_zeros(p,ns);

    add_matrices(arr1,arr1,temp1,ns,0+r1,0+c1,0+r1,ns+c1,0,0); // add(a,b)
    strassen(temp1,arr2,p,ns,0,0,r2+ns,c2+ns,0,0,cross_over_point); // P2 = (A+B)*H - temp1 h
    // Use P2 for result matrix
    add_matrices(p,res,res,ns,0,0,r3,c3+ns,r3,c3+ns);
    subtract_matrices(res,p,res,ns,r3,c3,0,0,r3,c3);
    populate_matrix_zeros(p,ns);

    add_matrices(arr1,arr1,temp1,ns,ns+r1,0+c1,ns+r1,ns+c1,0,0); // add(c,d)
    strassen(temp1,arr2,p,ns,0,0,r2,c2,0,0,cross_over_point); // P3 = (c+d)*e - temp1 e
    // Use P3 for result matrix
    add_matrices(p,res,res,ns,0,0,r3+ns,c3,r3+ns,c3);
    subtract_matrices(res,p,res,ns,r3+ns,c3+ns,0,0,r3+ns,c3+ns);
    populate_matrix_zeros(p,ns);

    subtract_matrices(arr2,arr2,temp1,ns,ns+r2,0+c2,0+r2,0+c2,0,0);
    strassen(arr1,temp1,p,ns,r1+ns,c1+ns,0,0,0,0,cross_over_point); // P4 - d temp1
    // Use P4 for result matrix
    add_matrices(p,res,res,ns,0,0,r3,c3,r3,c3);
    add_matrices(p,res,res,ns,0,0,r3+ns,c3,r3+ns,c3);
    populate_matrix_zeros(p,ns);

    add_matrices(arr1,arr1,temp1,ns,0+r1,0+c1,ns+r1,ns+c1,0,0);
    add_matrices(arr2,arr2,temp2,ns,0+r2,0+c2,ns+r2,ns+c2,0,0);
    strassen(temp1,temp2,p,ns,0,0,0,0,0,0,cross_over_point); // P5 = a*d+e*h
    // Use P5 for result matrix
    add_matrices(p,res,res,ns,0,0,r3,c3,r3,c3);
    add_matrices(p,res,res,ns,0,0,r3+ns,c3+ns,r3+ns,c3+ns);
    populate_matrix_zeros(p,ns);

    subtract_matrices(arr1,arr1,temp1,ns,0+r1,ns+c1,ns+r1,ns+c1,0,0);
    add_matrices(arr2,arr2,temp2,ns,ns+r2,0+c2,ns+r2,ns+c2,0,0);
    strassen(temp1,temp2,p,ns,0,0,0,0,0,0,cross_over_point); // P6
    // Use P6 for result matrix
    add_matrices(p,res,res,ns,0,0,r3,c3,r3,c3);
    populate_matrix_zeros(p,ns);

    subtract_matrices(arr1,arr1,temp1,ns,ns+r1,0+c1,0+r1,0+c1,0,0);
    add_matrices(arr2,arr2,temp2,ns,0+r2,0+c2,0+r2,ns+c2,0,0);
    strassen(temp1,temp2,p,ns,0,0,0,0,0,0,cross_over_point); // P7
    // Use P7 for result matrix
    add_matrices(p,res,res,ns,0,0,r3+ns,c3+ns,r3+ns,c3+ns);

    /**
    • AE +BG = P4 +P5 +P6 - P2
    • AF +BH = P1 +P2
    • CE +DG = P3 +P4
    • CF +DH = P1 −P3 +P5 +P7
    */
}

// Run strassen on input matrices of any size
void runStrassen(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int n, int cross_over_point)
{  
// Check to see if these matrices have dimensions of a power of 2. If not,
// the matrices must be resized and padded with zeroes to meet this criteria.
/**
 * TODO: is there a way to handle indexing so that we only need to add 1 row/column of zeros for odd-numbered
 *      matrices for each call to strassen()?
 * To avoid padding by zeros, we could pass two indices tracking the last non-bogus row and the last non-bogus column
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
// TODO: add function that pads an extra row/column to a matrix if the matrix has odd rows/columns

    int s = getNextPowerOf2(n);
    
    if (false) {
    //if (s!=n) {
        padMatrix(arr1,n,s); // Resize each matrix to the nearest power of 2
        padMatrix(arr2,n,s);
        padMatrix(res,n,s);
    }
    strassen(arr1,arr2,res,s,0,0,0,0,0,0,cross_over_point);
    return;

    /**
    //vector< vector<long long int> > arr1Resized(s, vector<long long int> (s, 0));
    //vector< vector<long long int> > arr2Resized(s, vector<long long int> (s, 0));
    //vector< vector<long long int> > resResized(s, vector<long long int> (s, 0));

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
    }*/
}

// ---------------------------Extra Verification/Utility Code Below-------------------------------

bool verify_strassen_equals_naive_multiplication(vector< vector<long long int> > &naive_result, vector< vector<long long int> > &strassen_result, int n) {
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
    int matrix_sizes_to_test = 8;
    int n_values[12] = {16,32,64,128,256,512,1024,2048};
    
    for (int i = 0; i < matrix_sizes_to_test; ++i) {
        int n = n_values[i];
        vector< vector<long long int> > arr1(n, vector<long long int> (n,0));
        vector< vector<long long int> > arr2(n, vector<long long int> (n,0));
        populate_matrix_values(arr1, n, n);
        populate_matrix_values(arr2, n, n);
        vector< vector<long long int> > naive_result(n, vector<long long int> (n, 0));
        vector< vector<long long int> > strassen_result(n, vector<long long int> (n, 0));
        
        chrono::steady_clock::time_point beginN = chrono::steady_clock::now();
        naive_matrix_multiplication(arr1,arr2,naive_result,n,0,0,0,0,0,0);
        chrono::steady_clock::time_point endN = chrono::steady_clock::now();
        cout << "Naive multiplication -  " << n << " by " << n << " matrix: " << chrono::duration_cast<chrono::milliseconds>(endN - beginN).count() << "[ms]" << endl;
        
        chrono::steady_clock::time_point beginS = chrono::steady_clock::now();
        runStrassen(arr1,arr2,strassen_result, n, cross_over_point);
        chrono::steady_clock::time_point endS = chrono::steady_clock::now();
        cout << "Strassen -  " << n << " by " << n << " matrix: " << chrono::duration_cast<chrono::milliseconds>(endS - beginS).count() << "[ms]" << endl;
        // Assertion fails if result of strassen multiplication doesn't equal result of naive multiplication
        assert(verify_strassen_equals_naive_multiplication(strassen_result,naive_result,n));
    }
}

void measure_crossover_point() {
    int crossover_values_to_test = 23;
    //int n_values[12] = {16,32,64,128,256,512,1024,2048};
    int n_values[23] = {16,32,36,40,44,48,52,56,60,64,68,72,96,112,128,192,256,384,512,768,1024,1536,2048};

    int n = 2048;
    vector< vector<long long int> > arr1(n, vector<long long int> (n,0));
    vector< vector<long long int> > arr2(n, vector<long long int> (n,0));
    vector< vector<long long int> > strassen_result(n, vector<long long int> (n, 0));
    populate_matrix_values(arr1, n, n);
    populate_matrix_values(arr2, n, n);

    ofstream outdata;
    outdata.open("crossover_point_times2048.dat", ios_base::app); // opens the file
    if( !outdata ) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    for (int i = 0; i < crossover_values_to_test; ++i) {
        outdata << n_values[i] << endl;
    }
    outdata << endl;
    int numtrials = 10;
    
    for (int i = 0; i < crossover_values_to_test; ++i) {
        int crossover_value = n_values[i];
        chrono::steady_clock::time_point beginS = chrono::steady_clock::now();
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        runStrassen(arr1,arr2,strassen_result, n, crossover_value);
        chrono::steady_clock::time_point endS = chrono::steady_clock::now();
        // Output time collected over 5 trials
        int time = chrono::duration_cast<chrono::milliseconds>(endS - beginS).count()/numtrials;
        outdata << time << endl;
        cout << "Strassen crossover point: " << crossover_value << ": " << time << "[ms]" << endl;
    }
    outdata << endl;

    outdata.close();
}

void read_in_matrix_values(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2,  int n, std::string input_file) {
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

// -------------------------------------Task 3-------------------------------------------
//populate_matrix_values_p
// TODO: debug
float get_num_triangles(vector< vector<long long int> > &matrix, int s) {
    float num_triangles = 0;
    for (int i = 0; i < s; ++i) {
        num_triangles += matrix[i][i];
    }
    
    return num_triangles/6;
}

void calc_triangles() {
    int n = 1024;
    float p_values[5] = {.01,.02,.03,.04,.05};
    float expected_triangles[5] = {178.4,1427.5,4817.7,11419.7,22304.1};
    vector< vector<long long int> > arr(n, vector<long long int> (n, 0));
    vector< vector<long long int> > cubed(n, vector<long long int> (n, 0));
    for (int i = 0; i < 5; ++i) {
        float p = p_values[i];
        populate_matrix_values_p(arr,n,p);
        //print_matrix(arr,5,5);
        naive_matrix_multiplication(arr,arr,cubed,n,0,0,0,0,0,0);
        naive_matrix_multiplication(arr,cubed,cubed,n,0,0,0,0,0,0);
        print_matrix(cubed,5,5);
        float expected_num_triangles = expected_triangles[i]; //(1024 choose 3) * p^3
        float num_triangles = get_num_triangles(cubed,n);
        cout << p << ": actual num triangles = " << num_triangles << ", expected num triangles = " << expected_num_triangles << endl;
    }
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

    vector< vector<long long int> > arr1(n, vector<long long int> (n, 0));
    vector< vector<long long int> > arr2(n, vector<long long int> (n, 0));
    vector< vector<long long int> > naive_result(n, vector<long long int> (n, 0));
    vector< vector<long long int> > strassen_result(n, vector<long long int> (n, 0));

    if (flag == 2) {
        populate_matrix_values(arr1, n, n); // Random 0's and 1's
        populate_matrix_values(arr2, n, n);
    } else {
        read_in_matrix_values(arr1,arr2,n,input_file);
    }

    
    
    // Run matrix multiplication 
    naive_matrix_multiplication(arr1,arr2,naive_result, n,0,0,0,0,0,0);
    runStrassen(arr1,arr2,strassen_result,n,cross_over_point);
    
    if (flag == 0) {
        output_values_along_diagonal(strassen_result,n);
    } else if (flag == 1) {
        // Compare runtimes of strassen and ours
        measure_multiplication_time();
    } else if (flag == 4) {
        measure_crossover_point();
    } else if (flag == 5) {
        calc_triangles();
    } else {
        print_matrix(arr1,n,n);
        print_matrix(arr2,n,n);
        print_matrix(strassen_result,n,n);
        print_matrix(naive_result,n,n);
        assert(verify_strassen_equals_naive_multiplication(strassen_result,naive_result,n));
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
