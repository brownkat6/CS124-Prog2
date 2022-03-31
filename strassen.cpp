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

int cross_over_point = 68;

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
        for (int j = i+1; j < s; ++j) {
            int v = (rand()%100 < p*100) ? 1 : 0;
            arr[i][j] = v;
            arr[j][i] = v;
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
    for (int i = 0; i < s; ++i) {
        cout << matrix[i][i] << endl;
    }
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

//---------------------------Strassens----------------------------------
void naive_matrix_multiplication(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int n, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Multiply arr1 and arr2 and store result
    // Note: We use i k j ordering to speedup multiplication
    int size1 = arr1.size();
    int size2 = arr2.size();
    for(int i = 0; i < min(size1-r1,n); ++i) 
        for (int k = 0; k < min(size1-c1,min(size2-r2,n)); ++k) 
            for(int j = 0; j < min(size2-c2,n); ++j) 
                res[i+r3][j+c3] += arr1[i+r1][k+c1] * arr2[k+r2][j+c2];
}

void naive_matrix_multiplication_no_padding(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int n) {
    // Multiply arr1 and arr2 and store result
    // Note: We use i k j ordering to speedup multiplication
    for(int i = 0; i < n; ++i) 
        for (int k = 0; k < n; ++k) 
            for(int j = 0; j < n; ++j) 
                res[i][j] += arr1[i][k] * arr2[k][j];
}

void add_matrices(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int s, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Add arr1 and arr2 and store result
    int size1 = arr1.size();
    int size2 = arr2.size();
    int sizeRes = res.size();
    for(int i = 0; i < min(sizeRes-r3,s); ++i) {
        for(int j = 0; j < min(sizeRes-c3,s); ++j) {
            res[i+r3][j+c3] = ((i+r1<size1 && j+c1<size1) ? arr1[i+r1][j+c1] : 0) + ((i+r2<size2 && j+c2<size2) ? arr2[i+r2][j+c2] : 0);
        }
    }
}

void subtract_matrices(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int s, int r1, int c1, int r2, int c2, int r3, int c3) {
    // Add arr1 and arr2 and store result
    int size1 = arr1.size();
    int size2 = arr2.size();
    int sizeRes = res.size();
    for(int i = 0; i < min(sizeRes-r3,s); ++i) {
        for(int j = 0; j < min(sizeRes-c3,s); ++j) {
            res[i+r3][j+c3] = ((i+r1<size1 && j+c1<size1) ? arr1[i+r1][j+c1] : 0) - ((i+r2<size2 && j+c2<size2) ? arr2[i+r2][j+c2] : 0);
        }
    }
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
}

// Run strassen on input matrices of any size
void runStrassen(vector< vector<long long int> > &arr1, vector< vector<long long int> > &arr2, vector< vector<long long int> > &res, int n, int cross_over_point)
{  
    int s = getNextPowerOf2(n);
    strassen(arr1,arr2,res,s,0,0,0,0,0,0,cross_over_point);
    return;
}

// ---------------------------Extra Verification/Utility Code Below-------------------------------

bool verify_matrices_are_equal(vector< vector<long long int> > &naive_result, vector< vector<long long int> > &strassen_result, int n) {
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (naive_result[i][j] != strassen_result[i][j]) {
                cout << "Assertion result: false" << endl;
                return false;
            }
        }
    }
    cout << "Assertion result: true" << endl;
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
        assert(verify_matrices_are_equal(strassen_result,naive_result,n));
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
float get_num_triangles(vector< vector<long long int> > &matrix, int s) {
    float num_triangles = 0;
    for (int i = 0; i < s; ++i) {
        num_triangles += static_cast< float >(matrix[i][i]);
    }
    
    return num_triangles/6;
}

void calc_triangles() {
    int n = 1024;
    float p_values[5] = {.01,.02,.03,.04,.05};
    float expected_triangles[5] = {178.4,1427.5,4817.7,11419.7,22304.1};
    vector< vector<long long int> > arr(n, vector<long long int> (n, 0));
    for (int i = 0; i < 5; ++i) {
        for (int trial = 0; trial < 5; ++trial) {
            float p = p_values[i];
            populate_matrix_values_p(arr,n,p);
            vector< vector<long long int> > squared(n, vector<long long int> (n, 0));
            vector< vector<long long int> > cubed(n, vector<long long int> (n, 0));
            runStrassen(arr,arr,squared,n,cross_over_point);
            runStrassen(arr,squared,cubed,n,cross_over_point);
            float expected_num_triangles = expected_triangles[i]; //(1024 choose 3) * p^3
            float num_triangles = get_num_triangles(cubed,n);
            cout << p << ": actual num triangles = " << num_triangles << ", expected num triangles = " << expected_num_triangles << endl;
        }
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
    vector< vector<long long int> > naive_result_no_pad(n, vector<long long int> (n, 0));
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
    naive_matrix_multiplication_no_padding(arr1,arr2,naive_result_no_pad,n);
    
    
    if (flag == 0) {
        output_values_along_diagonal(strassen_result,n);
        //assert(verify_matrices_are_equal(strassen_result,naive_result,n));
    } else if (flag == 1) {
        // Compare runtimes of strassen and ours
        measure_multiplication_time();
    } else if (flag == 4) {
        measure_crossover_point();
    } else if (flag == 5) {
        calc_triangles();
    } else {
        //output_values_along_diagonal(strassen_result,n);
        assert(verify_matrices_are_equal(strassen_result,naive_result,n));
        assert(verify_matrices_are_equal(naive_result_no_pad,naive_result,n));
        cout << "assertion passed" << endl;
    }    

    
    return 0;

    
}
