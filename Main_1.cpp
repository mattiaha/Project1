#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <time.h>
using namespace std;
using namespace arma;
ofstream ofile;
// function to find the right-hand solution;
inline double f(double x){return 100.0 * exp(-10.0 * x);
}
// function to find the analytic result;
inline double exact(double x) { return 1.0 - (1 - exp(-10)) * x - exp(-10 * x);
}

void spec(int, string);
void gen(int , string);
void LU_decomp(int, string); 

int main(int argc, char** argv)
{
	
	char* filename;
	// file name must be declared in command window ;
	if (argc <= 1) {
		cout << "Bad Usage:" << argv[0] <<
			"read also file name" << endl;
		exit(1);
	}
	else {
		filename = argv[1];
		int funcname; int exponent;
		// choose which function to run;
		cout << "Write in choice of function: 1 (special), 2 (general) or 3 (LU decomposition)\n";
		cin >> funcname;
		// choose maximum value of n;
		printf("Read in from screen max n-value for power of 10^n \n");
		cin >> exponent;
		if (funcname == 1) {
			spec(exponent, filename);
		}
		else if (funcname == 2) {
			gen(exponent, filename);
			
		}
		else if (funcname == 3) {
			LU_decomp(exponent, filename);
		}
		else {
			cout << "Bad Usage\n" <<
				"write in choice of function: 1 (special), 2 (general)  or 3 (LU_decomp)\n";
			cout << filename << endl;
				
			exit(2);

		}
	}

	return 0;
}

void spec(int exponent, string filename) {
	for (int i = 1; i <= exponent; i++) {
		int n = (int)pow(10.0, i);
		//add choice of function and max n to name of file;
		string fileout = to_string(i);
		string argument1 = "special_";
		string argument2 = filename;
		
		fileout.append(argument1);
		fileout.append(argument2);
		// declare step-size;
		double h = 1.0 / (n);
		double hh = h * h;
		//we declare the vectors to be used by using the Armadillo library for C++;
		vec d= vec(n + 1); vec b = vec(n + 1);	vec v = vec(n + 1); vec x = vec(n + 1);
		// we only have n-1 elements in RelativeError, as the max value will be 0 if the other elements are negative;
		vec RelativeError = vec(n - 1);
		d[0] = d[n] = 2.0; v[0] = v[n] = 0.0;
		
		clock_t start, finish;
		// The computer writes down the time at the start of the algorithm;
		start = clock();
		// fill vector d;
		for (int i = 1; i < n; i++) d[i] = (i + 1.0) / ((double)i);
		//fill vectors x and b;
		for (int i = 0; i <= n; i++) {
			x[i] = i * h;
			b[i] = hh * f(i*h);
			
		}
		
		// Forward substitution;
		for (int i = 2; i < n; i++) { b[i] = b[i] + b[i - 1] / d[i - 1]; }
		// Backward substitution;
		v[n - 1] = b[n-1] / d[n - 1];
		for (int i = n - 2; i > 0; i--) v[i] = (b[i] + v[i + 1]) / d[i];
		// The computer writes down the time when the algorithm has finished;
		finish = clock();
		ofile.open(fileout);
		// we write out the values to our selected file;
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile <<"    x:          approx:          exact:            relative error(log10)" << endl;
		for (int i = 1; i <= n-1; i++) {

			double xval = x[i];
			
			RelativeError[i-1] = log10(fabs((exact(xval) - v[i]) / exact(xval)));
			ofile << setw(15) << setprecision(8) << xval;
			ofile << setw(15) << setprecision(8) << v[i];
			ofile << setw(15) << setprecision(8) << exact(xval);
			ofile << setw(15) << setprecision(8) << RelativeError[i-1] << endl;
			
		}
		
		// write out current n, max error, mean error and time to the command window we run our code in;
		cout << "n:";
		cout << log10(n) << endl;
		cout<< "Max error (log10): ";
		cout << setw(15) << setprecision(8) << max(RelativeError) << endl;
		cout << "Mean error (log10)";
		cout << setw(15) << setprecision(8) << mean(RelativeError) << endl;
		cout << "Time spent:";
		cout << setw(15) << setprecision(8) << ((finish - start) / double(CLOCKS_PER_SEC)) << endl;
		ofile.close();
		
		
	}

}
void gen(int exponent, string filename) {
	// we only want to look at matrices up to the size of 10^4 x 10^4;
	if (exponent > 4) {
		cout << "n too high, maximum n allowed is 4" << endl;
		exit(3);
	}
	else {
		for (int j = 1; j <= exponent; j++) {
			int n = (int)pow(10.0, j);
			string fileout = to_string(j);
			string argument1 = "general_";
			string argument2 = filename;
			// Add power of 10 and choice of function to the file name;
			fileout.append(argument1);
			fileout.append(argument2);
			// We use the Armadillo library to generate an n x n matrix;
			mat A = zeros<mat>(n, n);
			// For a simple code we fill in the first and last row before the loop;
			A(0, 0) = 2.0; A(0, 1) = -1.0;
			A(n - 1, n - 1) = 2.0; A(n - 1, n - 2) = -1.0;
			
			for (int i = 1; i < n - 1; i++) {
				// Here we set the diagonal elements of the matrix to the desired value;
				A(i, i) = 2.0;
				A(i, i + 1) = -1.0;
				A(i, i - 1) = -1.0;
			}
			// declare vectors
			vec a = vec(n); vec b = vec(n); vec c = vec(n);
			vec u = vec(n); vec x = vec(n); vec ba = vec(n); vec s = vec(n); vec r = vec(n);
			// we have n-2 elements in RelativeError to prevent the first and last elements from throwing off mean and max values;
			vec RelativeError = vec(n-2);
			// declare step-size;
			double h = 1.0 / (n);
			double hh = h * h;
			clock_t start, finish;
			// start the algorithm for a general tridiagonal matrix;
			start = clock();
			
			// fill vectors x and ba(right-hand vector);
			for (int i = 0; i < n; i++) {
				x[i] = i * h;
				ba[i] = hh * f(i * h);

			}
			
			// we declare the end points before the loop;
			a[0] = 0; b[0] = A(0, 0); c[0] = A(0, 1);
			a[n - 1] = A(n - 1, n - 2); b[n - 1] = A(n - 1, n - 1); c[n - 1] = 0;
			// fill three vectors with the tridiagonal elements of A;
			for (int i = 1; i < n - 1; i++) {
				b[i] = A(i, i);
				a[i] = A(i, i - 1);
				c[i] = A(i, i + 1);
			}
			// declare start values of s and r;
			s[0] = b[0]; r[0] = hh * f(0);
			// forward substitution;
			for (int i = 1; i < n; i++) {
				s[i] = b[i] - (a[i] * c[i - 1]) / s[i - 1];
				r[i] = ba[i] - (a[i] * r[i - 1]) / s[i - 1];
			}
			
			// backward substitution;
			u[n - 1] = r[n - 1] / b[n - 1];
			for (int i = n - 1; i > 0; i--) {
				u[i - 1] = (r[i - 1] - c[i - 1] * u[i]) / s[i-1];
			}
			// The computer writes down the time when the algorithm has finished;
			finish = clock();

			ofile.open(fileout);
			// write data to file;
			ofile << setiosflags(ios::showpoint | ios::uppercase);
			ofile << "    x:          approx:          exact:            relative error(log10)" << endl;
			// ignoring the first element and last element of the arrays;
			for (int i = 1; i < n -1; i++) {
				double xval = x[i];
				// calculate the error between exact and analytic solution;
				RelativeError[i-1] = log10(fabs((exact(xval) - u[i]) / exact(xval)));
				ofile << setw(15) << setprecision(8) << xval;
				ofile << setw(15) << setprecision(8) << u[i];
				ofile << setw(15) << setprecision(8) << exact(xval);
				ofile << setw(15) << setprecision(8) << RelativeError[i-1] << endl;
			}

			// show current power of 10;
			cout << "n:";
			cout << log10(n) << endl;
			// write to terminal the biggest error for current power of 10;
			cout << "Max error (log10): ";
			cout << setw(15) << setprecision(8) << max(RelativeError) << endl;
			cout << "Mean error (log10)";
			cout << setw(15) << setprecision(8) << mean(RelativeError) << endl;
			// write to terminal time of algorithm for current power of 10;
			cout << "Time spent:";
			cout << setw(15) << setprecision(8) << ((finish - start) / double(CLOCKS_PER_SEC)) << endl;
			ofile.close();
			}

	}
}


void LU_decomp(int exponent, string filename) {
	// we only want to look at matrices up to the size of 10^4 x 10^4;
	if (exponent > 4) {
		cout << "n too high, maximum n allowed is 4" << endl;
		exit(3);
	}
	else {
		for (int j = 1; j <= exponent; j++) {
			int n = (int)pow(10.0, j);
			string fileout = to_string(j);
			string argument1 = "L_U";
			string argument2 = filename;
			// Add power of 10 and choice of function to the file name;
			fileout.append(argument1);
			fileout.append(argument2);
			// We use the Armadillo library to generate an n x n matrix;
			mat A = zeros<mat>(n, n);
			// For a simple code we fill in the first and last row before the loop;
			A(0, 0) = 2.0; A(0, 1) = -1.0;
			A(n - 1, n - 1) = 2.0; A(n - 1, n - 2) = -1.0;

			for (int i = 1; i < n - 1; i++) {
				// Here we set the diagonal elements of the matrix to the desired value;
				A(i, i) = 2.0;
				A(i, i + 1) = -1.0;
				A(i, i - 1) = -1.0;
			}
			// declare vectors;
			vec x = vec(n); vec b = vec(n); vec y = vec(n); vec sol = vec(n);
			vec RelativeError = vec(n-2);
			// declare step-size;
			double h = 1.0 / (n);
			double hh = h * h;
			for (int i = 0; i < n; i++) {
				x[i] = i * h;
				b[i] = hh * f(i * h);

			}
			//declare matrices L and U;
			mat L, U;
			clock_t start, finish;
			// start the algorithm for a general tridiagonal matrix;
			start = clock();
			// use armadillos function for LU-decomposition;
			lu(L, U, A);
			// use armadillos function to solve linear equations of the form Ax = b;
			y = solve(L, b);
			sol = solve(U, y);
			// the algorithm for finding the solution ends here;
			finish = clock();
			ofile.open(fileout);
			// write data to file;
			ofile << setiosflags(ios::showpoint | ios::uppercase);
			ofile << "    x:          approx:          exact:            relative error(log10)" << endl;
			// ignoring the first element of the arrays as this throws the relative error off;
			for (int i = 1; i <= n - 1; i++) {
				double xval = x[i];
				// calculate the error between exact and analytic solution;
				RelativeError[i-1] = log10(fabs((exact(xval) - sol[i]) / exact(xval)));
				ofile << setw(15) << setprecision(8) << xval;
				ofile << setw(15) << setprecision(8) << sol[i];
				ofile << setw(15) << setprecision(8) << exact(xval);
				ofile << setw(15) << setprecision(8) << RelativeError[i-1] << endl;
			}
			// The computer writes down the time when the algorithm has finished;

			// show current power of 10;
			cout << "n:";
			cout << log10(n) << endl;
			// write to terminal the biggest error for current power of 10;
			cout << "Max error (log10): ";
			cout << setw(15) << setprecision(8) << max(RelativeError) << endl;
			// write to terminal the mean error (log10);
			cout << "Mean error (log10)";
			cout << setw(15) << setprecision(8) << mean(RelativeError) << endl;
			// write to terminal time of algorithm for current power of 10;
			cout << "Time spent:";
			cout << setw(15) << setprecision(8) << ((finish - start) / double(CLOCKS_PER_SEC)) << endl;
			ofile.close();
		}
	}
}