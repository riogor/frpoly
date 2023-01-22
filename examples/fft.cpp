#include "fft.hpp"

#include <assert.h>
#include <ctime>
#include <iostream>

using namespace std;

int main()
{
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);
	cerr.tie(0);

	int n = 1000;

	frpoly::vll a(n), b(n);

	srand(time(NULL));
	for (int i = 0; i < n; ++i)
		a[i] = rand() % 5000, b[i] = rand() % 5000;

	frpoly::vll res_stupid;
	frpoly::multiply_polynomials_stupid(a, b, res_stupid);

	frpoly::vll res_fft;
	frpoly::multiply_polynomials_fft(a, b, res_fft);

	for (int i = 0; i < res_stupid.size(); ++i)
		assert(res_fft[i] == res_stupid[i]);

	return 0;
}