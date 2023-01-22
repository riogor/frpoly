#ifndef FRPOLY__FFT_HPP
#define FRPOLY__FFT_HPP

#include "pch.hpp"

#include "helpers.hpp"

namespace frpoly
{
using base = std::complex<double>;
using vbase = std::vector<base>;

#ifndef FRPOLY__FFT_OPTIMIZED

void calculate_roots(unsigned int n, vbase &roots)
{
	roots.resize(n);
	for (int i = 0; i < n / 2; ++i)
	{
		double angle = (2.0 * M_PI * i) / n;
		roots[i + n / 2] = base(cos(angle), sin(angle));
	}

	for (int i = n / 2 - 1; i >= 1; --i)
		roots[i] = roots[i * 2];
}

#else

void calculate_roots(unsigned int n, vbase &roots)
{
	roots.resize(n);
	roots[1] = base(1, 0);

	for (int k = 1; k < (32 - std::countl_zero(n)) - 1; ++k)
	{
		double angle = 2 * M_PI / (1 << (k + 1));
		base w(cos(angle), sin(angle));

		for (int i = (1 << (k - 1)); i < (1 << k); ++i)
		{
			roots[i << 1] = roots[i];
			roots[(i << 1) + 1] = roots[i] * w;
		}
	}
}

#endif

vbase calculate_roots(unsigned int n)
{
	vbase roots(n);

	calculate_roots(n, roots);

	return roots;
}

#ifndef FRPOLY__FFT_ITERATIVE

void fft(vbase &v, bool invert = false)
{
	int n = v.size();

	if (n == 1)
		return;

	vbase v0(n / 2), v1(n / 2);

	for (int i = 0; i < n / 2; ++i)
	{
		v0[i] = v[i * 2];
		v1[i] = v[i * 2 + 1];
	}

	fft(v0, invert);
	fft(v1, invert);

	double angle = 2.0 * M_PI / n * (invert ? -1 : 1);
	base w(1), wn(cos(angle), sin(angle));

	for (int i = 0; i < n / 2; ++i)
	{
		v[i] = v0[i] + w * v1[i];
		v[i + n / 2] = v0[i] - w * v1[i];

		if (invert)
			v[i] /= 2.0, v[i + n / 2] /= 2.0;

		w *= wn;
	}
}

void fft_inverted(vbase &v)
{
	fft(v, true);
}

void multiply_polynomials_fft(const vll &lhs, const vll &rhs, vll &result)
{
	vbase base_lhs(lhs.begin(), lhs.end()), base_rhs(rhs.begin(), rhs.end());

	int n = 1;
	while (n < std::max(lhs.size(), rhs.size()))
		n <<= 1;
	n <<= 1;

	base_lhs.resize(n), base_rhs.resize(n);

	fft(base_lhs), fft(base_rhs);
	for (int i = 0; i < n; ++i)
		base_lhs[i] *= base_rhs[i];

	fft_inverted(base_lhs);

	result.resize(n);
	for (int i = 0; i < n; ++i)
		result[i] = static_cast<ll>(base_lhs[i].real() + 0.5); // 0.5 for precision
}

#else

void fft(vbase &v, const vbase &roots, const vint &rev)
{
	int n = v.size();

	for (int i = 0; i < n; ++i)
		if (i < rev[i])
			std::swap(v[i], v[rev[i]]);

	for (int len = 2; len <= n; len <<= 1)
	{
		for (int i = 0; i < n; i += len)
		{
			for (int j = 0; j < len / 2; ++j)
			{
				base a = v[i + j], b = v[i + j + len / 2] * roots[j + len / 2];

				v[i + j] = a + b;
				v[i + j + len / 2] = a - b;
			}
		}
	}
}

void multiply_polynomials_fft(const vll &lhs, const vll &rhs, vll &result, const unsigned int n, const vbase &roots,
                              const vint &reverse)
{
	vbase base_lhs(lhs.begin(), lhs.end()), base_rhs(rhs.begin(), rhs.end());

	base_lhs.resize(n), base_rhs.resize(n);

	fft(base_lhs, roots, reverse), fft(base_rhs, roots, reverse);
	for (int i = 0; i < n; ++i)
		base_lhs[i] *= base_rhs[i];
	std::reverse(base_lhs.begin() + 1, base_lhs.end());

	fft(base_lhs, roots, reverse);
	for (int i = 0; i < n; ++i)
		base_lhs[i] /= n;

	result.resize(n);
	for (int i = 0; i < n; ++i)
		result[i] = static_cast<ll>(base_lhs[i].real() + 0.5); // 0.5 for precision
}

void multiply_polynomials_fft(const vll &lhs, const vll &rhs, vll &result)
{
	int n = 1;
	while (n < std::max(lhs.size(), rhs.size()))
		n <<= 1;
	n <<= 1;

	vbase roots = calculate_roots(n);
	vint reverse = calculate_reverse(n);

	multiply_polynomials_fft(lhs, rhs, result, n, roots, reverse);
}

#endif
} // namespace frpoly

#endif // FRPOLY__FFTP_HPP
