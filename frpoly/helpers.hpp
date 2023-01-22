#ifndef FRPOLY__HELPERS_HPP
#define FRPOLY__HELPERS_HPP

#include "pch.hpp"

namespace frpoly
{
ll binpow_mod(ll num, ll power, ll mod)
{
	ll result = 1;

	while (power != 0)
	{
		if (power & 1ll)
			result = (result * num) % mod, --power;
		else
			num = (num * num) % mod, power >>= 1;
	}

	return result;
}

ll find_group_gen(ll p)
{
	vll factors{};
	ll phi = p - 1, n = phi;
	for (ll i = 2; i * i <= n; ++i)
	{
		if (n % i == 0)
		{
			factors.push_back(i);

			while (n % i == 0)
				n /= i;
		}
	}

	if (n > 1)
		factors.push_back(n);

	for (ll gen = 2; gen <= p; ++gen)
	{
		bool is_gen_good = true;

		for (int i = 0; i < factors.size() && is_gen_good; ++i)
			is_gen_good &= (binpow_mod(gen, phi / factors[i], p) != 1);

		if (is_gen_good)
			return gen;
	}

	return -1;
}

vint calculate_reverse(unsigned int n) // зеркалим биты от 0
{
	vint reverse(n);

	int log2_n = 32 - __builtin_clz(n) - 1;

	for (int num = 0; num < n; ++num)
		for (int i = 0; i < log2_n; ++i)
			if (num & (1 << i))
				reverse[num] |= 1 << (log2_n - 1 - i);

	return reverse;
}

template <typename T>
void multiply_polynomials_stupid(const std::vector<T> &lhs, const std::vector<T> &rhs, std::vector<T> &result)
{
	result.resize(lhs.size() + rhs.size());

	for (int i = 0; i < lhs.size(); ++i)
		for (int j = 0; j < rhs.size(); ++j)
			result[i + j] += lhs[i] * rhs[j];
}

} // namespace frpoly

#endif // FRPOLY__HELPERS_HPP