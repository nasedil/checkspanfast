#ifndef SLAE_HPP_INCLUDED
#define SLAE_HPP_INCLUDED

#include <vector>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

using namespace std;

//=======================================================

int power(int a, int b)
{
    int result = a;
    while (--b)
        result *= a;
    return result;
}

//=======================================================

template <size_t N>
bool gauss_solve(vector<bitset<N>> a, bitset<N> b)
{
	int r = 0;
    for (int i = 0; i < a.size(); ++i) // for all columns
    {
        // first we exchange rows if needed, so that a[i][i] = 1
        if (!a[i][r])
        {
            int k = r+1;
            while ((k < N) && (!a[i][k]))
                ++k;
			if (k >= N)
				continue;
			else
			{
				bool tmp = b[k];
				b[k] = b[r];
				b[r] = tmp;
				for (int l = i; l < a.size(); ++l)
				{
					tmp = a[l][k];
					a[l][k] = a[l][r];
					a[l][r] = tmp;
				}
			}
        }
        // second, we make all coefficients under a[i][i] equal to 0
        for (int j = r+1; j < b.size(); ++j) // for all rows (16)
        {
            if (a[i][j])
            {
                for (int k = i; k < a.size(); ++k)
                {
                    a[k][j] = (a[k][j] != a[k][r]);
                }
                b[j] = (b[j] != b[r]);
            }
        }
        ++r;
    }
    // we chech all i from 7 to 15 for zero coefficients
    for (int i = r; i < b.size(); ++i)
    {
        if (b[i])
            return false;
    }
    return true;
}

//=======================================================

class Gauss_Presolve_Data
{
public:
    size_t n;
    int r;
    vector<int> rows_s;
    vector<int> rows_d;
    vector<bool> rows_o;
};

template <size_t N>
void gauss_presolve(vector<bitset<N>> a, Gauss_Presolve_Data & p)
{
	p.n = a.size();
    p.rows_s.clear();
    p.rows_d.clear();
    p.rows_o.clear();
    int r = 0;
    for (int i = 0; i < a.size(); ++i) // for all columns
    {
        // first we exchange rows if needed, so that a[i][i] = 1
        if (!a[i][r])
        {
            int k = r+1;
            while ((k < N) && (!a[i][k]))
                ++k;
			if (k >= N)
				continue;
			else
			{
				bool tmp;
				p.rows_s.push_back(r);
				p.rows_d.push_back(k);
				p.rows_o.push_back(true);
				for (int l = i; l < a.size(); ++l)
				{
					tmp = a[l][k];
					a[l][k] = a[l][r];
					a[l][r] = tmp;
				}
			}
        }
        // second, we make all coefficients under a[i][i] equal to 0
        for (int j = r+1; j < N; ++j) // for all rows (16)
        {
            if (a[i][j])
            {
                for (int k = i; k < a.size(); ++k)
                {
                    a[k][j] = (a[k][j] != a[k][r]);
                }
                p.rows_s.push_back(r);
                p.rows_d.push_back(j);
                p.rows_o.push_back(false);
            }
        }
        ++r;
    }
    p.r = r;
}

template <size_t N>
bool gauss_solve(Gauss_Presolve_Data & p, bitset<N> b)
{
    for (int i = 0; i < p.rows_o.size(); ++i)
    {
    	if (p.rows_o[i])
		{
			bool tmp = b[p.rows_s[i]];
			b[p.rows_s[i]] = b[p.rows_d[i]];
			b[p.rows_d[i]] = tmp;
		}
		else
		{
			b[p.rows_d[i]] = (b[p.rows_d[i]] != b[p.rows_s[i]]);
		}
    }
    for (int i = p.r; i < b.size(); ++i)
    {
        if (b[i])
            return false;
    }
    return true;
}

//=======================================================

template <size_t N>
bool binary_solve(vector<bitset<N>> a, bitset<N> b)
{
    uint_least64_t counter;
    uint_least64_t limit = power(2,a.size());
    for (counter = 1; counter < limit; ++counter)
    {
        boost::dynamic_bitset<> x(a.size(),counter);
        bitset<N> r(0);
        bool is_good = true;
        // we multiply a by x
        for (int i = 0; i < b.size(); ++i)
        {
            for (int j = 0; j < a.size(); ++j)
            {
                if (x[j])
                    r[i] = (r[i] != a[j][i]);
            }
            if (r[i] != b[i])
            {
                is_good = false;
                break;
            }
        }
        if (is_good)
            return true;
    }
    return false;
}

//=======================================================

template <size_t N>
bool binary_solve_explore(vector<bitset<N>> * a, bitset<N> * b, int depth, bitset<N> * current)
{
	//cout << "\n" << depth;
	if (depth == -1)
	{
		//cout << "\nGo up";
		for (int k = 0; k < N; ++k)
			if (current->test(k) != b->test(k))
				return false;
		return true;
	}
	else
	{
		bool result = binary_solve_explore(a, b, depth-1, current);
		if (result)
			return true;
		bitset<N> * current_add = new bitset<N>();
		for (int k = 0; k < N; ++k)
			(*current_add)[k] = current->test(k) != (*a)[depth][k];
		result = binary_solve_explore(a, b, depth-1, current_add);
		delete current_add;
		return result;
	}
}

template <size_t N>
bool binary_solve_recursive(vector<bitset<N>> a, bitset<N> b)
{
    uint_least64_t depth = a.size()-1;
    bitset<N> * current = new bitset<N>();
    current->reset();
	bool result = binary_solve_explore(&a, &b, depth, current);
    delete current;
    return result;
}

//=======================================================

#endif // SLAE_HPP_INCLUDED
