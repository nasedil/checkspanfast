#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

//=======================================================

class Timewatch
{
public:
	clock_t timestamp;
	Timewatch();
	float watch();
};

//=======================================================

Timewatch::Timewatch()
{
	timestamp = clock();
}

float Timewatch::watch()
{
	float elapsed_secs = float(clock() - timestamp) / CLOCKS_PER_SEC;
	timestamp = clock();
	return elapsed_secs;
}

//=======================================================

#endif // UTILS_HPP_INCLUDED
