#include "hahnseries.h"
#include <iostream>
#include "real.h"

using namespace std;

int main() {
	/////////////////// LEVEL 0 ///////////////////////
	using Hahn0 = HahnSeries<Real, Real>;
	
	Hahn0 h1(5.0, 1.0);
	Hahn0 h2(1.0, 2.0);
	Hahn0 h3(-1.0, 3.0/2);
	
	// using Hahn = HahnSeries<float, HahnSeries<float, float>>;

	cout << "h1 = " << h1 << endl;
	cout << "h2 = " << h2 << endl;
	cout << "h3 = "	<< h3 << endl;
	cout << "-h2 = " << -h2 << endl;
	cout << "h2 - h2 = " << h2 - h2 << endl;

	cout << "h1 + h2 = " << h1 + h2 << endl;
	cout << "h1 * h2 = " << h1 * h2 << endl;
	cout << "h1 * h3 = " << h1 * h3 << endl;

	cout << "(h1 - h2)^2 = " << (h1 - h2) * (h1 - h2) << endl;
	
	cout << boolalpha << "h1 < h2 = " << (h1 < h2) << endl;
	cout << boolalpha << "h1 < h3 = " << (h1 < h3) << endl;

	cout << boolalpha << "(1ω^4 -10ω^3) < (h1 - h2)^2 = "
	     << (HahnSeries<Real, Real>(1.0, 4.0) - HahnSeries<Real, Real>(10.0, 3.0)
		    < (h1 - h2) * (h1 - h2))
		 << endl;

	cout << "3.5 * h1 = " << (3.5 * h1) << endl;
	cout << "h1 * 3.5 = " << (h1 * Real(3.5)) << endl;

	/////////////////// LEVEL 1 ///////////////////////
	using Hahn1 = HahnSeries<Real, Hahn0>;
	
	Hahn1 h4(5.0, {1.0, 1.0});

	cout << "h4 = " << h4 << endl;

	return 0;
}

