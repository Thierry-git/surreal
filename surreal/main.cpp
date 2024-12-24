#include "hahnseries.h"
#include <iostream>

using namespace std;

int main() {
	HahnSeries<float, float> h1(5, 1);
	HahnSeries<float, float> h2(1, 2);
	HahnSeries<float, float> h3(-1, 3.0/2);
	
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
	     << (HahnSeries<float, float>(1, 4) - HahnSeries<float, float>(10, 3)
		    < (h1 - h2) * (h1 - h2))
		 << endl;

	cout << "3.5 * h1 = " << (3.5f * h1) << endl;
	cout << "h1 * 3.5 = " << (h1 * 3.5f) << endl;

	return 0;
}

