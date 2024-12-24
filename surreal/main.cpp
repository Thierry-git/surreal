#include "hahnseries.h"
#include <iostream>

using namespace std;

int main() {
	HahnSeries<float, float> h1(5, 1);
	HahnSeries<float, float> h2(1, 2);
	HahnSeries<float, float> h3(-1, 2);

	using Hahn = HahnSeries<float, HahnSeries<float, float>>;

	return 0;
}

