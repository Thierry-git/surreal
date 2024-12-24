#pragma once

#include <concepts>

#include <map>

#include <ostream>
#include <iomanip>




template <typename G>
concept additiveCommutativeGroup = requires (G a, G b, G c) {
	// Closure of addition
	{ a + b } -> std::same_as<G>;

	// Identity
	G(0) + a == a;

	// Inverses
	{ -a } -> std::same_as<G>;
	a - b == a + (-b);
	a - a == G(0);

	// Associativity
	a + (b + c) == (a + b) + c;

	// Commutativity
	a + b == b + a;
};

template <typename G>
concept totallyOrderedAdditiveCommutativeGroup =
	additiveCommutativeGroup<G> &&
	std::totally_ordered<G>;

template <typename K>
concept field = additiveCommutativeGroup<K>
	&& requires (K x, K y, K z) {
	// Closure of multiplication
	{ x * y } -> std::same_as<K>;

	// Identity
	K(1) * x == x;

	// Inverses
	{ x / y } -> std::same_as<K>;
	x * (K(1) / y) == x / y;
	x != K(0) ? x / x == K(1) : true; // Exist for x =/= 0

	// Associativity
	x * (y * z) == (x * y) * z;

	// Commutativity
	x * y == y * x;

	// Distributivity
	x * (y + z) == x * y + x * z;
};

template <typename T>
concept streamable = requires (std::ostream & os, const T & t) {
	{ os << t } -> std::same_as<std::ostream&>;
};


// TODO: enforce max size
template <field K, totallyOrderedAdditiveCommutativeGroup G>
class HahnSeries {
private:
	using Numeral = std::pair<G, K>;
	using CantorNormalForm = std::map<G, K>;

	// Member attributes
	CantorNormalForm cnf_;

public:
	// CONSTRUCTORS

	// 0 ∈ K[[ω^G]]
	HahnSeries()
	: HahnSeries(0, 0) {}
	// K c___> K[[ω^G]]
	// k |---> k
	HahnSeries(K k)
	: HahnSeries(k, 0) {}
	// K x G c___> K[[ω^G]]
	// (k,g) |---> kω^g
	HahnSeries(K k, G g);


	// ARITHMETIC

	//    2w^3 + 1w^1
	//+  -1w^3 + 4w^0
	//=   1w^3 + 1w^1 + 4w^0
	HahnSeries operator+(const HahnSeries& other) const {
		auto selfNumeral = cbegin();
		auto otherNumeral = other.cbegin();
		const auto selfEnd = cend();
		const auto otherEnd = other.cend();

		HahnSeries result;

		while (selfNumeral != selfEnd && otherNumeral != otherEnd) {
			if (selfNumeral->first > otherNumeral->first) {
				result.cnf_.insert(*selfNumeral);
				++selfNumeral;
			} else if (otherNumeral->first > selfNumeral->first) {
				result.cnf_.insert(*otherNumeral);
				++otherNumeral;
			} else { // same exponent
				K coefficient = selfNumeral->second + otherNumeral->second;
				if (coefficient != K(0))
					result.cnf_.insert({selfNumeral->first, coefficient});
				++selfNumeral;
				++otherNumeral;
			}
		}
		while (selfNumeral != selfEnd) {
			result.cnf_.insert(*selfNumeral);
			++selfNumeral;
		}
		while (otherNumeral != otherEnd) {
			result.cnf_.insert(*otherNumeral);
			++otherNumeral;
		}

		return result;
	}


	// DISPLAY

	friend std::ostream& operator<<(std::ostream& os, const HahnSeries<K, G>& hs)
		requires streamable<K> && streamable<G> {
		// hs = 0
		if (hs.cnf_.empty())
			return os << K(0);
		
		// hs ∈ K[[ω^G]] - K
		auto numeral = hs.cbegin();
		auto end = hs.cend();

		if (numeral != end) {
			hs.streamOutCoefficient(os, numeral->second)
			   << "ω^" << numeral->first;
			++numeral;
		}

		for (; numeral != end; ++numeral)
		{
			os << " + ";
			hs.streamOutCoefficient(os, numeral->second)
			   << "ω^" << numeral->first;
		}
		
		return os;
	}


	// ITERATION
	
	// Order is flipped to correspond to exploring higher numerals first
	auto begin()        { return cnf_.rbegin();  }
	auto cbegin() const { return cnf_.crbegin(); }
	auto end()          { return cnf_.rend();    }
	auto cend()   const { return cnf_.crend();   }


private:
	// Display helper
	static std::ostream& streamOutCoefficient(std::ostream& os, const K& k)
		requires streamable<K> {
		if constexpr (std::is_floating_point_v<K>)
			return os << std::fixed << std::setprecision(2) << k;
		return os << k;
	}
};


//////////////////////////////////// IMPLEMENTATION ////////////////////////////////////

template <field K, totallyOrderedAdditiveCommutativeGroup G>
HahnSeries<K, G>::HahnSeries(K k, G g) {
	if (k != K(0)) cnf_.emplace(g, k);
	// else an empty cnf_ represents 0
}

