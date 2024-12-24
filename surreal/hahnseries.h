#pragma once

#include <concepts>

#include <map>

#include <ostream>
#include <iomanip>

//////////////////////////////////// CONCEPTS ////////////////////////////////////

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
	HahnSeries(K k, G g){
		if (k != K(0)) cnf_.emplace(g, k);
		// else an empty cnf_ represents 0
	}


	// COMPARISON

	bool operator==(const HahnSeries& other) const {
		return cnf_ == other.cnf_;
	}

	std::strong_ordering operator<=>(const HahnSeries& other) const {
		auto selfNumeral = cbegin();
		auto otherNumeral = other.cbegin();
		const auto selfEnd = cend();
		const auto otherEnd = other.cend();

		while (selfNumeral != selfEnd && otherNumeral != otherEnd) {
			if (selfNumeral->first != otherNumeral->first) {
				// Different exponents - check which is larger and the sign of coefficient
				bool less = selfNumeral->first > otherNumeral->first
					? selfNumeral->second < K(0)
					: otherNumeral->second > K(0);
				return less ? std::strong_ordering::less
				            : std::strong_ordering::greater;
			}
			if (selfNumeral->second != otherNumeral->second) {
				if (selfNumeral->second < otherNumeral->second)
					return std::strong_ordering::less;
				else
					return std::strong_ordering::greater;
			}
			++selfNumeral;
			++otherNumeral;
		}

		// If we get here, one or both series are exhausted
		if (selfNumeral == selfEnd) {
			if (otherNumeral == otherEnd)
				return std::strong_ordering::equal;
			else
				return otherNumeral->second > K(0)
					? std::strong_ordering::less
					: std::strong_ordering::greater;
		}
		return selfNumeral->second < K(0)
			? std::strong_ordering::less
			: std::strong_ordering::greater;
	}


	// ARITHMETIC

	HahnSeries& operator+=(const HahnSeries& other) {
		CantorNormalForm result;
		auto selfNumeral = cbegin();
		auto otherNumeral = other.cbegin();
		const auto selfEnd = cend();
		const auto otherEnd = other.cend();
		
		while (selfNumeral != selfEnd && otherNumeral != otherEnd) {
			if (selfNumeral->first > otherNumeral->first) {
				result.insert(result.cbegin(), *selfNumeral);
				++selfNumeral;
			} else if (selfNumeral->first < otherNumeral->first) {
				result.insert(result.cbegin(), *otherNumeral);
				++otherNumeral;
			} else { // same exponent
				K coefficient = selfNumeral->second + otherNumeral->second;
				if (coefficient != K(0))
					result.emplace_hint(result.cbegin(), selfNumeral->first, coefficient);
				++selfNumeral;
				++otherNumeral;
			}
		}
		
		// Copy remaining terms
		while (selfNumeral != selfEnd) {
			result.insert(result.cbegin(), *selfNumeral);
			++selfNumeral;
		}
		while (otherNumeral != otherEnd) {
			result.insert(result.cbegin(), *otherNumeral);
			++otherNumeral;
		}
		
		// Overwrite cnf_
		cnf_.swap(result);
		return *this;
	}

	HahnSeries operator+(const HahnSeries& other) const {
		HahnSeries result = *this;
		result += other;
		return result;
	}

	HahnSeries operator-() const {
		HahnSeries result;
		for (const auto& numeral : cnf_) {
			// Know to insert numerals at the front
			// since they stay in the same order
			// after negation
			result.cnf_.emplace_hint(result.cnf_.cbegin(),
				numeral.first, -numeral.second);
		}
		return result;
	}

	HahnSeries& operator-=(const HahnSeries& other) {
		return *this += -other;
	}

	HahnSeries operator-(const HahnSeries& other) const {
		return *this + (-other);
	}

	// Performs f*g by distributing each term of f over g
	// That is, result = Σ(fi * g) where fi are terms of f
	HahnSeries& operator*=(const HahnSeries& other) {
		// O(1) check to avoid O(n) map copy when result must be zero
		if (cnf_.empty() || other.cnf_.empty()) {
			cnf_.clear();
			return *this;
		}

		// Store the original terms and ensure cnf_ is empty for accumulation
		CantorNormalForm original = cnf_;
		cnf_.clear();
		
		for (const auto& selfNumeral : original) {
			// Create a HahnSeries representing self_term * other
			HahnSeries termProduct;
			for (const auto& otherNumeral : other.cnf_) {
				K newCoeff = selfNumeral.second * otherNumeral.second;
				if (newCoeff != K(0)) {
					G newExp = selfNumeral.first + otherNumeral.first;
					termProduct.cnf_.emplace_hint(termProduct.cnf_.cbegin(), newExp, newCoeff);
				}
			}
			*this += termProduct;
		}
		return *this;
	}

	HahnSeries operator*(const HahnSeries& other) const {
		HahnSeries result = *this;
		result *= other;
		return result;
	}
	
	// K x K[[ω^G]] ---> K[[ω^G]]
	// (k,hs) |---> kω^0 * hs
	template<typename G2, typename K2>
	friend HahnSeries<G2, K2> operator*(const K2& scalar, const HahnSeries<G2, K2>& series) {
		return HahnSeries<G2, K2>(scalar) * series;
	}
	/*
	 * Note: K[[ω^G]] x K[[ω^G]] ---> K[[ω^G]] is
	 * defined by implicit casting of K to K[[ω^G]]
	 * and yield the same result as the above.
	*/ 

	// TODO: division

	// Note: Exponentiation is not defined at the level
	// of the Hahn series field.


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
	typename CantorNormalForm::reverse_iterator       begin()        { return cnf_.rbegin();  }
	typename CantorNormalForm::const_reverse_iterator cbegin() const { return cnf_.crbegin(); }
	typename CantorNormalForm::reverse_iterator       end()          { return cnf_.rend();    }
	typename CantorNormalForm::const_reverse_iterator cend()   const { return cnf_.crend();   }


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
