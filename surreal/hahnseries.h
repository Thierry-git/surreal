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


//////////////////////////////////// CLASS ////////////////////////////////////
/**
 * @brief Represents a Hahn series in the Hahn series field K[[ω^G]].
 * 
 * A Hahn series is a formal power series with well-ordered support,
 * where coefficients are from a field K and exponents from a totally ordered
 * abelian group G.
 * 
 * @tparam K The coefficient field type.
 * @tparam G The exponent group type.
 */
template <field K, totallyOrderedAdditiveCommutativeGroup G>
class HahnSeries {
private:
	using Numeral = std::pair<G, K>;
	using CantorNormalForm = std::map<G, K>;

	/** @brief Internal representation as map from exponents to coefficients */
	CantorNormalForm cnf_;

public:
	// CONSTRUCTORS

	/**
	 * @brief Constructs the zero series 0 ∈ K[[ω^G]].
	 */
	HahnSeries()
	: HahnSeries(0, 0) {}
	/**
	 * @brief Embeds an element of K into K[[ω^G]] as kω^0.
	 * 
	 * @param k The coefficient to embed
	 */
	HahnSeries(K k)
	: HahnSeries(k, 0) {}
	/**
	 * @brief Constructs a monomial kω^g in K[[ω^G]].
	 * 
	 * @param k The coefficient.
	 * @param g The exponent.
	 */
	HahnSeries(K k, G g){
		if (k != K(0)) cnf_.emplace(g, k);
		// else an empty cnf_ represents 0
	}


	// COMPARISON

	/**
	 * @brief Compares two Hahn series lexicographically by leading term
	 * 
	 * @param other The series to compare with
	 * @return std::strong_ordering The ordering relationship between the series
	 */
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
			if (selfNumeral->second != otherNumeral->second)
				return selfNumeral->second <=> otherNumeral->second;
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

	/**
	 * @brief Checks if two Hahn series are equal.
	 * 
	 * @param other The series to compare with.
	 * @return true If the series are equal.
	 * @return false If the series are not equal.
	 */
	bool operator==(const HahnSeries& other) const {
		return cnf_ == other.cnf_;
	}


	// ARITHMETIC

	/**
	 * @brief Adds another Hahn series to this one.
	 * 
	 * @param other The series to add.
	 * @return HahnSeries& Reference to this series after addition.
	 */
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

	/**
	 * @brief Adds two Hahn series.
	 * 
	 * @param lhs The first series.
	 * @param rhs The second series.
	 * @return HahnSeries The sum of the two series.
	 */
	friend HahnSeries operator+(const HahnSeries& lhs, const HahnSeries& rhs) {
		return HahnSeries(lhs) += rhs;
	}

	/**
	 * @brief Negates this series.
	 * 
	 * @return HahnSeries The negated series.
	 */
	HahnSeries operator-() const {
		HahnSeries result;
		for (const auto& numeral : cnf_) {
			result.cnf_.emplace_hint(result.cnf_.cbegin(), numeral.first, -numeral.second);
		}
		return result;
	}

	/**
	 * @brief Subtracts another Hahn series from this one.
	 * 
	 * @param other The series to subtract.
	 * @return HahnSeries& Reference to this series after subtraction.
	 */
	HahnSeries& operator-=(const HahnSeries& other) {
		return *this += -other;
	}

	/**
	 * @brief Subtracts two Hahn series.
	 * 
	 * @param lhs The first series.
	 * @param rhs The second series.
	 * @return HahnSeries The difference of the two series.
	 */
	friend HahnSeries operator-(const HahnSeries& lhs, const HahnSeries& rhs) {
		return HahnSeries(lhs) -= rhs;
	}

	/**
	 * @brief Multiplies this series by another one.
	 * 
	 * Performs f*g by distributing each term of f over g.
	 * That is, result = Σ(fi * g) where fi are terms of f.
	 * 
	 * @param other The series to multiply by.
	 * @return HahnSeries& Reference to this series after multiplication.
	 */
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

	/**
	 * @brief Multiplies two Hahn series.
	 * 
	 * @param lhs The first series.
	 * @param rhs The second series.
	 * @return HahnSeries The product of the two series.
	 */
	friend HahnSeries operator*(const HahnSeries& lhs, const HahnSeries& rhs) {
		HahnSeries result = lhs;
		return result *= rhs;
	}

	/**
	 * @brief Multiplies this series by a scalar.
	 * 
	 * @param scalar The scalar to multiply by.
	 * @return HahnSeries& Reference to this series after multiplication.
	 */
	HahnSeries& operator*=(const K& scalar) {
		for (auto& numeral : cnf_)
			numeral.second *= scalar;
		return *this;
	}

	/**
	 * @brief Multiplies a Hahn series by a scalar.
	 * 
	 * @param scalar The scalar to multiply by.
	 * @param hs The Hahn series to multiply.
	 * @return HahnSeries The product of the two series.
	 */
	friend HahnSeries operator*(const K& scalar, const HahnSeries& hs) {
		return HahnSeries(scalar) * hs;
	}


	// DISPLAY

	/**
	 * @brief Outputs the Hahn series to an output stream.
	 * 
	 * @param os The output stream.
	 * @param hs The Hahn series to output.
	 * @return std::ostream& The output stream.
	 */
	friend std::ostream& operator<<(std::ostream& os, const HahnSeries<K, G>& hs)
		requires streamable<K> && streamable<G> {
		// hs = 0
		if (hs.cnf_.empty())
			return os << K(0);
		
		// hs ∈ K[[ω^G]] - K
		auto numeral = hs.cbegin();
		auto end = hs.cend();

		if (numeral != end) {
			hs.streamOutCoefficient(os, numeral->second, numeral->first != G(0));
			if (numeral->first != G(0)) os << "ω";
			if (numeral->first != G(1)) os << "^" << numeral->first;
			++numeral;
		}

		for (; numeral != end; ++numeral)
		{
			os << " + ";
			hs.streamOutCoefficient(os, numeral->second, numeral->first != G(0));
			if (numeral->first != G(0)) os << "ω";
			if (numeral->first != G(1)) os << "^" << numeral->first;
		}
		
		return os;
	}

	// ITERATION
	
	/**
	 * @brief Returns a reverse iterator to the beginning (highest exponent).
	 * @return Iterator to the highest exponent term.
	 */
	typename CantorNormalForm::reverse_iterator begin() { return cnf_.rbegin(); }

	/**
	 * @brief Returns a const reverse iterator to the beginning.
	 * @return Const iterator to the highest exponent term.
	 */
	typename CantorNormalForm::const_reverse_iterator cbegin() const { return cnf_.crbegin(); }

	/**
	 * @brief Returns a reverse iterator to the end (lowest exponent)
	 * @return Iterator to the lowest exponent term
	 */
	typename CantorNormalForm::reverse_iterator end() { return cnf_.rend(); }

	/**
	 * @brief Returns a const reverse iterator to the end
	 * @return Const iterator to the lowest exponent term
	 */
	typename CantorNormalForm::const_reverse_iterator cend() const { return cnf_.crend(); }

private:
	// Display helper
	static std::ostream& streamOutCoefficient(std::ostream& os, const K& k, bool exponentNonZero)
		requires streamable<K> {
		if constexpr (std::is_floating_point_v<K>)
			return os << std::fixed << std::setprecision(2) << k;
		if (k == K(1))
			return exponentNonZero ? os : os << k;
		if (k == K(-1))
			return exponentNonZero ? os << "-" : os << k;
		return os << k;
	}
};

