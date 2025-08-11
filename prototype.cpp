// hahn_series_prototype.hpp — header-only prototype
// Minimal, eager Hahn-series core with exact Rational<K> and generic exponent G.
// Focus: clean invariants, efficient merges, leading-term API, scalar & series ops.
// Note: For the zero series, leading_exponent() is undefined (throws),
//       but leading_coefficient() returns 0 by convention.
// Build demo: define HAHN_SERIES_DEMO in one TU and compile.

// #pragma once

#include <map>
#include <numeric>
#include <compare>
#include <cstdint>
#include <ostream>
#include <stdexcept>

// -------------------- Small exact Rational (int64) --------------------
class Rational {
    std::int64_t num_ {0};
    std::int64_t den_ {1}; // > 0 always
    static std::int64_t sgn(std::int64_t x) { return (x>0) - (x<0); }
    static std::int64_t gcd64(std::int64_t a, std::int64_t b) {
        if (a<0) a = -a; if (b<0) b = -b;
        while (b){ auto t = a % b; a = b; b = t; } return a ? a : 1;
    }
    void normalize() {
        if (den_ == 0) throw std::domain_error("Rational: zero denominator");
        if (num_ == 0) { den_ = 1; return; }
        auto g = gcd64(num_, den_);
        num_ /= g; den_ /= g;
        // make denominator positive
        if (den_ < 0) { num_ = -num_; den_ = -den_; }
    }
public:
    constexpr Rational() = default;
    constexpr Rational(std::int64_t n) : num_(n), den_(1) {}
    Rational(std::int64_t n, std::int64_t d) : num_(n), den_(d) { normalize(); }

    std::int64_t num() const { return num_; }
    std::int64_t den() const { return den_; }

    // arithmetic
    Rational& operator+=(const Rational& r) {
        num_ = num_*r.den_ + r.num_*den_;
        den_ = den_*r.den_;
        normalize();
        return *this;
    }
    Rational& operator-=(const Rational& r) {
        num_ = num_*r.den_ - r.num_*den_;
        den_ = den_*r.den_;
        normalize();
        return *this;
    }
    Rational& operator*=(const Rational& r) {
        num_ *= r.num_; den_ *= r.den_;
        normalize();
        return *this;
    }
    Rational& operator/=(const Rational& r) {
        if (r.num_ == 0) throw std::domain_error("Rational: division by 0");
        num_ *= r.den_; den_ *= r.num_;
        normalize();
        return *this;
    }

    friend Rational operator+(Rational a, const Rational& b) { return a += b; }
    friend Rational operator-(Rational a, const Rational& b) { return a -= b; }
    friend Rational operator*(Rational a, const Rational& b) { return a *= b; }
    friend Rational operator/(Rational a, const Rational& b) { return a /= b; }
    Rational operator-() const { return Rational(-num_, den_); }

    // comparisons
    friend bool operator==(const Rational& a, const Rational& b) {
        return a.num_ == b.num_ && a.den_ == b.den_;
    }

    friend std::strong_ordering operator<=>(const Rational& a, const Rational& b) {
        // compare a.num/a.den ? b.num/b.den => cross-multiply (sign-safe)
        __int128 lhs = static_cast<__int128>(a.num_) * b.den_;
        __int128 rhs = static_cast<__int128>(b.num_) * a.den_;
        if (lhs < rhs) return std::strong_ordering::less;
        if (lhs > rhs) return std::strong_ordering::greater;
        return std::strong_ordering::equal;
    }

    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        if (r.den_ == 1) return os << r.num_;
        return os << r.num_ << '/' << r.den_;
    }
};

// Traits for zero/one (customize if needed for other fields)
template<class T> struct field_traits {
    static constexpr T zero() { return T(0); }
    static constexpr T one() { return T(1); }
};

template<> struct field_traits<Rational> {
    static constexpr Rational zero() { return Rational(0); }
    static constexpr Rational one() { return Rational(1); }
};

// Simple syntactic concepts (don’t try to assert algebraic laws at compile time)
#include <type_traits>

template<class T>
concept Ordered = requires (const T& a, const T& b){ { a <=> b } -> std::same_as<std::strong_ordering>; };

template<class T>
concept Addable = requires (T a, T b) { a+b; a-b; };

template<class T>
concept Multipliable = requires (T a, T b) { a*b; a/b; };

// -------------------- HahnSeries<K,G> --------------------
// Invariants: (1) no zero coefficients, (2) unique exponents, (3) map ordered ascending by G

template <class K, class G, class Compare = std::less<G>>
requires Ordered<G> && Addable<G> && Addable<K> && Multipliable<K>
class HahnSeries {
public:
    using coeff_type   = K;
    using exp_type     = G;
    using map_type     = std::map<G,K,Compare>;
    using const_iterator = typename map_type::const_iterator;
    using const_riterator = typename map_type::const_reverse_iterator;

private:
    map_type terms_; // ascending by exponent

    static constexpr K kzero() { return field_traits<K>::zero(); }

public:
    // ---- ctors ----
    HahnSeries() = default; // zero
    explicit HahnSeries(const K& k) { if (k != kzero()) terms_.emplace(G{}, k); }
    HahnSeries(const K& k, const G& g) { if (k != kzero()) terms_.emplace(g, k); }

    bool empty() const noexcept { return terms_.empty(); }
    bool is_zero() const noexcept { return terms_.empty(); }
    std::size_t size() const noexcept { return terms_.size(); }

    // leading term API (highest exponent)
    const_riterator rbegin() const { return terms_.rbegin(); }
    const_riterator rend()   const { return terms_.rend(); }

    // For zero series: exponent is undefined; coefficient is defined as 0.
    const G& leading_exponent() const { if (empty()) throw std::logic_error("leading_exponent of 0"); return terms_.rbegin()->first; }
    K leading_coefficient() const { return empty() ? kzero() : terms_.rbegin()->second; }

    // ---- equality ----
    friend bool operator==(const HahnSeries& a, const HahnSeries& b) { return a.terms_ == b.terms_; }

    // ---- comparison (lex by leading term) ----
    friend std::strong_ordering operator<=>(const HahnSeries& a, const HahnSeries& b) {
        auto i = a.terms_.rbegin(), e = a.terms_.rend();
        auto j = b.terms_.rbegin(), f = b.terms_.rend();
        const auto zero = field_traits<K>::zero();

        while (i!=e || j!=f) {
            // a has the current largest exponent
            if (i!=e && (j==f || i->first > j->first)){
                if (i->second > zero) return std::strong_ordering::greater;
                if (i->second < zero) return std::strong_ordering::less;
                ++i; continue;
            }
            // b has the current largest exponent
            if (j!=f && (i==e || j->first > i->first)){
                if (j->second > zero) return std::strong_ordering::less;
                if (j->second < zero) return std::strong_ordering::greater;
                ++j; continue;
            }
            // same exponent → compare coefficients
            if (auto ord = (i->second <=> j->second); ord != 0) return ord;
            ++i; ++j;
        }
        return std::strong_ordering::equal;
    }

    // ---- addition ----
    HahnSeries& operator+=(const HahnSeries& o) {
        map_type out;
        auto i = terms_.begin(), e = terms_.end();
        auto j = o.terms_.begin(), f = o.terms_.end();
        while (i!=e && j!=f) {
            if (i->first < j->first) { out.emplace_hint(out.end(), i->first, i->second); ++i; }
            else if (j->first < i->first) { out.emplace_hint(out.end(), j->first, j->second); ++j; }
            else { auto c = i->second + j->second; if (c != kzero()) out.emplace_hint(out.end(), i->first, c); ++i; ++j; }
        }
        for (; i!=e; ++i) out.emplace_hint(out.end(), i->first, i->second);
        for (; j!=f; ++j) out.emplace_hint(out.end(), j->first, j->second);
        terms_.swap(out);
        return *this;
    }
    friend HahnSeries operator+(HahnSeries a, const HahnSeries& b) { return a += b; }

    // ---- negation / subtraction ----
    HahnSeries operator-() const {
        HahnSeries r; for (auto&& [g,c] : terms_) r.terms_.emplace(g, -c); return r;
    }
    HahnSeries& operator-=(const HahnSeries& o) { return *this += (-o); }
    friend HahnSeries operator-(HahnSeries a, const HahnSeries& b) { return a -= b; }

    // ---- scalar multiply ----
    HahnSeries& operator*=(const K& s) {
        if (s == kzero()) { terms_.clear(); return *this; }
        for (auto it = terms_.begin(); it != terms_.end(); ) {
            it->second *= s;
            if (it->second == kzero()) it = terms_.erase(it); else ++it;
        }
        return *this;
    }
    friend HahnSeries operator*(const K& s, HahnSeries h) { h *= s; return h; }
    friend HahnSeries operator*(HahnSeries h, const K& s) { h *= s; return h; }

    // ---- series multiply (naive, combines on the fly) ----
    HahnSeries& operator*=(const HahnSeries& o) {
        if (terms_.empty() || o.terms_.empty()) { terms_.clear(); return *this; }
        map_type out;
        for (auto&& [g1,c1] : terms_) {
            for (auto&& [g2,c2] : o.terms_) {
                K c = c1 * c2; if (c == kzero()) continue;
                G g = g1 + g2;
                auto it = out.lower_bound(g);
                if (it != out.end() && !(out.key_comp()(g, it->first))) {
                    it->second += c;
                    if (it->second == kzero()) out.erase(it);
                } else {
                    out.emplace_hint(it, g, c);
                }
            }
        }
        terms_.swap(out);
        return *this;
    }
    friend HahnSeries operator*(HahnSeries a, const HahnSeries& b) { return a *= b; }

    // ---- streaming ----
    friend std::ostream& operator<<(std::ostream& os, const HahnSeries& h) {
        if (h.terms_.empty()) return os << field_traits<K>::zero();
        bool first = true;
        for (auto it = h.terms_.rbegin(); it != h.terms_.rend(); ++it) {
            const auto& g = it->first; const auto& c = it->second;
            if (!first) os << " + ";
            first = false;
            os << c << "*ω^" << g;
        }
        return os;
    }
};

// -------------------- Demo --------------------
#ifdef HAHN_SERIES_DEMO
#include <iostream>
int main(){
    using K = Rational;  // exact coefficients
    using G = Rational;  // simple totally ordered additive group
    using H = HahnSeries<K,G>;

    H a(K(5), G(1));      // 5*ω^1
    H b(K(1), G(2));      // 1*ω^2
    H c(K(-1), G(3));     // -1*ω^3

    std::cout << "a      = " << a << "\n";
    std::cout << "b      = " << b << "\n";
    std::cout << "c      = " << c << "\n";
    std::cout << "a+b    = " << (a+b) << "\n";
    std::cout << "a*b    = " << (a*b) << "\n";
    std::cout << "a*c    = " << (a*c) << "\n";
    std::cout << "3*a    = " << (K(3)*a) << "\n";
    std::cout << "(a+b)² = " << (a+b)*(a+b) << "\n";
    std::cout << "a²+b²  = " << a*a + b*b << "\n";
    std::cout << "a²+b²-(a+b)² = " << a*a + b*b - (a+b)*(a+b) << "\n";
    H e;
    std::cout << "0: " << e << "\n";
    e += H(K(1));
    for (int i = 0; i < 5; i++) { e *= a+b+c; }
    std::cout << "e = (a-c+b)⁵ = " << e << "\n";

    auto cmp = (a <=> b);
    std::cout << "(a < b) = " << std::boolalpha << (cmp == std::strong_ordering::less) << "\n";

    // Nested exponents prototype: use H as exponent (compile-time heavy; practical builds should add truncation policies)
    using H1 = HahnSeries<K,H>; // coefficients in K, exponents in H
    H1 d(K(7), a); // 7*ω^a
    std::cout << "d      = " << d << "\n";

    // TODO : include H0 in H1, Hn in H(n+1) ...
    /*
    cmp = (e <=> d);
    std::cout << "(e < d) =" << std::boolalpha << (cmp == std::strong_ordering::less) << "\n";
    */
}
#endif
