#include <compare>
#include <cmath>
#include <ostream>

struct Real {
    float val_;

    constexpr Real(float v) : val_(v) {}

    // Extend to NaN
    std::strong_ordering operator<=>(const Real& other) const {
        if (std::isnan(val_) && std::isnan(other.val_)) return std::strong_ordering::equal;
        if (std::isnan(val_)) return std::strong_ordering::less;
        if (std::isnan(other.val_)) return std::strong_ordering::greater;
        
        if (val_ < other.val_) return std::strong_ordering::less;
        if (val_ > other.val_) return std::strong_ordering::greater;
        return std::strong_ordering::equal;
    }

    bool operator==(const Real& other) const {
        // Extend to NaN
        if (std::isnan(val_) || std::isnan(other.val_)) return true;
        return val_ == other.val_;
    }

    Real& operator+=(const Real& other) { val_ += other.val_; return *this; }
    Real& operator-=(const Real& other) { val_ -= other.val_; return *this; }
    Real& operator*=(const Real& other) { val_ *= other.val_; return *this; }
    Real& operator/=(const Real& other) { val_ /= other.val_; return *this; }

    Real operator-() const { return Real(-val_); }

    Real operator+(const Real& other) const { return Real(*this) += other; }
    Real operator-(const Real& other) const { return Real(*this) -= other; }
    Real operator*(const Real& other) const { return Real(*this) *= other; }
    Real operator/(const Real& other) const { return Real(*this) /= other; }

    friend std::ostream& operator<<(std::ostream& os, const Real& r) {
        return os << r.val_;
    }
};
