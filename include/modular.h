#include <iostream>
#include <istream>
#include <memory>
#include <ostream>

class Modular {
    using T = int64_t;

  public:
    Modular() : value_(0) {}
    Modular(T value) : value_(value) { normalize(); }
    Modular(const Modular& other) { *this = other; }
    Modular& operator=(const Modular& other) {
        value_ = other.value_;

        return *this;
    }

    static void SetModulo(const T& mod_) { mod = mod_; }

    const T& GetValue() const { return value_; }

    friend bool operator==(const Modular& lhs, const Modular& rhs) {
        return lhs.value_ == rhs.value_;
    }
    friend bool operator!=(const Modular& lhs, const Modular& rhs) {
        return !(lhs == rhs);
    }
    friend bool operator<(const Modular& lhs, const Modular& rhs) {
        return lhs.value_ < rhs.value_;
    }
    friend bool operator>(const Modular& lhs, const Modular& rhs) {
        return lhs.value_ > rhs.value_;
    }

    Modular& operator+=(const Modular& other) {
        value_ += other.value_;
        if (value_ >= mod) {
            value_ -= mod;
        }

        return *this;
    }
    friend Modular operator+(const Modular& lhs, const Modular& rhs) {
        auto result(lhs);
        result += rhs;

        return result;
    }

    Modular& operator-=(const Modular& other) {
        value_ -= other.value_;
        if (value_ < 0) {
            value_ += mod;
        }

        return *this;
    }
    friend Modular operator-(const Modular& lhs, const Modular& rhs) {
        auto result(lhs);
        result -= rhs;

        return result;
    }

    Modular& operator*=(const Modular& other) {
        value_ = (value_ * other.value_) % mod;

        return *this;
    }
    friend Modular operator*(const Modular& lhs, const Modular& rhs) {
        auto result(lhs);
        result *= rhs;

        return result;
    }

    Modular GetInverse() const { return bin_pow(*this, mod - 2); }

    Modular& operator/=(const Modular& other) {
        *this = *this * other.GetInverse();

        return *this;
    }
    friend Modular operator/(const Modular& lhs, const Modular& rhs) {
        auto result(lhs);
        result /= rhs;

        return result;
    }

    friend std::istream& operator>>(std::istream& in, Modular& a) {
        in >> a.value_;
        a.normalize();

        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const Modular& a) {
        return out << a.value_;
    }

  private:
    static Modular bin_pow(Modular a, size_t n) {
        Modular result(1);
        while (n > 0) {
            if (n & 1) {
                result *= a;
            }
            a *= a;
            n >>= 1;
        }

        return result;
    }

    void normalize() { value_ = (value_ % mod + mod) % mod; }

    T value_;

    inline static T mod = 2;
};