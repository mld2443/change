#include <algorithm>    // max, min
#include <chrono>       // duration_cast, nanoseconds, now, steady_clock, time_point
#include <cmath>        // ceil, log2, log10, lround, powf, sin
#include <complex>      // complex
#include <immintrin.h>  // _tzcnt_u64
#include <iomanip>      // put_money, showbase
#include <iostream>     // cout, endl, locale
#include <map>          // map
#include <numeric>      // reduce
#include <string>       // string
#include <tuple>        // tuple
#include <utility>      // move
#include <vector>       // size_t, vector


using namespace std;

///////////
// TIMER //
///////////
class Timer {
public:
    Timer(const string &name): name(name), start(chrono::steady_clock::now()) { cout << "Start \"" << name << "\"." << endl; }
    ~Timer() {
        const auto end = chrono::steady_clock::now();
        const uint64_t elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        const uint32_t digits = static_cast<uint32_t>(log10(elapsed));
        const uint32_t scale = min(digits / 3u, 3u);
        const float trimmed = static_cast<float>(elapsed) / powf(1000.f, static_cast<float>(scale));
        static const char *unitScales[] = {"ns", "us", "ms", "s"};
        const char *units = unitScales[scale];

        cout << "-Stop \"" << name << "\", " << trimmed << units << " elapsed." << endl;
    }

private:
    const string name;
    const chrono::time_point<chrono::steady_clock> start;
};


////////////////////////////
// Fast Fourier Transform //
////////////////////////////
namespace {
    constexpr double  PI = 3.14159265358979323846264338327;

    template <typename T>
    void butterflyShuffleInPlace(vector<T> &data) {
        if (data.size() > 2uz) {
            // Resize data to the nearest power of 2
            const size_t logN = static_cast<size_t>(ceil(log2((double)data.size())));
            const size_t N = 1uz << logN;
            data.resize(N);

            // Shuffle the first half with the second
            for (size_t i = 0uz, i_rcp = 0uz; i < data.size(); ++i) {
                if (i < i_rcp)
                    swap(data[i], data[i_rcp]);

                size_t pow2 = data.size() >> 1uz;
                while (pow2 >= 1uz && i_rcp >= pow2) {
                    i_rcp -= pow2;
                    pow2 >>= 1uz;
                }
                i_rcp += pow2;
            }
        }
    }

    static vector<complex<double>> fwdTwiddles;
    static vector<complex<double>> invTwiddles;
    static vector<complex<double>> *twiddles = nullptr;

    void prepareTwiddles(size_t N, bool inverted) {
        twiddles = (inverted ? &invTwiddles : &fwdTwiddles);
        if (N <= twiddles->size())
            return;

        twiddles->resize(N);

        // Numerically stable way to walk the roots of unity.
        // -2*sin(0.5*theta)^2 + i*sin(theta) == e^(i*theta) - 1.0
        const double theta = (inverted ? -1.0 : 1.0) * PI / static_cast<double>(N);
        const double sqrcos = sin(0.5 * theta);
        const complex principle = { -2.0 * sqrcos * sqrcos, sin(theta) };
        complex omega = { 1.0, 0.0 };

        for(auto &twiddle : *twiddles) {
            twiddle = omega;
            omega += omega * principle;
        }
    }

    const complex<double>& getTwiddle(size_t N, size_t index) {
        const size_t factor = twiddles->size() / N;
        return (*twiddles)[factor * index];
    }

    // Adapted from Numerical Recipes
    void fastFourierTransform(vector<complex<double>>& data, bool inverted) {
        butterflyShuffleInPlace(data);

        // Precompute the twiddles since many can be re-used.
        const size_t halfN = data.size() >> 1uz;
        prepareTwiddles(halfN, inverted);

        // The Danielson-Lanczos routine.
        // Evaluation loop executed logN - 1 times.
        for (size_t evalStep = 1uz, twiddleStep = halfN; evalStep < data.size(); evalStep <<= 1uz, twiddleStep >>= 1uz) {
            // Inner loops combine to execute N times.
            for (size_t m = 0uz; m < evalStep; ++m) {
                for (size_t i = m; i < data.size(); i += (evalStep << 1uz)) {
                    const size_t j = i + evalStep;
                    const complex w = getTwiddle(halfN, m * twiddleStep);
                    const complex evaluation = data[j] * w;
                    data[j] = data[i] - evaluation;
                    data[i] = data[i] + evaluation;
                }
            }
        }
    }

    template <typename T>
    vector<complex<double>> fastFourierTransform(const vector<T> &data, size_t N, bool inverted = false) {
        vector<complex<double>> complexData(data.cbegin(), data.cend());
        complexData.resize(N);

        fastFourierTransform(complexData, inverted);

        return complexData;
    }

    template <typename T>
    vector<complex<double>> fastFourierTransform(const vector<T> &data, bool inverted = false) {
        vector<complex<double>> complexData(data.cbegin(), data.cend());

        fastFourierTransform(complexData, inverted);

        return complexData;
    }
}


////////////////
// Polynomial //
////////////////
// std::vector inline multiply
template <typename T>
static inline vector<T> operator*(const vector<T> &v1, const vector<T> &v2) {
    const size_t limit = min(v1.size(), v2.size());
    vector<T> result;
    result.reserve(limit);

    for (size_t i = 0uz; i < limit; ++i)
        result.push_back(v1[i] * v2[i]);

    return result;
}

struct DenominationPolynomial {
    size_t stride, limit;

    size_t size() const { return limit - ((limit - 1uz) % stride); }
};

template <typename CoeffType>
class Polynomial {
private:
    static vector<CoeffType>& dropZeros(vector<CoeffType>& vec) {
        while (!vec.empty() && vec.back() == CoeffType(0))
            vec.pop_back();
        return vec;
    }
    static vector<CoeffType> expand(const DenominationPolynomial &d) {
        vector<CoeffType> expansion(d.size());
        for (size_t i = 0uz; i < expansion.size(); i += d.stride)
            expansion[i] = CoeffType(1);
        return expansion;
    }
    [[nodiscard]]
    static vector<CoeffType> fftMultiply(const Polynomial &p1, const Polynomial &p2) {
        // Guaranteed size (1 + degree) of the resulting polynomial
        const size_t N = p1.m_coefficients.size() + p2.m_coefficients.size() - 1uz;

        // Convert polynomials over to complex numbers in sacrificial lists
        vector<complex<double>> pointValue1(p1.m_coefficients.begin(), p1.m_coefficients.end());
        vector<complex<double>> pointValue2(p2.m_coefficients.begin(), p2.m_coefficients.end());
        vector<CoeffType> result;

        // Grow lists to the new polynomial size
        pointValue1.resize(N);
        pointValue2.resize(N);
        result.reserve(N);

        // FFT polynomial evaluation
        fastFourierTransform(pointValue1);
        fastFourierTransform(pointValue2);

        // Multiply the complex point value pairs
        pointValue1 = move(pointValue1 * pointValue2);

        // FFT polynomial interpolation
        fastFourierTransform(pointValue1, true);

        // Discard the imaginary component and normalize the result
        for (size_t i = 0uz; i < N; ++i)
            result.push_back(CoeffType(lround(pointValue1[i].real() / ((double)pointValue1.size()))));

        return result;
    }
    [[nodiscard]]
    static vector<CoeffType> dumbMultiply(const Polynomial &p1, const DenominationPolynomial &d) {
        const size_t N = p1.m_coefficients.size() + d.size() - 1uz;
        vector<CoeffType> result(N);

        const size_t stop = d.size();
        for (size_t i = 0uz; i < stop; i += d.stride)
            for (size_t j = 0uz; j < p1.m_coefficients.size(); ++j)
                result[i + j] += p1.m_coefficients[j];

        return result;
    }

public:
    Polynomial(vector<CoeffType> &&coeffs) : m_coefficients(dropZeros(coeffs)) {}
    Polynomial(const DenominationPolynomial &d) : m_coefficients(expand(d)) {}

    vector<CoeffType>& getCoefficients() { return m_coefficients; }

    Polynomial<CoeffType> operator+(const Polynomial<CoeffType> &poly) const {
        const size_t coeffSize = max(m_coefficients.size(), poly.m_coefficients.size());
        vector<CoeffType> result;
        result.reserve(coeffSize);

        for (size_t i = 0uz; i < coeffSize; ++i)
            result.push_back((i <      m_coefficients.size() ?      m_coefficients[i] : CoeffType(0)) +
                             (i < poly.m_coefficients.size() ? poly.m_coefficients[i] : CoeffType(0)));

        return {move(result)};
    }
    [[nodiscard]]
    Polynomial<CoeffType> operator*(const Polynomial<CoeffType> &poly) const { return move(fftMultiply(*this, poly)); }
    [[nodiscard]]
    Polynomial<CoeffType> operator*(const DenominationPolynomial &d) const { return move(dumbMultiply(*this, d)); }

    Polynomial<CoeffType>& operator*=(const Polynomial<CoeffType> &poly) {
        m_coefficients = fftMultiply(*this, poly);
        return *this;
    }
    Polynomial<CoeffType>& operator*=(const DenominationPolynomial &d) {
        m_coefficients = dumbMultiply(*this, d);
        return *this;
    }

    template <typename T>
    friend Polynomial<T> multiplyAll(const vector<Polynomial<T>>& polys) {
        // Base cases
        if (polys.empty()) return { vector<T>{ T(1) } };
        if (polys.size() == 1uz) return polys[0uz];

        // Precompute the size of the resulting polynomial
        size_t numCombinedCoeffs = 1uz;
        for (const Polynomial<T> poly : polys)
            numCombinedCoeffs += poly.m_coefficients.size() - 1uz;

        // Convert the first polynomial to point-value form.
        auto pointValueProduct = fastFourierTransform(polys[0].m_coefficients, numCombinedCoeffs);

        // Transform remaining polynomials into pv form and inline multiply.
        for (size_t i = 1uz; i < polys.size(); ++i)
            pointValueProduct = move(pointValueProduct * fastFourierTransform(polys[i].m_coefficients, numCombinedCoeffs));

        // Convert back to coefficient form.
        fastFourierTransform(pointValueProduct, true);

        // Discard the imaginary component and normalize the result
        vector<T> result;
        result.reserve(numCombinedCoeffs);
        for (size_t i = 0uz; i < numCombinedCoeffs; ++i)
            result.push_back(T(lround(pointValueProduct[i].real() / ((double)pointValueProduct.size()))));

        return {move(result)};
    }

private:
    vector<CoeffType> m_coefficients;
};


//////////////
// CURRENCY //
//////////////
namespace {
    uint64_t nCr(uint64_t n, uint64_t r) {
        if (n < r) return 0ul;
        if (r == n) return 1ul;

        if (r << 1ul > n)
            r = n - r;

        uint64_t accum = 1ul, numer = n - r, denom = 0ul;

        while (denom < r) {
            ++numer;
            ++denom;
            accum = numer * accum / denom;
        }

        return accum;
    }

    uint64_t GCD(uint64_t u, uint64_t v) {
        if (u == 1ul || v == 1ul) return 1ul;
        if (u == 0ul) return v;
        if (v == 0ul) return u;

        const uint64_t i = _tzcnt_u64(u); u >>= i;
        const uint64_t j = _tzcnt_u64(v); v >>= j;
        const uint64_t k = min(i, j);

        while (true) {
            if (u > v)
                swap(u, v);

            v -= u;

            if (v <= 0ul)
                return u << k;

            v >>= _tzcnt_u64(v);
        }
    }
    uint64_t GCDv(vector<uint64_t>::const_iterator begin, vector<uint64_t>::const_iterator end) {
        return reduce(begin, end, 0ul, GCD);
    }

    uint64_t LCM(uint64_t u, uint64_t v) {
        const uint64_t gcd = GCD(u, v);

        if (gcd == 0ul)
            return 0ul;
        return u / gcd * v;
    }
    uint64_t LCMv(const vector<uint64_t> &v) {
        return reduce(v.begin(), v.end(), 1ul, LCM);
    }

    template <typename T>
    vector<T> remap(const vector<T> &v, T divisor) {
        if (v.empty()) return {};

        vector<T> result;
        result.reserve(v.size());

        if (*v.cbegin() == T(1))
            result.push_back(T(1));

        for (size_t i = (*v.cbegin() == T(1)); i < v.size(); ++i)
            result.push_back(v[i] / divisor);

        return result;
    }

    template <typename T>
    vector<T> denomReduceFFTFold(const vector<T>& denoms, T lcm) {
        if (denoms.empty()) return {T(0)};

        Polynomial<T> p(DenominationPolynomial(denoms[0], lcm));
        for (size_t i = 1uz; i < denoms.size(); ++i)
            p *= Polynomial<T>(DenominationPolynomial(denoms[i], lcm));

        return p.getCoefficients();
    }

    template <typename T>
    vector<T> denomReduceFFTCombined(const vector<T>& denoms, T lcm) {
        if (denoms.empty()) return {T(0)};

        vector<Polynomial<T>> polys;
        for (size_t i = 0uz; i < denoms.size(); ++i)
            polys.push_back(Polynomial<T>(DenominationPolynomial(denoms[i], lcm)));

        return multiplyAll(polys).getCoefficients();
    }

    template <typename T>
    vector<T> denomReduceDumbFold(const vector<T>& denoms, T lcm) {
        if (denoms.empty())
            return {T(0)};

        Polynomial<T> p(DenominationPolynomial(denoms[0], lcm));
        for (size_t i = 1uz; i < denoms.size(); ++i)
            p *= DenominationPolynomial(denoms[i], lcm);

        return p.getCoefficients();
    }
}

class Currency {
public:
    Currency(const vector<uint64_t> &&d)
    : m_denominations(d)
    , m_hasUnitValue(!m_denominations.empty() && m_denominations.front() == 1ul)
    , m_gcd(GCDv(m_denominations.cbegin() + (m_hasUnitValue ? 1ul : 0ul), m_denominations.cend()))
    , m_scaledDenominations(remap(m_denominations, m_gcd))
    , m_lcm(LCMv(m_denominations))
// Change below to any of the three "denomReduce*" functions above.
    , m_coefficients(denomReduceDumbFold(m_scaledDenominations, m_lcm / m_gcd))
    , m_memo() {}

    uint64_t analyticWaysToGiveChange(uint64_t amount) const {
        if (!m_hasUnitValue && amount % m_gcd != 0ul)
            return 0ul;

        const uint64_t r = m_denominations.size() - 1ul;
        const uint64_t limit = (m_coefficients.size()) * m_gcd;
        uint64_t count = 0ul;
        for (uint64_t exponent = amount % m_lcm; exponent < limit && exponent <= amount; exponent += m_lcm)
            count += m_coefficients[exponent / m_gcd] * nCr(((amount - exponent) / m_lcm) + r, r);

        return count;
    }
    uint64_t dynamicWaysToGiveChange(uint64_t amount) const {
        return dynamicWaysToGiveChange(amount, m_denominations.size() - 1ul);
    }

    void outputSizes() const {
        cout << "memo size: " << m_memo.size() * sizeof(*m_memo.begin()) << " B, coeffs size: " << m_coefficients.size() * sizeof(uint64_t) << " B" << endl;
    }

private:
    uint64_t dynamicWaysToGiveChange(uint64_t amount, uint64_t index) const {
        while (m_denominations[index] > amount && index > 0ul) --index;

        if (index == 0ul)
            return amount % m_denominations[0ul] == 0ul;
        if (auto search = m_memo.find({ amount, index }); search != m_memo.end())
            return search->second;

        uint64_t count = 0ul;
        for (uint64_t step = 0ul; step <= amount; step += m_denominations[index])
            count += dynamicWaysToGiveChange(amount - step, index - 1ul);

        m_memo[{ amount, index }] = count;
        return count;
    }

public:
    const vector<uint64_t>                  m_denominations;
    const bool                               m_hasUnitValue;
    const uint64_t                                    m_gcd;
    const vector<uint64_t>            m_scaledDenominations;
    const uint64_t                                    m_lcm;
    const vector<uint64_t>                   m_coefficients;
private:
    mutable map<tuple<uint64_t, uint64_t>, uint64_t> m_memo;
};


//////////
// MAIN //
//////////
int main() {
    Timer *tp = new Timer{"euro currency constructor"};
    //Currency uscoins({ 1, 5, 10, 25, 50, 1'00 });
    //Currency commonDollar({ 1, 5, 10, 25, 1'00, 5'00, 10'00, 20'00, 50'00, 100'00 });
    //Currency fullDollar({ 1, 5, 10, 25, 50, 1'00, 2'00, 5'00, 10'00, 20'00, 50'00, 100'00 });
    Currency euro({ 1, 2, 5, 10, 20, 50, 1'00, 2'00, 5'00, 10'00, 20'00, 50'00, 100'00, 200'00, 500'00 });
    //Currency primes({ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41 });
    //Currency powers({ 1<<0, 1<<1, 1<<2, 1<<3, 1<<4, 1<<5, 1<<6, 1<<7, 1<<8, 1<<9, 1<<10, 1<<11, 1<<12 });
    //Currency _5_7_13_23({5, 7, 13, 23});
    delete tp;
    tp = nullptr;

    Currency *selected = &euro;

    string loc = "en_US.UTF-8";
    cout << "Switching to locale '" << loc << "'" << endl;
    cout.imbue(locale(loc));

    const uint64_t bigNumber = 399'99ul;
    {
        Timer _{"Analytic"};
        cout << showbase << put_money(static_cast<long double>(bigNumber)) << " -> " << selected->analyticWaysToGiveChange(bigNumber) << endl;
    }
    {
        Timer _{"Dynamic"};
        cout << showbase << put_money(static_cast<long double>(bigNumber)) << " -> " << selected->dynamicWaysToGiveChange(bigNumber) << endl;
    }
    selected->outputSizes();
    cout.imbue(locale());
    return 0;
}
