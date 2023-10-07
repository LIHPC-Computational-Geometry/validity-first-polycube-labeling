#include <geogram/basic/assert.h> // for geo_assert_not_reached

#include <algorithm>    // for std::min() and std::max()
#include <limits>       // for std::numeric_limits<>::min() and max()
#include <cmath>        // for std::std::pow()

class BasicStats {
public:

    BasicStats() {
        reset();
    }

    void reset() {
        min_ = std::numeric_limits<double>::max();
        max_ = std::numeric_limits<double>::min();
        sum_ = 0.0;
        avg_ = 0.0;
        variance_ = 0.0;
        count_ = 0;
    }

    double min() {
        if(count_ == 0) { geo_assert_not_reached; }
        return min_;
    }

    double max() {
        if(count_ == 0) { geo_assert_not_reached; }
        return max_;
    }

    double sum() {
        return sum_;
    }

    unsigned int count() {
        return count_;
    }

    double avg() {
        if(count_ == 0) { geo_assert_not_reached; }
        return avg_;
    }

    double variance() {
        if(count_ == 0) { geo_assert_not_reached; }
        return variance_;
    }

    double sd() {
        if(count_ == 0) { geo_assert_not_reached; }
        return std::sqrt(variance_);
    }

    void insert(double value) {
        count_++;
        double avg_prev = avg_;
        min_ = std::min(min_,value);
        max_ = std::max(max_,value);
        sum_ += value;
        avg_ = sum_ / (double) count_;
        // recursive formula for variance, thank you Did https://math.stackexchange.com/a/375022
        variance_ = variance_ + std::pow(avg_prev,2.0) - std::pow(avg_,2.0) + (( std::pow(value,2.0) - variance_ - std::pow(avg_prev,2.0) ) / count_);
    }

private:

    double min_;
    double max_;
    double sum_;
    double avg_; // needed for the standard deviation, so store it
    double variance_;
    unsigned int count_;
};