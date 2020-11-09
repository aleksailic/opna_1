#define CATCH_CONFIG_RUNNER

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>

#include <boost/program_options.hpp>
#include <boost/rational.hpp>
#include <boost/optional/optional.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>

#include <catch.hpp>
#include <prettyprint.hpp>
#include <fort.hpp>


static const char *GetProgramName(const char *path) {
    const char *last = path;
    while (*path++) {
        last = (*path == '/') ? path : last;
    }
    return last;
}

TEST_CASE("GetProgramName") {
    constexpr auto kExamplePath1 = "/home/ilic/opna-1/opna_1";
    constexpr auto kExamplePath2 = "/opna_1";
    constexpr auto kExpectedName = "/opna_1";

    REQUIRE(strcmp(GetProgramName(kExamplePath1), kExpectedName) == 0);
    REQUIRE(strcmp(GetProgramName(kExamplePath2), kExpectedName) == 0);
}

namespace mp = boost::multiprecision;

using IntType = long long;
using FloatType = mp::cpp_dec_float_100;
using Fraction = boost::rational<IntType>;
using DecimalNumber = std::tuple<IntType, FloatType>;

std::ostream &operator<<(std::ostream &os, const DecimalNumber &number) {
    return os << '[' << std::get<0>(number) << ',' << std::get<1>(number) << ']';
}

boost::optional<std::string> SafeSubstring(const std::string &string, size_t pos, size_t len = std::string::npos) {
    if (pos >= string.length())
        return boost::none;
    else
        return string.substr(pos, len);
}

DecimalNumber ParseDecimalNumber(std::string number) {
    auto dot_index = number.find('.');
    if (dot_index == std::string::npos) {
        return std::make_tuple(std::stol(number), 0);
    }
    std::string a = number.substr(0, dot_index);
    std::string d = "0." + SafeSubstring(number, dot_index + 1).get_value_or("0");

    return std::make_tuple(std::stol(a), FloatType(d));
}

DecimalNumber ParseDecimalNumber(FloatType number) {
    long a = std::floor(static_cast<long double>(number));
    auto d = number - a;
    return std::make_tuple(a, d);
}

TEST_CASE("ParseDecimalNumber") {
    REQUIRE(ParseDecimalNumber("3.14") == std::make_tuple(3, FloatType("0.14")));
    REQUIRE(ParseDecimalNumber(1.5L) == std::make_tuple(1, 0.5L));
    REQUIRE(ParseDecimalNumber("125") == std::make_tuple(125, 0.0));
    REQUIRE_THROWS(ParseDecimalNumber("NaN"));
}

std::vector<DecimalNumber> ContinuedFractionIndices(std::vector<DecimalNumber> &indices, size_t total_iterations) {
    if (indices.size() == 0) {
        throw std::invalid_argument("Indices must not be empty");
    }

    for (int i = indices.size(); i < total_iterations; i++) {
        auto[a, d] = indices.back();
        indices.push_back(ParseDecimalNumber(FloatType(1) / d));
    }

    return indices;
}

std::vector<DecimalNumber> ContinuedFractionIndices(DecimalNumber number, size_t iterations) {
    std::vector<DecimalNumber> indices;
    indices.reserve(iterations);

    indices.push_back(number);

    return ContinuedFractionIndices(indices, iterations);
}

TEST_CASE("ContinuedFractionIndices") {
    const auto kExampleNumber = boost::math::constants::pi<FloatType>();
    const std::vector<long> kExpectedResults = {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84,
                                                2, 1, 1, 15, 3, 13, 1, 4, 2, 6, 6, 99, 1, 2, 2, 6, 3, 5, 1, 1, 6, 8, 1,
                                                7, 1, 2, 3, 7, 1, 2, 1, 1, 12, 1, 1, 1, 3, 1, 1, 8, 1, 1, 2, 1, 6, 1, 1,
                                                5, 2, 2, 3, 1, 2, 4, 4, 16, 1, 161, 45, 1, 22, 1, 2, 2, 1};
    auto indices = ContinuedFractionIndices(ParseDecimalNumber(kExampleNumber), kExpectedResults.size());
    REQUIRE(std::equal(kExpectedResults.begin(), kExpectedResults.end(), indices.begin(), [](auto lhs, auto rhs) {
        INFO("Checking number " << lhs << " and " << rhs);
        auto comparison = lhs == std::get<0>(rhs);
        CHECK(comparison);
        return comparison;
    }));
}

Fraction ContinuedFractionEvaluator(std::vector<DecimalNumber> indices) {
    return std::accumulate(std::next(indices.rbegin()),
                           indices.rend(),
                           Fraction(std::get<0>(indices.back()), 1),
                           [](Fraction accumulator, DecimalNumber number) {
                               return 1 / accumulator + std::get<0>(number);
                           });
}

TEST_CASE("ContinuedFractionEvaluator") {
    const auto kExampleNumber = boost::math::constants::pi<FloatType>();
    const std::vector<Fraction> kExpectedResults = {{3,      1},
                                                    {22,     7},
                                                    {333,    106},
                                                    {355,    133},
                                                    {103993, 33102}};
    constexpr auto kDegree = 10;

    std::vector<Fraction> results(kDegree);
    std::generate(results.begin(), results.end(), [n = 1, kExampleNumber]() mutable {
        auto indices = ContinuedFractionIndices(ParseDecimalNumber(kExampleNumber), n++);
        return ContinuedFractionEvaluator(indices);
    });

    std::cout << results;
}

struct EvaluatedIteration {
    std::vector<DecimalNumber> indices;
    Fraction fraction;
    FloatType evaluated_fraction;
};

void PrintEvaluationTable(const std::vector<EvaluatedIteration> &iterations,
                          size_t precision = std::numeric_limits<FloatType>::digits10) {
    fort::char_table table;
    /* Change border style */
    table.set_border_style(FT_DOUBLE2_STYLE);

    table << fort::header
          << "Iteration" << "Indices" << "Fraction" << "Evaluated fraction" << fort::endr;

    int i = 1;
    for (const auto &iteration : iterations) {
        table << i++;

        std::ostringstream oss;
        oss << '[';
        for (int j = 0; j < iteration.indices.size(); j++) {
            oss << std::get<0>(iteration.indices[j]);
            if (j < iteration.indices.size() - 1)
                oss << ',';
        }
        oss << ']';
        table << oss.str();
        table << iteration.fraction;

        oss.str("");
        oss.clear();
        oss.precision(precision);
        oss << iteration.evaluated_fraction;

        table << oss.str();
        table << fort::endr;
    }
    std::cout << table.to_string();
}

std::vector<EvaluatedIteration> EvaluateDecimalNumber(DecimalNumber number, int iterations) {
    std::vector<DecimalNumber> indices = ContinuedFractionIndices(number, 1);
    std::vector<EvaluatedIteration> evaluated_iterations;
    evaluated_iterations.push_back(
            {indices, ContinuedFractionEvaluator(indices), ContinuedFractionEvaluator(indices).numerator()}
    );

    for (int i = 2; i <= iterations; i++) {
        ContinuedFractionIndices(indices, i);
        auto fraction = ContinuedFractionEvaluator(indices);
        FloatType evaluated_fraction = FloatType(fraction.numerator()) / FloatType(fraction.denominator());
        evaluated_iterations.push_back({indices, fraction, evaluated_fraction});
    }

    return evaluated_iterations;
}

TEST_CASE("EvaluateDecimalNumber") {
    const auto kExampleNumber = boost::math::constants::pi<FloatType>();
    constexpr auto kIterations = 20;
    EvaluateDecimalNumber(ParseDecimalNumber(kExampleNumber), kIterations);
}

int main(int argc, const char *argv[]) {
    namespace po = boost::program_options;

    std::string number;
    size_t iterations;
    size_t precision;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("iterations", po::value<size_t>(&iterations)->default_value(20), "number of iterations (default 20)")
            ("p,print", "print evaluation table")
            ("precision", po::value<size_t>(&precision)->default_value(std::numeric_limits<FloatType>::digits10),
             "how many decimals should result have")
            ("t,test", "run built-in tests")
            ("number", po::value<std::string>(&number)->required(), "number to be parsed");

    po::positional_options_description pos;
    pos.add("number", 1);

    auto print_help = [argv, &desc]() {
        std::cout << "Usage: ." << GetProgramName(argv[0]) << "[OPTIONS] <NUMBER>\n";
        std::cout << desc;
    };

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
        po::notify(vm);
    } catch (const po::error &e) {
        print_help();
        return 1;
    }

    if (vm.count("test")) {
        Catch::Session().run(argc, argv);
        return 0;
    }

    auto evaluated_iterations = EvaluateDecimalNumber(ParseDecimalNumber(number), iterations);
    if (vm.count("p")) {
        PrintEvaluationTable(evaluated_iterations, precision);
    } else {
        std::cout.precision(precision);
        std::cout << evaluated_iterations.back().evaluated_fraction << std::endl;
    }


}