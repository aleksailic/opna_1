/* MIT License
 *
 * Copyright (c) 2020 Aleksa Ilic <aleksa.d.ilic@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#define CATCH_CONFIG_RUNNER

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <regex>

#include <boost/program_options.hpp>
#include <boost/rational.hpp>
#include <boost/optional/optional.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/constants/constants.hpp>

#ifdef TESTING_ENABLED
#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <prettyprint.hpp>
#endif

#include <fort.hpp>

static constexpr auto kVersion = "v0.1.0";
static constexpr auto kProgramName = "OPNA_1: Continual Fraction Generator";
static constexpr auto kUnderlineType = '=';

// -- CHANGE THESE TYPES FOR MORE PRECISION
namespace mp = boost::multiprecision;
using IntType = mp::int1024_t;
using FloatType = mp::number<mp::cpp_dec_float<256>>;
using Fraction = boost::rational<IntType>;
using DecimalNumber = std::tuple<IntType, FloatType>;

struct Config {
    size_t precision;
    size_t iterations;
};
static constexpr Config kDefaultConfig{14, 14};

std::ostream &operator<<(std::ostream &os, const DecimalNumber &number) {
    return os << '(' << std::get<0>(number) << ',' << std::get<1>(number) << ')';
}

boost::optional<std::string> SafeSubstring(const std::string &string, size_t pos, size_t len = std::string::npos) {
    if (pos >= string.length())
        return boost::none;
    else
        return string.substr(pos, len);
}

DecimalNumber ParseDecimalNumber(FloatType number) {
    auto a = mp::floor(number);
    auto d = number - a;
    return std::make_tuple(a.convert_to<IntType>(), d);
}

DecimalNumber ParseDecimalNumber(std::string number) {
    if (number == "pi") {
        return ParseDecimalNumber(boost::math::constants::pi<FloatType>());
    } else if (number == "phi") {
        return ParseDecimalNumber(boost::math::constants::phi<FloatType>());
    } else if (number == "e") {
        return ParseDecimalNumber(boost::math::constants::e<FloatType>());
    }

    auto dot_index = number.find('.');
    if (dot_index == std::string::npos) {
        return std::make_tuple(std::stol(number), 0);
    }
    std::string a = number.substr(0, dot_index);
    std::string d = "0." + SafeSubstring(number, dot_index + 1).get_value_or("0");

    return std::make_tuple(std::stoll(a), FloatType(d));
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

Fraction ContinuedFractionEvaluator(std::vector<DecimalNumber> indices) {
    return std::accumulate(std::next(indices.rbegin()),
                           indices.rend(),
                           Fraction(std::get<0>(indices.back()), 1),
                           [](Fraction accumulator, DecimalNumber number) {
                               return 1 / accumulator + std::get<0>(number);
                           });
}

enum class EvaluatedType {
    I, II
};
struct EvaluatedIteration {
    std::vector<DecimalNumber> indices;
    Fraction fraction;
    FloatType evaluated_fraction;
    FloatType diff;
    EvaluatedType type;
};

EvaluatedIteration
EvaluateIteration(FloatType number_searched, const std::vector<DecimalNumber> &indices, EvaluatedType type) {
    Fraction fraction = ContinuedFractionEvaluator(indices);
    if (fraction < 0) {
        throw std::overflow_error(
                "Underlying type overflow. Compile with larger integer and/or floating point types.");
    }
    auto evaluated_fraction = FloatType(fraction.numerator()) / FloatType(fraction.denominator());
    auto diff = number_searched - evaluated_fraction;
    return EvaluatedIteration{indices, fraction, evaluated_fraction, diff, type};
}

void PrintEvaluatedIteration(fort::char_table &os, const EvaluatedIteration &iteration, size_t precision) {
    std::ostringstream oss;
    oss << '[';
    for (int j = 0; j < iteration.indices.size(); j++) {
        oss << std::get<0>(iteration.indices[j]);
        if (j < iteration.indices.size() - 1)
            oss << ',';
    }
    oss << ']';
    os << oss.str();

    oss.str("");
    oss.clear();
    oss << iteration.fraction;
    if (iteration.type == EvaluatedType::I)
        oss << '*';
    os << oss.str();

    os << iteration.evaluated_fraction.str(precision, std::ios::fixed);

    os << std::regex_replace(iteration.diff.str(precision, std::ios::scientific), std::regex("e"), " * 10^");
}

void PrintEvaluationTable(const std::vector<EvaluatedIteration> &iterations,
                          size_t precision = std::numeric_limits<FloatType>::digits10) {
    fort::char_table table;
    table.set_border_style(FT_NICE_STYLE);

    table << fort::header
          << "Iteration" << "Indices" << "Fraction" << "Evaluated fraction" << "Difference" << fort::endr;

    int i = 1;
    for (const auto &iteration : iterations) {
        table << i++;
        PrintEvaluatedIteration(table, iteration, precision);
        table << fort::endr;
        table[i - 1][4].set_cell_text_align(fort::text_align::right);
    }

    std::cout << table.to_string();
}

std::vector<EvaluatedIteration>
FindInBetweenApproximations(FloatType number_searched, std::vector<DecimalNumber> indices, FloatType diff) {
    auto from = IntType(1);
    auto to = std::get<0>(indices.back());
    indices.pop_back();

    std::vector<EvaluatedIteration> iterations;
    for (auto i = from; i < to; i++) {
        indices.push_back(std::make_tuple(i, 0));
        auto candidate = EvaluateIteration(number_searched, indices, EvaluatedType::II);
        if (mp::abs(candidate.diff) < mp::abs(diff)) {
            iterations.push_back(candidate);
        }
        indices.pop_back();
    }
    return iterations;
}

std::vector<EvaluatedIteration>
EvaluateDecimalNumber(DecimalNumber number, size_t iterations, bool find_in_between = false) {
    const auto number_searched = FloatType(std::get<0>(number)) + FloatType(std::get<1>(number));
    std::vector<DecimalNumber> indices = ContinuedFractionIndices(number, 1);
    std::vector<EvaluatedIteration> evaluated_iterations;

    auto fraction = ContinuedFractionEvaluator(indices);
    auto evaluated_fraction = FloatType(fraction.numerator()) / FloatType(fraction.denominator());
    evaluated_iterations.push_back(
            {indices, fraction, evaluated_fraction, number_searched - evaluated_fraction, EvaluatedType::I});

    try {
        for (int i = 2; i <= iterations; i++) {
            ContinuedFractionIndices(indices, i);
            auto evaluated_iteration = EvaluateIteration(number_searched, indices, EvaluatedType::I);
            if (mp::abs(evaluated_iteration.diff) > mp::abs(evaluated_iterations.back().diff)) {
                throw std::logic_error(
                        "More iterations result in bigger deviation. Compile with larger integer and/or floating point types.");
            }
            if (find_in_between) {
                auto in_between_aproximations = FindInBetweenApproximations(number_searched, indices,
                                                                            evaluated_iterations.back().diff);
                evaluated_iterations.insert(evaluated_iterations.end(), in_between_aproximations.begin(),
                                            in_between_aproximations.end());
            }
            evaluated_iterations.push_back(evaluated_iteration);
        }
    } catch (const std::range_error &err) {
        std::cerr << err.what() << " Stopped at iteration: " << indices.size() << std::endl;
    } catch (const std::overflow_error &err) {
        std::cerr << err.what() << " Stopped at iteration: " << indices.size() << std::endl;
    } catch (const std::logic_error &err) {
        std::cerr << err.what() << " Stopped at iteration: " << indices.size() << std::endl;
    }

    return evaluated_iterations;
}

static const char *GetProgramName(const char *path) {
    const char *last = path;
    while (*path++) {
        last = (*path == '/' || *path == '\\') ? path : last;
    }
    return last;
}

#ifdef TESTING_ENABLED
TEST_CASE("GetProgramName") {
    constexpr auto kExamplePath1 = "/home/ilic/opna-1/opna_1";
    constexpr auto kExamplePath2 = "C:\\Users\\My Documents\\opna_1";
    constexpr auto kExpectedName1 = "/opna_1";
    constexpr auto kExpectedName2 = "\\opna_1";

    REQUIRE(strcmp(GetProgramName(kExamplePath1), kExpectedName1) == 0);
    REQUIRE(strcmp(GetProgramName(kExamplePath2), kExpectedName2) == 0);
}
TEST_CASE("ContinuedFractionEvaluator") {
    const auto kExampleNumber = boost::math::constants::pi<FloatType>();
    const std::vector<Fraction> kExpectedResults = {{3,        1},
                                                    {22,       7},
                                                    {333,      106},
                                                    {355,      113},
                                                    {103993,   33102},
                                                    {104348,   33215},
                                                    {208341,   66317},
                                                    {312689,   99532},
                                                    {833719,   265381},
                                                    {1146408,  364913},
                                                    {4272943,  1360120},
                                                    {5419351,  1725033},
                                                    {80143857, 25510582}};
    constexpr auto kDegree = 13;
    std::vector<Fraction> results(kDegree);
    std::generate(results.begin(), results.end(), [n = 1, kExampleNumber]() mutable {
        auto indices = ContinuedFractionIndices(ParseDecimalNumber(kExampleNumber), n++);
        return ContinuedFractionEvaluator(indices);
    });

    REQUIRE(std::equal(kExpectedResults.begin(), kExpectedResults.end(), results.begin(), [](auto lhs, auto rhs) {
        INFO("Checking fraction " << lhs << " and " << rhs);
        CHECK(lhs == rhs);
        return lhs == rhs;
    }));
}
TEST_CASE("ContinuedFractionIndices") {
    const auto kExampleNumber = boost::math::constants::pi<FloatType>();
    const std::vector<IntType> kExpectedResults = {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1,
                                                   84,
                                                   2, 1, 1, 15, 3, 13, 1, 4, 2, 6, 6, 99, 1, 2, 2, 6, 3, 5, 1, 1, 6, 8,
                                                   1,
                                                   7, 1, 2, 3, 7, 1, 2, 1, 1, 12, 1, 1, 1, 3, 1, 1, 8, 1, 1, 2, 1, 6, 1,
                                                   1,
                                                   5, 2, 2, 3, 1, 2, 4, 4, 16, 1, 161, 45, 1, 22, 1, 2, 2, 1};
    auto indices = ContinuedFractionIndices(ParseDecimalNumber(kExampleNumber), kExpectedResults.size());
    REQUIRE(std::equal(kExpectedResults.begin(), kExpectedResults.end(), indices.begin(), [](auto lhs, auto rhs) {
        INFO("Checking number " << lhs << " and " << rhs);
        auto comparison = lhs == std::get<0>(rhs);
        CHECK(comparison);
        return comparison;
    }));
}
TEST_CASE("ParseDecimalNumber") {
    REQUIRE(ParseDecimalNumber("3.14") == std::make_tuple(3, FloatType("0.14")));
    REQUIRE(ParseDecimalNumber(1.5L) == std::make_tuple(1, 0.5L));
    REQUIRE(ParseDecimalNumber("125") == std::make_tuple(125, 0.0));
    REQUIRE_THROWS(ParseDecimalNumber("NaN"));
}
#else
int main(int argc, const char *argv[]) {
    namespace po = boost::program_options;

    Config config = kDefaultConfig;

    auto program_name = GetProgramName(argv[0]);
    std::string number;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "print help message")
            ("version", "print version information")
            ("examples", "show examples")
            ("table", "print evaluation table of every iteration")
            ("inbetween", "show best in-between continual fraction approximations also")
            ("iterations", po::value<size_t>(&config.iterations)->default_value(kDefaultConfig.iterations),
             "number of iterations")
            ("precision", po::value<size_t>(&config.precision)->default_value(kDefaultConfig.precision),
             "how many decimals should result have")
            ("number", po::value<std::string>(&number), "number to be parsed");;

    po::positional_options_description pos_desc;
    pos_desc.add("number", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "Usage: ." << program_name << " [OPTIONS] <NUMBER|CONSTANT> \n"
                  << desc << '\n'
                  << "Allowed constants: pi, phi, e \n";
        return 0;
    }

    if (vm.count("version")) {
        std::ostringstream name_header;
        name_header << kProgramName << ' ' << kVersion;

        std::string underline(name_header.str().size(), kUnderlineType);

        std::cout << name_header.str() << "\n"
                  << underline << "\n"
                  << "Integer type w/: " << std::numeric_limits<IntType>::digits << "bits\n"
                  << "Floating type w/: " << std::numeric_limits<FloatType>::digits << "bits\n";
        return 0;
    }

    if (vm.count("examples")) {
        fort::char_table table;
        table.set_border_style(FT_EMPTY_STYLE);

        const std::vector<std::tuple<std::string, std::string>> examples = {
                std::make_tuple("pi",
                                "Evaluates pi with default (" + std::to_string(config.iterations) + ") iterations"),
                std::make_tuple("pi --table --iterations 5 --precision 10",
                                "Print evaluation table with 5 iterations and up to 10 decimal places for pi expansion")
        };

        table << fort::header
              << "" << "Command" << "Description" << fort::endr;

        for (size_t i = 0; i < examples.size(); i++) {
            const auto&[command, description] = examples[i];
            table << i + 1 << std::string(".") + program_name + " " + command << description << fort::endr;
        }

        std::cout << table.to_string();
        return 0;
    }

    if (!vm.count("number")) {
        std::cerr << "You need to specify a number. Use --help for usage information \n";
        return 1;
    }

    auto evaluated_iterations = EvaluateDecimalNumber(ParseDecimalNumber(number), config.iterations, vm.count("inbetween"));
    if (vm.count("table")) {
        PrintEvaluationTable(evaluated_iterations, config.precision);
    } else {
        const auto &eval = evaluated_iterations.back();
        std::ostringstream oss;
        oss.precision(config.precision);
        oss << eval.evaluated_fraction;

        fort::char_table table;
        table.set_border_style(FT_EMPTY_STYLE);
        table << fort::header
              << "Iterations" << "Indices" << "Fraction" << "Evaluated Fraction" << "Difference" << fort::endr;
        table << evaluated_iterations.size();
        PrintEvaluatedIteration(table, eval, config.precision);
        table << fort::endr;

        std::cout << table.to_string();
    }
}
#endif