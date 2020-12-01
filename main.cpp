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

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <tuple>
#include <regex>

#include <boost/algorithm/string.hpp>
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
#else
#endif

#include <fort.hpp>

static constexpr auto kVersion = "v0.2.0";
static constexpr auto kProgramName = "OPNA_1: Continued Fraction Generator";
static constexpr auto kUnderlineType = '=';

// -- CHANGE THESE TYPES FOR MORE PRECISION
namespace mp = boost::multiprecision;
using IntType = mp::int1024_t;
using FloatType = mp::number<mp::cpp_dec_float<256>>;
using Fraction = boost::rational<IntType>;
using DecimalNumber = std::tuple<IntType, FloatType>;

enum class Field {
    ITERATION, INDICES, FRACTION, EVALUATED_FRACTION, DIFFERENCE
};

constexpr const char *FieldToString(Field field) {
    switch (field) {
        case Field::ITERATION:
            return "Iteration";
        case Field::INDICES:
            return "Indices";
        case Field::FRACTION:
            return "Fraction";
        case Field::EVALUATED_FRACTION:
            return "Evaluated fraction";
        case Field::DIFFERENCE:
            return "Difference";
    }
}

struct Config {
    size_t precision = 14;
    size_t iterations = 14;
    IntType max_denominator = std::numeric_limits<IntType>::max();
    bool find_in_between = false;
    bool headerless = false;
    const struct ft_border_style *table_style = FT_NICE_STYLE;
    std::vector<Field> displayed_fields = {Field::ITERATION, Field::INDICES,
                                           Field::FRACTION, Field::EVALUATED_FRACTION,
                                           Field::DIFFERENCE};
};
static const Config kDefaultConfig;

/// Pretty prints DecimalNumber to stream
std::ostream &operator<<(std::ostream &os, const DecimalNumber &number) {
    return os << '(' << std::get<0>(number) << ',' << std::get<1>(number) << ')';
}

/// Performs substring operation with range check without throwing errors
/// \param string String to be processed
/// \param position Position from which to extract
/// \param length How many characters are to be processed
/// \return Extracted substring if range check passes or boost::none
boost::optional<std::string>
SafeSubstring(const std::string &string, size_t position, size_t length = std::string::npos) {
    if (position >= string.length())
        return boost::none;
    else
        return string.substr(position, length);
}

/// Splits decimal number into whole and fraction part
/// \param numer Number to be parsed
/// \return Tuple containing whole and fraction
DecimalNumber ParseDecimalNumber(FloatType number) {
    auto a = mp::floor(number);
    auto d = number - a;
    return std::make_tuple(a.convert_to<IntType>(), d);
}

/// Split decimal number into whole and fraction part
/// \param number Number to be parsed
/// \return Tuple containing whole and fraction
DecimalNumber ParseDecimalNumber(std::string number) {
    if (number == "pi") {
        return ParseDecimalNumber(boost::math::constants::pi<FloatType>());
    } else if (number == "phi") {
        return ParseDecimalNumber(boost::math::constants::phi<FloatType>());
    } else if (number == "e") {
        return ParseDecimalNumber(boost::math::constants::e<FloatType>());
    } else if (number == "catalan") {
        return ParseDecimalNumber(boost::math::constants::catalan<FloatType>());
    }

    auto dot_index = number.find('.');
    if (dot_index == std::string::npos) {
        return std::make_tuple(std::stol(number), 0);
    }
    std::string a = number.substr(0, dot_index);
    std::string d = "0." + SafeSubstring(number, dot_index + 1).get_value_or("0");

    return std::make_tuple(std::stoll(a), FloatType(d));
}

/// Generates continued fraction indices
/// \param indices Generated indices from previous iterations
/// \param total_iterations Total number of iterations to run (max indices size)
/// \return Generated indices
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

/// Generates continued fraction indices
/// \param number Decimal number for which to calculate indices
/// \param iterations Number of iterations
/// \return Generated indices
std::vector<DecimalNumber> ContinuedFractionIndices(DecimalNumber number, size_t iterations) {
    std::vector<DecimalNumber> indices;
    indices.reserve(iterations);

    indices.push_back(number);
    return ContinuedFractionIndices(indices, iterations);
}

/// Evaluate continued fraction indices to get fractional representation
/// \param indices Continued fraction indices
/// \return Evaluated fraction
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

/// Evaluates single iteration from indices and populates all relevant fields
/// \param number_searched Number for which we are deducing optimal indices
/// \param indices Generated continued fraction indices
/// \param type Type of iteration (First order/Second order)
/// \return
EvaluatedIteration
EvaluateIteration(FloatType number_searched, const std::vector<DecimalNumber> &indices,
                  EvaluatedType type) {
    Fraction fraction = ContinuedFractionEvaluator(indices);
    if (fraction < 0) {
        throw std::overflow_error(
                "Underlying type overflow. Compile with larger integer and/or floating point types.");
    }
    auto evaluated_fraction = FloatType(fraction.numerator()) / FloatType(fraction.denominator());
    auto diff = number_searched - evaluated_fraction;
    return EvaluatedIteration{indices, fraction, evaluated_fraction, diff, type};
}

void PrintEvaluatedIteration(fort::char_table &table, const EvaluatedIteration &iteration, const Config &config) {
    for (const auto &header : config.displayed_fields) {
        std::ostringstream oss;
        switch (header) {
            case Field::ITERATION:
                table << iteration.indices.size();
                break;
            case Field::INDICES:
                oss << '[';
                for (int j = 0; j < iteration.indices.size(); j++) {
                    oss << std::get<0>(iteration.indices[j]);
                    if (j < iteration.indices.size() - 1)
                        oss << ',';
                }
                oss << ']';
                table << oss.str();
                break;
            case Field::FRACTION:
                oss << iteration.fraction;
                if (iteration.type == EvaluatedType::I)
                    oss << '*';
                table << oss.str();
                break;
            case Field::EVALUATED_FRACTION:
                table << iteration.evaluated_fraction.str(config.precision, std::ios::fixed);
                break;
            case Field::DIFFERENCE:
                table << std::regex_replace(iteration.diff.str(config.precision, std::ios::scientific), std::regex("e"),
                                            " * 10^");
                break;
        }
    }
    table << fort::endr;
}

void PrintEvaluationTableHeader(fort::char_table &table, const Config &config) {
    if (!config.headerless) {
        table << fort::header;
        for (const auto &field : config.displayed_fields) {
            table << FieldToString(field);
        }
        table << fort::endr;
    }
}

void PrintEvaluationTable(const std::vector<EvaluatedIteration> &iterations, const Config &config) {
    fort::char_table table;
    table.set_border_style(config.table_style);

    auto diff_it = std::find(config.displayed_fields.begin(), config.displayed_fields.end(), Field::DIFFERENCE);
    PrintEvaluationTableHeader(table, config);

    size_t counter = 1;
    for (const auto &iteration : iterations) {
        PrintEvaluatedIteration(table, iteration, config);
        if(diff_it != config.displayed_fields.end()){
            table[counter++][std::distance(config.displayed_fields.begin(), diff_it)].set_cell_text_align(fort::text_align::right);
        }
    }

    std::cout << table.to_string();
}

/// Find indices of in-between approximations (Indices of second order)
/// \param number_searched Number for which we are deducing second order indices
/// \param indices Calculated continued fraction indices of the first order
/// \param diff Evaluated difference of indices of the first order
/// \return Vector of evaluated iterations of approximations of the second order
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

/// Processes decimal number with given configuration. Throws overflow error.
/// \param number Decimal number to be processed
/// \param config Main program configuration
/// \return Vector of evaluated iterations
std::vector<EvaluatedIteration>
ProcessDecimalNumber(DecimalNumber number, const Config &config) {
    const auto number_searched = FloatType(std::get<0>(number)) + FloatType(std::get<1>(number));
    std::vector<DecimalNumber> indices = ContinuedFractionIndices(number, 1);
    std::vector<EvaluatedIteration> evaluated_iterations;

    auto fraction = ContinuedFractionEvaluator(indices);
    auto evaluated_fraction = FloatType(fraction.numerator()) / FloatType(fraction.denominator());
    evaluated_iterations.push_back(
            {indices, fraction, evaluated_fraction, number_searched - evaluated_fraction, EvaluatedType::I});

    try {
        for (size_t iteration_number = 2; iteration_number <= config.iterations; iteration_number++) {
            ContinuedFractionIndices(indices, iteration_number);
            auto evaluated_iteration = EvaluateIteration(number_searched, indices, EvaluatedType::I);
            if (mp::abs(evaluated_iteration.diff) > mp::abs(evaluated_iterations.back().diff)) {
                throw std::overflow_error(
                        "More iterations result in bigger deviation. Compile with larger integer and/or floating point types.");
            }
            if (config.find_in_between) {
                auto in_between_aproximations = FindInBetweenApproximations(number_searched, indices,
                                                                            evaluated_iterations.back().diff);
                evaluated_iterations.insert(evaluated_iterations.end(), in_between_aproximations.begin(),
                                            in_between_aproximations.end());
            }
            evaluated_iterations.push_back(evaluated_iteration);
            if (evaluated_iterations.back().fraction.denominator() > config.max_denominator) {
                evaluated_iterations.erase(std::lower_bound(evaluated_iterations.begin(), evaluated_iterations.end(),
                                                            config.max_denominator,
                                                            [](const auto &iteration, const auto &max_denominator) {
                                                                return iteration.fraction.denominator() <=
                                                                       max_denominator;
                                                            }), evaluated_iterations.end());
                break;
            }
        }
    } catch (const std::overflow_error &err) {
        std::cerr << err.what() << " Stopped at iteration: " << indices.size() << std::endl;
    }

    std::sort(evaluated_iterations.begin(), evaluated_iterations.end(),
              [](auto &lhs, auto &rhs) {
                  return mp::abs(lhs.diff) >= mp::abs(rhs.diff);
              });

    return evaluated_iterations;
}

/// Retrieves program's base name from the exec path
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
    std::string table_style;
    std::string fields;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "print help message")
            ("version", "print version information")
            ("examples", "show examples")
            ("table", "print evaluation table of every iteration")
            ("table_style", po::value<std::string>(&table_style), "Specify table style: nice|double|simple|empty")
            ("headerless", "Do not display header when showing results")
            ("fields", po::value<std::string>(&fields),
             "Comma delimited list of fields to be shown: iter,ind,frac,eval,diff")
            ("maxdenominator", po::value<IntType>(&config.max_denominator),
             "maximum denominator up to which to iterate")
            ("inbetween", "show best in-between continual fraction approximations as well")
            ("iterations", po::value<size_t>(&config.iterations)->default_value(kDefaultConfig.iterations),
             "number of iterations")
            ("precision", po::value<size_t>(&config.precision)->default_value(kDefaultConfig.precision),
             "how many decimals should result have")
            ("number", po::value<std::string>(&number), "number to be parsed");

    po::positional_options_description pos_desc;
    pos_desc.add("number", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "Usage: ." << program_name << " [OPTIONS] <NUMBER|CONSTANT> \n"
                  << desc << '\n'
                  << "Allowed constants: pi, phi, e, catalan \n";
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
                                "Print evaluation table with 5 iterations and up to 10 decimal places for pi expansion"),
                std::make_tuple("phi --table --maxdenominator 400 --inbetween",
                                "Print evaluation table for phi where maximum fraction approximation denominator is less than or equal to 400"),
                std::make_tuple("phi --table --inbetween",
                                "Print evaluation table for phi with default (" + std::to_string(config.iterations) +
                                ") iterations and also find best in-between approximations")
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

    if (vm.count("inbetween")) {
        config.find_in_between = true;
    }

    if (vm.count("headerless")) {
        config.headerless = true;
    }

    if (vm.count("table_style")) {
        if (boost::iequals(table_style, "nice")) {
            config.table_style = FT_NICE_STYLE;
        } else if (boost::iequals(table_style, "simple")) {
            config.table_style = FT_SIMPLE_STYLE;
        } else if (boost::iequals(table_style, "double")) {
            config.table_style = FT_DOUBLE2_STYLE;
        } else if (boost::iequals(table_style, "empty")) {
            config.table_style = FT_EMPTY_STYLE;
        } else {
            std::cerr << "Incorrect table style supplied: " << table_style << ". Use --help for usage information \n";
            return 1;
        }
    }

    if (vm.count("fields")) {
        config.displayed_fields.clear();
        std::vector<std::string> field_vector;
        boost::split(field_vector, fields, [](char c) { return c == ','; });
        for (const auto &field:field_vector) {
            if (boost::iequals(field, "iter")) {
                config.displayed_fields.push_back(Field::ITERATION);
            } else if (boost::iequals(field, "ind")) {
                config.displayed_fields.push_back(Field::INDICES);
            } else if (boost::iequals(field, "frac")) {
                config.displayed_fields.push_back(Field::FRACTION);
            } else if (boost::iequals(field, "eval")) {
                config.displayed_fields.push_back(Field::EVALUATED_FRACTION);
            } else if (boost::iequals(field, "diff")) {
                config.displayed_fields.push_back(Field::DIFFERENCE);
            } else {
                std::cerr << "Incorrect field name supplied: " << field << ". Use --help for usage information \n";
                return 1;
            }
        }
    }

    auto evaluated_iterations = ProcessDecimalNumber(ParseDecimalNumber(number), config);
    if (vm.count("table")) {
        PrintEvaluationTable(evaluated_iterations, config);
    } else {
        const auto &eval = evaluated_iterations.back();

        fort::char_table table;
        table.set_border_style(config.table_style);
        PrintEvaluationTableHeader(table, config);
        PrintEvaluatedIteration(table, eval, config);

        std::cout << table.to_string();
    }
}

#endif