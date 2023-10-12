#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include <random>
#include <iterator>
#include <cmath>

#include "../include/genetics.hpp"

typedef float domain_t;
typedef std::vector<domain_t> polynomial_t;

domain_t interpret(polynomial_t polynomial, domain_t variable) {
	if (polynomial.size() == 0) return static_cast<domain_t>(0);
	domain_t result = polynomial[polynomial.size()-1];
	for (size_t i = 1; i < polynomial.size(); ++i) {
		result = result * variable + polynomial[polynomial.size()-i-1];
	}
	return result;
}

std::string listing(polynomial_t polynomial, std::string varname = "x") {
	std::string result = "";
	for (size_t i = 0; i < polynomial.size(); ++i) {
		//if (polynomial[polynomial.size()-i-1] == 0) continue;
		if (result == "") {
			result = std::to_string(polynomial[polynomial.size()-i-1]) + (i < polynomial.size()-1 ? varname : "") + (polynomial.size() >= 3 && i < polynomial.size() - 2 ? std::string("^") + std::to_string(polynomial.size()-i-1) : "");
		} else {
			if (polynomial[polynomial.size()-i-1] >= 0) {
				result += " + " + std::to_string(polynomial[polynomial.size()-i-1]);
			} else {
				result += " - " + std::to_string(-polynomial[polynomial.size()-i-1]);
			}
			result += (i < polynomial.size()-1 ? varname : "") + (i < polynomial.size() - 2 ? std::string("^") + std::to_string(polynomial.size()-i-1) : "");
		}
	}
	if (result == "") return "0";
	else return result.substr(0, result.size()-3);
}

double difference_with(polynomial_t const& one, polynomial_t const& another) {
	double result = 0, monomial_difference;
	size_t min_size = std::min(one.size(), another.size()),
		   max_size = std::max(one.size(), another.size());

	for (size_t i = 0; i < min_size; ++i) {
		// MAGIC!
		/*
		monomial_difference = -(std::min(one[i], another[i]) - std::max(one[i], another[i]))
							  / std::max(std::max(one[i], another[i]), -std::min(one[i], another[i]));
		*/

		// Less magic
		domain_t coefficient_distance = std::max(one[i], another[i]) - std::min(one[i], another[i]);
		domain_t maximum_absolute_coefficient = std::max(std::abs(one[i]), std::abs(another[i]));

		monomial_difference = coefficient_distance / maximum_absolute_coefficient; // With normalization
		//monomial_difference = coefficient_distance; // Without normalization

		result += monomial_difference * i; // Greater power - greater importance
		//result += monomial_difference; // Plain
	}

	// MORE MAGIC!!!
	// Why arithmetic progression? It's like monomial_difference for all other coefficients is 1.
	result += static_cast<double>(min_size + max_size - 1) * (max_size - min_size) / 2;
	//result += max_size - min_size;

	result /= (static_cast<double>(max_size - 1) * max_size / 2); // Normilize result (with greatest progression member)
	//result /= max_size; // Normilize result (with polynomial degree)

	return result;
};

class PolynomialCost {
	private:
		domain_t inaccuracy;
		size_t size;

		friend std::ostream& operator<<(std::ostream& os, PolynomialCost const& c);

	public:
		friend void swap(PolynomialCost& one, PolynomialCost& another) {
			using std::swap;
			swap(one.inaccuracy, another.inaccuracy);
			swap(one.size, another.size);
		}
		PolynomialCost() : inaccuracy(.0), size(0) { }
		PolynomialCost(domain_t inaccuracy, size_t size) : inaccuracy(inaccuracy), size(size) { }
		PolynomialCost(PolynomialCost const& another) : inaccuracy(another.inaccuracy), size(another.size) { }
		PolynomialCost(PolynomialCost&& another) {
			swap(*this, another);
		}
		PolynomialCost& operator=(PolynomialCost another) {
			swap(*this, another);
			return *this;
		}

		bool operator==(PolynomialCost const& another) const {
			return (inaccuracy == another.inaccuracy && size == another.size);
		}
		bool operator!=(PolynomialCost const& another) const {
			return (!(*this == another));
		}

		bool operator<(PolynomialCost const& another) const {
			return ((inaccuracy < another.inaccuracy) || (inaccuracy == another.inaccuracy && size < another.size));
		}
		bool operator>=(PolynomialCost const& another) const {
			return (!(*this < another));
		}
		bool operator>(PolynomialCost const& another) const {
			return (*this >= another && *this != another);
		}
		bool operator<=(PolynomialCost const& another) const {
			return (!(*this > another));
		}

		size_t get_size() const {
			return size;
		}
		domain_t get_inaccuracy() const {
			return inaccuracy;
		}
};

std::ostream& operator<<(std::ostream& os, PolynomialCost const& c) {
	return (os << "(" << c.inaccuracy << "; " << c.size << ")");
}

class PolyFitness : public genetics::IFitness<polynomial_t, PolynomialCost> {
public:
	PolyFitness(std::vector<std::pair<domain_t, domain_t>> target) : target(target) { }
	~PolyFitness() override = default;

	void cost(genetics::Generation<polynomial_t> const& generation, genetics::GenerationsCosts<PolynomialCost>& costs) override {
		auto& polynomials = std::get<0>(generation);
		auto& ages = std::get<1>(generation);

		domain_t inaccuracy, output;

		for (size_t index = 0; index < polynomials.size(); ++index) {
			inaccuracy = 0.f;
			for (auto&& testcase : target) {
				output = interpret(polynomials[index], testcase.first);
				//temp = (std::max(testcase.second, temp) - std::min(testcase.second, temp));
				inaccuracy += (testcase.second - output) * (testcase.second - output);
			}
			costs.emplace_back(inaccuracy, polynomials[index].size());
		}
	}

private:
	std::vector<std::pair<domain_t, domain_t>> target;
};

class PolyCrossover : public genetics::ICrossover<polynomial_t, PolynomialCost> {
public:
	~PolyCrossover() override = default;

	bool is_symmetric() { return false; } // Is crossover symmetric?
	
	polynomial_t cross(genetics::Generation<polynomial_t> const& generation, genetics::GenerationsCosts<PolynomialCost> const& costs, std::size_t parent1, std::size_t parent2) override {
		auto& polynomials = std::get<0>(generation);
		auto& ages = std::get<1>(generation);


		auto& one = polynomials[parent1];
		auto& another = polynomials[parent2];

		// Generate new poly
		size_t first_joint = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, one.size()),
			   second_joint = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, another.size());
		polynomial_t new_polynomial(first_joint + another.size() - second_joint);
		std::copy_n(std::begin(one), first_joint, std::begin(new_polynomial));
		std::copy_n(std::next(std::begin(another), second_joint), another.size() - second_joint, std::next(std::begin(new_polynomial), first_joint));

		// How similar parents are?
		/*
		domain_t cumulative_sum = 0;
		size_t const min_size = std::min(one.size(), another.size());

		for (size_t i = 0; i < min_size; ++i) {
			cumulative_sum += std::abs(one[i] - another[i]);
		}

		double const similarity = (10.0 * min_size - cumulative_sum) / (10.0 * min_size);
		*/
		double const similarity = difference_with(one, another);

		// Will new poly mutate and how?
		double constexpr mutation_probability = 0.5;

		// Note: yep, with similarity > 0.5, there always will be a mutation
		if (genetics::RandomGenerator::get_instance().get_random_float<double>(0.0, 1.0) > mutation_probability + similarity
			|| new_polynomial.size() == 0) {
			return new_polynomial;
		}

		// There can be 6 types of mutation: changing coefficient value, (inserting/losing) a coefficient, cutting (head/tail) and "knockoff"
		std::size_t choice = genetics::RandomGenerator::get_instance().get_random_int<std::size_t>(0, 5);
		switch (choice) {
		case 0: { // change coefficient value
			domain_t update = genetics::RandomGenerator::get_instance().get_random_float<domain_t>(-3.f, 3.f);
			for (auto& coefficient : new_polynomial) {
				coefficient *= update;
			}
			break;
		}
		case 1: { // insert new coefficients/monomials
			size_t how_much = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, new_polynomial.size() / 4);
			size_t cell;
			for (size_t i = 0; i < how_much; ++i) {
				cell = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, new_polynomial.size() - 1);
				new_polynomial.insert(std::next(std::begin(new_polynomial), cell),
				                      genetics::RandomGenerator::get_instance().get_random_float<domain_t>(-1.f, 1.f));
			}
			break;
		}
		case 2: { // remove/zero coefficients/monimials
			size_t how_much = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, new_polynomial.size() / 4);
			size_t cell;
			for (size_t i = 0; i < how_much; ++i) {
				cell = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, new_polynomial.size()-1);
				new_polynomial[cell] = 0.f;
			}
			break;
		}
		case 3: { // cut poly head
			size_t cell = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, new_polynomial.size() - 1);
			new_polynomial.erase(std::begin(new_polynomial), std::next(std::begin(new_polynomial), cell));
			break;
		}
		case 4: { // cut poly tail
			size_t cell = genetics::RandomGenerator::get_instance().get_random_int<size_t>(0, new_polynomial.size() - 1);
			new_polynomial.erase(std::next(std::begin(new_polynomial), cell + 1), std::end(new_polynomial));
			break;
		}
		case 5: // "knock-off"
			for (auto& coefficient : new_polynomial) {
				coefficient = std::max(-0.25f, std::min(coefficient, 0.25f));
			}
			break;
		}

		return new_polynomial;
	}
};

class PolySelection : public genetics::ISelection<polynomial_t, PolynomialCost> {
public:
	PolySelection(std::size_t gen_survivors, std::size_t max_generations_til_end) : gen_survivors(gen_survivors), max_generations_til_end(max_generations_til_end) { }
	~PolySelection() override = default;

	// Number of specimens selected from each generation.
	std::size_t survivors() override { return gen_survivors; }
	void set_survivors(std::size_t new_val) { gen_survivors = new_val; }

	// Limit of generations for one evolve call
	std::optional<std::size_t> max_generations() override { return max_generations_til_end; }
	void set_max_generations(std::size_t new_val) { max_generations_til_end = new_val; }

	// Selection will be performed only every generations_till_eliminaion generations
	std::size_t generations_till_eliminaion() override { return generations_per_phase; }
	void set_generations_till_elimination(std::size_t new_val) { generations_per_phase = new_val; }

private:
	std::size_t gen_survivors;
	std::size_t max_generations_til_end;
	std::size_t generations_per_phase = 1;
};


int main() {
	std::vector<std::pair<domain_t, domain_t>> target;
	size_t test_pairs, generation_cap, survivors;
	std::pair<domain_t, domain_t> buffer;

	std::cout << "How many test pairs? " << std::flush;
	std::cin >> test_pairs;

	std::cout << "Please, enter tests in format: <source value> <target value>" << std::endl;
	for (size_t i = 0; i < test_pairs; ++i) {
		std::cin >> buffer.first >> buffer.second;
		target.push_back(buffer);
	}

	std::cout << "How many generations this civilization should persist? " << std::flush;
	std::cin >> generation_cap;
	std::cout << "How many survivors may live each generation? " << std::flush;
	std::cin >> survivors;

	std::vector<polynomial_t> initial_polys;
	for (size_t i = 0; i < survivors; ++i) {
		polynomial_t pbuffer;
		std::size_t const guessed_length = genetics::RandomGenerator::get_instance().get_random_int<std::size_t>(0, 10);
		for (size_t j = 0; j < guessed_length; ++j) {
			pbuffer.push_back(genetics::RandomGenerator::get_instance().get_random_float<domain_t>(-10.f, 10.f));
		}
		initial_polys.push_back(std::move(pbuffer));
	}

	auto the_nonglitch = genetics::new_generation(initial_polys);

	auto fitness = std::make_shared<PolyFitness>(target);
	auto crossover = std::make_shared<PolyCrossover>();
	auto selection = std::make_shared<PolySelection>(survivors, generation_cap); // TODO: NOT gencap. I want to see progress

	genetics::Environment<polynomial_t, PolynomialCost> world(fitness, crossover, selection);

	//std::size_t const gens_till_death = 1;
	std::size_t const gens_till_death = 3; // XXX: VERY HEAVY. Initial survivors should be calculated carefully. Even "13" is big enough. "20" won't fit in 32 GB.
	                                       //      But! This gives an unparalleled variety of specimen to algo.
	selection->set_generations_till_elimination(gens_till_death);

	genetics::GenerationsCosts<PolynomialCost> how_fit;

	char ans;
	std::cout << "\nWould you like to see the first settlers (y/N)? " << std::flush;
	std::cin.get();
	std::cin.get(ans);
	if (ans == 'y' || ans == 'Y') {
		fitness->cost(the_nonglitch, how_fit);
		auto& specimens = std::get<0>(the_nonglitch);
		for (size_t i = 0; i < specimens.size(); ++i)
			std::cout << listing(specimens.at(i)) << "; fitness = " << how_fit[i] << std::endl;
	}
	std::cout << std::endl;

	std::size_t matings = survivors * (survivors - 1);
	for (std::size_t i = 1; i < gens_till_death; ++i) matings *= matings - 1;

	// Tries to write a "." every time when another 1 hundred million matings occured.
	if (matings > 100'000'000) {
		double how_many_gens_for_10k_matings = 100'000'000.0 / matings;

		std::size_t new_gen_cap = std::max(static_cast<std::size_t>(how_many_gens_for_10k_matings), gens_till_death);
		double generation_pool = 0.0;

		std::cout << "Wait for world history to happen " << std::flush;

		selection->set_max_generations(new_gen_cap);

		size_t full_iterations = generation_cap / new_gen_cap;

		for (size_t i = 0; i < full_iterations; ++i) {
			the_nonglitch = world.evolve(std::move(the_nonglitch));
			generation_pool += new_gen_cap;
			while (generation_pool >= how_many_gens_for_10k_matings) {
				std::cout << ".";
				generation_pool -= how_many_gens_for_10k_matings;
			}
			std::cout << std::flush;
		}

		// Last gens
		selection->set_max_generations(generation_cap - full_iterations * new_gen_cap);
		the_nonglitch = world.evolve(std::move(the_nonglitch));
		std::cout << std::endl << std::endl;
	} else {
		the_nonglitch = world.evolve(std::move(the_nonglitch));
	}

	how_fit.clear();
	fitness->cost(the_nonglitch, how_fit);

	if (how_fit[0].get_inaccuracy() > 0.001) {
		std::cout << "Civilization of the Nonglitch fell, unable to match your goal." << std::endl << std::endl;
	} else {
		std::cout << "The Nonglitch ascended." << std::endl << std::endl;
	}

	auto& specimens = std::get<0>(the_nonglitch);
	auto& ages = std::get<1>(the_nonglitch);

	std::cout << "Best match:" << std::endl;
	for (auto& testcase : target) {
		std::cout << '\t' << testcase.first << " -> " << interpret(specimens.at(0), testcase.first) << std::endl;
	}
	std::cout << "That fits like: " << how_fit[0] << std::endl;
	std::cout << "Author: " << listing(specimens.at(0)) << std::endl;
	std::cout << "Age: " << ages.value()[0] << std::endl;

	std::cout << std::endl << "Other last survivors:" << std::endl;
	for (size_t i = 1; i < std::min(specimens.size(), 20lu); ++i) {
		std::cout << listing(specimens.at(i)) << ";\n\tfitness = " << how_fit[i] << ";\n\tage = " << ages.value()[i] << ";" << std::endl;
	}
}
