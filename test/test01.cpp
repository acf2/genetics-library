#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iterator>

#include "../include/genetics.hpp"

using namespace std;

typedef long long int gene_t;
typedef std::size_t cost_t;

class Fitness : public genetics::IFitness<gene_t, cost_t> {
public:
	~Fitness() override = default;

	void cost(genetics::Generation<gene_t> const& generation, genetics::GenerationsCosts<cost_t>& costs) override {
		auto& specimens = std::get<genetics::SPECIMENS_ID>(generation);
		costs.resize(specimens.size());
		std::transform(begin(specimens), end(specimens), begin(costs), [](gene_t const& specimen) -> cost_t {
			gene_t constexpr ideal = 100;
			return static_cast<cost_t>(std::max(specimen, ideal) - std::min(specimen, ideal));
		});
	}
};

class Crossover : public genetics::ICrossover<gene_t, cost_t> {
public:
	~Crossover() override = default;

	bool does_commute() const override { return true; }

	gene_t cross(genetics::Generation<gene_t> const& generation, genetics::GenerationsCosts<cost_t> const& costs, std::size_t parent1, std::size_t parent2) {
		auto& specimens = std::get<genetics::SPECIMENS_ID>(generation);
		gene_t new_specimen = (specimens[parent1] + specimens[parent2]) / 2;

		// Mutation
		if (genetics::RandomGenerator::get_instance().get_random_float<float>(0.f, 1.f) < 0.3f) {
			if (genetics::RandomGenerator::get_instance().get_random_int<int>(0, 1)) {
				--new_specimen;
			} else {
				++new_specimen;
			}
		}	

		return new_specimen;
	}
};

class Selection : public genetics::ISelection<gene_t, cost_t> {
public:
	~Selection() override = default;

	std::size_t survivors() const override { return 100; };

	std::optional<std::size_t> max_generations() const override { return 10; }
};

int main(int argc, char const *argv[]) {
	genetics::Environment<gene_t, cost_t> env(std::make_shared<Fitness>(), std::make_shared<Crossover>(), std::make_shared<Selection>());

	auto first_gen = genetics::new_generation<gene_t>({1, 888, 42, 0xDEADCAFE});

	auto result = env.evolve(first_gen);

	auto& specimens0 = std::get<genetics::SPECIMENS_ID>(first_gen);
	auto& specimens1 = std::get<genetics::SPECIMENS_ID>(result);

	cout << "GENERATIONS: " << std::get<genetics::GENERATION_COUNT_ID>(result) << endl;
	for (auto i = specimens0.begin(); i != specimens0.end(); ++i) {
		cout << *i << ' ';
	}
	cout << endl;

	for (auto i = specimens1.begin(); i != specimens1.end(); ++i) {
		cout << *i << ' ';
	}
	cout << endl;

	return 0;
}
