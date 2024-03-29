#ifndef __GENETICS_H
#define __GENETICS_H

// C++20 is expected

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <optional>
#include <memory>

#include <map>
#include <functional>
#include <algorithm>
#include <random>
#include <cassert>
#include <utility>

namespace genetics {

// Stolen & modified from here: https://stackoverflow.com/a/17074810
// About permutations and The Great Fukken Confusion: https://en.wikiversity.org/wiki/Permutation_notation
//   hint: it's clear and obvious only until you try to reason with someone else about it.
template <typename T, typename Compare>
std::vector<std::size_t> sort_to_permutation(const std::vector<T>& vec, Compare&& compare) {
	std::vector<std::size_t> passive_permutation(vec.size()), active_permutation(vec.size());
	std::iota(passive_permutation.begin(), passive_permutation.end(), 0);
	std::sort(passive_permutation.begin(), passive_permutation.end(),
	          [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
	for (std::size_t i = 0; i < passive_permutation.size(); ++i) {
		active_permutation[passive_permutation[i]] = i;
	}
	return active_permutation;
}

// I can't believe C++ does not have something like this out-of-the-box
template <typename T>
void apply_permutation_in_place(std::vector<T>& vec, std::vector<std::size_t> permutation) {
	std::size_t current_place = 0;
	while (current_place < vec.size()) {
		if (permutation[current_place] == current_place) {
			++current_place;
			continue;
		}
		std::swap(vec[current_place], vec[permutation[current_place]]);
		std::swap(permutation[current_place], permutation[permutation[current_place]]);
	}
}


class RandomGenerator {
private:
	std::mt19937 engine;

	RandomGenerator() {
		std::random_device rd;
		engine.seed(rd());
	}
	RandomGenerator(RandomGenerator const&) = delete;
	RandomGenerator& operator=(RandomGenerator&) = delete;
	~RandomGenerator() { }

public:
	static RandomGenerator& get_instance() {
		static RandomGenerator instance;
		return instance;
	}

	template <typename T>
	T get_random_float(T from, T to) {
		return std::uniform_real_distribution<T>(from, to)(engine);
	}

	template <typename T>
	T get_random_int(T from, T to) {
		return std::uniform_int_distribution<T>(from, to)(engine);
	}
};

std::size_t constexpr SPECIMENS_ID = 0,
                      GENERATION_COUNT_ID = 1,
                      AGE_ID = 2;

template <typename Genome>
using Generation = std::tuple<std::vector<Genome>,
                              std::size_t,
                              std::optional<std::vector<std::size_t>>>;

template <typename Genome>
Generation<Genome> new_generation(std::vector<Genome> specimens) {
	std::size_t const specimens_size = specimens.size();
	return std::make_tuple(std::move(specimens), 0, std::vector<std::size_t>(specimens_size, 0));
}

template <typename Cost>
using GenerationsCosts = std::vector<Cost>;

template <typename Genome, typename Cost>
class IFitness {
public:
	virtual ~IFitness() = default;

	// All-at-once computation, because from framework POV it parallelizes poorly
	// Also, must be computed "in place": the most complicated part for framework user
	// Framework guarantees that size of costs will be 0, and it's capacity will be somewhat close to optimal
	virtual void cost(Generation<Genome> const& generation, GenerationsCosts<Cost>& costs) = 0;
};

/*
 *	Could be controlled by Genome class.
 *	Must be *very* fast, because it is expected to be done synchronously
 *	in sequential manner.
 */
template <typename Genome, typename Cost>
class ICrossover {
public:
	virtual ~ICrossover() = default;

	virtual bool does_commute() const = 0; // Is crossover symmetric?
	virtual std::size_t default_offspring_amount() const { return 1; } // How much offspring should a pair have by default?

	// How to scale offspring amount due to high fitness? - This is a must, if one wants to implement specimen aging
	//     ...after all, this is how real evolution works
	virtual std::size_t offspring_amount(Generation<Genome> const& generation, GenerationsCosts<Cost> const& costs, std::size_t parent1, std::size_t parent2) { return 1; }
	
	// What mutations will offspring have? Probably this should be decided/done by crossover operation
	
	virtual Genome cross(Generation<Genome> const& generation, GenerationsCosts<Cost> const& costs, std::size_t parent1, std::size_t parent2) = 0;
};

template <typename Genome, typename Cost>
class ISelection {
public:
	virtual ~ISelection() = default;

	// Number of specimens selected from each generation.
	virtual std::size_t survivors() const = 0;

	// Check for some specimen with good enough fitness
	virtual bool is_good_enough(Genome const& specimen, Cost const& cost) { return false; }
	
	// Limit of generations for one evolve call
	virtual std::optional<std::size_t> max_generations() const { return std::nullopt; }

	// Selection will be performed only every generations_till_eliminaion generations
	virtual std::size_t generations_till_eliminaion() const { return 1; }
};

// XXX: Probably should be distributed among crossover and selection
/*
template <typename Genome, typename Cost>
class ISimilarity {
public:
	virtual ~ISimilarity() = default;

	// After choosing specimen with the best fitness, add penalties for being similar to already chosen pool
	//virtual bool discourage_similarity() { return false; }
};
*/

template <typename Genome, typename Cost>
class Environment {
public:
	Environment(std::shared_ptr<IFitness<Genome, Cost>> fitness,
	            std::shared_ptr<ICrossover<Genome, Cost>> crossover,
	            std::shared_ptr<ISelection<Genome, Cost>> selection,
	            //std::shared_ptr<ISimilarity<Genome, Cost>> similarity,
	            std::size_t number_of_threads = 1)
	  : fitness(fitness)
	  , crossover(crossover)
	  , selection(selection)
	  //, similarity(similarity)
	  , number_of_threads(number_of_threads) { }

	// Main function
	Generation<Genome> evolve(Generation<Genome> generation) {
		auto& specimens = std::get<SPECIMENS_ID>(generation);
		auto& generation_count = std::get<GENERATION_COUNT_ID>(generation);
		auto& ages = std::get<AGE_ID>(generation);
		GenerationsCosts<Cost> costs;

		std::size_t approximate_size_of_generation_container = compute_approximate_size_of_generation_container();

		specimens.reserve(approximate_size_of_generation_container);
		if (ages) {
			ages.value().reserve(approximate_size_of_generation_container);
		}
		costs.reserve(approximate_size_of_generation_container);

		// Compute costs first time
		compute_fitness(generation, costs);

		auto starting_generation = generation_count;
		auto const max_generations = selection->max_generations();
		for (;;) {
			compute_crossover(generation, costs);
			++generation_count;

			if (generation_count % selection->generations_till_eliminaion() == 0) {
				compute_fitness(generation, costs);
				eliminate_losers(generation, costs);

				// NOTE: Why? Because first specimen must be the best, after sorting
				if (selection->is_good_enough(specimens[0], costs[0])) {
					break;
				}
			}

			if (max_generations && (generation_count - starting_generation) >= max_generations.value()) {
				break;
			}
		}

		return generation;
	}

private:
	std::size_t number_of_threads;

	std::shared_ptr<IFitness<Genome, Cost>> fitness;
	std::shared_ptr<ICrossover<Genome, Cost>> crossover;
	std::shared_ptr<ISelection<Genome, Cost>> selection;
	//std::shared_ptr<ISimilarity<Genome, Cost>> similarity;

	// Calculate approximate size of generation container (to minimize reallocations)
	// offspring_amount is not considered, but it should converge to optimal capacity quickly (provided that resize() of vector does not change capacity)
	std::size_t compute_approximate_size_of_generation_container() const {
		std::size_t generation_size_estimation = selection->survivors();

		// For each generation, when elimination did not occur, we must cumulatively multiply required space
		// It will be something akin to survivors ^ (2 ^ generations): increases VERY fast
		// Each generation:
		//   - All specimen crossover with all specimens (n * n)
		//   - But! Specimens do not crossover with themselves (- n)
		//   - And if crossover commutes, we don't need half of that (/ 2)
		for (std::size_t i = 0; i < selection->generations_till_eliminaion(); ++i) {
			std::size_t new_specimens = generation_size_estimation
			                            * generation_size_estimation
			                            - generation_size_estimation;
			if (crossover->does_commute()) new_specimens /= 2;
			new_specimens *= crossover->default_offspring_amount();

			generation_size_estimation += new_specimens;
		}
		return generation_size_estimation;
	}

	// modifies generation
	void compute_crossover(Generation<Genome>& generation, GenerationsCosts<Cost> const& costs) {
		// Sequential for now, but it should not be heavy
		// TODO: Redesign for any number of threads
		// NOTES:
		//   - Commutative crossover should not use the same specimens twice: use only half of the matrix
		//   - Default offspring amount controls how many times crossover for this particular pair is called
		//   - Offspring amount shows how much more offspring should a pair have, if it's fitness is high (close to zero)
		//   - Each crossover phase must increase age of specimen (because it is used when computing fitness)
		//     In short? Parents age, when they give birth to offspring
		std::size_t const default_offspring = crossover->default_offspring_amount();

		auto& specimens = std::get<SPECIMENS_ID>(generation);
		auto& ages = std::get<AGE_ID>(generation);

		std::size_t const specimens_old_size = specimens.size();

		if (ages) {
			for (std::size_t i = 0; i < ages.value().size(); ++i) {
				++ages.value()[i];
			}
		}

		for (std::size_t i = 0; i < specimens_old_size; ++i) {
			std::size_t starting_specimen = crossover->does_commute() ? 0 : i+1;

			for (std::size_t j = i+1; j < specimens_old_size; ++j) {
				if (i == j) continue;

				std::size_t const current_offspring_amount = crossover->offspring_amount(generation, costs, i, j) * default_offspring;

				for (std::size_t count = 0; count < current_offspring_amount; ++count) {
					specimens.emplace_back(crossover->cross(generation, costs, i, j)); // One cannot be sure if Genome is heavy or not
					if (ages) {
						ages.value().emplace_back(0);
					}
				}
			}
		}
	}

	// modifies costs
	// Sequential! User MUST add parallelization here by himself.
	// It's the most heavy function of all, but it depends on all specimens at once
	// thus, it's not possible to parallelize this in framework
	void compute_fitness(Generation<Genome> const& generation, GenerationsCosts<Cost>& costs) {
		costs.resize(0);
		fitness->cost(generation, costs);
	}

	// modifies generation
	void eliminate_losers(Generation<Genome>& generation, GenerationsCosts<Cost> const& costs) {
		// Sequential for now, but it should not be heavy
		// TODO: Redesign for any number of threads
		// NOTES:
		//   - All "winners" should be moved/swapped to first positions
		//   - Then container should be resized to target size
		auto& specimens = std::get<SPECIMENS_ID>(generation);
		auto& ages = std::get<AGE_ID>(generation);

		auto&& sort_permutation = sort_to_permutation(costs, [](Cost const& one, Cost const& another) { return one < another; });

		apply_permutation_in_place(specimens, sort_permutation);
		if (specimens.size() > selection->survivors()) {
			specimens.resize(selection->survivors());
		}
		if (ages) {
			apply_permutation_in_place(ages.value(), sort_permutation);
			if (ages.value().size() > selection->survivors()) {
				ages.value().resize(selection->survivors());
			}
		}
	}
};

} // namespace genetics

#endif //__GENETICS_H

