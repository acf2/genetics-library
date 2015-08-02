#ifndef __GNTC_H__
#define __GNTC_H__

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <string>
#include <functional>
#include <map>
#include <cassert>

#include <iostream>

/*
 *	Простейшая, наивная реализация генетического алгоритма на шаблоне.
 *	Перед использованием надо инициализировать ГПСЧ. 
 */

namespace Gntc {

	// Ошибки понадобятся для расширения
	class Error {
	public:
		virtual std::string what() = 0;
	};

	class UnexpectedOrder : public Error {
	private:
		std::string msg;
	public:
		UnexpectedOrder(std::string S) : msg(S) { }
		std::string what() override {
			return("Unexpected order of launch in " + msg);
		}
	};
	// --------

	inline size_t min_field(size_t parents) {
		return(((parents+1)*parents)/2);
	}

	inline float random_chance() {
		return(float(rand()%1000)/1000.0);
	}

	template< typename Dna >
	class World {
	private:
		
		size_t max_generations; // Если 0, то неограниченно
		size_t required_fitness; // Если 0, то вплоть до идеального представителя вида

		// Здесь переменные, хранящие функции cross и fitness
		std::function< Dna(Dna const&, Dna const&) > cross;
		std::function< void(Dna&) > mutate;
		std::function< size_t(Dna const&) > fitness;
	
		std::vector< Dna > generation0;

		bool working;
		std::vector< Dna > best_cache;
		size_t gennum_cache;

	public:
		/*
		 *	num_generation0 - размер нулевого поколения
		 *	generation0 - массив с представителями нулевого поколения
		 */
		World(
			size_t num_generation0_ini,
			Dna* generation0_ini,
			std::function< Dna(Dna const&, Dna const&) > cross_ini,
			std::function< void(Dna&) > mutate_ini,
			std::function< size_t(Dna const&) > fitness_ini,
			size_t required_fitness_ini
		) : cross(cross_ini),
			mutate(mutate_ini),
			fitness(fitness_ini),
			required_fitness(required_fitness_ini)
		{
			generation0.assign(generation0_ini, generation0_ini + num_generation0_ini);
			best_cache.clear();
			gennum_cache = 0;
			working = 0;
		}
		~World() = default;

		Dna ChaseDream(
			size_t parents_size = 0,
			float mutation_probability = 0.05,
			bool symmetric = 1,
			size_t max_generations = 0
			size_t 
		) {
			bool genbased_algo = max_generations != 0;
			
			if (parents_size == 0)
				parents_size = generation0.size();
			
			std::vector< Dna > Parents(generation0);
			std::multimap< size_t, Dna > Field;

			size_t generation = 0;
			size_t mutual_fitness = size_t(-1);
			for (size_t i = 0; i < Parents.size(); ++i)
				if (mutual_fitness > fitness(Parents[i]))
					mutual_fitness = fitness(Parents[i]);

			// Цикл скрещивания
			// (генетический алгоритм ⇒ поколение < максимум поколений) и (приспособленность >= необходимой приспособленности)
			Dna crossbuffer;
			typename std::multimap< size_t, Dna >::iterator iter;
			
			while ((!genbased_algo || generation < max_generations) && (mutual_fitness > required_fitness)) {
				Field.clear();

				//Скрещивание
				for (size_t i = 0; i < Parents.size(); ++i) {
					for (size_t j = (symmetric ? i+1 : 0); j < Parents.size(); ++j) {
						crossbuffer = cross(Parents[i], Parents[j]);
						if ( random_chance() <= mutation_probability )
							mutate(crossbuffer);
						Field.emplace(fitness(crossbuffer), crossbuffer);
					}
				}

				iter = Field.begin();
				mutual_fitness = iter->first; //Лучший
				++generation;
				for (size_t i = 0; i < Parents.size(); ++i) {
					assert(iter != Field.end());
					Parents[i] = iter->second;
					++iter;
				}

			}
			best_cache = Parents;
			gennum_cache = generation;
			return Parents[0];
		}
	
		size_t acc_gennum() {
			return gennum_cache;
		}
	};

}

#endif //__GNTC_H__
