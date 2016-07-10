#ifndef __GENETICS_H
#define __GENETICS_H

#include <cstdlib>
#include <cstddef>
#include <vector>
//#include <string>
#include <map>
#include <functional>
#include <algorithm>
#include <chrono>
#include <random>
#include <cassert>

//#include <iostream>

/*
 *	Простейшая, наивная реализация генетического алгоритма на шаблоне.
 *	Перед использованием надо инициализировать ГПСЧ.
 */

namespace Genetics {

	inline float random_chance() {
		return(static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
	}

	/*
	 *	Класс популяции служит для снижения сложности восприятия работы с особями.
	 *	Проблема в том, что номер поколения и прочие данные порою должны быть связаны.
	 */
	template<typename DNA>
	class Population {
	private:
		std::vector<DNA> specimens;
		size_t generation_number;

	public:
		Population() : generation_number(0) {}
		Population(Population const& another) : specimens(another.specimens), generation_number(another.generation_number) { }
		Population(std::vector<DNA> specimens, size_t generation_number) : specimens(specimens), generation_number(generation_number) {}
		~Population() { }

		Population& operator=(Population const& another) {
			specimens = another.specimens;
			generation_number = another.generation_number;
		}

		// Не лучшая идея, но мне нравится считать разыменование популяции вектором из особей
		std::vector<DNA>& operator*() {
			return specimens;
		}

		std::vector<DNA>* operator->() {
			return &specimens;
		}

		size_t get_generation_number() const {
			return generation_number;
		}
	};

	template<typename Fitness_t, typename DNA>
	class World {
	private:
		std::function<DNA(DNA const&, DNA const&)> crossover;
		bool symmetric_crossover;
		std::function<void(DNA*)> mutate;
		std::function<Fitness_t(DNA const&)> fitness;

	public:
		World(
			std::function<DNA(DNA const&, DNA const&)> crossover,
			std::function<void(DNA*)> mutate,
			std::function<Fitness_t(DNA const&)> fitness,
			bool symmetric_crossover = 1
		) : crossover(crossover),
			mutate(mutate),
			fitness(fitness),
			symmetric_crossover(symmetric_crossover) { }
		~World() = default;

		Population Evolve(
			// Предтечи
			Population generation0,
			// При каком уровне приспособленности остановить выполнение
			// Ноль, или эквивалентный ему уровень - вплоть до идеального представителя вида
			Fitness_t required_fitness,
			// Количество особей, что отбираются для порождения следующего поколения
			// При нуле принимает значение количества особей в первом поколении
			size_t parents_size = 0,
			// Вероятность возникновения мутации (0 <= mutation_probability <= 1)
			// Значения больше 1 приравниваются к 1, меньше 0 приравниваются к 0
			float mutation_probability = 0.05,
			// На сколько поколений можно максимально ЕЩЁ продвинуться вперёд. Ноль - без ограничений
			size_t max_generations = 0
		) {
			bool is_generation_based = max_generations != 0;
			if (parents_size == 0)
				parents_size = generation0.size();

			//std::vector<DNA> Parents(generation0);
			std::vector<DNA> Parents(parents_size);
			std::copy_n(generation0->begin(), parents_size, Parents.begin());

			std::multimap<Fitness_t, DNA> Field;
			size_t generation_number = 0;
			Fitness_t achieved_fitness = Parents[0];
			for (size_t i = 1; i < Parents.size(); ++i)
				if (fitness(Parents[i]) < achieved_fitness)
					achieved_fitness = fitness(Parents[i]);

			// Цикл скрещивания
			// (генетический алгоритм ⇒ поколение < максимум поколений) и (необходимая приспособленность < достигнутой приспособленности)
			DNA crossbuffer;
			typename std::multimap<Fitness_t, DNA>::iterator iter;

			while ((!is_generation_based || generation_number < max_generations) && (required_fitness < achieved_fitness)) {
				Field.clear();

				//Скрещивание
				for (size_t i = 0; i < Parents.size(); ++i) {
					for (size_t j = (symmetric ? i+1 : 0); j < Parents.size(); ++j) {
						crossbuffer = crossover(Parents[i], Parents[j]);
						if (random_chance() <= mutation_probability)
							mutate(&crossbuffer);
						Field.emplace(fitness(crossbuffer), crossbuffer);
					}
				}

				iter = Field.begin();
				achieved_fitness = iter->first; //Лучший
				++generation_number;
				for (size_t i = 0; i < Parents.size(); ++i) {
					assert(iter != Field.end());
					Parents[i] = iter->second;
					++iter;
				}

			}
			Population result(std::vector<DNA>(Field.size()), generation0.get_generation_number() + generation_number);
			return Parents[0];
		}

		size_t acc_gennum() {
			return gennum_cache;
		}
	};

}

#endif //__GENETICS_H
