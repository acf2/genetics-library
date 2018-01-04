#ifndef __GENETICS_H
#define __GENETICS_H

#include <cstdlib>
#include <cstddef>
#include <vector>

#include <map>
#include <functional>
#include <algorithm>
#include <random>
#include <cassert>
#include <utility>


namespace Genetics {

	// Достаточно синглтона Мэйерса
	// Когда буду делать более навороченную штуку, можно будет подумать и о многопоточном
	// ...
	// ... он и так многопоточный?..
	class RandomFloatGenerator {
	private:
		std::mt19937 engine;
		std::uniform_real_distribution<float> distribution;

		RandomFloatGenerator() : distribution(0.0, 1.0) {
			std::random_device rd;
			engine.seed(rd());
		}
		RandomFloatGenerator(RandomFloatGenerator const&) = delete;
		RandomFloatGenerator& operator=(RandomFloatGenerator&) = delete;
		~RandomFloatGenerator() { }
	public:
		static RandomFloatGenerator& get_instance() {
			static RandomFloatGenerator instance;
			return instance;
		}

		float get_random_float() {
			return distribution(engine);
		}
	};

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
		friend void swap(Population& one, Population& another) {
			using std::swap;
			swap(one.specimens, another.specimens);
			swap(one.generation_number, another.generation_number);
		}
		Population() : generation_number(0) {}
		Population(std::vector<DNA> specimens, size_t generation_number = 0) : specimens(specimens), generation_number(generation_number) {}
		Population(Population const& another) : specimens(another.specimens), generation_number(another.generation_number) { }
		Population(Population&& another) {
			swap(*this, another);
		}
		~Population() { }

		Population& operator=(Population another) {
			swap(*this, another);
			return *this;
		}

		Population& operator=(Population&& another) {
			swap(*this, another);
			return *this;
		}

		// Не лучшая идея, но мне нравится считать разыменование популяции вектором из особей
		std::vector<DNA>& operator*() {
			return specimens;
		}

		std::vector<DNA>* operator->() {
			return &specimens;
		}

		size_t generation() const {
			return generation_number;
		}
	};

	template<typename DNA, typename Fitness_t>
	class World {
	private:
		std::function<Fitness_t(DNA const&)> fitness;
		std::function<DNA(DNA const&, DNA const&)> crossover;
		bool symmetric_crossover;
		std::function<void(DNA*)> mutate;

	public:
		static const size_t unlimited = 0;

		World(
			std::function<Fitness_t(DNA const&)> fitness,
			std::function<DNA(DNA const&, DNA const&)> crossover,
			bool symmetric_crossover,
			std::function<void(DNA*)> mutate
		) : fitness(fitness),
			crossover(crossover),
			symmetric_crossover(symmetric_crossover),
			mutate(mutate) { }
		~World() = default;

		Population<DNA> evolve(
			// Первые родители
			Population<DNA> generation0,
			// При каком уровне приспособленности остановить выполнение
			// Ноль, или эквивалентный ему уровень - вплоть до идеального представителя вида
			Fitness_t required_fitness,
			// Вероятность возникновения мутации (0 <= mutation_probability <= 1)
			// Значения больше 1 приравниваются к 1, меньше 0 приравниваются к 0
			float mutation_probability = 0.05,
			// Количество особей, что отбираются для порождения следующего поколения
			// При нуле принимает значение количества особей в первом поколении
			size_t parents_size = unlimited,
			// На сколько поколений можно максимально ЕЩЁ продвинуться вперёд. Ноль - без ограничений
			size_t max_generations = unlimited
		) {
			bool is_generation_based = max_generations != 0;
			size_t active_parents;
			if (parents_size == 0)
				parents_size = generation0->size();

			std::vector<DNA> Parents(parents_size);
			active_parents = std::min(parents_size, generation0->size());
			std::copy_n(generation0->begin(), active_parents, Parents.begin());

			std::multimap<Fitness_t, DNA> Field;
			size_t generation_number = 0;
			Fitness_t achieved_fitness = fitness(Parents[0]);
			for (size_t i = 1; i < active_parents; ++i)
				if (fitness(Parents[i]) < achieved_fitness)
					achieved_fitness = fitness(Parents[i]);

			// Цикл скрещивания
			// пока (генетический алгоритм ⇒ поколение < максимум поколений) и (необходимая приспособленность < достигнутой приспособленности)
			// Выход в середине тела цикла
			DNA crossbuffer;
			size_t crossbuffer_fitness;
			RandomFloatGenerator& generator = RandomFloatGenerator::get_instance();

			for (;;) {
				Field.clear();

				//Скрещивание
				for (size_t i = 0; i < active_parents; ++i) {
					for (size_t j = (symmetric_crossover ? i+1 : 0); j < active_parents; ++j) {
						crossbuffer = std::move(crossover(Parents[i], Parents[j]));
						if (generator.get_random_float() <= mutation_probability)
							mutate(&crossbuffer);
						crossbuffer_fitness = fitness(crossbuffer);
						Field.emplace(crossbuffer_fitness, std::move(crossbuffer));
					}
				}

				achieved_fitness = Field.begin()->first; //Лучший
				++generation_number;

				if ((is_generation_based && generation_number >= max_generations) || !(required_fitness < achieved_fitness))
					break;
				
				size_t parent_number = 0;
				for (typename std::multimap<Fitness_t, DNA>::iterator iter = Field.begin(); parent_number < Parents.size() && iter != Field.end(); ++parent_number, ++iter) {
					Parents[parent_number] = std::move(iter->second);
				}
				active_parents = parent_number;

			}
			Population<DNA> result(
				std::vector<DNA>(Field.size()), 
				generation0.generation() + generation_number);
			std::transform(
				Field.begin(), 
				Field.end(), 
				result->begin(), 
				[](std::pair<Fitness_t, DNA> p) -> DNA { return std::move(p.second); });
			return result;
		}
	};
}

#endif //__GENETICS_H
