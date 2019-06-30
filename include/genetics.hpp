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

	const size_t unlimited = 0;

	template<typename DNA, typename Fitness_t>
	struct WorldSettings {
		std::function<std::vector<Fitness_t>(std::vector<DNA> const&)> fitness;
		std::function<DNA(DNA const&, DNA const&)> crossover;
		bool is_crossover_symmetric = 1;
		std::function<DNA(DNA&&)> mutate;
		std::function<bool(DNA const&, DNA const&)> is_similar_to = [](DNA const& one, DNA const& another) -> bool { return 0; };

		// 0 <= mutation_probability <= 1
		// if > 1 then 1; if < 0 then 0
		float mutation_probability;
		// Number of specimens selected from each generation.
		// (zero or "unlimited" for unlimited version, HIGHLY not recommended)
		size_t survivors;
		// Selection will be performed every <death_rate> generations.
		size_t death_rate = 1;
	};

	template<typename DNA, typename Fitness_t>
	Population<DNA> evolve(
		WorldSettings<DNA, Fitness_t> settings,
		Population<DNA> generation0,
		std::function<bool(Fitness_t)> good_enough, // Check for target fitness
		size_t max_generations = unlimited
	) {
		bool is_generation_based = max_generations != unlimited;
		decltype(settings.mutation_probability) mutation_probability = std::max(0.f, std::min(1.f, settings.mutation_probability));

		// field это базовый контейнер, с которым будет работать evolve
		std::multimap<Fitness_t, DNA> field;
		// newbies нужна для обновления field (см. ниже)
		std::vector<DNA> newbies;
		/*
		 *	Из-за природы функции fitness (которая может зависеть и от текущей популяции), после
		 *	получения всех особей некоторого поколения необходимо пересчитать очки приспособленности
		 *	Для более-менее оптимального обращения с памятью при заметных размерах типа DNA и дороговизны вызова fitness
		 *	будет использован следующий алгоритм:
		 *		0. Имеется старое поле field
		 *		1. Заполняется вектор vector<DNA> newbies с помощью кроссовера (c std::move, поэтому копирования нет)
		 *		2. Все родители из field добавляются в newbies с помощью std::move
		 *		3. field очищается clear
		 *		4. Пересчёт
		 *			4.1 Для всех newbies одновременно вычисляется значение fitness
		 *			    (поэтому его производительность зависит от пользовательского коллбэка)
		 *			4.2 Каждый из newbies заносится в новое field c помощью std::move вместе с его значением fitness
		 *		5. Получено новое поле field
		 *
		 *	Функция update_field и будет выполнять этот алгоритм после заполнения вектора newbies.
		 *	Она принимает на вход старое field и массив newbies, возвращает новое filed и опустошает newbies (с помощью clear)
		 */
		std::function<void()> update_field = [&settings, &field, &newbies]() -> void {
			newbies.reserve(newbies.size() + field.size());
			for (auto&& elem : field) {
				newbies.push_back(std::move(elem.second));
			}
			std::vector<Fitness_t> newbies_fitnesses = settings.fitness(newbies);

			field.clear();
			for (size_t i = 0; i < newbies.size(); ++i) {
				field.emplace(std::move(newbies_fitnesses[i]), std::move(newbies[i]));
			}
			newbies.clear();
		};

		// Вот так это будет происходить
		newbies = *generation0;
		update_field();

		size_t generation_number = generation0.generation();
		size_t end_generation = generation_number + max_generations;
		Fitness_t achieved_fitness = field.begin()->first;

		// Цикл скрещивания
		// пока (генетический алгоритм ⇒ поколение < максимум поколений) и (необходимая приспособленность < достигнутой приспособленности)
		// Выход в середине тела цикла
		RandomFloatGenerator& generator = RandomFloatGenerator::get_instance();

		while ((!is_generation_based || generation_number < end_generation) && !good_enough(achieved_fitness)) {
			// Crossover
			for (auto first_specimen = field.begin(); first_specimen != field.end(); ++first_specimen) {
				for (auto second_specimen = (settings.is_crossover_symmetric ? std::next(first_specimen) : field.begin()); second_specimen != field.end(); ++second_specimen) {
					newbies.emplace_back(std::move(settings.crossover(first_specimen->second, second_specimen->second)));
					if (generator.get_random_float() <= mutation_probability)
						newbies.back() = settings.mutate(std::move(newbies.back()));
				}
			}
			update_field();

			++generation_number;

			// Selection will be done in the following way:
			//   1) Select the most fit first and pull it from the pool
			//   2) Discard every specimen that "is similar to" previously selected
			//   3) Go to 1
			if (settings.survivors != unlimited && field.size() > settings.survivors && generation_number % settings.death_rate == 0) {
				std::multimap<Fitness_t, DNA> old_field = std::move(field);
				field.clear();

				size_t i;
				for (i = 0; i < settings.survivors && !old_field.empty(); ++i) {
					auto&& newcomer = *std::begin(old_field);

					field.insert(newcomer);

					// XXX: remove_if does not work with maps!
					for (auto iter = std::begin(old_field); iter != std::end(old_field);) {
						if (settings.is_similar_to(iter->second, newcomer.second)) {
							iter = old_field.erase(iter);
						} else {
							++iter;
						}
					}
				}

				// If there is not enough survivors, just clone the best one
				for (; i < settings.survivors; ++i) {
					assert(!field.empty());
					field.insert(*std::begin(field));
				}
			}

			achieved_fitness = field.begin()->first;
		}

		Population<DNA> result(
			std::vector<DNA>(field.size()),
			generation_number);
		std::transform(field.begin(), field.end(), result->begin(),
					   [](std::pair<Fitness_t, DNA> p) -> DNA {
						   return std::move(p.second);
					   });
		return result;
	}
}

#endif //__GENETICS_H
