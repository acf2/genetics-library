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
	class World {
	private:
		std::function<Fitness_t(std::vector<DNA> const&, size_t)> fitness;
		std::function<DNA(DNA const&, DNA const&)> crossover;
		bool symmetric_crossover;
		std::function<void(DNA&)> mutate;

	public:
		friend void swap(World& one, World& another) {
			using std::swap;
			swap(one.fitness, another.fitness);
			swap(one.crossover, another.crossover);
			swap(one.symmetric_crossover, another.symmetric_crossover);
			swap(one.mutate, another.mutate);
		}

		World(
			std::function<Fitness_t(std::vector<DNA> const&, size_t)> fitness,
			std::function<DNA(DNA const&, DNA const&)> crossover,
			bool symmetric_crossover,
			std::function<void(DNA&)> mutate
		) : fitness(fitness),
			crossover(crossover),
			symmetric_crossover(symmetric_crossover),
			mutate(mutate) { }

		World(World const& another)
		  : fitness(another.fitness),
			crossover(another.crossover),
			symmetric_crossover(another.symmetric_crossover),
			mutate(another.mutate) { }

		World(World&& another) {
			swap(*this, another);
		}
		~World() = default;

		World& operator=(World another) {
			swap(*this, another);
			return *this;
		}

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
			// При нуле - неограничено
			size_t survivors = unlimited,
			// На сколько поколений можно максимально ЕЩЁ продвинуться вперёд. Ноль - без ограничений
			size_t max_generations = unlimited,
			// Через сколько поколений отделять выживших от мертвецов
			size_t death_rate = 1
		) {
			bool is_generation_based = max_generations != unlimited;

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
			 *			4.1 Для каждого из newbies один раз вычисляется значение fitness
			 *			4.2 Каждый из newbies заносится в новое field c помощью std::move вместе с его значением fitness
			 *		5. Получено новое поле field
			 *
			 *	Функция update_field и будет выполнять этот алгоритм после заполнения вектора newbies.
			 *	Она принимает на вход старое field и массив newbies, возвращает новое filed и опустошает newbies (с помощью clear)
			 */
			std::function<void()> update_field = [this, &field, &newbies]() -> void {
				newbies.reserve(newbies.size() + field.size());
				for (auto& elem : field) {
					newbies.push_back(std::move(elem.second));
				}
				std::vector<Fitness_t> newbies_fitnesses;
				newbies_fitnesses.reserve(newbies.size());
				for (size_t i = 0; i < newbies.size(); ++i) {
					newbies_fitnesses.push_back(std::move(fitness(newbies, i)));
				}

				field.clear();
				for (size_t i = 0; i < newbies.size(); ++i) {
					field.emplace(std::move(newbies_fitnesses[i]), std::move(newbies[i]));
				}
				newbies.clear();
			};

			// Вот так это будет происходить
			newbies = *generation0;
			update_field();

			size_t generation_number = 0;
			Fitness_t achieved_fitness = field.begin()->first;

			// Цикл скрещивания
			// пока (генетический алгоритм ⇒ поколение < максимум поколений) и (необходимая приспособленность < достигнутой приспособленности)
			// Выход в середине тела цикла
			RandomFloatGenerator& generator = RandomFloatGenerator::get_instance();

			while ((!is_generation_based || generation_number < max_generations) && (required_fitness < achieved_fitness)) {
				//Скрещивание
				for (auto iter1 = field.begin(); iter1 != field.end(); ++iter1) {
					for (auto iter2 = (symmetric_crossover ? std::next(iter1) : field.begin()); iter2 != field.end(); ++iter2) {
						newbies.push_back(std::move(crossover(iter1->second, iter2->second)));
						if (generator.get_random_float() <= mutation_probability)
							mutate(newbies.back());
					}
				}
				update_field();

				if (survivors != unlimited && field.size() > survivors && generation_number % death_rate == 0)
					field.erase(std::next(field.begin(), survivors), field.end());

				achieved_fitness = field.begin()->first; //Лучший
				++generation_number;
			}

			Population<DNA> result(
				std::vector<DNA>(field.size()),
				generation0.generation() + generation_number);
			std::transform(
				field.begin(),
				field.end(),
				result->begin(),
				[](std::pair<Fitness_t, DNA> p) -> DNA { return std::move(p.second); });
			return result;
		}
	};
}

#endif //__GENETICS_H
