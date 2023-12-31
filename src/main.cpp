#include "util.hpp"

#include <iostream>
#include <conio.h>
#include <string>

Matrix matrix;
Algorithms algo;

void mainMenu() {
	std::cout << "\n==== MENU GLOWNE ====\n";
	std::cout << "1.Wczytaj z pliku - " << algo.getBenchmarkFile() << "\n";
	std::cout << "2.Ustaw czas wykonywania [s] (kryterium stopu) - "
		<< algo.getMaxExecutionTime().count() << "\n";
	std::cout << "3.Ustaw rozmiar populacji poczatkowej - "
		<< algo.getStartingPopulationSize() << "\n";
	std::cout << "4.Ustaw wspolczynnik mutacji - "
		<< algo.getMutationConstant() << "\n";
	std::cout << "5.Ustaw wspolczynnik krzyzowania - "
		<< algo.getCrossoverConstant() << "\n";
	std::cout << "6.Wybierz metode krzyzowania - " <<
		(algo.getCurrentCrossoverType() ? "EAX" : "Order Crossover") << "\n";
	std::cout << "7.Wybierz metode mutacji - ";
	switch (algo.getCurrentMutationType()) {
	case 0:
		std::cout << "Odwrotna kolejnosc\n";
		break;
	case 1:
		std::cout << "Zamiana miejsc\n";
		break;
	case 2:
		std::cout << "Wstaw w miejsce\n";
		break;
	default:
		break;
	}
	std::cout << "8.Uruchom algorytm\n";
	std::cout << "0.Wyjdz\n";
	std::cout << "Podaj opcje:";
}

void mutationMenu() {
	clear();

	char option;
	do {
		std::cout << "\n==== MUTACJA ====\n";
		std::cout << "1.Odwrotna kolejnosc\n";
		std::cout << "2.Zamiana miejsc\n";
		std::cout << "3.Wstaw w miejsce\n";
		// std::cout << "4.Wstaw podsciezke - tylko wyzarzanie\n";
		std::cout << "0.Powrot\n";
		std::cout << "Podaj opcje:";
		option = _getche();
		std::cout << std::endl;

		switch (option) {
		case '1':
			algo.setMutationType(INVERSE);
			return;

		case '2':
			algo.setMutationType(SWAP);
			return;

		case '3':
			algo.setMutationType(INSERT);
			return;
			/*
		case '4':
			algo.setMutationType(INSERT_SUB);
			return;*/
		}
	} while (option != '0');
}

void crossoverMenu() {
	clear();

	char option;
	do {
		std::cout << "\n==== KRZYZOWANIE ====\n";
		std::cout << "1.Order Crossover (OX)\n";
		std::cout << "2.(EAX)\n";
		std::cout << "0.Powrot\n";
		std::cout << "Podaj opcje:";
		option = _getche();
		std::cout << std::endl;

		switch (option) {
		case '1':
			algo.setCrossoverType(0);
			return;

		case '2':
			algo.setCrossoverType(1);
			return;
		}
	} while (option != '0');
}

int main() {
	char option;
	std::string fileName;
	int value;
	double valueD;

	algo.initRandom();

	do {
		mainMenu();
		option = _getche();
		std::cout << std::endl;

		switch (option) {
		case '1':
			std::cout << " Podaj nazwe zbioru:";
			std::cin >> fileName;
			if (matrix.loadFromFile(fileName)) {
				matrix.display();
				algo.setBenchmarkFile(fileName);
			}
			else algo.setBenchmarkFile("BRAK");
			
			break;

		case '2':
			std::cout << "Podaj dlugosc wykonania (s):";
			std::cin >> value;
			algo.setStopCriterium(value);
			clear();
			break;

		case '3':
			std::cout << "Podaj rozmiar populacji startowej: ";
			std::cin >> value;
			algo.setStartingPopulationSize(value);
			clear();
			break;

		case '4':
			std::cout << "Podaj wspolczynnik mutacji:";
			std::cin >> valueD;
			algo.setMutationConstant(valueD);
			clear();
			break;

		case '5':
			std::cout << "Podaj wspolczynnik krzyzowania:";
			std::cin >> valueD;
			algo.setCrossoverConstant(valueD);
			clear();
			break;

		case '6':
			crossoverMenu();
			break;
		case '7':
			mutationMenu();
			break;
		case '8':
			algo.geneticAlgorithm((Matrix*)&matrix);
			clear();
			algo.displayResults();
			break;
		}
	} while (option != '0');

	return 0;
}