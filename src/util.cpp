#include "util.hpp"
#include <iostream>
#include <string>
#include <numeric>

void Matrix::loadFromFile(std::string fileName) {
	int n = fileName.find(".txt");

	if (n != std::string::npos
		&& fileName.length() - n == 4) {
		loadFromTXT(fileName);
		return;
	}

	n = fileName.find(".atsp");

	if (n != std::string::npos
		&& fileName.length() - n == 5) {
		loadFromATSP(fileName);
	}
	else
		std::cout << "Nie mozna otworzyc pliku!\n";
}

void Matrix::loadFromTXT(std::string fileName) {
	std::fstream file;
	file.open(fileName, std::ios::in);

	if (file.good()) {
		// je?li istnieje poprzednia matryca, czy?cimy
		if (mat.size()) mat.clear();

		int dimension = 0, cross = 0, valueInMatrix = 0;
		std::string stringTemp;

		file >> dimension;

		this->size = dimension;
		this->mat.reserve(dimension);
		std::vector<std::vector<int>>::iterator matIter = mat.begin();

		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				file >> valueInMatrix;

				// pierwszy wiersz wymaga specjalnego traktowania - musimy zrobi? push_back vectorów
				if (mat.size() < size) {
					std::vector<int> tempVec;
					tempVec.push_back(!mat.size() ? 0 : valueInMatrix);
					mat.push_back(tempVec);

					matIter = mat.begin();
				}
				// reszta normalny pushback
				else {
					if (i == j) (*matIter).push_back(0);
					else (*matIter).push_back(valueInMatrix);

					matIter++;
					if (matIter == mat.end()) matIter = mat.begin();
				}
			}
		}

		file.close();
	}
	else std::cout << "Plik nie zostal otworzony!\n";
}

void Matrix::loadFromATSP(std::string fileName) {
	std::fstream file;
	file.open(fileName, std::ios::in);

	if (file.good()) {
		// je?li istnieje poprzednia matryca, czy?cimy
		if (mat.size()) mat.clear();

		int dimension = 0, cross = 0, valueInMatrix = 0;
		std::string stringTemp;

		// Znalezienie wielko?ci matrycy
		do {
			file >> stringTemp;
		} while (stringTemp != "DIMENSION:");

		file >> dimension;

		this->size = dimension;
		this->mat.reserve(dimension);
		std::vector<std::vector<int>>::iterator matIter = mat.begin();

		do {
			file >> stringTemp;
		} while (stringTemp != "EDGE_WEIGHT_SECTION");

		// Dotarcie do warto?ci zapisanych w matrycy, podgl?d warto?ci zapisanych na przek?tnych (?)
		file >> stringTemp;
		cross = std::stoi(stringTemp);

		while (stringTemp != "EOF") {
			valueInMatrix = std::stoi(stringTemp);

			// pierwszy wiersz wymaga specjalnego traktowania - musimy zrobi? push_back vectorów
			if (mat.size() < size) {
				std::vector<int> tempVec;
				tempVec.push_back(!mat.size() ? 0 : valueInMatrix);
				mat.push_back(tempVec);

				matIter = mat.begin();
			}
			// reszta normalny pushback
			else {
				if (valueInMatrix == cross || valueInMatrix == 0) (*matIter).push_back(0);
				else (*matIter).push_back(valueInMatrix);

				matIter++;
				if (matIter == mat.end()) matIter = mat.begin();
			}

			file >> stringTemp;
		}
		file.close();
	}
	else std::cout << "Plik nie zostal otworzony!\n";
}

void Matrix::display() {
	std::vector<std::vector<int>>::iterator matIter;
	std::vector<int>::iterator matIterInner;

	for (int i = 0; i < mat.size(); i++) {
		for (matIter = mat.begin(); matIter != mat.end(); matIter++) {
			matIterInner = (*matIter).begin();
			std::advance(matIterInner, i);

			int spacesToAdd = 2;
			if (*matIterInner) {
				spacesToAdd -= (int)log10((double)(*matIterInner));
			}

			std::cout << (*matIterInner) << " ";
			while (spacesToAdd > 0) {
				std::cout << " ";
				spacesToAdd--;
			}
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

void Algorithms::displayResults() {
	std::cout << "\nDlugosc sciezki: " << pathLength << "\n";
	std::cout << "Kolejnosc wierzcholkow:\n0 ";
	for (auto a : this->vertexOrder) std::cout << a << " ";
	std::cout << "0\nCzas trwania algorytmu: " << runningTime.count() << "s\n";
}

void Algorithms::setStopCriterium(int value) {
	maxExecutionTime = std::chrono::seconds(value);
}

void Algorithms::setMutationConstant(double value) {
	mutationConstant = value;
}

/*
void Algorithms::setCrossoverConstant(double value) {
	crossoverConstant = value;
}*/

void Algorithms::setStartingPopulationSize(int value) {
	if (value < 3) {
		startingPopulationSize = 3;
		std::cout << "Wartosc musi byc minimum 3!\n";
	}
	else 
		startingPopulationSize = value;
}

void Algorithms::setMutationType(MutationType type) {
	currentMutationType = type;
}

void Algorithms::initRandom() {
	gen.seed(rd());
}

void Algorithms::benchmark(Matrix* matrix) {
	int maxIter = 10;
	int time = 120;

	matrix->loadFromFile("ftv47.atsp");
	setStopCriterium(time);
	benchmarkFile = "TS_47.txt";
	for (int j = 0; j < 3; j++) {
		currentMutationType = (MutationType)j;
		for (int i = 0; i < maxIter; i++) {
			//tabuSearchBenchmark(matrix, i);
		}
	}
	benchmarkFile = "SA_47.txt";
	//coolingConstant = 0.75;
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < maxIter; i++) {
			//simulatedAnnealingBenchmark(matrix, i);
		}
		//coolingConstant -= 0.25;
	}

	matrix->loadFromFile("ftv170.atsp");
	setStopCriterium(time * 2);
	benchmarkFile = "TS_170.txt";
	for (int j = 0; j < 3; j++) {
		currentMutationType = (MutationType)j;
		for (int i = 0; i < maxIter; i++) {
			//tabuSearchBenchmark(matrix, i);
		}
	}
	benchmarkFile = "SA_170.txt";
	//coolingConstant = 0.75;
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < maxIter; i++) {
		//	simulatedAnnealingBenchmark(matrix, i);
		}
		//coolingConstant -= 0.25;
	}

	matrix->loadFromFile("rbg403.atsp");
	setStopCriterium(time * 3);
	benchmarkFile = "TS_403.txt";
	for (int j = 0; j < 3; j++) {
		currentMutationType = (MutationType)j;
		for (int i = 0; i < maxIter; i++) {
			//tabuSearchBenchmark(matrix, i);
		}
	}
	benchmarkFile = "SA_403.txt";
	//coolingConstant = 0.75;
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < maxIter; i++) {
			//simulatedAnnealingBenchmark(matrix, i);
		}
		//coolingConstant -= 0.25;
	}
}

std::tuple<int, int> Algorithms::generateRandomTwoPositions(int lowerBound, int higherBound, bool correctOrder) {
	// generate two positions
	std::uniform_int_distribution<> distribution(lowerBound, higherBound);

	int indexOne = distribution(gen), indexTwo = distribution(gen);

	while (indexOne == indexTwo)
		indexTwo = distribution(gen);

	if (!correctOrder)
		return std::make_tuple(indexOne, indexTwo);

	if (indexTwo < indexOne) {
		int temp = indexOne;
		indexOne = indexTwo;
		indexTwo = temp;
	}

	return std::make_tuple(indexOne, indexTwo);
}

std::tuple<std::vector<short>, int> Algorithms::generateInitialSolution(Matrix* matrix) {
	// Greedy method
	std::vector<short> possibleVertices, returnVector;
	int matrixSize = matrix->size, returnLength = 0;
	returnVector.reserve(matrixSize);
	for (int i = 1; i < matrixSize; i++)
		possibleVertices.push_back(i);

	int currentVertex = 0;

	while (possibleVertices.size()) {
		int value = INT_MAX, lowestIndex = 0;
		// znajdz najkrotsza mozliwa sciezke
		for (int i = 1; i < matrixSize; i++) {
			if (matrix->mat[i][currentVertex] < value
				&& (std::find(possibleVertices.begin(), possibleVertices.end(), i) != std::end(possibleVertices))) {
				value = matrix->mat[i][currentVertex];
				lowestIndex = i;
			}
		}

		// dodaj do generowanego rozwiazania
		returnVector.push_back(lowestIndex);
		returnLength += value;

		// zmien "obecny" wierzcholek
		currentVertex = lowestIndex;

		// usun z mozliwych do wybrania wierzcholkow
		possibleVertices.erase(
			std::find(possibleVertices.begin(), possibleVertices.end(), lowestIndex)
		);
	}

	returnLength += matrix->mat[0][currentVertex];

	return std::make_tuple(returnVector, returnLength);
}

std::tuple<std::vector<short>, int> Algorithms::generateNewSolution(Matrix* matrix, int notAllowedSecondary) {
	// Greedy method
	std::vector<short> possibleVertices, returnVector;
	int matrixSize = matrix->size, returnLength = INT_MAX;
	returnVector.reserve(matrixSize);
	for (int i = 1; i < matrixSize; i++)
		possibleVertices.push_back(i);

	std::uniform_int_distribution<> u(1, matrixSize - 1);

	// zmien "obecny" wierzcholek
	int currentVertex = notAllowedSecondary;
	do {
		currentVertex = u(rd);
	} while (currentVertex == notAllowedSecondary);
	returnVector.push_back(currentVertex);
	returnLength = matrix->mat[currentVertex][0];

	// usun z mozliwych do wybrania wierzcholkow
	possibleVertices.erase(
		std::find(possibleVertices.begin(), possibleVertices.end(), currentVertex)
	);

	while (possibleVertices.size()) {
		int value = INT_MAX, lowestIndex = 0;
		// znajdz najkrotsza mozliwa sciezke
		for (int i = 1; i < matrixSize; i++) {
			if (matrix->mat[i][currentVertex] < value
				&& (std::find(possibleVertices.begin(), possibleVertices.end(), i) != std::end(possibleVertices))) {
				value = matrix->mat[i][currentVertex];
				lowestIndex = i;
			}
		}

		// dodaj do generowanego rozwiazania
		returnVector.push_back(lowestIndex);
		returnLength += value;

		// zmien "obecny" wierzcholek
		currentVertex = lowestIndex;

		// usun z mozliwych do wybrania wierzcholkow
		possibleVertices.erase(
			std::find(possibleVertices.begin(), possibleVertices.end(), lowestIndex)
		);
	}

	returnLength += matrix->mat[0][currentVertex];

	return std::make_tuple(returnVector, returnLength);
}

std::tuple<std::vector<short>, int> Algorithms::generateNewSolutionV(Matrix* matrix, const short startVertex) {
	// Greedy method
	if (!startVertex)
		return generateInitialSolution(matrix);

	int matrixSize = matrix->size;
	std::vector<short> localPossibleVertices(matrixSize - 1), returnVector;
	int returnLength = INT_MAX;
	returnVector.reserve(matrixSize);

	// usun z mozliwych do wybrania wierzcholkow
	std::iota(localPossibleVertices.begin(), localPossibleVertices.end(), 1);
	localPossibleVertices.erase(
		std::find(localPossibleVertices.begin(), localPossibleVertices.end(), startVertex)
	);
	returnVector.push_back(startVertex);
	returnLength = matrix->mat[startVertex][0];

	// zmien "obecny" wierzcholek
	int currentVertex = 0;
	while (localPossibleVertices.size()) {
		int value = INT_MAX, lowestIndex = 1;
		// znajdz najkrotsza mozliwa sciezke
		for (int i = 1; i < matrixSize; i++) {
			if (matrix->mat[i][currentVertex] < value
				&& (std::find(localPossibleVertices.begin(), localPossibleVertices.end(), i) != std::end(localPossibleVertices))) {
				value = matrix->mat[i][currentVertex];
				lowestIndex = i;
			}
		}

		// dodaj do generowanego rozwiazania
		returnVector.push_back(lowestIndex);
		returnLength += value;

		// zmien "obecny" wierzcholek
		currentVertex = lowestIndex;

		// usun z mozliwych do wybrania wierzcholkow
		localPossibleVertices.erase(
			std::find(localPossibleVertices.begin(), localPossibleVertices.end(), lowestIndex)
		);
	}

	returnLength += matrix->mat[0][currentVertex];

	return std::make_tuple(returnVector, returnLength);
}

std::vector<short> Algorithms::inverse(std::vector<short>* currentOrder, int firstPosition, int secondPosition) {
	std::vector<short> returnVector(*currentOrder);
	int vectorSize = currentOrder->size();

	std::tuple<int, int> t;
	// generate two positions
	if (!firstPosition && !secondPosition)
		t = generateRandomTwoPositions(0, vectorSize - 1);
	else t = std::make_tuple(firstPosition, secondPosition);

	std::vector<short>::iterator low = returnVector.begin(), high = returnVector.begin();
	std::advance(low, std::get<0>(t));
	std::advance(high, std::get<1>(t) + 1);
	std::reverse(low, high);

	return returnVector;
}

std::vector<short> Algorithms::swap(std::vector<short>* currentOrder, int firstPosition, int secondPosition) {
	std::vector<short> returnVector(*currentOrder);
	int vectorSize = currentOrder->size();

	std::tuple<int, int> t;
	// generate two positions
	if (!firstPosition && !secondPosition)
		t = generateRandomTwoPositions(0, vectorSize - 1);
	else t = std::make_tuple(firstPosition, secondPosition);

	// swap positions
	std::vector<short>::iterator low = returnVector.begin(), high = returnVector.begin();

	std::advance(low, std::get<0>(t));
	std::advance(high, std::get<1>(t));
	std::swap(low, high);

	return returnVector;
}

std::vector<short> Algorithms::insert(std::vector<short>* currentOrder, int firstPosition, int secondPosition) {
	std::vector<short> returnVector(*currentOrder);
	int vectorSize = currentOrder->size();

	std::tuple<int, int> t;
	// generate two positions
	if (!firstPosition && !secondPosition)
		t = generateRandomTwoPositions(0, vectorSize - 1, false);
	else t = std::make_tuple(firstPosition, secondPosition);

	// insert
	if (std::get<0>(t) > std::get<1>(t)) {
		std::rotate(returnVector.rend() - std::get<0>(t) - 1,
			returnVector.rend() - std::get<0>(t),
			returnVector.rend() - std::get<1>(t));
	}
	else {
		std::rotate(returnVector.begin() + std::get<0>(t),
			returnVector.begin() + std::get<0>(t) + 1,
			returnVector.begin() + std::get<1>(t) + 1);
	}

	return returnVector;
}

std::vector<short> Algorithms::insertSub(std::vector<short>* currentOrder) {
	std::vector<short> returnVector(*currentOrder);
	int vectorSize = currentOrder->size();
	std::vector<short>::iterator beginCurrent = currentOrder->begin();

	// generate two positions
	std::tuple<int, int> t = generateRandomTwoPositions(0, vectorSize - 1);
	std::uniform_int_distribution<> distribution(0, vectorSize - 1);
	int indexThree = distribution(gen);

	// insert group
	if (indexThree > std::get<1>(t)) {
		if (std::get<0>(t) == 0) {
			returnVector.insert(returnVector.begin(), beginCurrent + std::get<1>(t) + 1, beginCurrent + indexThree + 1);
			returnVector.insert(returnVector.begin() + (indexThree - std::get<1>(t)), beginCurrent + std::get<0>(t), beginCurrent + std::get<1>(t) + 1);
		}
		else {
			returnVector.insert(returnVector.begin(), beginCurrent, beginCurrent + std::get<0>(t));
			returnVector.insert(returnVector.end(), beginCurrent + std::get<1>(t) + 1, beginCurrent + indexThree + 1);
			returnVector.insert(returnVector.end(), beginCurrent + std::get<0>(t), beginCurrent + std::get<1>(t) + 1);
		}
	}
	else if (indexThree < std::get<0>(t)) {
		returnVector.insert(returnVector.begin(), beginCurrent, beginCurrent + indexThree);
		returnVector.insert(returnVector.end(), beginCurrent + std::get<0>(t), beginCurrent + std::get<1>(t) + 1);
		returnVector.insert(returnVector.end(), beginCurrent + indexThree, beginCurrent + std::get<0>(t));
		returnVector.insert(returnVector.end(), beginCurrent + std::get<1>(t) + 1, currentOrder->end());
	}
	else {
		returnVector.insert(returnVector.begin(), beginCurrent, beginCurrent + std::get<0>(t));
		returnVector.insert(returnVector.end(), beginCurrent + std::get<1>(t) + 1, beginCurrent + std::min(std::get<1>(t) + (std::get<1>(t) - std::get<0>(t)), (int)currentOrder->size()));
		returnVector.insert(returnVector.end(), beginCurrent + std::get<0>(t), beginCurrent + std::get<1>(t) + 1);
	}

	return returnVector;
}

QueueData Algorithms::getNewOrder(std::vector<short>* currentOrder, int anchorOne, int anchorTwo, std::vector<std::vector<int>>* matrix) {
	QueueData tempData;
	short previousVertex = 0;

	switch (currentMutationType)
	{
	case INVERSE:
		tempData.pathOrder = inverse(currentOrder, anchorOne, anchorTwo);
		break;
	case SWAP:
		tempData.pathOrder = swap(currentOrder, anchorOne, anchorTwo);
		break;
	case INSERT:
		tempData.pathOrder = insert(currentOrder, anchorOne, anchorTwo);
		break;
	default:
		break;
	}

	tempData.pathLength = 0;

	for (int k = 0; k < tempData.pathOrder.size(); k++) {
		tempData.pathLength += (*matrix)[tempData.pathOrder[k]][previousVertex];
		previousVertex = tempData.pathOrder[k];
	}

	tempData.pathLength += (*matrix)[0][previousVertex];
	tempData.anchorOne = tempData.pathOrder[anchorOne];
	tempData.anchorTwo = tempData.pathOrder[anchorTwo];

	return tempData;
}

int Algorithms::calculateCandidate(std::vector<short>* candidateOrder, Matrix* matrix) {
	int pathLength = 0, previousVector = 0;
	std::vector<std::vector<int>>* matrixStart = &(matrix->mat);
	std::vector<short>::iterator ending = candidateOrder->end(), beginning = candidateOrder->begin();

	for (auto iter = beginning; iter != ending; iter++) {
		pathLength += (*matrixStart)[*iter][previousVector];
		previousVector = *iter;
	}

	pathLength += (*matrixStart)[0][previousVector];

	return pathLength;
}

void clear() {
#if defined _WIN32
	system("cls");
#elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
	system("clear");
#endif
}