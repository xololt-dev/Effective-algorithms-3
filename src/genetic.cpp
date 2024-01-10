#include "util.hpp"

#include <numeric>
#include <thread>

// https://dl.acm.org/doi/pdf/10.1145/3071178.3071305
// populacja - ocena - selekcja - pula rodzicielska - operacje genetyczne - subpopulacja - ocena - sukcesja - ...

void Algorithms::geneticAlgorithm(Matrix* matrix, int benchmarkIteration) {
	if (!matrix->size) {
		std::cout << "Nie wczytano matrycy!\n";
		return;
	}
	
	if (!benchmarkIteration) {
		if (currentCrossoverType)
			geneticEAX(matrix);
		else
			geneticOX(matrix);
	}
	else {
		if (currentCrossoverType)
			geneticEAXBench(matrix, benchmarkIteration);
		else
			geneticOXBench(matrix, benchmarkIteration);
	}	
}

void Algorithms::geneticOX(Matrix* matrix) {
	int iterations = 0;
	// Priority queue func
	struct {
		bool operator()(const PathData x, const PathData y) const
			{ return x.pathLength < y.pathLength; };
	} compareS;

	pathLength = INT_MAX;

	// Starting population generation
	std::vector<PathData> population = generateStartingPopulation(matrix);
	std::sort(population.begin(), population.end(), compareS);

	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
	std::vector<double> lowerBounds;
	lowerBounds.reserve(population.size());
	std::vector<std::tuple<int, int>> parents;
	parents.reserve(population.size());
	std::vector<PathData> childrenData;
	childrenData.reserve(population.size());

	while (std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - start) < maxExecutionTime) {

		// Update best solution if applicable
		if (population[0].pathLength < pathLength) {
			pathLength = population[0].pathLength;
			vertexOrder = population[0].pathOrder;
			runningTime = std::chrono::duration_cast<std::chrono::microseconds>
				(std::chrono::steady_clock::now() - start);

			for (short& s : population[0].pathOrder) {
				std::cout << s << " ";
			}
			std::cout << "\n";
		}

		// Get lower bounds
		lowerBounds = getVertexLowerBounds(population.size());

		// Generate pairs
		parents = generateParents(&lowerBounds, population.size());

		// Dump gen. to vector

		// Genetic operations
		std::uniform_real_distribution<> distribution(0.0, 1.0);
		for (std::tuple<int, int>& t : parents) {
			if (distribution(gen) > crossoverConstant)
				continue;

			// Crossover
			PathData child = generateChildOX(matrix, &population[std::get<0>(t)].pathOrder,
				&population[std::get<1>(t)].pathOrder);
			
			child.pathLength = calculateCandidate(&child.pathOrder, matrix);
			// Mutate
			if (distribution(gen) <= mutationConstant) {
				if (child.pathLength < pathLength)
					childrenData.push_back(PathData(child));

				std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
				child = getNewOrder(&child.pathOrder, std::get<0>(iT),
					std::get<1>(iT), &(matrix->mat));
			}

			childrenData.emplace_back(std::move(child));
		}

		// Combine to get new generation
		population = getNewRandomGeneration(matrix, &population, &childrenData);

		std::sort(population.begin(), population.end(), compareS);

		iterations++;
		lowerBounds.clear();
		parents.clear();
		childrenData.clear();
	}

	std::cout << iterations << "\n";
}

void Algorithms::geneticEAX(Matrix* matrix) {
	int iterations = 0;
	// Priority queue func
	struct {
		bool operator()(const PathData x, const PathData y) const
		{
			return x.pathLength < y.pathLength;
		};
	} compareS;

	pathLength = INT_MAX;

	// Starting population generation
	std::vector<PathData> population = generateStartingPopulation(matrix);

	std::sort(population.begin(), population.end(), compareS);

	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
	std::vector<double> lowerBounds;
	lowerBounds.reserve(population.size());
	std::vector<std::tuple<int, int>> parents;
	parents.reserve(population.size());
	std::vector<PathData> childrenData;
	childrenData.reserve(population.size());

	while (std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - start) < maxExecutionTime) {

		// Update best solution if applicable
		if (population[0].pathLength < pathLength) {
			pathLength = population[0].pathLength;
			vertexOrder = population[0].pathOrder;
			runningTime = std::chrono::duration_cast<std::chrono::microseconds>
				(std::chrono::steady_clock::now() - start);

			for (short& s : population[0].pathOrder) {
				std::cout << s << " ";
			}
			std::cout << "\n";
		}

		// Get lower bounds
		lowerBounds = getVertexLowerBounds(population.size());

		// Generate pairs
		parents = generateParents(&lowerBounds, population.size());

		// Generate children

		// Genetic operations
		std::uniform_real_distribution<> distribution(0.0, 1.0);
		total = 0,
		repeat = 0;
		for (std::tuple<int, int>& t : parents) {
			if (distribution(gen) > crossoverConstant)
				continue;
			total++;
			// Crossover
			PathData child = generateChildEAX(matrix, &population[std::get<0>(t)].pathOrder,
				&population[std::get<1>(t)].pathOrder);

			// Mutate
			if (distribution(gen) <= mutationConstant) {
				std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
				child = getNewOrder(&child.pathOrder, std::get<0>(iT),
					std::get<1>(iT), &(matrix->mat));
			}
			else child.pathLength = calculateCandidate(&child.pathOrder, matrix);

			childrenData.emplace_back(std::move(child));
		}

		// Combine to get new generation
		population = getNewRandomGeneration(matrix, &population, &childrenData);

		std::sort(population.begin(), population.end(), compareS);
		iterations++;
		lowerBounds.clear();
		parents.clear();
		childrenData.clear();
	}
	std::cout << iterations << "\n";
}

void Algorithms::geneticOXBench(Matrix* matrix, int iteration) {
	// Priority queue func
	struct {
		bool operator()(const PathData x, const PathData y) const
		{
			return x.pathLength < y.pathLength;
		};
	} compareS;

	pathLength = INT_MAX;

	// Starting population generation
	std::vector<PathData> population = generateStartingPopulation(matrix);
	std::sort(population.begin(), population.end(), compareS);

	std::fstream file;
	file.open(benchmarkFile, std::fstream::app | std::fstream::out);

	if (file.good()) {
		std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

		while (std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::steady_clock::now() - start) < maxExecutionTime) {

			// Update best solution if applicable
			if (population[0].pathLength < pathLength) {
				pathLength = population[0].pathLength;
				vertexOrder = population[0].pathOrder;
				runningTime = std::chrono::duration_cast<std::chrono::microseconds>
					(std::chrono::steady_clock::now() - start);

				file << iteration << ";" << currentCrossoverType << ";" << currentMutationType <<
					";" << startingPopulationSize << ";" << startingPopulationRandom << ";" <<
					mutationConstant << ";" << crossoverConstant << ";" << 
					std::chrono::duration_cast<std::chrono::microseconds>(runningTime).count() << ";" << pathLength << ";\n";
			}

			// Get lower bounds
			std::vector<double> lowerBounds = getVertexLowerBounds(population.size());

			// Generate pairs
			std::vector<std::tuple<int, int>> parents =
				generateParents(&lowerBounds, population.size());

			// Dump gen. to vector
			std::vector<PathData> childrenData;
			childrenData.reserve(population.size());

			// Genetic operations
			std::uniform_real_distribution<> distribution(0.0, 1.0);
			for (std::tuple<int, int> t : parents) {
				if (distribution(gen) > crossoverConstant)
					continue;

				// Crossover
				PathData child = generateChildOX(matrix, &population[std::get<0>(t)].pathOrder,
					&population[std::get<1>(t)].pathOrder);

				child.pathLength = calculateCandidate(&child.pathOrder, matrix);
				// Mutate
				if (distribution(gen) <= mutationConstant) {
					if (child.pathLength < pathLength)
						childrenData.push_back(PathData(child));

					std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
					child = getNewOrder(&child.pathOrder, std::get<0>(iT),
						std::get<1>(iT), &(matrix->mat));
				}

				childrenData.emplace_back(std::move(child));
			}

			// Combine to get new generation
			population = getNewRandomGeneration(matrix, &population, &childrenData);

			std::sort(population.begin(), population.end(), compareS);
		}

		file << iteration << ";" << currentCrossoverType << ";" << currentMutationType <<
			";" << startingPopulationSize << ";" << startingPopulationRandom << ";" <<
			mutationConstant << ";" << crossoverConstant << ";" << 
			std::chrono::duration_cast<std::chrono::microseconds>(maxExecutionTime).count() << ";" << pathLength << ";0-";
		for (auto& a : vertexOrder) file << a << "-";
		file << "0\n";

		file.close();
	}
	else std::cout << "File not opened!\n";
}

void Algorithms::geneticEAXBench(Matrix* matrix, int iteration) {
	// Priority queue func
	struct {
		bool operator()(const PathData x, const PathData y) const
		{
			return x.pathLength < y.pathLength;
		};
	} compareS;

	pathLength = INT_MAX;

	// Starting population generation
	std::vector<PathData> population = generateStartingPopulation(matrix);

	std::sort(population.begin(), population.end(), compareS);

	std::fstream file;
	file.open(benchmarkFile, std::fstream::app | std::fstream::out);

	if (file.good()) {
		std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

		while (std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::steady_clock::now() - start) < maxExecutionTime) {

			// Update best solution if applicable
			if (population[0].pathLength < pathLength) {
				pathLength = population[0].pathLength;
				vertexOrder = population[0].pathOrder;
				runningTime = std::chrono::duration_cast<std::chrono::microseconds>
					(std::chrono::steady_clock::now() - start);

				file << iteration << ";" << currentCrossoverType << ";" << currentMutationType <<
					";" << startingPopulationSize << ";" << startingPopulationRandom << ";" <<
					mutationConstant << ";" << crossoverConstant << ";" << 
					std::chrono::duration_cast<std::chrono::microseconds>(runningTime).count() << ";" << pathLength << ";\n";
			}

			// Get lower bounds
			std::vector<double> lowerBounds = getVertexLowerBounds(population.size());

			// Generate pairs
			std::vector<std::tuple<int, int>> parents =
				generateParents(&lowerBounds, population.size());

			// Generate children
			std::vector<PathData> childrenData;
			childrenData.reserve(population.size());

			// Genetic operations
			std::uniform_real_distribution<> distribution(0.0, 1.0);
			for (std::tuple<int, int>& t : parents) {
				if (distribution(gen) > crossoverConstant)
					continue;

				// Crossover
				PathData child = generateChildEAX(matrix, &population[std::get<0>(t)].pathOrder,
					&population[std::get<1>(t)].pathOrder);

				// Mutate
				if (distribution(gen) <= mutationConstant) {
					std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
					child = getNewOrder(&child.pathOrder, std::get<0>(iT),
						std::get<1>(iT), &(matrix->mat));
				}
				else child.pathLength = calculateCandidate(&child.pathOrder, matrix);

				childrenData.emplace_back(std::move(child));
			}

			// Combine to get new generation
			population = getNewRandomGeneration(matrix, &population, &childrenData);

			std::sort(population.begin(), population.end(), compareS);
		}

		file << iteration << ";" << currentCrossoverType << ";" << currentMutationType <<
			";" << startingPopulationSize << ";" << startingPopulationRandom << ";" <<
			mutationConstant << ";" << crossoverConstant << ";" << 
			std::chrono::duration_cast<std::chrono::microseconds>(maxExecutionTime).count() << ";" << pathLength << ";0-";
		for (auto& a : vertexOrder) file << a << "-";
		file << "0\n";

		file.close();
	}
	else std::cout << "File not opened!\n";	
}

std::vector<PathData> Algorithms::generateStartingPopulation(Matrix* matrix) {
	std::vector<PathData> returnVec;
	returnVec.reserve(startingPopulationSize);
	std::vector<short> leftVertices(matrix->size - 1);
	std::iota(leftVertices.begin(), leftVertices.end(), 1);

	if (startingPopulationRandom == true) {
		generateRandomStartingPopulation(matrix, &leftVertices, &returnVec);
	}
	else {
		generateGreedyStartingPopulation(matrix, &leftVertices, &returnVec);

		if (startingPopulationSize >= matrix->size) {
			std::vector<short> leftVertices(matrix->size - 1);
			std::iota(leftVertices.begin(), leftVertices.end(), 1);

			generateRandomStartingPopulation(matrix, &leftVertices, &returnVec);
		}
	}
	
	return returnVec;
}

void Algorithms::generateRandomStartingPopulation(Matrix* matrix, std::vector<short>* leftVertices, std::vector<PathData>* returnVec) {
	while (returnVec->size() < startingPopulationSize) {
		PathData tempData;
		std::shuffle(leftVertices->begin(), leftVertices->end(), gen);

		short previousVertex = 0;
		tempData.pathLength = 0;
		tempData.pathOrder = *leftVertices;

		for (int k = 0; k < tempData.pathOrder.size(); k++) {
			tempData.pathLength += matrix->mat[tempData.pathOrder[k]][previousVertex];
			previousVertex = tempData.pathOrder[k];
		}

		tempData.pathLength += matrix->mat[0][previousVertex];

		returnVec->emplace_back(std::move(tempData));
	}
}

void Algorithms::generateGreedyStartingPopulation(Matrix* matrix, std::vector<short>* leftVertices, std::vector<PathData>* returnVec) {
	int cap = 0;
	PathData tempData;

	std::tuple<std::vector<short>, int> t = generateInitialSolution(matrix);
	short firstOne = std::get<0>(t)[0];
	leftVertices->erase(
		std::find(leftVertices->begin(), leftVertices->end(), firstOne)
	);

	tempData.pathOrder = std::get<0>(t);
	tempData.pathLength = std::get<1>(t);
	returnVec->emplace_back(std::move(tempData));
	
	if (startingPopulationSize < matrix->size)
		cap = startingPopulationSize;
	else cap = matrix->size - 1;

	for (int i = 1; i < cap; i++) {
		PathData tempDataTwo;

		std::shuffle(leftVertices->begin(), leftVertices->end(), gen);

		t = generateNewSolutionV(matrix, (*leftVertices)[0]);

		firstOne = std::get<0>(t)[0];
		leftVertices->erase(
			std::find(leftVertices->begin(), leftVertices->end(), firstOne)
		);

		tempDataTwo.pathOrder = std::get<0>(t);
		tempDataTwo.pathLength = std::get<1>(t);
		returnVec->emplace_back(std::move(tempDataTwo));
	}
}

std::vector<double> Algorithms::getVertexLowerBounds(int vectorSize) {
	std::vector<double> vertexLowerBound(vectorSize, 0.0);
	
	for (int i = 1; i < vertexLowerBound.size(); i++) {
		double current = (1.0 - (double) i / (double) vectorSize);
		vertexLowerBound[i] = current + vertexLowerBound[i - 1];
	}

	double sum = vertexLowerBound[vertexLowerBound.size() - 1];

	for (int i = 1; i < vertexLowerBound.size(); i++) {
		vertexLowerBound[i] /= sum;
	}

	return vertexLowerBound;
}

std::vector<std::tuple<int, int>> Algorithms::generateParents(std::vector<double>* boundsVector, int vectorSize) {
	std::vector<std::tuple<int, int>> pairsVector;
	pairsVector.reserve(boundsVector->size());
	// Maybe saved "wrong" generations to a vector for later usage?
	std::vector<int> generatedRandoms;
	generatedRandoms.reserve(4);

	// Generate number(s) [0, 1]
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	int childrenGenerated = 0;

	while (childrenGenerated < vectorSize) {
		double generated;
		int firstCandidate = -1;
		int secondCandidate = 0;
		std::vector<double>::iterator boundIter = boundsVector->begin();

		// Check if we generated something we can take
		if (!generatedRandoms.empty()) {
			firstCandidate = generatedRandoms[0];

			generatedRandoms.erase(
				std::find(generatedRandoms.begin(), generatedRandoms.end(), firstCandidate)
			);
		}
		else {
			generated = distribution(gen);
			
			while (*boundIter < generated
				&& boundIter != boundsVector->end()) {
				firstCandidate++;
				boundIter++;
			}
			// firstCandidate--;
		}

		// Check if we generated something we can take
		if (!generatedRandoms.empty()) {
			for (int& vertex : generatedRandoms) {
				if (firstCandidate != vertex) {
					secondCandidate = vertex;
					generatedRandoms.erase(
						std::find(generatedRandoms.begin(), generatedRandoms.end(), secondCandidate)
					);
					break;
				}
			}
		}
		
		// If we did, save and end this iteration early
		if (secondCandidate != firstCandidate) 
			pairsVector.push_back(std::tuple<int, int>(firstCandidate, secondCandidate));
		
		// Otherwise, go standard
		else {
			while (true) {
				secondCandidate = -1;
				generated = distribution(gen);
				boundIter = boundsVector->begin();

				while (*boundIter < generated
					&& boundIter != boundsVector->end()) {
					secondCandidate++;
					boundIter++;
				}
				// secondCandidate--;

				// Found candidate, save and exit
				if (secondCandidate != firstCandidate) {
					pairsVector.push_back(std::tuple<int, int>(firstCandidate, secondCandidate));
					break;
				}
				// Otherwise save for later
				else generatedRandoms.push_back(secondCandidate);
			}
		}

		childrenGenerated++;
	}

	return pairsVector;
}

PathData Algorithms::generateChildOX(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	PathData generatedChild;
	
	std::tuple<int, int> t = generateRandomTwoPositions(0, matrix->size - 2);

	std::vector<short>::iterator pathIterF =
		firstParent->begin();
	std::vector<short>::iterator pathIterB =
		firstParent->begin();

	std::advance(pathIterF, std::get<0>(t));
	std::advance(pathIterB, std::get<1>(t));

	// Copy fragment
	std::vector<short> segmentFromParent(std::get<1>(t) - std::get<0>(t) + 1);
	segmentFromParent.reserve(firstParent->size());
	std::copy(pathIterF, pathIterB + 1, segmentFromParent.begin());

	// From last place till end
	pathIterB = secondParent->begin();
	std::advance(pathIterB, std::get<1>(t) + 1);

	int buffer = 0;
	while (pathIterB != secondParent->end()) {
		if (std::find(segmentFromParent.begin(), segmentFromParent.end(), *pathIterB)
			== segmentFromParent.end()) {
			segmentFromParent.push_back(*pathIterB);
		}
		else buffer++;

		pathIterB++;
	}

	// From begin to first place copied
	pathIterB = secondParent->begin();
	std::advance(pathIterB, std::get<0>(t) + buffer);
	pathIterF = secondParent->begin();
	std::vector<short> childPath;
	childPath.reserve(firstParent->size());

	while (pathIterF != pathIterB) {
		if (std::find(segmentFromParent.begin(), segmentFromParent.end(), *pathIterF)
			== segmentFromParent.end()) {
			childPath.push_back(*pathIterF);
		}
		else pathIterB++;

		pathIterF++;
	}

	// Connect two vectors
	childPath.insert(childPath.end(), segmentFromParent.begin(), segmentFromParent.end());
	
	generatedChild.pathOrder = std::move(childPath);

	return generatedChild;
}

PathData Algorithms::generateChildEAX(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	PathData generatedChild;

	if (*firstParent == *secondParent) {
		repeat++;
		std::uniform_real_distribution<> distribution(0.0, 1.0);
		// Mutate
		if (distribution(gen) <= mutationConstant * (repeat + 1) / total) {
			std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
			generatedChild = getNewOrder(firstParent, std::get<0>(iT),
				std::get<1>(iT), &(matrix->mat));
		}
		else generatedChild.pathOrder = std::vector<short>(firstParent->begin(), firstParent->end());
		
		return generatedChild;
	}

	int matrixSize = matrix->size;
	generatedChild.pathOrder.reserve(matrixSize - 1);
	// Edge tables
	EdgeTable edgeTable(matrixSize - 1);

	// Fill edge table
	updateTable(matrix, &edgeTable, firstParent, secondParent);

	EdgeTable referenceCopy = edgeTable;

	// Vertices left
	std::vector<short> verticesLeft(matrixSize - 1);
	std::iota(verticesLeft.begin(), verticesLeft.end(), 1);
	std::shuffle(verticesLeft.begin(), verticesLeft.end(), gen);

	// Current vertex
	std::uniform_int_distribution<> distribution(0, 1);
	short current = 0;
	if (distribution(gen))
		current = (*firstParent)[0];
	else current = (*secondParent)[0];

	while (!verticesLeft.empty()) {
		// Add vertex
		generatedChild.pathOrder.push_back(current);

		// Delete current from verticesLeft
		auto iterLeft = std::find(verticesLeft.begin(), verticesLeft.end(), current);

		if (iterLeft != verticesLeft.end())
			verticesLeft.erase(iterLeft);

		// Find occurences
		std::vector<short> toDelete(referenceCopy.singleEdge[current - 1].begin(),
			referenceCopy.singleEdge[current - 1].end());
		
		toDelete.resize(referenceCopy.singleEdge[current - 1].size() + 
			referenceCopy.doubleEdge[current - 1].size());

		std::copy_backward(referenceCopy.doubleEdge[current - 1].begin(),
			referenceCopy.doubleEdge[current - 1].end(), toDelete.end());

		// Delete occurences
		deleteOccurences(&edgeTable, current, &toDelete);
		
		// If not found next == -1
		short next = getNext(&edgeTable, current, generatedChild.pathOrder[0]);

		if (next == -1) {
			if (!verticesLeft.empty()) 
				current = verticesLeft[0];
			else break;
		}
		else current = next;
	}
	// std::cout << "\n";
	return generatedChild;
}

void Algorithms::updateTable(Matrix* matrix, EdgeTable* edgeTable, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	int matrixSize = matrix->size;
	if (matrixSize < 250) {
		for (int i = 1; i <= matrixSize - 1; i++) {
			updateTableFind(edgeTable, firstParent, i);
			updateTableFind(edgeTable, secondParent, i);
		}
	}
	else updateTableBig(edgeTable, firstParent, secondParent);
}

void Algorithms::updateTableFind(EdgeTable* edgeTable, std::vector<short>* parent, int vertexToFind) {
	std::vector<short>::iterator iter =
		std::find(parent->begin(), parent->end(), vertexToFind);

	// Check if double connection
		// Needs to be changed to checking singleEdge vectors, since "same connection" can be in different places of vectors
		// Maybe only look forward and for end() look at front()

	short firstValue = vertexToFind - 1, // -1 to correct for vector offset
	      secondValue = (*iter == parent->back() ? parent->front() - 1 : (*(iter + 1) - 1));

	bool first = edgeTable->isInSingle(firstValue, secondValue + 1);
	bool second = edgeTable->isInSingle(secondValue, firstValue + 1);

	// Exists in single
	if (first || second) {
		// Delete from single...
		edgeTable->singleEdge[firstValue].erase(
			std::find(edgeTable->singleEdge[firstValue].begin(),
				edgeTable->singleEdge[firstValue].end(),
				secondValue + 1)
		);
		edgeTable->singleEdge[secondValue].erase(
			std::find(edgeTable->singleEdge[secondValue].begin(),
				edgeTable->singleEdge[secondValue].end(),
				firstValue + 1)
		);
		// ... move to double
		edgeTable->doubleEdge[firstValue].push_back(secondValue + 1);
		edgeTable->doubleEdge[secondValue].push_back(firstValue + 1);
	}
	// Doesn't exist in single, add to single
	else {
		edgeTable->singleEdge[firstValue].push_back(secondValue + 1);
		edgeTable->singleEdge[secondValue].push_back(firstValue + 1);
	}
}

// Better for bigger matrices
void Algorithms::updateTableBig(EdgeTable* edgeTable, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	// Iterate through firstParent vertices
	for (auto iter = firstParent->begin(); iter != firstParent->end(); iter++) {
		short secondValue = (*iter == firstParent->back() ? firstParent->front() : *(iter + 1));

		edgeTable->singleEdge[*iter - 1].push_back(secondValue);
		edgeTable->singleEdge[secondValue - 1].push_back(*iter);
	}

	// Iterate through secondParent vertices
	for (auto iter = secondParent->begin(); iter != secondParent->end(); iter++) {
		short secondValue = (*iter == secondParent->back() ? secondParent->front() : *(iter + 1));

		bool first = edgeTable->isInSingle(*iter - 1, secondValue);
		bool second = edgeTable->isInSingle(secondValue - 1, *iter);

		// Exists in single
		if (first || second) {
			// Delete from single...
			edgeTable->singleEdge[*iter - 1].erase(
				std::find(edgeTable->singleEdge[*iter - 1].begin(),
					edgeTable->singleEdge[*iter - 1].end(),
					secondValue)
			);
			edgeTable->singleEdge[secondValue - 1].erase(
				std::find(edgeTable->singleEdge[secondValue - 1].begin(),
					edgeTable->singleEdge[secondValue - 1].end(),
					*iter)
			);
			// ... move to double
			edgeTable->doubleEdge[*iter - 1].push_back(secondValue);
			edgeTable->doubleEdge[secondValue - 1].push_back(*iter);
		}
		// Doesn't exist in single, add to single
		else {
			edgeTable->singleEdge[*iter - 1].push_back(secondValue);
			edgeTable->singleEdge[secondValue - 1].push_back(*iter);
		}
	}
}

void Algorithms::deleteOccurences(EdgeTable* edgeTable, short currentVertex, std::vector<short>* vertexToDelete) {
	std::vector<short>::iterator iterOccur;

	for (short& s : *vertexToDelete) {
		iterOccur = std::find(edgeTable->singleEdge[s - 1].begin(),
			edgeTable->singleEdge[s - 1].end(), currentVertex);

		// If found
		if (iterOccur != edgeTable->singleEdge[s - 1].end())
			edgeTable->singleEdge[s - 1].erase(iterOccur);

		iterOccur = std::find(edgeTable->doubleEdge[s - 1].begin(),
			edgeTable->doubleEdge[s - 1].end(), currentVertex);

		// If found
		if (iterOccur != edgeTable->doubleEdge[s - 1].end()) 
			edgeTable->doubleEdge[s - 1].erase(iterOccur);
	}
}

short Algorithms::getNext(EdgeTable* edgeTable, short current, short currentFallback) {
	short next = 0;
	int tableSize = INT_MAX;
	std::vector<short>::iterator temp = edgeTable->doubleEdge[current - 1].begin();

	// Double edges exist
	if (!edgeTable->doubleEdge[current - 1].empty()) {
		while (temp != edgeTable->doubleEdge[current - 1].end()) {
			// If table size is lower than previous, set as next
			if (edgeTable->singleEdge[*temp - 1].size() +
				edgeTable->doubleEdge[*temp - 1].size() < INT_MAX) {
				tableSize = edgeTable->singleEdge[*temp - 1].size() +
					edgeTable->doubleEdge[*temp - 1].size();

				next = *temp;
			}
			temp++;
		}

		// Delete occurences from local
		edgeTable->doubleEdge[current - 1].erase(
			std::find(edgeTable->doubleEdge[current - 1].begin(),
			edgeTable->doubleEdge[current - 1].end(), next)
		);

		return next;
	}

	// Single edges exist
	if (!edgeTable->singleEdge[current - 1].empty()) {
		temp = edgeTable->singleEdge[current - 1].begin();

		while (temp != edgeTable->singleEdge[current - 1].end()) {
			// If table size is lower than previous, set as next
			if (edgeTable->singleEdge[*temp - 1].size() +
				edgeTable->doubleEdge[*temp - 1].size() < INT_MAX) {
				tableSize = edgeTable->singleEdge[*temp - 1].size() +
					edgeTable->doubleEdge[*temp - 1].size();

				next = *temp;
			}
			temp++;
		}

		// Delete occurences from local
		edgeTable->singleEdge[current - 1].erase(
			std::find(edgeTable->singleEdge[current - 1].begin(),
				edgeTable->singleEdge[current - 1].end(), next)
		);

		return next;
	}

	// Check other side
	if (currentFallback) 
		return getNext(edgeTable, currentFallback, 0);

	return -1;
}

std::vector<PathData> Algorithms::getNewRandomGeneration(Matrix* matrix, std::vector<PathData>* parents, std::vector<PathData>* children) {
	if (children->empty())
		return *parents;
	
	int pSize = parents->size();
	std::vector<PathData> newGeneration;//(parents->size() / 2);
	newGeneration.reserve(pSize);

	auto compare = [](PathData x, PathData y)
		{ return x.pathLength < y.pathLength; };
	// Only use when you wanna be elite about children
	std::partial_sort(children->begin(), children->begin() + 1, children->end(), compare);

	std::vector<int> possibleIndexCombined(pSize + children->size() - 2);

	// Doesn't include position of elite child
	std::iota(possibleIndexCombined.begin(), possibleIndexCombined.begin() + pSize - 1, 1);
	std::iota(possibleIndexCombined.begin() + pSize - 1, possibleIndexCombined.end(), pSize + 1);
		
	std::shuffle(possibleIndexCombined.begin(), possibleIndexCombined.end(), gen);
	
	newGeneration.emplace_back(std::move((*parents)[0]));
	newGeneration.emplace_back(std::move((*children)[0]));

	for (int i = 0; i < startingPopulationSize - 2; i++) {
		if (possibleIndexCombined[i] < pSize)
			newGeneration.emplace_back(std::move((*parents)[possibleIndexCombined[i]]));
		else
			newGeneration.emplace_back(std::move((*children)[possibleIndexCombined[i] - pSize]));
	}

	return newGeneration;
}