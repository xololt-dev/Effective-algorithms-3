#include "util.hpp"

#include <numeric>
#include <thread>

// https://dl.acm.org/doi/pdf/10.1145/3071178.3071305
// populacja - ocena - selekcja - pula rodzicielska - operacje genetyczne - subpopulacja - ocena - sukcesja - ...

void Algorithms::geneticAlgorithm(Matrix* matrix) {
	if (!matrix->size) {
		std::cout << "Nie wczytano matrycy!\n";
		return;
	}
	
	if (currentCrossoverType)
		geneticEAX(matrix);
	else
		geneticOX(matrix);
}

// TODO: check if swapping to vector also helps here
void Algorithms::geneticOX(Matrix* matrix) {
	// Priority queue func
	struct {
		bool operator()(const QueueData x, const QueueData y) const
			{ return x.pathLength > y.pathLength; };
	} compareS;

	pathLength = INT_MAX;

	// Starting population generation
	std::vector<QueueData> qD = generateStartingPopulation(matrix);
	std::priority_queue<QueueData, std::vector<QueueData>, decltype(compareS)>
		queue(qD.begin(), qD.end());

	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

	while (std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - start) < maxExecutionTime) {

		// Update best solution if applicable
		if (queue.top().pathLength < pathLength) {
			pathLength = queue.top().pathLength;
			vertexOrder = queue.top().pathOrder;
			runningTime = std::chrono::duration_cast<std::chrono::microseconds>
				(std::chrono::steady_clock::now() - start);

			for (short s : queue.top().pathOrder) {
				std::cout << s << " ";
			}
			std::cout << "\n";
		}

		// Get lower bounds
		std::vector<double> lowerBounds = getVertexLowerBounds(queue.size());

		// Generate pairs
		std::vector<std::tuple<short, short>> parents = 
			generateParents(&lowerBounds, queue.size());

		// Dump gen. to vector
		std::vector<QueueData> parentsData;
		std::vector<QueueData> childrenData;
		while (!queue.empty()) {
			parentsData.push_back(queue.top());
			queue.pop();
		}

		// Genetic operations
		std::uniform_real_distribution<> distribution(0.0, 1.0);
		for (std::tuple<short, short> t : parents) {
			if (distribution(gen) > crossoverConstant)
				continue;

			// Crossover
			QueueData child = generateChildOX(matrix, &parentsData[std::get<0>(t)].pathOrder,
				&parentsData[std::get<1>(t)].pathOrder);
			
			// Mutate
			if (distribution(gen) <= mutationConstant) {
				std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
				child = getNewOrder(&child.pathOrder, std::get<0>(iT),
					std::get<1>(iT), &(matrix->mat));
			}
			else child.pathLength = calculateCandidate(&child.pathOrder, matrix);

			childrenData.push_back(child);
		}

		// Combine to get new generation
		std::vector<QueueData> newGen = getNewRandomGeneration(matrix,
			&parentsData, &childrenData);

		queue = std::priority_queue<QueueData, std::vector<QueueData>, 
			decltype(compareS)> (newGen.begin(), newGen.end());
	}
}

void Algorithms::geneticEAX(Matrix* matrix) {
	int iterations = 0;
	// Priority queue func
	struct {
		bool operator()(const QueueData x, const QueueData y) const
		{
			return x.pathLength < y.pathLength;
		};
	} compareS;

	pathLength = INT_MAX;

	// Starting population generation
	std::vector<QueueData> population = generateStartingPopulation(matrix);

	std::sort(population.begin(), population.end(), compareS);

	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

	while (std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - start) < maxExecutionTime) {

		// Update best solution if applicable
		if (population[0].pathLength < pathLength) {
			pathLength = population[0].pathLength;
			vertexOrder = population[0].pathOrder;
			runningTime = std::chrono::duration_cast<std::chrono::microseconds>
				(std::chrono::steady_clock::now() - start);

			for (short s : population[0].pathOrder) {
				std::cout << s << " ";
			}
			std::cout << "\n";
		}

		// Get lower bounds
		std::vector<double> lowerBounds = getVertexLowerBounds(population.size());

		// Generate pairs
		std::vector<std::tuple<short, short>> parents =
			generateParents(&lowerBounds, population.size());

		// Dump gen. to vector
		// std::vector<QueueData> parentsData(qD.begin(), qD.end());
		//parents.reserve(queue.size());
		std::vector<QueueData> childrenData;

		// Genetic operations
		std::uniform_real_distribution<> distribution(0.0, 1.0);
		for (std::tuple<short, short> t : parents) {
			if (distribution(gen) > crossoverConstant)
				continue;

			// Crossover
			// Can generate too big children
			QueueData child = generateChildEAX(matrix, &population[std::get<0>(t)].pathOrder,
				&population[std::get<1>(t)].pathOrder);

			// Mutate
			if (distribution(gen) <= mutationConstant) {
				std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
				child = getNewOrder(&child.pathOrder, std::get<0>(iT),
					std::get<1>(iT), &(matrix->mat));
			}
			else child.pathLength = calculateCandidate(&child.pathOrder, matrix);

			childrenData.push_back(child);
		}

		// Combine to get new generation
		population = getNewRandomGeneration(matrix, &population, &childrenData);

		std::sort(population.begin(), population.end(), compareS);
		iterations++;
	}
	std::cout << iterations << "\n";
}

std::vector<QueueData> Algorithms::generateStartingPopulation(Matrix* matrix) {
	std::vector<QueueData> returnVec;
	std::vector<short> leftVertices(matrix->size - 1);
	std::iota(leftVertices.begin(), leftVertices.end(), 1);
	QueueData tempData;

	if (startingPopulationRandom == true) {
		for (int i = 0; i < startingPopulationSize; i++) {
			std::shuffle(leftVertices.begin() + 1, leftVertices.end(), gen);

			short previousVertex = 0;
			tempData.pathLength = 0;
			tempData.pathOrder = leftVertices;

			for (int k = 0; k < tempData.pathOrder.size(); k++) {
				tempData.pathLength += matrix->mat[tempData.pathOrder[k]][previousVertex];
				previousVertex = tempData.pathOrder[k];
			}

			tempData.pathLength += matrix->mat[0][previousVertex];

			returnVec.push_back(tempData);
		}
	}
	else {
		std::tuple<std::vector<short>, int> t = generateInitialSolution(matrix);
		short firstOne = std::get<0>(t)[0];
		leftVertices.erase(
			std::find(leftVertices.begin(), leftVertices.end(), firstOne)
		);

		tempData.pathOrder = std::get<0>(t);
		tempData.pathLength = std::get<1>(t);
		tempData.anchorOne = 0;
		tempData.anchorTwo = 0;
		returnVec.push_back(tempData);

		if (startingPopulationSize < matrix->size) {
			for (int i = 1; i < startingPopulationSize; i++) {
				std::shuffle(leftVertices.begin() + 1, leftVertices.end(), gen);

				t = generateNewSolutionV(matrix, leftVertices[0]);

				firstOne = std::get<0>(t)[0];
				leftVertices.erase(
					std::find(leftVertices.begin(), leftVertices.end(), firstOne)
				);

				tempData.pathOrder = std::get<0>(t);
				tempData.pathLength = std::get<1>(t);
				tempData.anchorOne = 0;
				tempData.anchorTwo = 0;
				returnVec.push_back(tempData);
			}
		}
		else {
			for (int i = 2; i < matrix->size; i++) {
				std::shuffle(leftVertices.begin() + 1, leftVertices.end(), gen);

				t = generateNewSolutionV(matrix, leftVertices[0]);

				firstOne = std::get<0>(t)[0];
				leftVertices.erase(
					std::find(leftVertices.begin(), leftVertices.end(), firstOne)
				);

				tempData.pathOrder = std::get<0>(t);
				tempData.pathLength = std::get<1>(t);
				tempData.anchorOne = 0;
				tempData.anchorTwo = 0;
				returnVec.push_back(tempData);
			}

			std::vector<short> leftVertices(matrix->size - 1);
			std::iota(leftVertices.begin(), leftVertices.end(), 1);

			while (returnVec.size() < startingPopulationSize) {
				std::shuffle(leftVertices.begin() + 1, leftVertices.end(), gen);

				short previousVertex = 0;
				tempData.pathLength = 0;
				tempData.pathOrder = leftVertices;

				for (int k = 0; k < tempData.pathOrder.size(); k++) {
					tempData.pathLength += matrix->mat[tempData.pathOrder[k]][previousVertex];
					previousVertex = tempData.pathOrder[k];
				}

				tempData.pathLength += matrix->mat[0][previousVertex];

				returnVec.push_back(tempData);
			}
		}		
	}
	
	return returnVec;
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

std::vector<std::tuple<short, short>> Algorithms::generateParents(std::vector<double>* boundsVector, int vectorSize) {
	std::vector<std::tuple<short, short>> pairsVector;
	// Maybe saved "wrong" generations to a vector for later usage?
	std::vector<short> generatedRandoms;

	// Generate number(s) [0, 1]
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	int childrenGenerated = 0;

	while (childrenGenerated < vectorSize) {
		double generated;
		int firstCandidate = 0;
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
			firstCandidate--;
		}

		// Check if we generated something we can take
		if (!generatedRandoms.empty()) {
			for (short vertex : generatedRandoms) {
				if (firstCandidate != vertex) {
					secondCandidate = vertex;
					generatedRandoms.erase(
						std::find(generatedRandoms.begin(), generatedRandoms.end(), secondCandidate)
					);
					break;
				}
			}
			// If we did, save and end this iteration early
			if (secondCandidate != firstCandidate) {
				pairsVector.push_back(std::tuple<short, short>(firstCandidate, secondCandidate));
				childrenGenerated++;
				continue;
			}
		}
		
		/*
		do {
			secondCandidate = 0;
			generatedSecond = distribution(gen);
			boundIter = boundsVector->begin();

			while (*boundIter < generatedSecond
				&& boundIter != boundsVector->end()) {
				secondCandidate++;
				boundIter++;
			}
			secondCandidate--;
		} while (secondCandidate == firstCandidate);
		*/
		// Otherwise, go standard
		while (true) {
			secondCandidate = 0;
			generated = distribution(gen);
			boundIter = boundsVector->begin();

			while (*boundIter < generated
				&& boundIter != boundsVector->end()) {
				secondCandidate++;
				boundIter++;
			}
			secondCandidate--;
			// Found candidate, save and exit
			if (secondCandidate != firstCandidate) {
				pairsVector.push_back(std::tuple<short, short>(firstCandidate, secondCandidate));
				break;
			}
			// Otherwise save for later
			else generatedRandoms.push_back(secondCandidate);
		}

		childrenGenerated++;
	}

	return pairsVector;
}

QueueData Algorithms::generateChildOX(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	QueueData generatedChild;
	
	std::tuple<int, int> t = generateRandomTwoPositions(0, matrix->size - 2);

	std::vector<short>::iterator pathIterF =
		firstParent->begin();
	std::vector<short>::iterator pathIterB =
		firstParent->begin();

	std::advance(pathIterF, std::get<0>(t));
	std::advance(pathIterB, std::get<1>(t));

	// Copy fragment
	std::vector<short> segmentFromParent(std::get<1>(t) - std::get<0>(t) + 1);
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
	
	generatedChild.pathOrder = childPath;

	return generatedChild;
}

QueueData Algorithms::generateChildEAX(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	QueueData generatedChild;

	if (*firstParent == *secondParent) {
		generatedChild.pathOrder = *firstParent;
		return generatedChild;
	}
		
	int matrixSize = matrix->size;
	// Edge tables
	EdgeTable edgeTable(matrixSize - 1);

	// Fill edge table
	for (int i = 1; i <= firstParent->size(); i++) {
		updateTable(&edgeTable, firstParent, i);
		updateTable(&edgeTable, secondParent, i);
	}

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

		short next = -1;
		// Find occurences
		// values with offset -1
		std::vector<short> occurences = findOccurences(&edgeTable, &next, current);

		// Delete occurences
		// If not found next == -1
		next = getNext(&edgeTable, current, generatedChild.pathOrder[0]);
		if (next == -1 && !verticesLeft.empty())
			current = verticesLeft[0];
		else current = next;
	}

	return generatedChild;
}

void Algorithms::updateTable(EdgeTable* edgeTable, std::vector<short>* parent, int vertexToFind) {
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

std::vector<short> Algorithms::findOccurences(EdgeTable* edgeTable, short* next, short currentVertex) {
	std::vector<short> occurences;
	std::vector<short>::iterator iterOccur;

	for (int i = 0; i < edgeTable->singleEdge.size(); i++) {
		// Skip self check
		if (i + 1 == currentVertex)
			continue;

		// Max amount of occurences is 4, so discontinue search when max size reached
		if (occurences.size() == 4)
			break;

		iterOccur = std::find(edgeTable->singleEdge[i].begin(), edgeTable->singleEdge[i].end(), currentVertex);
		// If found
		if (iterOccur != edgeTable->singleEdge[i].end()) {
			occurences.push_back(i);
			edgeTable->singleEdge[i].erase(iterOccur);
			continue;
		}

		iterOccur = std::find(edgeTable->doubleEdge[i].begin(), edgeTable->doubleEdge[i].end(), currentVertex);
		// If found
		if (iterOccur != edgeTable->doubleEdge[i].end()) {
			occurences.push_back(i);
			*next = i + 1;
			edgeTable->doubleEdge[i].erase(iterOccur);
		}
		// All possibilities of move in occurences
	}

	return occurences;
}

short Algorithms::getNext(EdgeTable* edgeTable, short current, short currentFallback) {
	short next = 0,
		  tableSize = SHRT_MAX;
	std::vector<short>::iterator temp = edgeTable->doubleEdge[current - 1].begin();

	// Double edges exist
	if (!edgeTable->doubleEdge[current - 1].empty()) {
		while (temp != edgeTable->doubleEdge[current - 1].end()) {
			// If table size is lower than previous, set as next
			if (edgeTable->singleEdge[*temp - 1].size() +
				edgeTable->doubleEdge[*temp - 1].size() < SHRT_MAX) {
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
				edgeTable->doubleEdge[*temp - 1].size() < SHRT_MAX) {
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
	if (currentFallback) {
		return getNext(edgeTable, currentFallback, 0);
	}

	return -1;
}

std::vector<QueueData> Algorithms::getNewRandomGeneration(Matrix* matrix, std::vector<QueueData>* parents, std::vector<QueueData>* children/*, void* structure*/) {
	std::vector<QueueData> newGeneration;

	if (children->empty()) {
		for (int i = 0; i < parents->size(); i++)
			newGeneration.push_back((*parents)[i]);
		return newGeneration;
	}

	std::vector<short> possibleIndexCombined(parents->size() + children->size() - 2);

	std::iota(possibleIndexCombined.begin(), possibleIndexCombined.end(), 1);
	
	if (currentCrossoverType) {
		auto compare = [](QueueData x, QueueData y)
		{ return x.pathLength < y.pathLength; };
		// Only use when you wanna be elite about children
		std::sort(children->begin(), children->end(), compare);
	}
	else {
		auto compare = [](QueueData x, QueueData y)
		{ return x.pathLength > y.pathLength; };
		// Only use when you wanna be elite about children
		std::sort(children->begin(), children->end(), compare);
	}
	
	std::vector<short>::iterator iter = 
		std::find(possibleIndexCombined.begin(), possibleIndexCombined.end(), parents->size());
	
	if (iter != possibleIndexCombined.end())
		possibleIndexCombined.erase(iter);

	std::shuffle(possibleIndexCombined.begin(), possibleIndexCombined.end(), gen);

	newGeneration.push_back((*parents)[0]);
	newGeneration.push_back((*children)[0]);

	for (int i = 0; i < startingPopulationSize - 2; i++) {
		if (possibleIndexCombined[i] < parents->size())
			newGeneration.push_back((*parents)[possibleIndexCombined[i]]);
		else
			newGeneration.push_back((*children)[possibleIndexCombined[i] - parents->size()]);
	}

	return newGeneration;
}