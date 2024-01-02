#include "util.hpp"

#include <numeric>

// https://dl.acm.org/doi/pdf/10.1145/3071178.3071305
// populacja - ocena - selekcja - pula rodzicielska - operacje genetyczne - subpopulacja - ocena - sukcesja - ...

void Algorithms::geneticAlgorithm(Matrix* matrix) {
	if (!matrix->size) {
		std::cout << "Nie wczytano matrycy!\n";
		return;
	}

	// geneticOX(matrix);
	newGeneticOX(matrix);
}

void Algorithms::geneticOX(Matrix* matrix) {
	// Priority queue
	struct {
		bool operator()(const QueueData x, const QueueData y) const 
		{ return x.pathLength > y.pathLength; };
	} compareS;

	// Starting population generation
	std::vector<QueueData> qD = generateStartingPopulation(matrix);
	std::priority_queue<QueueData, std::vector<QueueData>, decltype(compareS)> 
		queue(qD.begin(), qD.end());
	// std::priority_queue<QueueData, std::vector<QueueData>,
		//decltype(compare)> queue(compare);
	
	std::vector<short> currentSolution = queue.top().pathOrder;
	int currentLength = queue.top().pathLength;
	int currentPopulationSize = queue.size();
	int iterationsDone = 1;
	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

	while (std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - start) < maxExecutionTime) {
		
		currentPopulationSize = queue.size();

		// Grading
		if (queue.top().pathLength < currentLength) {
			currentLength = queue.top().pathLength;
			currentSolution = queue.top().pathOrder;
			runningTime = std::chrono::duration_cast<std::chrono::microseconds>
				(std::chrono::steady_clock::now() - start);
		}

			// Certain bounds

		// Selection
			// Generate possibilities
		std::vector<double> probabilites(queue.size());
		std::vector<double> vertexLowerBound(queue.size(), 0.0);
		double sum = 0.0;
		for (int i = 0; i < probabilites.size(); i++) {
			double current = (1.0 - (double)(i + 1) / (double)currentPopulationSize);
			probabilites[i] = current;
			sum += current;
		}

		for (int i = 0; i < probabilites.size(); i++) {
			probabilites[i] /= sum;
		}

		for (int i = 1; i < vertexLowerBound.size(); i++) {
			vertexLowerBound[i] = vertexLowerBound[i - 1] + probabilites[i - 1];
		}

		std::vector<QueueData> queueCandidates;
		std::vector<QueueData> generatedChildren;
		while (!queue.empty()) {
			queueCandidates.push_back(queue.top());
			queue.pop();
		}
		// Generate number(s) [0, 1]
		std::uniform_real_distribution<> distribution(0.0, 1.0);

		int childrenGenerated = 0;
		while (childrenGenerated < currentPopulationSize) {
			double generated = distribution(gen);
			int firstCandidate = 0;
			int secondCandidate = 0;
			std::vector<double>::iterator boundIter = vertexLowerBound.begin();
			while (*boundIter < generated 
				&& boundIter != vertexLowerBound.end()) {
				firstCandidate++;
				boundIter++;
			}
			firstCandidate--;
			
			double generatedSecond;
			do {
				secondCandidate = 0;
				generatedSecond = distribution(gen);
				boundIter = vertexLowerBound.begin();

				while (*boundIter < generatedSecond
					&& boundIter != vertexLowerBound.end()) {
					secondCandidate++;
					boundIter++;
				}
				secondCandidate--;
			} while (secondCandidate == firstCandidate);

			// two diff parents
		// Genetic operations
			// OX
			std::tuple<int, int> t = generateRandomTwoPositions(0, matrix->size - 2);//currentPopulationSize);

			std::vector<short>::iterator pathIterF = 
				queueCandidates[firstCandidate].pathOrder.begin();
			std::vector<short>::iterator pathIterB =
				queueCandidates[firstCandidate].pathOrder.begin();
			
			std::advance(pathIterF, std::get<0>(t));
			std::advance(pathIterB, std::get<1>(t));

			// Copy fragment
			std::vector<short> segmentFromParent(std::get<1>(t) - std::get<0>(t) + 1);
			std::copy(pathIterF, pathIterB + 1, segmentFromParent.begin());

			// From last place till end
			pathIterB = queueCandidates[secondCandidate].pathOrder.begin();
			std::advance(pathIterB, std::get<1>(t) + 1);
			
			int buffer = 0;
			while (pathIterB != queueCandidates[secondCandidate].pathOrder.end()) {
				if (std::find(segmentFromParent.begin(), segmentFromParent.end(), *pathIterB)
					== segmentFromParent.end()) {
					segmentFromParent.push_back(*pathIterB);
				}
				else buffer++;

				pathIterB++;
			}

			// From begin to first place copied
			pathIterB = queueCandidates[secondCandidate].pathOrder.begin();
			std::advance(pathIterB, std::get<0>(t) + buffer);
			pathIterF = queueCandidates[secondCandidate].pathOrder.begin();
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

			QueueData d;
			d.pathOrder = childPath;

			generatedChildren.push_back(d);

			childrenGenerated++;
		}
			
		// Subpopulation
			// Update population size (aka do we need to generate more)
		//if (distribution(gen) > crossoverConstant)
			//continue;
		// currentPopulationSize = startingPopulationSize + log2(iterationsDone);
		int previousGenNumber = (currentPopulationSize + 2 - 1) / 2; // ceil int

		std::vector<QueueData>::iterator queueDataIter = queueCandidates.begin();
		std::advance(queueDataIter, previousGenNumber);
		queueCandidates.erase(queueDataIter, queueCandidates.end());

		queue = std::priority_queue<QueueData, std::vector<QueueData>, decltype(compareS)>
			(queueCandidates.begin(), queueCandidates.end());

		std::priority_queue<QueueData, std::vector<QueueData>, decltype(compareS)>
			queueChildren;

		for (int i = 0; i < generatedChildren.size(); i++) {
			double randomNum = distribution(gen);
			if (randomNum > crossoverConstant)
				continue;

			QueueData d;

			if (randomNum <= mutationConstant) {
				std::tuple<int, int> iT = generateRandomTwoPositions(0, matrix->size - 2);
				d = getNewOrder(&(generatedChildren[i].pathOrder), std::get<0>(iT), std::get<1>(iT), &(matrix->mat));
			}
			else {
				d.pathLength = calculateCandidate(&(generatedChildren[i].pathOrder), matrix);
				d.pathOrder = generatedChildren[i].pathOrder;
			}

			queueChildren.push(d);
		}

		// Not enough children, fill with previous
		if (queueChildren.size() + previousGenNumber < currentPopulationSize) {
			queue = queueChildren;

			// Fill
			for (int i = 0; i < queueCandidates.size(); i++)
				queue.push(queueCandidates[i]);
		}
		// Enough children - standard distrib.
		else {
			std::vector<QueueData>::iterator queueDataIter = queueCandidates.begin();
			std::advance(queueDataIter, previousGenNumber);
			queueCandidates.erase(queueDataIter, queueCandidates.end());

			queue = std::priority_queue<QueueData, std::vector<QueueData>, decltype(compareS)>
				(queueCandidates.begin(), queueCandidates.end());

			while (!queueChildren.empty()) {
				queue.push(queueChildren.top());
				queueChildren.pop();
			}
			/*
			for (int i = 0; i < currentPopulationSize - previousGenNumber; i++) {
				queue.push(queueChildren.top());
				queueChildren.pop();
			}
			*/
		}

		if (queue.top().pathLength < currentLength) {
			currentLength = queue.top().pathLength;
			currentSolution = queue.top().pathOrder;
			runningTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
		}

		iterationsDone++;
	}

	pathLength = currentLength;
	vertexOrder = currentSolution;
}

void Algorithms::newGeneticOX(Matrix* matrix) {
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
				&parentsData[std::get<0>(t)].pathOrder);
			
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

// TODO EAX
void Algorithms::geneticEAX(Matrix* matrix) {
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
				&parentsData[std::get<0>(t)].pathOrder);

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
			decltype(compareS)>(newGen.begin(), newGen.end());
	}
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

std::vector<std::vector<short>> Algorithms::createESet(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent) {
	std::vector<std::vector<short>> eSet;

	// If parents ==, exit
	if (*firstParent == *secondParent)
		return;

	// Selected parent to "take edge from"
	int parent = 0;
	// Current vertex
	short current = 0;

	// Choose vertex randomly to start sub cycle process
	std::uniform_int_distribution<> distribution(1, firstParent->size());
	current = distribution(gen);

	std::vector<short> possibleVertices(firstParent->size());
	std::iota(possibleVertices.begin(), possibleVertices.end(), 1);
	// delete current from it

	// Take turns taking edges/vertices from parents
	do {
		if (!parent) {
			parent = 1;
		}
		else {
			parent = 0;
		}

	} while (true);
	
	// Check if E-Set has been achieved

	return eSet;
}

std::vector<QueueData> Algorithms::getNewRandomGeneration(Matrix* matrix, std::vector<QueueData>* parents, std::vector<QueueData>* children) {
	std::vector<QueueData> newGeneration;

	if (children->empty()) {
		for (int i = 0; i < parents->size(); i++)
			newGeneration.push_back((*parents)[i]);
		return newGeneration;
	}

	std::vector<short> possibleIndexCombined(parents->size() + children->size() - 2);

	std::iota(possibleIndexCombined.begin(), possibleIndexCombined.end(), 1);
	
	auto compare = [](QueueData x, QueueData y)
	{ return x.pathLength > y.pathLength; };

	// Only use when you wanna be elite about children
	std::sort(children->begin(), children->end(), compare);
	
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