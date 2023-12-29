#include "util.hpp"

#include <numeric>

// https://dl.acm.org/doi/pdf/10.1145/3071178.3071305
// populacja - ocena - selekcja - pula rodzicielska - operacje genetyczne - subpopulacja - ocena - sukcesja - ...

void Algorithms::geneticAlgorithm(Matrix* matrix) {
	geneticOX(matrix);
}

void Algorithms::geneticOX(Matrix* matrix) {
	// Priority queue
	struct {
		bool operator()(const QueueData x, const QueueData y) const 
		{ return x.pathLength > y.pathLength; };
	} compareS;

	auto compare = [](QueueData x, QueueData y)
	{ return x.pathLength > y.pathLength; };

	// Starting population generation
	std::vector<QueueData> qD = generateStartingPopulation(matrix);
	std::priority_queue<QueueData, std::vector<QueueData>, decltype(compareS)> queue(qD.begin(), qD.end());
	// std::priority_queue<QueueData, std::vector<QueueData>,
		//decltype(compare)> queue(compare);
	
	std::vector<short> currentSolution = queue.top().pathOrder;
	int currentLength = queue.top().pathLength;
	int currentPopulationSize = queue.size();
	int iterationsDone = 1;
	std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

	while (std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - start) < maxExecutionTime) {
	
		// Grading
		if (queue.top().pathLength < currentLength) {
			currentLength = queue.top().pathLength;
			currentSolution = queue.top().pathOrder;
		}

			// Certain bounds

		// Selection
			// Generate possibilities
		std::vector<double> probabilites(queue.size());
		std::vector<double> vertexLowerBound(queue.size());
		double sum = 0.0;
		for (int i = 0; i < probabilites.size(); i++) {
			vertexLowerBound.push_back(sum);
			double current = (1.0 - (double)(i + 1) / (double)currentPopulationSize);
			probabilites[i] = current;
			sum += current;
		}

		for (int i = 0; i < probabilites.size(); i++) {
			probabilites[i] /= sum;
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
			int firstCandidate = -1;
			int secondCandidate = -1;
			std::vector<double>::iterator boundIter = vertexLowerBound.begin();
			while (*boundIter < generated 
				&& boundIter != vertexLowerBound.end()) {
				firstCandidate++;
				boundIter++;
			}
			
			double generatedSecond;
			do {
				secondCandidate = -1;
				generatedSecond = distribution(gen);

				while (*boundIter < generated 
					&& boundIter != vertexLowerBound.end()) {
					firstCandidate++;
					boundIter++;
				}
			} while (secondCandidate == firstCandidate);

			// two diff parents
		// Genetic operations
			// OX
			std::tuple<int, int> t = generateRandomTwoPositions(0, queue.size());

			std::vector<short>::iterator pathIterF = 
				queueCandidates[firstCandidate].pathOrder.begin();
			std::vector<short>::iterator pathIterB =
				queueCandidates[firstCandidate].pathOrder.begin();
			
			std::advance(pathIterF, std::get<0>(t));
			std::advance(pathIterB, std::get<1>(t));

			// Copy fragment
			std::vector<short> segmentFromParent(pathIterF, pathIterB);

			// From last place till end
			pathIterB = queueCandidates[secondCandidate].pathOrder.begin();
			std::advance(pathIterB, std::get<1>(t) + 1);

			while (pathIterB != queueCandidates[secondCandidate].pathOrder.end()) {
				if (std::find(segmentFromParent.begin(), segmentFromParent.end(), *pathIterB) 
					== segmentFromParent.end()) {
					segmentFromParent.push_back(*pathIterB);
				}

				pathIterB++;
			}

			// From begin to first place copied
			pathIterB = queueCandidates[secondCandidate].pathOrder.begin();
			std::advance(pathIterB, std::get<0>(t));
			pathIterF = queueCandidates[secondCandidate].pathOrder.begin();
			std::vector<short> childPath;

			while (pathIterF != pathIterB) {
				if (std::find(segmentFromParent.begin(), segmentFromParent.end(), *pathIterF)
					== segmentFromParent.end()) {
					childPath.push_back(*pathIterF);
				}

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
		// generatedChildren

		// Grading 2

			// Update population size (aka do we need to generate more)
		currentPopulationSize = startingPopulationSize + log2(iterationsDone);

		iterationsDone++;
	}
}

std::vector<QueueData> Algorithms::generateStartingPopulation(Matrix* matrix) {
	std::vector<QueueData> returnVec;
	std::vector<short> leftVertices;
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
		returnVec.push_back(tempData);

		for (int i = 1; i < startingPopulationSize; i++) {
			std::shuffle(leftVertices.begin() + 1, leftVertices.end(), gen);

			t = generateNewSolutionV(matrix, leftVertices[0]);

			firstOne = std::get<0>(t)[0];
			leftVertices.erase(
				std::find(leftVertices.begin(), leftVertices.end(), firstOne)
			);

			tempData.pathOrder = std::get<0>(t);
			tempData.pathLength = std::get<1>(t);
			returnVec.push_back(tempData);
		}
	}
	
	return returnVec;
}
