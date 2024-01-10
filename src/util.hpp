#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <iostream>
#include <tuple>
#include <queue>

enum MutationType {
	INVERSE,
	SWAP,
	INSERT,
	INSERT_SUB
};

class Matrix {
public:
	int size = 0;
	// wypelniane wierszami
	std::vector<std::vector<int>> mat;

	bool loadFromFile(std::string fileName);
	void display();

private:
	void loadFromTXT(std::string fileName);
	void loadFromATSP(std::string fileName);
};

struct PathData {
	std::vector<short> pathOrder;
	int pathLength;

	PathData() {}

	~PathData() {}

	PathData(PathData& other) {		//copy constructor
		pathOrder = other.pathOrder;
		pathLength = other.pathLength;
	}
	
	PathData& operator=(PathData other) {
		pathOrder = other.pathOrder;
		pathLength = other.pathLength;
		/*
		std::swap(pathOrder, other.pathOrder);
		std::swap(pathLength, other.pathLength);
		std::swap(anchorOne, other.anchorOne);
		std::swap(anchorTwo, other.anchorTwo);
		*/
		return *this;
	}
	
	PathData(PathData&& qd) noexcept :
		pathOrder(std::move(qd.pathOrder)),       // explicit move of a member of class type
		pathLength(std::exchange(qd.pathLength, 0)) // explicit move of a member of non-class type
	{}
};

struct EdgeTable{
	std::vector<std::vector<short>> singleEdge,
									doubleEdge;

	bool isInSingle(int valueOne, int valueTwo) {
		if (singleEdge[valueOne].empty())
			return false;

		// Not found
		if (std::find(singleEdge[valueOne].begin(),
			singleEdge[valueOne].end(),
			valueTwo) == singleEdge[valueOne].end())
			return false;

		return true;
	};

	void display() {
		for (int i = 0; i < singleEdge.size(); i++) {
			if (!singleEdge[i].empty() || !doubleEdge[i].empty()) {
				std::cout << i + 1 << ": ";
				for (short s : singleEdge[i])
					std::cout << s << " ";
				for (short s : doubleEdge[i])
					std::cout << s << "+ ";
				std::cout << "\n";
			}
		}
		std::cout << "\n";
	}

	EdgeTable() {}

	EdgeTable(int matrixSize) :
		singleEdge(std::vector<std::vector<short>>(matrixSize)),
		doubleEdge(std::vector<std::vector<short>>(matrixSize)) {
		for (std::vector<short>& v : singleEdge) {
			v.reserve(4);
		}
		for (std::vector<short>& v : doubleEdge) {
			v.reserve(4);
		}
	}

	~EdgeTable() {}
};

class Algorithms {
public:
	void geneticAlgorithm(Matrix* matrix, int benchmark = 0);
	void benchmark(Matrix* matrix);
	void displayResults();

	void initRandom();

	void setBenchmarkFile(std::string fileName) { benchmarkFile = fileName; };
	void setStartingPopulationSize(int value);
	void setStartingPopulationRandom() 
		{ startingPopulationRandom = !startingPopulationRandom; };
	void setStopCriterium(int value);
	void setMutationConstant(double value);
	void setCrossoverConstant(double value) { crossoverConstant = value; };
	void setCrossoverType(int type)	
		{ type ? currentCrossoverType = 1 : currentCrossoverType = 0; };
	void setMutationType(MutationType type);

	std::string getBenchmarkFile() { return benchmarkFile; };
	std::chrono::duration<double> getMaxExecutionTime() { return maxExecutionTime; };
	int getStartingPopulationSize() { return startingPopulationSize; };
	bool getStartingPopulationRandom() { return startingPopulationRandom; };
	double getMutationConstant() { return mutationConstant; };
	double getCrossoverConstant() { return crossoverConstant; };
	int getCurrentCrossoverType() { return currentCrossoverType; };
	MutationType getCurrentMutationType() { return currentMutationType; };

private:
	std::string benchmarkFile = "BRAK";
	int pathLength = INT_MAX;
	std::vector<short> vertexOrder = {};

	std::chrono::duration<double> maxExecutionTime = std::chrono::seconds(1);
	std::chrono::duration<double> runningTime;

	int startingPopulationSize = 3; // min 3
	bool startingPopulationRandom = false; // true = random, false = greedy path
	double mutationConstant = 0.01f;
	double crossoverConstant = 0.8f;
	int currentCrossoverType = 0; // 0 OX, 1 EAX
	MutationType currentMutationType = INSERT;

	double repeat = 0;
	double total = 0;

	// Random
	std::random_device rd;
	std::mt19937 gen;
	std::tuple<int, int> generateRandomTwoPositions(int lowerBound, int higherBound, bool correctOrder = 1);

	std::tuple<std::vector<short>, int> generateInitialSolution(Matrix* matrix);
	std::tuple<std::vector<short>, int> generateNewSolution(Matrix* matrix, int notAllowedSecondary);
	std::tuple<std::vector<short>, int> generateNewSolutionV(Matrix* matrix, const short startVertex = 0);

	// Genetic
	void geneticOX(Matrix* matrix);
	void geneticEAX(Matrix* matrix);

	void geneticOXBench(Matrix* matrix, int iteration);
	void geneticEAXBench(Matrix* matrix, int iteration);

	std::vector<PathData> generateStartingPopulation(Matrix* matrix);
	void generateRandomStartingPopulation(Matrix* matrix, std::vector<short>* leftVertices, std::vector<PathData>* returnVec);
	void generateGreedyStartingPopulation(Matrix* matrix, std::vector<short>* leftVertices, std::vector<PathData>* returnVec);
	std::vector<double> getVertexLowerBounds(int vectorSize);
	std::vector<std::tuple<int, int>> generateParents(std::vector<double>* boundsVector, int vectorSize);
	// OX
	PathData generateChildOX(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent);
	// EAX
	PathData generateChildEAX(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent);
	
	void updateTable(Matrix* matrix, EdgeTable* edgeTable, std::vector<short>* firstParent, std::vector<short>* secondParent);
	void updateTableFind(EdgeTable* edgeTable, std::vector<short>* parent, int vertexToFind);
	void updateTableBig(EdgeTable* edgeTable, std::vector<short>* firstParent, std::vector<short>* secondParent);
	void deleteOccurences(EdgeTable* edgeTable, short currentVertex, std::vector<short>* vertexToDelete);
	short getNext(EdgeTable* edgeTable, short current, short currentFallback);

	std::vector<PathData> getNewRandomGeneration(Matrix* matrix, std::vector<PathData>* parents, std::vector<PathData>* children);

	PathData getNewOrder(std::vector<short>* currentOrder, int anchorOne, int anchorTwo, std::vector<std::vector<int>>* matrix);

	std::vector<short> inverse(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> swap(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> insert(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> insertSub(std::vector<short>* currentOrder);

	int calculateCandidate(std::vector<short>* candidateOrder, Matrix* matrix);
};

void clear();