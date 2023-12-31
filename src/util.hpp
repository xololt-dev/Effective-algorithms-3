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

struct QueueData {
	std::vector<short> pathOrder;
	int pathLength;
	int anchorOne;
	int anchorTwo;
};

class Algorithms {
public:
	void geneticAlgorithm(Matrix* matrix);
	void benchmark(Matrix* matrix);
	void displayResults();

	void initRandom();

	void setBenchmarkFile(std::string fileName) { benchmarkFile = fileName; };
	void setStartingPopulationSize(int value);
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

	// Random
	std::random_device rd;
	std::mt19937 gen;
	std::tuple<int, int> generateRandomTwoPositions(int lowerBound, int higherBound, bool correctOrder = 1);

	std::tuple<std::vector<short>, int> generateInitialSolution(Matrix* matrix);
	std::tuple<std::vector<short>, int> generateNewSolution(Matrix* matrix, int notAllowedSecondary);
	std::tuple<std::vector<short>, int> generateNewSolutionV(Matrix* matrix, const short startVertex = 0);

	// Genetic
	void geneticOX(Matrix* matrix);
	void newGeneticOX(Matrix* matrix);
	// void geneticEAX(Matrix* matrix);

	std::vector<QueueData> generateStartingPopulation(Matrix* matrix);
	std::vector<double> getVertexLowerBounds(int vectorSize);
	std::vector<std::tuple<short, short>> generateParents(std::vector<double>* boundsVector, int vectorSize);
	QueueData generateChild(Matrix* matrix, std::vector<short>* firstParent, std::vector<short>* secondParent);
	std::vector<QueueData> getNewRandomGeneration(Matrix* matrix, std::vector<QueueData>* parents, std::vector<QueueData>* children);

	QueueData getNewOrder(std::vector<short>* currentOrder, int anchorOne, int anchorTwo, std::vector<std::vector<int>>* matrix);

	std::vector<short> inverse(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> swap(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> insert(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> insertSub(std::vector<short>* currentOrder);

	int calculateCandidate(std::vector<short>* candidateOrder, Matrix* matrix);
};

void clear();