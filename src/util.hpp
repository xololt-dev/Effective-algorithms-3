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

	void loadFromFile(std::string fileName);
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
	void setStopCriterium(int value);
	void setMutationConstant(double value);
	void setCrossoverConstant(double value);
	void setStartingPopulationSize(int value);
	void setMutationType(MutationType type);

	void initRandom();

	void geneticAlgorithm(Matrix* matrix);
	
	void benchmark(Matrix* matrix);

	void displayResults();

private:
	std::string benchmarkFile;
	int pathLength;
	std::vector<short> vertexOrder;

	std::chrono::duration<double> maxExecutionTime = std::chrono::seconds(30);
	std::chrono::duration<double> runningTime;

	int startingPopulationSize = 10;
	bool startingPopulationRandom = false; // true = random, false = greedy path
	double mutationConstant = 0.99f;
	double crossoverConstant = 0.99f;
	MutationType currentMutationType = INVERSE;

	// Random
	std::random_device rd;
	std::mt19937 gen;
	std::tuple<int, int> generateRandomTwoPositions(int lowerBound, int higherBound, bool correctOrder = 1);

	std::tuple<std::vector<short>, int> generateInitialSolution(Matrix* matrix);
	std::tuple<std::vector<short>, int> generateNewSolution(Matrix* matrix, int notAllowedSecondary);
	std::tuple<std::vector<short>, int> generateNewSolutionV(Matrix* matrix, const short startVertex = 0);

	// Genetic
	void geneticOX(Matrix* matrix);
	// void geneticEAX(Matrix* matrix);
	std::vector<QueueData> generateStartingPopulation(Matrix* matrix);
	
	// TS
	QueueData generateBestMove(std::vector<short>* currentOrder, int bestLength, std::vector<std::vector<int>>* matrix);
	QueueData getNewOrder(std::vector<short>* currentOrder, int anchorOne, int anchorTwo, std::vector<std::vector<int>>* matrix);

	// SA
	int getPathDelta(Matrix* matrix);

	std::vector<short> generateRandomCandidate(std::vector<short>* currentOrder);
	std::vector<short> inverse(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> swap(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> insert(std::vector<short>* currentOrder, int firstPosition = 0, int secondPosition = 0);
	std::vector<short> insertSub(std::vector<short>* currentOrder);

	bool changeSolutions(int candidatePath, int currentPath, double currentTemp);
	int calculateCandidate(std::vector<short>* candidateOrder, Matrix* matrix);
};

void clear();