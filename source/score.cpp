/*F1Score is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  F1Score is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <list>
#include <windows.h>
#include <direct.h>

using namespace std;

std::vector<std::set<unsigned int>*>* partitionA;
std::vector<std::set<unsigned int>*>* partitionB;

#define	LOWER_THRESHOLD 2 

/*
* score.cpp用法
* 聚类结果A 聚类结果B 输出文件(注意是相对路径和文件名)
* test_data\\clusterA.txt test_data\\clusterB.txt data_out
*/
/** Parses the sets of a file **/
//从文件读取community
/* 输入文件格式
* 一行代表一个community
1 2 3 4 5
6 7 8 9 10
*/
std::vector<std::set<unsigned int>*>* ParseSets(std::ifstream& file) {
	std::list<std::set<unsigned int>*>* auxSets = new std::list<std::set<unsigned int>*>();
	std::string line;
	std::string delimiter(" ");
	// 每个line中存有一个community
	while (std::getline(file, line)) {
		std::set<unsigned int>* set = new std::set<unsigned int>();
		size_t nextToken = 0;
		std::istringstream stream(line);
		unsigned int node;
		while (stream >> node) {
			set->insert(node);
		}

		auxSets->push_back(set);
	}

	std::vector<std::set<unsigned int>*>* returnSet = new std::vector<std::set<unsigned int>*>();
	std::list<std::set<unsigned int>*>::iterator it;
	unsigned int i = 0;
	for (it = auxSets->begin(); it != auxSets->end(); ++it, ++i) {
		returnSet->push_back(*it);
	}
	return returnSet;
}

/*
计算两个community的F1-score
*/
double F1Score(std::set<unsigned int>& com1, std::set<unsigned int>& com2, double* precision, double* recall) {
	*precision = 0.0f;
	*recall = 0.0f;

	//computing intersection
	unsigned int counter = 0;
	if (com1.size() < com2.size()) {
		for (std::set<unsigned int>::iterator iterCom1 = com1.begin(); iterCom1 != com1.end(); iterCom1++) {
			if (com2.find(*iterCom1) != com2.end()) {
				counter++;
			}
		}
	}
	else {
		for (std::set<unsigned int>::iterator iterCom2 = com2.begin(); iterCom2 != com2.end(); iterCom2++) {
			if (com1.find(*iterCom2) != com1.end()) {
				counter++;
			}
		}
	}
	*precision = counter / (double)com1.size();
	*recall = counter / (double)com2.size();
	if (*precision + *recall > 0.0f) {
		return 2 * (*precision) * (*recall) / ((*precision) + (*recall));
	}
	return 0.0f;
}

static double Average(double array[], int size) {
	double accum = 0.0;
	for (int i = 0; i < size; ++i) {
		accum += array[i];
	}
	return accum / size;
}

// 格式，<输入文件1><输入文件2><输出文件1><输出文件2>
int socre(int argc, char** argv)
{
	if (argc < 4) {
		std::cout << "Wrong number of arguments. Usage: F1Score <partitionA> <partitionB> <datafilename>" << std::endl;
		exit(0);
	}

	char pwd[1000];
	_getcwd(pwd, 1000);
	cout << pwd << endl;

	std::ifstream inputFileA;
	char* tmp = argv[1];
	inputFileA.open(argv[1]);
	if (!inputFileA.is_open()) {
		std::cout << "PARTITION FILE A NOT FOUND" << std::endl;
		exit(1);
	}
	std::ifstream inputFileB;
	inputFileB.open(argv[2]);
	if (!inputFileB.is_open()) {
		std::cout << "PARTITION FILE B NOT FOUND" << std::endl;
		exit(1);
	}

	std::cout << "Parsing Input Files" << std::endl;
	partitionA = ParseSets(inputFileA);
	inputFileA.close();

	partitionB = ParseSets(inputFileB);
	inputFileB.close();

	/*
	* nodePartitionA
	* 这里记录每个node节点所在分区
	* 比如node 1: 2 3 5 6(这里指的是community的序号)
	*/
	std::map<unsigned int, std::set<unsigned int> > nodePartitionA;
	std::map<unsigned int, std::set<unsigned int> > nodePartitionB;

	std::cout << "Partition A size: " << partitionA->size() << std::endl;
	std::cout << "Partition B size: " << partitionB->size() << std::endl;
	double maxSetsF1Scores = 0.0f;
	double maxTagsF1Scores = 0.0f;

	std::cout << "Creating Indexes of Partition A" << std::endl;
	for (unsigned int i = 0; i < partitionA->size(); i++) {
		std::set<unsigned int>* community = (*partitionA)[i];
		for (auto it = community->begin(); it != community->end(); it++) {
			unsigned int node = *it;
			auto it2 = nodePartitionA.find(node);
			nodePartitionA[node].insert(i);
		}
	}

	std::cout << "Creating Indexes of Partition B" << std::endl;
	for (unsigned int i = 0; i < partitionB->size(); i++) {
		std::set<unsigned int>* community = (*partitionB)[i];
		for (auto it = community->begin(); it != community->end(); it++) {
			unsigned int node = *it;
			auto it2 = nodePartitionB.find(node);
			nodePartitionB[node].insert(i);
		}
	}

	double* f1ScorePartitionA = new double[partitionA->size()];
	double* f1ScorePartitionB = new double[partitionB->size()];

	double* f1ScorePrecisionPartitionA = new double[partitionA->size()];
	double* f1ScorePrecisionPartitionB = new double[partitionB->size()];

	double* f1ScoreRecallPartitionA = new double[partitionA->size()];
	double* f1ScoreRecallPartitionB = new double[partitionB->size()];

	int* bestIdPartitionA = new int[partitionA->size()];
	int* bestIdPartitionB = new int[partitionB->size()];

	double* precisionPartitionA = new double[partitionA->size()];
	double* precisionPartitionB = new double[partitionB->size()];

	double* recallPartitionA = new double[partitionA->size()];
	double* recallPartitionB = new double[partitionB->size()];

	int* coveragePartitionA = new int[partitionA->size()];
	int* coveragePartitionB = new int[partitionB->size()];


	/** Initializing arrays to 0 **/
	std::memset(f1ScorePartitionA, 0, sizeof(double) * partitionA->size());
	std::memset(f1ScorePrecisionPartitionA, 0, sizeof(double) * partitionA->size());
	std::memset(f1ScoreRecallPartitionA, 0, sizeof(double) * partitionA->size());
	std::memset(precisionPartitionA, 0, sizeof(double) * partitionA->size());
	std::memset(recallPartitionA, 0, sizeof(double) * partitionA->size());
	std::memset(coveragePartitionA, 0, sizeof(int) * partitionA->size());
	std::memset(bestIdPartitionA, 0, sizeof(int) * partitionA->size());

	std::memset(f1ScorePartitionB, 0, sizeof(double) * partitionB->size());
	std::memset(f1ScorePrecisionPartitionB, 0, sizeof(double) * partitionB->size());
	std::memset(f1ScoreRecallPartitionB, 0, sizeof(double) * partitionB->size());
	std::memset(precisionPartitionB, 0, sizeof(double) * partitionB->size());
	std::memset(recallPartitionB, 0, sizeof(double) * partitionB->size());
	std::memset(coveragePartitionB, 0, sizeof(int) * partitionB->size());
	std::memset(bestIdPartitionB, 0, sizeof(int) * partitionB->size());

	std::cout << "Computing F1Score of partition A" << std::endl;
	for (unsigned int i = 0; i < partitionA->size(); i++) {
		std::set<unsigned int>* community = (*partitionA)[i];
		std::set<unsigned int> communities;

		// Selecting candidate communities to compare with.
		// coveragePartitionA
		for (auto it = community->begin(); it != community->end(); it++) {
			// it: partitionA[i][iter]--某个节点
			auto it2 = nodePartitionB.find(*it);
			if (it2 != nodePartitionB.end()) {
				// it2: nodePartition[it]--包含顶点it的全部commmunity，
				for (auto it3 = (*it2).second.begin(); it3 != (*it2).second.end(); it3++) {
					communities.insert(*it3);
				}
				coveragePartitionA[i]++;
			}
		}

		// communities存放于当前迭代的community[i]有交集的所有commnity序号
		for (auto it = communities.begin(); it != communities.end(); it++) {
			std::set<unsigned int>* community2 = (*partitionB)[*it];
			double precision;
			double recall;
			// 记录最好的f1-score相关信息
			double f1Score = F1Score(*community, *community2, &precision, &recall);
			if (f1Score >= f1ScorePartitionA[i]) {
				f1ScorePartitionA[i] = f1Score;
				f1ScorePrecisionPartitionA[i] = precision;
				f1ScoreRecallPartitionA[i] = recall;
				bestIdPartitionA[i] = *it;
			}
			// 记录最好的precision和recall相关信息
			precisionPartitionA[i] = precision > precisionPartitionA[i] ? precision : precisionPartitionA[i];
			recallPartitionA[i] = recall > recallPartitionA[i] ? recall : recallPartitionA[i];
		}
	}

	std::cout << "Computing F1Score of partition B" << std::endl;
	for (unsigned int i = 0; i < partitionB->size(); i++) {
		std::set<unsigned int>* community = (*partitionB)[i];
		std::set<unsigned int> communities;

		// Selecting candidate communities to compare with.
		for (auto it = community->begin(); it != community->end(); it++) {
			auto it2 = nodePartitionA.find(*it);
			if (it2 != nodePartitionA.end()) {
				for (auto it3 = (*it2).second.begin(); it3 != (*it2).second.end(); it3++) {
					communities.insert(*it3);
				}
				coveragePartitionB[i]++;
			}
		}

		for (auto it = communities.begin(); it != communities.end(); it++) {
			std::set<unsigned int>* community2 = (*partitionA)[*it];
			double precision;
			double recall;
			double f1Score = F1Score(*community, *community2, &precision, &recall);
			if (f1Score >= f1ScorePartitionB[i]) {
				f1ScorePartitionB[i] = f1Score;
				f1ScorePrecisionPartitionB[i] = precision;
				f1ScoreRecallPartitionB[i] = recall;
				bestIdPartitionB[i] = *it;
			}
			precisionPartitionB[i] = precision > precisionPartitionB[i] ? precision : precisionPartitionB[i];
			recallPartitionB[i] = recall > recallPartitionB[i] ? recall : recallPartitionB[i];
		}
	}

	std::cout << "Average precision partition A: " << Average(precisionPartitionA, partitionA->size()) << std::endl;
	std::cout << "Average recall partition A: " << Average(recallPartitionA, partitionA->size()) << std::endl;

	std::cout << "Average precision partition B: " << Average(precisionPartitionB, partitionB->size()) << std::endl;
	std::cout << "Average recall partition B: " << Average(recallPartitionB, partitionB->size()) << std::endl;

	// 这里计算的F1-score就是average F1-score
	std::cout << "F1Score: " << (Average(f1ScorePartitionA, partitionA->size()) + Average(f1ScorePartitionB, partitionB->size())) / 2 << std::endl;

	std::ofstream outputFileA;
	outputFileA.open(std::string(argv[3]).append("-partA").c_str() + string(".txt"));
	outputFileA << "id|size|precision|recall|f1score|f1scorePrecision|f1scoreRecall|coverage|id2|other_coverage\n";
	for (int i = 0; i < partitionA->size(); ++i) {
		std::set<unsigned int>* community = (*partitionA)[i];
		outputFileA << i << "|" << community->size() << "|" << precisionPartitionA[i] << "|" << recallPartitionA[i] << "|" << f1ScorePartitionA[i] << "|" << f1ScorePrecisionPartitionA[i] << "|" << f1ScoreRecallPartitionA[i] << "|" << coveragePartitionA[i] / (double)community->size() << "|" << bestIdPartitionA[i] << "|" << coveragePartitionB[bestIdPartitionA[i]] / (double)(*partitionB)[bestIdPartitionA[i]]->size() << "\n";
	}
	outputFileA.close();

	std::ofstream outputFileB;
	outputFileB.open(std::string(argv[3]).append("-partB").c_str() + string(".txt"));
	outputFileB << "id|size|precision|recall|f1score|f1scorePrecision|f1scoreRecall|coverage|id2|other_coverage\n";
	for (int i = 0; i < partitionB->size(); ++i) {
		std::set<unsigned int>* community = (*partitionB)[i];
		outputFileB << i << "|" << community->size() << "|" << precisionPartitionB[i] << "|" << recallPartitionB[i] << "|" << f1ScorePartitionB[i] << "|" << f1ScorePrecisionPartitionB[i] << "|" << f1ScoreRecallPartitionB[i] << "|" << coveragePartitionB[i] / (double)community->size() << "|" << bestIdPartitionB[i] << "|" << coveragePartitionA[bestIdPartitionB[i]] / (double)(*partitionA)[bestIdPartitionB[i]]->size() << "\n";
	}
	outputFileB.close();

	delete[] f1ScorePartitionA;
	delete[] f1ScorePartitionB;

	delete[] f1ScorePrecisionPartitionA;
	delete[] f1ScorePrecisionPartitionB;

	delete[] f1ScoreRecallPartitionA;
	delete[] f1ScoreRecallPartitionB;

	delete[] precisionPartitionA;
	delete[] precisionPartitionB;

	delete[] recallPartitionA;
	delete[] recallPartitionB;

	delete[] coveragePartitionA;
	delete[] coveragePartitionB;

	delete[] bestIdPartitionA;
	delete[] bestIdPartitionB;

	return 0;
}