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
#include <time.h>
#include "utils.h"

class Edge
{
public:
	Edge(uint32_t u, uint32_t v)
		: first(u), second(v)
	{
	}
	bool operator < (const Edge& edge) const
	{
		if (edge.first == this->first)
			return this->second < edge.second;
		else
			return this->first < edge.first;
	}
	Edge edge_reverse() { return Edge(this->second, this->first); }

	//Ϊ�˷��������Ϊpublic
	uint32_t first;
	uint32_t second;
};

using namespace std;

std::vector<std::set<Edge>*>* partitionA;
std::vector<std::set<Edge>*>* partitionB;

#define	LOWER_THRESHOLD 2 

/*
* score.cpp�÷�
* ������A ������B ����ļ�(ע�������·�����ļ���)
* test_data\\clusterA.txt test_data\\clusterB.txt data_out
*/
/** Parses the sets of a file **/
//���ļ���ȡcommunity
/* �����ļ���ʽ
* һ�д���һ��community
1 2 3 4 5
6 7 8 9 10
*/
std::vector<std::set<Edge>*>* ParseSets(std::ifstream& file) {
	std::list<std::set<Edge>*>* auxSets = new std::list<std::set<Edge>*>();
	std::string line;
	std::string delimiter(" ");
	// ÿ��line�д���һ��community
	while (std::getline(file, line)) {
		std::set<Edge>* set = new std::set<Edge>();
		size_t nextToken = 0;
		std::istringstream stream(line);
		unsigned int sourceNode, targetNode;
		while (stream >> sourceNode) {
			stream >> targetNode;
			set->insert(Edge(sourceNode, targetNode));
		}
		if(!set->empty())
			auxSets->push_back(set);
	}

	std::vector<std::set<Edge>*>* returnSet = new std::vector<std::set<Edge>*>();
	std::list<std::set<Edge>*>::iterator it;
	unsigned int i = 0;
	for (it = auxSets->begin(); it != auxSets->end(); ++it, ++i) {
		returnSet->push_back(*it);
	}
	return returnSet;
}

/*
��������community��F1-score
*/
double F1Score(std::set<Edge>& com1, std::set<Edge>& com2, double* precision, double* recall) {
	*precision = 0.0f;
	*recall = 0.0f;

	//computing intersection
	unsigned int counter = 0;
	if (com1.size() < com2.size()) {
		for (std::set<Edge>::iterator iterCom1 = com1.begin(); iterCom1 != com1.end(); iterCom1++) {
			// �������ֿ���,��Ҫreverse
			if (com2.find(*iterCom1) != com2.end() || com2.find( Edge(iterCom1->second, iterCom1->first)) != com2.end()) {
				counter++;
			}
		}
	}
	else {
		for (std::set<Edge>::iterator iterCom2 = com2.begin(); iterCom2 != com2.end(); iterCom2++) {
			if (com1.find(*iterCom2) != com1.end() || com1.find( Edge(iterCom2->second, iterCom2->first)) != com1.end()) {
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

// ��ʽ��<�����ļ�1><�����ļ�2><����ļ�1><����ļ�2>
int main(int argc, char** argv)
{
	long unsigned initTime = StartClock();
	if (argc < 4) {
		std::cout << "Wrong number of arguments. Usage: F1Score <partitionA> <partitionB> <datafilename>" << std::endl;
		exit(0);
	}

	/*char pwd[1000];
	_getcwd(pwd, 1000);
	cout << pwd << endl;*/

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

	// �߷�������������
	/*
	* nodePartitionA
	* �����¼ÿ��node�ڵ����ڷ���
	* ����node 1: 2 3 5 6(����ָ����community�����)
	*/
	std::map<Edge, std::set<unsigned int> > nodePartitionA;
	std::map<Edge, std::set<unsigned int> > nodePartitionB;

	std::cout << "Partition A size: " << partitionA->size() << std::endl;
	std::cout << "Partition B size: " << partitionB->size() << std::endl;
	double maxSetsF1Scores = 0.0f;
	double maxTagsF1Scores = 0.0f;

	std::cout << "Creating Indexes of Partition A" << std::endl;
	for (unsigned int i = 0; i < partitionA->size(); i++) {
		std::set<Edge>* community = (*partitionA)[i];
		for (auto it = community->begin(); it != community->end(); it++) {
			// Edge edge = *it;
			// auto it2 = nodePartitionA.find(node);
			nodePartitionA[*it].insert(i);
		}
	}

	std::cout << "Creating Indexes of Partition B" << std::endl;
	for (unsigned int i = 0; i < partitionB->size(); i++) {
		std::set<Edge>* community = (*partitionB)[i];
		for (auto it = community->begin(); it != community->end(); it++) {
			// unsigned int node = *it;
			// auto it2 = nodePartitionB.find(node);
			nodePartitionB[*it].insert(i);
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
		std::set<Edge>* community = (*partitionA)[i];
		std::set<unsigned int> communities; // ÿ��communityA��communityB�а����˶��ٸ�

		// Selecting candidate communities to compare with.
		// coveragePartitionA
		for (auto it = community->begin(); it != community->end(); it++) {
			// it: partitionA[i][iter]--ĳ���ڵ�
			auto it2 = nodePartitionB.find(*it);
			if (it2 != nodePartitionB.end()) {
				// it2: nodePartition[it]--��������it��ȫ��commmunity��
				for (auto it3 = (*it2).second.begin(); it3 != (*it2).second.end(); it3++) {
					communities.insert(*it3);
				}
				coveragePartitionA[i]++;
			}
		}

		// communities����ڵ�ǰ������community[i]�н���������commnity���
		for (auto it = communities.begin(); it != communities.end(); it++) {
			std::set<Edge>* community2 = (*partitionB)[*it];
			double precision;
			double recall;
			// ��¼��õ�f1-score�����Ϣ
			double f1Score = F1Score(*community, *community2, &precision, &recall);
			if (f1Score >= f1ScorePartitionA[i]) {
				f1ScorePartitionA[i] = f1Score;
				f1ScorePrecisionPartitionA[i] = precision;
				f1ScoreRecallPartitionA[i] = recall;
				bestIdPartitionA[i] = *it;
			}
			// ��¼��õ�precision��recall�����Ϣ
			precisionPartitionA[i] = precision > precisionPartitionA[i] ? precision : precisionPartitionA[i];
			recallPartitionA[i] = recall > recallPartitionA[i] ? recall : recallPartitionA[i];
		}
	}

	std::cout << "Computing F1Score of partition B" << std::endl;
	for (unsigned int i = 0; i < partitionB->size(); i++) {
		std::set<Edge>* community = (*partitionB)[i];
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
			std::set<Edge>* community2 = (*partitionA)[*it];
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

	// ��������F1-score����average F1-score
	std::cout << "F1Score: " << (Average(f1ScorePartitionA, partitionA->size()) + Average(f1ScorePartitionB, partitionB->size())) / 2 << std::endl;

	std::ofstream outputFileA;
	outputFileA.open(std::string(argv[3]).append("-partA").c_str() + string(".txt"));
	outputFileA << "id|size|precision|recall|f1score|f1scorePrecision|f1scoreRecall|coverage|id2|other_coverage\n";
	for (int i = 0; i < partitionA->size(); ++i) {
		std::set<Edge>* community = (*partitionA)[i];
		outputFileA << i << "|" << community->size() << "|" << precisionPartitionA[i] << "|" << recallPartitionA[i] << "|" << f1ScorePartitionA[i] << "|" << f1ScorePrecisionPartitionA[i] << "|" << f1ScoreRecallPartitionA[i] << "|" << coveragePartitionA[i] / (double)community->size() << "|" << bestIdPartitionA[i] << "|" << coveragePartitionB[bestIdPartitionA[i]] / (double)(*partitionB)[bestIdPartitionA[i]]->size() << "\n";
	}
	outputFileA.close();

	std::ofstream outputFileB;
	outputFileB.open(std::string(argv[3]).append("-partB").c_str() + string(".txt"));
	outputFileB << "id|size|precision|recall|f1score|f1scorePrecision|f1scoreRecall|coverage|id2|other_coverage\n";
	for (int i = 0; i < partitionB->size(); ++i) {
		std::set<Edge>* community = (*partitionB)[i];
		outputFileB << i << "|" << community->size() << "|" << precisionPartitionB[i] << "|" << recallPartitionB[i] << "|" << f1ScorePartitionB[i] << "|" << f1ScorePrecisionPartitionB[i] << "|" << f1ScoreRecallPartitionB[i] << "|" << coveragePartitionB[i] / (double)community->size() << "|" << bestIdPartitionB[i] << "|" << coveragePartitionA[bestIdPartitionB[i]] / (double)(*partitionA)[bestIdPartitionB[i]]->size() << "\n";
	}
	outputFileB.close();
	long unsigned spentTime = StopClock(initTime);
	printf("%-32s %lu ms\n", "Compute time:", spentTime);

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