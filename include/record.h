#pragma once
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <time.h>
#include <fstream>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include "types.h"
#include "utils.h"
#include <stdlib.h>
#include <direct.h>
#include <iostream>
#include <vector>
#include <queue>
using namespace std;

// ����node-community��¼
struct NodeCommunityRecord
{
    uint32_t node_id;
    uint32_t com_id;
    set<uint32_t> neighbor;
    NodeCommunityRecord() {} // Ĭ�Ϲ��캯��
    NodeCommunityRecord(uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor) {
        this->node_id = node_id;
        this->com_id = com_id;
        this->neighbor = neighbor;
    }
    bool operator<(const NodeCommunityRecord& other) const {
        return neighbor.size() < other.neighbor.size();//�󶥶�
    }
    // ����vector<NodeCommunityRecord>��find����
    bool operator == (uint32_t com_id) {
        return (this->com_id == com_id);
    }
};

// �����ھӣ�community��¼,�����㷨2
struct NodeNeighborCommunity
{
    set<uint32_t>neighbor;
    set<uint32_t>community;
    uint32_t degree;
    NodeNeighborCommunity() { degree = 0; } // Ĭ�Ϲ��캯��
    void update(uint32_t adjNode, uint32_t com_id) {
        neighbor.insert(adjNode);
        community.insert(com_id);
        degree++;
    }
};

class Record
{
public:
    Record(std::vector< std::vector<std::multiset<NodeCommunityRecord>>>& nodeCommunityList,
        std::vector< std::vector< vector<Edge> > >& communityEdgeList,
        std::vector< uint32_t >& maxVolumeList);//���캯��
    bool IsAvailable(uint32_t threshold_id, uint32_t node_id); //�ö����Ƿ��пɷ���(�����������������)��community
    void edgeInsert(uint32_t threshold_id, uint32_t com_id, Edge edge); //�߲��뵽community��
    void edgeCreateNewCommunity(uint32_t threshold_id, Edge edge); //ĳ���½�һ��community
    void insertNodeComRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor);  //����һ��Node-Community�ļ�¼
    void updateNodeComRecord_AddNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode);  //����Node-Community�ļ�¼��������µ��Ǵ󶥶ѵļ�¼
    void updateNodeComRecord_AddMulNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor);  //����Node-Community�ļ�¼��һ������Ӷ������
    void updateNodeComRecord_DeleteNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode);  //����Node-Community�ļ�¼������ɾ��ĳһ���ھ�
    NodeCommunityRecord getTopRecord(uint32_t threshold_id, uint32_t node_id);
    multiset<NodeCommunityRecord>::iterator findRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id);
    void edgeTrans(uint32_t threshold_id, uint32_t min_id, uint32_t max_id); // ��ת�ƺ���
    void print(const char* fileName, uint32_t threshold_id, bool removeSingleton); // ��ӡ���

    // ��Ա����
    std::vector< std::vector<std::multiset<NodeCommunityRecord>>> nodeCommunityList;
    std::vector< std::vector< vector<Edge>>> communityEdgeList;
    std::vector< uint32_t > maxVolumeList;
};

// �㷨2����ʵ��
// �����ڵ��ھӵ�ͬʱ��Ϊÿ��community�������
void updateCommunityScore(uint32_t threshold_id, vector<vector<NodeNeighborCommunity>>& record, map<int, double>& com_score, uint32_t node_id);

// �㷨3���ƶȵļ���
double cal_sim(uint32_t threshold_id, vector<vector<map<uint32_t, set<uint32_t>>>>& record, uint32_t sourceNode, uint32_t targetNode, vector<uint32_t>& vertexDegree, map< uint32_t, set<uint32_t>>::iterator src_it, map< uint32_t, set<uint32_t>>::iterator tar_it);