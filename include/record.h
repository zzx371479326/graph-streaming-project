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

// 定义node-community记录
struct NodeCommunityRecord
{
    uint32_t node_id;
    uint32_t com_id;
    set<uint32_t> neighbor;
    NodeCommunityRecord() {} // 默认构造函数
    NodeCommunityRecord(uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor) {
        this->node_id = node_id;
        this->com_id = com_id;
        this->neighbor = neighbor;
    }
    bool operator<(const NodeCommunityRecord& other) const {
        return neighbor.size() < other.neighbor.size();//大顶堆
    }
    // 用于vector<NodeCommunityRecord>的find函数
    bool operator == (uint32_t com_id) {
        return (this->com_id == com_id);
    }
};

// 定义邻居，community记录,用于算法2
struct NodeNeighborCommunity
{
    set<uint32_t>neighbor;
    set<uint32_t>community;
    uint32_t degree;
    NodeNeighborCommunity() { degree = 0; } // 默认构造函数
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
        std::vector< uint32_t >& maxVolumeList);//构造函数
    bool IsAvailable(uint32_t threshold_id, uint32_t node_id); //该顶点是否有可分配(不超过最大容量限制)的community
    void edgeInsert(uint32_t threshold_id, uint32_t com_id, Edge edge); //边插入到community中
    void edgeCreateNewCommunity(uint32_t threshold_id, Edge edge); //某边新建一个community
    void insertNodeComRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor);  //插入一条Node-Community的记录
    void updateNodeComRecord_AddNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode);  //更新Node-Community的记录，这里更新的是大顶堆的记录
    void updateNodeComRecord_AddMulNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor);  //更新Node-Community的记录，一次性添加多个顶点
    void updateNodeComRecord_DeleteNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode);  //更新Node-Community的记录，这里删除某一个邻居
    NodeCommunityRecord getTopRecord(uint32_t threshold_id, uint32_t node_id);
    multiset<NodeCommunityRecord>::iterator findRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id);
    void edgeTrans(uint32_t threshold_id, uint32_t min_id, uint32_t max_id); // 边转移函数
    void print(const char* fileName, uint32_t threshold_id, bool removeSingleton); // 打印结果

    // 成员变量
    std::vector< std::vector<std::multiset<NodeCommunityRecord>>> nodeCommunityList;
    std::vector< std::vector< vector<Edge>>> communityEdgeList;
    std::vector< uint32_t > maxVolumeList;
};

// 算法2函数实现
// 遍历节点邻居的同时，为每个community计算分数
void updateCommunityScore(uint32_t threshold_id, vector<vector<NodeNeighborCommunity>>& record, map<int, double>& com_score, uint32_t node_id);

// 算法3相似度的计算
double cal_sim(uint32_t threshold_id, vector<vector<map<uint32_t, set<uint32_t>>>>& record, uint32_t sourceNode, uint32_t targetNode, vector<uint32_t>& vertexDegree, map< uint32_t, set<uint32_t>>::iterator src_it, map< uint32_t, set<uint32_t>>::iterator tar_it);