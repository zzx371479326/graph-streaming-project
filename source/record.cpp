#include"record.h"

Record::Record(std::vector< std::vector<std::multiset<NodeCommunityRecord>>>& nodeCommunityList,
    std::vector< std::vector< vector<Edge> > >& communityEdgeList,
    std::vector< uint32_t >& maxVolumeList){
    this->nodeCommunityList = nodeCommunityList;
    this->communityEdgeList = communityEdgeList;
    this->maxVolumeList = maxVolumeList;
}

bool Record::IsAvailable(uint32_t threshold_id, uint32_t node_id) {
    // 如果该节点从未被分配过
    if (nodeCommunityList[threshold_id][node_id].size() == 0)
        return false;
    // 如果该节点所在全部community容量满了
    multiset<NodeCommunityRecord>::iterator it = nodeCommunityList[threshold_id][node_id].begin();
    while (it != nodeCommunityList[threshold_id][node_id].end()) {
        if (communityEdgeList[threshold_id][it->com_id].size() >= maxVolumeList[threshold_id] || communityEdgeList[threshold_id][it->com_id].size() <= 0) {
            it = nodeCommunityList[threshold_id][node_id].erase(it);
            if (it == nodeCommunityList[threshold_id][node_id].end())
                break;
            continue;
        }
        else
            it++;
    }
    //for(;it != nodeCommunityList[threshold_id][node_id].end();it++){
    //    // 如果扫描的时候该分区满了，则将该分区的信息从列表中删除
    //    if (communityEdgeList[threshold_id][it->com_id].size() >= maxVolumeList[threshold_id] || communityEdgeList[threshold_id][it->com_id].size() <= 0) {
    //        it = nodeCommunityList[threshold_id][node_id].erase(it);
    //        if (it == nodeCommunityList[threshold_id][node_id].end())
    //            break;
    //        continue;
    //    }
    //}
    // 如果扫描完还有community的记录，则说明仍然有可用的分区
    if (nodeCommunityList[threshold_id][node_id].size() == 0)
        return false;
    return true;
}

void Record::edgeInsert(uint32_t threshold_id, uint32_t com_id, Edge edge) {
    communityEdgeList[threshold_id][com_id].push_back(edge);
}

void Record::edgeCreateNewCommunity(uint32_t threshold_id, Edge edge) {
    // 新建一个commmunity
    communityEdgeList[threshold_id].push_back(vector<Edge>{edge});
    uint32_t com_id = communityEdgeList[threshold_id].size() - 1; // 标号从0开始
    // 记录信息
    insertNodeComRecord(threshold_id, edge.first, com_id, set<uint32_t>{edge.second});
    insertNodeComRecord(threshold_id, edge.second, com_id, set<uint32_t>{edge.first});
}

void Record::insertNodeComRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor) {
    nodeCommunityList[threshold_id][node_id].insert(NodeCommunityRecord(node_id, com_id, neighbor));
}

void Record::updateNodeComRecord_AddNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode) {
    // 向node_id的记录中添加otherNode
    // 必须先删除，再insert，不能直接修改
    multiset<NodeCommunityRecord>::iterator it;
    it = findRecord(threshold_id, node_id, com_id);
    // 判断是否有com_id的记录
    if (it != nodeCommunityList[threshold_id][node_id].end()) {
        NodeCommunityRecord tmp = *it;
        tmp.neighbor.insert(otherNode);
        nodeCommunityList[threshold_id][node_id].erase(it);
        nodeCommunityList[threshold_id][node_id].insert(tmp);
    }
    // 没有的话直接创建新的即可
    else {
        insertNodeComRecord(threshold_id, node_id, com_id, set<uint32_t>{otherNode});
    }
}

void Record::updateNodeComRecord_AddMulNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor) {
    // 向node_id的记录中添加otherNode
    // 必须先删除，再insert，不能直接修改
    multiset<NodeCommunityRecord>::iterator it;
    it = findRecord(threshold_id, node_id, com_id);
    // 判断是否有com_id的记录
    if (it != nodeCommunityList[threshold_id][node_id].end()) {
        NodeCommunityRecord tmp = *it;
        tmp.neighbor.insert(neighbor.begin(), neighbor.end());
        nodeCommunityList[threshold_id][node_id].erase(it);
        nodeCommunityList[threshold_id][node_id].insert(tmp);
    }
    // 没有的话直接创建新的即可
    else {
        insertNodeComRecord(threshold_id, node_id, com_id, neighbor);
    }
}

void Record::updateNodeComRecord_DeleteNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode) {
    multiset<NodeCommunityRecord>::iterator it = findRecord(threshold_id, node_id, com_id);
    // 以下语句出现了bug
    if (it == nodeCommunityList[threshold_id][node_id].end())
        return;
    // 如果neighbor只有一个，则直接删除该记录
    NodeCommunityRecord tmp = *it;
    nodeCommunityList[threshold_id][node_id].erase(it);
    if (tmp.neighbor.size() > 1){
        tmp.neighbor.erase(otherNode);
        nodeCommunityList[threshold_id][node_id].insert(tmp);
    }
}

NodeCommunityRecord Record::getTopRecord(uint32_t threshold_id, uint32_t node_id) {
    // 返回最大容量的community
    return *nodeCommunityList[threshold_id][node_id].begin();
}

multiset<NodeCommunityRecord>::iterator Record::findRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id){
    multiset<NodeCommunityRecord>::iterator it;
    for (it = nodeCommunityList[threshold_id][node_id].begin(); it != nodeCommunityList[threshold_id][node_id].end(); it++)
        if (it->com_id == com_id)
            break;
    return it;
}

void Record::edgeTrans(uint32_t threshold_id, uint32_t min_id, uint32_t max_id) {
    // 找出各自community
    uint32_t min_com = nodeCommunityList[threshold_id][min_id].begin()->com_id;
    uint32_t max_com = nodeCommunityList[threshold_id][max_id].begin()->com_id;
    // 对communityEdgeList进行遍历
    vector<Edge>::iterator it;
    for (it = communityEdgeList[threshold_id][min_com].begin(); it != communityEdgeList[threshold_id][min_com].end(); it++) {
        // 如果该边含有min_id，则执行移除和插入操作
        if (it->first == min_id || it->second == min_id) {
            uint32_t otherNode = it->first == min_id ? it->second : it->first;
            // 先执行边的插入删除操作
            edgeInsert(threshold_id, max_com, *it);
            it = communityEdgeList[threshold_id][min_com].erase(it); // 边删除操作
            // 更新ohter节点的记录
            // otherNode记录处理...
            updateNodeComRecord_DeleteNode(threshold_id, otherNode, min_com, min_id); // 删除记录(有Bug--边 (763， 465909)会无法删除记录，因为记录已经被删除了？？
            updateNodeComRecord_AddNode(threshold_id, otherNode, max_com, min_id); // 添加记录
            if (it == communityEdgeList[threshold_id][min_com].end())
                break;
            continue;
        }
    }
    // min_id记录处理...
    // 找出min_node的min_com记录指针
    multiset<NodeCommunityRecord>::iterator min_com_it = findRecord(threshold_id, min_id, min_com);
    // 向max_com的记录中添加
    updateNodeComRecord_AddMulNode(threshold_id, min_id, max_com, min_com_it->neighbor);
    // 删除min_id的min_com记录
    nodeCommunityList[threshold_id][min_id].erase(min_com_it);
}

void Record::print(const char* fileName, uint32_t threshold_id, bool removeSingleton) {
    std::ofstream outFile;
    outFile.open(fileName);
    uint32_t nbCommunities = 0;
    for (auto kv : communityEdgeList[threshold_id]) {
        if (!removeSingleton || kv.size() > 1) {
            vector<Edge>::iterator it2 = kv.begin();
            while (it2 != kv.end()) {
                outFile << it2->first << ' ' << it2->second << ' ';
                it2++;
            }
            outFile << endl;
        }
    }
    printf("%-32s %i\n", "Number of Communities:", communityEdgeList[threshold_id].size());
    outFile.close();
}

/*
// 算法2函数实现
// 遍历节点邻居的同时，为每个community计算分数
void updateCommunityScore(uint32_t threshold_id, vector<vector<NodeNeighborCommunity>>& record, map<int, double>& com_score, uint32_t node_id) {
    // 遍历某个节点的所有邻居
    for (set<uint32_t>::iterator it_node = record[threshold_id][node_id].neighbor.begin(); it_node != record[threshold_id][node_id].neighbor.end(); it_node++) {
        // 遍历该节点的所有community
        for (set<uint32_t>::iterator it_com = record[threshold_id][*it_node].community.begin(); it_com != record[threshold_id][*it_node].community.end(); it_com++) {
            map<int, double>::iterator com_score_it = com_score.find(*it_com);
            if (com_score_it != com_score.end()) {
                com_score[*it_com] += 1 / (record[threshold_id][*it_node].degree + 1);
            }
            else
                com_score[*it_com] = 1 / (record[threshold_id][*it_node].degree + 1);
        }
    }
}

// 算法3函数实现
// 算法3相似度的计算
double cal_sim(uint32_t threshold_id, vector<vector<map<uint32_t, set<uint32_t>>>>& record, uint32_t sourceNode, uint32_t targetNode, vector<uint32_t>& vertexDegree, map< uint32_t, set<uint32_t>>::iterator src_it, map< uint32_t, set<uint32_t>>::iterator tar_it) {
    uint32_t com_id = src_it->first;
    double sim = 0;
    // 提取出target在该com的邻居
    set<uint32_t>tar_neighbor = tar_it->second;
    // 提取出source在该com的邻居
    set<uint32_t>src_neighbor = src_it->second;
    // 先计算src的所有邻居
    set<uint32_t>::iterator set_it;
    for (set_it = src_neighbor.begin(); set_it != src_neighbor.end(); set_it++) {
        // 找出该src_neighbor在该community的邻边集合
        set<uint32_t>src_neighbor_set = record[threshold_id][*set_it][com_id];
        int result_size = tar_neighbor.size() > src_neighbor_set.size() ? tar_neighbor.size() : src_neighbor_set.size();
        vector<uint32_t>insection_set(result_size); //结果集要用vector存放,并且初始化大小
        set_intersection(tar_neighbor.begin(), tar_neighbor.end(), src_neighbor_set.begin(), src_neighbor_set.end(), insection_set.begin());
        sim += (insection_set.size()) / (vertexDegree[targetNode] + vertexDegree[*set_it]);
    }
    // 再计算tar所有邻居
    for (set_it = tar_neighbor.begin(); set_it != tar_neighbor.end(); set_it++) {
        // 找出该src_neighbor在该community的邻边集合
        set<uint32_t>tar_neighbor_set = record[threshold_id][*set_it][com_id];
        int result_size = src_neighbor.size() > tar_neighbor_set.size() ? src_neighbor.size() : tar_neighbor_set.size();
        vector<uint32_t>insection_set(result_size); //结果集要用vector存放,并且初始化大小
        set_intersection(src_neighbor.begin(), src_neighbor.end(), tar_neighbor_set.begin(), tar_neighbor_set.end(), insection_set.begin());
        sim += (insection_set.size()) / (vertexDegree[sourceNode] + vertexDegree[*set_it]);
    }
    return sim;
}
*/