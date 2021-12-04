#include"record.h"

Record::Record(std::vector< std::vector<std::multiset<NodeCommunityRecord>>>& nodeCommunityList,
    std::vector< std::vector< vector<Edge> > >& communityEdgeList,
    std::vector< uint32_t >& maxVolumeList){
    this->nodeCommunityList = nodeCommunityList;
    this->communityEdgeList = communityEdgeList;
    this->maxVolumeList = maxVolumeList;
}

bool Record::IsAvailable(uint32_t threshold_id, uint32_t node_id) {
    // ����ýڵ��δ�������
    if (nodeCommunityList[threshold_id][node_id].size() == 0)
        return false;
    // ����ýڵ�����ȫ��community��������
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
    //    // ���ɨ���ʱ��÷������ˣ��򽫸÷�������Ϣ���б���ɾ��
    //    if (communityEdgeList[threshold_id][it->com_id].size() >= maxVolumeList[threshold_id] || communityEdgeList[threshold_id][it->com_id].size() <= 0) {
    //        it = nodeCommunityList[threshold_id][node_id].erase(it);
    //        if (it == nodeCommunityList[threshold_id][node_id].end())
    //            break;
    //        continue;
    //    }
    //}
    // ���ɨ���껹��community�ļ�¼����˵����Ȼ�п��õķ���
    if (nodeCommunityList[threshold_id][node_id].size() == 0)
        return false;
    return true;
}

void Record::edgeInsert(uint32_t threshold_id, uint32_t com_id, Edge edge) {
    communityEdgeList[threshold_id][com_id].push_back(edge);
}

void Record::edgeCreateNewCommunity(uint32_t threshold_id, Edge edge) {
    // �½�һ��commmunity
    communityEdgeList[threshold_id].push_back(vector<Edge>{edge});
    uint32_t com_id = communityEdgeList[threshold_id].size() - 1; // ��Ŵ�0��ʼ
    // ��¼��Ϣ
    insertNodeComRecord(threshold_id, edge.first, com_id, set<uint32_t>{edge.second});
    insertNodeComRecord(threshold_id, edge.second, com_id, set<uint32_t>{edge.first});
}

void Record::insertNodeComRecord(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor) {
    nodeCommunityList[threshold_id][node_id].insert(NodeCommunityRecord(node_id, com_id, neighbor));
}

void Record::updateNodeComRecord_AddNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode) {
    // ��node_id�ļ�¼�����otherNode
    // ������ɾ������insert������ֱ���޸�
    multiset<NodeCommunityRecord>::iterator it;
    it = findRecord(threshold_id, node_id, com_id);
    // �ж��Ƿ���com_id�ļ�¼
    if (it != nodeCommunityList[threshold_id][node_id].end()) {
        NodeCommunityRecord tmp = *it;
        tmp.neighbor.insert(otherNode);
        nodeCommunityList[threshold_id][node_id].erase(it);
        nodeCommunityList[threshold_id][node_id].insert(tmp);
    }
    // û�еĻ�ֱ�Ӵ����µļ���
    else {
        insertNodeComRecord(threshold_id, node_id, com_id, set<uint32_t>{otherNode});
    }
}

void Record::updateNodeComRecord_AddMulNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, set<uint32_t> neighbor) {
    // ��node_id�ļ�¼�����otherNode
    // ������ɾ������insert������ֱ���޸�
    multiset<NodeCommunityRecord>::iterator it;
    it = findRecord(threshold_id, node_id, com_id);
    // �ж��Ƿ���com_id�ļ�¼
    if (it != nodeCommunityList[threshold_id][node_id].end()) {
        NodeCommunityRecord tmp = *it;
        tmp.neighbor.insert(neighbor.begin(), neighbor.end());
        nodeCommunityList[threshold_id][node_id].erase(it);
        nodeCommunityList[threshold_id][node_id].insert(tmp);
    }
    // û�еĻ�ֱ�Ӵ����µļ���
    else {
        insertNodeComRecord(threshold_id, node_id, com_id, neighbor);
    }
}

void Record::updateNodeComRecord_DeleteNode(uint32_t threshold_id, uint32_t node_id, uint32_t com_id, uint32_t otherNode) {
    multiset<NodeCommunityRecord>::iterator it = findRecord(threshold_id, node_id, com_id);
    // ������������bug
    if (it == nodeCommunityList[threshold_id][node_id].end())
        return;
    // ���neighborֻ��һ������ֱ��ɾ���ü�¼
    NodeCommunityRecord tmp = *it;
    nodeCommunityList[threshold_id][node_id].erase(it);
    if (tmp.neighbor.size() > 1){
        tmp.neighbor.erase(otherNode);
        nodeCommunityList[threshold_id][node_id].insert(tmp);
    }
}

NodeCommunityRecord Record::getTopRecord(uint32_t threshold_id, uint32_t node_id) {
    // �������������community
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
    // �ҳ�����community
    uint32_t min_com = nodeCommunityList[threshold_id][min_id].begin()->com_id;
    uint32_t max_com = nodeCommunityList[threshold_id][max_id].begin()->com_id;
    // ��communityEdgeList���б���
    vector<Edge>::iterator it;
    for (it = communityEdgeList[threshold_id][min_com].begin(); it != communityEdgeList[threshold_id][min_com].end(); it++) {
        // ����ñߺ���min_id����ִ���Ƴ��Ͳ������
        if (it->first == min_id || it->second == min_id) {
            uint32_t otherNode = it->first == min_id ? it->second : it->first;
            // ��ִ�бߵĲ���ɾ������
            edgeInsert(threshold_id, max_com, *it);
            it = communityEdgeList[threshold_id][min_com].erase(it); // ��ɾ������
            // ����ohter�ڵ�ļ�¼
            // otherNode��¼����...
            updateNodeComRecord_DeleteNode(threshold_id, otherNode, min_com, min_id); // ɾ����¼(��Bug--�� (763�� 465909)���޷�ɾ����¼����Ϊ��¼�Ѿ���ɾ���ˣ���
            updateNodeComRecord_AddNode(threshold_id, otherNode, max_com, min_id); // ��Ӽ�¼
            if (it == communityEdgeList[threshold_id][min_com].end())
                break;
            continue;
        }
    }
    // min_id��¼����...
    // �ҳ�min_node��min_com��¼ָ��
    multiset<NodeCommunityRecord>::iterator min_com_it = findRecord(threshold_id, min_id, min_com);
    // ��max_com�ļ�¼�����
    updateNodeComRecord_AddMulNode(threshold_id, min_id, max_com, min_com_it->neighbor);
    // ɾ��min_id��min_com��¼
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
// �㷨2����ʵ��
// �����ڵ��ھӵ�ͬʱ��Ϊÿ��community�������
void updateCommunityScore(uint32_t threshold_id, vector<vector<NodeNeighborCommunity>>& record, map<int, double>& com_score, uint32_t node_id) {
    // ����ĳ���ڵ�������ھ�
    for (set<uint32_t>::iterator it_node = record[threshold_id][node_id].neighbor.begin(); it_node != record[threshold_id][node_id].neighbor.end(); it_node++) {
        // �����ýڵ������community
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

// �㷨3����ʵ��
// �㷨3���ƶȵļ���
double cal_sim(uint32_t threshold_id, vector<vector<map<uint32_t, set<uint32_t>>>>& record, uint32_t sourceNode, uint32_t targetNode, vector<uint32_t>& vertexDegree, map< uint32_t, set<uint32_t>>::iterator src_it, map< uint32_t, set<uint32_t>>::iterator tar_it) {
    uint32_t com_id = src_it->first;
    double sim = 0;
    // ��ȡ��target�ڸ�com���ھ�
    set<uint32_t>tar_neighbor = tar_it->second;
    // ��ȡ��source�ڸ�com���ھ�
    set<uint32_t>src_neighbor = src_it->second;
    // �ȼ���src�������ھ�
    set<uint32_t>::iterator set_it;
    for (set_it = src_neighbor.begin(); set_it != src_neighbor.end(); set_it++) {
        // �ҳ���src_neighbor�ڸ�community���ڱ߼���
        set<uint32_t>src_neighbor_set = record[threshold_id][*set_it][com_id];
        int result_size = tar_neighbor.size() > src_neighbor_set.size() ? tar_neighbor.size() : src_neighbor_set.size();
        vector<uint32_t>insection_set(result_size); //�����Ҫ��vector���,���ҳ�ʼ����С
        set_intersection(tar_neighbor.begin(), tar_neighbor.end(), src_neighbor_set.begin(), src_neighbor_set.end(), insection_set.begin());
        sim += (insection_set.size()) / (vertexDegree[targetNode] + vertexDegree[*set_it]);
    }
    // �ټ���tar�����ھ�
    for (set_it = tar_neighbor.begin(); set_it != tar_neighbor.end(); set_it++) {
        // �ҳ���src_neighbor�ڸ�community���ڱ߼���
        set<uint32_t>tar_neighbor_set = record[threshold_id][*set_it][com_id];
        int result_size = src_neighbor.size() > tar_neighbor_set.size() ? src_neighbor.size() : tar_neighbor_set.size();
        vector<uint32_t>insection_set(result_size); //�����Ҫ��vector���,���ҳ�ʼ����С
        set_intersection(src_neighbor.begin(), src_neighbor.end(), tar_neighbor_set.begin(), tar_neighbor_set.end(), insection_set.begin());
        sim += (insection_set.size()) / (vertexDegree[sourceNode] + vertexDegree[*set_it]);
    }
    return sim;
}
*/