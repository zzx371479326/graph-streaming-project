#include"record.h"

/*
* edgeList: �洢��ͼ
* nodeCommunityList: ÿ�����������Ķ�Ԫ��<community,volume>�������ж������ʼ��Ϊÿ�������community�ǿ�
* communityEdgeList: ÿ��community�������ı�
* maxVolumeList: �������
* communityVolumeList: ÿ��community������
*/
int EdgeModAlgo(const std::vector< Edge >& edgeList,
    Record& record,
    uint32_t condition,
    uint32_t randSeed = 0
) {
    /*
    // Random shuffle
    std::srand(randSeed); // use current time as seed for random generator
    std::vector< Edge > shuffledEdgeList(edgeList);
    std::random_shuffle(shuffledEdgeList.begin(), shuffledEdgeList.end(), myrandom);
    */
    // Aggregation
    std::vector< uint32_t > nextCommunityIdList(record.maxVolumeList.size(), 1);
    // for ѭ������vector
    for (Edge edge : edgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        // cout << sourceNode << ',' << targetNode << ' '; //������
        // ��ÿһ��threshold�����б���
        for (uint32_t i = 0; i < record.maxVolumeList.size(); i++) {
            // ��ǰ����ֵ
            uint32_t maxVolume = record.maxVolumeList[i];
            // Case1:
            // �ñߵ����������δ����community
            bool src_available = record.IsAvailable(i, sourceNode);
            bool tar_abailable = record.IsAvailable(i, targetNode);
            if (!src_available && !tar_abailable)
            {
                // �½�һ��commmunity
                record.edgeCreateNewCommunity(i, edge);
            }
            // Case2:
            // ���ֻ��һ�˶��������community����ñ߷��������community
            else if (!(src_available && tar_abailable))
            {
                Node newNode = !record.IsAvailable(i, sourceNode) ? sourceNode : targetNode;
                Node oldNode = record.IsAvailable(i, sourceNode) ? sourceNode : targetNode;
                // �õ���ӦoldNode�ļ�¼
                NodeCommunityRecord oldNodeTopRecord = record.getTopRecord(i, oldNode);
                // �߼��ϸ���
                record.edgeInsert(i, oldNodeTopRecord.com_id, edge);
                // ��¼����
                record.insertNodeComRecord(i, newNode, oldNodeTopRecord.com_id, set<uint32_t>{newNode}); // �����½ڵ��¼
                record.updateNodeComRecord_AddNode(i, oldNode, 0, newNode); // ����oldNode��˵ļ�¼
            }
            // Case3:
            // ������˶�����п��õ�community
            else
            {
                // cout << "node merge..." << endl;//����
                // ����С�Ľڵ�ϲ���������Ľڵ�
                // ------------------
                // �����ҳ�����Ҫ�ϲ���community�����
                uint32_t src_com = record.getTopRecord(i, sourceNode).com_id;
                uint32_t tar_com = record.getTopRecord(i, targetNode).com_id;
                // �����ͬһ��communty�����úϲ�
                if (src_com == tar_com) {
                    // �߼��ϸ���
                    record.edgeInsert(i, src_com, edge);
                    // ��¼����
                    record.updateNodeComRecord_AddNode(i, sourceNode, 0, targetNode);
                    record.updateNodeComRecord_AddNode(i, targetNode, 0, sourceNode);
                    continue;
                }
                // �����ҳ�community volume��С�ߺͽϴ���
                uint32_t min_id = record.communityEdgeList[i][src_com].size() < record.communityEdgeList[i][tar_com].size() ? sourceNode : targetNode;
                uint32_t max_id = min_id == sourceNode ? targetNode : sourceNode;
                // ����ת��
                // ����ط���bug?
                if (sourceNode == 763 && targetNode == 465909) {
                    int a = 0;
                }
                record.edgeTrans(i, min_id, max_id);
                // ���ñ������max_com��
                record.edgeInsert(i, record.nodeCommunityList[i][max_id].begin()->com_id, edge);
            }
        }
    }
    return 0;
}

/*
// �㷨2ʵ�֣�Modularity��Greedy�㷨
int EdgeModGreedyAlgo(const std::vector< Edge >& edgeList,
    vector<vector<NodeNeighborCommunity>>& record,
    vector< std::vector< vector<Edge>>> communityEdgeList, // �������
    std::vector< uint32_t >& maxVolumeList,
    uint32_t condition,
    uint32_t randSeed = 0
) {
    // Random shuffle
    std::srand(randSeed); // use current time as seed for random generator
    std::vector< Edge > shuffledEdgeList(edgeList);
    std::random_shuffle(shuffledEdgeList.begin(), shuffledEdgeList.end(), myrandom);
    // Aggregation
    std::vector< uint32_t > nextCommunityIdList(maxVolumeList.size(), 1);
    // ��¼��
    ; // ÿ�������community-�ڱ���Ϣ

    // for ѭ������vector
    for (Edge edge : shuffledEdgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        // ��ÿһ��threshold�����б���
        for (uint32_t i = 0; i < maxVolumeList.size(); i++) {
            // ��ǰ����ֵ
            uint32_t maxVolume = maxVolumeList[i];
            // ���ü�¼ֵ
            double argmax = -1;
            // Case1:
            // �ñߵ����������δ����community
            if (record[i][sourceNode].degree == 0 && record[i][sourceNode].degree == 0)
            {
                // �½�һ��commmunity
                communityEdgeList[i].push_back(vector<Edge>{edge});
                int com_id = communityEdgeList.size() - 1;// ��Ŵ�0��ʼ
                record[i][sourceNode].update(targetNode, com_id);
                record[i][targetNode].update(sourceNode, com_id);
            }
            // Case��©: community size�������
            // Case2:
            // �ж����Ѿ�������community
            else {
                // ���Ҷ�����ڱ߼���
                map<int, double>com_score;
                // ����source���ھ�
                updateCommunityScore(i, record, com_score, sourceNode);
                // ����target���ھ�
                updateCommunityScore(i, record, com_score, targetNode);
                // ��map����ѡ����������community_id
                int max_id = -1;
                double max_score = -1;
                for (map<int, double>::iterator it = com_score.begin(); it != com_score.end(); it++) {
                    // �����community������������max���򲻼�������
                    if (communityEdgeList[i][it->first].size() >= maxVolume)
                        continue;
                    if (it->second > max_score) {
                        max_id = it->first;
                        max_score = it->second;
                    }
                }
                // �ҵ������community֮�󣬽��ñ߲��뵽��community����
                communityEdgeList[i][max_id].push_back(edge);
                record[i][sourceNode].update(targetNode, max_id);
                record[i][targetNode].update(sourceNode, max_id);
            }
        }
    }
    return 0;
}



// �㷨3ʵ�֣���������ṹGreedy�㷨
int EdgeNetGreedyAlgo(const std::vector< Edge >& edgeList,
    vector<vector<map<uint32_t, set<uint32_t>>>>& record, //��¼����
    vector< std::vector< vector<Edge>>> communityEdgeList, // �������
    vector<uint32_t>vertexDegree, // �����ܶ���
    std::vector< uint32_t >& maxVolumeList,
    uint32_t condition,
    uint32_t randSeed = 0
) {
    // Random shuffle
    std::srand(randSeed); // use current time as seed for random generator
    std::vector< Edge > shuffledEdgeList(edgeList);
    std::random_shuffle(shuffledEdgeList.begin(), shuffledEdgeList.end(), myrandom);
    // Aggregation
    std::vector< uint32_t > nextCommunityIdList(maxVolumeList.size(), 1);
    // ��¼��
    ; // ÿ�������community-�ڱ���Ϣ

    // for ѭ������vector
    for (Edge edge : shuffledEdgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        // ��ÿһ��threshold�����б���
        for (uint32_t i = 0; i < maxVolumeList.size(); i++) {
            // ��ǰ����ֵ
            uint32_t maxVolume = maxVolumeList[i];
            // Case1:
            // ����ñߵ����������¼��Ϊ�գ����½�community
            if (record[i][sourceNode].size() == 0 && record[i][targetNode].size() == 0) {
                // �½�һ��commmunity
                communityEdgeList[i].push_back(vector<Edge>{edge});
                int com_id = communityEdgeList.size() - 1;// ��Ŵ�0��ʼ
                record[i][sourceNode].insert(pair<uint32_t, set<uint32_t>>(com_id, set<uint32_t>{targetNode}));
                record[i][targetNode].insert(pair<uint32_t, set<uint32_t>>(com_id, set<uint32_t>{sourceNode}));
            }
            // �������һ������
            else if (1) {

            }
            else {
                // ǰ��: ����i�Ͷ���j��Ҫ��ͬһ��community���м�¼
                // ����м�¼����������Ӧ��¼��community��Ѱ�����ƶ�֮������community
                // ��������i,j��community setȻ���ÿһ����ͬ��community���б���
                double sim_sum = 0;
                int max_id = -1;
                for (map< uint32_t, set<uint32_t>>::iterator src_it = record[i][sourceNode].begin(); src_it != record[i][sourceNode].end(); src_it++) {
                    // ��������Ѿ����ˣ�����Ҫ������community
                    if(communityEdgeList[i][src_it->first].size() >= maxVolume)
                        break;
                    map< uint32_t, set<uint32_t>>::iterator tar_it = record[i][targetNode].find(src_it->first);
                    if (tar_it == record[i][targetNode].end())
                        break;
                    // ���ƥ��ɹ�
                    double sim = 0;
                    // �ȼ����붥��i�����ıߵ����ƶ�
                    sim += cal_sim(i, record, sourceNode, targetNode, vertexDegree, src_it, tar_it);
                    if (sim > sim_sum)
                    {
                        max_id = src_it->first;
                        sim_sum = sim;
                    }
                }
                // ���ƥ��ʧ�ܣ���i,jû�й�ͬ�Ķ���
            }
        }
        vertexDegree[sourceNode]++;
        vertexDegree[targetNode]++;
    }
    return 0;
}

void AlgoPrint(const char* fileName,
    vector< vector<Edge>>& communityEdgeList,
    bool removeSingleton) // singleton:�Ƿ�ȥ�����������community
{
    std::ofstream outFile;
    outFile.open(fileName);
    uint32_t nbCommunities = 0;
    // ����community�е�ÿһ��Ԫ��cv
    for (auto kv : communityEdgeList) {
        if (!removeSingleton || kv.size() > 1) {
            std::vector<Edge>::iterator it2 = kv.second.begin();
            Node nodeId;
            while (true) {
                nodeId = *it2;
                ++it2;
                if (it2 != kv.second.end()) {
                    outFile << nodeId << " ";
                }
                else {
                    break;
                }
            }
            outFile << nodeId << std::endl;
            nbCommunities++;
        }
    }
    printf("%-32s %i\n", "Number of Communities:", nbCommunities);
    outFile.close();
}

*/