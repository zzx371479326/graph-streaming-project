#include"record.h"

/*
* edgeList: 存储的图
* nodeCommunityList: 每个顶点所属的二元组<community,volume>，可能有多个，初始化为每个顶点的community是空
* communityEdgeList: 每个community所包含的边
* maxVolumeList: 最大容量
* communityVolumeList: 每个community的容量
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
    // for 循环遍历vector
    for (Edge edge : edgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        // cout << sourceNode << ',' << targetNode << ' '; //调试用
        // 对每一个threshold都进行遍历
        for (uint32_t i = 0; i < record.maxVolumeList.size(); i++) {
            // 当前的阈值
            uint32_t maxVolume = record.maxVolumeList[i];
            // Case1:
            // 该边的两个顶点均未分配community
            bool src_available = record.IsAvailable(i, sourceNode);
            bool tar_abailable = record.IsAvailable(i, targetNode);
            if (!src_available && !tar_abailable)
            {
                // 新建一个commmunity
                record.edgeCreateNewCommunity(i, edge);
            }
            // Case2:
            // 如果只有一端顶点分配了community，则该边分配给最大的community
            else if (!(src_available && tar_abailable))
            {
                Node newNode = !record.IsAvailable(i, sourceNode) ? sourceNode : targetNode;
                Node oldNode = record.IsAvailable(i, sourceNode) ? sourceNode : targetNode;
                // 得到相应oldNode的记录
                NodeCommunityRecord oldNodeTopRecord = record.getTopRecord(i, oldNode);
                // 边集合更新
                record.edgeInsert(i, oldNodeTopRecord.com_id, edge);
                // 记录更新
                record.insertNodeComRecord(i, newNode, oldNodeTopRecord.com_id, set<uint32_t>{newNode}); // 更新新节点记录
                record.updateNodeComRecord_AddNode(i, oldNode, 0, newNode); // 更新oldNode最顶端的记录
            }
            // Case3:
            // 如果两端顶点均有可用的community
            else
            {
                // cout << "node merge..." << endl;//测试
                // 容量小的节点合并到容量大的节点
                // ------------------
                // 首先找出各自要合并的community的序号
                uint32_t src_com = record.getTopRecord(i, sourceNode).com_id;
                uint32_t tar_com = record.getTopRecord(i, targetNode).com_id;
                // 如果是同一个communty，则不用合并
                if (src_com == tar_com) {
                    // 边集合更新
                    record.edgeInsert(i, src_com, edge);
                    // 记录更新
                    record.updateNodeComRecord_AddNode(i, sourceNode, 0, targetNode);
                    record.updateNodeComRecord_AddNode(i, targetNode, 0, sourceNode);
                    continue;
                }
                // 否则找出community volume较小者和较大者
                uint32_t min_id = record.communityEdgeList[i][src_com].size() < record.communityEdgeList[i][tar_com].size() ? sourceNode : targetNode;
                uint32_t max_id = min_id == sourceNode ? targetNode : sourceNode;
                // 进行转移
                // 这个地方有bug?
                if (sourceNode == 763 && targetNode == 465909) {
                    int a = 0;
                }
                record.edgeTrans(i, min_id, max_id);
                // 将该边添加至max_com中
                record.edgeInsert(i, record.nodeCommunityList[i][max_id].begin()->com_id, edge);
            }
        }
    }
    return 0;
}

/*
// 算法2实现：Modularity的Greedy算法
int EdgeModGreedyAlgo(const std::vector< Edge >& edgeList,
    vector<vector<NodeNeighborCommunity>>& record,
    vector< std::vector< vector<Edge>>> communityEdgeList, // 结果集合
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
    // 记录表
    ; // 每个顶点的community-邻边信息

    // for 循环遍历vector
    for (Edge edge : shuffledEdgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        // 对每一个threshold都进行遍历
        for (uint32_t i = 0; i < maxVolumeList.size(); i++) {
            // 当前的阈值
            uint32_t maxVolume = maxVolumeList[i];
            // 设置记录值
            double argmax = -1;
            // Case1:
            // 该边的两个顶点均未分配community
            if (record[i][sourceNode].degree == 0 && record[i][sourceNode].degree == 0)
            {
                // 新建一个commmunity
                communityEdgeList[i].push_back(vector<Edge>{edge});
                int com_id = communityEdgeList.size() - 1;// 标号从0开始
                record[i][sourceNode].update(targetNode, com_id);
                record[i][targetNode].update(sourceNode, com_id);
            }
            // Case遗漏: community size满的情况
            // Case2:
            // 有顶点已经分配了community
            else {
                // 查找顶点的邻边集合
                map<int, double>com_score;
                // 遍历source的邻居
                updateCommunityScore(i, record, com_score, sourceNode);
                // 遍历target的邻居
                updateCommunityScore(i, record, com_score, targetNode);
                // 在map里面选出分数最大的community_id
                int max_id = -1;
                double max_score = -1;
                for (map<int, double>::iterator it = com_score.begin(); it != com_score.end(); it++) {
                    // 如果该community的容量超过了max，则不计入其中
                    if (communityEdgeList[i][it->first].size() >= maxVolume)
                        continue;
                    if (it->second > max_score) {
                        max_id = it->first;
                        max_score = it->second;
                    }
                }
                // 找到了最佳community之后，将该边插入到该community里面
                communityEdgeList[i][max_id].push_back(edge);
                record[i][sourceNode].update(targetNode, max_id);
                record[i][targetNode].update(sourceNode, max_id);
            }
        }
    }
    return 0;
}



// 算法3实现：基于网络结构Greedy算法
int EdgeNetGreedyAlgo(const std::vector< Edge >& edgeList,
    vector<vector<map<uint32_t, set<uint32_t>>>>& record, //记录集合
    vector< std::vector< vector<Edge>>> communityEdgeList, // 结果集合
    vector<uint32_t>vertexDegree, // 顶点总度数
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
    // 记录表
    ; // 每个顶点的community-邻边信息

    // for 循环遍历vector
    for (Edge edge : shuffledEdgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        // 对每一个threshold都进行遍历
        for (uint32_t i = 0; i < maxVolumeList.size(); i++) {
            // 当前的阈值
            uint32_t maxVolume = maxVolumeList[i];
            // Case1:
            // 如果该边的两个顶点记录均为空，则新建community
            if (record[i][sourceNode].size() == 0 && record[i][targetNode].size() == 0) {
                // 新建一个commmunity
                communityEdgeList[i].push_back(vector<Edge>{edge});
                int com_id = communityEdgeList.size() - 1;// 标号从0开始
                record[i][sourceNode].insert(pair<uint32_t, set<uint32_t>>(com_id, set<uint32_t>{targetNode}));
                record[i][targetNode].insert(pair<uint32_t, set<uint32_t>>(com_id, set<uint32_t>{sourceNode}));
            }
            // 如果其中一个顶点
            else if (1) {

            }
            else {
                // 前提: 顶点i和顶点j都要在同一个community中有记录
                // 如果有记录，则在有相应记录的community上寻找相似度之和最大的community
                // 创建顶点i,j的community set然后对每一个共同的community进行遍历
                double sim_sum = 0;
                int max_id = -1;
                for (map< uint32_t, set<uint32_t>>::iterator src_it = record[i][sourceNode].begin(); src_it != record[i][sourceNode].end(); src_it++) {
                    // 如果容量已经满了，则需要跳过该community
                    if(communityEdgeList[i][src_it->first].size() >= maxVolume)
                        break;
                    map< uint32_t, set<uint32_t>>::iterator tar_it = record[i][targetNode].find(src_it->first);
                    if (tar_it == record[i][targetNode].end())
                        break;
                    // 如果匹配成功
                    double sim = 0;
                    // 先计算与顶点i相连的边的相似度
                    sim += cal_sim(i, record, sourceNode, targetNode, vertexDegree, src_it, tar_it);
                    if (sim > sim_sum)
                    {
                        max_id = src_it->first;
                        sim_sum = sim;
                    }
                }
                // 如果匹配失败，则i,j没有共同的顶点
            }
        }
        vertexDegree[sourceNode]++;
        vertexDegree[targetNode]++;
    }
    return 0;
}

void AlgoPrint(const char* fileName,
    vector< vector<Edge>>& communityEdgeList,
    bool removeSingleton) // singleton:是否去除单个顶点的community
{
    std::ofstream outFile;
    outFile.open(fileName);
    uint32_t nbCommunities = 0;
    // 对于community中的每一个元素cv
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