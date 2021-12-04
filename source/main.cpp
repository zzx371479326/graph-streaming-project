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
#include "record.h"
using namespace std;


#define CHECK_ARGUMENT_STRING(index, option,variable,setVariable) \
    if( strcmp(argv[index],option) == 0 ){ \
            setVariable = true; \
            if( (index+1) < argc ) { \
                variable = argv[index+1]; \
            } else { \
                printf( "Invalid options.\n" ); \
                return 1;\
            }\
        }

#define CHECK_ARGUMENT_FLOAT(index, option,variable,setVariable) \
    if( strcmp(argv[index],option) == 0 ){ \
            setVariable = true; \
            if( (index+1) < argc ) { \
                variable = atof(argv[index+1]); \
            } else { \
                printf( "Invalid options.\n" ); \
                return 1;\
            }\
        }

#define CHECK_ARGUMENT_INT(index, option,variable,setVariable) \
    if( strcmp(argv[index],option) == 0 ){ \
            setVariable = true; \
            if( (index+1) < argc ) { \
                variable = atoi(argv[index+1]); \
            } else { \
                printf( "Invalid options.\n" ); \
                return 1;\
            }\
        }

#define CHECK_FLAG(index, option,setVariable) \
    if( strcmp(argv[index],option) == 0 ){ \
            setVariable = true; \
        }

int StreamComAlgo(const std::vector< Edge >& edgeList,
                 std::vector< std::vector< uint32_t > >& nodeCommunityList,
                 std::vector< uint32_t >& nodeDegree,
                 std::vector< std::vector< uint32_t > >& communityVolumeList,
                 std::vector< uint32_t >& maxVolumeList,
                 uint32_t condition,
                 uint32_t randSeed=0
                 ) {
    // Random shuffle
    std::srand(randSeed); // use current time as seed for random generator
    std::vector< Edge > shuffledEdgeList (edgeList);
    std::random_shuffle(shuffledEdgeList.begin(), shuffledEdgeList.end(), myrandom);
    // Aggregation
    std::vector< uint32_t > nextCommunityIdList (maxVolumeList.size(), 1);
    // for 循环遍历vector
    for (Edge edge : shuffledEdgeList) {
        Node sourceNode = edge.first;
        Node targetNode = edge.second;
        nodeDegree[sourceNode]++;
        nodeDegree[targetNode]++;
        // 对每一个threshold都进行遍历
        for (uint32_t i = 0; i < maxVolumeList.size(); i++) {
            // 当前的阈值
            uint32_t maxVolume = maxVolumeList[i];
            // nextCommunityIdList中存储了当前一共分配了多少个community
            if(nodeCommunityList[i][sourceNode] == 0) {
                nodeCommunityList[i][sourceNode] = nextCommunityIdList[i];
                nextCommunityIdList[i]++;
            }
            if(nodeCommunityList[i][targetNode] == 0) {
                nodeCommunityList[i][targetNode] = nextCommunityIdList[i];
                nextCommunityIdList[i]++;
            }
            uint32_t sourceCommunityVolume = ++communityVolumeList[i][nodeCommunityList[i][sourceNode]];
            uint32_t targetCommunityVolume = ++communityVolumeList[i][nodeCommunityList[i][targetNode]];
            if(((condition == 0) && ((sourceCommunityVolume <= maxVolume) && (targetCommunityVolume <= maxVolume)))
            || ((condition == 1) && ((sourceCommunityVolume <= maxVolume) || (targetCommunityVolume <= maxVolume)))) {
                if(sourceCommunityVolume < targetCommunityVolume) {
                    communityVolumeList[i][nodeCommunityList[i][sourceNode]] -= nodeDegree[sourceNode];
                    communityVolumeList[i][nodeCommunityList[i][targetNode]] += nodeDegree[sourceNode];
                    nodeCommunityList[i][sourceNode] = nodeCommunityList[i][targetNode];
                } else {
                    communityVolumeList[i][nodeCommunityList[i][targetNode]] -= nodeDegree[targetNode];
                    communityVolumeList[i][nodeCommunityList[i][sourceNode]] += nodeDegree[targetNode];
                    nodeCommunityList[i][targetNode] = nodeCommunityList[i][sourceNode];
                }
            }
        }
    }
    return 0;
}

/*
用法示例：
-f com-amazon.ungraph.txt --skip 4 -o output.txt
*/

static void PrintUsage() {
    printf("Usage: streamcom <flags>\n");
    printf("Availaible flags:\n");
    printf("\t-algo [algo name] : vertexClu or edgeClu.\n");
    printf("\t-f [graph file name] : Specifies the graph file.\n");
    printf("\t--skip [number of lines] : Specifies the number of lines to skip in the graph file.\n");
    printf("\t-o [output file name] : Specifies the output file for communities.\n");
    printf("\t--vmax-start [maximum volume range start] : Specifies the maximum volume for the aggregation phase, beginning of the range.\n");
    printf("\t--vmax-end [maximum volume range end] : Specifies the maximum volume for the aggregation phase, end of the range.\n");
    printf("\t-c [condition] : Specifies the aggregation condition for the aggregation phase: 0 is AND and 1 is OR.\n");
    printf("\t--seed [random seed]: Specifies the random seed if the edges need to be shuffle.\n");
    printf("\t--niter [number of iteration]: Specifies the number of iteration of the algorithm.\n");
}

int EdgeModAlgo(const std::vector< Edge >& edgeList,
    Record& record,
    uint32_t condition,
    uint32_t randSeed = 0
);

void vertexcluster(
    uint32_t nIter,
    uint32_t volumeThresholdStart,
    uint32_t volumeThresholdEnd,
    uint32_t condition,
    int randomSeed,
    Node maxNodeId,
    long unsigned totalTime,
    long unsigned initTime,
    long unsigned spentTime,
    long unsigned loadingTime,
    long unsigned algorithmTime,
    vector< Edge >& edgeList,
    char* outputFileName
) 
{
    // 以下部分全部注释
    //======================================================================

    for (uint32_t iter = 0; iter < nIter; iter ++) {
        // 不同的volume对应不同的结果，因此每个volume都需要跑一次
        std::vector< uint32_t > volumeThresholdList;
        for (uint32_t volumeThreshold = volumeThresholdStart; volumeThreshold <= volumeThresholdEnd; volumeThreshold++) {
            volumeThresholdList.push_back(volumeThreshold);
        }
        // 序号从零开始
        std::vector< uint32_t > nodeDegree (maxNodeId + 1);
        // 每个顶点对应的community序号；volumeThresholdList.size()的意思是可能有多组实验数据
        std::vector< std::vector< uint32_t > > nodeCommunityList (volumeThresholdList.size(), std::vector< uint32_t > (maxNodeId + 1, 0));

        std::vector< std::vector< uint32_t > > communityVolumeList (volumeThresholdList.size(), std::vector< uint32_t > (maxNodeId + 1, 0));

        //=================== ALGORITHM  =======================================
        printf("Start algorithm...\n");
        initTime = StartClock();
        StreamComAlgo(edgeList, nodeCommunityList, nodeDegree, communityVolumeList, volumeThresholdList, condition, randomSeed * (1 + iter));
        spentTime = StopClock(initTime);
        algorithmTime += spentTime;
        totalTime += spentTime;
        printf("%-32s %lu ms\n", "Algorithm time:", spentTime);
        //======================================================================

        //======================== PRINT RESULTS ===============================
        initTime = StartClock();
        for (uint32_t i = 0; i < volumeThresholdList.size(); i++) {
            std::map< uint32_t, std::set< Node > > communities;
            GetCommunities(nodeCommunityList[i], maxNodeId, communities);
            std::string outputFileName_i = outputFileName;
            outputFileName_i += "_" + std::to_string(iter) + "_" + std::to_string(volumeThresholdList[i]);
            printf("%-32s %i\t", "Volume theshold: ", volumeThresholdList[i]);
            PrintPartition(outputFileName_i.c_str(), communities, false);
        }
        spentTime = StopClock(initTime);
        totalTime += spentTime;
        printf("%-32s %lu ms\n", "Print partition time:", spentTime);
        //======================================================================
    }

    printf("\n");
    printf("*******************************************************\n");
    printf("%-32s %-10lu ms\n", "Loading time:", loadingTime);
    double avgAlgorithmTime = ((float) algorithmTime) / ((float) nIter);
    printf("%-32s %0.2f ms\n", "Average algorithm time:", avgAlgorithmTime);
    printf("%-32s %-10lu ms\n", "Total execution time:", totalTime);
    printf("*******************************************************\n");
    return;
    
}
void edgecluster(
    uint32_t nIter,
    uint32_t volumeThresholdStart,
    uint32_t volumeThresholdEnd,
    uint32_t condition,
    int randomSeed,
    Node maxNodeId,
    long unsigned totalTime,
    long unsigned initTime,
    long unsigned spentTime,
    long unsigned loadingTime,
    long unsigned algorithmTime,
    vector< Edge >& edgeList,
    char* outputFileName
    )
{
    // 一下是调试部分
    for (uint32_t iter = 0; iter < nIter; iter++) {
        // 不同的volume对应不同的结果，因此每个volume都需要跑一次
        std::vector< uint32_t > volumeThresholdList;
        for (uint32_t volumeThreshold = volumeThresholdStart; volumeThreshold <= volumeThresholdEnd; volumeThreshold++) {
            volumeThresholdList.push_back(volumeThreshold);
        }
        // 以下部分数据初始化------------------------
        // 序号从零开始
        std::vector< uint32_t > nodeDegree(maxNodeId + 1);
        // 每个顶点对应的community序号；volumeThresholdList.size()的意思是可能有多组实验数据
        // 多重vector的初始化: 最后一维只需要有数字即可
        vector< std::vector<std::multiset<NodeCommunityRecord>>> nodeCommunityList(volumeThresholdList.size(), vector<multiset<NodeCommunityRecord>>(maxNodeId + 1));
        vector< std::vector< vector<Edge> > > communityEdgeList(volumeThresholdList.size());
        std::vector< std::vector< uint32_t > > communityVolumeList(volumeThresholdList.size());
        Record record(nodeCommunityList, communityEdgeList, volumeThresholdList);
        //=================== ALGORITHM  =======================================
        printf("Start algorithm...\n");
        initTime = StartClock();
        EdgeModAlgo(edgeList, record, condition, randomSeed * (1 + iter));
        spentTime = StopClock(initTime);
        algorithmTime += spentTime;
        totalTime += spentTime;
        printf("%-32s %lu ms\n", "Algorithm time:", spentTime);
        ////======================================================================
        //// 
        ////======================== PRINT RESULTS ===============================
        initTime = StartClock();
        for (uint32_t i = 0; i < volumeThresholdList.size(); i++) {
            std::string outputFileName_i = outputFileName;
            outputFileName_i += "_" + std::to_string(iter) + "_" + std::to_string(volumeThresholdList[i]) + ".txt";
            printf("%-32s %i\t\n", "Volume theshold: ", volumeThresholdList[i]);
            record.print(outputFileName_i.c_str(), i, false);
        }
        spentTime = StopClock(initTime);
        totalTime += spentTime;
        printf("%-32s %lu ms\n", "Print partition time:", spentTime);
        //======================================================================
    }
    return ;
}

int main(int argc, char ** argv) {
    bool algoSet = false;
    bool graphFileNameSet = false;
    bool outputFileNameSet = false;
    bool volumeThresholdStartSet = false;
    bool volumeThresholdEndSet = false;
    bool conditionSet = false;
    bool randomSeedSet = false;
    bool nIterSet = false;
    bool nbLinesToSkipSet = false;
    char* algoName = NULL;
    char * graphFileName = NULL;
    char * outputFileName = NULL;
    uint32_t volumeThresholdStart = 3;
    uint32_t volumeThresholdEnd = 3;
    uint32_t condition = 0;
    int randomSeed = 0;
    uint32_t nIter = 1;
    uint32_t nbLinesToSkip = 0;

    //输出当前的工作路径，方便调试
    char* buffer = _getcwd(NULL, 0);
    cout << buffer << endl;

    for(int i = 1; i < argc; i++) {
        CHECK_ARGUMENT_STRING(i, "-algo", algoName, algoSet);
        CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet);
        CHECK_ARGUMENT_INT(i, "--skip", nbLinesToSkip, nbLinesToSkipSet);
        CHECK_ARGUMENT_STRING(i, "-o", outputFileName, outputFileNameSet)
        CHECK_ARGUMENT_INT(i, "--vmax-start", volumeThresholdStart, volumeThresholdStartSet);
        CHECK_ARGUMENT_INT(i, "--vmax-end", volumeThresholdEnd, volumeThresholdEndSet);
        CHECK_ARGUMENT_INT(i, "-c", condition, conditionSet);
        CHECK_ARGUMENT_INT(i, "--seed", randomSeed, randomSeedSet);
        CHECK_ARGUMENT_INT(i, "--niter", nIter, nIterSet);
    }
    if (!algoSet) {
        printf("Algorithm name not set\n");
        PrintUsage();
        return 1;
    }

    if (!graphFileNameSet) {
        printf("Graph filename not set\n");
        PrintUsage();
        return 1;
    }

    if (!outputFileNameSet) {
        printf("Output filename not set\n");
        PrintUsage();
        return 1;
    }

    if (!randomSeedSet) {
        randomSeed = std::time(0);
    }

    if (!volumeThresholdEndSet) {
        volumeThresholdEnd = volumeThresholdStart;
    }

    printf("%-32s %i\n", "Volume threshold (Range Start):", volumeThresholdStart);
    printf("%-32s %i\n", "Volume threshold (Range End):", volumeThresholdEnd);

    printf("%-32s %i\n", "Random Seed:", randomSeed);

    long unsigned totalTime = 0,
                  initTime = 0,
                  spentTime = 0,
                  loadingTime = 0,
                  algorithmTime = 0;

    // Allocating list for edges
    std::vector< Edge > edgeList;

    //==================== LOAD THE GRAPH ==================================
    initTime = StartClock();
    printf("%-32s %s\n", "Graph file name:", graphFileName);
    printf("%-32s %i\n", "Number of lines to skip:", nbLinesToSkip);
    Node maxNodeId; //int
    LoadGraph(graphFileName, edgeList, maxNodeId, nbLinesToSkip); //加载图
    spentTime = StopClock(initTime);
    loadingTime = spentTime;
    totalTime += spentTime;
    printf("%-32s %lu\n", "Nb of edges:", edgeList.size());
    printf("%-32s %lu ms\n", "Loading time:", loadingTime);
    if (strcmp(algoName, "vertexClu") == 0)
        vertexcluster(nIter, volumeThresholdStart, volumeThresholdEnd, condition, randomSeed, maxNodeId, totalTime, initTime, spentTime, loadingTime, algorithmTime, edgeList, outputFileName);
    else if (strcmp(algoName, "edgeClu") == 0)
        edgecluster(nIter, volumeThresholdStart, volumeThresholdEnd, condition, randomSeed, maxNodeId, totalTime, initTime, spentTime, loadingTime, algorithmTime, edgeList, outputFileName);
}
