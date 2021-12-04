#ifndef DYCOLA_UTILS_H
#define DYCOLA_UTILS_H

#define WIN32

#include <stdlib.h>
#include <map>
#include <vector>
#include <set>
#include <utility>
#include <time.h>
#include <fstream>
#include "types.h"
#include <iostream>

//Ê±¼äº¯Êý
#include <time.h>
#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif
#ifdef WIN32
int gettimeofday(struct timeval* tp, void* tzp);
#endif

int myrandom (int i);
long unsigned StartClock();
long unsigned StopClock(long unsigned initTime);
int LoadGraph(char * graphFileName,
              std::vector< Edge >& edgeList,
              Node& maxNodeId,
              uint32_t nbLinesToSkip=0);
int BuildNeighborhoods(std::vector< Edge >& edgeList, std::vector< NodeSet >& nodeNeighbors);
int PrintPartition(const char* fileName,
                   std::map< uint32_t, std::set< Node > >& communities,
                   bool removeSingleton=false);
int PrintStats(const char* fileName,
               uint32_t nbCommunities,
               unsigned long algoTime);
int PrintStats(const char* fileName,
               uint32_t nbRuns,
               std::vector< uint32_t > nbCommunities,
               std::vector< unsigned long > executionTimes);
int GetCommunities(const std::vector< uint32_t > nodeCommunity,
                   Node maxNodeId,
                   std::map< uint32_t, std::set< Node > >& communities);

#endif
