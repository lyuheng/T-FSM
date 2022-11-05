#pragma once

#include <atomic>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <stack>
#include <list>
#include <fstream>

#include "constack.h"
#include "conque.h"
#include "conmap.h"
#include "rwlock.h"
#include "taskprogmap.h"

#include "grami.h"

using namespace std;

#define WAIT_TIME_WHEN_IDLE 100000

// enable printing all found frequent patterns
#define VERBOSE_MODE

// enable optimized version of subgraph matching methods
#define OPTIMIZED_MATCH

// enable timeout search
#define TIMEOUT_THRESHOLD 0.1

// enable time threshold for decomposing tasks for timeout task
#define DECOMPOSE_TIME_THRESHOLD 1

// Set max vertex sizes of mined patterns
#define MAX_PATTERN_SIZE 32



#define SEQNO_LIMIT 9223372036854775807

GraMi grami;

int activeQ_num(0);
rwlock activeQ_lock;
void *global_activeQ_list;
int activeQ_list_capacity = 40;


void *global_data_stack;

size_t MINI_BATCH_NUM = 800;

int RT_THRESHOLD_FOR_REFILL = 800;

atomic<int> qid(0);

size_t num_compers = 32;

condition_variable cv_go;
mutex mtx_go;
// Protected by mtx_go above
bool ready_go = true;

atomic<int> global_num_idle(0);
atomic<bool> global_end_label(false);

// results counter
vector<ui> results_counter;

TaskProgMap global_prog_map;

ofstream * fout;
mutex fmtx;

// use for push down pruning
unordered_map<subPattern, VtxSetVec, subPatternHashCode> cache;
rwlock cache_mtx;

