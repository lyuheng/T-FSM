#pragma once

#include "global.h"
#include "grami.h"
#include "gmatch.h"
#include "decompose.h"

#include <iostream>
#include <string>
#include <deque>
#include <queue>
#include <numeric>
#include <unistd.h>

struct Task;
struct TimeOutTask;
struct TimeOutTask_Container;
struct task_container;

typedef constack<task_container *> DataStack;
typedef deque<Task *> TaskQ;
typedef list<TimeOutTask_Container *> TimeOutTaskList;
typedef list<task_container *> Qlist;


typedef conmap<string, VtxSetVec>::bucket bucket;
typedef conmap<string, VtxSetVec>::bucket_map bucket_map;

// 3 return values
#define MATCH_FOUND -1
#define MATCH_NOT_FOUND -2
#define TIMEOUT -3

using namespace std;
using namespace std::chrono;

struct Task
{
    VertexID u; // query vertex id
    VertexID v; // a vertex in domain[u]
    ui v_idx; // v's index in domain[u]

    Task(VertexID u_, VertexID v_, ui v_idx_): u(u_), v(v_), v_idx(v_idx_) {}
    Task(Task * task)
    {
        u = task->u;
        v = task->v;
        v_idx = task->v_idx;
    }
};

struct TimeOutTask
{
    Task task;
    MatchingStatus * ms;
    task_prog * prog; // creation and deletion are done externally

    bool skip;

	TaskID get_id()
    {
		return prog->tid;
	}

	int get_qid()
    {
		return prog->query_id;
	}

    TimeOutTask(Task *task_, ui size): task(task_) 
    {
        ms = new MatchingStatus(size);
        ms->vid = task.v;
        ms->v_idx = task.v_idx;

        skip = false; // set false as default
    }
    ~TimeOutTask()
    {
        // delete task;
        delete ms;
    }
};

struct TimeOutTask_Container
{
    VertexID vid;
    conque<TimeOutTask *> task_queue; // since it's conque, no lock here

    TimeOutTask_Container(int vid_): vid(vid_) {}
    
    ~TimeOutTask_Container()
    {
        TimeOutTask *totask;
        unordered_set<Task *> deleted; // prevent multiple deletion 
        while(task_queue.dequeue(totask))
        {
            delete totask;
        }
        assert(task_queue.empty()); // make sure task_queue is empty
    }
};

struct task_container
{
    int qid;
    vector<FILE *> fout_map;

    TaskQ Q_domain;
    mutex Q_domain_mtx;

    // =========== Add by Timeout List ===========
    TimeOutTaskList L_timeout;
    rwlock L_timeout_mtx;

    atomic<int> domain_front; // these two front VertexIDs indicate which is greater, determining priority of 2 queues.
                              // update after pop_front(.)

    // ============================================

    Pattern *pattern;

    // Progress status
    // *** for round-robin refill
    VertexID next_vq; // 0->1->2->3
    vector<ui> next_domainPos; // next_domainPos[v_q] = next position in domain(v_q) to refill with
    // *** for deciding the end condition of each v_q 
    vector<ui> domain_finished; // domain_finished[v_q] = how many domain elements of v_q are processed
                                // if domain_finished[v_q] == pattern->get_cands()[v_q].size(), v_q is done
    ui domain_done; // counter for candidate-finished v_q's 

    VtxSetVec domain_matches; // domain_matches[v_q] = matched domain elements, could be set by v_q'
    mutex * domain_matches_mtx; // domain_matches_mtx[v_q] protects domain_matches[v_q]


    Edges ***edge_matrix; // matching index, shared by all tasks
                          // initialized when query is added to activeQ_list

    mutex refill_mtx;  // only 1 comper performs refill at a time
                       // protects next_vq, next_domainPos

    mutex end_mtx; // protects domain_finished

    mutex domain_done_mtx; // protects domain_done

    mutex * non_cand_mtx;

    atomic<bool> frequent_tag;
    
    bool * vq_stops_refill; // vq_stops_refill[v_q] = false, if domain_matches[vq] < support
                                  // vq_stops_refill[v_q] = true, if domain_matches[vq] >= support

    rwlock * vq_stops_refill_mtx; // mutexes protect vq_stops_refill


    ui comper_counter;
    mutex comper_counter_mtx;

    atomic<bool> normal_exit;

    vector<vector<subPattern> > maps;


    vector<VertexID *> all_matching_order;
    vector<VertexID **> all_bn;
    vector<ui *> all_bn_count;


    task_container(int id)
    {
        qid = id;
    }

    void init() // called when query is added to activeQ_list
    {   
        int size = pattern->size();

        // Temporarily initialize the size
        for(ui i = 0; i < size; ++i)
        {
            TimeOutTask_Container * new_tc = new TimeOutTask_Container(i);
            L_timeout.push_back(new_tc);
        }
        
        domain_front = -1;
    
        frequent_tag = true;

        normal_exit = false;

        comper_counter = 0;

        next_vq = 0;

        domain_done = 0;

        next_domainPos.assign(size, 0);
    
        domain_finished.assign(size, 0);

        vq_stops_refill = new bool[size];
        memset(vq_stops_refill, false, sizeof(bool)*size);
        

        domain_matches.resize(size);

        domain_matches_mtx = new mutex[size];
        non_cand_mtx = new mutex[size];

        vq_stops_refill_mtx = new rwlock[size];
        
        edge_matrix = new Edges **[pattern->size()];
        for (ui i = 0; i < pattern->size(); ++i) {
            edge_matrix[i] = new Edges *[pattern->size()];
        }

        for(ui i = 0; i < pattern->size(); ++i) {
            for(ui j = 0; j < pattern->size(); ++j) {
                edge_matrix[i][j] = NULL;
            }
        }

        all_matching_order.resize(size, NULL);
        all_bn.resize(size, NULL);
        all_bn_count.resize(size, NULL);
    }

    ~task_container()
    {   
        delete[] domain_matches_mtx;
        delete[] non_cand_mtx;
        delete[] vq_stops_refill_mtx;

        delete[] vq_stops_refill;

        for (ui i = 0; i < pattern->size(); ++i) {
            for (ui j = 0; j < pattern->size(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;

        for (ui i = 0; i < pattern->size(); ++i)
        {
            if (all_bn[i] != NULL)
            {
                for (ui j = 0; j < pattern->size(); ++j)
                {
                    delete[] all_bn[i][j]; // FIXME: 
                }
            }
            delete[] all_bn[i];
            delete[] all_bn_count[i];
            delete[] all_matching_order[i];
        }

        delete pattern; // created by extend(.)


        // delete remaining tasks in Q_domain
        for(auto it = Q_domain.begin(); it != Q_domain.end(); ++it)
        {
            delete *it;
        }
        // delete remaining timeout tasks in L_timeout
        for(auto it = L_timeout.begin(); it != L_timeout.end(); ++it)
        {
            delete *it; // delete all space inside L_timeout
        }
    }

    inline bool nothing_to_refill()
    {
        return next_vq == pattern->size();
    }

    inline bool domain_finish(VertexID u)
    {
        return domain_finished[u] == pattern->get_cands()[u].size();
    }

    inline void domain_increment(VertexID u)
    {
        domain_finished[u]++;
    }

    inline bool domain_early_exit(VertexID u)
    {   
        // if # matches + # non-completed < nsupport, then early exit
        // # non-completed = # total - # finished
        domain_matches_mtx[u].lock();
        bool exit = domain_matches[u].size() + pattern->get_cands()[u].size() - domain_finished[u] < grami.nsupport_;
        domain_matches_mtx[u].unlock();

        return exit;
    }
};

class Comper
{
public:
    GMatchEngine gmatch_engine;
    GMatchEngine exist_engine;

    Decompose decomposer;

    int thread_id;

    thread main_thread;

    DataStack & data_stack = *(DataStack *)global_data_stack;

    Qlist & activeQ_list = *(Qlist *)global_activeQ_list;

    task_container* tc;

    int cur_qid;
    TaskID cur_tid;

    ui seqno;

    Comper()
    {
        seqno = 0;
    }

    ~Comper()
    {
        main_thread.join();
    }

    TaskID get_next_taskID()
	{
		TaskID id = thread_id;
		id = (id << 48); // first 16-bit is thread_id
		id += seqno;
		seqno++;
		assert(seqno < SEQNO_LIMIT);
		return id;
	}

    void add_timeout_task(TimeOutTask *totask, bool is_root_task) 
    {   
        VertexID vid = totask->task.u;

        if(is_root_task)
        {
			totask->prog = new task_prog(get_next_taskID(), -1, cur_qid);
		}
        else
        {
			totask->prog = new task_prog(get_next_taskID(), cur_tid, cur_qid);
			task_prog* parent_prog = global_prog_map.get(cur_tid);
			parent_prog->t_counter++; //no need lock since it happens before "completed -> ture"
		}
        global_prog_map.insert(totask->prog->tid, totask->prog);

        tc->L_timeout_mtx.rdlock();

        if(!tc->L_timeout.empty())
        {
            for(auto it = tc->L_timeout.begin(); it != tc->L_timeout.end(); ++it)
            {
                if((*it)->vid == vid)
                {
                    (*it)->task_queue.enqueue(totask);
                    tc->L_timeout_mtx.unlock();
                    return;
                }

                if((*it)->vid > vid) break;
            }
        }
        tc->L_timeout_mtx.unlock();

        // can't find one container's qid == id
        // that means we have to add a new container at the end of list
    
        // 1. allocate a new TimeOutTask_Container
        // TimeOutTask_Container * new_tc = new TimeOutTask_Container(vid);
        // new_tc->task_queue.enqueue(totask);

        // 2. upgrade lock to writelock, push into L_timeout
        // tc->L_timeout_mtx.wrlock();
        // tc->L_timeout.push_back(new_tc);
        // tc->L_timeout_mtx.unlock();
    }

    int get_timeout_min()
    {
        int ret;

        tc->L_timeout_mtx.rdlock();
        if(!tc->L_timeout.empty())
        {
            for(auto it = tc->L_timeout.begin(); it != tc->L_timeout.end(); ++it)
            {
                if(!(*it)->task_queue.empty())
                {
                    ret = (*it)->vid;
                    tc->L_timeout_mtx.unlock();
                    return ret;
                }
            }
        }
        tc->L_timeout_mtx.unlock();
        return -1;
    }

    bool get_timeout_task(TimeOutTask * &totask)
    {
        tc->L_timeout_mtx.rdlock();

        if(!tc->L_timeout.empty())
        {
            for(auto it = tc->L_timeout.begin(); it != tc->L_timeout.end(); ++it)
            {
                if(!(*it)->task_queue.empty())
                {
                    bool succ = (*it)->task_queue.dequeue(totask);
                    tc->L_timeout_mtx.unlock();
                    return succ;
                }
            }
        }
        tc->L_timeout_mtx.unlock();
        return false;
    }

    int compute(Task *task)
    {   
        Pattern *pattern = tc->pattern;

        cur_qid = tc->qid;

        gmatch_engine.set(&grami.pruned_graph, pattern);
        // gmatch_engine.generateGQLQueryPlan(task->u);
        // gmatch_engine.generateBN();
        gmatch_engine.matching_order = tc->all_matching_order[task->u]; // common
        gmatch_engine.bn = tc->all_bn[task->u]; // common
        gmatch_engine.bn_count = tc->all_bn_count[task->u]; // common
        gmatch_engine.edge_matrix = tc->edge_matrix; // common

        vector<VertexID> embedding(pattern->size());
        embedding[task->u] = task->v; // u_id -> v_id

        TimeOutTask * totask = new TimeOutTask(task, pattern->size());
        int value = gmatch_engine.execute(task->u, task->v_idx, embedding, *(totask->ms)); // use TimeOut

        
        if (value == MATCH_NOT_FOUND)
        {
            tc->non_cand_mtx[task->u].lock();
            pattern->non_candidates[task->u].insert(task->v);
            tc->non_cand_mtx[task->u].unlock();
            delete task;
            delete totask;
        }
        else if(value == MATCH_FOUND)
        {
            for (ui j = 0; j < pattern->size(); ++j)
            {
                // unordered_set<VertexID> & auto_nodes = pattern->autos[idx];
                // for (auto j: auto_nodes)
                {
                    VertexID matched_vertex = embedding[j];
                    
                    tc->domain_matches_mtx[j].lock();
                    tc->domain_matches[j].insert(matched_vertex);

                    // ========== Early Termination ==========
                    tc->vq_stops_refill_mtx[j].rdlock();
                    if(!tc->vq_stops_refill[j])
                    {
                        if(tc->domain_matches[j].size() >= grami.nsupport_)
                        {
                            tc->vq_stops_refill_mtx[j].unlock();
                            tc->vq_stops_refill_mtx[j].wrlock(); // upgrade rdlock to wrlock
                            tc->vq_stops_refill[j] = true; // prevent repetitve increment
                            tc->vq_stops_refill_mtx[j].unlock();

                            refill_check(j);
                        }
                        else tc->vq_stops_refill_mtx[j].unlock();
                    }
                    else tc->vq_stops_refill_mtx[j].unlock();
                    // ========== Early Termination done ==========

                    tc->domain_matches_mtx[j].unlock();
                }
            }
            delete task;
            delete totask;
        }
        else // value == TIMEOUT
        {
            tc->vq_stops_refill_mtx[task->u].rdlock();

            // only if domain of A is not done, plus current pattern is still frequent so far
            if(!tc->vq_stops_refill[task->u] && tc->frequent_tag) 
            {
                // throw current task into L_timeout
                tc->vq_stops_refill_mtx[task->u].unlock();

                totask->ms->depth = 1;
                totask->ms->embedding.resize(pattern->size());
                totask->ms->embedding[task->u] = task->v;

                totask->ms->idx_embedding.resize(pattern->size());
                totask->ms->idx_embedding[task->u] = task->v_idx;
                totask->ms->visited.insert(task->v);
                totask->skip = true; 
 
                add_timeout_task(totask, true);

#ifdef VERBOSE_MODE
                cout << "add TIMEOUT successfully!!!!!!!!!!" << endl;
#endif
            }
            else
            {
                tc->vq_stops_refill_mtx[task->u].unlock();
                delete totask;
            }
            delete task;
        }

        // gmatch_engine.clear_match();
        gmatch_engine.reset();

        return value;
    }

    /** Obsolete, not used !!!!!!!!! */
    int compute_timeout(Task *task, MatchingStatus &ms)
    {
        Pattern *pattern = tc->pattern;

        gmatch_engine.set(&grami.pruned_graph, pattern);
        gmatch_engine.generateGQLQueryPlan(task->u);
        gmatch_engine.generateBN();
        gmatch_engine.edge_matrix = tc->edge_matrix; // common

        vector<VertexID> embedding(pattern->size());
        embedding[task->u] = task->v; // u_id -> v_id
        
        int value = gmatch_engine.resumeSearch(task->u, ms, embedding); // use TimeOut

        if (value == MATCH_NOT_FOUND)
        {
            tc->non_cand_mtx[task->u].lock();
            pattern->non_candidates[task->u].insert(task->v);
            tc->non_cand_mtx[task->u].unlock();
        }
        else if(value == MATCH_FOUND)
        {
            for (ui j = 0; j < pattern->size(); ++j)
            {
                VertexID matched_vertex = embedding[j];
                
                tc->domain_matches_mtx[j].lock();
                tc->domain_matches[j].insert(matched_vertex);

                // ========== Early Termination ==========
                tc->vq_stops_refill_mtx[j].rdlock();
                if(!tc->vq_stops_refill[j])
                {
                    if(tc->domain_matches[j].size() >= grami.nsupport_)
                    {
                        tc->vq_stops_refill_mtx[j].unlock();
                        tc->vq_stops_refill_mtx[j].wrlock(); // upgrade rdlock to wrlock
                        tc->vq_stops_refill[j] = true; // prevent repetitve increment
                        tc->vq_stops_refill_mtx[j].unlock();

                        refill_check(j);
                    }
                    else tc->vq_stops_refill_mtx[j].unlock();
                }
                else tc->vq_stops_refill_mtx[j].unlock();
                // ========== Early Termination done ==========

                tc->domain_matches_mtx[j].unlock();
            }
        }

        gmatch_engine.clear_match();
        gmatch_engine.reset();

        return value;
    }

    bool backtrack(task_prog * prog, bool &tag)
    {   
        bool finish_all = false;

        prog->prog_lock.lock();
		prog->completed = true;
        bool last_tag = prog->decompose_tag; // TODO
		if(prog->t_counter == 0) 
        {
			prog->prog_lock.unlock();
			bool flag = true;
			while(flag)
			{
				TaskID parent_id = prog->parent_tid;

				global_prog_map.erase(prog->tid);
				delete prog;

				if(parent_id == -1)
				{
					flag = false;
                    finish_all = true;
                    tag = last_tag;
				}
                else
                {
					prog = global_prog_map.get(parent_id);
					prog->prog_lock.lock();
					prog->t_counter--;
					flag = (prog->completed && prog->t_counter == 0);

                    prog->decompose_tag = prog->decompose_tag || last_tag;
                    last_tag = prog->decompose_tag;

					prog->prog_lock.unlock();
				}
			}
		}
		else prog->prog_lock.unlock();

        return finish_all;
    }


#ifdef OPTIMIZED_MATCH

    /** Iterative style of dfs_decompose(.). No fast forward temporarily */
    void dfs_decompose_iterative(ui depth, Pattern *pattern, vector<Domain> &candidates, VertexID *matching_order, 
                    vector<VertexID> & embedding, vector<VertexID>& idx_embedding, unordered_set<VertexID> &visited, 
                    VertexID **bn, ui *bn_count, ui &num_decomposed_task, bool &found_result, TimeOutTask *totaskOri, 
                    Task *task, bool &proceed, ui ffdepth, vector<ui> &counter)
    {
        for (int i = 0; i < depth; ++i)
        {
            gmatch_engine.mvisited_array[embedding[matching_order[i]]] = true;
        }

        ui cur_depth = depth;
        ui max_depth = pattern->size();

        if (cur_depth == 0)
        {
            VertexID start_vertex = matching_order[cur_depth];
            counter[cur_depth] = 0; //#### TODO:
            gmatch_engine.midx_count[cur_depth] = pattern->get_cands()[start_vertex].size();

            for (ui i = 0; i < gmatch_engine.midx_count[cur_depth]; ++i) 
            {
                gmatch_engine.mvalid_candidate_idx[cur_depth][i] = i;
            }
        }
        else
        {
            counter[cur_depth] = 0; //#### TODO: No fast forward temporarily
            gmatch_engine.generateValidCandidateIndex(cur_depth, idx_embedding);
        }

        while(true)
        {
            while (counter[cur_depth] < gmatch_engine.midx_count[cur_depth])
            {
                if(found_result) break; // if found one result, we don't search any further.
                if(!proceed) break;

                ui valid_idx = gmatch_engine.mvalid_candidate_idx[cur_depth][counter[cur_depth]];
                ui u = matching_order[cur_depth];
                ui v = candidates[u].candidate[valid_idx];

                if (gmatch_engine.mvisited_array[v])
                {
                    counter[cur_depth] += 1;
                    continue;
                }

                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                gmatch_engine.mvisited_array[v] = true;
                counter[cur_depth] += 1;

                if (cur_depth == max_depth - 1)
                { 
                    found_result = true;
                    gmatch_engine.mvisited_array[v] = false;
                }
                else if (gmatch_engine.countElaspedTime() < TIMEOUT_THRESHOLD)
                {
                    cur_depth += 1;
                    counter[cur_depth] = 0;
                    gmatch_engine.generateValidCandidateIndex(cur_depth, idx_embedding);
                }
                else
                {
                    tc->vq_stops_refill_mtx[task->u].rdlock();
                    // only if domain of A is not done, plus current pattern is still frequent so far
                    if(!tc->vq_stops_refill[task->u] && tc->frequent_tag) 
                    {
                        // throw current task into L_timeout
                        tc->vq_stops_refill_mtx[task->u].unlock();

                        // decompose task
                        TimeOutTask * totask = new TimeOutTask(task, pattern->size());
                        totask->ms->depth = cur_depth+1;
                        totask->ms->embedding = embedding;
                        totask->ms->idx_embedding = idx_embedding;
                        totask->ms->visited = visited; // TODO: 
                        add_timeout_task(totask, false);

#ifdef VERBOSE_MODE
                        // cout << "add TIMEOUT again successfully!!!!!!!!!!" << endl;
#endif

                        num_decomposed_task ++;
                    }
                    else
                    {
                        tc->vq_stops_refill_mtx[task->u].unlock();

                        totaskOri->prog->prog_lock.lock();
                        totaskOri->prog->decompose_tag = true;
                        totaskOri->prog->prog_lock.unlock();
                        proceed = false;
                    }
                    gmatch_engine.mvisited_array[v] = false;
                }
            }
            cur_depth -= 1;
            if (cur_depth < depth)
                break;
            else
                gmatch_engine.mvisited_array[embedding[matching_order[cur_depth]]] = false;
        }

        for (int i = 0; i < depth; ++i)
        {
            gmatch_engine.mvisited_array[embedding[matching_order[i]]] = false;
        }
    }

#endif

    void dfs_decompose(ui depth, Pattern *pattern, vector<Domain> &candidates, VertexID *matching_order, 
                    vector<VertexID> & embedding, vector<VertexID>& idx_embedding, unordered_set<VertexID> &visited, 
                    VertexID **bn, ui *bn_count, ui &num_decomposed_task, bool &found_result, TimeOutTask *totaskOri, Task *task, bool &proceed, 
                    ui ffdepth, vector<ui> &counter, bool &skip)
    {

        struct timeb cur_time;
		double drun_time;

        if(depth == pattern->size()) 
        {
            found_result = true;
        }
        else
        {
            VertexID u = matching_order[depth];
            vector<ui> valid_candidate_idx;
            ui valid_candidate_cnt;
            if(depth == 0)
            {
                valid_candidate_cnt = pattern->get_cands()[u].size();
                for(ui i = 0; i < valid_candidate_cnt; ++i)
                {
                    valid_candidate_idx.push_back(i);
                }
            }
            else
            {
                gmatch_engine.generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                            bn[depth], bn_count[depth]);
            }

            ui i = skip? counter[depth]: 0;
            // ui i = 0;
            for( ; i < valid_candidate_cnt; ++i)
            {
                if(skip)
                {
                    if(depth == ffdepth-1)
                    {
                        // cout << "jump successfully!!!!!!!!!!" << endl;
                        skip = false;
                    }
                }

                ui valid_idx = valid_candidate_idx[i];
                VertexID v = pattern->get_cands()[u].candidate[valid_idx];

                if(visited.find(v) != visited.end()) continue;
                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited.insert(v);

                ftime(&cur_time);
                drun_time = cur_time.time-gmatch_engine.gtime_start.time+(double)(cur_time.millitm-gmatch_engine.gtime_start.millitm)/1000;

                if (skip) // Fast forward to recover last traversal place
                {
                    dfs_decompose(depth+1, pattern, candidates, matching_order, embedding, idx_embedding, 
                                visited, bn, bn_count, num_decomposed_task, found_result, totaskOri, task, proceed, ffdepth, counter, skip);
                }
                else if (drun_time < DECOMPOSE_TIME_THRESHOLD) 
                {
                    dfs_decompose(depth+1, pattern, candidates, matching_order, embedding, idx_embedding, 
                                        visited, bn, bn_count, num_decomposed_task, found_result, totaskOri, task, proceed, ffdepth, counter, skip);

                    if(found_result) return; // if found one result, we don't search any further.
                    if(!proceed) return;
                }
                else
                {
                    tc->vq_stops_refill_mtx[task->u].rdlock();
                    // only if domain of A is not done, plus current pattern is still frequent so far
                    if(!tc->vq_stops_refill[task->u] && tc->frequent_tag) 
                    {
                        // throw current task into L_timeout
                        tc->vq_stops_refill_mtx[task->u].unlock();

                        // decompose task
                        TimeOutTask * totask = new TimeOutTask(task, pattern->size());
                        totask->ms->depth = depth+1;
                        totask->ms->embedding = embedding;
                        totask->ms->idx_embedding = idx_embedding;
                        totask->ms->visited = visited;
                        add_timeout_task(totask, false);
#ifdef VERBOSE_MODE
                        // cout << "add TIMEOUT again successfully!!!!!!!!!!" << endl;
#endif

                        num_decomposed_task ++;
                    }
                    else
                    {
                        tc->vq_stops_refill_mtx[task->u].unlock();

                        totaskOri->prog->prog_lock.lock();
                        totaskOri->prog->decompose_tag = true;
                        totaskOri->prog->prog_lock.unlock();

                        proceed = false;
                        return;
                    }
                }
                visited.erase(v);
            }
        }
    }

    /**
     * If all decomposed tasks under (task->u, task->v) are done, return true
     * Otherwise, return false
     */
    bool compute_decomposed(TimeOutTask *totask, MatchingStatus &ms) 
    {   
        Pattern *pattern = tc->pattern;
        Task *task = &(totask->task);

        tc->domain_matches_mtx[task->u].lock(); // first lock match, then lock vq_stop_refill
        tc->vq_stops_refill_mtx[task->u].rdlock();
        if(tc->vq_stops_refill[task->u] || !tc->frequent_tag || tc->domain_matches[task->u].find(task->v) != tc->domain_matches[task->u].end())
        {
            totask->prog->prog_lock.lock();
            totask->prog->decompose_tag = true;
            totask->prog->prog_lock.unlock();

            tc->vq_stops_refill_mtx[task->u].unlock();
            tc->domain_matches_mtx[task->u].unlock();
        }
        else
        {
            tc->vq_stops_refill_mtx[task->u].unlock();
            tc->domain_matches_mtx[task->u].unlock();

            cur_qid = totask->get_qid();
            cur_tid = totask->get_id();

            gmatch_engine.set(&grami.pruned_graph, pattern);
            // gmatch_engine.generateGQLQueryPlan(task->u); 
            // gmatch_engine.generateBN();
            gmatch_engine.matching_order = tc->all_matching_order[task->u]; // common
            gmatch_engine.bn = tc->all_bn[task->u]; // common
            gmatch_engine.bn_count = tc->all_bn_count[task->u]; // common
            gmatch_engine.edge_matrix = tc->edge_matrix; // common

            // DFS trackers
            ui num_decomposed_task = 0;
            bool found_result = false;
            bool proceed = true;
            
            // Start searching ...
            ftime(&gmatch_engine.gtime_start);

#ifdef OPTIMIZED_MATCH
            dfs_decompose_iterative(ms.depth, pattern, pattern->get_cands(), gmatch_engine.matching_order, ms.embedding, ms.idx_embedding, ms.visited, 
                            gmatch_engine.bn, gmatch_engine.bn_count, num_decomposed_task, found_result, totask, task, proceed, ms.ffdepth, ms.counter);
#else
            dfs_decompose(ms.depth, pattern, pattern->get_cands(), gmatch_engine.matching_order, ms.embedding, ms.idx_embedding, ms.visited, 
                            gmatch_engine.bn, gmatch_engine.bn_count, num_decomposed_task, found_result, totask, task, proceed, ms.ffdepth, ms.counter, 
                            totask->skip); // with divided-and-conquer, not in gmatch.h
#endif


            if(found_result)
            {
                // find a result
                for (ui j = 0; j < pattern->size(); ++j)
                {
                    // unordered_set<VertexID> & auto_nodes = pattern->autos[idx];
                    // for (auto j: auto_nodes)
                    {
                        VertexID matched_vertex = ms.embedding[j];
                        
                        tc->domain_matches_mtx[j].lock();
                        tc->domain_matches[j].insert(matched_vertex);

                        // ========== Early Termination ==========
                        
                        tc->vq_stops_refill_mtx[j].rdlock();
                        if(!tc->vq_stops_refill[j])
                        {
                            if(tc->domain_matches[j].size() >= grami.nsupport_)
                            {
                                tc->vq_stops_refill_mtx[j].unlock();
                                tc->vq_stops_refill_mtx[j].wrlock(); // upgrade rdlock to wrlock
                                tc->vq_stops_refill[j] = true; // prevent repetitve increment
                                tc->vq_stops_refill_mtx[j].unlock();

                                refill_check(j);
                            }
                            else tc->vq_stops_refill_mtx[j].unlock();
                        }
                        else tc->vq_stops_refill_mtx[j].unlock();

                        // ========== Early Termination done ==========

                        tc->domain_matches_mtx[j].unlock();
                    }
                }
            }
            // gmatch_engine.clear_match();
            gmatch_engine.reset();
        }

        bool check;
        bool all_finished = backtrack(totask->prog, check);

        if(all_finished && !check)
        {
            tc->domain_matches_mtx[task->u].lock();
            if (tc->domain_matches[task->u].find(task->v) == tc->domain_matches[task->u].end())
            {
                tc->domain_matches_mtx[task->u].unlock();
                tc->non_cand_mtx[task->u].lock();
                pattern->non_candidates[task->u].insert(task->v);
                tc->non_cand_mtx[task->u].unlock();
            }
            else tc->domain_matches_mtx[task->u].unlock();
        }
        
        return all_finished;
    }

    void delete_tc()  // called when q is finished
    {
        // delete tc
        task_container *tc_delete;
        activeQ_lock.wrlock();
        auto it = activeQ_list.begin();

        int delete_qid = tc->qid;

        for(; it != activeQ_list.end(); ++it)
        {
            tc_delete = *it;
            if(tc_delete->qid == delete_qid) 
            {
                break;
            }
        }
        assert(it != activeQ_list.end());
        activeQ_list.erase(it);
        activeQ_num--;

        activeQ_lock.unlock();

        delete tc;
    }

    void remove_non_cands()
    {
        Pattern *pattern = tc->pattern;
        for(ui i = 0; i < pattern->size(); ++i) 
        {   
            auto & cands = pattern->get_cands()[i].candidate;
            auto & non_cands = pattern->non_candidates[i];
            vector<VertexID> new_cands;
            for(auto it = cands.begin(); it !=cands.end(); ++it)
            {
                if(non_cands.find(*it) == non_cands.end())
                {
                    new_cands.push_back(*it);
                }
            }
            cands.swap(new_cands);
        }
    }

    bool check_freq_and_extend_express(task_container * tc_)
    {
        Pattern *pattern = tc_->pattern;

        results_counter[thread_id]++;
        // cout << "##### " << pattern->size() << " ########" << endl;

#ifdef VERBOSE_MODE
        cout << "#####" << endl;
        cout << pattern->toString();
        cout << "#####" << endl;
#endif 

        // cout << pattern->toString() << endl;

        subPattern key;
        // pattern->toCache(key);
        key.pattern.vertices_ = pattern->vertices_;
        key.pattern.edge2vertex = pattern->edge2vertex;
        cache_mtx.wrlock();
        cache[key].swap(pattern->non_candidates);
        cache_mtx.unlock();

        PatternPVec ext_pattern_vec;
        grami.extend((*pattern), ext_pattern_vec);
        pattern->prog->children_cnt = ext_pattern_vec.size();

        // special case, otherwise memory leak!!
        if(pattern->prog->children_cnt == 0)
        {
            delete pattern->prog;
        }

        for (auto it = ext_pattern_vec.begin(); it != ext_pattern_vec.end(); ++it)
        {
            // add into datastack
            task_container *new_tc = new task_container(qid++);

            Pattern * child_pattern = *it;
            child_pattern->parent_prog = pattern->prog;
            child_pattern->non_candidates.resize(child_pattern->size());
            new_tc->pattern = child_pattern;
            
            data_stack.enstack(new_tc);
        }
        return true;
    }   

    // will remove non-candidates from domain
    bool check_freq_and_extend()  // called when q is finished normally (not early termination)
    {   
        Pattern *pattern = tc->pattern;

        // get min frequency among v_q's
        // int freq = tc->domain_matches[0].size();

        // for (ui i = 1; i < pattern->size(); ++i)
        // {
        //     if (freq > tc->domain_matches[i].size())
        //     {
        //         freq = tc->domain_matches[i].size();
        //     }
        // }

        // // extend when frequent
        // if(freq >= grami.nsupport_)
        
        {
            results_counter[thread_id]++;

            // cout << "##### " << pattern->size() << " ########" << endl;

            
            // if (pattern->size()*(pattern->size()-1) /2 == pattern->get_nedges() && pattern->distinct_labels())
            // {
#ifdef VERBOSE_MODE
                cout << "~~~~~" << endl;
                cout << pattern->toString();
                
                // cout << tc->domain_matches[0].size() << " " <<  tc->domain_matches[1].size() << " " << tc->domain_matches[2].size() << endl;
                cout << "~~~~~" << endl;
#endif
            // }
            

            // cout << pattern->toString() << endl;

            remove_non_cands();

            // ===== push down pruning ======

            // string code = pattern->toString();

            subPattern key;
            // pattern->toCache(key);
            key.pattern.vertices_ = pattern->vertices_;
            key.pattern.edge2vertex = pattern->edge2vertex;
            cache_mtx.wrlock();
            cache[key].swap(pattern->non_candidates);
            cache_mtx.unlock();

            // ===== push down pruning done ======

            PatternPVec ext_pattern_vec;
            grami.extend((*pattern), ext_pattern_vec);
            pattern->prog->children_cnt = ext_pattern_vec.size();

            // special case, otherwise memory leak!!
            if(pattern->prog->children_cnt == 0)
            {
                delete pattern->prog;
            }

            for (auto it = ext_pattern_vec.begin(); it != ext_pattern_vec.end(); ++it)
            {
                // add into datastack
                task_container *new_tc = new task_container(qid++);

                Pattern * child_pattern = *it;
                child_pattern->parent_prog = tc->pattern->prog;
                child_pattern->non_candidates.resize(child_pattern->size());
                new_tc->pattern = child_pattern;
                
                data_stack.enstack(new_tc);
            }
            return true;

        } // the last child will delete tc->pattern->prog, otherwise
        // else // if not frequent, there's no need to keep progress and itself
        // {
        //     return false;
        // }
    }

    bool refill_Q_adjust(bool unlock)
    {
        bool exit = false;
        while(!exit)
        {
            tc->vq_stops_refill_mtx[tc->next_vq].rdlock();
            if(tc->vq_stops_refill[tc->next_vq])
            {
                tc->vq_stops_refill_mtx[tc->next_vq].unlock();

                tc->next_vq ++;
                if(tc->next_vq == tc->pattern->size()) // if last one, exit
                {
                    if(unlock) tc->refill_mtx.unlock();
                    return false;
                }
            }
            else
            {
                tc->vq_stops_refill_mtx[tc->next_vq].unlock();
                exit = true;
            }
        }

        return true;
    }

    bool refill_Q()
    {
        if(!tc->frequent_tag) return false;
    
        vector<Task *> tmp_vector;
        Pattern * pattern = tc->pattern;

        tc->refill_mtx.lock();

        if(!tc->nothing_to_refill())
        {

            if(!refill_Q_adjust(true))
            {
                return false;
            }

            while(tmp_vector.size() < MINI_BATCH_NUM)
            {
                VertexID & idx = tc->next_domainPos[tc->next_vq];
                Domain & vq_candidates = pattern->get_cands()[tc->next_vq];
                
                Task * t = new Task(tc->next_vq, vq_candidates[idx], idx);
                tmp_vector.push_back(t);

                // update
                idx++;
                if(idx == vq_candidates.size())
                {
                    tc->next_vq ++;
                    if(tc->nothing_to_refill())
                        break;
                    if(!refill_Q_adjust(false))
                        break;
                }
            }

            tc->refill_mtx.unlock();

            // spawn tasks
            tc->Q_domain_mtx.lock();
            tc->Q_domain.insert(tc->Q_domain.end(), tmp_vector.begin(), tmp_vector.end());

            // ========== update domain_front ===========
            
            tc->domain_front = tc->Q_domain.front()->u;

            // ========== update domain_front done ===========

            tc->Q_domain_mtx.unlock();

            return true;
        }
        else 
        {
            tc->refill_mtx.unlock();
            return false;
        }
    }

    void refill_check(VertexID id)
    {
        // first do not think about conque's deletion
        // ====== delete timeoutQ ======

        tc->L_timeout_mtx.wrlock();
        auto it = tc->L_timeout.begin();
        for( ; it != tc->L_timeout.end(); ++it)
        {
            if((*it)->vid == id)
            {
                delete *it;
                tc->L_timeout.erase(it);
                break;
            }
        }
        tc->L_timeout_mtx.unlock();

        // ==================================

        tc->domain_done_mtx.lock();
        tc->domain_done ++;
    
        if (tc->domain_done == tc->pattern->size())
        {
            tc->normal_exit = true;

            // fout[thread_id] << "QID: " << tc->qid << " set normal_exit to true" << endl;
        }
        tc->domain_done_mtx.unlock();
    }

    void postprocess(VertexID u)
    {   
        tc->end_mtx.lock();

        if(tc->frequent_tag)
        {
            tc->domain_increment(u);
            
            // if(tc->domain_finished[u] == tc->pattern->get_cands()[u].size() && tc->domain_matches[u].size() >= grami.nsupport_) 
            // {
            //     tc->vq_stops_refill_mtx[u].wrlock(); // upgrade rdlock to wrlock
            //     tc->vq_stops_refill[u] = true; // prevent repetitve increment
            //     tc->vq_stops_refill_mtx[u].unlock();

            //     refill_check(u);
            // }

            // ======= update tc->frequent_tag =======
            if(tc->domain_early_exit(u))
            {
                tc->frequent_tag = false; // after refill and compute!!!

                // fout[thread_id] << "QID: " << tc->qid << " set frequent_tag to false" << endl;
            }
            // ======= update tc->frequent_tag done ======= 
        }
        tc->end_mtx.unlock();
    }

    bool get_and_process_tasks()
    {   
        Task *task = NULL;
        TimeOutTask *totask = NULL;
        bool succ;

        bool block = false;

        activeQ_lock.rdlock();

        if(!activeQ_list.empty())
        {
            bool refilled = false;
            for(auto it = activeQ_list.begin(); it != activeQ_list.end(); ++it)
            {
                tc = *it;

                if(tc->comper_counter_mtx.try_lock()) // prevent deadlock
                {   
                    if((tc->normal_exit || !tc->frequent_tag) && tc->comper_counter >= 1)
                    {
                        // won't interfere tc's deletion
                        tc->comper_counter_mtx.unlock();
                        continue;
                    }
                    else
                    {
                        tc->comper_counter ++;
                        tc->comper_counter_mtx.unlock();
                    }
                }
                else 
                {
                    block = true;
                    continue;
                }

                bool task_obtained = false;

                if(tc->frequent_tag && !tc->normal_exit)
                {
                    // 1. nothing in Q_timeout
                    // 2. Q_domain hasn't done yet
                    int timeout_front = get_timeout_min();
                    if(timeout_front == -1 || (tc->domain_front != -1 && tc->domain_front <= timeout_front)) // < is very unlikely, <= ensures priority of Q_domain.
                    {
                        if(tc->Q_domain_mtx.try_lock())
                        {
                            if(tc->Q_domain.size() < RT_THRESHOLD_FOR_REFILL)
                            {
                                tc->Q_domain_mtx.unlock();

                                refilled = refill_Q();
                            }
                            else
                            {
                                tc->Q_domain_mtx.unlock();
                            }

                            if(tc->Q_domain_mtx.try_lock())
                            {
                                if(!tc->Q_domain.empty())
                                {
                                    task = tc->Q_domain.front();
                                    tc->Q_domain.pop_front();

                                    // ========== update domain_front ===========
                                    if(!tc->Q_domain.empty())
                                    {
                                        tc->domain_front = tc->Q_domain.front()->u;
                                    }
                                    else
                                    {
                                        tc->domain_front = -1; // -1 means Q_domain is temporarily empty
                                    }
                                    // ========== update domain_front done ===========

                                    tc->Q_domain_mtx.unlock();

                                    activeQ_lock.unlock();

                                    task_obtained = true;
                                    
                                    VertexID u = task->u;

                                    // ====== Early Termination =======
                                    // bool skip_compute = false;
                                    // tc->domain_matches_mtx[u].lock();

                                    // if(tc->domain_matches[u].find(v) != tc->domain_matches[u].end())
                                    // {
                                    //     skip_compute = true;
                                    // }
                                    // tc->domain_matches_mtx[u].unlock();
                                    // // ====== Early Termination done =======

                                    // if(!skip_compute)

                                    int value = compute(task);
                                    
                                    if(value != TIMEOUT)
                                    {
                                        postprocess(u);
                                    }

                                }
                                else
                                {
                                    tc->Q_domain_mtx.unlock();
                                }
                            }
                        }
                    }
                    else // tc->domain_front > timeout_front
                    {

                        TimeOutTask * totask;
                        succ = get_timeout_task(totask);
        
                        if(succ) 
                        {
                            activeQ_lock.unlock();

                            Task * task = &(totask->task);
                            task_obtained = true;

                            VertexID u = task->u;

                            // =========== Decompose pruning ============
                            
                            bool existence = true;
                            
                            //**                     
                            // if totask is on first level
                            if(totask->skip)
                            {
                                vector<vector<subPattern> > maps = tc->maps;
                                for(ui i = 0; i < maps.size(); ++i)
                                {
                                    vector<subPattern> &edge_removed = maps[i];
                                    for(auto it = edge_removed.begin(); it != edge_removed.end(); ++it)
                                    {
                                        Pattern & key = (*it).pattern;
                                        vector<VertexID> & subgraph_mappings = (*it).mapping;
                                        ui subgraph_size = key.size();
                                        int map_index = grami.search_mapping(subgraph_mappings, u);
                                        if(map_index == -1) continue;
                                        vector<VertexID> embedding(subgraph_size);
                                        embedding[map_index] = task->v_idx;

                                        key.prog = new PatternProgress;

                                        key.get_cands().resize(subgraph_size);
                                        for (ui l = 0; l < subgraph_size; ++l)
                                        {
                                            key.get_cands()[l].candidate = tc->pattern->get_cands()[subgraph_mappings[l]].candidate; // deep copy
                                        }
                                        exist_engine.set(&grami.pruned_graph, &key);
                                        existence = exist_engine.searchExistence(map_index, task->v_idx, embedding);
                                        exist_engine.clear_all();

                                        delete key.prog;

                                        if (!existence)
                                        {
                                            tc->non_cand_mtx[task->u].lock();
                                            tc->pattern->non_candidates[u].insert(task->v);
                                            tc->non_cand_mtx[task->u].unlock();

                                            postprocess(u);
                                            delete totask;
                                            
                                            break;
                                        }
                                    }
                                    if (!existence) break;
                                }
                            }
                            //*/
                            // =========== Decompose pruning done ============
                            
                            if(existence)
                            {
                                if(compute_decomposed(totask, *(totask->ms)))
                                {
                                    // if all timeout tasks under (u,v) are done
                                    postprocess(u);

                                    // delete totask->task;
                                }

                                delete totask;
                            }
                        }
                    }
                } // end if

                // ======= delete =======
                tc->comper_counter_mtx.lock();
                tc->comper_counter--;
                if(tc->comper_counter == 0)
                {
                    if(tc->normal_exit) // activeQ is unlock
                    {
                        if(!task_obtained) activeQ_lock.unlock();
                        if(!check_freq_and_extend())
                        {
                            cout << "check_freq_and_extend() #########!!!!!!!!!!" << "QID: " << tc->qid<< endl;
                            // fout[thread_id] << "delete qid = " << tc->qid << "'s prog" << endl;
                            
                            delete tc->pattern->prog;
                        }

                        // fout[thread_id] << "delete QID: " << tc->qid << " in normal exit" << endl;

                        delete_tc();
                        return true;
                    }
                    else if(!tc->frequent_tag)
                    {
                        if(!task_obtained) activeQ_lock.unlock();
                        // fout[thread_id] << "delete qid = " << tc->qid << "'s prog" << endl;

                        delete tc->pattern->prog;

                        // fout[thread_id] << "delete QID: " << tc->qid << " in frequent_tag" << endl;

                        delete_tc();
                        return true;
                    }
                    else tc->comper_counter_mtx.unlock();
                }
                else tc->comper_counter_mtx.unlock();
                
                if(task_obtained) // activeQ is unlock
                {
                    return true;
                }
                else if(refilled)
                {
                    activeQ_lock.unlock();
                    return true;
                }
            }
        }
        activeQ_lock.unlock();

        // reaching here means we don't get any task 

        succ = false;
        if(!data_stack.empty())
        {
            activeQ_lock.wrlock();
            if(activeQ_num < activeQ_list_capacity)
            {
                activeQ_num++;
                activeQ_lock.unlock();

                task_container *tc_new;
                succ = data_stack.destack(tc_new); // another comper may have taken tc from data_stack

                if(succ)
                {
                    tc_new->init();

                    // ====== compute automorphisms =======
                    // tc_new->pattern->autos.resize(tc_new->pattern->size());
                    // Graph dataGraphAuto;
                    // dataGraphAuto.toGraph(*(tc_new->pattern));
                    
                    // gmatch_engine.set(&dataGraphAuto, tc_new->pattern);
                    // gmatch_engine.executeAuto(tc_new->pattern->autos);
                    // gmatch_engine.clear_all();
                    // tc_new->pattern->clear_candidate();
                    // ====== compute automorphisms done =======


                    // ====== push down pruning ==========

                    
                    decomposer.decompose(*(tc_new->pattern), tc_new->maps);

                    for (ui i = 0; i < tc_new->maps.size(); ++i)
                    {
                        vector<subPattern> &edge_removed = tc_new->maps[i];
                        for (auto it = edge_removed.begin(); it != edge_removed.end(); ++it)
                        {
                            subPattern & key = (*it);

                            cache_mtx.rdlock();
                            auto pos = cache.find(key);
                            if (pos != cache.end())
                            {
                                vector<VertexID> &subgraph_mappings = (*it).mapping;
                                VtxSetVec &node_non_candidates = pos->second;

                                cache_mtx.unlock();

                                for (ui j = 0; j < subgraph_mappings.size(); ++j)
                                {
                                    unordered_set<VertexID> &non_cans = tc_new->pattern->non_candidates[subgraph_mappings[j]];
                                    unordered_set<VertexID> &non_cans_check = node_non_candidates[j];
                                    for (auto it2 = non_cans_check.begin(); it2 != non_cans_check.end(); ++it2)
                                    {
                                        non_cans.insert(*it2);
                                    }
                                }
                            }
                            else cache_mtx.unlock();
                        }
                    }
                
                    // ====== push down pruning done ==========

                    // ============= unique label pruning ====================

                    
                    if(tc_new->pattern->distinct_labels() && tc_new->pattern->is_acyclic())
                    {
                        gmatch_engine.set(&grami.pruned_graph, tc_new->pattern);
                        bool is_freq = gmatch_engine.filterToConsistency(grami.nsupport_);
                        gmatch_engine.reset();

                        PatternProgress * pattern_prog = tc_new->pattern->parent_prog;
                        if(pattern_prog != NULL)
                        {
                            pattern_prog->children_mtx.lock();
                            pattern_prog->children_cnt--;
                            if(pattern_prog->children_cnt == 0)
                            {
                                delete pattern_prog;
                            }
                            else
                            {
                                pattern_prog->children_mtx.unlock();
                            }
                        }

                        if(is_freq)
                        {
                            check_freq_and_extend_express(tc_new);
                            delete tc_new;
                        }
                        else
                        {
                            delete tc_new->pattern->prog;
                            delete tc_new;
                        }
                        activeQ_lock.wrlock();
                        activeQ_num--;
                        activeQ_lock.unlock();
                    }
                    // ============= unique label pruning done ====================
                    else
                    {
                        gmatch_engine.set(&grami.pruned_graph, tc_new->pattern);

                        bool keep = gmatch_engine.DPisoFilter(false, grami.nsupport_); // degree-based pruning

                        // delete parent pattern
                        PatternProgress * pattern_prog = tc_new->pattern->parent_prog;
                        if(pattern_prog != NULL)
                        {
                            pattern_prog->children_mtx.lock();
                            pattern_prog->children_cnt--;
                            if(pattern_prog->children_cnt == 0)
                            {
                                delete pattern_prog;
                            }
                            else
                            {
                                pattern_prog->children_mtx.unlock();
                            }
                        }

                        if(keep)
                        {
                            gmatch_engine.buildTable(tc_new->edge_matrix);
                            for (ui i = 0; i < tc_new->pattern->size(); ++i)
                            {
                                gmatch_engine.generateGQLQueryPlan(i, tc_new->all_matching_order[i]);
                                gmatch_engine.generateBN(tc_new->all_bn_count[i], tc_new->all_bn[i], tc_new->all_matching_order[i]);
                            }
                            gmatch_engine.reset();

                            activeQ_lock.wrlock();
                            activeQ_list.push_back(tc_new);
                            activeQ_lock.unlock();
                        }
                        else // current pattern is not frequent because at least one vertex has domain size < support 
                        {
                            activeQ_lock.wrlock();
                            activeQ_num--;
                            activeQ_lock.unlock();

                            delete tc_new->pattern->prog;
                            delete tc_new;
                        }
                    }
                }
                else // get nothing from data_stack
                {
                    activeQ_lock.wrlock();
                    activeQ_num--;
                    activeQ_lock.unlock();
                }
            }
            else
            {
                activeQ_lock.unlock();
            }
        }
        return succ || block;
    }

    void run()
    {
        while (global_end_label == false) // Otherwise, thread terminates
        {
            // Process task or batch of tasks
            bool task_found = get_and_process_tasks();

            // Means that the Queues are empty
            if (!task_found)
            {
                unique_lock<mutex> lck(mtx_go);
                ready_go = false;
                global_num_idle++;
                while (!ready_go)
                {
                    cv_go.wait(lck);
                }
                global_num_idle--;
            }
        }
    }

    void start(int thread_id)
    {
        this->thread_id = thread_id;
        main_thread = thread(&Comper::run, this);
    }
};