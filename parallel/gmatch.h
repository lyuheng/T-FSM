#pragma once

#include "graph.h"
#include "leapfrogjoin.h"
#include "intersection/computesetintersection.h"

#include "global.h"

#include <algorithm>
#include <sys/timeb.h>

#define INVALID_VERTEX_ID 1000000000


struct MatchingStatus
{
    bool init;
    VertexID vid;
    ui v_idx;
    ui depth;

    ui ffdepth; // fast-forward detph

    vector<ui> counter;

    // matching context 
    // used for decomposed tasks 
    vector<VertexID> embedding;
    vector<ui> idx_embedding;
    unordered_set<VertexID> visited;

    MatchingStatus(ui size): init(false), ffdepth(1) 
    {
        counter.assign(size, 0);
    }
};

class GMatchEngine
{
public:
    Pattern *query_graph;
    Graph *data_graph;

    // =====================================

    VertexID *bfs_order, *matching_order;
    TreeNode *tree;

    Edges ***edge_matrix;

    VertexID **bn;
    ui *bn_count;
    
    // ====================================
    // Added from iterative-style backtrack

#ifdef OPTIMIZED_MATCH
    ui midx_count[MAX_PATTERN_SIZE];
    ui *mtemp_buffer;
    ui *mvalid_candidate_idx[MAX_PATTERN_SIZE];
    bool *mvisited_array;
#endif

    // =========== add from GraMi ===========

    timeb gtime_start;

    // ======================================

    GMatchEngine()
    {
        bfs_order = NULL;
        matching_order = NULL;
        tree = NULL;

#ifdef OPTIMIZED_MATCH

        mtemp_buffer = new ui[grami.pruned_graph.maxLabelFreq];

        for (ui i = 0; i < MAX_PATTERN_SIZE; ++i)
        {
            mvalid_candidate_idx[i] = new ui[grami.pruned_graph.maxLabelFreq];
        }
        mvisited_array = new bool[grami.pruned_graph.size()];
        memset(mvisited_array, false, sizeof(bool)*grami.pruned_graph.size());
#endif
    }
    
    ~GMatchEngine()
    {
#ifdef OPTIMIZED_MATCH
        if (mtemp_buffer != NULL)
            delete[] mtemp_buffer;
        for (ui i = 0; i < MAX_PATTERN_SIZE; ++i)
        {
            delete[] mvalid_candidate_idx[i];
        }

        delete[] mvisited_array;
#endif
    }

    void set(Graph *data_graph_, Pattern *query_graph_) 
    {
        data_graph = data_graph_;
        query_graph = query_graph_;
    }

    /** find one result, with timeout */
    int execute(VertexID u, ui v_idx, vector<VertexID> &embedding, MatchingStatus &ms);

    /** find one result, without timeout */
    int execute(VertexID u, ui v_idx, vector<VertexID> &embedding);

    /** find all results, without timeout */
    void executeAuto(vector<unordered_set<VertexID>> &results);

    // ========== TimeOut search ==========

    bool searchExistence(VertexID u, ui v_idx, vector<VertexID> &embedding);
    int resumeSearch(VertexID u, MatchingStatus &ms, vector<VertexID> &embedding);

    // ========== Filtering ==============

    bool LDF_pruning(bool automorphism);
    bool DPisoFilter(bool automorphism, ui support = 0);
    void generateDPisoFilterPlan();
    VertexID selectDPisoStartVertex();
    void bfsTraversal(VertexID root_vertex);
    bool pruneCandidates(VertexID query_vertex, VertexID *pivot_vertices, ui pivot_vertices_count, ui *flag, ui *updated_flag);
    void compactCandidates();
    bool isCandidateSetValid();
    bool isCandidateSetValid(ui support);

    bool filterToConsistency(ui support = 0);

    // ========== Ordering ===============

    void generateGQLQueryPlan(int start_vertex);
    void generateGQLQueryPlan(int start_vertex, ui *&matching_order_);
    VertexID selectGQLStartVertex();
    void updateValidVertices(VertexID query_vertex, vector<bool> &visited, vector<bool> &adjacent);
    void generateBN();
    void generateBN(ui *&bn_count_, ui **&bn_, ui *&matching_order_);
    void buildTable();

    void buildTable(Edges ***&edge_matrix_);

    // =========== Backtracking ================

    /** If finding one result, return true; otherwise return false. Use timeout and record status */
    int backtrack(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding, unordered_set<VertexID> &visited,
                    vector<ui> &counter, MatchingStatus &ms);
    
    /** If finding one result, return true; otherwise return false. w/o using timeout and recording status */
    int backtrack(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                                unordered_set<VertexID> &visited, vector<ui> &counter);

#ifdef OPTIMIZED_MATCH
    /** Iterative-style backtrack until finding one result, with timeout and recording status */
    int LFTJ(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                                unordered_set<VertexID> &visited, vector<ui> &counter, MatchingStatus &ms);
#endif

    /** Find all results */
    void backtrackAll(ui depth, VertexID *embedding, ui *idx_embedding, unordered_set<VertexID> &visited, 
                        VtxSetVec &results);
    
    /** Resume backtracking after TimeOut */
    int backtrackResume(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                        unordered_set<VertexID> &visited, vector<ui> &counter, MatchingStatus &ms, bool &skip);

    /** Find local candidates via intersection */
    void generateValidCandidateIndex(VertexID u, ui *idx_embedding, ui &valid_candidate_cnt, 
                                    vector<ui> &valid_candidate_index, VertexID *bn, ui bn_cnt);

    /** Find local candidates via intersection */
    void generateValidCandidateIndex(VertexID u, vector<ui> &idx_embedding, ui &valid_candidate_cnt, 
                                    vector<ui> &valid_candidate_index, VertexID *bn, ui bn_cnt);

#ifdef OPTIMIZED_MATCH
    /** Find local candidates via intersection */
    void generateValidCandidateIndex(ui depth, vector<ui> &idx_embedding);
#endif

    // ============ Utility ====================
    double countElaspedTime();

#ifdef OPTIMIZED_MATCH
    void resetVisitedArray(ui depth, vector<VertexID> &embedding);
#endif

    // ============ Clear ======================
    void clear_all();
    void clear_match();
    void clear_filter();
    void clear_table();
    void reset();

};

/**
 * Match u(query_v) to v(data_v) with timeout
 */
int GMatchEngine::execute(VertexID u, ui v_idx, vector<VertexID> &embedding, MatchingStatus &ms)
{   
    unordered_set<VertexID> visited;
    vector<ui> idx_embedding(query_graph->size());
    idx_embedding[u] = v_idx;
    visited.insert(embedding[u]);
    vector<ui> counter(query_graph->size(), 0);

    ftime(&gtime_start);  // record starting time

#ifdef OPTIMIZED_MATCH
    int value = LFTJ(1, embedding, idx_embedding, visited, counter, ms);
#else
    int value = backtrack(1, embedding, idx_embedding, visited, counter, ms);
#endif

    return value;
}

/**
 * Match u(query_v) to v(data_v) without timeout
 */
int GMatchEngine::execute(VertexID u, ui v_idx, vector<VertexID> &embedding)
{   
    unordered_set<VertexID> visited;
    vector<ui> idx_embedding(query_graph->size());
    idx_embedding[u] = v_idx;
    visited.insert(embedding[u]);
    vector<ui> counter(query_graph->size(), 0);

    int value = backtrack(1, embedding, idx_embedding, visited, counter);

    return value;
}

void GMatchEngine::executeAuto(vector<unordered_set<VertexID>> &results)
{
    DPisoFilter(true);
    generateGQLQueryPlan(-1);
    generateBN();
    buildTable();

    unordered_set<VertexID> visited;
    VertexID *embedding = new VertexID[query_graph->size()];
    ui *idx_embedding = new ui[query_graph->size()];
    backtrackAll(0, embedding, idx_embedding, visited, results);
    
    delete[] embedding;
    delete[] idx_embedding;
}

bool GMatchEngine::searchExistence(VertexID u, ui v_idx, vector<VertexID> &embedding)
{   
    generateDPisoFilterPlan();
    
    generateGQLQueryPlan(u);
    generateBN();
    buildTable();

    unordered_set<VertexID> visited;
    vector<ui> idx_embedding(query_graph->size());
    idx_embedding[u] = v_idx;
    visited.insert(embedding[u]);
    vector<ui> counter(query_graph->size(), 0);
    MatchingStatus ms(query_graph->size());

    ftime(&gtime_start);  // record starting time
    
    int value = backtrack(1, embedding, idx_embedding, visited, counter, ms);

    clear_filter();

    if(value == -2)
        return false;
    else 
        return true;
}

int GMatchEngine::resumeSearch(VertexID u, MatchingStatus &ms, vector<VertexID> &embedding)
{   
    vector<ui> idx_embedding(query_graph->size());
    idx_embedding[u] = ms.v_idx;
    unordered_set<VertexID> visited;
    visited.insert(embedding[u]);

    bool skip = true;
    int value = backtrackResume(1, embedding, idx_embedding, visited, ms.counter, ms, skip);

    return value;
}

/** if(parent == NULL) 
 *      nodesByLabel
 *  else
 *      copy parent candidate
*/
bool GMatchEngine::LDF_pruning(bool automorphism)
{   
    query_graph->get_cands().resize(query_graph->size());

    unordered_map<vLabel, int> *nlf = new unordered_map<vLabel, int>[query_graph->size()];

    query_graph->build_nlf(nlf);

    if(!automorphism) {
        for(ui i = 0; i < query_graph->size(); ++i) {
            vLabel vlabel = query_graph->get_vertex_label(i);
            int degree = query_graph->get_vertex_degree(i);
            
            if(query_graph->parent_prog == NULL) {
                vector<VertexID> &data_vertices = data_graph->nodesByLabel[vlabel]; 

                for(ui j = 0; j < data_vertices.size(); ++j) {
                    VertexID data_vertex = data_vertices[j];
                    if (data_graph->get_vertex_degree(data_vertex) >= degree) { // deg_G >= deg_Q
                        
                        bool pass_nlf = true;
                        for(auto iter = nlf[i].begin(); iter != nlf[i].end(); ++iter)
                        {   
                            if(data_graph->nlf[data_vertex].find(iter->first) == data_graph->nlf[data_vertex].end()) pass_nlf = false;
                            else
                            {
                                if(data_graph->nlf[data_vertex][iter->first] < iter->second) pass_nlf = false;
                            }
                        }
                        if(pass_nlf) query_graph->get_cands()[i].candidate.push_back(data_vertex);
                    }
                }
            }
            else {
                if(i < query_graph->get_parent_cands().size()) { // while parent pattern doesn't have this vertex, go to else.
                    for(auto it = query_graph->get_parent_cands()[i].candidate.begin(); 
                            it != query_graph->get_parent_cands()[i].candidate.end(); ++it) {

                        if(data_graph->get_vertex_degree(*it) >= degree) {

                            bool pass_nlf = true;
                            for(auto iter = nlf[i].begin(); iter != nlf[i].end(); ++iter)
                            {
                                if(data_graph->nlf[*it].find(iter->first) == data_graph->nlf[*it].end()) pass_nlf = false;
                                else
                                {
                                    if(data_graph->nlf[*it][iter->first] < iter->second) pass_nlf = false;
                                }
                            }
                            
                            if(pass_nlf)
                            {
                                if(query_graph->non_candidates[i].find(*it) == query_graph->non_candidates[i].end()) { // push-down pruning
                                    query_graph->get_cands()[i].candidate.push_back(*it);
                                }
                            }
                        }
                    }
                }
                else {
                    vector<VertexID> &data_vertices = data_graph->nodesByLabel[vlabel]; 

                    for(ui j = 0; j < data_vertices.size(); ++j) {
                        VertexID data_vertex = data_vertices[j];
                        if (data_graph->get_vertex_degree(data_vertex) >= degree) { // deg_G >= deg_Q
                            
                            bool pass_nlf = true;
                            for(auto iter = nlf[i].begin(); iter != nlf[i].end(); ++iter)
                            {   
                                if(data_graph->nlf[data_vertex].find(iter->first) == data_graph->nlf[data_vertex].end()) pass_nlf = false;
                                else
                                {
                                    if(data_graph->nlf[data_vertex][iter->first] < iter->second) pass_nlf = false;
                                }
                            }
                            
                            if(pass_nlf)
                            {
                                if(query_graph->non_candidates[i].find(data_vertex) == query_graph->non_candidates[i].end()) { // push-down pruning
                                    query_graph->get_cands()[i].candidate.push_back(data_vertex);
                                }
                            }
                        }
                    }
                }
            }
            if(query_graph->get_cands()[i].size() == 0) {

                delete[] nlf; 
                return false;
            }
        }
    } 
    else {
        for(ui i = 0; i < query_graph->size(); ++i) {
            vLabel vlabel = query_graph->get_p_vertex(i)->label;
            int degree = query_graph->get_vertex_degree(i);
            
            vector<VertexID> &data_vertices = data_graph->nodesByLabel[vlabel]; 

            for(ui j = 0; j < data_vertices.size(); ++j) {
                VertexID data_vertex = data_vertices[j];
                if (data_graph->get_vertex_degree(data_vertex) >= degree) { // deg_G >= deg_Q
                    query_graph->get_cands()[i].candidate.push_back(data_vertex);
                }
            }
            if(query_graph->get_cands()[i].size() == 0) {
                return false;
            }
        }
    }
    delete[] nlf; 

    return true;
}

bool GMatchEngine::DPisoFilter(bool automorphism, ui support)
{   

    if (!LDF_pruning(automorphism))
        return false;

    generateDPisoFilterPlan();

    ui query_vertices_num = query_graph->size();
    ui* updated_flag = new ui[data_graph->size()];
    ui* flag = new ui[data_graph->size()];
    std::fill(flag, flag + data_graph->size(), 0);


    // The number of refinement is k. According to the original paper, we set k as 3.
    // for(ui k = 0; k < 3; ++k) {
    //     if(k % 2 == 0) {
    //         for(int i = 1; i < query_vertices_num; ++i) {
    //             VertexID query_vertex = bfs_order[i];
    //             TreeNode& node = tree[query_vertex];
    //             pruneCandidates(query_vertex, node.bn_, node.bn_count_, flag, updated_flag);
    //         }
    //     }
    //     else {
    //         for(int i = query_vertices_num - 2; i >= 0; --i) {
    //             VertexID query_vertex = bfs_order[i];
    //             TreeNode& node = tree[query_vertex];
    //             pruneCandidates(query_vertex, node.fn_, node.fn_count_, flag, updated_flag);
    //         }
    //     }
    // }

    bool shrunk = true;
    while(shrunk) {
        for(int i = 1; i < query_vertices_num; ++i) {
            VertexID query_vertex = bfs_order[i];
            TreeNode& node = tree[query_vertex];
            shrunk = pruneCandidates(query_vertex, node.bn_, node.bn_count_, flag, updated_flag);
        }
    
        for(int i = query_vertices_num - 2; i >= 0; --i) {
            VertexID query_vertex = bfs_order[i];
            TreeNode& node = tree[query_vertex];
            shrunk = shrunk | pruneCandidates(query_vertex, node.fn_, node.fn_count_, flag, updated_flag);
        }
    }

    clear_filter();

    compactCandidates();

    delete[] updated_flag;
    delete[] flag;

    return isCandidateSetValid(support);
}

bool GMatchEngine::filterToConsistency(ui support)
{
    if (!LDF_pruning(false))
        return false;

    generateDPisoFilterPlan();

    ui query_vertices_num = query_graph->size();
    ui* updated_flag = new ui[data_graph->size()];
    ui* flag = new ui[data_graph->size()];
    std::fill(flag, flag + data_graph->size(), 0);

    bool shrunk = true;
    while(shrunk) {
        for(int i = 1; i < query_vertices_num; ++i) {
            VertexID query_vertex = bfs_order[i];
            TreeNode& node = tree[query_vertex];
            shrunk = pruneCandidates(query_vertex, node.bn_, node.bn_count_, flag, updated_flag);
        }
    
        for(int i = query_vertices_num - 2; i >= 0; --i) {
            VertexID query_vertex = bfs_order[i];
            TreeNode& node = tree[query_vertex];
            shrunk = shrunk | pruneCandidates(query_vertex, node.fn_, node.fn_count_, flag, updated_flag);
        }
    }

    clear_filter();
    
    compactCandidates();

    delete[] updated_flag;
    delete[] flag;

    return isCandidateSetValid(support);
}

void GMatchEngine::generateDPisoFilterPlan()
{
    VertexID start_vertex = selectDPisoStartVertex();
    bfsTraversal(start_vertex);

    ui query_vertices_num = query_graph->size();
    vector<int> order_index(query_vertices_num);
    for(ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = bfs_order[i];
        order_index[query_vertex] = i;
    }

    for(ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = bfs_order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count = query_graph->get_vertex_degree(u);
        for(ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = query_graph->get_p_vertex(u)->edges[j].to;
            if (order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
            }
            else {
                tree[u].fn_[tree[u].fn_count_++] = u_nbr;
            }
        }
    }
}

VertexID GMatchEngine::selectDPisoStartVertex() 
{
    double min_score = data_graph->size();
    VertexID start_vertex = 0;

    for(ui i = 0; i < query_graph->size(); ++i) {
        ui degree = query_graph->get_vertex_degree(i);
        if (degree <= 1)
            continue;

        double cur_score = query_graph->get_cands()[i].size() / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }
    return start_vertex;
}

void GMatchEngine::bfsTraversal(VertexID root_vertex)
{
    ui vertex_num = query_graph->size();

    queue<VertexID> bfs_queue;
    vector<bool> visited(vertex_num, false);

    tree = new TreeNode[vertex_num];
    for(ui i = 0; i < vertex_num; ++i) {
        tree[i].initialize(vertex_num);
    }
    bfs_order = new VertexID[vertex_num];

    ui visited_vertex_count = 0;
    bfs_queue.push(root_vertex);
    visited[root_vertex] = true;
    tree[root_vertex].level_ = 0;
    tree[root_vertex].id_ = root_vertex;

    while(!bfs_queue.empty()) {
        VertexID u = bfs_queue.front();
        bfs_queue.pop();
        bfs_order[visited_vertex_count++] = u;

        ui u_nbrs_count = query_graph->get_vertex_degree(u);
        for(ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = query_graph->get_p_vertex(u)->edges[i].to;

            if (!visited[u_nbr]) {
                bfs_queue.push(u_nbr);
                visited[u_nbr] = true;
                tree[u_nbr].id_ = u_nbr;
                tree[u_nbr].parent_ = u;
                tree[u_nbr].level_ = tree[u].level_ + 1; // parent.level+1
                tree[u].children_[tree[u].children_count_++] = u_nbr;
            }
        }
    }
}

bool GMatchEngine::pruneCandidates(VertexID query_vertex, VertexID *pivot_vertices, ui pivot_vertices_count, 
                                    ui *flag, ui *updated_flag)
{   
    bool shrunk = false;

    vLabel query_vertex_label = query_graph->get_p_vertex(query_vertex)->label;
    ui query_vertex_degree = query_graph->get_vertex_degree(query_vertex);

    ui count = 0;
    ui updated_flag_count = 0;
    for(ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];

        eLabel edge_label = query_graph->get_edge_label(query_vertex, pivot_vertex);

        for(ui j = 0; j < query_graph->get_cands()[pivot_vertex].size(); ++j) {
            VertexID v = query_graph->get_cands()[pivot_vertex].candidate[j];

            if(v == INVALID_VERTEX_ID)
                continue;

            for(ui k = 0; k < data_graph->get_vertex_degree(v); ++k) {
                VertexID v_nbr = data_graph->get_p_vertex(v)->edges[k].to;
                vLabel v_nbr_label = data_graph->get_p_vertex(v_nbr)->label;
                ui v_nbr_degree = data_graph->get_vertex_degree(v_nbr);

                eLabel v_nbr_edge_label = data_graph->get_p_vertex(v)->edges[k].label;

                if (flag[v_nbr] == count && edge_label == v_nbr_edge_label  // new constraint of edge label
                        && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
                    flag[v_nbr] += 1;

                    if (count == 0) {
                        updated_flag[updated_flag_count++] = v_nbr;
                    }
                }
            }
        }
        count += 1;
    }

    for(ui i = 0; i < query_graph->get_cands()[query_vertex].size(); ++i) {
        VertexID v = query_graph->get_cands()[query_vertex].candidate[i];
        if (v == INVALID_VERTEX_ID)
            continue;

        if (flag[v] != count) {
            query_graph->get_cands()[query_vertex].candidate[i] = INVALID_VERTEX_ID;
            shrunk = true;
        }
    }

    for(ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
    return shrunk;
}

void GMatchEngine::compactCandidates() // TODO erase?
{
    for(ui i = 0; i < query_graph->size(); ++i) {
        VertexID query_vertex = i;

        for(auto it = query_graph->get_cands()[query_vertex].candidate.begin(); 
                    it != query_graph->get_cands()[query_vertex].candidate.end(); ) {
            if (*it == INVALID_VERTEX_ID) {
                it = query_graph->get_cands()[query_vertex].candidate.erase(it);  
            } else {
                ++it;
            }
        }
    }
}

bool GMatchEngine::isCandidateSetValid() {
    for(ui i = 0; i < query_graph->size(); ++i) {
        if(query_graph->get_cands()[i].size() == 0)
            return false;
    }
    return true;
}

bool GMatchEngine::isCandidateSetValid(ui support) {
    for(ui i = 0; i < query_graph->size(); ++i) {
        if(query_graph->get_cands()[i].size() < support)
            return false;
    }
    return true;
}

void GMatchEngine::generateGQLQueryPlan(int start_vertex, ui *&matching_order_)
{
    /**
     * Select the vertex v such that (1) v is adjacent to the selected vertices; and (2) v has the minimum number of candidates.
     */
    vector<bool> visited_vertices(query_graph->size(), false);
    vector<bool> adjacent_vertices(query_graph->size(), false);
    matching_order_ = new ui[query_graph->size()];

    if(start_vertex == -1)
        start_vertex = selectGQLStartVertex();

    matching_order_[0] = start_vertex;
    updateValidVertices(start_vertex, visited_vertices, adjacent_vertices);

    for(ui i = 1; i < query_graph->size(); ++i) {
        VertexID next_vertex;
        ui min_value = data_graph->size() + 1;
        for(ui j = 0; j < query_graph->size(); ++j) {
            VertexID cur_vertex = j;

            if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                if (query_graph->get_cands()[cur_vertex].size() < min_value) {
                    min_value = query_graph->get_cands()[cur_vertex].size();
                    next_vertex = cur_vertex;
                }
                else if (query_graph->get_cands()[cur_vertex].size() == min_value && 
                        query_graph->get_vertex_degree(cur_vertex) > query_graph->get_vertex_degree(next_vertex)) {
                    next_vertex = cur_vertex;
                }
            }
        }
        updateValidVertices(next_vertex, visited_vertices, adjacent_vertices);
        matching_order_[i] = next_vertex;
    }
}

void GMatchEngine::generateGQLQueryPlan(int start_vertex)
{
    /**
     * Select the vertex v such that (1) v is adjacent to the selected vertices; and (2) v has the minimum number of candidates.
     */
    vector<bool> visited_vertices(query_graph->size(), false);
    vector<bool> adjacent_vertices(query_graph->size(), false);
    matching_order = new ui[query_graph->size()];

    if(start_vertex == -1)
        start_vertex = selectGQLStartVertex();

    matching_order[0] = start_vertex;
    updateValidVertices(start_vertex, visited_vertices, adjacent_vertices);

    for(ui i = 1; i < query_graph->size(); ++i) {
        VertexID next_vertex;
        ui min_value = data_graph->size() + 1;
        for(ui j = 0; j < query_graph->size(); ++j) {
            VertexID cur_vertex = j;

            if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                if (query_graph->get_cands()[cur_vertex].size() < min_value) {
                    min_value = query_graph->get_cands()[cur_vertex].size();
                    next_vertex = cur_vertex;
                }
                else if (query_graph->get_cands()[cur_vertex].size() == min_value && 
                        query_graph->get_vertex_degree(cur_vertex) > query_graph->get_vertex_degree(next_vertex)) {
                    next_vertex = cur_vertex;
                }
            }
        }
        updateValidVertices(next_vertex, visited_vertices, adjacent_vertices);
        matching_order[i] = next_vertex;
    }
}

VertexID GMatchEngine::selectGQLStartVertex()
{
    /**
     * Select the vertex with the minimum number of candidates as the start vertex.
     * Tie Handling:
     *  1. degree
     *  2. label id
     */

    VertexID start_vertex = 0;

    for(ui i = 1; i < query_graph->size(); ++i) {
        VertexID cur_vertex = i;

        if (query_graph->get_cands()[cur_vertex].size() < query_graph->get_cands()[start_vertex].size()) {
            start_vertex = cur_vertex;
        }
        else if (query_graph->get_cands()[cur_vertex].size() == query_graph->get_cands()[start_vertex].size()
                && query_graph->get_vertex_degree(cur_vertex) > query_graph->get_vertex_degree(start_vertex)) {
            start_vertex = cur_vertex;
        }
    }
    return start_vertex;
}

void GMatchEngine::updateValidVertices(VertexID query_vertex, vector<bool> &visited,
                                        vector<bool> &adjacent) 
{   
    visited[query_vertex] = true;
    ui nbr_cnt = query_graph->get_vertex_degree(query_vertex);

    for(ui i = 0; i < nbr_cnt; ++i) {
        ui nbr = query_graph->get_p_vertex(query_vertex)->edges[i].to;
        adjacent[nbr] = true;
    }
}

void GMatchEngine::generateBN(ui *&bn_count_, ui **&bn_, ui *&matching_order_) 
{
    ui query_vertices_num = query_graph->size();
    bn_count_ = new ui[query_vertices_num];
    std::fill(bn_count_, bn_count_ + query_vertices_num, 0);
    bn_ = new ui *[query_vertices_num];
    for(ui i = 0; i < query_vertices_num; ++i) {
        bn_[i] = new ui[query_vertices_num];
    }

    vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[matching_order_[0]] = true;
    for(ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = matching_order_[i];

        ui nbrs_cnt = query_graph->get_vertex_degree(vertex);
        for(ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = query_graph->get_p_vertex(vertex)->edges[j].to;

            if(visited_vertices[nbr]) {
                bn_[i][bn_count_[i]++] = nbr;
            }
        }
        visited_vertices[vertex] = true;
    }
}

void GMatchEngine::generateBN() 
{
    ui query_vertices_num = query_graph->size();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for(ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
    }

    vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[matching_order[0]] = true;
    for(ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = matching_order[i];

        ui nbrs_cnt = query_graph->get_vertex_degree(vertex);
        for(ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = query_graph->get_p_vertex(vertex)->edges[j].to;

            if(visited_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }
        visited_vertices[vertex] = true;
    }
}

void GMatchEngine::buildTable()
{   
    ui query_vertices_num = query_graph->size();
    ui* flag = new ui[data_graph->size()];
    ui* updated_flag = new ui[data_graph->size()];
    std::fill(flag, flag + data_graph->size(), 0);

    edge_matrix = new Edges **[query_graph->size()];
    for (ui i = 0; i < query_graph->size(); ++i) {
        edge_matrix[i] = new Edges *[query_graph->size()];
    }

    for(ui i = 0; i < query_vertices_num; ++i) {
        for(ui j = 0; j < query_vertices_num; ++j) {
            edge_matrix[i][j] = NULL;
        }
    }

    vector<VertexID> build_table_order(query_vertices_num);
    for(ui i = 0; i < query_vertices_num; ++i) {
        build_table_order[i] = i;
    }


    std::sort(build_table_order.begin(), build_table_order.end(), [=](VertexID l, VertexID r) {
        if (query_graph->get_vertex_degree(l) == query_graph->get_vertex_degree(r)) {
            return l < r;
        }
        return query_graph->get_vertex_degree(l) > query_graph->get_vertex_degree(r);
    });

    vector<ui> temp_edges(data_graph->get_nedges() * 2);

    for(auto u : build_table_order) {
        ui u_nbrs_count = query_graph->get_vertex_degree(u);

        ui updated_flag_count = 0;

        for(ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = query_graph->get_p_vertex(u)->edges[i].to;

            if (edge_matrix[u][u_nbr] != NULL)
                continue;

            eLabel edge_label = query_graph->get_edge_label(u_nbr, u);

            if (updated_flag_count == 0) { // record u's CS in the first round

                for(ui j = 0; j < query_graph->get_cands()[u].size(); ++j) {
                    VertexID v = query_graph->get_cands()[u].candidate[j];
                    flag[v] = j + 1; // v's ID -> index in u's CS
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_nbr][u] = new Edges;
            edge_matrix[u_nbr][u]->vertex_count_ = query_graph->get_cands()[u_nbr].size();
            edge_matrix[u_nbr][u]->offset_ = new ui[query_graph->get_cands()[u_nbr].size() + 1];

            edge_matrix[u][u_nbr] = new Edges;
            edge_matrix[u][u_nbr]->vertex_count_ = query_graph->get_cands()[u].size();
            edge_matrix[u][u_nbr]->offset_ = new ui[query_graph->get_cands()[u].size() + 1];
            std::fill(edge_matrix[u][u_nbr]->offset_, edge_matrix[u][u_nbr]->offset_ + query_graph->get_cands()[u].size() + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            for(ui j = 0; j < query_graph->get_cands()[u_nbr].size(); ++j) {
                VertexID v = query_graph->get_cands()[u_nbr].candidate[j];
                edge_matrix[u_nbr][u]->offset_[j] = local_edge_count;

                ui v_nbrs_count = data_graph->get_vertex_degree(v);

                ui local_degree = 0;

                for(ui k = 0; k < v_nbrs_count; ++k) {
                    VertexID v_nbr = data_graph->get_p_vertex(v)->edges[k].to; // u's neighors's candidate's neighbor (in G)

                    eLabel v_nbr_edge_label = data_graph->get_edge_label(v, v_nbr);
                    
                    // v_nbr is in u's CS && constraints of edge label
                    if(flag[v_nbr] != 0 && v_nbr_edge_label == edge_label) {
                        ui position = flag[v_nbr] - 1; // get index in u's CS
                        temp_edges[local_edge_count++] = position;
                        edge_matrix[u][u_nbr]->offset_[position + 1] += 1; // if u matches to position(th), u_nbr can match to j(th)
                        local_degree += 1;
                    }
                }

                if(local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_nbr][u]->offset_[query_graph->get_cands()[u_nbr].size()] = local_edge_count;
            edge_matrix[u_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_nbr][u]->edge_ = new ui[local_edge_count];
            std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_nbr][u]->edge_);

            edge_matrix[u][u_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_nbr]->edge_ = new ui[local_edge_count];

            local_max_degree = 0;
            for(ui j = 1; j <= query_graph->get_cands()[u].size(); ++j) {
                if (edge_matrix[u][u_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_nbr]->offset_[j];
                }
                edge_matrix[u][u_nbr]->offset_[j] += edge_matrix[u][u_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_nbr]->max_degree_ = local_max_degree;

            for(ui j = 0; j < query_graph->get_cands()[u_nbr].size(); ++j) {
                ui begin = j; // index in u_nbr
                for(ui k = edge_matrix[u_nbr][u]->offset_[begin]; k < edge_matrix[u_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix[u_nbr][u]->edge_[k]; // index in u

                    edge_matrix[u][u_nbr]->edge_[edge_matrix[u][u_nbr]->offset_[end]++] = begin;
                }
            }

            for(ui j = query_graph->get_cands()[u].size(); j >= 1; --j) { // reset edge_matrix[u][u_nbr]->offset_
                edge_matrix[u][u_nbr]->offset_[j] = edge_matrix[u][u_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_nbr]->offset_[0] = 0;
        }

        for(ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }

    delete[] flag;
    delete[] updated_flag;
}

void GMatchEngine::buildTable(Edges ***&edge_matrix_)
{       
    ui query_vertices_num = query_graph->size();
    ui* flag = new ui[data_graph->size()];
    ui* updated_flag = new ui[data_graph->size()];
    std::fill(flag, flag + data_graph->size(), 0);


    vector<VertexID> build_table_order(query_vertices_num);
    for(ui i = 0; i < query_vertices_num; ++i) {
        build_table_order[i] = i;
    }




    std::sort(build_table_order.begin(), build_table_order.end(), [=](VertexID l, VertexID r) {
        if (query_graph->get_vertex_degree(l) == query_graph->get_vertex_degree(r)) {
            return l < r;
        }
        return query_graph->get_vertex_degree(l) > query_graph->get_vertex_degree(r);
    });

    vector<ui> temp_edges(data_graph->get_nedges() * 2);

            
    for(auto u : build_table_order) {
        ui u_nbrs_count = query_graph->get_vertex_degree(u);

        ui updated_flag_count = 0;

        for(ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = query_graph->get_p_vertex(u)->edges[i].to;

            if (edge_matrix_[u][u_nbr] != NULL)
                continue;

            eLabel edge_label = query_graph->get_edge_label(u_nbr, u);

            if (updated_flag_count == 0) { // record u's CS in the first round

                for(ui j = 0; j < query_graph->get_cands()[u].size(); ++j) {
                    VertexID v = query_graph->get_cands()[u].candidate[j];
                    flag[v] = j + 1; // v's ID -> index in u's CS
                    updated_flag[updated_flag_count++] = v;
                }
            }


            edge_matrix_[u_nbr][u] = new Edges;
            edge_matrix_[u_nbr][u]->vertex_count_ = query_graph->get_cands()[u_nbr].size();
            edge_matrix_[u_nbr][u]->offset_ = new ui[query_graph->get_cands()[u_nbr].size() + 1];

            edge_matrix_[u][u_nbr] = new Edges;
            edge_matrix_[u][u_nbr]->vertex_count_ = query_graph->get_cands()[u].size();
            edge_matrix_[u][u_nbr]->offset_ = new ui[query_graph->get_cands()[u].size() + 1];
            std::fill(edge_matrix_[u][u_nbr]->offset_, edge_matrix_[u][u_nbr]->offset_ + query_graph->get_cands()[u].size() + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            for(ui j = 0; j < query_graph->get_cands()[u_nbr].size(); ++j) {
                VertexID v = query_graph->get_cands()[u_nbr].candidate[j];
                edge_matrix_[u_nbr][u]->offset_[j] = local_edge_count;

                ui v_nbrs_count = data_graph->get_vertex_degree(v);

                ui local_degree = 0;

                for(ui k = 0; k < v_nbrs_count; ++k) {
                    VertexID v_nbr = data_graph->get_p_vertex(v)->edges[k].to; // u's neighors's candidate's neighbor (in G)

                    eLabel v_nbr_edge_label = data_graph->get_edge_label(v, v_nbr);
                    
                    // v_nbr is in u's CS && constraints of edge label
                    if(flag[v_nbr] != 0 && v_nbr_edge_label == edge_label) {
                        ui position = flag[v_nbr] - 1; // get index in u's CS
                        temp_edges[local_edge_count++] = position;
                        edge_matrix_[u][u_nbr]->offset_[position + 1] += 1; // if u matches to position(th), u_nbr can match to j(th)
                        local_degree += 1;
                    }
                }

                if(local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }


            edge_matrix_[u_nbr][u]->offset_[query_graph->get_cands()[u_nbr].size()] = local_edge_count;
            edge_matrix_[u_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix_[u_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix_[u_nbr][u]->edge_ = new ui[local_edge_count];
            std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix_[u_nbr][u]->edge_);

            edge_matrix_[u][u_nbr]->edge_count_ = local_edge_count;
            edge_matrix_[u][u_nbr]->edge_ = new ui[local_edge_count];

            

            local_max_degree = 0;
            for(ui j = 1; j <= query_graph->get_cands()[u].size(); ++j) {
                if (edge_matrix_[u][u_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix_[u][u_nbr]->offset_[j];
                }
                edge_matrix_[u][u_nbr]->offset_[j] += edge_matrix_[u][u_nbr]->offset_[j - 1];
            }

            edge_matrix_[u][u_nbr]->max_degree_ = local_max_degree;


            for(ui j = 0; j < query_graph->get_cands()[u_nbr].size(); ++j) {
                ui begin = j; // index in u_nbr
                for(ui k = edge_matrix_[u_nbr][u]->offset_[begin]; k < edge_matrix_[u_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix_[u_nbr][u]->edge_[k]; // index in u

                    edge_matrix_[u][u_nbr]->edge_[edge_matrix_[u][u_nbr]->offset_[end]++] = begin;
                }
            }

            for(ui j = query_graph->get_cands()[u].size(); j >= 1; --j) { // reset edge_matrix[u][u_nbr]->offset_
                edge_matrix_[u][u_nbr]->offset_[j] = edge_matrix_[u][u_nbr]->offset_[j - 1];
            }
            edge_matrix_[u][u_nbr]->offset_[0] = 0;
        }


        for(ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }

    delete[] flag;
    delete[] updated_flag;
}
/** 
 * backtrack until finding one result, without timeout
 */
int GMatchEngine::backtrack(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                                unordered_set<VertexID> &visited, vector<ui> &counter)
{
    if(depth == query_graph->size()) {
        return -1;
    }
    else {
        VertexID u = matching_order[depth];
        vector<ui> valid_candidate_idx;
        ui valid_candidate_cnt;
        if(depth == 0) {
            valid_candidate_cnt = query_graph->get_cands()[u].size();
            for(ui i = 0; i < valid_candidate_cnt; ++i) {
                valid_candidate_idx.push_back(i);
            }
        } else {
            generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                        bn[depth], bn_count[depth]);
        }
        for(ui i = 0; i < valid_candidate_cnt; ++i) {

            counter[depth] = i; // used as stack counter

            ui valid_idx = valid_candidate_idx[i];
            VertexID v = query_graph->get_cands()[u].candidate[valid_idx];
            if(visited.find(v) != visited.end()) {
                continue;
            }
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited.insert(v);

            int value = backtrack(depth+1, embedding, idx_embedding, visited, counter);

            if(value == -1) return -1;
            visited.erase(v);
        }
    }
    return -2;
}


double GMatchEngine::countElaspedTime()
{
    struct timeb cur_time;
    ftime(&cur_time);
    return cur_time.time-gtime_start.time+(double)(cur_time.millitm-gtime_start.millitm)/1000;
}


#ifdef OPTIMIZED_MATCH

void GMatchEngine::resetVisitedArray(ui depth, vector<VertexID> &embedding)
{
    for (ui i = 0; i <= depth; ++i)
    {
        mvisited_array[embedding[matching_order[i]]] = false;
    }
}

/**
 * iterative-style backtrack until finding one result, with timeout and recording status,
 * and only execute once.
 */
int GMatchEngine::LFTJ(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                                unordered_set<VertexID> &visited, vector<ui> &counter, MatchingStatus &ms)
{
    mvisited_array[embedding[matching_order[0]]] = true;

    ui cur_depth = depth;
    ui max_depth = query_graph->size();

    if (cur_depth == 0)
    {
        VertexID start_vertex = matching_order[cur_depth];
        counter[cur_depth] = 0; //#### TODO:
        midx_count[cur_depth] = query_graph->get_cands()[start_vertex].size();

        for (ui i = 0; i < midx_count[cur_depth]; ++i) 
        {
            mvalid_candidate_idx[cur_depth][i] = i;
        }
    }
    else
    {
        counter[cur_depth] = 0; //#### TODO: 
        generateValidCandidateIndex(cur_depth, idx_embedding);
    }

    while (true)
    {
        while (counter[cur_depth] < midx_count[cur_depth])
        {
            ui valid_idx = mvalid_candidate_idx[cur_depth][counter[cur_depth]];
            ui u = matching_order[cur_depth];
            ui v = query_graph->get_cands()[u].candidate[valid_idx];

            // if(visited.find(v) != visited.end()) 
            if (mvisited_array[v])
            {
                counter[cur_depth] += 1;
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;

            mvisited_array[v] = true;

            counter[cur_depth] += 1;
            if (cur_depth == max_depth - 1) 
            { 
                resetVisitedArray(cur_depth, embedding); 
                return -1;
            }

            if (countElaspedTime() < TIMEOUT_THRESHOLD)
            {
                cur_depth += 1;
                counter[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding);
            }
            else
            {
                if(!ms.init)
                {
                    ms.init = true;
                    ms.ffdepth = cur_depth+1;
                    ms.counter = counter;
                }
                mvisited_array[v] = false;
                resetVisitedArray(cur_depth, embedding); 
                return -3;
            }
        }
        cur_depth -= 1;
        if (cur_depth < depth)
            break;
        else
        {
            mvisited_array[embedding[matching_order[cur_depth]]] = false;
        }
    }

    mvisited_array[embedding[matching_order[0]]] = false;
    return -2;
}
#endif


/** 
 * backtrack until finding one result, with timeout and recording status
 */
int GMatchEngine::backtrack(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                                unordered_set<VertexID> &visited, vector<ui> &counter, MatchingStatus &ms)
{   
    if(depth == query_graph->size()) {
        return -1;
    }

    struct timeb cur_time;
    double drun_time;

    ftime(&cur_time);
    drun_time = cur_time.time-gtime_start.time+(double)(cur_time.millitm-gtime_start.millitm)/1000;
    if(drun_time >= TIMEOUT_THRESHOLD) {
        if(!ms.init) {
            ms.init = true;
            ms.ffdepth = depth;
            ms.counter = counter;
        }
        return -3;
    }
    else {
        VertexID u = matching_order[depth];
        vector<ui> valid_candidate_idx;
        ui valid_candidate_cnt;
        if(depth == 0) {
            valid_candidate_cnt = query_graph->get_cands()[u].size();
            for(ui i = 0; i < valid_candidate_cnt; ++i) {
                valid_candidate_idx.push_back(i);
            }
        } else {
            generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                        bn[depth], bn_count[depth]);
        }
        for(ui i = 0; i < valid_candidate_cnt; ++i) {

            counter[depth] = i; // used as stack counter

            int valid_idx = valid_candidate_idx[i];
            int v = query_graph->get_cands()[u].candidate[valid_idx];
            if(visited.find(v) != visited.end()) {
                continue;
            }
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited.insert(v);

            int value = backtrack(depth+1, embedding, idx_embedding, visited, counter, ms);

            if(value == -1) return -1;
            else if(value == -3) return -3;

            visited.erase(v);
        }
        counter[depth] = 0;
    }
    return -2;
}

int GMatchEngine::backtrackResume(ui depth, vector<VertexID> &embedding, vector<ui> &idx_embedding,   
                                unordered_set<VertexID> &visited, vector<ui> &counter, MatchingStatus &ms, bool &skip)
{   
    if(depth == query_graph->size()) {
        return -1;
    }
    else {
        VertexID u = matching_order[depth];
        vector<ui> valid_candidate_idx;
        ui valid_candidate_cnt;
        if(depth == 0) {
            valid_candidate_cnt = query_graph->get_cands()[u].size();
            for(ui i = 0; i < valid_candidate_cnt; ++i) {
                valid_candidate_idx.push_back(i);
            }
        } else {
            generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                        bn[depth], bn_count[depth]);
        }
        ui i = skip? counter[depth]: 0;
        for(; i < valid_candidate_cnt; ++i) {

            if(skip && depth == ms.depth) skip = false;

            int valid_idx = valid_candidate_idx[i];
            int v = query_graph->get_cands()[u].candidate[valid_idx];
            if(visited.find(v) != visited.end()) {
                continue;
            }
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited.insert(v);

            int value = backtrackResume(depth+1, embedding, idx_embedding, visited, counter, ms, skip);

            if(value == -1) return -1;

            visited.erase(v);
        }
    }
    return -2;
}

/**
 * If automorphisms, find all matching results.
 */
void GMatchEngine::backtrackAll(ui depth, VertexID *embedding, ui *idx_embedding,   
                                unordered_set<VertexID> &visited, VtxSetVec &results)
{
    if(depth == query_graph->size()) {
        for(ui i = 0; i < query_graph->size(); ++i) {
            results[i].insert(embedding[i]);
        }
    } else {
        VertexID u = matching_order[depth];
        vector<ui> valid_candidate_idx;
        ui valid_candidate_cnt;
        if(depth == 0) { // # BN = 0
            valid_candidate_cnt = query_graph->get_cands()[u].size();
            for(ui i = 0; i < valid_candidate_cnt; ++i) {
                valid_candidate_idx.push_back(i);
            }
        } else {
            generateValidCandidateIndex(u, idx_embedding, valid_candidate_cnt, valid_candidate_idx,
                                        bn[depth], bn_count[depth]);
        }
        for(ui i = 0; i < valid_candidate_cnt; ++i) {
            int valid_idx = valid_candidate_idx[i];
            int v = query_graph->get_cands()[u].candidate[valid_idx];
            if(visited.find(v) != visited.end()) {
                continue;
            }
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited.insert(v);

            backtrackAll(depth+1, embedding, idx_embedding, visited, results);

            visited.erase(v);
        }
    }
}

void GMatchEngine::generateValidCandidateIndex(VertexID u, ui *idx_embedding, ui &valid_candidate_cnt, 
                                                vector<ui> &valid_candidate_index, VertexID *bn, ui bn_cnt)
{
    if(bn_cnt == 1) {
        VertexID current_bn = bn[0];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];

        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];
        valid_candidate_cnt = current_candidates_count;
        for(ui i = 0; i < current_candidates_count; ++i) {
            valid_candidate_index.push_back(current_candidates[i]);
        }

    } else {
        vector<vector<ui> > vecs;
        for(ui i = 0; i < bn_cnt; ++i) {
            VertexID current_bn = bn[i];

            Edges &current_edge = *edge_matrix[current_bn][u];
            ui current_index_id = idx_embedding[current_bn];

            ui current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
            ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

            vector<ui> vec(current_candidates, current_candidates+current_candidates_count);
            

            vecs.push_back(vec);
        }

        valid_candidate_index = leapfrogJoin(vecs);
        valid_candidate_cnt = valid_candidate_index.size();
    }
}

void GMatchEngine::generateValidCandidateIndex(VertexID u, vector<ui> &idx_embedding, ui &valid_candidate_cnt, 
                                                vector<ui> &valid_candidate_index, VertexID *bn, ui bn_cnt)
{
    if(bn_cnt == 1) {
        VertexID current_bn = bn[0];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];

        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];
        valid_candidate_cnt = current_candidates_count;
        for(ui i = 0; i < current_candidates_count; ++i) {
            valid_candidate_index.push_back(current_candidates[i]);
        }

    } else {
        vector<vector<ui> > vecs;
        for(ui i = 0; i < bn_cnt; ++i) {
            VertexID current_bn = bn[i];

            Edges &current_edge = *edge_matrix[current_bn][u];
            ui current_index_id = idx_embedding[current_bn];

            ui current_candidates_count =
                current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
            ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

            vector<ui> vec(current_candidates, current_candidates+current_candidates_count);
            
            vecs.push_back(vec);
        }

        valid_candidate_index = leapfrogJoin(vecs);
        valid_candidate_cnt = valid_candidate_index.size();
    }
}


#ifdef OPTIMIZED_MATCH

void GMatchEngine::generateValidCandidateIndex(ui depth, vector<ui> &idx_embedding)
{
    VertexID u = matching_order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;

    Edges& previous_edge = *edge_matrix[previous_bn][u];

    valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    memcpy(mvalid_candidate_idx[depth], previous_candidates, valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_count[depth]; ++i) {
        
        VertexID current_bn = bn[depth][i];

        Edges& current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];
        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        if (current_candidates_count < valid_candidates_count)
            ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, mvalid_candidate_idx[depth], valid_candidates_count,
                        mtemp_buffer, temp_count);
        else
            ComputeSetIntersection::ComputeCandidates(mvalid_candidate_idx[depth], valid_candidates_count, current_candidates, current_candidates_count,
                        mtemp_buffer, temp_count);
        
        for(int i = 0; i < temp_count; ++i)
        {
            mvalid_candidate_idx[depth][i] = mtemp_buffer[i];
        }
        valid_candidates_count = temp_count;
    }
    midx_count[depth] = valid_candidates_count;
}

#endif


/**
 * ================================= below are clear methods ================================
 */

void GMatchEngine::clear_all()
{   
    // clear_filter();
    clear_table();
    clear_match();
    reset();
}

void GMatchEngine::clear_match()
{
    for(ui i = 0; i < query_graph->size(); ++i) {
        delete[] bn[i];
    }
    delete[] bn;
    delete[] bn_count;
    delete[] matching_order;
}

void GMatchEngine::clear_filter()
{
    delete[] tree;
    delete[] bfs_order;
}

void GMatchEngine::clear_table()
{
    for (ui i = 0; i < query_graph->size(); ++i) {
        for (ui j = 0; j < query_graph->size(); ++j) {
            delete edge_matrix[i][j];
        }
        delete[] edge_matrix[i];
    }
    delete[] edge_matrix;
}

void GMatchEngine::reset()
{
    query_graph = NULL;
    data_graph = NULL;
}