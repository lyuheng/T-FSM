#pragma once

#include "graph.h"
#include "gmatch.h"
#include "decompose.h"
#include "pretty_print.h"
#include "canonical.h"
#include "setting.h"

#include <map>
#include <fstream>

using std::ofstream;

// typedef map<dfs_code_t, Pattern, cmp> PatternMap;

// typedef unordered_map<string, Pattern> PatternMap;

typedef unordered_map<dfs_code_t, Pattern, DfsCodeHashCode> PatternMap; // only used for initial_pattern_set
typedef vector<Pattern> PatternVec;

class GraMi
{
public:
    int nsupport_;
    Graph pruned_graph; // pruned graph after deleting infrequent vertice and edges
    // PatternMap init_pattern_map;
    PatternMap init_pattern_map;

    unordered_map<vLabel, unordered_set<vLabel>> neighbor_labels;

    GMatchEngine gmatch_engine; // subgraph matching & find automorphisms
    GMatchEngine exist_engine;  // check existence

    Decompose decomposer;

    unordered_map<string, VectorOfSet> cache; // TODO trie later

    ofstream fout;

    ui total_num;

    GraMi(int nsupport) : nsupport_(nsupport), pruned_graph(nsupport)
    {
        char file[20];
        strcpy(file, "./GraMi_results.txt");
        fout.open(file);
        total_num = 0;
    }

    ~GraMi()
    {
        fout.close();
    }

    void initialize();
    void extend(Pattern &pattern, PatternVec &ext_pattern_map);
    bool isCan(Pattern &pattern);
    int frequency(Pattern &pattern);

    void project();
    void mine_subgraph(Pattern &pattern); // recursive function

    void find_existence(Pattern &pattern, VectorOfSet &found_results);
    int get_domain_minsize(Pattern &pattern);

    VertexID has_been_computed(VectorOfSet &autos, vector<bool> &preComputed, VertexID index);

    int search_mapping(vector<VertexID> &mapping, VertexID u);

    void report(Pattern &pattern);

    void verify(Pattern &pattern);
};

void GraMi::initialize()
{
    for (ui i = 0; i < pruned_graph.size(); ++i)
    {
        const vertex_t *vertex = pruned_graph.get_p_vertex(i);
        vLabel firstlabel = vertex->label;
        int degree = vertex->edges.size();
        for (ui j = 0; j < degree; j++)
        {
            VertexID nbrid = vertex->edges[j].to;
            eLabel edgelabel = vertex->edges[j].label;
            vLabel nbrlabel = pruned_graph.get_vertex_label(nbrid);

            if(firstlabel <= nbrlabel) // Partial Pruning
            {
                dfs_code_t dfscode(0, 1, firstlabel, edgelabel, nbrlabel);
                Pattern pattern;
                pattern.init(dfscode, NULL);

                if (init_pattern_map.find(dfscode) == init_pattern_map.end())
                {
                    init_pattern_map[dfscode] = pattern;
                }
            }
        }
    }

    for (auto it = init_pattern_map.begin(); it != init_pattern_map.end(); )
    {
        auto &pattern = it->second;
        if (frequency(pattern) < nsupport_)
        {
            it = init_pattern_map.erase(it);
        }
        else
        {
            it++;
        }
    }

    for (auto &pa : init_pattern_map)
    {
        auto &pattern = pa.second;
        dfs_code_t &last_edge = pattern.dfscodes[0];
        vLabel from_label = last_edge.from_label;
        vLabel to_label = last_edge.to_label;

        neighbor_labels[from_label].insert(to_label);
        neighbor_labels[to_label].insert(from_label);
    }

    cout << "Initial Pattern Map size: " << init_pattern_map.size() << endl;
}

/** rightmost extension */
void GraMi::extend(Pattern &pattern, PatternVec &ext_pattern_vec)
{
    VertexID last_id = pattern.size() - 1;
    vertex_t *last_vertex = pattern.get_p_vertex(last_id);
    vLabel min_label = pattern.dfscodes[0].from_label;
    // edge_t &last_edge = pattern.dfscodes[pattern.right_most_path[0]];
    unordered_set<eLabel> &freq_edge_labels = pruned_graph.freqEdgeLabels;

    // find extension of last node, forward
    if(Settings::maxNumNodes == -1 || pattern.size() < Settings::maxNumNodes)
    {
        for (vLabel nbr_label : neighbor_labels[last_vertex->label])
        {
            for (eLabel edge_label : freq_edge_labels)
            {
                dfs_code_t dfs_code(last_id, last_id + 1, last_vertex->label, edge_label, nbr_label);

                if(nbr_label >= min_label)
                {
                    ext_pattern_vec.resize(ext_pattern_vec.size() + 1);
                    ext_pattern_vec.back().copy(pattern);
                    ext_pattern_vec.back().extend(dfs_code);
                }
            }
        }
    }
    
    // find extension of last node, backward
    for (ui i = pattern.right_most_path.size(); i > 1; --i)
    {
        VertexID cur_id = pattern.dfscodes[pattern.right_most_path[i - 1]].from;
        vertex_t *cur_vertex = pattern.get_p_vertex(cur_id);
        vLabel cur_nbr_label = pattern.dfscodes[pattern.right_most_path[i - 1]].to_label;
        eLabel cur_edge_label = pattern.dfscodes[pattern.right_most_path[i - 1]].edge_label;

        if (neighbor_labels[last_vertex->label].find(cur_vertex->label) != neighbor_labels[last_vertex->label].end())
        {
            bool is_neighbor = false;
            for (ui j = 0; j < last_vertex->edges.size(); ++j)
            {
                vertex_t *nbr_vertex = pattern.get_p_vertex(last_vertex->edges[j].to);
                if (cur_vertex->id == nbr_vertex->id)
                {
                    is_neighbor = true;
                    break;
                }
            }
            if (is_neighbor)
                continue;
            for (eLabel edge_label : freq_edge_labels)
            {
                dfs_code_t dfs_code(last_id, cur_vertex->id, last_vertex->label, edge_label, cur_vertex->label);
                if(edge_label > cur_edge_label || (edge_label == cur_edge_label && last_vertex->label >= cur_nbr_label))
                { 
                    ext_pattern_vec.resize(ext_pattern_vec.size() + 1);
                    ext_pattern_vec.back().copy(pattern);
                    ext_pattern_vec.back().extend(dfs_code);
                }
            }
        }
    }

    // find extension on rightmost path, forward
    if(Settings::maxNumNodes == -1 || pattern.size() < Settings::maxNumNodes)
    {
        for (ui i = 0; i < pattern.right_most_path.size(); ++i)
        {
            int idx = pattern.right_most_path[i];
            int cur_vertex_id = pattern.dfscodes[idx].from;
            vertex_t *cur_vertex = pattern.get_p_vertex(cur_vertex_id);
            
            vLabel cur_nbr_label = pattern.dfscodes[idx].to_label;
            eLabel cur_edge_label = pattern.dfscodes[idx].edge_label;

            for (vLabel nbr_label : neighbor_labels[cur_vertex->label])
            {
                for (eLabel edge_label : freq_edge_labels)
                {
                    dfs_code_t dfs_code(cur_vertex->id, last_id + 1, cur_vertex->label, edge_label, nbr_label);

                    if(cur_edge_label < edge_label || (cur_edge_label == edge_label && cur_nbr_label <= nbr_label))
                    {
                        ext_pattern_vec.resize(ext_pattern_vec.size() + 1);
                        ext_pattern_vec.back().copy(pattern);
                        ext_pattern_vec.back().extend(dfs_code);
                    }
                }
            }
        }
    }
}

void GraMi::project()
{   
    initialize();
    for (auto it = init_pattern_map.begin(); it != init_pattern_map.end(); ++it)
    {   
        mine_subgraph(it->second);
    }
    cout << "[RESULT] # Frequent Patterns: " << total_num << endl;

    delete[] pruned_graph.nlf;
}

void GraMi::mine_subgraph(Pattern &pattern)
{
    if (!isCan(pattern))
    {
        return; 
    }
    total_num++;
    report(pattern);
    // verify(pattern);
    PatternVec ext_pattern_vec;
    extend(pattern, ext_pattern_vec);

    for (auto it = ext_pattern_vec.begin(); it != ext_pattern_vec.end(); ++it)
    {
        int nsupport = frequency(*it);

        if (nsupport < nsupport_)
        {   
            continue;
        }
        mine_subgraph(*it);
    }
}


int GraMi::frequency(Pattern &pattern)
{
    if (pattern.size() == 2)
    {
        vLabel labelA = pattern.get_vertex_label(0);
        vLabel labelB = pattern.get_vertex_label(1);
        eLabel edgeLabel = pattern.get_p_vertex(0)->edges[0].label;

        PairVertices pv;

        if (labelA <= labelB)
            pv.set(labelA, edgeLabel, labelB);
        else
            pv.set(labelB, edgeLabel, labelA);

        if (pruned_graph.hashedEdges.find(pv) == pruned_graph.hashedEdges.end())
            return 0;
        else
        {
            EdgeCandidate & edge_cand = pruned_graph.hashedEdges[pv];
            return std::min(edge_cand.candA.size(), edge_cand.candB.size());
        }
    }

    // push-down pruning, decompose
    pattern.non_candidates.resize(pattern.size());

    vector<vector<subPattern> > maps;
    decomposer.decompose(pattern, maps);
    for (ui i = 0; i < maps.size(); ++i)
    {
        vector<subPattern> &edge_removed = maps[i];
        for (auto it = edge_removed.begin(); it != edge_removed.end(); ++it)
        {
            Pattern key = (*it).pattern;
            if (cache.find(key.toString()) != cache.end())
            {
                vector<VertexID> &subgraph_mappings = (*it).mapping;
                VectorOfSet &node_non_candidates = cache[key.toString()];
                for (ui j = 0; j < subgraph_mappings.size(); ++j)
                {
                    unordered_set<VertexID> &non_cans = pattern.non_candidates[subgraph_mappings[j]];
                    unordered_set<VertexID> &non_cans_check = node_non_candidates[j];
                    for (auto it2 = non_cans_check.begin(); it2 != non_cans_check.end(); ++it2)
                    {
                        non_cans.insert(*it2);
                    }
                }
            }
        }
    }

    VectorOfSet found_results(pattern.size());

    find_existence(pattern, found_results);

    // since lazy search, this frequence may not be accurate
    int freq = found_results[0].size();
    for (ui i = 1; i < pattern.size(); ++i)
    {
        if (freq > found_results[i].size())
        {
            freq = found_results[i].size();
        }
    }

    string code = pattern.toString();
    cache[code] = pattern.non_candidates; // deep copy

    pattern.non_candidates.clear();

    return freq;
}

int GraMi::get_domain_minsize(Pattern &pattern)
{
    int min = pattern.candidates[0].size();
    for (ui i = 0; i < pattern.size(); ++i)
    {
        if (min > pattern.candidates[i].size())
        {
            min = pattern.candidates[i].size();
        }
    }
    return min;
}

void GraMi::find_existence(Pattern &pattern, VectorOfSet &found_results)
{ 
    if (pattern.size() == 2 && pattern.get_vertex_label(0) != pattern.get_vertex_label(1))
    {   
        // special case
        // Do arc pruning until reaching consistency
        gmatch_engine.set(&pruned_graph, &pattern);
        gmatch_engine.filterToConsistency(nsupport_);
        // gmatch_engine.clear_filter();
        gmatch_engine.reset();

        // clone
        for (ui i = 0; i < pattern.size(); ++i)
        {
            found_results[i].insert(pattern.candidates[i].candidate.begin(), pattern.candidates[i].candidate.end());
        }
        return;
    }

    if (pattern.distinct_labels() && pattern.is_acyclic())
    {   
        // special case
        // Do arc pruning until reaching consistency
        gmatch_engine.set(&pruned_graph, &pattern);
        gmatch_engine.filterToConsistency(nsupport_);
        // gmatch_engine.clear_filter();
        gmatch_engine.reset();

        // clone
        for (ui i = 0; i < pattern.size(); ++i)
        {
            found_results[i].insert(pattern.candidates[i].candidate.begin(), pattern.candidates[i].candidate.end());
        }
        return;
    }

    // Find automorphisms
    vector<bool> pre_computed(pattern.size(), false);
    VectorOfSet autos(pattern.size());

    // convert a pattern to a Graph

    Graph dataGraphAuto;
    dataGraphAuto.toGraph(pattern);
    
    // find automorphisms
    gmatch_engine.set(&dataGraphAuto, &pattern);
    gmatch_engine.executeAuto(autos);
    gmatch_engine.clear_all();
    pattern.clear_candidate();

    // Now do the search
    gmatch_engine.set(&pruned_graph, &pattern);
    if(!gmatch_engine.DPisoFilter(false, nsupport_)) // degree-based pruning.
    {
        gmatch_engine.reset();
        return;
    }
    
    gmatch_engine.buildTable();

    vector<MatchingStatus> status_vector;

    // auto pre_non_candidates = pattern.non_candidates; // deep copy

    for (int i = pattern.size() - 1; i >= 0; --i)
    {
        bool search = true;
        VertexID pre_idx = has_been_computed(autos, pre_computed, i);
        if (i != pre_idx)
        {
            search = false;
            found_results[i] = found_results[pre_idx]; // deep copy
        }
        if (search)
        {
            ui idx = 0;
            gmatch_engine.generateGQLQueryPlan(i);
            gmatch_engine.generateBN();
            status_vector.clear();

            for (auto it = pattern.candidates[i].candidate.begin(); it != pattern.candidates[i].candidate.end(); ++it, ++idx)
            {
                if (found_results[i].find(*it) != found_results[i].end())
                {
                    continue;
                }

                vector<VertexID> embedding(pattern.size());
                embedding[i] = *it; // u_id -> v_id
                MatchingStatus ms(pattern.size());
                int value = gmatch_engine.execute(i, idx, embedding, ms);

                if (value == -3)
                {
                    ms.vid = *it;
                    ms.v_idx = idx;
                    status_vector.push_back(ms);

                    cout << "ADD INTO TIMEOUT QUEUE..." << endl;
                }
                else if (value == -2)
                {
                    unordered_set<VertexID> &auto_nodes = autos[i];
                    for (auto auto_node : auto_nodes)
                    {
                        pattern.non_candidates[auto_node].insert(*it);

                        // # fail > # total - support
                        // if(pattern.non_candidates[auto_node].size() - pre_non_candidates[auto_node].size() > pattern.candidates[auto_node].size() - nsupport_)
                        // {
                        //     gmatch_engine.clear_table();
                        //     gmatch_engine.clear_match();
                        //     gmatch_engine.reset();
                        //     return;
                        // }
                    }
                }
                else if (value == -1)
                {
                    for (ui j = 0; j < pattern.size(); ++j)
                    {
                        VertexID matched_vertex = embedding[j];
                        unordered_set<VertexID> &auto_nodes = autos[j];
                        for (auto auto_node : auto_nodes)
                        {
                            found_results[auto_node].insert(matched_vertex);
                        }
                    }
                    if (found_results[i].size() >= nsupport_)
                    {
                        break; // move to the next variable
                    }
                }
            }

            // into TimeOut search
            if (found_results[i].size() + status_vector.size() < nsupport_)
            {   // quick check
                gmatch_engine.clear_table();
                gmatch_engine.clear_match();
                gmatch_engine.reset();
                return;
            }

            // decompose
            vector<vector<subPattern> > maps;
            decomposer.decompose(pattern, maps);

            for (ui j = 0; j < status_vector.size(); ++j)
            {
                if (found_results[i].size() + status_vector.size() - j < nsupport_)
                {   // quick check
                    gmatch_engine.clear_table();
                    gmatch_engine.clear_match();
                    gmatch_engine.reset();
                    return;
                }

                bool existence = true;
                MatchingStatus &ms = status_vector[j];


                for (ui k = 0; k < maps.size(); ++k)
                {
                    vector<subPattern> &edge_removed = maps[k];
                    for (auto it = edge_removed.begin(); it != edge_removed.end(); ++it)
                    {
                        Pattern &key = (*it).pattern;
                        vector<VertexID> &subgraph_mappings = (*it).mapping;
                        ui subgraph_size = key.size();

                        // search mappings of i
                        int map_index = search_mapping(subgraph_mappings, i);

                        if (map_index == -1)
                            continue;

                        vector<VertexID> embedding(subgraph_size);
                        embedding[map_index] = ms.vid;

                        key.candidates.resize(subgraph_size);

                        for (ui l = 0; l < subgraph_size; ++l)
                        {
                            key.candidates[l].candidate = pattern.candidates[subgraph_mappings[l]].candidate; // deep copy
                        }

                        exist_engine.set(&pruned_graph, &key);
                        existence = exist_engine.searchExistence(map_index, ms.v_idx, embedding);
                        exist_engine.clear_all();

                        if (!existence)
                        {
                            pattern.non_candidates[i].insert(ms.vid);

                            break; // ms.vid fails
                        }
                    }
                    if (!existence)
                        break;
                }

                if (!existence)
                    continue;

                // continue the search
                vector<VertexID> embedding(pattern.size());
                embedding[i] = ms.vid;
                int value = gmatch_engine.resumeSearch(i, ms, embedding);

                if (value == -2)
                {
                    unordered_set<VertexID> &auto_nodes = autos[i];
                    for (auto auto_node : auto_nodes)
                    {
                        pattern.non_candidates[auto_node].insert(ms.vid);

                        // # fail > # total - support
                        // if(pattern.non_candidates[auto_node].size() - pre_non_candidates[auto_node].size() > pattern.candidates[auto_node].size() - nsupport_)
                        // {
                        //     gmatch_engine.clear_table();
                        //     gmatch_engine.clear_match();
                        //     gmatch_engine.reset();
                        //     return;
                        // }
                    }
                }
                else
                { // value == -1
                    for (ui j = 0; j < pattern.size(); ++j)
                    {
                        VertexID matched_vertex = embedding[j];
                        unordered_set<VertexID> &auto_nodes = autos[j];
                        for (auto auto_node : auto_nodes)
                        {
                            found_results[auto_node].insert(matched_vertex);
                        }
                    }
                    if (found_results[i].size() >= nsupport_)
                    {
                        break; // move to the next variable
                    }
                }
            }
            gmatch_engine.clear_match();
        }
        if (found_results[i].size() < nsupport_)
        {
            return;
        }

        pre_computed[i] = 1;
    }

    gmatch_engine.clear_table();
    gmatch_engine.reset();

    // remove non_candidates from candidates
    for(ui i = 0; i < pattern.size(); ++i) 
    {
        auto & cands = pattern.candidates[i].candidate;
        auto & non_cands = pattern.non_candidates[i];
        vector<VertexID> new_cands;
        for(auto it = cands.begin(); it != cands.end(); ++it)
        {
            if(non_cands.find(*it) == non_cands.end())
            {
                new_cands.push_back(*it);
            }
        }
        cands.swap(new_cands);
    }
}

VertexID GraMi::has_been_computed(VectorOfSet &autos, vector<bool> &preComputed, VertexID index)
{   
    unordered_set<VertexID> &auto_nodes = autos[index];
    for (auto it = auto_nodes.begin(); it != auto_nodes.end(); ++it)
    {
        if (preComputed[*it])
            return *it;
    }
    return index; //else return the same index
}

int GraMi::search_mapping(vector<VertexID> &mapping, VertexID u)
{
    for (ui i = 0; i < mapping.size(); ++i)
    {
        if (mapping[i] == u)
            return i;
    }
    return -1;
}

void GraMi::report(Pattern &pattern)
{   
    fout << "=====================================" << std::endl;
    for (ui i = 0; i < pattern.size(); ++i)
    {
        vertex_t *vertex = pattern.get_p_vertex(i);
        fout << "v " << vertex->id << " " << vertex->label << std::endl;
    }
    for (ui i = 0; i < pattern.dfscodes.size(); ++i)
    {
        fout << "e " << pattern.dfscodes[i].from << " " << pattern.dfscodes[i].to
             << " " << pattern.dfscodes[i].edge_label << std::endl;
    }
}

bool GraMi::isCan(Pattern &pattern)
{
	Graph MIN_GRAPH; //@wenwen: reuse Graph
	Pattern MIN_PATTERN; //@wenwen: reuse Pattern
	DfsCodes& codes = pattern.dfscodes;
	if (codes.size() == 1)
		return (true);

	MIN_GRAPH.toGraph(pattern);
	MIN_PATTERN.dfscodes.clear ();

	RootForwardEdge_Projected root;
	EdgeList	edges; //@??? edge list?

	for (VertexID from = 0; from < MIN_GRAPH.size() ; ++from) 
    {
		if (get_forward_root (MIN_GRAPH, (MIN_GRAPH.vertices_)[from], edges))
        {
            for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
            {
                RootForwardEdge edge((MIN_GRAPH.vertices_)[from].label,(*it)->label, (MIN_GRAPH.vertices_)[(*it)->to].label);
                root[edge].push(*it);
            }
        }
    }
	RootForwardProjected_iter fromlabel = root.begin();
	//the first edge in root is the minimum edge, so push it into MIN_PATTERN
	RootForwardEdge edge = fromlabel->first;

	MIN_PATTERN.dfscodes.emplace_back (0, 1, edge.fromLabel, edge.elabel, edge.toLabel);

	return project_is_min(fromlabel->second, MIN_GRAPH, MIN_PATTERN, pattern);
}


void GraMi::verify(Pattern &pattern)
{   
    for (ui i = 0; i < pattern.size(); ++i)
    {
        vertex_t *vertex = pattern.get_p_vertex(i);
        cout << vertex->label << " ";
    }

    cout << endl;
}