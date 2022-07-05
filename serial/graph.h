#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <queue>
#include <algorithm>

#include "types.h"

#define FILE_MAX_LINE 1024
#define INVALID_EDGE_LABEL 100000000

using std::cout;
using std::endl;
using std::queue;
using std::string;
using std::unordered_map;
using std::map;
using std::unordered_set;
using std::vector;

class Pattern;

struct edge_t
{
    edge_t() : from(0), label(0), to(0), id(0) {}

    edge_t(VertexID from, eLabel label, VertexID to, EdgeID id) : from(from), label(label), to(to), id(id) {}

    VertexID from;
    eLabel label;
    VertexID to;
    EdgeID id;
};

typedef vector<edge_t> EdgeType;

struct vertex_t
{
    vertex_t() : id(0), label(0), edges(0) {}
    vertex_t(VertexID id, vLabel label) : id(id), label(label), edges(0) {}

    VertexID id;
    vLabel label;
    EdgeType edges;
};

typedef vector<vertex_t> Vertices;

struct PairVertices
{
    vLabel from_label, to_label;
    eLabel edge_label;

    bool operator==(const PairVertices &pv) const
    {
        return (from_label == pv.from_label) && (to_label == pv.to_label) &&
               (edge_label == pv.edge_label);
    }

    void set(vLabel from_label_, eLabel edge_label_, vLabel to_label_)
    {
        from_label = from_label_;
        edge_label = edge_label_;
        to_label = to_label_;
    }
};

class PairVerticesHashCode 
{
public:
    size_t operator()(const PairVertices& pv) const
    {
        size_t seed = 0;
        hash_combine(seed, pv.from_label);
        hash_combine(seed, pv.edge_label);
        hash_combine(seed, pv.to_label);
        return seed;
    }
};

/**
 * @brief represent candidates of a single edge in G.
 */
struct EdgeCandidate
{
    unordered_set<VertexID> candA, candB;
};

class Graph
{
public:
    int nsupport_;

    VertexID vertex_id;
    EdgeID edge_id;

    int nedges_;
    Vertices vertices_;

    unordered_map<vLabel, vector<VertexID>> nodesByLabel; // nodesByLabel[vLabel] = list of v-ids with vLabel

    unordered_set<eLabel> freqEdgeLabels;

    unordered_map<PairVertices, EdgeCandidate, PairVerticesHashCode> hashedEdges;

    unordered_map<vLabel, int> *nlf;


    Graph() : nedges_(0)
    {
        vertices_.resize(0);
    }

    Graph(int support) : nedges_(0), nsupport_(support)
    {
        vertices_.resize(0);
    }

    int size() const
    {
        return vertices_.size();
    }

    void resize(int s)
    {
        vertices_.resize(s);
    }

    void set_nedges(int nedges)
    {
        nedges_ = nedges;
    }

    int get_nedges() const
    {
        return nedges_;
    }

    void set_vertices(const Vertices &vertices)
    {
        vertices_ = vertices;
    }

    Vertices *get_p_vertices()
    {
        return &vertices_;
    }

    const Vertices *get_p_vertices() const
    {
        return &vertices_;
    }

    vertex_t get_vertex(VertexID index)
    {
        return vertices_[index];
    }

    const vertex_t get_vertex(VertexID index) const
    {
        return vertices_[index];
    }

    vertex_t *get_p_vertex(VertexID index)
    {
        return &vertices_[index];
    }

    const vertex_t *get_p_vertex(VertexID index) const
    {
        return &vertices_[index];
    }

    ui get_vertex_degree(VertexID index)
    {
        return get_p_vertex(index)->edges.size();
    }

    eLabel get_vertex_label(VertexID index)
    {
        return get_p_vertex(index)->label;
    }

    void clear()
    {
        nedges_ = 0;
        vertices_.clear();
    }

    void load_graph(Graph &pruned_graph, const string &input_file, const string &separator);
    void read_input(const string &input_file, const string &separator, vector<vector<string>> &input);
    void find_frequent_labels(unordered_map<eLabel, int> &edgesByLabel,
                              unordered_set<vLabel> &freqNodeLabels);

    void construct_graph(vector<vector<string>> &input, unordered_map<eLabel, int> &edgesByLabel);
    void construct_freq_graph(Graph &pruned_graph, vector<vector<string>> &input,
                              unordered_map<eLabel, int> &edgesByLabel,
                              unordered_set<vLabel> &freqNodeLabels);

    eLabel get_edge_label(VertexID u, VertexID v);

    void build_nlf();

    void toGraph(Pattern &pattern);
};

// DFS Code
struct dfs_code_t
{
    dfs_code_t() : from(0), to(0), from_label(0), edge_label(0), to_label(0) {}

    dfs_code_t(VertexID from, VertexID to, vLabel from_label, eLabel edge_label, vLabel to_label) : from(from), to(to),
                                                                                                    from_label(from_label), edge_label(edge_label), to_label(to_label) {}

    dfs_code_t(const dfs_code_t &other)
    {
        this->from = other.from;
        this->to = other.to;
        this->from_label = other.from_label;
        this->edge_label = other.edge_label;
        this->to_label = other.to_label;
    }

    bool operator!=(const dfs_code_t &t) const
    {
        return (from != t.from) || (to != t.to) ||
               (from_label != t.from_label) || (edge_label != t.edge_label) ||
               (to_label != t.to_label);
    }

    bool operator==(const dfs_code_t &t) const
    {
        return (from == t.from) && (to == t.to) &&
               (from_label == t.from_label) && (edge_label == t.edge_label) &&
               (to_label == t.to_label);
    }

    /** return true if this < other; otherwise return false, c.f. gSpan paper Page10. */
    bool compareTo(const dfs_code_t &other) const
    {

        if (this->from > this->to)
        { // backward
            if (other.from < other.to)
            { // forward
                return true;
            }
            else
            { // both backward
                if (this->to < other.to)
                    return true;
                else if (this->to == other.to)
                {
                    if (this->edge_label < other.edge_label)
                        return true;
                }
            }
        }
        else
        { // forward
            if (other.from < other.to)
            { // both forward
                if (other.from < this->from)
                    return true;
                else if (other.from == this->from)
                {
                    if (this->from_label < other.from_label)
                        return true;
                    else if (this->from_label == other.from_label)
                    {
                        if (this->edge_label < other.edge_label)
                            return true;
                        else if (this->edge_label == other.edge_label)
                        {
                            if (this->to_label < other.to_label)
                                return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    string toString() 
    {
        string txt = std::to_string(from) + "_" + std::to_string(to) + "_" + std::to_string(from_label) + "_"
                        + std::to_string(edge_label) + "_" + std::to_string(to_label);
        return txt;
    } 

    void print()
    {
        cout << from << " " << to << " " << from_label << " " << edge_label << " " << to_label << endl;
    }

    VertexID from;
    VertexID to;
    vLabel from_label;
    eLabel edge_label;
    vLabel to_label;
};

class DfsCodeHashCode 
{
public:
    size_t operator()(const dfs_code_t& code) const
    {   
        size_t seed = 0;
        hash_combine(seed, code.from);
        hash_combine(seed, code.to);
        hash_combine(seed, code.from_label);
        hash_combine(seed, code.edge_label);
        hash_combine(seed, code.to_label);
        return seed;
    }
};

typedef vector<dfs_code_t> DfsCodes;

struct Domain
{
    vector<VertexID> candidate; // ascending order by default, no duplicates

    int size()
    {
        return candidate.size();
    }

    void sort()
    {
        std::sort(candidate.begin(), candidate.end());
    }

    int search(VertexID vid)
    { // binary search
        int l = 0, r = candidate.size() - 1, m;
        while (l <= r)
        {
            m = l + (r - l) / 2;
            if (candidate[m] == vid)
                return m;
            else if (candidate[m] < vid)
                l = m + 1;
            else
                r = m - 1;
        }
        return -1;
    }

    void shrink()
    {
        if (candidate.size() < candidate.capacity() / 2)
        {
            candidate.shrink_to_fit();
        }
    }
};

typedef vector<unordered_set<VertexID>> VectorOfSet;

class Pattern
{
public:
    Pattern *parent;  // parent pointer

    VertexID vertex_id; // vertex ID counter
    EdgeID edge_id;     // edge ID counter
    Vertices vertices_; // pattern graph

    DfsCodes dfscodes;           // dfscodes[] = <from, to, from_lb, edge_lb, to_lb>
    vector<int> right_most_path; // indices inside dfscodes, from bottom to top

    vector<Domain> candidates; // candidates[qryVID] = candidate vertices in qryVID’s domain

    VectorOfSet non_candidates; // invalid assignment, used for push-down pruning, will be added into cache

    unordered_map<EdgeID, edge_t> edge2vertex; // EdgeID to edge entity

    Pattern() {}

    void init(dfs_code_t &dfscode, Pattern *parent_)
    {
        add_edge(dfscode);
        dfscodes.push_back(dfscode);
        right_most_path.push_back(0);
        parent = NULL;
    }

    /** candidates wouldn't be assigned in deep copy */
    // Pattern &operator=(Pattern &p)
    // {
    //     this->parent = &p;
    //     this->children_cnt = 0;
    //     this->vertex_id = p.vertex_id;
    //     this->edge_id = p.edge_id;
    //     this->vertices_ = p.vertices_;
    //     this->dfscodes = p.dfscodes;
    //     this->right_most_path = p.right_most_path;
    //     this->edge2vertex = p.edge2vertex;
    //     return *this;
    // }

    void copy(Pattern &p)
    {
        this->parent = &p;
        this->vertex_id = p.vertex_id;
        this->edge_id = p.edge_id;
        this->vertices_ = p.vertices_;
        this->dfscodes = p.dfscodes;
        this->right_most_path = p.right_most_path;
        this->edge2vertex = p.edge2vertex;
    }

    int size() const
    {
        return vertices_.size();
    }

    void resize(int s)
    {
        vertices_.resize(s);
    }

    ui get_nedges() const
    {
        return edge_id;
    }

    vertex_t *get_p_vertex(VertexID index)
    {
        return &vertices_[index];
    }

    const vertex_t *get_p_vertex(VertexID index) const
    {
        return &vertices_[index];
    }

    ui get_vertex_degree(VertexID index)
    {
        return get_p_vertex(index)->edges.size();
    }

    eLabel get_vertex_label(VertexID index)
    {
        return get_p_vertex(index)->label;
    }

    void insert(EdgeType &edges, VertexID from, eLabel edge_label, VertexID to, ui edgeId);

    void add_edge(dfs_code_t &dfscode);
    void buildRMPath();
    void extend(dfs_code_t &dfscode);

    bool is_acyclic();
    bool distinct_labels();
    eLabel get_edge_label(VertexID u, VertexID v);

    void build_nlf(unordered_map<vLabel, int> *nlf);

    string toString() const;

    void clear_candidate()
    {
        candidates.clear();
    }
};


// ########################################################################################
// =========================== below are Graph methods ====================================
// ########################################################################################

void Graph::load_graph(Graph &pruned_graph, const string &input_file, const string &separator)
{
    vector<vector<string>> input;
    unordered_map<eLabel, int> edgesByLabel;
    unordered_set<vLabel> freqNodeLabels;
    read_input(input_file, separator, input);
    construct_graph(input, edgesByLabel);
    find_frequent_labels(edgesByLabel, freqNodeLabels);
    construct_freq_graph(pruned_graph, input, edgesByLabel, freqNodeLabels);

    pruned_graph.nlf = new unordered_map<vLabel, int>[pruned_graph.size()];
    pruned_graph.build_nlf();
}

void Graph::read_input(const string &input_file, const string &separator, vector<vector<string>> &input)
{
    std::ifstream fin(input_file.c_str());
    char line[FILE_MAX_LINE];

    if (!fin.is_open())
    {
        std::cout << "Can not open the graph file " << input_file << " ." << std::endl;
        exit(-1);
    }

    int num_line = 0;
    while (fin.getline(line, FILE_MAX_LINE))
    {
        char *pch = NULL;
        pch = strtok(line, separator.c_str());
        input.resize(num_line + 1);
        while (pch != NULL)
        {
            input[num_line].emplace_back(pch);
            pch = strtok(NULL, separator.c_str());
        }
        ++num_line;
    }
    fin.close();
}

void Graph::construct_graph(vector<vector<string>> &input, unordered_map<eLabel, int> &edgesByLabel)
{
    int edge_id = 0;
    for (ui i = 0; i < input.size(); ++i)
    {

        if (input[i][0] == "t" || input[i][0] == "#")
            continue; // by SHIXUAN INPUT

        else if (input[i][0] == "v")
        {
            int id = atoi(input[i][1].c_str());
            vLabel label = atoi(input[i][2].c_str());
            vertices_.emplace_back(id, label);
        }
        else if (input[i][0] == "e")
        {
            int from = atoi(input[i][1].c_str());
            int to = atoi(input[i][2].c_str());
            eLabel label = atoi(input[i][3].c_str()); // no edge label by default
            // Add an undirected edge
            // Forward direction edge
            vertices_[from].edges.emplace_back(from, label, to, edge_id);
            // Backward direction edge
            vertices_[to].edges.emplace_back(to, label, from, edge_id);
            ++edge_id;
        }
        else
        {
            std::cout << "Reading input error!" << std::endl;
        }
    }

    set_nedges(edge_id);

    for (ui i = 0; i < size(); ++i)
    {

        // std::sort(vertices_[i].edges.begin(), vertices_[i].edges.end(), [](const edge_t &a, const edge_t &b)
        // { 
        //     return a.to < b.to;  // ascending order 
        // });

        int degree = get_vertex_degree(i);
        vLabel label = vertices_[i].label;
        if (nodesByLabel.find(label) == nodesByLabel.end())
        {
            nodesByLabel[label] = vector<VertexID>();
        }
        nodesByLabel[label].push_back(i);

        for (ui j = 0; j < degree; ++j)
        {
            eLabel edgeLabel = vertices_[i].edges[j].label;
            if (edgesByLabel.find(edgeLabel) == edgesByLabel.end())
            {
                edgesByLabel[edgeLabel] = 0;
            }
            edgesByLabel[edgeLabel]++;
        }
    }
}

void Graph::find_frequent_labels(unordered_map<eLabel, int> &edgesByLabel,
                                 unordered_set<vLabel> &freqNodeLabels)
{
    for (auto &pa : nodesByLabel)
    {
        vLabel label = pa.first;
        int freq = pa.second.size();
        if (freq >= nsupport_)
        {
            freqNodeLabels.insert(label);
        }
    }
    for (auto &pa : edgesByLabel)
    {
        eLabel label = pa.first;
        int freq = pa.second;
        if (freq / 2 >= nsupport_) // undirected graph
        { 
            freqEdgeLabels.insert(label);
        }
    }
}

void Graph::construct_freq_graph(Graph &pruned_graph, vector<vector<string>> &input,
                                 unordered_map<eLabel, int> &edgesByLabel,
                                 unordered_set<vLabel> &freqNodeLabels)
{
    edge_id = 0;
    vertex_id = 0;
    vector<vLabel> labels;
    unordered_map<VertexID, VertexID> id_map;

    for (ui i = 0; i < input.size(); ++i)
    {

        if (input[i][0] == "t" || input[i][0] == "#")
            continue;

        else if (input[i][0] == "v")
        {
            VertexID id = atoi(input[i][1].c_str());
            vLabel label = atoi(input[i][2].c_str());
            labels.emplace_back(label);

            if (freqNodeLabels.find(label) != freqNodeLabels.end())
            {
                pruned_graph.vertices_.emplace_back(vertex_id, label);
                id_map[id] = vertex_id; // id_map[old_id] = new_id
                vertex_id++;
            }
        }
        else if (input[i][0] == "e")
        {
            VertexID from = atoi(input[i][1].c_str());
            VertexID to = atoi(input[i][2].c_str());
            eLabel label = atoi(input[i][3].c_str()); // SHIXUAN INPUT;
            vLabel label_from = labels[from];
            vLabel label_to = labels[to];

            if (freqNodeLabels.find(label_from) != freqNodeLabels.end() &&
                freqNodeLabels.find(label_to) != freqNodeLabels.end() &&
                freqEdgeLabels.find(label) != freqEdgeLabels.end())
            {

                pruned_graph.vertices_[id_map[from]].edges.emplace_back(id_map[from], label, id_map[to], edge_id);
                pruned_graph.vertices_[id_map[to]].edges.emplace_back(id_map[to], label, id_map[from], edge_id);
                ++edge_id;
            }

            PairVertices pv;
            if(label_from < label_to)
            {   
                pv.set(label_from, label, label_to);
                pruned_graph.hashedEdges[pv].candA.insert(id_map[from]);
                pruned_graph.hashedEdges[pv].candB.insert(id_map[to]);
            }
            else if(label_from == label_to)
            {
                pv.set(label_from, label, label_to);
                pruned_graph.hashedEdges[pv].candA.insert(id_map[from]);
                pruned_graph.hashedEdges[pv].candA.insert(id_map[to]);
                pruned_graph.hashedEdges[pv].candB.insert(id_map[to]);
                pruned_graph.hashedEdges[pv].candB.insert(id_map[from]);
            }
            else 
            {
                pv.set(label_to, label, label_from);
                pruned_graph.hashedEdges[pv].candA.insert(id_map[to]);
                pruned_graph.hashedEdges[pv].candB.insert(id_map[from]);
            }
        }
        else
        {
            std::cout << "Reading input error!" << std::endl;
        }
    }
    pruned_graph.set_nedges(edge_id);

    for (ui i = 0; i < pruned_graph.size(); ++i)
    {
        vLabel label = pruned_graph.vertices_[i].label;
        if (pruned_graph.nodesByLabel.find(label) == pruned_graph.nodesByLabel.end())
        {
            pruned_graph.nodesByLabel[label] = vector<VertexID>();
        }
        pruned_graph.nodesByLabel[label].push_back(i);
    }

    pruned_graph.freqEdgeLabels = freqEdgeLabels;

    // prune infrequent hasdedEdges
    for (auto it = pruned_graph.hashedEdges.begin(); it != pruned_graph.hashedEdges.end();)
    {
        if (it->second.candA.size() < nsupport_ || it->second.candB.size() < nsupport_)
        {
            it = pruned_graph.hashedEdges.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

eLabel Graph::get_edge_label(VertexID u, VertexID v)
{
    // int l = 0, r = get_p_vertex(u)->edges.size()-1, m;
    // while(l <= r) {
    //     m = l + (r - l)/2;
    //     if(get_p_vertex(u)->edges[m].to == v)
    //         return get_p_vertex(u)->edges[m].label;
    //     else if(get_p_vertex(u)->edges[m].to < v)
    //         l = m + 1;
    //     else
    //         r = m - 1;
    // }
    // std::cout << "Cannot find edge between "<< u << " and " << v << std::endl;
    // return INVALID_EDGE_LABEL;

    ui degree = get_vertex_degree(u);
    for (ui i = 0; i < degree; ++i)
    {
        if (get_p_vertex(u)->edges[i].to == v)
        {
            return get_p_vertex(u)->edges[i].label;
        }
    }
    std::cout << "Cannot find edge between " << u << " and " << v << std::endl;
    return INVALID_EDGE_LABEL;
}

void Graph::toGraph(Pattern &pattern)
{
    vertices_ = pattern.vertices_; // deep copy

    // nodesByLabel
    for (ui i = 0; i < pattern.size(); ++i)
    {
        vLabel label = pattern.get_vertex_label(i);
        if (nodesByLabel.find(label) == nodesByLabel.end())
        {
            nodesByLabel[label] = vector<VertexID>();
        }
        nodesByLabel[label].push_back(i);
    }
    // set_nedges
    set_nedges(pattern.edge_id);
}

void Graph::build_nlf()
{
    for (ui i = 0; i < size(); ++i)
    {
        ui degree = get_p_vertex(i)->edges.size();
        for (ui j = 0; j < degree; ++j)
        {
            VertexID nbr = get_p_vertex(i)->edges[j].to;
            vLabel label = get_vertex_label(nbr);
            if (nlf[i].find(label) == nlf[i].end())
            {
                nlf[i][label] = 0;
            }
            nlf[i][label] += 1;
        }
    }
}

// ########################################################################################
// =========================== below are Pattern methods ==================================
// ########################################################################################

void Pattern::insert(EdgeType &edges, VertexID from, eLabel edge_label, VertexID to, ui edgeId)
{
    auto it = edges.begin();
    for ( ; it != edges.end(); ++it)
    {
        if(it->to > to)
        {
            edges.emplace(it, from, edge_label, to, edgeId);
            return;
        }
    }
    edges.emplace(it, from, edge_label, to, edgeId);
}

void Pattern::buildRMPath()
{
    right_most_path.clear();
    int prev_id = -1;
    for (ui i = dfscodes.size(); i > 0; --i)
    {
        if (dfscodes[i - 1].from < dfscodes[i - 1].to &&
            (right_most_path.empty() || prev_id == dfscodes[i - 1].to))
        {
            prev_id = dfscodes[i - 1].from;
            right_most_path.push_back(i - 1); // reverse order
        }
    }
}

void Pattern::extend(dfs_code_t &dfscode)
{
    add_edge(dfscode);
    dfscodes.push_back(dfscode);
    buildRMPath();
}

void Pattern::add_edge(dfs_code_t &dfscode)
{
    if (size() == 0)
    {
        resize(2);
        vertices_[0].id = 0;
        vertices_[0].label = dfscode.from_label;
        vertices_[1].id = 1;
        vertices_[1].label = dfscode.to_label;
        vertices_[0].edges.emplace_back(0, dfscode.edge_label, 1, 0);
        vertices_[1].edges.emplace_back(1, dfscode.edge_label, 0, 0);
        edge2vertex[0] = vertices_[0].edges.back();
        edge_id = 1;
        vertex_id = 2;
    }
    else
    {
        VertexID from_id = dfscode.from;
        VertexID to_id = dfscode.to;
        if (from_id < vertex_id && to_id < vertex_id)
        {
            // backward edge
            vertices_[from_id].edges.emplace_back(from_id, dfscode.edge_label, to_id, edge_id);
            vertices_[to_id].edges.emplace_back(to_id, dfscode.edge_label, from_id, edge_id);
            edge2vertex[edge_id] = vertices_[from_id].edges.back();
            edge_id++;
        }
        else
        {
            // forward edge
            vertices_.emplace_back(to_id, dfscode.to_label);
            vertices_[from_id].edges.emplace_back(from_id, dfscode.edge_label, to_id, edge_id);
            vertices_[to_id].edges.emplace_back(to_id, dfscode.edge_label, from_id, edge_id);
            edge2vertex[edge_id] = vertices_[from_id].edges.back();
            vertex_id++;
            edge_id++;
        }
    }
}

void Pattern::build_nlf(unordered_map<vLabel, int> *nlf)
{
    for (ui i = 0; i < size(); ++i)
    {
        int degree = get_p_vertex(i)->edges.size();
        for (ui j = 0; j < degree; ++j)
        {
            int nbr = get_p_vertex(i)->edges[j].to;
            vLabel label = get_vertex_label(nbr);
            if (nlf[i].find(label) == nlf[i].end())
            {
                nlf[i][label] = 0;
            }
            nlf[i][label] += 1;
        }
    }
}

bool Pattern::distinct_labels()
{
    unordered_set<vLabel> set;
    for (ui i = 0; i < size(); i++)
    {
        vLabel vlabel = get_vertex_label(i);
        if (set.find(vlabel) == set.end())
        {
            set.insert(vlabel);
        }
        else
        {
            return false;
        }
    }
    return true;
}

bool Pattern::is_acyclic()
{
    int cnt;
    queue<VertexID> q;
    unordered_set<VertexID> visited;
    q.push(0);

    while (!q.empty())
    {
        VertexID cur = q.front();
        q.pop();
        visited.insert(cur);
        cnt = 0;
        int degree = get_p_vertex(cur)->edges.size();
        for (ui i = 0; i < degree; ++i)
        {
            VertexID nbr = get_p_vertex(cur)->edges[i].to;
            if (visited.find(nbr) != visited.end())
            {
                cnt++;
                if (cnt > 1)
                {
                    return false;
                }
            }
            else
            {
                q.push(nbr);
            }
        }
    }
    return true;
}

eLabel Pattern::get_edge_label(VertexID u, VertexID v)
{
    // int l = 0, r = get_p_vertex(u)->edges.size()-1, m;
    // while(l <= r) {
    //     m = l + (r - l)/2;
    //     if(get_p_vertex(u)->edges[m].to == v)
    //         return get_p_vertex(u)->edges[m].label;
    //     else if(get_p_vertex(u)->edges[m].to < v)
    //         l = m + 1;
    //     else
    //         r = m - 1;
    // }
    // std::cout << "Cannot find edge between "<< u << " and " << v << std::endl;
    // return INVALID_EDGE_LABEL;

    ui degree = get_vertex_degree(u);
    for (ui i = 0; i < degree; ++i)
    {
        if (get_p_vertex(u)->edges[i].to == v)
        {
            return get_p_vertex(u)->edges[i].label;
        }
    }
    std::cout << "Cannot find edge between " << u << " and " << v << std::endl;
    return INVALID_EDGE_LABEL;
}

string Pattern::toString() const
{
    string key = "";
    for (ui i = 0; i < size(); ++i)
    {
        key += "v " + std::to_string(get_p_vertex(i)->id) + ' ' + std::to_string(get_p_vertex(i)->label) + '\n';
    }

    for(ui i = 0; i < get_nedges(); ++i)
    {   
        if(edge2vertex.find(i) == edge2vertex.end()) 
            cout << "EdgeID " << i << " is not existed..."<< endl;
        else
            key += "e " + std::to_string(edge2vertex.at(i).from) + ' ' + std::to_string(edge2vertex.at(i).to) + ' ' + std::to_string(edge2vertex.at(i).label) + '\n';
    }
    return key;
}