#pragma once

#include "graph.h"

#include <deque>
using std::deque;


class Decompose
{
public:
    void decompose(Pattern &pattern, vector<vector<subPattern> > &mappings);

    int indexOf(vector<VertexID> &vec, VertexID u);

    void insert(EdgeType &vec, VertexID from, eLabel edge_label, VertexID to, ui edgeId);
};

void Decompose::decompose(Pattern &pattern, vector<vector<subPattern> > &mappings)
{
    vector<bool> color;
    deque<VertexID> DFSstack;

    // iterate over each edge
    for (ui i = 0; i < pattern.get_nedges(); ++i)
    {
        color.assign(pattern.size(), false);
        EdgeID currentEdge = i;

        VertexID idA = pattern.edge2vertex[currentEdge].from;
        VertexID idB = pattern.edge2vertex[currentEdge].to;

        vector<subPattern> currentEdgeMappings;

        for (ui j = 0; j < pattern.size(); ++j)
        {
            if (color[j])
                continue;

            subPattern sub_pattern;
            ui edgeId = 0;

            DFSstack.clear();
            DFSstack.push_front(j);

            sub_pattern.pattern.vertices_.emplace_back(0, pattern.get_vertex_label(j));

            sub_pattern.mapping.push_back(j);

            while (!DFSstack.empty())
            {
                VertexID currentNode = DFSstack.front();
                DFSstack.pop_front();

                if (color[currentNode])
                    continue;
                color[currentNode] = true;

                ui degree = pattern.get_vertex_degree(currentNode);

                for (ui k = 0; k < degree; ++k)
                {
                    VertexID nbr = pattern.get_p_vertex(currentNode)->edges[k].to;
                    if ((currentNode == idA && nbr == idB) || (currentNode == idB && nbr == idA))
                        continue;
                    if (color[nbr])
                        continue;

                    VertexID currentNode_id = indexOf(sub_pattern.mapping, currentNode);
                    int nbr_id = indexOf(sub_pattern.mapping, nbr);
                    if (nbr_id == -1)
                    {
                        nbr_id = sub_pattern.mapping.size();
                        sub_pattern.pattern.vertices_.emplace_back(nbr_id, pattern.get_vertex_label(nbr));
                        sub_pattern.mapping.push_back(nbr);
                    }
                    eLabel edge_label = pattern.get_p_vertex(currentNode)->edges[k].label;
                    insert(sub_pattern.pattern.vertices_[currentNode_id].edges, currentNode_id, edge_label, nbr_id, edgeId);
                    insert(sub_pattern.pattern.vertices_[nbr_id].edges, nbr_id, edge_label, currentNode_id, edgeId);
                    sub_pattern.pattern.edge2vertex.emplace_back(currentNode_id, edge_label, nbr_id, edgeId);

                    DFSstack.push_front(nbr);
                    edgeId ++;
                }
            }
            if (sub_pattern.pattern.vertices_.size() > 1 && indexOf(sub_pattern.mapping, pattern.size() - 1) >= 0) 
            {
                currentEdgeMappings.push_back(sub_pattern);
            }
        }
        mappings.push_back(currentEdgeMappings);
    }
}

/**
 * return u's index inside vec, return -1 if not existed.
 */
int Decompose::indexOf(vector<VertexID> &vec, VertexID u)
{
    for (int i = 0; i < vec.size(); ++i)
    {
        if (vec[i] == u)
            return i;
    }
    return -1;
}

void Decompose::insert(EdgeType &vec, VertexID from, eLabel edge_label, VertexID to, ui edgeId)
{   
    auto it = vec.begin();
    for ( ; it != vec.end(); ++it)
    {
        if(it->to > to)
        {
            vec.emplace(it, from, edge_label, to, edgeId);
            return;
        } 
    }
    vec.emplace(it, from, edge_label, to, edgeId);
}