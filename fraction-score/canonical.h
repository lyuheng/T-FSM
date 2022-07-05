#pragma once
#include "graph.h"

//forward edge, used for generate children
struct RootForwardEdge
{
	vLabel fromLabel;
	eLabel elabel;
	vLabel toLabel;

	RootForwardEdge(vLabel from, eLabel e, vLabel to): elabel(e), toLabel(to), fromLabel(from) {};

	friend bool operator<(const RootForwardEdge& e1, const RootForwardEdge& e2)
	{
	    if (e1.fromLabel != e2.fromLabel)
			return e1.fromLabel < e2.fromLabel; 
		else 
		{
			if (e1.elabel != e2.elabel) 
				return e1.elabel < e2.elabel;
			else 
				return e1.toLabel < e2.toLabel;
				
		}
	}

	friend bool operator == (const RootForwardEdge &e1, const RootForwardEdge &e2)
	{
		return (e1.fromLabel == e2.fromLabel && e1.toLabel == e2.toLabel && e1.elabel == e2.elabel);
	}
};

struct ForwardEdge
{
	eLabel elabel;
	vLabel toLabel;

	ForwardEdge(eLabel e, vLabel to): elabel(e), toLabel(to) {};

	friend bool operator<(const ForwardEdge& e1, const ForwardEdge& e2)
	{
		if(e1.elabel < e2.elabel) return true;
	    return e1.toLabel < e2.toLabel;
	}

	friend bool operator == (const ForwardEdge &e1, const ForwardEdge &e2)
	{
		return (e1.toLabel == e2.toLabel && e1.elabel == e2.elabel);
	}
};

typedef edge_t Edge;

class GspanProjTrans
{
public:
	vector<Edge *> pgraph; //projected graph, elements are edges in a data graph that constitute the right-most path
	GspanProjTrans() {};
	GspanProjTrans(Edge* _e, GspanProjTrans* _prev)
	{
		if(_prev !=0)
			pgraph = _prev->pgraph; //deep copy
		pgraph.push_back(_e);
	};
};

//pdb to certain pattern
class GspanProjDB: public vector<GspanProjTrans> 
{
public:
	//int sup; // pdb may contain multiple matched instances of the same transaction
	
	//@wenwen: push for root patern's DB
	void push (Edge *edge)
	{
		resize (size() + 1);
		GspanProjTrans &d = (*this)[size()-1];
		d.pgraph.push_back(edge);
	}
	
	//@wenwen: push for other patern's DB
	void push (Edge *edge, GspanProjTrans &prev)
	{
		resize (size() + 1);
		GspanProjTrans &d = (*this)[size()-1];
		//if(prev!=0)
		d.pgraph = prev.pgraph; //@wenwen: deep copy
		d.pgraph.push_back(edge);
	}
	/*
	//support counting. Each graph will support current pattern only once
	unsigned int support ()
	{
		unsigned int size = 0;
		set<unsigned int> visited;

		for (vector<GspanProjTrans>::iterator cur = begin(); cur != end(); ++cur) {
			int tid = cur->tid;
			if (visited.find(tid) == visited.end()) {
				visited.insert(tid);
				++size;
			}
		}

		sup = size;
		return size;
	}
	*/
};


//store all visited vertex and edge. History will be built with a GspanProj and its corresponding GspanTrans object by build()
//used to avoid expanding an already visited edges in the transaction when growing projected transaction into the pdbs of the children
class History: public vector<Edge *> 
{
private:
	vector<bool> edge;
	vector<bool> vertex;

public:
	bool hasEdge   (ui id) { return edge[id]; }
	bool hasVertex (ui id) { return vertex[id]; }
	void build (Graph &graph, GspanProjTrans &e)
	{
		// first build history
		clear ();
		edge.clear ();
		edge.resize (graph.nedges_); //@wenwen: nedges???
		vertex.clear ();
		vertex.resize (graph.vertices_.size());
		//assert(&e); //@wenwen: test
		{
			for(ui i = 0; i<e.pgraph.size(); i++){
				push_back (e.pgraph[i]);	// this line eats 8% of overall instructions(!)
				edge[e.pgraph[i]->id] = vertex[e.pgraph[i]->from] = vertex[e.pgraph[i]->to] = true;
			}
		}
	}
	History() {};
	History (Graph &g, GspanProjTrans &p) { build (g, p); }
};


typedef map<RootForwardEdge, GspanProjDB>	RootForwardEdge_Projected; // to store all forward edges expansion and their projected database, only use in first round.
typedef RootForwardEdge_Projected::iterator	RootForwardProjected_iter;
typedef map<ForwardEdge, GspanProjDB>	ForwardEdge_Projected; // to store all forward edges expansion and their projected database, only use in first round.
typedef ForwardEdge_Projected::iterator	ForwardProjected_iter;
typedef eLabel BackEdge;
typedef map<BackEdge, GspanProjDB>		BackEdge_Projected;// to store all forward edges expansion and their projected database.
typedef BackEdge_Projected::iterator	BackProjected_iter;

typedef vector<Edge *> EdgeList;

/* get_forward_root ()
 * 	get all forward edges that start from v
*/
bool get_forward_root (Graph &g, vertex_t &v, EdgeList &result)
{
	result.clear ();
	for (EdgeType::iterator it = v.edges.begin(); it != v.edges.end(); ++it) {
		if (v.label <= (g.vertices_)[it->to].label)
			result.push_back (&(*it));
	}
	return (! result.empty());
}

/* get_backward (graph, e1, e2, history);
*	@e1: the edges in the right-most path but e2
*	@e2: the right-most edge.

*    get all back edges that begin from the right-most vertex
 */
Edge *get_backward (Graph &graph, Edge* e1, Edge* e2, History& history)
{
	if (e1 == e2)
		return 0;

	for (EdgeType::iterator it = (graph.vertices_)[e2->to].edges.begin() ;
		it != (graph.vertices_)[e2->to].edges.end() ; ++it)
	{
		if (history.hasEdge (it->id))
			continue;
   	   	/*only extend to the vertex that its edge in the right-most path has a smaller edge label than the extended edge
   	   	or their edge's label are the same but dest vertex's label of the earlier is smaller to avoid redundancy.*/
		if ( (it->to == e1->from) && //which means it is a back edge
        (
            (e1->label < it->label) ||
            (
                (e1->label == it->label) &&
                ((graph.vertices_)[e1->to].label <= (graph.vertices_)[e2->to].label)
            )
        )
        )
		{
			return &(*it);
		}
	}

	return 0;
}

/* get_forward_pure ()
 * @e: the last edge in the rightmost path
 * @result: the output
 *
 * get all forward edges that begin from the rightest vertex
*/
bool get_forward_pure (Graph &graph, Edge *e, vLabel minlabel, History& history, EdgeList &result)
{
	result.clear ();
	/* Walk all edges leaving from vertex e->to.
	 */
	for (EdgeType::iterator it = (graph.vertices_)[e->to].edges.begin();
		it != (graph.vertices_)[e->to].edges.end(); ++it)
	{
		/* -e-> [e->to] -it-> [it->to]
		 */
		/*only extend to the vertex which has a larger label than the begin vertex to avoid redundancy.*/
		if (minlabel > (graph.vertices_)[it->to].label || history.hasVertex (it->to))
			continue;

		result.push_back (&(*it));
	}

	return (! result.empty());
}

/* get_forward_pure ()
 * 	@e: the edge in the rightmost path
 * 	@result: the output
 *
 * 	get all forward edges that start from vertex in the right-most path but the right-most vertex
*/
bool get_forward_rmpath (Graph &graph, Edge *e, vLabel minlabel, History& history, EdgeList &result){
	result.clear ();

	vLabel tolabel = (graph.vertices_)[e->to].label;

	for (EdgeType::iterator it = (graph.vertices_)[e->from].edges.begin() ;
		it != (graph.vertices_)[e->from].edges.end() ; ++it)
	{
		vLabel tolabel2 = (graph.vertices_)[it->to].label;
		/*only extend to edges that source vertex's edge in the right-most path has a smaller edge label than the extended edge or they are
		 * 	 equal but the earlier's dest vertex label is smaller to avoid redundancy.*/
		if (e->to == it->to || minlabel > tolabel2 || history.hasVertex (it->to))
			continue;

		if (e->label < it->label || (e->label == it->label && tolabel <= tolabel2))
			result.push_back (&(*it));
	}

	return (! result.empty());
};

bool project_is_min (GspanProjDB &projected, Graph &MIN_GRAPH, Pattern &MIN_PATTERN, Pattern &pat)
{
	MIN_PATTERN.buildRMPath();
	const vector<int>& rmpath = MIN_PATTERN.right_most_path; //@wenwen: contain the last nodes???
	vLabel minlabel = MIN_PATTERN.dfscodes[0].from_label;
	VertexID maxtoc = MIN_PATTERN.dfscodes[rmpath[0]].to;

	/*===========================begin===========================
	*enum all back edge, forward edge and find the minimum edge. Then store the minimum edge into MIN_PATTERN
	*/
	//first, back edge
	{
		BackEdge_Projected  root; //@wenwen: root[elabel]
		bool flg = false;
		VertexID new_to = 0;
		vLabel new_tolabel = 0;

		for (ui i = rmpath.size()-1; ! flg  && i >= 1; --i) {
			for (ui n = 0; n < projected.size(); ++n) {
				GspanProjTrans &cur = projected[n];
				History history (MIN_GRAPH, cur);
					Edge *e = get_backward (MIN_GRAPH, history[rmpath[i]], history[rmpath[0]], history);
					if (e) {
						//BackEdge edge(MIN_PATTERN.dfscodes[rmpath[i]].from,e->label, MIN_PATTERN.dfscodes[rmpath[i]].from_label);
						root[e->label].push (e, cur);
						new_to = MIN_PATTERN.dfscodes[rmpath[i]].from;
						new_tolabel= MIN_PATTERN.dfscodes[rmpath[i]].from_label;
						flg = true;
					}
				}
			}

			if (flg) {
				BackProjected_iter to = root.begin();
				BackEdge first_elabel = to->first;
				//MIN_PATTERN.dfscodes.emplace_back (maxtoc, edge.to, -1, edge.elabel, -1);
				MIN_PATTERN.dfscodes.emplace_back (maxtoc, new_to, MIN_PATTERN.dfscodes[rmpath[0]].to_label , first_elabel, new_tolabel);
				if (pat.dfscodes[MIN_PATTERN.dfscodes.size()-1] != MIN_PATTERN.dfscodes[MIN_PATTERN.dfscodes.size()-1]) return false;
				return project_is_min (to->second, MIN_GRAPH, MIN_PATTERN, pat);
			}
		}

		//then forward edge
		{
			bool flg = false;
			VertexID new_from = 0;
			vLabel new_fromlabel = 0;
			ForwardEdge_Projected root;
			EdgeList edges;

			for (ui n = 0; n < projected.size(); ++n) {
				GspanProjTrans &cur = projected[n];
				History history (MIN_GRAPH, cur);
				if (get_forward_pure (MIN_GRAPH, history[rmpath[0]], minlabel, history, edges)) {
					flg = true;
					for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it){
						ForwardEdge edge( (*it)->label, (MIN_GRAPH.vertices_)[(*it)->to].label);
						new_from = maxtoc;
						new_fromlabel = (MIN_GRAPH.vertices_)[(*it)->from].label;
						root[edge].push(*it, cur);
					}
				}
			}

			for (ui i = 0; ! flg && i < (int)rmpath.size(); ++i) {
				for (ui n = 0; n < projected.size(); ++n) {
					GspanProjTrans &cur = projected[n];
					History history (MIN_GRAPH, cur);
					if (get_forward_rmpath (MIN_GRAPH, history[rmpath[i]], minlabel, history, edges)) {
						flg = true;
						for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it){
							ForwardEdge edge((*it)->label,(MIN_GRAPH.vertices_)[(*it)->to].label);
							new_from = MIN_PATTERN.dfscodes[rmpath[i]].from;
							new_fromlabel = (MIN_GRAPH.vertices_)[(*it)->from].label;
							root[edge].push(*it, cur);
						}
					}
				}
			}

			if (flg) {
				ForwardProjected_iter from  = root.begin();
				ForwardEdge edge = from->first;
				//MIN_PATTERN.dfscodes.emplace_back (edge.from, maxtoc + 1, -1, edge.elabel, edge.toLabel);
				MIN_PATTERN.dfscodes.emplace_back (new_from, maxtoc + 1, new_fromlabel, edge.elabel, edge.toLabel);
				if (pat.dfscodes[MIN_PATTERN.dfscodes.size()-1] != MIN_PATTERN.dfscodes [MIN_PATTERN.dfscodes.size()-1]) return false;
				return project_is_min (from->second, MIN_GRAPH, MIN_PATTERN, pat);
			}
		}
		/*===========================end===========================*/

		return true;
	}