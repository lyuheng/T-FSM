#include "grami.h"
#include "setting.h"

int Settings::maxNumNodes;

int main(int argc, char *argv[])
{   
    string fileName;
    int support;
    //load graph file
	char * argfilename = getCmdOption(argv, argv + argc, "-file");
	if(argfilename)
	{
		fileName = string(argfilename);
	}

    //get user-given support threshold
	char * argSupport = getCmdOption(argv, argv + argc, "-freq");
	if(argSupport)
	{
		support = atoi(argSupport);
	}

    //parameter to set the maximum subgraph size (in terms of the number of nodes)
	char * argMaxNodes = getCmdOption(argv, argv + argc, "-maxNodes");
	if(argMaxNodes)
	{
		Settings::maxNumNodes = atoi(argMaxNodes);
	}
	else
		Settings::maxNumNodes = -1;


    GraMi grami(support);
    Graph graph(support);

    graph.load_graph(grami.pruned_graph, fileName, " "); // separated by space

    grami.project();

    cout << "done" << endl;

    return 0;
}