#include "grami.h"
#include "worker.h"
#include "setting.h"

#include <numeric>

int Settings::maxNumNodes;

using namespace std::chrono;

int main(int argc, char *argv[])
{   
    string fileName;
    int support, thread_num;
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

    //parameter to set the maximum subgraph size (in terms of the number of vertices)
	char * argMaxNodes = getCmdOption(argv, argv + argc, "-maxNodes");
	if(argMaxNodes)
	{
		Settings::maxNumNodes = atoi(argMaxNodes);
	}
	else
		Settings::maxNumNodes = -1;

    //get user-given number of threads
    char * argThreads = getCmdOption(argv, argv + argc, "-thread");
	if(argThreads)
	{
		thread_num = atoi(argThreads);
	}

    auto time1 = steady_clock::now();

    grami.nsupport_ = support;
    grami.pruned_graph.nsupport_ = support;
   
    Worker worker(thread_num);

    worker.load_data(support, fileName);

    auto time2 = steady_clock::now();

    cout << "[TIME] Load Graph time: " << (float)duration_cast<milliseconds>(time2 - time1).count()/1000 << " s" << endl;

    worker.run();

    auto time3 = steady_clock::now();
    cout << "[TIME] Mining time: " << (float)duration_cast<milliseconds>(time3 - time2).count()/1000 << " s" << endl;

    cout << "[TIME] Total time: " << (float)duration_cast<milliseconds>(time3 - time1).count()/1000 << " s" << endl;

    cout << "[INFO] # Frequent Patterns: " << std::accumulate(results_counter.begin(), results_counter.end(), 0) << endl;

    return 0;
}
