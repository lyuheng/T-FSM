#include "comper.h"
#include <unistd.h>

class Worker
{  
public:
    Comper *compers = nullptr;
    
    DataStack *data_stack;

    Qlist *activeQ_list;

    Worker(int comper_num)
    {
        global_data_stack = data_stack = new DataStack;
        global_activeQ_list = activeQ_list = new Qlist;

        num_compers = comper_num;
        global_end_label = false;

        results_counter.assign(comper_num, 0);

        // fout = new ofstream[32];
        // for(int i=0; i<32; i++)
        // {
        //     char file[200];
        //     sprintf(file, "./log_%d", i);
        //     fout[i].open(file);
        // }
    }

    ~Worker()
    {
        if (compers)
            delete[] compers;

        delete activeQ_list;
        delete data_stack;

        delete[] grami.pruned_graph.nlf;

        // for(ui i=0; i<32; i++)
        //     fout[i].close();
    }


    void load_data(int support, const string &file_path)
    { 
        Graph graph(support);
        grami.nsupport_ = support;

        graph.load_graph(grami.pruned_graph, file_path, " ", num_compers); // separated by space

        initialize_pattern();
    }

    void initialize_pattern()
    {
        results_counter[0] += grami.initialize();

        vector<task_container *> sep_results;

        // vector<PatternMap::iterator> parallel_vec;

        // for (auto it = grami.init_pattern_map.begin(); it != grami.init_pattern_map.end(); ++it)
        // {
        //     parallel_vec.push_back(it);
        // }

// #pragma omp parallel for schedule(dynamic, 1) num_threads(num_compers)
        for(auto it = grami.init_pattern_map.begin(); it != grami.init_pattern_map.end(); ++it)
        // for (ui i = 0; i < parallel_vec.size(); ++i)
        {
            // auto it = parallel_vec[i];
            // int tid = omp_get_thread_num();

            PatternPVec ext_pattern_vec;

            grami.extend(*(it->second), ext_pattern_vec);

            for(auto it2 = ext_pattern_vec.begin(); it2 != ext_pattern_vec.end(); ++it2) 
            {
                task_container *new_tc = new task_container(qid++);
                new_tc->pattern = *it2;
                new_tc->pattern->non_candidates.resize(new_tc->pattern->size());

                sep_results.push_back(new_tc);
            }

            delete it->second->prog;
            delete it->second;
        }
    
        // for(ui i = 0; i < num_compers; ++i)
        // {
        //     total += sep_results[i].size();
        //     data_stack->enstack(sep_results[i]);
        // }
        data_stack->enstack(sep_results);
        
        cout << "In Initialization, size of data_stack: " << sep_results.size() << endl;
    }

    void run()
    {   
        // create compers
        compers = new Comper[num_compers];
        for (int i = 0; i < num_compers; i++)
        {
            compers[i].start(i);
        }

        while (global_end_label == false)
        {   
            // Avoid busy-checking
            usleep(WAIT_TIME_WHEN_IDLE);

            activeQ_lock.rdlock();
            if(activeQ_num > 0)
            {
                // cout << "Branch A  " << activeQ_num << endl;
            	activeQ_lock.unlock();
            	mtx_go.lock();
				ready_go = true;
				// Release threads to compute tasks
				cv_go.notify_all();
				mtx_go.unlock();
            }
            else
            {
                cout << "Branch B" << endl;
            	activeQ_lock.unlock();
            	if(!data_stack->empty())
                {   
					mtx_go.lock();
					ready_go = true;
					// Release threads to compute tasks
					cv_go.notify_all();
					mtx_go.unlock();
				}
				else 
                {
                    cout << "Branch C" << endl;
					mtx_go.lock();

                    cout << "### " << global_num_idle << " " << num_compers << endl;
					if (global_num_idle == num_compers)
					{
						// every thread is waiting, guaranteed since mtx_go is locked
						// Since we are in else-branch, Qreg must be empty
						cout << "global_num_idle: " << global_num_idle << endl;
						global_end_label = true;
						ready_go = true;
						// Release threads
						cv_go.notify_all();
					}
					// else, some threads are still processing tasks, check in next round
					mtx_go.unlock();
				}
            }
        }
    }
};