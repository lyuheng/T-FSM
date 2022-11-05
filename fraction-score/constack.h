//########################################################################
//## Copyright 2022 Da Yan http://www.cs.uab.edu/yanda
//##
//## Licensed under the Apache License, Version 2.0 (the "License");
//## you may not use this file except in compliance with the License.
//## You may obtain a copy of the License at
//##
//## //http://www.apache.org/licenses/LICENSE-2.0
//##
//## Unless required by applicable law or agreed to in writing, software
//## distributed under the License is distributed on an "AS IS" BASIS,
//## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//## See the License for the specific language governing permissions and
//## limitations under the License.
//########################################################################

#ifndef CONSTACK_H
#define CONSTACK_H

#include <mutex>
#include <stack>
#include <vector>
#include <iostream>

using namespace std;

template <typename T> 
class constack
{
public:
    stack<T> s;
    mutex mtx;

    void enstack(T value) 
    {
        mtx.lock();
        s.push(value);
        mtx.unlock();
    }

    void enstack(vector<T> & val_vec)
    {
        mtx.lock();
        for(T &el: val_vec)
        {
            s.push(el);
        }
        mtx.unlock();
    }

    bool destack(T & to_get) 
    {
        mtx.lock();
        if(!s.empty())
        {
            to_get = s.top();
            s.pop();
            mtx.unlock();
            return true;
        }
        else
        {
            mtx.unlock();
            return false;
        }
    }

    bool empty()
    {
    	mtx.lock();
        bool ret = s.empty();
        mtx.unlock();
        return ret;
    }
};

#endif
