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

#ifndef CONMAP_H
#define CONMAP_H

#define CONMAP_BUCKET_NUM 10000 //should be proportional to the number of threads on a machine

//idea: 2-level hashing
//1. id % CONMAP_BUCKET_NUM -> bucket_index
//2. bucket[bucket_index] -> give id, get content

//now, we can dump zero-cache, since GC can directly scan buckets one by one

#include <vector>
#include <ext/hash_map>
#include <ext/hash_set>
#include <mutex>
//#include <map>

#define hash_map __gnu_cxx::hash_map
#define hash_set __gnu_cxx::hash_set
using namespace std;

namespace __gnu_cxx {
template <>
struct hash<long long> {
    size_t operator()(long long val) const
    {
        return (size_t) val;
    }
};
}

template <typename K, typename V> struct conmap_bucket
{
	typedef hash_map<K, V> KVMap;
//	typedef map<K, V> KVMap;
	mutex mtx;
	KVMap bucket;

	inline void lock()
	{
		mtx.lock();
	}

	inline void unlock()
	{
		mtx.unlock();
	}

	KVMap & get_map()
	{
		return bucket;
	}

	//returns true if inserted
	//false if an entry with this key already exists
	bool insert(K key, V val)
	{
		auto ret = bucket.insert(
			std::pair<K, V>(key, val)
		);
		return ret.second;
	}

	//returns whether deletion is successful
	bool erase(K key)
	{
		size_t num_erased = bucket.erase(key);
		return (num_erased == 1);
	}
};

template <typename K, typename V> struct conmap
{
public:
	typedef conmap_bucket<K, V> bucket;
	typedef hash_map<K, V> bucket_map;
	bucket* buckets;

	conmap()
	{
		buckets = new bucket[CONMAP_BUCKET_NUM];
	}

	bucket & get_bucket(K key)
	{
		return buckets[key % CONMAP_BUCKET_NUM];
	}

	bucket & pos(size_t pos)
	{
		return buckets[pos];
	}

	~conmap()
	{
		delete[] buckets;
	}
};

#endif
