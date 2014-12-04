#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <utility>

#include <libcuckoo/cuckoohash_config.h> // for SLOT_PER_BUCKET
#include <libcuckoo/cuckoohash_map.hh>
//#include <libcuckoo/city_hasher.hh>
#include "test_util.cc"

typedef uint32_t KeyType;
typedef uint32_t ValType;
typedef std::pair<KeyType, ValType> KVPair;

size_t start_power = 10;

size_t thread_num = sysconf(_SC_NPROCESSORS_ONLN);
size_t numkeys = (1U << 10) * SLOT_PER_BUCKET;

bool use_strings = false;

class InsertFindEnvironment {
public:
    InsertFindEnvironment() : smalltable(2 << start_power) {
        // Sets up the random number generator
        uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::cout << "seed = " << seed << std::endl;
        std::uniform_int_distribution<ValType> v_dist(std::numeric_limits<ValType>::min(), std::numeric_limits<ValType>::max());
        std::mt19937_64 gen(seed);

        keys.resize(numkeys);
        vals.resize(numkeys);
        nonkeys.resize(numkeys);

        // Inserting elements into the table
        for (size_t i = 0; i < numkeys; i++) {
            keys[i] = i;
            vals[i] = v_dist(gen);
            //std::cout << "Inserting key #" << i << std::endl;
            //EXPECT_TRUE(smalltable.insert(keys[i], vals[i]));
        }

        std::vector<std::thread> threads;
        size_t keys_per_thread = numkeys / thread_num;
        for (size_t i = 0; i < thread_num - 1; i++) {
            threads.emplace_back(insert_thread_with_val<KeyType, ValType>, 
                                 std::ref(smalltable), keys.begin()+i*keys_per_thread, 
                                 keys.begin()+(i+1)*keys_per_thread,
                                 vals.begin()+i*keys_per_thread, 
                                 vals.begin()+(i+1)*keys_per_thread);
        }

        threads.emplace_back(insert_thread_with_val<KeyType, ValType>, 
                                 std::ref(smalltable), keys.begin()+(thread_num-1)*keys_per_thread, 
                                 keys.begin()+numkeys,
                                 vals.begin()+(thread_num-1)*keys_per_thread, 
                                 vals.begin()+numkeys);
        
        for (size_t i = 0; i < threads.size(); i++) {
            threads[i].join();
        }
        std::cout << "Finished joining all threads, size of table is" << smalltable.size() << std::endl;
        // Fills up nonkeys with keys that aren't in the table
        std::uniform_int_distribution<KeyType> k_dist(std::numeric_limits<KeyType>::min(), std::numeric_limits<KeyType>::max());
        for (size_t i = 0; i < numkeys; i++) {
            KeyType k;
            do {
                k = k_dist(gen);
            } while (k < numkeys);
            nonkeys[i] = k;
        }
    }

    cuckoohash_map<KeyType, ValType, CityHasher<KeyType> > smalltable;
    std::vector<KeyType> keys;
    std::vector<ValType> vals;
    std::vector<KeyType> nonkeys;
};

InsertFindEnvironment* env;

// Makes sure that we can find all the keys with their matching values
// in the small and big tables
void FindKeysInTables() {
    ASSERT_EQ(env->smalltable.size(), numkeys);
    ValType retval;
    for (size_t i = 0; i < numkeys; i++) {
        EXPECT_TRUE(env->smalltable.find(env->keys[i], retval));
        EXPECT_EQ(retval, env->vals[i]);
    }
}

// Makes sure than none of the nonkeys are in either table
void FindNonkeysInTables() {
    ValType retval;
    for (size_t i = 0; i < numkeys; i++) {
        EXPECT_FALSE(env->smalltable.find(env->nonkeys[i], retval));
    }
}

int main(int argc, char **argv) {
    const char* args[] = {"--start-power", "--num-keys", "--thread-num"};
    size_t* arg_vars[] = {&start_power, &numkeys, &thread_num};
    const char* arg_help[] = {"The power argument given to the hashtable during initialization",
                              "The total number of keys to insert",
                              "The number of threads to insert with"};
    const char* flags[] = {"--use-strings"};
    bool* flag_vars[] = {&use_strings};
    const char* flag_help[] = {"If set, the key type of the map will be std::string"};
    parse_flags(argc, argv, "A benchmark for inserts", args, arg_vars, arg_help,
                sizeof(args)/sizeof(const char*), flags, flag_vars, flag_help,
                sizeof(flags)/sizeof(const char*));

    env = new InsertFindEnvironment;

    std::cout << "Running FindKeysInTables" << std::endl;
    FindKeysInTables();
    std::cout << "Running FindNonkeysInTables" << std::endl;
    FindNonkeysInTables();
}
