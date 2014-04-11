/* Tests the # of retries for finds with one insert thread running for a
 * specific load range in a partially-filled table */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <stdint.h>
#include <sys/time.h>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>

#include <libcuckoo/cuckoohash_config.h> // for SLOT_PER_BUCKET
#include <libcuckoo/cuckoohash_map.hh>
//#include <libcuckoo/city_hasher.hh>
#include "test_util.cc"

//typedef uint32_t KeyType;
typedef std::string KeyType2;
typedef uint32_t ValType;

typedef uint64_t KeyType;

// The power argument passed to the hashtable constructor. This can be
// set with the command line flag --power.
size_t power = 22;
// The number of threads spawned for inserts. This can be set with the
// command line flag --thread-num-inserts
size_t thread_num_inserts = sysconf(_SC_NPROCESSORS_ONLN);
// The number of threads spawned for finds. This can be set with the 
// command line flag --thread-num-finds
size_t thread_num_finds = sysconf(_SC_NPROCESSORS_ONLN);
// The load factor to fill the table up to before testing throughput.
// This can be set with the command line flag --begin-load.
size_t begin_load = 0;
// The maximum load factor to fill the table up to when testing
// throughput. This can be set with the command line flag
// --end-load.
size_t end_load = 90;
// The seed which the random number generator uses. This can be set
// with the command line flag --seed
size_t seed = 0;

// Whether to use strings as the key
bool use_strings = false;

template <class KType>
void insert_thread(cuckoohash_map<KType, ValType, CityHasher<KType> >& table,
                   typename std::vector<KType>::iterator begin,
                   typename std::vector<KType>::iterator end) {
    for (;begin != end; begin++) {
        if (table.hashpower() > power) {
            std::cerr << "Expansion triggered" << std::endl;
            exit(1);
        }
        //std::cout << "Inserting" << *begin << std::endl;
        ASSERT_TRUE(table.insert(*begin, 0));

    }
}


template <class KType>
class AllEnvironment {
public: 
    AllEnvironment()
        : numkeys((1U << power) * SLOT_PER_BUCKET), table(numkeys), keys(numkeys), nonkeys(numkeys), finished(false), num_reads(){
            // Sets up the random number generator
            if (seed == 0) {
                seed = std::chrono::system_clock::now().time_since_epoch().count();
            }
            //std::cout << "seed = " << seed << std::endl;
            gen.seed(seed);

            // We fill the keys array with integers between numkeys and
            // 2*numkeys, shuffled randomly
            keys[0] = numkeys;
            for (size_t i = 1; i < numkeys; i++) {
                const size_t swapind = gen() % i;
                keys[i] = keys[swapind];
                keys[swapind] = generateKey<KType>(i+numkeys);
            }

            // Fills up nonkeys with keys that aren't in the table
            for (size_t i = 0; i < numkeys; i++) {
                size_t ind = gen();
                do {
                    ind = gen();
                } while (ind <= 2*numkeys);
                nonkeys[i] = generateKey<KType>(ind);
            }

            // prefill table to given begin_load factor
            std::vector<std::thread> threads;
            size_t keys_per_thread = numkeys *(begin_load/100.0)/thread_num_inserts;
            for(size_t i = 0; i < thread_num_inserts; i++) {
                threads.emplace_back(insert_thread<KType>, std::ref(table), keys.begin() + i*keys_per_thread, keys.begin() + (i+1)*keys_per_thread);
            }
            for(size_t i = 0; i < threads.size(); i++) {
                threads[i].join();
            }
            init_size = table.size();
            assert(init_size == keys_per_thread*thread_num_inserts);
            //std::cout << "Initial load factor: " << init_size << std::endl;
        }

    size_t numkeys;
    cuckoohash_map<KType, ValType, CityHasher<KType> > table;
    std::vector<KType> keys;
    std::vector<KType> nonkeys;
    std::mt19937_64 gen;
    std::atomic<bool> finished;
    cacheint num_reads;
    size_t init_size;
};

template <class KType>
void RetryTest(AllEnvironment<KType> *env) {
    std::vector<std::thread> insert_threads;
    std::vector<std::thread> find_threads;
    size_t keys_per_insert_thread = env->numkeys * ((end_load-begin_load) / 100.0) / thread_num_inserts;
    
    for (size_t i = 0; i < thread_num_inserts; i++) {
        insert_threads.emplace_back(insert_thread<KType>, std::ref(env->table), 
                                                   env->keys.begin() + env->init_size + i*keys_per_insert_thread,
                                                   env->keys.begin() + env->init_size + (i+1)*keys_per_insert_thread);
    }

    size_t keys_per_find_thread = env->numkeys / thread_num_finds;
    for (size_t i = 0; i < thread_num_finds; i++) {
        find_threads.emplace_back(read_thread<KType, ValType>, std::ref(env->table),
                                 env->nonkeys.begin() + i*keys_per_find_thread,
                                 env->nonkeys.begin() + (i+1)*keys_per_find_thread,
                                 std::ref(env->num_reads),
                                 false,
                                 std::ref(env->finished));
    }
    for (size_t i = 0; i < insert_threads.size(); i++) {
        insert_threads[i].join();
    }
    env->finished.store(true);
    for (size_t i = 0; i < find_threads.size(); i++) {
        find_threads[i].join();
    }
    
    size_t num_inserts = env->table.size() - env->init_size;
    assert(num_inserts == keys_per_insert_thread*thread_num_inserts);
    size_t num_read_retries = env->table.number_retries();
    size_t num_reads = env->num_reads.num.load();
    // Reports the results
    /*std::cout << "----------Results----------" << std::endl;
    std::cout << "Final load factor:\t" << env->table.load_factor() << std::endl;
    std::cout << "Number of inserts:\t" << num_inserts << std::endl;
    std::cout << "Number of reads:\t" << num_reads << std::endl;
    std::cout << "Number of read retries:\t" << num_read_retries << std::endl;
    std::cout << "Avg # of retries: " << (double) num_read_retries / num_reads << std::endl;
    std::cout << "Percent inserts: " << (double) num_inserts/ (num_inserts + num_reads) << std::endl;*/
    std::cout << (begin_load + end_load)/2 << " , " << num_inserts << " , " << num_reads << " , " << num_read_retries << " , " <<
    thread_num_finds << " , " << thread_num_inserts << " , " <<
    (double) num_read_retries / num_reads << " , " << (double) num_inserts/ (num_inserts + num_reads) << std::endl;
}

int main(int argc, char** argv) {
    const char* args[] = {"--power", "--thread-num-inserts", "--thread-num-finds", "--begin-load",
     "--end-load", "--seed"};
    size_t* arg_vars[] = {&power, &thread_num_inserts, &thread_num_finds, &begin_load, &end_load, &seed};
    const char* arg_help[] = {"The power argument given to the hashtable during initialization",
                              "The number of threads to spawn for inserts",
                              "The number of threads to spawn for finds",
                              "The load factor to fill the table up to before testing",
                              "The maximum load factor to fill the table up to when testing",
                              "The seed used by the random number generator"};
    const char* flags[] = {"--use-strings"};
    bool* flag_vars[] = {&use_strings};
    const char* flag_help[] = {"If set, the key type of the map will be std::string"};
    parse_flags(argc, argv, "Testing nmber of retries", args, arg_vars, arg_help,
                sizeof(args)/sizeof(const char*), flags, flag_vars, flag_help,
                sizeof(flags)/sizeof(const char*));

    if (begin_load >= 100) {
        std::cerr << "--begin-load must be between 0 and 99" << std::endl;
        exit(1);
    } else if (begin_load >= end_load) {
        std::cerr << "--end-load must be greater than --begin-load" << std::endl;
        exit(1);
    }

    std::cout << "Size of cacheint" << sizeof(std::atomic<size_t>) << std::endl;

    if (use_strings) {
        auto *env = new AllEnvironment<KeyType2>;
        RetryTest(env);
        delete env;
    } else {
        auto *env = new AllEnvironment<KeyType>;
        RetryTest(env);
        delete env;
    }
}
