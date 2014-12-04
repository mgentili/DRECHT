/* Tests the throughput (queries/sec) of only inserts between a
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
#include <stdlib.h>
typedef uint64_t KeyType;
typedef std::string KeyType2;
typedef uint32_t ValType;

// The power argument passed to the hashtable constructor. This can be
// set with the command line flag --power.
size_t power = 22;
// The number of threads spawned for inserts. This can be set with the
// command line flag --thread-num
size_t thread_num = sysconf(_SC_NPROCESSORS_ONLN);
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

double setup_time = 0;
Json::Value info;
Json::StyledWriter writer;

// Inserts the keys in the given range (with value 0), exiting if there is an expansion
template <class KType>
void insert_thread(cuckoohash_map<KType, ValType, CityHasher<KType> >& table,
                   typename std::vector<KType>::iterator begin,
                   typename std::vector<KType>::iterator end) {
    timeval t1, t2;
    for (;begin != end; begin++) {
        /*if (table.hashpower() > power) {
            std::cerr << "Expansion triggered" << std::endl;
            exit(1);
        }*/
        //table.size();
        gettimeofday(&t1, NULL);
        ASSERT_TRUE(table.insert(*begin, 0));
        gettimeofday(&t2, NULL);
        long int elapsed_time = (t2.tv_sec - t1.tv_sec) * 1000000; // sec to us
        elapsed_time += (t2.tv_usec - t1.tv_usec); // us
        printf("%ld, %ld\n", t1.tv_sec*1000000 + t1.tv_usec, elapsed_time);
        //std::cout <<  << ',' << elapsed_time << std::endl;
    }
}

template <class KType>
class InsertEnvironment {
public:
    // We allocate the vectors with the total amount of space in the
    // table, which is bucket_count() * SLOT_PER_BUCKET
    InsertEnvironment()
        : tablesize((1U << power) * SLOT_PER_BUCKET ), numkeys( tablesize* (1 + (int) end_load/100)), keys(numkeys) {
        // Sets up the random number generator
        if (seed == 0) {
            seed = std::chrono::system_clock::now().time_since_epoch().count();
        }
        //std::cout << "seed = " << seed << std::endl;
        gen.seed(seed);

        timeval t1, t2;
        gettimeofday(&t1, NULL);
        table.initialize(tablesize);
        gettimeofday(&t2, NULL);
        double elapsed_time = (t2.tv_sec - t1.tv_sec) * 1000.0; // sec to ms
        elapsed_time += (t2.tv_usec - t1.tv_usec) / 1000.0; // us to ms
        setup_time = elapsed_time/1000;
        info["setup_time"] = setup_time;
        // We fill the keys array with integers between numkeys and
        // 2*numkeys, shuffled randomly
        keys[0] = numkeys;
        for (size_t i = 1; i < numkeys; i++) {
            const size_t swapind = gen() % i;
            keys[i] = keys[swapind];
            keys[swapind] = generateKey<KType>(i+numkeys);
        }

        /*for (size_t i = 0; i < numkeys; i++) {
                std::cout << "Key#" << i << "is " << keys[i] << std::endl;
        }*/

        // We prefill the table to begin_load with thread_num threads,
        // giving each thread enough keys to insert
        std::vector<std::thread> threads;
        size_t keys_per_thread = tablesize * (begin_load / 100.0) / thread_num;
        for (size_t i = 0; i < thread_num; i++) {
            threads.emplace_back(insert_thread<KType>, std::ref(table), keys.begin()+i*keys_per_thread, keys.begin()+(i+1)*keys_per_thread);
        }
        for (size_t i = 0; i < threads.size(); i++) {
            threads[i].join();
        }

        init_size = table.size();
        ASSERT_TRUE(init_size == keys_per_thread * thread_num);

        info["test_type"] = "insert_throughput";
        info["initial_table_size"] = (int) tablesize;
        info["num_keys"] = (int) numkeys;
        info["begin_load"] = table.load_factor();
        //std::cout << "Table with capacity " << numkeys << " prefilled to a load factor of " << table.load_factor() << std::endl;
        //std::cout << numkeys << ", " << table.load_factor();
    }

    size_t tablesize;
    size_t numkeys;
    cuckoohash_map<KType, ValType, CityHasher<KType> > table;
    std::vector<KType> keys;
    std::mt19937_64 gen;
    size_t init_size;
};

template <class KType>
void InsertThroughputTest(InsertEnvironment<KType> *env) {
    std::vector<std::thread> threads;
    size_t keys_per_thread = env->tablesize * ((end_load-begin_load) / 100.0) / thread_num;
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    for (size_t i = 0; i < thread_num; i++) {
        threads.emplace_back(insert_thread<KType>, std::ref(env->table), env->keys.begin()+(i*keys_per_thread)+env->init_size, env->keys.begin()+((i+1)*keys_per_thread)+env->init_size);
    }
    for (size_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }
    gettimeofday(&t2, NULL);
    double elapsed_time = (t2.tv_sec - t1.tv_sec) * 1000.0; // sec to ms
    elapsed_time += (t2.tv_usec - t1.tv_usec) / 1000.0; // us to ms
    size_t num_inserts = env->table.size() - env->init_size;
    // Reports the results
    /*std::cout << "----------Results----------" << std::endl;
    std::cout << "Final load factor:\t" << env->table.load_factor() << std::endl;
    std::cout << "Number of inserts:\t" << num_inserts << std::endl;
    std::cout << "Time elapsed:\t" << elapsed_time/1000 << " seconds" << std::endl;
    std::cout << "Throughput: " << std::fixed << (double)num_inserts / (elapsed_time/1000) << " inserts/sec" << std::endl;
    std::cout << ", " << env->table.load_factor();
    std::cout << ", " << num_inserts;
    std::cout << ", " << elapsed_time/1000;
    std::cout << ", " << (double)num_inserts / (elapsed_time/1000) << std::endl;*/
    info["end_load"] = env->table.load_factor();
    info["num_inserts"] = (int) num_inserts;
    info["final_table_size"] = (int) ((1U << env->table.hashpower()) * SLOT_PER_BUCKET);
    info["insert_time"] = elapsed_time/1000;
    info["throughput"] = (double)num_inserts / (elapsed_time/1000); 
    double total_time = elapsed_time/1000 + setup_time;
    info["throughput with setup"] = (double)num_inserts/total_time;
    std::cout << writer.write( info ) << std::endl;
}

int main(int argc, char** argv) {
    const char* args[] = {"--power", "--thread-num", "--begin-load", "--end-load", "--seed"};
    size_t* arg_vars[] = {&power, &thread_num, &begin_load, &end_load, &seed};
    const char* arg_help[] = {"The power argument given to the hashtable during initialization",
                              "The number of threads to spawn for each type of operation",
                              "The load factor to fill the table up to before testing throughput",
                              "The maximum load factor to fill the table up to when testing throughput",
                              "The seed used by the random number generator"};
    const char* flags[] = {"--use-strings"};
    bool* flag_vars[] = {&use_strings};
    const char* flag_help[] = {"If set, the key type of the map will be std::string"};
    parse_flags(argc, argv, "A benchmark for inserts", args, arg_vars, arg_help,
                sizeof(args)/sizeof(const char*), flags, flag_vars, flag_help,
                sizeof(flags)/sizeof(const char*));

    if (begin_load >= 100) {
        std::cerr << "--begin-load must be between 0 and 99" << std::endl;
        exit(1);
    } else if (begin_load >= end_load) {
        std::cerr << "--end-load must be greater than --begin-load" << std::endl;
        exit(1);
    }
    if (use_strings) {
        auto *env = new InsertEnvironment<KeyType2>;
        InsertThroughputTest(env);
        delete env;
    } else {
        auto *env = new InsertEnvironment<KeyType>;
        InsertThroughputTest(env);
        delete env;
    }
}
