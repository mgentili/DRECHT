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

typedef uint64_t KeyType;
typedef std::string KeyType2;
typedef uint32_t ValType;

// The power argument passed to the hashtable constructor. This can be
// set with the command line flag --power.
size_t power = 22;
// The number of threads spawned for upserts. This can be set with the
// command line flag --thread-num
size_t thread_num = sysconf(_SC_NPROCESSORS_ONLN);
// The load factor to fill the table up to before testing throughput.
// This can be set with the command line flag --begin-load.
size_t begin_load = 40;
// The seed which the random number generator uses. This can be set
// with the command line flag --seed
size_t seed = 0;

//percentage of workload that are inserts. This can be set with the 
//command line flag --percent-inserts. Right now, the rest are finds.
size_t percent_inserts = 10;

//percentage of finds in the table. This can be set with
//the command line flag --percent-find-in
size_t percent_find_in = 40;

size_t percent_find_out = 40;
size_t percent_deletes = 10;

// Whether to use strings as the key
bool use_strings = false;

size_t test_len = 10;

std::atomic<bool> finished(false);

// Inserts the keys in the given range (with value 0), exiting if there is an expansion
template <class KType>
void insert_thread(cuckoohash_map<KType, ValType, CityHasher<KType> >& table,
                   typename std::vector<KType>::iterator begin,
                   typename std::vector<KType>::iterator end) {
    for (;begin != end; begin++) {
        ASSERT_TRUE(table.insert(*begin, 0));
    }
}

// Inserts keys with probability percent_inserts
// Finds keys in table with probability percent_find_in
// Finds keys not in table with probability 1-percent_inserts-percent_find_in
// Doesn't support multithreading yet
template <class KType>
void mixed_thread(cuckoohash_map<KType, ValType, CityHasher<KType> >& table,
                   typename std::vector<KType>& keys,
                   typename std::vector<KType>& nonkeys,
                   size_t initpos,
                   size_t endpos
                   ) {

    ValType retval;
    size_t insert_pos = initpos;
    size_t delete_pos = initpos;
    size_t find_pos = 0;
    size_t optype_ind = 0;

    while (true) {
        if (finished.load(std::memory_order_acquire)) {
                return;
        }
        if( optype_ind < percent_inserts) {
            ASSERT_TRUE(table.insert(keys[insert_pos], 0));
            insert_pos++;
        } else if( optype_ind < percent_inserts + percent_find_in) {
            ASSERT_TRUE(table.find(keys[delete_pos + find_pos], retval));
            find_pos = (find_pos + 1) % percent_inserts;
        } else if( optype_ind < percent_inserts + percent_find_in + percent_deletes) {
            ASSERT_TRUE(table.erase(keys[delete_pos]));
            delete_pos++;
        } else {
            ASSERT_TRUE(!table.find(keys[delete_pos + find_pos], retval));
            find_pos = (find_pos + 1) % percent_inserts;
        }
        optype_ind++;
        if( optype_ind % 100 == 0) {
            optype_ind = 0;
        }
    }
    //table.list();
    //std::cout << sizeof(KType) << std::endl;
    //std::cout << sizeof(ValType) << std::endl;
}

template <class KType>
class MixedEnvironment {
public:
    MixedEnvironment()
        : numkeys((1U << power) * SLOT_PER_BUCKET), 
        keys(numkeys), nonkeys(numkeys) {

        // Sets up the random number generator
        seed = std::chrono::system_clock::now().time_since_epoch().count();
        //std::cout << "seed = " << seed << std::endl;
        gen.seed(seed);

        table.initialize(numkeys);
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

        // We prefill the table to begin_load with thread_num threads,
        // giving each thread enough keys to upsert
        std::vector<std::thread> threads;
        size_t keys_per_thread = numkeys * (begin_load / 100.0) / thread_num;
        for (size_t i = 0; i < thread_num; i++) {
            threads.emplace_back(insert_thread<KType>, std::ref(table), keys.begin()+i*keys_per_thread, keys.begin()+(i+1)*keys_per_thread);
        }
        for (size_t i = 0; i < threads.size(); i++) {
            threads[i].join();
        }

        init_size = table.size();
        ASSERT_TRUE(init_size == keys_per_thread * thread_num);
        //std::cout << "Table with capacity " << numkeys << " prefilled to a load factor of " << table.load_factor() << std::endl;
        std::cout << numkeys << ", " << table.load_factor();
    }

    size_t numkeys;
    cuckoohash_map<KType, ValType, CityHasher<KType> > table;
    std::vector<KType> keys;
    std::vector<KType> nonkeys;
    std::mt19937_64 gen;
    size_t init_size;
    size_t init_updates;
};

template <class KType>
void MixedTest(MixedEnvironment<KType> *env) {
    std::vector<std::thread> threads;
    size_t keys_per_thread = env->numkeys * ((100-begin_load) / 100.0) / thread_num;
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    for (size_t i = 0; i < thread_num; i++) {
        threads.emplace_back(mixed_thread<KType>, std::ref(env->table), 
            std::ref(env->keys),
            std::ref(env->nonkeys),
            env->init_size + i*keys_per_thread,
            env->init_size + (i+1)*keys_per_thread);
    }
    sleep(test_len);
    finished.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }
    gettimeofday(&t2, NULL);
    double elapsed_time = (t2.tv_sec - t1.tv_sec) * 1000.0; // sec to ms
    elapsed_time += (t2.tv_usec - t1.tv_usec) / 1000.0; // us to ms
    size_t num_inserts = env->table.num_inserts() - env->init_size;
    size_t num_deletes = env->table.num_deletes();
    size_t num_finds_in = num_inserts * percent_find_in / ((double) percent_inserts);
    size_t num_finds_out = num_inserts * percent_find_out / ((double) percent_inserts);
    double throughput = (num_inserts + num_deletes + num_finds_in + num_finds_out)/(elapsed_time/1000);
    // Reports the results
    /*std::cout << "----------Results----------" << std::endl;
    std::cout << "Final load factor:\t" << env->table.load_factor() << std::endl;
    std::cout << "Number of updates:\t" << num_updates << std::endl;
    std::cout << "Number of inserts:\t" << num_inserts << std::endl;
    std::cout << "Number of finds in:\t" << num_finds_in << std::endl;
    std::cout << "Number of finds out:\t" << num_finds_out << std::endl;
    std::cout << "Time elapsed:\t" << elapsed_time/1000 << " seconds" << std::endl;
    std::cout << "Throughput: " << std::fixed 
        << (double)(num_upserts + num_finds_in + num_finds_out)/ (elapsed_time/1000) 
        << " operations/sec" << std::endl;*/
    std::cout << ", " << env->table.load_factor();
    std::cout << ", " << num_inserts;
    std::cout << ", " << num_deletes;
    std::cout << ", " << num_finds_in;
    std::cout << ", " << num_finds_out;
    std::cout << ", " << elapsed_time/1000;
    std::cout << ", " << throughput;
    //std::cout << ", " << 
    //  (double)(num_upserts + num_finds_in + num_finds_out)/(elapsed_time/1000) << std::endl;
}

int main(int argc, char** argv) {
    const char* args[] = {"--power", "--thread-num", "--begin-load", "--seed", "--percent-inserts", "--percent-find-in", "--percent-find-out", "--percent-deletes", "--time"};
    size_t* arg_vars[] = {&power, &thread_num, &begin_load, &seed, &percent_inserts, &percent_find_in, &percent_find_out, &percent_deletes, &test_len};
    const char* arg_help[] = {"The power argument given to the hashtable during initialization",
                              "The number of threads to spawn for each type of operation",
                              "The load factor to fill the table up to before testing throughput",
                              "The seed used by the random number generator",
                              "The percent of inserts in the workload",
                              "The percent of finds that are in the table",
                              "The percent of finds that are not in the table",
                              "The percent of deletes in the workload",
                              "How long to run the test"};

    const char* flags[] = {"--use-strings"};
    bool* flag_vars[] = {&use_strings};
    const char* flag_help[] = {"If set, the key type of the map will be std::string"};
    parse_flags(argc, argv, "A benchmark for upserts", args, arg_vars, arg_help,
                sizeof(args)/sizeof(const char*), flags, flag_vars, flag_help,
                sizeof(flags)/sizeof(const char*));

    if (begin_load >= 100) {
        std::cerr << "--begin-load must be between 0 and 99" << std::endl;
        exit(1);
    } 

    if (percent_inserts >= 100 || percent_inserts < 0) {
        std::cerr << "--percent-inserts must be between 0 and 99" << std::endl;
        exit(1);
    } else if (percent_find_in >= 100 || percent_find_in < 0) {
        std::cerr << "--percent-find-in must be between 0 and 99" << std::endl;
        exit(1);
    } else if (percent_inserts + percent_find_in > 100) {
        std::cerr << "the sum of inserts + finds_in must be <= 100" << std::endl;
        exit(1);
    }

    if (use_strings) {
        auto *env = new MixedEnvironment<KeyType2>;
        MixedTest(env);
        delete env;
    } else {
        auto *env = new MixedEnvironment<KeyType>;
        MixedTest(env);
        delete env;
    }
}
