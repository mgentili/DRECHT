/* Tests the throughput (queries/sec) of inserts and finds between a
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
// The number of threads spawned for upserts. This can be set with the
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

//percentage of workload that are inserts. This can be set with the 
//command line flag --percent-inserts. Right now, the rest are finds.
size_t percent_inserts = 10;

//percentage of finds that are in the table. This can be set with
//the command line flag --percent-find-in
size_t percent_find_in = 80;

// Whether to use strings as the key
bool use_strings = false;

size_t num_finds_in = 0;
size_t num_finds_out = 0;

//enum table_call {INSERT_CALL, FIND_IN_CALL, FIND_OUT_CALL};

template <class KType>
struct insert_thread_args {
    cuckoohash_map<KType, ValType, CityHasher<KType> >& table;
    typename std::vector<KType>::iterator begin;
    typename std::vector<KType>::iterator end;
};

/*template <class KType>
struct insert_find_thread_args {
    cuckoohash_map<KType, ValType, CityHasher<KType> >& table;
    typename std::vector<KType>& keys;
    typename std::vector<KType>& nonkeys;
    cacheint* reads;
    bool in_table;
    std::atomic<bool>* finished;
};*/


// Inserts the keys in the given range (with value 0), exiting if there is an expansion
template <class KType>
void insert_thread(insert_thread_args<KType> it_args) {
    cuckoohash_map<KType, ValType, CityHasher<KType> >& table = it_args.table;
    auto begin = it_args.begin;
    auto end = it_args.end;

    for (;begin != end; begin++) {
        if (table.hashpower() > power) {
            std::cerr << "Expansion triggered" << std::endl;
            exit(1);
        }
        //ASSERT_TRUE(table.upsert(*begin, generateKey<ValType>(gen()) ));
        ASSERT_TRUE(table.insert(*begin, 0));
    }
    //table.list();
}

// Inserts keys with probability percent_inserts
// Finds keys in table with probability percent_find_in
// Finds keys not in table with probability 1-percent_inserts-percent_find_in
// Doesn't support multithreading yet
template <class KType>
void insertfind_thread(cuckoohash_map<KType, ValType, CityHasher<KType> >& table,
                   typename std::vector<KType>& keys,
                   typename std::vector<KType>& nonkeys,
                   size_t table_initial_size,
                   size_t keys_to_add,
                   size_t numkeys
                   ) {

    ValType retval;
    size_t endpos = table_initial_size + keys_to_add;
    size_t initpos = table_initial_size;
    size_t nonkey_pos = 0;
    size_t tablecalls_pos = 0;
    size_t inkey_pos = initpos - 1; //start looking for keys from back forwards
    size_t percent_inserts_finds_in = percent_inserts + percent_find_in;
    //size_t inkey_pos = 0;
    while(initpos != endpos) {
        if (table.hashpower() > power) {
            std::cerr << "Expansion triggered" << std::endl;
            exit(1);
        }
        tablecalls_pos %= 100;

        if(tablecalls_pos < percent_inserts) { //we should insert
            //std::cout << "Inserting" << keys[initpos] << std::endl;
            ASSERT_TRUE(table.insert(keys[initpos], 0));
            initpos++;
        } else if(tablecalls_pos < percent_inserts_finds_in ){ //we should find a key inside table
                /*std::cout << "Finding with inkey_pos:" << inkey_pos 
                    << "initpos: " << initpos << "key:" << keys[inkey_pos % initpos] 
                    << ", which should be in table"<< std::endl;*/
                ASSERT_TRUE(table.find(keys[inkey_pos % initpos], retval));
                inkey_pos--;
                //inkey_pos++;
                num_finds_in++;
        } else {
                //std::cout << "Finding" << nonkeys[nonkey_pos % numkeys] << "which should not be in table"<< std::endl;
                ASSERT_TRUE(!(table.find(nonkeys[nonkey_pos % numkeys], retval)));
                nonkey_pos++;
                num_finds_out++;
        }
        tablecalls_pos++;

        //usleep(100000);
    }
    //table.list();
    //std::cout << sizeof(KType) << std::endl;
    //std::cout << sizeof(ValType) << std::endl;
}

template <class KType>
class InsertFindEnvironment {
public:
    // We allocate the vectors with the total amount of space in the
    // table, which is bucket_count() * SLOT_PER_BUCKET
    InsertFindEnvironment()
        : numkeys((1U << power) * SLOT_PER_BUCKET), table(numkeys), 
        keys(numkeys), nonkeys(numkeys) {

        // Sets up the random number generator
        seed = std::chrono::system_clock::now().time_since_epoch().count();
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

        // We prefill the table to begin_load with thread_num threads,
        // giving each thread enough keys to upsert
        std::vector<std::thread> threads;
        size_t keys_per_thread = numkeys * (begin_load / 100.0) / thread_num;
        for (size_t i = 0; i < thread_num; i++) {
            threads.emplace_back(insert_thread<KType>, insert_thread_args<KType>{std::ref(table), keys.begin()+i*keys_per_thread, keys.begin()+(i+1)*keys_per_thread});
        }
        for (size_t i = 0; i < threads.size(); i++) {
            threads[i].join();
        }
        init_size = table.size();
        init_updates = table.number_updates();
        ASSERT_TRUE((init_size + init_updates) == keys_per_thread * thread_num);
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
void InsertFindTest(InsertFindEnvironment<KType> *env) {
    std::vector<std::thread> threads;
    size_t keys_per_thread = env->numkeys * ((end_load-begin_load) / 100.0) / thread_num;
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    for (size_t i = 0; i < thread_num; i++) {
        threads.emplace_back(insertfind_thread<KType>, std::ref(env->table), 
            std::ref(env->keys),
            std::ref(env->nonkeys),
            env->init_size,
            keys_per_thread,
            env->numkeys);
    }
    for (size_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }
    gettimeofday(&t2, NULL);
    double elapsed_time = (t2.tv_sec - t1.tv_sec) * 1000.0; // sec to ms
    elapsed_time += (t2.tv_usec - t1.tv_usec) / 1000.0; // us to ms
    size_t num_updates = env->table.number_updates() - env->init_updates;
    size_t num_inserts = env->table.size() - env->init_size;
    size_t num_upserts = num_updates + num_inserts;
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
    std::cout << ", " << num_updates;
    std::cout << ", " << num_inserts;
    std::cout << ", " << num_finds_in;
    std::cout << ", " << num_finds_out;
    std::cout << ", " << elapsed_time/1000;
    std::cout << ", " << 
        (double)(num_upserts + num_finds_in + num_finds_out)/(elapsed_time/1000) << std::endl;
}

int main(int argc, char** argv) {
    const char* args[] = {"--power", "--thread-num", "--begin-load",
     "--end-load", "--seed", "--percent-inserts", "--percent-find-in"};
    size_t* arg_vars[] = {&power, &thread_num, &begin_load, &end_load, &seed, &percent_inserts, &percent_find_in};
    const char* arg_help[] = {"The power argument given to the hashtable during initialization",
                              "The number of threads to spawn for each type of operation",
                              "The load factor to fill the table up to before testing throughput",
                              "The maximum load factor to fill the table up to when testing throughput",
                              "The seed used by the random number generator",
                              "The percent of inserts in the workload",
                              "The percent of finds that are in the table"};
    const char* flags[] = {"--use-strings"};
    bool* flag_vars[] = {&use_strings};
    const char* flag_help[] = {"If set, the key type of the map will be std::string"};
    parse_flags(argc, argv, "A benchmark for upserts", args, arg_vars, arg_help,
                sizeof(args)/sizeof(const char*), flags, flag_vars, flag_help,
                sizeof(flags)/sizeof(const char*));

    if (begin_load >= 100) {
        std::cerr << "--begin-load must be between 0 and 99" << std::endl;
        exit(1);
    } else if (begin_load >= end_load) {
        std::cerr << "--end-load must be greater than --begin-load" << std::endl;
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
        auto *env = new InsertFindEnvironment<KeyType2>;
        InsertFindTest(env);
        delete env;
    } else {
        auto *env = new InsertFindEnvironment<KeyType>;
        InsertFindTest(env);
        delete env;
    }
}
