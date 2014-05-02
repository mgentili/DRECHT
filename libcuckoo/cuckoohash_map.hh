#ifndef _CUCKOOHASH_MAP_HH
#define _CUCKOOHASH_MAP_HH

#include <atomic>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <list>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>

#include "cuckoohash_config.h"
#include "cuckoohash_util.h"

/*! cuckoohash_map is the hash table class. */
template <class Key, class T, class Hash = std::hash<Key>, class Pred = std::equal_to<Key> >
class cuckoohash_map {

    // Structs and functions used internally

    typedef enum {
        ok = 0,
        failure = 1,
        failure_key_not_found = 2,
        failure_key_duplicated = 3,
        failure_space_not_enough = 4,
        failure_function_not_supported = 5,
        failure_table_full = 6,
        failure_under_expansion = 7,
        failure_key_moved = 8,
        failure_already_migrating_all = 9
    } cuckoo_status;

    /* This is a hazard pointer, used to indicate which version of the
     * TableInfo is currently being used in the thread. Since
     * cuckoohash_map operations can run simultaneously in different
     * threads, this variable is thread local. Note that this variable
     * can be safely shared between different cuckoohash_map
     * instances, since multiple operations cannot occur
     * simultaneously in one thread. The hazard pointer variable
     * points to a pointer inside a global list of pointers, that each
     * map checks before deleting any old TableInfo pointers. */
    static __thread void** hazard_pointer;

    /* A GlobalHazardPointerList stores a list of pointers that cannot
     * be deleted by an expansion thread. Each thread gets its own
     * node in the list, whose data pointer it can modify without
     * contention. */
    class GlobalHazardPointerList {
        std::list<void*> hp_;
        std::mutex lock_;
    public:
        /* new_hazard_pointer creates and returns a new hazard pointer for
         * a thread. */
        void** new_hazard_pointer() {
            lock_.lock();
            hp_.emplace_back(nullptr);
            void** ret = &hp_.back();
            lock_.unlock();
            return ret;
        }

        /* delete_unused scans the list of hazard pointers, deleting
         * any pointers in old_pointers that aren't in this list.
         * If it does delete a pointer in old_pointers, it deletes
         * that node from the list. */
        template <class Ptr>
        void delete_unused(std::list<Ptr*>& old_pointers) {
            lock_.lock();
            auto it = old_pointers.begin();
            while (it != old_pointers.end()) {
                std::cout << "Hazard pointer" << *it << std::endl;
                bool deleteable = true;
                for (auto hpit = hp_.cbegin(); hpit != hp_.cend(); hpit++) {
                    if (*hpit == *it) {
                        deleteable = false;
                        break;
                    }
                }
                if (deleteable) {
                    std::cout << "Deleting hazard pointer" << *it << std::endl;
                    LIBCUCKOO_DBG("deleting %p\n", *it);
                    delete *it;
                    it = old_pointers.erase(it);
                } else {
                    it++;
                }
            }
            lock_.unlock();
        }
    };

    // As long as the thread_local hazard_pointer is static, which
    // means each template instantiation of a cuckoohash_map class
    // gets its own per-thread hazard pointer, then each template
    // instantiation of a cuckoohash_map class can get its own
    // global_hazard_pointers list, since different template
    // instantiations won't interfere with each other.
    static GlobalHazardPointerList global_hazard_pointers;

    /* check_hazard_pointer should be called before any public method
     * that loads a table snapshot. It checks that the thread local
     * hazard pointer pointer is not null, and gets a new pointer if
     * it is null. */
    static inline void check_hazard_pointer() {
        if (hazard_pointer == nullptr) {
            hazard_pointer = global_hazard_pointers.new_hazard_pointer();
        }
    }

    /* Once a function is finished with a version of the table, it
     * calls unset_hazard_pointer so that the pointer can be freed if
     * it needs to. */
    static inline void unset_hazard_pointer() {
        *hazard_pointer = nullptr;
    }

    /* counterid stores the per-thread counter index of each
     * thread. */
    static __thread int counterid;

    /* check_counterid checks if the counterid has already been
     * determined. If not, it assigns a counterid to the current
     * thread by picking a random core. This should be called at the
     * beginning of any function that changes the number of elements
     * in the table. */
    static inline void check_counterid() {
        if (counterid < 0) {
            counterid = rand() % kNumCores;
        }
    }

    /* reserve_calc takes in a parameter specifying a certain number
     * of slots for a table and returns the smallest hashpower that
     * will hold n elements. */
    size_t reserve_calc(size_t n) {
        size_t new_hashpower = (size_t)ceil(log2((double)n / (double)SLOT_PER_BUCKET));
        assert(n <= hashsize(new_hashpower) * SLOT_PER_BUCKET);
        return new_hashpower;
    }

public:
    //! key_type is the type of keys.
    typedef Key               key_type;
    //! value_type is the type of key-value pairs.
    typedef std::pair<Key, T> value_type;
    //! mapped_type is the type of values.
    typedef T                 mapped_type;
    //! hasher is the type of the hash function.
    typedef Hash              hasher;
    //! key_equal is the type of the equality predicate.
    typedef Pred              key_equal;

    /*! The constructor creates a new hash table with enough space for
     * \p n elements. If the constructor fails, it will throw an
     * exception. */
    explicit cuckoohash_map(size_t n = DEFAULT_SIZE) {
        cuckoo_init(reserve_calc(n));
    }

    /*! The destructor deletes any remaining table pointers managed by
     * the hash table, also destroying all remaining elements in the
     * table. */
    ~cuckoohash_map() {
        TableInfo *ti = table_info.load();
        if (ti != nullptr) {
            delete ti;
        }
        for (auto it = old_table_infos.cbegin(); it != old_table_infos.cend(); it++) {
            delete *it;
        }
    }

    /*! TODO: clear removes all the elements in the hash table, calling
     *  their destructors. */
    void clear() {
        check_hazard_pointer();
        expansion_lock.lock();
        TableInfo *ti = snapshot_and_lock_write_all();
        assert(ti != nullptr);
        cuckoo_clear(ti);
        unlock_write_all(ti);
        expansion_lock.unlock();
        unset_hazard_pointer();
    }

    void force_migration() {
        check_hazard_pointer();
        TableInfo *ti_old, *ti_new;
        snapshot_old(ti_old);
        snapshot_new(ti_new);
        try_migrate_all(ti_old, ti_new, 1);
    }

    /*! size returns the number of items currently in the hash table.
     * Since it doesn't lock the table, elements can be inserted
     * during the computation, so the result may not necessarily be
     * exact. */
    size_t size() {
        check_hazard_pointer();
        TableInfo *ti_old, *ti_new;
        snapshot_old(ti_old);
        snapshot_new(ti_new);
        size_t s1 = cuckoo_size(ti_old);
        size_t s2 = cuckoo_size(ti_new);
        unset_hazard_pointer();
        return s1+s2;
    }

    size_t number_retries() {
        check_hazard_pointer();
        const TableInfo *ti = snapshot_table_nolock();
        size_t retries = 0;
        for (size_t i = 0; i < ti->num_retries.size(); i++) {
            //std::cout << "Num retries" << ti->num_retries[i].num.load();
            retries += ti->num_retries[i].num.load();
        }
        unset_hazard_pointer();
        return retries;
    }

    /*! empty returns true if the table is empty. */
    bool empty() {
        return size() == 0;
    }

    /* undergoing_expansion returns true if there are currently two tables in use */
    bool undergoing_expansion() {
        TableInfo* ti_new;
        snapshot_new(ti_new);
        return (ti_new == nullptr);
    }
    /*! hashpower returns the hashpower of the table, which is log<SUB>2</SUB>(the
     * number of buckets). */
    size_t hashpower() {
        check_hazard_pointer();
        TableInfo* ti = snapshot_table_nolock();
        const size_t hashpower = ti->hashpower_;
        unset_hazard_pointer();
        return hashpower;
    }

    /*! bucket_count returns the number of buckets in the table. */
    size_t bucket_count() {
        check_hazard_pointer();
        TableInfo *ti = snapshot_table_nolock();
        size_t buckets = hashsize(ti->hashpower_);
        unset_hazard_pointer();
        return buckets;
    }

    /*! load_factor returns the ratio of the number of items in the
     * table to the total number of available slots in the table. */
    float load_factor() {
        check_hazard_pointer();
        const TableInfo *ti = snapshot_table_nolock();
        const float lf = cuckoo_loadfactor(ti);
        unset_hazard_pointer();
        return lf;
    }

    /*! find searches through the table for \p key, and stores
     * the associated value it finds in \p val. */
    bool find(const key_type& key, mapped_type& val) {
        check_hazard_pointer();
        //check_counterid(); //TODO: Remove when not counting reads
        size_t hv = hashed_key(key);
        
        TableInfo *ti;
        cuckoo_status res;

        snapshot_old(ti);
        std::cout << "Starting find of key" << key << std::endl;
        res = find_one(ti, hv, key, val);

        // couldn't find key in bucket, and one of the buckets was moved to new table
        if (res == failure_key_moved) {
            snapshot_new(ti);
            assert(ti != nullptr);
            res = find_one(ti, hv, key, val);
            assert(res == ok || res == failure_key_not_found);
        }

        unset_hazard_pointer(); //only need to keep track of ptr to old table b/c new never will be deleted until in old pos
        return (res == ok);
    }

    /*! This version of find does the same thing as the two-argument
     * version, except it returns the value it finds, throwing an \p
     * std::out_of_range exception if the key isn't in the table. */
    mapped_type find(const key_type& key) {
        mapped_type val;
        bool done = find(key, val);

        if (done) {
            return val;
        } else {
            throw std::out_of_range("key not found in table");
        }
    }

    /*! insert puts the given key-value pair into the table. It first
     * checks that \p key isn't already in the table, since the table
     * doesn't support duplicate keys. If the table is out of space,
     * insert will automatically expand until it can succeed. Note
     * that expansion can throw an exception, which insert will
     * propagate. */
    bool insert(const key_type& key, const mapped_type& val) {
        check_hazard_pointer();
        check_counterid();
        size_t hv = hashed_key(key);
        TableInfo *ti_old, *ti_new;
        size_t i1_o, i2_o, i1_n, i2_n;
        cuckoo_status res;

        snapshot_old(ti_old);
        std::cout << "Starting insert of key" << key << std::endl;
        // lock and don't unlock
        res = insert_one(ti_old, hv, key, val, i1_o, i2_o);

        snapshot_new(ti_new); //need to check this to see if expansion in progress
        std::cout << "In insert with old pointer" << ti_old << "new pointer" << ti_new << std::endl;
        // expansion in progress, so moving buckets over
        if (ti_new != ti_old && ti_new != nullptr) {
            std::cout << "trying to migrate buckets for insert!" << std::endl;
            //try to move both buckets for now
            //migrate bucket shouldn't ever fail (e.g. need to resize)
            //because load factor should be sufficiently low
            //assumes both old buckets are locked
            try_migrate_bucket(ti_old, ti_new, i1_o);
            try_migrate_bucket(ti_old, ti_new, i2_o);

            // if sufficient number of buckets moved, start a thread that starts from beginning of table to end
            // trying to move each bucket
            std::cout << "Percent moved buckets is" << double(count_migrated_buckets(ti_old))/hashsize(ti_old->hashpower_) << std::endl;
            if (double(count_migrated_buckets(ti_old))/hashsize(ti_old->hashpower_) > MIGRATE_THRESHOLD) {
                try_migrate_all(ti_old, ti_new, 1);
            }

        }

        // This is triggered only if we couldn't find the key in either
        // old table bucket and one of the buckets was moved
        // We hold both old bucket locks at the start
        if (res == failure_key_moved) {
            snapshot_new(ti_new);
            assert(ti_new != nullptr);
            std::cout << "result is failure_key_moved! with new pointer" << ti_new << std::endl;
            res = insert_one(ti_new, hv, key, val, i1_n, i2_n);
            assert(res == failure_key_duplicated || res == ok);
            unlock_write_two(ti_new, i1_n, i2_n); //TODO: When are we unlocking?
        }
        std::cout << "Ending insert of key" << key << std::endl;
        unlock_write_two(ti_old, i1_o, i2_o); //TODO: When are we unlocking?
        std::cout << "Unsetting hazard pointer" << key << std::endl;
        unset_hazard_pointer();
        std::cout << "Returning!" << key << std::endl;
        return (res == ok);
    }
    
    /*! erase removes \p key and it's associated value from the table,
     * calling their destructors. If \p key is not there, it returns
     * false. */
    bool erase(const key_type& key) {
        check_hazard_pointer();
        check_counterid();
        size_t hv = hashed_key(key);
        TableInfo *ti_old, *ti_new;
        size_t i1_o, i2_o, i1_n, i2_n;
        cuckoo_status res;

        snapshot_old(ti_old);
        // lock and don't unlock
        res = delete_one(ti_old, hv, key, i1_o, i2_o);
        snapshot_new(ti_new);

        if (ti_new != ti_old && ti_new != nullptr) {
            std::cout << "We're migrating buckets in erase!" << std::endl;
            try_migrate_bucket(ti_old, ti_new, i1_o);
            try_migrate_bucket(ti_old, ti_new, i2_o);
        }

        if (res == failure_key_moved) {
            snapshot_new(ti_new);
            assert(ti_new != nullptr);
            res = delete_one(ti_new, hv, key, i1_n, i2_n);
            assert(res == ok || res == failure_key_not_found);
            unlock_write_two(ti_new, i1_n, i2_n);
        }
        unlock_write_two(ti_old, i1_o, i2_o);
        unset_hazard_pointer();
        return (res == ok);
    }

    /*! rehash will size the table using a hashpower of \p n. Note
     * that the number of buckets in the table will be 2<SUP>\p
     * n</SUP> after expansion, so the table will have 2<SUP>\p
     * n</SUP> &times; \ref SLOT_PER_BUCKET slots to store items in.
     * If \p n is not larger than the current hashpower, then the
     * function does nothing. It returns true if the table expansion
     * succeeded, and false otherwise. rehash can throw an exception
     * if the expansion fails to allocate enough memory for the larger
     * table. */
    bool rehash(size_t n) {
        check_hazard_pointer();
        TableInfo* ti = snapshot_table_nolock();
        if (n <= ti->hashpower_) {
            return false;
        }
        const cuckoo_status st = cuckoo_expand_simple(n);
        unset_hazard_pointer();
        return (st == ok);
    }

    /*! reserve will size the table to have enough slots for at least
     * \p n elements. If the table can already hold that many
     * elements, the function has no effect. Otherwise, the function
     * will expand the table to a hashpower sufficient to hold \p n
     * elements. It will return true if there was an expansion, and
     * false otherwise. reserve can throw an exception if the
     * expansion fails to allocate enough memory for the larger
     * table. */
    bool reserve(size_t n) {
        check_hazard_pointer();
        TableInfo* ti = snapshot_table_nolock();
        if (n <= hashsize(ti->hashpower_) * SLOT_PER_BUCKET) {
            return false;
        }
        const cuckoo_status st = cuckoo_expand_simple(reserve_calc(n));
        unset_hazard_pointer();
        return (st == ok);
    }

    /*! hash_function returns the hash function object used by the
     * table. */
    hasher hash_function() {
        return hashfn;
    }

    /*! key_eq returns the equality predicate object used by the
     * table. */
    key_equal key_eq() {
        return eqfn;
    }

private:

    /* cacheint is a cache-aligned atomic integer type. */
    struct cacheint {
        std::atomic<size_t> num;
        cacheint() {}
        cacheint(cacheint&& x) {
            num.store(x.num.load());
        }
    } __attribute__((aligned(64)));

    /* The Bucket type holds SLOT_PER_BUCKET keys and values, and a
     * occupied bitset, which indicates whether the slot at the given
     * bit index is in the table or not. It allows constructing and
     * destroying key-value pairs separate from allocating and
     * deallocating the memory. */
    struct Bucket {
        cacheint version;
        std::bitset<SLOT_PER_BUCKET> occupied;
        key_type keys[SLOT_PER_BUCKET];
        mapped_type vals[SLOT_PER_BUCKET];
        bool need_check_alternate;
        bool hasmigrated;

        void setKV(size_t pos, const key_type& k, const mapped_type& v) {
            occupied.set(pos);
            new (keys+pos) key_type(k);
            new (vals+pos) mapped_type(v);
        }

        void eraseKV(size_t pos) {
            occupied.reset(pos);
            (keys+pos)->~key_type();
            (vals+pos)->~mapped_type();
        }

        void clear() {
            for (size_t i = 0; i < SLOT_PER_BUCKET; i++) {
                if (occupied[i]) {
                    eraseKV(i);
                }
            }
            version.num.store(0);
        }
    };

    /* TableInfo contains the entire state of the hashtable. We
     * allocate one TableInfo pointer per hash table and store all of
     * the table memory in it, so that all the data can be atomically
     * swapped during expansion. */
    struct TableInfo {
        // 2**hashpower is the number of buckets
        size_t hashpower_;

        // unique pointer to the array of buckets
        Bucket* buckets_;

        // per-core counters for the number of inserts and deletes
        std::vector<cacheint> num_inserts;
        std::vector<cacheint> num_deletes;
        std::vector<cacheint> num_migrated_buckets;
        // counter for the number of find retries
        std::vector<cacheint> num_retries;

        /* The constructor allocates the memory for the table. For
         * buckets, it uses the bucket_allocator, so that we can free
         * memory independently of calling its destructor. It
         * allocates one cacheint for each core in num_inserts and
         * num_deletes. */
        TableInfo(const size_t hashtable_init) {
            buckets_ = nullptr;
            try {
                hashpower_ = hashtable_init;
                buckets_ = bucket_allocator.allocate(hashsize(hashpower_));
                for (size_t i = 0; i < hashsize(hashpower_); i++) {
                    bucket_allocator.construct(buckets_+i);
                }
                num_inserts.resize(kNumCores);
                num_deletes.resize(kNumCores);
                num_migrated_buckets.resize(kNumCores);
                num_retries.resize(kNumCores);

            } catch (const std::bad_alloc&) {
                if (buckets_ != nullptr) {
                    bucket_allocator.deallocate(buckets_, hashsize(hashpower_));
                }
                throw;
            }
        }

        ~TableInfo() {
            for (size_t i = 0; i < hashsize(hashpower_); i++) {
                buckets_[i].clear();
            }
            bucket_allocator.deallocate(buckets_, hashsize(hashpower_));
        }

    };
    std::atomic<TableInfo*> table_info;
    std::atomic<TableInfo*> new_table_info;

    /* old_table_infos holds pointers to old TableInfos that were
     * replaced during expansion. This keeps the memory alive for any
     * leftover operations, until they are deleted by the global
     * hazard pointer manager. */
    std::list<TableInfo*> old_table_infos;

    static hasher hashfn;
    static key_equal eqfn;
    static std::allocator<Bucket> bucket_allocator;


    /* find_one searches a specific table instance for the value corresponding to a given hash value 
     * It doesn't take any locks */
    cuckoo_status find_one(const TableInfo *ti, size_t hv, const key_type& key, mapped_type& val) {
        size_t i1, i2, v1_i, v2_i, v1_f, v2_f;
        cuckoo_status st;
        i1 = index_hash(ti, hv);
        i2 = alt_index(ti, hv, i1);
        do {
            get_version_two(ti, i1, i2, v1_i, v2_i);
            st = cuckoo_find(key, val, hv, ti, i1, i2);
            get_version_two(ti, i1, i2, v1_f, v2_f);
            //ti->num_retries[counterid].num.fetch_add(1, std::memory_order_relaxed);
        } while(!check_version_two(v1_i, v2_i, v1_f, v2_f));

        return st;
    }

    /* insert_one tries to insert a key-value pair into a specific table instance.
     * It will return a failure only if the key is already in the table.
     * Regardless, it starts with the buckets unlocked and ends with buckets locked 
     */
    cuckoo_status insert_one(TableInfo *ti, size_t hv, const key_type& key,
                            const mapped_type& val, size_t& i1, size_t& i2) {
        i1 = index_hash(ti, hv);
        i2 = alt_index(ti, hv, i1);
        //std::cout << "In insert loop for key" << key << std::endl;
        lock_write_two(ti, i1, i2);
        //TODO: Don't unlock '
        cuckoo_status res = cuckoo_insert(key, val, hv, ti, i1, i2);
        //std::cout << "We finished an insert loop for key" << key << std::endl;
        lock_write_two(ti, i1, i2);
        //TODO: If unsuccessful because of expansion, make sure that nothing inserted in between
        return res;
        //return cuckoo_insert(key, val, hv, ti, i1, i2);
    }

    cuckoo_status delete_one(TableInfo *ti, size_t hv, const key_type& key,
                             size_t& i1, size_t& i2) {
        i1 = index_hash(ti, hv);
        i2 = alt_index(ti, hv, i1);
        lock_write_two(ti, i1, i2);
        cuckoo_status res = cuckoo_delete(key, hv, ti, i1, i2);

        return res;
    }

    // assumes bucket in old table is locked the whole time
    // tries to migrate bucket, returning true on success, false on failure
    bool try_migrate_bucket(TableInfo* ti_old, TableInfo* ti_new, size_t old_bucket) {
        if (ti_old->buckets_[old_bucket].hasmigrated) {
            //std::cout << "Already migrated bucket " << old_bucket << std::endl;
            return false;
        }
        size_t i1, i2, hv;
        key_type key;
        mapped_type val;
        cuckoo_status res;
        (void) res;
        assert(!ti_old->buckets_[old_bucket].hasmigrated);
        for (size_t i = 0; i < SLOT_PER_BUCKET; i++) {
            if (!ti_old->buckets_[old_bucket].occupied[i]) {
                continue;
            }
            key = ti_old->buckets_[old_bucket].keys[i];
            val = ti_old->buckets_[old_bucket].vals[i];

            hv = hashed_key(key);
            i1 = index_hash(ti_new, hv);
            i2 = alt_index(ti_new, hv, i1);
            //std::cout << "Moving element key: " << key << "val" << val << "to pos" << i1 << "or" << i2 << std::endl;
            lock_write_two(ti_new, i1, i2);
            res = cuckoo_insert(key, val, hv, ti_new, i1, i2); // TODO: cuckoo_insert unlocks locks for now
            assert(res == ok); // cannot have inserted into new table before
            //unlock_write_two(ti_new, i1, i2);

            //to keep track of current number of elements in table
            ti_old->num_deletes[counterid].num.fetch_add(1, std::memory_order_relaxed);
        }

        ti_old->buckets_[old_bucket].hasmigrated = true;
        ti_old->num_migrated_buckets[counterid].num.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

    void migrate_bucket_range(TableInfo* ti_old, TableInfo* ti_new, size_t begin, size_t end) {
        check_hazard_pointer();
        check_counterid();
        assert(ti_old == table_info.load());
        assert(ti_new == new_table_info.load());
        for (size_t i = begin; i < end; i++) {
            lock_write(ti_old, i);
            std::cout << "In migrate_bucket_range. Trying to migrate bucket " << i << std::endl;
            try_migrate_bucket(ti_old, ti_new, i);
            unlock_write(ti_old, i);
        }

        std::cout << "migrate_bucket_range done!" << std::endl;

        if(count_migrated_buckets(ti_old) == hashsize(ti_old->hashpower_)) {
            cuckoo_status res = cuckoo_expand_end();
            if (res != ok) {
                std::cout << "Someone tried to call cuckoo_expand_end before me" << std::endl;
                return;
            }
            migrate_all_lock.unlock();
        }


    }

    cuckoo_status try_migrate_all(TableInfo* ti_old, TableInfo* ti_new, size_t threadnum) {
        if(!migrate_all_lock.try_lock()) {
            //std::cout << "Already migrating all" << std::endl;
            return failure_already_migrating_all;
        }

        std::cout << "In try_migrate_all and starting new threads" << std::endl;
        const size_t buckets_per_thread = hashsize(ti_old->hashpower_) / threadnum;
        std::vector<std::thread> migrate_threads(threadnum);
        for (size_t i = 0; i < threadnum-1; i++) {
            migrate_threads[i] = std::thread( &cuckoohash_map<Key, T, Hash>::migrate_bucket_range,
                this, ti_old, ti_new, 
                i*buckets_per_thread, (i+1)*buckets_per_thread);
        }
        // remaining buckets
        migrate_threads[threadnum-1] = std::thread(&cuckoohash_map<Key, T, Hash>::migrate_bucket_range,
            this, ti_old, ti_new, 
            (threadnum-1)*buckets_per_thread, hashsize(ti_old->hashpower_));

        for (size_t i = 0; i < threadnum; i++) {
            migrate_threads[i].detach();
        }
        return ok;
    }

    /* get_version gets the version for a given bucket index */
    static inline void get_version(const TableInfo *ti, const size_t i, size_t& v) {
        v = ti->buckets_[i].version.num.load();
    }

    static inline void get_version_two(const TableInfo *ti, const size_t i1, const size_t i2, 
                                       size_t& v1, size_t& v2) {
        v1 = ti->buckets_[i1].version.num.load();
        v2 = ti->buckets_[i2].version.num.load();
    }

    

    /* check_version makes sure that the final version is the same as the initial version, and 
     * that the initial version is not dirty
     */
    static inline bool check_version(const size_t v_i, const size_t v_f) {
        return (v_i==v_f && (v_i & W) == 0);
    }

    static inline bool check_version_two(const size_t v1_i, const size_t v2_i,
                                         const size_t v1_f, const size_t v2_f) {
        return check_version(v1_i,v1_f) && check_version(v2_i, v2_f);
    }

    /* lock locks a bucket based on lock_level privileges
     */
    static inline void lock(const TableInfo *ti, const size_t i, size_t lock_level) {
        size_t v;
        get_version(ti, i, v);
        /*if ((v & lock_level) != 0) {
            std::cout << "Starting lock of" << i << "with version" << v << "and lock level" << lock_level << std::endl;
        }*/
        while(true) {
            if( (v & lock_level) != 0) {
                get_version(ti, i, v);
                continue;
            }

            if(ti->buckets_[i].version.num.compare_exchange_weak(
                v, (size_t) (v | lock_level) , std::memory_order_release, std::memory_order_relaxed)) {
                //get_version(ti, i, v);
                //std::cout << "Ending lock of" << i << "with version" << v << "and lock level" << lock_level << std::endl;
                return;
            }
        }
    }

    static inline void unlock(const TableInfo *ti, const size_t i, const size_t lock_level) {
        /*if ((~ti->buckets_[i].version.num.load() & lock_level) != 0) {
            std::cout << "Using table with size" << ti->hashpower_ << std::endl;
            std::cout << "Starting unlock of" << i << "with version" << ti->buckets_[i].version.num.load() << "and lock level" << lock_level << std::endl;
        }*/
        ti->buckets_[i].version.num += (1<<2) - lock_level;
        //std::cout << "Ending unlock of" << i << "with version"<< ti->buckets_[i].version.num.load() << "and lock level" << lock_level << std::endl;
    }

    static inline void try_lock_write(const TableInfo *ti, const size_t i) {

    }
    /* lock_paths locks out concurrent writers, but allows concurrent readers
     */
    static inline void lock_path(const TableInfo *ti, const size_t i) {
        lock(ti, i, P);
    }

    static inline void unlock_path(const TableInfo *ti, const size_t i) {
        unlock(ti, i, P);
    }
    /* lock_write locks out concurrent writers and readers
     */
    static inline void lock_write(const TableInfo *ti, const size_t i) {
        lock(ti, i, P | W);
    }

    static inline void unlock_write(const TableInfo *ti, const size_t i) {
        unlock(ti, i, P | W);
    }
    
    static inline void lock_path_two(const TableInfo *ti, size_t i1, size_t i2) {
        if (i1 < i2) {
            lock_path(ti, i1);
            lock_path(ti, i2);
        } else if (i2 < i1) {
            lock_path(ti, i2);
            lock_path(ti, i1);
        } else {
            lock_path(ti, i1);
        }
    }

    static inline void path_unlock_two(const TableInfo *ti, size_t i1, size_t i2) {
        path_unlock(ti, i1);
        if (i1 != i2) {
            path_unlock(ti, i2);
        }
    }

    static inline void lock_write_two(const TableInfo *ti, size_t i1, size_t i2) {
        if (i1 < i2) {
            lock_write(ti, i1);
            lock_write(ti, i2);
        } else if (i2 < i1) {
            lock_write(ti, i2);
            lock_write(ti, i1);
        } else {
            lock_write(ti, i1);
        }
    }

    static inline void unlock_write_two(const TableInfo *ti, size_t i1, size_t i2) {
        unlock_write(ti, i1);
        if (i1 != i2) {
            unlock_write(ti, i2);
        }
    }

    static inline void upgrade_to_writelock(const TableInfo *ti, size_t i) {
        ti->buckets_[i].versions.num |= W;
    }

    static inline void upgrade_to_writelock_two(const TableInfo *ti, size_t i1, size_t i2) {
        upgrade_to_writelock(ti, i1);
        upgrade_to_writelock(ti, i2);
    }

    static inline void get_version_three(const TableInfo *ti, const size_t i1, const size_t i2,
                                         const size_t i3, size_t& v1, size_t& v2, size_t& v3) {
        v1 = ti->buckets_[i1].version.num.load();
        v2 = ti->buckets_[i2].version.num.load();
        v3 = ti->buckets_[i3].version.num.load();
    }

    static inline void lock_write_three(const TableInfo *ti, size_t i1,
                                  size_t i2, size_t i3) {
        // If any are the same, we just run lock_two
        if (i1 == i2) {
            lock_write_two(ti, i1, i3);
        } else if (i2 == i3) {
            lock_write_two(ti, i1, i3);
        } else if (i1 == i3) {
            lock_write_two(ti, i1, i2);
        } else {
            if (i1 < i2) {
                if (i2 < i3) {
                    lock_write(ti, i1);
                    lock_write(ti, i2);
                    lock_write(ti, i3);
                } else if (i1 < i3) {
                    lock_write(ti, i1);
                    lock_write(ti, i3);
                    lock_write(ti, i2);
                } else {
                    lock_write(ti, i3);
                    lock_write(ti, i1);
                    lock_write(ti, i2);
                }
            } else if (i2 < i3) {
                if (i1 < i3) {
                    lock_write(ti, i2);
                    lock_write(ti, i1);
                    lock_write(ti, i3);
                } else {
                    lock_write(ti, i2);
                    lock_write(ti, i3);
                    lock_write(ti, i1);
                }
            } else {
                lock_write(ti, i3);
                lock_write(ti, i2);
                lock_write(ti, i1);
            }
        }
    }

    /* unlock_three unlocks the three given buckets */
    static inline void unlock_write_three(const TableInfo *ti, size_t i1,
                                    size_t i2, size_t i3) {
        unlock_write(ti, i1);
        if (i2 != i1) {
            unlock_write(ti, i2);
        }
        if (i3 != i1 && i3 != i2) {
            unlock_write(ti, i3);
        }
    }

    void snapshot_old(TableInfo*& ti) {
    TryAcquire:
        ti = table_info.load();
        *hazard_pointer = static_cast<void*>(ti); //to keep track of which tables are currently in use
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
    }

    // don't need to because it can never be garbage collected before it becomes an "old" table pointer 
    void snapshot_new(TableInfo*& ti) {
        ti = new_table_info.load();
    }

    /* snapshot_and_get_version_two loads the table_info pointer and gets the versions
     * of the buckets associated with the given hash value. It stores the
     * table_info and the two bucket versions in reference variables.
     * Since the positions of the bucket versions depends on the number
     * of buckets in the table, the table_info pointer needs to be
     * grabbed first. */
    void snapshot_and_get_version_two(const size_t hv, TableInfo*& ti,
                               size_t& i1, size_t& i2, size_t& v1, size_t& v2) {
    TryAcquire:
        ti = table_info.load();
        *hazard_pointer = static_cast<void*>(ti);
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
        i1 = index_hash(ti, hv);
        i2 = alt_index(ti, hv, i1);
        get_version_two(ti, i1, i2, v1, v2);
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
    }

    void snapshot_and_lock_path_two(const size_t hv, TableInfo*& ti,
                               size_t& i1, size_t& i2) {
    TryAcquire:
        //std::cout << "In snapshot and inc version two" << std::endl;
        ti = table_info.load();
        *hazard_pointer = static_cast<void*>(ti);
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
        i1 = index_hash(ti, hv);
        i2 = alt_index(ti, hv, i1);
        lock_path_two(ti, i1, i2);
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
    }

    void snapshot_and_lock_write_two(const size_t hv, TableInfo*& ti,
                               size_t& i1, size_t& i2) {
    TryAcquire:
        //std::cout << "In snapshot and inc version two" << std::endl;
        ti = table_info.load();
        *hazard_pointer = static_cast<void*>(ti);
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
        i1 = index_hash(ti, hv);
        i2 = alt_index(ti, hv, i1);
        lock_write_two(ti, i1, i2);
        if (ti != table_info.load()) {
            goto TryAcquire;
        }
    }

    /* snapshot_and_lock_write_all increases the version of every counter to 1mod2
     * so that no inserts or finds will be able to be performed (effectively
     * locking all the buckets)
     */
    TableInfo *snapshot_and_lock_write_all() {
        assert(!expansion_lock.try_lock());
        TableInfo *ti = table_info.load();
        *hazard_pointer = static_cast<void*>(ti);
        for (size_t i = 0; i < hashsize(ti->hashpower_); i++) {
            lock_write(ti, i);
        }
        return ti;
    }

    /* unlock_write_all increases the version of every counter (from 1mod2 to 0mod2)
     * effectively releases all the "locks" on the buckets */
    inline void unlock_write_all(TableInfo *ti) {
        for (size_t i = 0; i < hashsize(ti->hashpower_); i++) {
            unlock_write(ti, i);
        }
    }

    /* snapshot_table_nolock loads the table info pointer and sets the
     * hazard pointer, whithout locking anything. There is a
     * possibility that after loading a snapshot and setting the
     * hazard pointer, an expansion runs and create a new version of
     * the table, leaving the old one for deletion. To deal with that,
     * we check that the table_info we loaded is the same as the
     * current one, and if it isn't, we try again. Whenever we check
     * if (ti != table_info.load()) after setting the hazard pointer,
     * there is an ABA issue, where the address of the new table_info
     * equals the address of a previously deleted one, however it
     * doesn't matter, since we would still be looking at the most
     * recent table_info in that case. */
    TableInfo* snapshot_table_nolock() {
        TableInfo *ti;
    TryAcquire:
        ti = table_info.load();
        *hazard_pointer = static_cast<void*>(ti);
        if (ti != table_info.load()) {
            std::cout << "Table info changed!" << std::endl;
            goto TryAcquire;
        }
        return ti;
    }

    static const size_t W = 1 << 1;
    static const size_t P = 1 << 1;

    // key size in bytes
    static const size_t kKeySize = sizeof(key_type);

    // value size in bytes
    static const size_t kValueSize = sizeof(mapped_type);

    // size of a bucket in bytes
    static const size_t kBucketSize = sizeof(Bucket);

    // number of cores on the machine
    static const size_t kNumCores;

    // The maximum number of cuckoo operations per insert. This must
    // be less than or equal to SLOT_PER_BUCKET^(MAX_BFS_DEPTH+1)
    static const size_t MAX_CUCKOO_COUNT = 500;

    // The maximum depth of a BFS path
    static const size_t MAX_BFS_DEPTH = 4;

    // the % of moved buckets above which migrate_all is called
    static constexpr double MIGRATE_THRESHOLD = 0.8;

    /* hashsize returns the number of buckets corresponding to a given
     * hashpower. */
    static inline size_t hashsize(const size_t hashpower) {
        return 1U << hashpower;
    }

    /* hashmask returns the bitmask for the buckets array
     * corresponding to a given hashpower. */
    static inline size_t hashmask(const size_t hashpower) {
        return hashsize(hashpower) - 1;
    }

    /* hashed_key hashes the given key. */
    static inline size_t hashed_key(const key_type &key) {
        return hashfn(key);
    }

    /* index_hash returns the first possible bucket that the given
     * hashed key could be. */
    static inline size_t index_hash(const TableInfo *ti, const size_t hv) {
        return hv & hashmask(ti->hashpower_);
    }

    /* alt_index returns the other possible bucket that the given
     * hashed key could be. It takes the first possible bucket as a
     * parameter. Note that this function will return the first
     * possible bucket if index is the second possible bucket, so
     * alt_index(ti, hv, alt_index(ti, hv, index_hash(ti, hv))) ==
     * index_hash(ti, hv). */
    static inline size_t alt_index(const TableInfo *ti,
                                   const size_t hv,
                                   const size_t index) {
        // ensure tag is nonzero for the multiply
        const size_t tag = (hv >> ti->hashpower_) + 1;
        /* 0x5bd1e995 is the hash constant from MurmurHash2 */
        return (index ^ (tag * 0x5bd1e995)) & hashmask(ti->hashpower_);
    }

    /* CuckooRecord holds one position in a cuckoo path. */
    typedef struct  {
        size_t bucket;
        size_t slot;
        key_type key;
    }  CuckooRecord;

    /* b_slot holds the information for a BFS path through the
     * table */
    struct b_slot {
        // The bucket of the last item in the path
        size_t bucket;
        // a compressed representation of the slots for each of the
        // buckets in the path.
        size_t pathcode;
        // static_assert(pow(SLOT_PER_BUCKET, MAX_BFS_DEPTH+1) <
        //               std::numeric_limits<decltype(pathcode)>::max(),
        //               "pathcode may not be large enough to encode a cuckoo path");
        // The 0-indexed position in the cuckoo path this slot
        // occupies
        int depth;
        b_slot() {}
        b_slot(const size_t b, const size_t p, const int d)
            : bucket(b), pathcode(p), depth(d) {}
    } __attribute__((__packed__));

    /* b_queue is the queue used to store b_slots for BFS cuckoo
     * hashing. */
    class b_queue {
        b_slot slots[MAX_CUCKOO_COUNT+1];
        size_t first;
        size_t last;

    public:
        b_queue() : first(0), last(0) {}


        void enqueue(b_slot x) {
            slots[last] = x;
            last = (last == MAX_CUCKOO_COUNT) ? 0 : last+1;
            assert(last != first);
        }

        b_slot dequeue() {
            assert(first != last);
            b_slot& x = slots[first];
            first = (first == MAX_CUCKOO_COUNT) ? 0 : first+1;
            return x;
        }

        bool not_full() {
            const size_t next = (last == MAX_CUCKOO_COUNT) ? 0 : last+1;
            return next != first;
        }
    } __attribute__((__packed__));

    /* slot_search searches for a cuckoo path using breadth-first
       search. It starts with the i1 and i2 buckets, and, until it finds
       a bucket with an empty slot, adds each slot of the bucket in the
       b_slot. If the queue runs out of space, it fails. */
    static b_slot slot_search(const TableInfo *ti, const size_t i1,
                              const size_t i2) {
        b_queue q;
        // The initial pathcode informs cuckoopath_search which bucket
        // the path starts on
        q.enqueue(b_slot(i1, 0, 0));
        q.enqueue(b_slot(i2, 1, 0));
        while (q.not_full()) {
            b_slot x = q.dequeue();
            // Picks a random slot to start from
            for (size_t slot = 0; slot < SLOT_PER_BUCKET && q.not_full(); slot++) {
                lock_path(ti, x.bucket);

                if (ti->buckets_[x.bucket].hasmigrated) {
                    unlock_path(ti, x.bucket);
                    return b_slot(0, 0, -1);
                }
                if (!ti->buckets_[x.bucket].occupied[slot]) {
                    // We can terminate the search here
                    x.pathcode = x.pathcode * SLOT_PER_BUCKET + slot;
                    unlock_path(ti, x.bucket);
                    return x;
                }
                // Create a new b_slot item, that represents the
                // bucket we would look at after searching x.bucket
                // for empty slots.
                const size_t hv = hashed_key(ti->buckets_[x.bucket].keys[slot]);
                unlock_path(ti, x.bucket);
                b_slot y(alt_index(ti, hv, x.bucket),
                         x.pathcode * SLOT_PER_BUCKET + slot, x.depth+1);

                // Check if any of the slots in the prospective bucket
                // are empty, and, if so, return that b_slot. We lock
                // the bucket so that no changes occur while
                // iterating.
                lock_path(ti, y.bucket);
                if (ti->buckets_[y.bucket].hasmigrated) {
                    unlock_path(ti, y.bucket);
                    return b_slot(0, 0, -1);
                }
                for (size_t j = 0; j < SLOT_PER_BUCKET; j++) {
                    if (!ti->buckets_[y.bucket].occupied.test(j)) {
                        y.pathcode = y.pathcode * SLOT_PER_BUCKET + j;
                        unlock_path(ti, y.bucket);
                        return y;
                    }
                }
                unlock_path(ti, y.bucket);

                // No empty slots were found, so we push this onto the
                // queue
                if (y.depth != static_cast<int>(MAX_BFS_DEPTH)) {
                    q.enqueue(y);
                }
            }
        }
        // We didn't find a short-enough cuckoo path, so the queue ran
        // out of space. Return a failure value.
        return b_slot(0, 0, -1);
    }

    /* cuckoopath_search finds a cuckoo path from one of the starting
     * buckets to an empty slot in another bucket. It returns the
     * depth of the discovered cuckoo path on success, and -1 on
     * failure. Since it doesn't take locks on the buckets it
     * searches, the data can change between this function and
     * cuckoopath_move. Thus cuckoopath_move checks that the data
     * matches the cuckoo path before changing it. */
    static int cuckoopath_search(const TableInfo *ti, CuckooRecord* cuckoo_path,
                                 const size_t i1, const size_t i2) {
        b_slot x = slot_search(ti, i1, i2);
        if (x.depth == -1) {
            return -1;
        }
        // Fill in the cuckoo path slots from the end to the beginning
        for (int i = x.depth; i >= 0; i--) {
            cuckoo_path[i].slot = x.pathcode % SLOT_PER_BUCKET;
            x.pathcode /= SLOT_PER_BUCKET;
        }
        /* Fill in the cuckoo_path buckets and keys from the beginning
         * to the end, using the final pathcode to figure out which
         * bucket the path starts on. Since data could have been
         * modified between slot_search and the computation of the
         * cuckoo path, this could be an invalid cuckoo_path. */
        CuckooRecord *curr = cuckoo_path;
        if (x.pathcode == 0) {
            curr->bucket = i1;
            lock_path(ti, curr->bucket);
            if (ti->buckets_[curr->bucket].hasmigrated) {
                    unlock_path(ti, curr->bucket);
                    return -1;
            }
            if (!ti->buckets_[curr->bucket].occupied[curr->slot]) {
                // We can terminate here
                unlock_path(ti, curr->bucket);
                return 0;
            }
            curr->key = ti->buckets_[curr->bucket].keys[curr->slot];
            unlock_path(ti, curr->bucket);
        } else {
            assert(x.pathcode == 1);
            curr->bucket = i2;
            lock_path(ti, curr->bucket);
            if (ti->buckets_[curr->bucket].hasmigrated) {
                    unlock_path(ti, curr->bucket);
                    return -1;
            }
            if (!ti->buckets_[curr->bucket].occupied[curr->slot]) {
                // We can terminate here
                unlock_path(ti, curr->bucket);
                return 0;
            }
            curr->key = ti->buckets_[curr->bucket].keys[curr->slot];
            unlock_path(ti, curr->bucket);
        }
        for (int i = 1; i <= x.depth; i++) {
            CuckooRecord *prev = curr++;
            const size_t prevhv = hashed_key(prev->key);
            assert(prev->bucket == index_hash(ti, prevhv) ||
                   prev->bucket == alt_index(ti, prevhv, index_hash(ti, prevhv)));
            // We get the bucket that this slot is on by computing the
            // alternate index of the previous bucket
            curr->bucket = alt_index(ti, prevhv, prev->bucket);
            lock_path(ti, curr->bucket);
            if (ti->buckets_[curr->bucket].hasmigrated) {
                    unlock_path(ti, curr->bucket);
                    return -1;
            }
            if (!ti->buckets_[curr->bucket].occupied[curr->slot]) {
                // We can terminate here
                unlock_path(ti, curr->bucket);
                return i;
            }
            curr->key = ti->buckets_[curr->bucket].keys[curr->slot];
            unlock_path(ti, curr->bucket);
        }
        return x.depth;
    }


    /* cuckoopath_move moves keys along the given cuckoo path in order
     * to make an empty slot in one of the buckets in cuckoo_insert.
     * Before the start of this function, the two insert-locked
     * buckets were unlocked in run_cuckoo. At the end of the
     * function, if the function returns true (success), then the last
     * bucket it looks at (which is either i1 or i2 in run_cuckoo)
     * remains locked. If the function is unsuccessful, then both
     * insert-locked buckets will be unlocked. */
    static bool cuckoopath_move(TableInfo *ti, CuckooRecord* cuckoo_path,
                                size_t depth, const size_t i1, const size_t i2) {

        if (depth == 0) {
            /* There is a chance that depth == 0, when
             * try_add_to_bucket sees i1 and i2 as full and
             * cuckoopath_search finds one empty. In this case, we
             * lock both buckets. If the bucket that cuckoopath_search
             * found empty isn't empty anymore, we unlock them and
             * return false. Otherwise, the bucket is empty and
             * insertable, so we hold the locks and return true. */
            const size_t bucket = cuckoo_path[0].bucket;
            assert(bucket == i1 || bucket == i2);
            lock_write_two(ti, i1, i2);
            if (ti->buckets_[i1].hasmigrated || ti->buckets_[i2].hasmigrated) {
                unlock_write_two(ti, i1, i2);
                return false;
            }
            if (!ti->buckets_[bucket].occupied[cuckoo_path[0].slot]) {
                return true;
            } else {
                unlock_write_two(ti, i1, i2);
                return false;
            }
        }

        while (depth > 0) {
            CuckooRecord *from = cuckoo_path + depth - 1;
            CuckooRecord *to   = cuckoo_path + depth;
            size_t fb = from->bucket;
            size_t fs = from->slot;
            size_t tb = to->bucket;
            size_t ts = to->slot;

            size_t ob = 0;
            if (depth == 1) {
                /* Even though we are only swapping out of i1 or i2,
                 * we have to lock both of them along with the slot we
                 * are swapping to, since at the end of this function,
                 * i1 and i2 must be locked. */
                ob = (fb == i1) ? i2 : i1;
                lock_write_three(ti, fb, tb, ob);
            } else {
                lock_write_two(ti, fb, tb);
            }

            if (ti->buckets_[fb].hasmigrated ||
                ti->buckets_[tb].hasmigrated ||
                ti->buckets_[ob].hasmigrated) {

                if (depth == 1) {
                    unlock_write_three(ti, fb, tb, ob);
                } else {
                    unlock_write_two(ti, fb, tb);
                }
                return false;
            }
            /* We plan to kick out fs, but let's check if it is still
             * there; there's a small chance we've gotten scooped by a
             * later cuckoo. If that happened, just... try again. Also
             * the slot we are filling in may have already been filled
             * in by another thread, or the slot we are moving from
             * may be empty, both of which invalidate the swap. */
            if (!eqfn(ti->buckets_[fb].keys[fs], from->key) ||
                ti->buckets_[tb].occupied[ts] ||
                !ti->buckets_[fb].occupied[fs]) {

                if (depth == 1) {
                    unlock_write_three(ti, fb, tb, ob);
                } else {
                    unlock_write_two(ti, fb, tb);
                }
                return false;
            }

            ti->buckets_[tb].setKV(ts, ti->buckets_[fb].keys[fs], ti->buckets_[fb].vals[fs]);
            ti->buckets_[fb].eraseKV(fs);
            ti->buckets_[fb].need_check_alternate = true;
            if (depth == 1) {
                // Don't unlock fb or ob, since they are needed in
                // cuckoo_insert. Only unlock tb if it doesn't unlock
                // the same bucket as fb or ob.
                if (tb != fb && tb != ob) {
                    unlock_write(ti, tb);
                }
            } else {
                unlock_write_two(ti, fb, tb);
            }
            depth--;
        }
        return true;
    }

    /* run_cuckoo performs cuckoo hashing on the table in an attempt
     * to free up a slot on either i1 or i2. On success, the bucket
     * and slot that was freed up is stored in insert_bucket and
     * insert_slot. In order to perform the search and the swaps, it
     * has to unlock both i1 and i2, which can lead to certain
     * concurrency issues, the details of which are explained in the
     * function. If run_cuckoo returns ok (success), then the slot it
     * freed up is still locked. Otherwise it is unlocked. */
    cuckoo_status run_cuckoo(TableInfo *ti, const size_t i1, const size_t i2,
                             size_t &insert_bucket, size_t &insert_slot) {

        CuckooRecord cuckoo_path[MAX_BFS_DEPTH+1];

        // We must unlock i1 and i2 here, so that cuckoopath_search
        // and cuckoopath_move can lock buckets as desired without
        // deadlock. cuckoopath_move has to look at either i1 or i2 as
        // its last slot, and it will lock both buckets and leave them
        // locked after finishing. This way, we know that if
        // cuckoopath_move succeeds, then the buckets needed for
        // insertion are still locked. If cuckoopath_move fails, the
        // buckets are unlocked and we try again. This unlocking does
        // present two problems. The first is that another insert on
        // the same key runs and, finding that the key isn't in the
        // table, inserts the key into the table. Then we insert the
        // key into the table, causing a duplication. To check for
        // this, we search i1 and i2 for the key we are trying to
        // insert before doing so (this is done in cuckoo_insert, and
        // requires that both i1 and i2 are locked). Another problem
        // is that an expansion runs and changes table_info, meaning
        // the cuckoopath_move and cuckoo_insert would have operated
        // on an old version of the table, so the insert would be
        // invalid. For this, we check that ti == table_info.load()
        // after cuckoopath_move, signaling to the outer insert to try
        // again if the comparison fails.
        unlock_write_two(ti, i1, i2);

        bool done = false;
        while (!done) {
            int depth = cuckoopath_search(ti, cuckoo_path, i1, i2);
            if (depth < 0) {
                break;
            }

            if (cuckoopath_move(ti, cuckoo_path, depth, i1, i2)) {
                insert_bucket = cuckoo_path[0].bucket;
                insert_slot = cuckoo_path[0].slot;
                assert(insert_bucket == i1 || insert_bucket == i2);
                assert(ti->buckets_[i1].version.num.load() % 2 == 1);
                assert(ti->buckets_[i2].version.num.load() % 2 == 1);
                assert(!ti->buckets_[insert_bucket].occupied[insert_slot]);
                done = true;
                break;
            }
        }

        if (!done) {
            return failure; //happens if path is too long
        } 
        //don't need below because we don't care if operating on old table,
        //will find new table anyways
        /*else if (ti != table_info.load()) {
            // Unlock i1 and i2 and signal to cuckoo_insert to try
            // again. Since we set the hazard pointer to be ti, this
            // check isn't susceptible to an ABA issue, since a new
            // pointer can't have the same address as ti.
            unlock_write_two(ti, i1, i2);
            return failure_under_expansion;
        }*/
        return ok;
    }

    /* cuckoo_insert tries to insert the given key-value pair into an
     * empty slot in i1 or i2, performing cuckoo hashing if necessary.
     * It expects the locks to be taken outside the function, but they
     * are released here, since different scenarios require different
     * handling of the locks. Before inserting, it checks that the key
     * isn't already in the table. cuckoo hashing presents multiple
     * concurrency issues, which are explained in the function. */
    cuckoo_status cuckoo_insert(const key_type &key, const mapped_type &val,
                                const size_t hv, TableInfo *ti,
                                const size_t i1, const size_t i2) {
        mapped_type oldval;
        int open1, open2;
        cuckoo_status res1, res2;
        
        res1 = try_add_to_bucket(ti, key, val, i1, open1);
        if (res1 == failure_key_duplicated) {
            unlock_write_two(ti, i1, i2);
            //std::cout << "Key duplicated" << key << " , " << val << " , " << hv << std::endl;
            return failure_key_duplicated;
        }

        res2 = try_add_to_bucket(ti, key, val, i2, open2);
        if (res2 == failure_key_duplicated) {
            unlock_write_two(ti, i1, i2);
            //std::cout << "Key duplicated" << key << " , " << val << " , " << hv << std::endl;
            return failure_key_duplicated;
        }

        if (res1 == failure_key_moved || res2 == failure_key_moved) {
            unlock_write_two(ti, i1, i2);
            std::cout << "Key moved" << key << " , " << val << " , " << hv << std::endl;
            return failure_key_moved;
        }

        if (open1 != -1) {
            add_to_bucket(ti, key, val, i1, open1);
            unlock_write_two(ti, i1, i2);
            return ok;
        }
        ti->buckets_[i1].need_check_alternate = true;

        //std::cout << "Have to check second bucket" << std::endl;
        if (open2 != -1) {
            add_to_bucket(ti, key, val, i2, open2);
            unlock_write_two(ti, i1, i2);
            return ok;
        }

        // we are unlucky, so let's perform cuckoo hashing
        //std::cout << "Have to run cuckoo hashing" << std::endl;
        size_t insert_bucket = 0;
        size_t insert_slot = 0;
        cuckoo_status st = run_cuckoo(ti, i1, i2, insert_bucket, insert_slot);
        if (st == ok) {
            //TODO: remove
            assert(ti->buckets_[i1].version.num.load() % 2 == 1);
            assert(ti->buckets_[i2].version.num.load() % 2 == 1);
            assert(!ti->buckets_[insert_bucket].occupied[insert_slot]);
            assert(insert_bucket == index_hash(ti, hv) || insert_bucket == alt_index(ti, hv, index_hash(ti, hv)));
            /* Since we unlocked the buckets during run_cuckoo,
             * another insert could have inserted the same key into
             * either i1 or i2, so we check for that before doing the
             * insert. */
            if (cuckoo_find(key, oldval, hv, ti, i1, i2) == ok) {
                unlock_write_two(ti, i1, i2);
                return failure_key_duplicated;
            }
            add_to_bucket(ti, key, val, insert_bucket, insert_slot);
            unlock_write_two(ti, i1, i2);
            return ok;
        }

        assert(st == failure);
        std::cout << "Hash table is full (hashpower = " << ti->hashpower_ << "with load factor = " << cuckoo_loadfactor(ti) << std::endl;
        LIBCUCKOO_DBG("hash table is full (hashpower = %zu, hash_items = %zu, load factor = %.2f), need to increase hashpower\n",
                      ti->hashpower_, cuckoo_size(ti), cuckoo_loadfactor(ti));

        //we will create a new table and set the new table pointer to it
        if (cuckoo_expand_start(ti->hashpower_+1) == failure_under_expansion) {
            std::cout << "Somebody swapped the table pointer before I did. Anyways, it's changed!" << std::endl;
            return failure_key_moved; //TODO: Change to a different error message?
        }

        return failure_key_moved;
    }

        /* try_read_from-bucket will search the bucket for the given key
     * and store the associated value if it finds it. */
    static cuckoo_status try_read_from_bucket(const TableInfo *ti, 
                                     const key_type &key, mapped_type &val,
                                     const size_t i) {
        if (ti->buckets_[i].hasmigrated) {
            return failure_key_moved;
        }

        for (size_t j = 0; j < SLOT_PER_BUCKET; j++) {
            if (!ti->buckets_[i].occupied[j]) {
                continue;
            }

            if (eqfn(key, ti->buckets_[i].keys[j])) {
                val = ti->buckets_[i].vals[j];
                return ok;
            }
        }
        return failure_key_not_found;
    }

    /* add_to_bucket will insert the given key-value pair into the
     * slot. */
    static inline void add_to_bucket(TableInfo *ti, 
                                     const key_type &key, const mapped_type &val,
                                     const size_t i, const size_t j) {
        assert(!ti->buckets_[i].occupied[j]);
        
        ti->buckets_[i].setKV(j, key, val);
        ti->num_inserts[counterid].num.fetch_add(1, std::memory_order_relaxed);
    }

    /* try_add_to_bucket will search the bucket and store the index of
     * an empty slot if it finds one, or -1 if it doesn't. Regardless,
     * it will search the entire bucket and return false if it finds
     * the key already in the table (duplicate key error) and true
     * otherwise. */
    static cuckoo_status try_add_to_bucket(TableInfo *ti, 
                                  const key_type &key, const mapped_type &val,
                                  const size_t i, int& j) {
        j = -1;
        bool found_empty = false;

        if (ti->buckets_[i].hasmigrated) {
            return failure_key_moved;
        }
        for (size_t k = 0; k < SLOT_PER_BUCKET; k++) {
            if (ti->buckets_[i].occupied[k]) {
                if (eqfn(key, ti->buckets_[i].keys[k])) {
                    return failure_key_duplicated;
                }
            } else {
                if (!found_empty) {
                    found_empty = true;
                    j = k;
                }
            }
        }
        return ok;
    }

    /* try_del_from_bucket will search the bucket for the given key,
     * and set the slot of the key to empty if it finds it. */
    static cuckoo_status try_del_from_bucket(TableInfo *ti, 
                                    const key_type &key, const size_t i) {
        if (ti->buckets_[i].hasmigrated) {
            return failure_key_moved;
        }

        for (size_t j = 0; j < SLOT_PER_BUCKET; j++) {
            if (!ti->buckets_[i].occupied[j]) {
                continue;
            }
            if (eqfn(ti->buckets_[i].keys[j], key)) {
                ti->buckets_[i].eraseKV(j);
                ti->num_deletes[counterid].num.fetch_add(1, std::memory_order_relaxed);
                return ok;
            }
        }
        return failure_key_not_found;
    }

    /* cuckoo_find searches the table for the given key and value,
     * storing the value in the val if it finds the key. It expects
     * the locks to be taken and released outside the function. */
    static cuckoo_status cuckoo_find(const key_type& key, mapped_type& val,
                                     const size_t hv, const TableInfo *ti,
                                     const size_t i1, const size_t i2) {
        cuckoo_status res1, res2;
        res1 = try_read_from_bucket(ti, key, val, i1);
        if (res1 == ok) {
            return ok;
        }
        res2 = try_read_from_bucket(ti, key, val, i2);
        if(res2 == ok) {
            return ok;
        }

        if (res1 == failure_key_moved || res2 == failure_key_moved) {
            return failure_key_moved;
        }

        return failure_key_not_found;
    }

    /* TODO: Change to support prefetching
     * cuckoo_delete searches the table for the given key and sets the
     * slot with that key to empty if it finds it. It expects the
     * locks to be taken and released outside the function. */
    cuckoo_status cuckoo_delete(const key_type &key, const size_t hv,
                                TableInfo *ti, const size_t i1,
                                const size_t i2) {
        cuckoo_status res1, res2;
        res1 = try_del_from_bucket(ti, key, i1);
        if (res1 == ok) {
            return ok;
        }
        res2 = try_del_from_bucket(ti, key, i2);
        if (res2 == ok) {
            return ok;
        }

        //couldn't find key in either bucket, and a bucket was moved
        if (res1 == failure_key_moved || res2 == failure_key_moved) {
            return failure_key_moved;
        }
        return failure_key_not_found;
    }

    /* cuckoo_init initializes the hashtable, given an initial
     * hashpower as the argument. */
    cuckoo_status cuckoo_init(const size_t hashtable_init) {
        table_info.store(new TableInfo(hashtable_init));
        cuckoo_clear(table_info.load());
        return ok;
    }

    /* count_migrated_buckets returns the number of migrated buckets in the given
     * table. */
    size_t count_migrated_buckets(const TableInfo *ti) {
        size_t num_migrated = 0;
        for (size_t i = 0; i < ti->num_inserts.size(); i++) {
            num_migrated += ti->num_migrated_buckets[i].num.load();
        }
        return num_migrated;
    }

    /* cuckoo_clear empties the table, calling the destructors of all
     * the elements it removes from the table. It assumes the locks
     * are taken as necessary. */
    cuckoo_status cuckoo_clear(TableInfo *ti) {
        for (size_t i = 0; i < hashsize(ti->hashpower_); i++) {
            ti->buckets_[i].clear();
        }
        for (size_t i = 0; i < ti->num_inserts.size(); i++) {
            ti->num_inserts[i].num.store(0);
            ti->num_deletes[i].num.store(0);
            ti->num_migrated_buckets[i].num.store(0);
            ti->num_retries[i].num.store(0);
        }
        return ok;
    }

    /* cuckoo_size returns the number of elements in the given
     * table. */
    size_t cuckoo_size(const TableInfo *ti) {
        if (ti == nullptr) {
            std::cout<< "New table doesn't exist yet" << std::endl;
            return 0;
        }
        size_t inserts = 0;
        size_t deletes = 0;
        for (size_t i = 0; i < ti->num_inserts.size(); i++) {
            inserts += ti->num_inserts[i].num.load();
            deletes += ti->num_deletes[i].num.load();
        }
        return inserts-deletes;
    }

    /* cuckoo_loadfactor returns the load factor of the given table. */
    float cuckoo_loadfactor(const TableInfo *ti) {
        return 1.0 * cuckoo_size(ti) / SLOT_PER_BUCKET / hashsize(ti->hashpower_);
    }

    /* insert_into_table is a helper function used by
     * cuckoo_expand_simple to fill up the new table. */
    static void insert_into_table(cuckoohash_map<Key, T, Hash>& new_map, const TableInfo *old_ti, size_t i, size_t end) {
        for (;i < end; i++) {
            for (size_t j = 0; j < SLOT_PER_BUCKET; j++) {
                if (old_ti->buckets_[i].occupied[j]) {
                    new_map.insert(old_ti->buckets_[i].keys[j], old_ti->buckets_[i].vals[j]);
                }
            }
        }
    }

    // expansion_lock serializes functions that call
    // snapshot_and_lock_all, thereby ensuring that multiple
    // expansions and iterator constructions cannot occur
    // simultaneously.
    std::mutex expansion_lock;
    std::mutex migrate_all_lock;

    /* cuckoo_expand_start tries to create a new table, succeeding as long as 
     * there is no ongoing expansion and no other thread has already created a new table
     * (new_table_pointer != nullptr) */
    cuckoo_status cuckoo_expand_start(size_t n) {
        //we only want to create a new table if there is no ongoing expansion already
        //Also accouunts for possibility of aomebody already swapping the table pointer
        if (new_table_info.load() != nullptr) {
            return failure_under_expansion;
        }
        TableInfo *ti = new TableInfo(n); //new, larger table
        TableInfo *expected = nullptr; 
        if (!new_table_info.compare_exchange_weak(expected, ti, 
                                                    std::memory_order_release,
                                                    std::memory_order_relaxed)) {
            delete ti;
            return failure_under_expansion;
        }

        return ok;
    }

    /* cuckoo_expand_end is triggered once all buckets from the old table have been moved over.
     * It tries to swap the new_table pointer to the old, succeeding as long as no other thread
     * has already done so. It then adds the old_table pointer to a list of pointers that will
     * be garbage collected, and sets the new table pointer to nullptr (so for a small period
     * of time we can have both old_table and new_table being the same) */
    cuckoo_status cuckoo_expand_end() {
        TableInfo* old_ti = table_info.load();
        TableInfo* new_ti = new_table_info.load();
        if (!table_info.compare_exchange_weak(old_ti, new_ti, 
                                                    std::memory_order_release,
                                                    std::memory_order_relaxed)) {
            std::cout << "Someone swapped the old table pointer before we could" << std::endl;
            return failure_under_expansion;
        }

        std::cout << "Finished swapping old_ti with new" << std::endl;
        // Rather than deleting ti now, we store it in
        // old_table_infos. The hazard pointer manager will delete it
        // if no other threads are using the pointer.
        old_table_infos.push_back(old_ti);
        global_hazard_pointers.delete_unused(old_table_infos);
        new_table_info.store(nullptr); //TODO: Need to CAS?
        std::cout << "Old pointer is " << old_ti << "New pointer is" << table_info.load() << std::endl;
        std::cout << "New table pointer is " << new_table_info.load() << std::endl;
        return ok;
    }

    /* cuckoo_expand_simple is a simpler version of expansion than
     * cuckoo_expand, which will double the size of the existing hash
     * table. It needs to take all the bucket locks, since no other
     * operations can change the table during expansion. If some other
     * thread is holding the expansion thread at the time, then it
     * will return failure_under_expansion. */
    cuckoo_status cuckoo_expand_simple(size_t n) {
        if (!expansion_lock.try_lock()) {
            unset_hazard_pointer();
            return failure_under_expansion;
        }

        TableInfo *ti = snapshot_and_lock_write_all();
        assert(ti != nullptr);
        if (n <= ti->hashpower_) {
            // Most likely another expansion ran before this one could
            // grab the locks
            unlock_write_all(ti);
            unset_hazard_pointer();
            expansion_lock.unlock();
            return failure_under_expansion;
        }

        try {
            // Creates a new hash table with hashpower n and adds all
            // the elements from the old buckets
            cuckoohash_map<Key, T, Hash> new_map(hashsize(n) * SLOT_PER_BUCKET);
            const size_t threadnum = kNumCores;
            const size_t buckets_per_thread = hashsize(ti->hashpower_) / threadnum;
            std::vector<std::thread> insertion_threads(threadnum);
            for (size_t i = 0; i < threadnum-1; i++) {
                insertion_threads[i] = std::thread(
                    insert_into_table, std::ref(new_map),
                    ti, i*buckets_per_thread, (i+1)*buckets_per_thread);
            }
            insertion_threads[threadnum-1] = std::thread(
                insert_into_table, std::ref(new_map), ti,
                (threadnum-1)*buckets_per_thread, hashsize(ti->hashpower_));
            for (size_t i = 0; i < threadnum; i++) {
                insertion_threads[i].join();
            }
            // Sets this table_info to new_map's. It then sets new_map's
            // table_info to nullptr, so that it doesn't get deleted when
            // new_map goes out of scope
            table_info.store(new_map.table_info.load());
            new_map.table_info.store(nullptr);
        } catch (const std::bad_alloc&) {
            // Unlocks resources and rethrows the exception
            unlock_write_all(ti);
            unset_hazard_pointer();
            expansion_lock.unlock();
            throw;
        }

        // Rather than deleting ti now, we store it in
        // old_table_infos. The hazard pointer manager will delete it
        // if no other threads are using the pointer.
        old_table_infos.push_back(ti);
        unlock_write_all(ti);
        unset_hazard_pointer();
        global_hazard_pointers.delete_unused(old_table_infos);
        expansion_lock.unlock();
        return ok;
    }

    // Iterator definitions
    friend class const_iterator;
    friend class iterator;
public:

    /*! A const_iterator is an iterator through the table that is
     * thread safe. For the duration of its existence, it takes all
     * the locks on the table it is given, thereby ensuring that no
     * other threads can modify the table while the iterator is in
     * use. Note that this also means that only one iterator can be
     * active on a table at one time and furthermore that all
     * operations on the table, except the \ref size, \ref empty, \ref
     * hashpower, \ref bucket_count, and \ref load_factor methods,
     * will stall until the iterator loses its lock. For this reason,
     * we suggest using the \ref snapshot_table method if possible,
     * since it is less error-prone. The iterator allows movement
     * forward and backward through the table as well as dereferencing
     * items in the table. It maintains the invariant that the
     * iterator is either an end iterator (which points past the end
     * of the table), or points to a filled slot. As soon as the
     * iterator looses its lock on the table, all dereference and
     * movement operations will throw an exception. */
    class const_iterator {
        /* The constructor locks the entire table, retrying until
         * snapshot_and_lock_all succeeds. Then it calculates end_pos
         * and begin_pos and sets index and slot to the beginning or
         * end of the table, based on the boolean argument. We keep
         * this constructor private (but expose it to the
         * cuckoohash_map class), since we don't want users calling
         * it. */
        const_iterator(cuckoohash_map<Key, T, Hash, Pred> *hm, bool is_end) {
            cuckoohash_map<Key, T, Hash, Pred>::check_hazard_pointer();
            hm_ = hm;
            hm->expansion_lock.lock();
            ti_ = hm_->snapshot_and_lock_write_all();
            assert(ti_ != nullptr);

            has_table_lock = true;

            index_ = slot_ = 0;

            set_end(end_pos.first, end_pos.second);
            set_begin(begin_pos.first, begin_pos.second);
            if (is_end) {
                index_ = end_pos.first;
                slot_ = end_pos.second;
            } else {
                index_ = begin_pos.first;
                slot_ = begin_pos.second;
            }
        }

        friend class cuckoohash_map<Key, T, Hash, Pred>;

    public:

        /*! This is an rvalue-reference constructor that takes the
          lock from \p it and copies its state. To create an iterator
          from scratch, call the \ref cbegin or \ref cend methods of
          cuckoohash_map. */
        const_iterator(const_iterator&& it) {
            if (this == &it) {
                return;
            }
            memcpy(this, &it, sizeof(const_iterator));
            it.has_table_lock = false;
        }

        /*! The assignment operator behaves identically to the
         * rvalue-reference constructor. */
        const_iterator* operator=(const_iterator&& it) {
            if (this == &it) {
                return this;
            }
            memcpy(this, &it, sizeof(const_iterator));
            it.has_table_lock = false;
            return this;
        }

        /*! release unlocks the table, thereby freeing it up for other
         * operations, but also invalidating all future operations
         * with this iterator. */
        void release() {
            if (has_table_lock) {
                hm_->unlock_write_all(ti_);
                hm_->expansion_lock.unlock();
                cuckoohash_map<Key, T, Hash, Pred>::unset_hazard_pointer();
                has_table_lock = false;
            }
        }

        /*! The destructor simply calls \ref release. */
        ~const_iterator() {
            release();
        }

        /*! is_end returns true if the iterator is at end_pos, which means
         * it is past the end of the table. */
        bool is_end() {
            return (index_ == end_pos.first && slot_ == end_pos.second);
        }

        /*! is_begin returns true if the iterator is at begin_pos, which
         * means it is at the first item in the table. */
        bool is_begin() {
            return (index_ == begin_pos.first && slot_ == begin_pos.second);
        }

    protected:
        // For the arrow dereference operator, we return a pointer to
        // a lightweight pair consisting of const references to the
        // key and value under the iterator.
        typedef std::pair<const Key&, const T&> ref_pair;
        // Since we can't initialize a ref_pair before knowing what it
        // points to, we use std::aligned_storage to reserve
        // unititialized space for the object, which we then construct
        // with placement new.
        typename std::aligned_storage<sizeof(ref_pair), std::alignment_of<ref_pair>::value>::type data;

    public:
        /*! The dereference operator returns a value_type copied from
         *  the key-value pair under the iterator. */
        value_type operator*() {
            check_lock();
            if (is_end()) {
                throw std::out_of_range(end_dereference);
            }
            assert(ti_->buckets_[index_].occupied[slot_]);
            return {ti_->buckets_[index_].keys[slot_], ti_->buckets_[index_].vals[slot_]};
        }

        /*! The arrow dereference operator returns a pointer to an
         *  internal std::pair which contains const references to the
         *  key and value under the iterator. */
        ref_pair* operator->() {
            check_lock();
            if (is_end()) {
                throw std::out_of_range(end_dereference);
            }
            assert(ti_->buckets_[index_].occupied[slot_]);
            ref_pair *data_ptr = static_cast<ref_pair*>(static_cast<void*>(&data));
            new (data_ptr) ref_pair(ti_->buckets_[index_].keys[slot_], ti_->buckets_[index_].vals[slot_]);
            return data_ptr;
        }

        /*! The prefix increment operator moves the iterator forwards
         * to the next nonempty slot. If it reaches the end of the
         * table, it becomes an end iterator. It throws an exception
         * if the iterator is already at the end of the table. */
        const_iterator* operator++() {
            check_lock();
            if (is_end()) {
                throw std::out_of_range(end_increment);
            }
            forward_filled_slot(index_, slot_);
            return this;
        }

        /*! The postfix increment operator behaves identically to the
         *  prefix increment operator. */
        const_iterator* operator++(int) {
            check_lock();
            if (is_end()) {
                throw std::out_of_range(end_increment);
            }
            forward_filled_slot(index_, slot_);
            return this;
        }

        /*! The prefix decrement operator moves the iterator backwards
         * to the previous nonempty slot. If we aren't at the
         * beginning, then the backward_filled_slot operation should
         * not fail. If we are, it throws an exception. */
        const_iterator* operator--() {
            check_lock();
            if (is_begin()) {
                throw std::out_of_range(begin_decrement);
            }
            backward_filled_slot(index_, slot_);
            return this;
        }

        /*! The postfix decrement operator behaves identically to the
         *  prefix decrement operator. */
        const_iterator* operator--(int) {
            check_lock();
            if (is_begin()) {
                throw std::out_of_range(begin_decrement);
            }
            backward_filled_slot(index_, slot_);
            return this;
        }

    protected:
        // A pointer to the associated hashmap
        cuckoohash_map<Key, T, Hash, Pred> *hm_;

        // The hashmap's table info
        typename cuckoohash_map<Key, T, Hash, Pred>::TableInfo *ti_;

        // Indicates whether the iterator has the table lock
        bool has_table_lock;

        // Stores the bucket and slot of the end iterator, which is one
        // past the end of the table. It is initialized during the
        // iterator's constructor.
        std::pair<size_t, size_t> end_pos;

        // Stotres the bucket and slot of the begin iterator, which is the
        // first filled position in the table. It is initialized during
        // the iterator's constructor. If the table is empty, it points
        // past the end of the table, to the same position as end_pos.
        std::pair<size_t, size_t> begin_pos;

        // The bucket index of the item being pointed to
        size_t index_;

        // The slot in the bucket of the item being pointed to
        size_t slot_;

        /* set_end sets the given index and slot to one past the last
         * position in the table. */
        void set_end(size_t& index, size_t& slot) {
            index = hm_->bucket_count();
            slot = 0;
        }

        /* set_begin sets the given pair to the position of the first
         * element in the table. */
        void set_begin(size_t& index, size_t& slot) {
            if (hm_->empty()) {
                set_end(index, slot);
            } else {
                index = slot = 0;
                // There must be a filled slot somewhere in the table
                if (!ti_->buckets_[index].occupied[slot]) {
                    forward_filled_slot(index, slot);
                    assert(!is_end());
                }
            }
        }

        /* forward_slot moves the given index and slot to the next
         * available slot in the forwards direction. It returns true if it
         * successfully advances, and false if it has reached the end of
         * the table, in which case it sets index and slot to end_pos. */
        bool forward_slot(size_t& index, size_t& slot) {
            if (slot < SLOT_PER_BUCKET-1) {
                ++slot;
                return true;
            } else if (index < hm_->bucket_count()-1) {
                ++index;
                slot = 0;
                return true;
            } else {
                set_end(index, slot);
                return false;
            }
        }

        /* backward_slot moves index and slot to the next available slot
         * in the backwards direction. It returns true if it successfully
         * advances, and false if it has reached the beginning of the
         * table, setting the index and slot back to begin_pos. */
        bool backward_slot(size_t& index, size_t& slot) {
            if (slot > 0) {
                --slot;
                return true;
            } else if (index > 0) {
                --index;
                slot = SLOT_PER_BUCKET-1;
                return true;
            } else {
                set_begin(index, slot);
                return false;
            }
        }

        /* forward_filled_slot moves index and slot to the next filled
         * slot. */
        bool forward_filled_slot(size_t& index, size_t& slot) {
            bool res = forward_slot(index, slot);
            if (!res) {
                return false;
            }
            while (!ti_->buckets_[index].occupied[slot]) {
                res = forward_slot(index, slot);
                if (!res) {
                    return false;
                }
            }
            return true;
        }

        /* backward_filled_slot moves index and slot to the previous
         * filled slot. */
        bool backward_filled_slot(size_t& index, size_t& slot) {
            bool res = backward_slot(index, slot);
            if (!res) {
                return false;
            }
            while (!ti_->buckets_[index].occupied[slot]) {
                res = backward_slot(index, slot);
                if (!res) {
                    return false;
                }
            }
            return true;
        }


        /* check_lock throws an exception if the iterator doesn't have a
         * lock. */
        void check_lock() {
            if (!has_table_lock) {
                throw std::runtime_error("Iterator does not have a lock on the table");
            }
        }

        // Other error messages
        static constexpr char end_dereference[] =
            "Cannot dereference: iterator points past the end of the table";
        static constexpr char end_increment[] =
            "Cannot increment: iterator points past the end of the table";
        static constexpr char begin_decrement[] =
            "Cannot decrement: iterator points to the beginning of the table";
    };


    /*! An iterator supports the same operations as the const_iterator
     *  and provides an additional \ref set_value method to allow
     *  changing values in the table. */
    class iterator : public const_iterator {
        /* This constructor does the same thing as the private
         * const_iterator one. */
        iterator(cuckoohash_map<Key, T, Hash, Pred> *hm, bool is_end)
            : const_iterator(hm, is_end) {}

        friend class cuckoohash_map<Key, T, Hash, Pred>;

    public:

        /*! This constructor is identical to the rvalue-reference
         *  constructor of const_iterator. */
        iterator(iterator&& it)
            : const_iterator(std::move(it)) {}

        /*! This constructor allows converting from a const_iterator
         *  to an iterator. */
        iterator(const_iterator&& it)
            : const_iterator(std::move(it)) {}

        /* The assignment operator behaves identically to the
         * rvalue-reference constructor. */
        iterator* operator=(iterator&& it) {
            if (this == &it) {
                return this;
            }
            memcpy(this, &it, sizeof(iterator));
            it.has_table_lock = false;
            return this;
        }

        /*! set_value sets the value pointed to by the iterator to \p
         * val. This involves modifying the hash table itself, but
         * since we have a lock on the table, we are okay. We are only
         * changing the value in the bucket, so the element will
         * retain it's position in the table. */
        void set_value(const mapped_type val) {
            this->check_lock();
            if (this->is_end()) {
                throw std::out_of_range(this->end_dereference);
            }
            assert(this->ti_->buckets_[this->index_].occupied[this->slot_]);
            this->ti_->buckets_[this->index_].vals[this->slot_] = val;
        }
    };

// Public iterator functions
public:
    /*! cbegin returns a const_iterator to the first filled slot in the
     * table. */
    const_iterator cbegin() {
        return const_iterator(this, false);
    }

    /*! cend returns a const_iterator set past the end of the table. */
    const_iterator cend() {
        return const_iterator(this, true);
    }

    /*! begin returns an iterator to the first filled slot in the
     * table. */
    iterator begin() {
        return iterator(this, false);
    }

    /*! end returns an iterator set past the end of the table. */
    iterator end() {
        return iterator(this, true);
    }

    /*! snapshot_table allocates a vector and, using a const_iterator
     * stores all the elements currently in the table. */
    std::vector<value_type> snapshot_table() {
        const_iterator it = cbegin();
        size_t table_size = size();
        std::vector<value_type> items(table_size);
        size_t ind = 0;
        while (!it.is_end()) {
            items[ind++] = *it;
            it++;
        }
        return items;
    }
};

// Initializing the static members
template <class Key, class T, class Hash, class Pred>
__thread void** cuckoohash_map<Key, T, Hash, Pred>::hazard_pointer = nullptr;

template <class Key, class T, class Hash, class Pred>
__thread int cuckoohash_map<Key, T, Hash, Pred>::counterid = -1;

template <class Key, class T, class Hash, class Pred>
typename cuckoohash_map<Key, T, Hash, Pred>::hasher
cuckoohash_map<Key, T, Hash, Pred>::hashfn;

template <class Key, class T, class Hash, class Pred>
typename cuckoohash_map<Key, T, Hash, Pred>::key_equal
cuckoohash_map<Key, T, Hash, Pred>::eqfn;

template <class Key, class T, class Hash, class Pred>
typename std::allocator<typename cuckoohash_map<Key, T, Hash, Pred>::Bucket>
cuckoohash_map<Key, T, Hash, Pred>::bucket_allocator;


template <class Key, class T, class Hash, class Pred>
typename cuckoohash_map<Key, T, Hash, Pred>::GlobalHazardPointerList
cuckoohash_map<Key, T, Hash, Pred>::global_hazard_pointers;

template <class Key, class T, class Hash, class Pred>
const size_t cuckoohash_map<Key, T, Hash, Pred>::kNumCores =
    std::thread::hardware_concurrency() == 0 ?
    sysconf(_SC_NPROCESSORS_ONLN) : std::thread::hardware_concurrency();

template <class Key, class T, class Hash, class Pred>
constexpr char cuckoohash_map<Key, T, Hash, Pred>::const_iterator::end_dereference[62];

template <class Key, class T, class Hash, class Pred>
constexpr char cuckoohash_map<Key, T, Hash, Pred>::const_iterator::end_increment[60];

template <class Key, class T, class Hash, class Pred>
constexpr char cuckoohash_map<Key, T, Hash, Pred>::const_iterator::begin_decrement[64];

#endif
