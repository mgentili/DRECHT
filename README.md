This is an implementation of a dynamically resizable, concurrent hashtable based heavily off of the [libcuckoo library](http://efficient.github.io/libcuckoo/). 

It supports concurrent inserts and deletes, and is dynamically resizable in the sense that upon failing to find a sufficiently short cuckoo path, a new table will be allocated, and the table will resume operation, all the while migrating buckets from the new table to the old one.

To generate performance graphs, navigate to the tests directory
$ cd tests 
$ python benchmarks.py -g -i -z -r

The -g flag tells it to generate graphs, the -i flag tells to it to generate insert data, the -z flag is for resizing data, and the -r flag if for read data. 
