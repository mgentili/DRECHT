LDFLAGS = -lpthread -L$(top_builddir)/libcuckoo -lcityhash
DEPS = $(top_builddir)/libcuckoo/libcityhash.la

noinst_PROGRAMS = test_insert_and_find.out stress_checked.out insert_throughput.out read_throughput.out insert_and_find_retries.out
TESTS = test_insert_and_find.out stress_checked.out insert_throughput.out read_throughput.out insert_and_find_retries.out

test_insert_and_find_out_SOURCES = test_insert_and_find.cc
stress_checked_out_SOURCES = stress_checked.cc
insert_throughput_out_SOURCES = insert_throughput.cc
read_throughput_out_SOURCES = read_throughput.cc
insert_and_find_retries_out_SOURCES = insert_and_find_retries.cc

AM_CPPFLAGS = -I$(top_srcdir)
