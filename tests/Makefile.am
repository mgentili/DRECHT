LDFLAGS = -lpthread -L$(top_builddir)/libcuckoo -lcityhash -ljsoncpp
DEPS = $(top_builddir)/libcuckoo/libcityhash.la

noinst_PROGRAMS = test_resizing.out test_insert_and_find.out stress_checked.out insert_throughput.out read_throughput.out mixed_throughput.out
TESTS = test_resizing.out test_insert_and_find.out stress_checked.out insert_throughput.out read_throughput.out mixed_throughput.out

test_insert_and_find_out_SOURCES = test_insert_and_find.cc
stress_checked_out_SOURCES = stress_checked.cc
insert_throughput_out_SOURCES = insert_throughput.cc
read_throughput_out_SOURCES = read_throughput.cc
test_resizing_out_SOURCE = test_resizing.cc
mixed_throughput_out_SOURCE = mixed_throughput.cc
AM_CPPFLAGS = -I$(top_srcdir)
