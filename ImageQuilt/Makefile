CC = g++


SERIAL_TARGET = ImageQuiltBin
SERIAL_SRCS = $(wildcard ./ImageQuilt/*.cpp)
SERIAL_OBJS = $(SERIAL_SRCS:.cpp=.o)

PARALLEL_TARGET_v1 = ImageQuiltBin_p1
PARALLEL_SRCS_v1 = $(wildcard ./ImageQuilt_p1/*.cpp)
PARALLEL_OBJS_v1 = $(PARALLEL_SRCS_v1:.cpp=.o)

CXXFLAGS = -std=c++11 -Wall -Iextlibs/include/ -IImageQuilt/

all: clean $(SERIAL_TARGET) $(PARALLEL_TARGET_v1)

$(SERIAL_TARGET): $(SERIAL_OBJS)
	$(CC) $(CXXFLAGS) $^ -o $@

$(PARALLEL_TARGET_v1): $(PARALLEL_OBJS_v1)
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(SERIAL_OBJS) $(SERIAL_TARGET) $(PARALLEL_OBJS_v1) $(PARALLEL_TARGET_v1) 
