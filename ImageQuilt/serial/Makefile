CC = g++


TARGET = ImageQuiltBin
SRCS = $(wildcard ./*.cpp)
OBJS = $(SRCS:.cpp=.o)

CXXFLAGS = -std=c++11 -O2 -Wall -isystem ../extlibs/include/ -I/

all: clean $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(OBJS) $(TARGET) $(PARALLEL_OBJS_v1)
