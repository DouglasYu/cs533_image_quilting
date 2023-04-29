CC = g++


TARGET = ImageQuiltBin
SRCS = $(wildcard ./*.cpp)
OBJS = $(SRCS:.cpp=.o)

CXXFLAGS = -std=c++20 -O2 -g -Wall -isystem ../extlibs/include/ -I/ -pthread

all: clean $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(OBJS) $(TARGET)
