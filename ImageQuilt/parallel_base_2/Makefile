CC = g++


TARGET = ImageQuiltBin
SRCS = $(wildcard ./*.cpp)
OBJS = $(SRCS:.cpp=.o)

CXXFLAGS = -std=c++20 -O0 -g -Wall -isystem ../extlibs/include/ -I/

all: clean $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(OBJS) $(TARGET)
