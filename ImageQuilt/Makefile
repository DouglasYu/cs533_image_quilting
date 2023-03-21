CC = g++


TARGET = ImageQuiltBin
SRCS = $(wildcard ./ImageQuilt/*.cpp)
OBJS = $(SRCS:.cpp=.o)

CXXFLAGS = -std=c++11 -Wall -Iextlibs/include/ -IImageQuilt/

all: clean $(TARGET)

$(TARGET): $(OBJS) $
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(OBJS) $(TARGET)
