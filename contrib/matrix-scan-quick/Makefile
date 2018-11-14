CXX      = g++ 
CXXFLAGS = -O3 -Wall
#CXXFLAGS = -O3 -Wall -stdlib=libc++
SRC     = $(wildcard *.cpp) 
OBJS    = $(SRC:.cpp=.o)
APP     = matrix-scan-quick

# compile: $(OBJS)
# 	$(CXX) $(OBJS) -lc++ -o $(APP)

compile: $(OBJS)
	$(CXX) $(OBJS) -o $(APP)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $<

%.o: %.c
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f *.o $(APP)

all: clean compile

