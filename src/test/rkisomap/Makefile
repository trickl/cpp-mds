CC = g++
LIBTOOL = libtool
CDEFS   = 
CFLAGS	= -Wall -g -I../../main -I/usr/include/boost
LFLAGS	= -L/home/tgee/projects/trickl/cpp-mds/lib -lmds -lgsl -llapack
SRCS	= main.cpp 

TARGET	= rkisomap
OBJS	= $(SRCS:.cpp=.o)

default : all
all	: $(TARGET)

clean	:
	rm -f $(OBJS)

clobber	: 
	rm $(TARGET)

$(TARGET): $(OBJS)
	$(LIBTOOL) --mode=link $(CC) $(OBJS) $(LFLAGS) -o $(TARGET)

.cpp.o  :
	$(LIBTOOL) --mode=compile $(CC) -c $(CDEFS) $(CFLAGS) $< -o $*.o

install :
	$(LIBTOOL) install -c $(TARGET) $(BINDIR)

