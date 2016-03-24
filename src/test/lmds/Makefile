CC = g++
LIBTOOL = libtool
CDEFS   = 
CFLAGS	= -Wall -g -I/home/Sleeper/projects/cluster/src -I/usr/include/boost-1_33_1
LFLAGS	= -L/home/Sleeper/projects/cluster/lib -lmds -lgsl -liosmanip -llapack
SRCS	= main.cpp 

TARGET	= lmds
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

