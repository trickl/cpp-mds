CC = g++
LIBTOOL = libtool
CDEFS   = 
CFLAGS	= -Wall -g -I/home/tgee/projects/c++/algo/src -I/usr/include/boost-1_40_0
LFLAGS	= -L/home/tgee/projects/c++/algo/lib -lmds -lgsl -liosmanip -llapack
SRCS	= main.cpp 

TARGET	= smacof_1
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

