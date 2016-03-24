CC = g++
LIBTOOL = libtool
CDEFS   = 
CFLAGS	= -Wall -g -std=c++0x -I/home/tgee/projects/c++/algo/src -I/usr/include/boost-1_40_0 -I/home/tgee/projects/external/libssrckdtree-1.0.0
LFLAGS	= -L/home/tgee/projects/c++/algo/lib -ldraw -lmds -liosmanip -lcolor -lgsl -llapack -lfltk_gl -lmysqlclient -lmysqlpp
SRCS	= main.cpp 

TARGET	= pdbmds
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

