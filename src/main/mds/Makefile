CC = g++
LIBDIR  = ~/projects/trickl/cpp-mds/lib
LIBTOOL = libtool
CDEFS   = -DNDEBUG
CFLAGS	= -Wall -march=athlon64 -g -std=c++0x -fpermissive  \
          -I/usr/include/boost  -I/usr/include/boost-sandbox \
          -I/usr/include/ietl
LFLAGS	= 
SRCS	= cmds.cpp \
          smds.cpp \
          fdmds.cpp \
          cha96mds.cpp \
          sqmds.cpp \
          sqpmds.cpp \
          smacof_1.cpp \
          isomap.cpp \
          rkisomap.cpp \
          stress.cpp \
          euclidean_distance.cpp \
          guttman.cpp \
          lmds.cpp
#          wmds.cpp 
#          lmds.cpp \
#          imds.cpp \
#          pdbmds.cpp
#          pmds.cpp \

TARGET	= libmds.la
.SUFFIXES : .lo
OBJS	= $(SRCS:.cpp=.lo)

default : all
all	: $(TARGET)

clean	:
	rm -f $(OBJS) $(TARGET)

clobber	: 
	rm $(TARGET)

depend  :
	makedepend -- $(CFLAGS) -- $(SRCS)

$(TARGET): $(OBJS)
	$(LIBTOOL) --mode=link $(CC) $(LFLAGS) -o $@ $(OBJS) -rpath $(LIBDIR)

.cpp.lo:
	$(LIBTOOL) --mode=compile $(CC) -c $(CDEFS) $(CFLAGS) $< -o $*.o 

install : $(TARGET)
	$(LIBTOOL) --mode=install cp $(TARGET) $(LIBDIR)
