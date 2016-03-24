CC = g++
LIBDIR  = ~/projects/c++/algo/lib
LIBTOOL = libtool
CDEFS   = -DNDEBUG
CFLAGS	= -Wall -march=athlon64 -g -std=c++0x  \
          -I/home/tgee/projects/c++/algo/src \
          -I/usr/include/boost-1_40_0  -I/usr/include/boost-sandbox \
          -I/usr/include/mysql -I/usr/include/mysql++ \
          -I/home/tgee/projects/external/libssrckdtree-1.0.0
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
          cpmds.cpp 
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
