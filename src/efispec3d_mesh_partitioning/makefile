
## Main config file
include ../../config/make.config

## Objects and sources
COBJS	= $(CSRCS:.c=.o)
CSRCS	= $(wildcard *.c)

all:	$(COBJS)
	$(CC) $(CFLAGS) $(COBJS) $(LIBS) -o $(EXEC_MESH)  

clean:
	$(RM) $(COBJS) $(EXEC_MESH)



