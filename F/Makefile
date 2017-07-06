FC=ftn
FFLAGS=-g

.SUFFIXES: .f90

DEPS=Makefile 

OBJ=ioutils.o  
SSEND=mpiio_ssend
BCAST=mpiio_bcast
MASTER_VEC=mpiio_master_vector
MASTER_SUB=mpiio_master_subarray
READ_VEC=mpiio_read_vector
READ_SUB=mpiio_read_subarray
READALL_VEC=mpiio_readall_vector
READALL_SUB=mpiio_readall_subarray
SUB_READ_VEC=mpiio_sub_read_vector
SUB_READ_SUB=mpiio_sub_read_subarray
SUB_READALL_VEC=mpiio_sub_readall_vector
SUB_READALL_SUB=mpiio_sub_readall_subarray
INDIVID=mpiio_individual

all: $(SSEND) $(BCAST) $(MASTER_VEC) $(MASTER_SUB) $(READ_VEC) $(READ_SUB) $(READALL_SUB) $(READALL_VEC) $(SUB_READ_VEC) $(SUB_READ_SUB) $(SUB_READALL_VEC) $(SUB_READALL_SUB) $(INDIVID)

.PHONY:all

%.o:%.f90 $(DEPS) 
	$(FC) -c -o $@ $< $(FFLAGS)

$(BCAST): $(OBJ) $(BCAST).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(SSEND): $(OBJ) $(SSEND).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(INDIVID): $(OBJ) $(INDIVID).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(MASTER_VEC): $(OBJ) $(MASTER_VEC).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(MASTER_SUB): $(OBJ) $(MASTER_SUB).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(READ_VEC): $(OBJ) $(READ_VEC).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(READ_SUB): $(OBJ) $(READ_SUB).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(READALL_VEC): $(OBJ) $(READALL_VEC).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(READALL_SUB): $(OBJ) $(READALL_SUB).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(SUB_READ_VEC): $(OBJ) $(SUB_READ_VEC).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(SUB_READALL_VEC): $(OBJ) $(SUB_READALL_VEC).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(SUB_READ_SUB): $(OBJ) $(SUB_READ_SUB).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

$(SUB_READALL_SUB): $(OBJ) $(SUB_READALL_SUB).o
	$(FC) -o $@ $@.o $(OBJ) $(FFLAGS)

clean:
	rm -fr *.o core* $(SSEND) $(BCAST) $(MASTER_VEC) $(MASTER_SUB) $(READ_VEC) $(READ_SUB) $(READALL_SUB) $(READALL_VEC) $(SUB_READ_VEC) $(SUB_READ_SUB) $(SUB_READALL_VEC) $(SUB_READALL_SUB) $(INDIVID)
