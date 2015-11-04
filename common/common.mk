SRC_DIR=$(realpath ../src)
LIB_SRC_DIR=${SRC_DIR}/lib
OBJ_DIR=$(realpath ../obj/${PETSC_ARCH})
BIN_DIR=$(realpath ../bin/${PETSC_ARCH})

LIB_SRCS=fd.c scale.c bps.c fem_inf.c angmoment.c oce1.c oce2.c bspline.c dvr.c mat.c 	pot.c
LIB_OBJS=${addprefix ${OBJ_DIR}/, ${patsubst %.c,%.o,${LIB_SRCS}}}

CFLAGS=-std=gnu99
LIBGSL=-lgsl

${OBJ_DIR}/%.o : ${SRC_DIR}/lib/%.c
	${PCC} -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} -o $@ $< 

my_clean: 
	rm -f ${OBJ_DIR}/*.o
	rm -f ${BIN_DIR}*.out

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}

${BIN_DIR}:
	mkdir -p ${BIN_DIR}

