# === PETSc Configuration === #

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CPPFLAGS += -I${PETSC_DIR}/${PETSC_ARCH}/include
CPPFLAGS += -DPETSC_USE_COMPLEX

CFLAGS   += -O2 -Wall -Wextra

LDLIBS   += -L${PETSC_ARCH}/lib

SRC = \
	main.c \
	reader.c \
	fc_factors.c \
	self_factors.c \
	self_writer.c \
	builder.c \
	assemble.c

OBJ = $(SRC:.c=.o)

TARGET = ham_builder

# === Build rules === #

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CLINKER) -o $@ $^ $(PETSC_LIB) $(LDLIBS)

# Optional test executable
test: test.o
	$(CLINKER) -o $@ $^ $(PETSC_LIB) $(LDLIBS)

# === Cleanup === #

clean:
	$(RM) $(OBJ) test.o $(TARGET) test a.out

# === Debug (uncomment to use) === #

# $(info [DEBUG] CPPFLAGS = $(CPPFLAGS))
# $(info [DEBUG] LDLIBS   = $(LDLIBS))

