CFLAGS = -Wuninitialized -g -fPIC


ODIR=obj
SRCDIR=src
BASEDIR=base
BASE_LIB_LINK=/home/cades/dev/QuaC/build/base/libbase.a 
EXAMPLESDIR=examples
EXAMPLES=$(basename $(notdir $(wildcard $(EXAMPLESDIR)/*.c)))
TESTDIR=tests
TESTS=$(basename $(notdir $(wildcard $(TESTDIR)/*test*.c)))
MPI_TESTS=$(addprefix mpi_,$(TESTS))
CFLAGS += -isystem $(SRCDIR) -isystem $(BASEDIR)

LIBQUAC=libQuaC.so

include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules

_DEPS = quantum_gates.h dm_utilities.h operators.h solver.h operators_p.h quac.h quac_p.h kron_p.h qasm_parser.h error_correction.h interface_xacc_ir.h
DEPS  = $(patsubst %,$(SRCDIR)/%,$(_DEPS))

_OBJ  = quac.o operators.o solver.o kron.o dm_utilities.o quantum_gates.o error_correction.o qasm_parser.o interface_xacc_ir.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_TEST_OBJ  = unity.o timedep_test.o imag_ham.o
TEST_OBJ = $(patsubst %,$(ODIR)/%,$(_TEST_OBJ))

_TEST_DEPS = tests.h
TEST_DEPS  = $(patsubst %,$(TESTDIR)/%,$(_TEST_DEPS))

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	@mkdir -p $(@D)
	${PETSC_COMPILE} -c -o $@ $< $(CFLAGS) ${PETSC_KSP_LIB} ${PETSC_CC_INCLUDES} ${BASE_LIB_LINK}

$(ODIR)/%.o: $(EXAMPLESDIR)/%.c $(DEPS)
	@mkdir -p $(@D)
	${PETSC_COMPILE} -c -o $@ $< $(CFLAGS) ${PETSC_KSP_LIB} ${PETSC_CC_INCLUDES}

$(ODIR)/%.o: $(TESTDIR)/%.c $(DEPS) $(TEST_DEPS)
	@mkdir -p $(@D)
	@${PETSC_COMPILE} -c -o $@ $< $(CFLAGS) ${PETSC_KSP_LIB} ${PETSC_CC_INCLUDES}

all: examples
shared_lib: ${LIBQUAC}

examples: clean_test $(EXAMPLES)

$(TESTS) : CFLAGS += -DUNIT_TEST
$(TESTS) : % : $(ODIR)/%.o $(OBJ) $(TEST_OBJ)
	${CLINKER} -o $@ $^ $(CFLAGS) ${PETSC_KSP_LIB}
	@echo 'running '$@
	@-./$@ -ts_adapt_type none > tmp_test_results
	-@grep FAIL tmp_test_results || true
	@cat tmp_test_results >> test_results
	@rm tmp_test_results

$(MPI_TESTS) : CFLAGS += -DUNIT_TEST
$(MPI_TESTS) : $(TESTS)
	@$(eval tmp=$(subst mpi_,,$@))
	@echo 'running '$@
	@-mpiexec -np 2 ./$(tmp) -ts_adapt_type none > tmp_test_results
	-@grep FAIL tmp_test_results || true
	@cat tmp_test_results >> test_results
	@rm tmp_test_results

.phony: clean_test test count_fails

count_fails:
	@echo "All failures listed below"
	@grep FAIL test_results || true

clean_test:
	rm -f $(TEST_OBJ)
	@rm -f test_results

test: clean_test $(TESTS) count_fails

mpi_test: clean_test $(MPI_TESTS) count_fails

$(EXAMPLES) : % : $(ODIR)/%.o $(OBJ)
	${CLINKER} -o $@ $^ $(CFLAGS) ${PETSC_KSP_LIB}

${LIBQUAC} : $(OBJ)
	@echo "Linking shared library: "
	${CLINKER} -shared -o ${LIBQUAC} $^ $(CFLAGS) ${PETSC_KSP_LIB} ${BASE_LIB_LINK}

.PHONY: clean

clean:
	rm -f $(ODIR)/*
	rm -f $(EXAMPLES)
	rm -f $(TESTS)
