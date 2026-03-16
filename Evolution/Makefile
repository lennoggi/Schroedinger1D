# Include and parse the option list, but don't do that when running 'make clean'
ifeq (, ${filter clean, ${MAKECMDGOALS}})
	include ${options}

	# Disable optimization if so requested by the user
	ifeq (${OPTIMIZE}, no)
		CXX_OPTIMIZE_FLAGS =
	else ifneq (${OPTIMIZE}, yes)
		$(error Set option OPTIMIZE to either yes or no in the optionlist (current setting: ${OPTIMIZE}))
	endif

	# Disable warnings if so requested by the user
	ifeq (${WARN}, no)
		CXX_WARN_FLAGS =
	else ifneq (${WARN}, yes)
		$(error Set option WARN to either yes or no in the optionlist (current setting: ${WARN}))
	endif

	# Disable debug mode if so requested by the user
	ifeq (${DEBUG}, no)
		CXX_DEBUG_FLAGS =
	else ifneq (${DEBUG}, yes)
		$(error Set option DEBUG to either yes or no in the optionlist (current setting: ${DEBUG}))
	endif
endif


# Objects to be built from regular C++ (not CUDA) source files
OBJ = Main.o Potentials.o Wave_functions.o
EXE = Schroedinger1D_exe


# Build all targets
build: $(OBJ)
	$(CXX) -o $(EXE) $(OBJ) ${LDFLAGS} -L$(HDF5_LIBS_DIR) $(HDF5_LIBS)

Main.o: Main.cc include/Check_parameters.hh include/Declare_functions.hh Parameters.hh
	$(CXX) $(CXXFLAGS) $(CXX_OPTIMIZE_FLAGS) $(CXX_WARN_FLAGS) $(CXX_DEBUG_FLAGS) -I$(HDF5_INC_DIR) -c Main.cc

Potentials.o: Potentials.cc include/Declare_functions.hh Parameters.hh
	$(CXX) $(CXXFLAGS) $(CXX_OPTIMIZE_FLAGS) $(CXX_WARN_FLAGS) $(CXX_DEBUG_FLAGS) -c Potentials.cc

Wave_functions.o: Wave_functions.cc include/Declare_functions.hh Parameters.hh
	$(CXX) $(CXXFLAGS) $(CXX_OPTIMIZE_FLAGS) $(CXX_WARN_FLAGS) $(CXX_DEBUG_FLAGS) -c Wave_functions.cc


# Remove the executable and all object files
# NOTE: icpx generates *.o.tmp files which should be removed when cleaning
.PHONY : clean
clean:
	${RM} ${EXE} ${OBJ} *.o.tmp
