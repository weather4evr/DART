# Makefile created by mkmf $Id$ 

include ../../../../build_templates/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys jedi_to_dart
fixsys: 
	@ (cd ../../../../assimilation_code/modules/utilities; ./fixsystem $(FC) )
module_map_utils.o: ../../../../models/wrf/module_map_utils.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../models/wrf/module_map_utils.f90
mpi_utilities_mod.o: ../../../../assimilation_code/modules/utilities/mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/mpi_utilities_mod.f90
obs_sequence_mod.o: ../../../../assimilation_code/modules/observations/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/observations/obs_sequence_mod.f90
constants.o: ../../../../observations/obs_converters/jedi2DART/enkf/constants.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/constants.f90
params.o: ../../../../observations/obs_converters/jedi2DART/params.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/params.f90
time_manager_mod.o: ../../../../assimilation_code/modules/utilities/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/time_manager_mod.f90
dart_time_io_mod.o: ../../../../assimilation_code/modules/io/dart_time_io_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o netcdf_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/io/dart_time_io_mod.f90
radinfo.o: ../../../../observations/obs_converters/jedi2DART/radinfo.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/radinfo.f90
dart_obs_seq_mod.o: ../../../../observations/obs_converters/jedi2DART/dart_obs_seq_mod.f90 kinds.o params.o mpisetup.o obs_utilities_mod.o obs_sequence_mod.o obs_def_mod.o types_mod.o obs_kind_mod.o location_mod.o time_manager_mod.o radinfo.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/dart_obs_seq_mod.f90
location_mod.o: ../../../../assimilation_code/location/threed_sphere/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o obs_kind_mod.o mpi_utilities_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/location/threed_sphere/location_mod.f90
sort_mod.o: ../../../../assimilation_code/modules/utilities/sort_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/sort_mod.f90
jedi_to_dart.o: ../../../../observations/obs_converters/jedi2DART/jedi_to_dart.f90 kinds.o mpisetup.o params.o radinfo.o mpi_readobs.o dart_obs_seq_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/jedi_to_dart.f90
mpisetup.o: ../../../../observations/obs_converters/jedi2DART/mpisetup.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/mpisetup.f90
obs_utilities_mod.o: ../../../../observations/obs_converters/utilities/obs_utilities_mod.f90 types_mod.o utilities_mod.o netcdf_utilities_mod.o obs_def_mod.o obs_sequence_mod.o time_manager_mod.o location_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/utilities/obs_utilities_mod.f90
kinds.o: ../../../../observations/obs_converters/jedi2DART/enkf/kinds.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/kinds.f90
random_seq_mod.o: ../../../../assimilation_code/modules/utilities/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/random_seq_mod.f90
types_mod.o: ../../../../assimilation_code/modules/utilities/types_mod.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/types_mod.f90
default_model_mod.o: ../../../../models/utilities/default_model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o netcdf_utilities_mod.o ensemble_manager_mod.o dart_time_io_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../models/utilities/default_model_mod.f90
no_cray_win_mod.o: ../../../../assimilation_code/modules/utilities/no_cray_win_mod.f90 mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/no_cray_win_mod.f90
state_structure_mod.o: ../../../../assimilation_code/modules/io/state_structure_mod.f90 utilities_mod.o obs_kind_mod.o types_mod.o netcdf_utilities_mod.o sort_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/io/state_structure_mod.f90
obs_def_utilities_mod.o: ../../../../observations/forward_operators/obs_def_utilities_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/forward_operators/obs_def_utilities_mod.f90
model_mod.o: ../../../../models/wrf/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o netcdf_utilities_mod.o mpi_utilities_mod.o random_seq_mod.o obs_kind_mod.o ensemble_manager_mod.o sort_mod.o distributed_state_mod.o default_model_mod.o state_structure_mod.o module_map_utils.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../models/wrf/model_mod.f90
readsatobs.o: ../../../../observations/obs_converters/jedi2DART/enkf/readsatobs.f90 kinds.o params.o netcdf_mod.o mpisetup.o jedi_time_routines_mod.o constants.o radinfo.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/readsatobs.f90
jedi_time_routines_mod.o: ../../../../observations/obs_converters/jedi2DART/enkf/jedi_time_routines_mod.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/jedi_time_routines_mod.f90
readconvobs.o: ../../../../observations/obs_converters/jedi2DART/enkf/readconvobs.f90 kinds.o constants.o params.o netcdf_mod.o mpisetup.o jedi_time_routines_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/readconvobs.f90
netcdf_utilities_mod.o: ../../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90
location_io_mod.o: ../../../../assimilation_code/location/utilities/location_io_mod.f90 types_mod.o netcdf_utilities_mod.o location_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/location/utilities/location_io_mod.f90
ensemble_manager_mod.o: ../../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o sort_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90
readozobs.o: ../../../../observations/obs_converters/jedi2DART/enkf/readozobs.f90 kinds.o params.o mpisetup.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/readozobs.f90
assim_model_mod.o: ../../../../assimilation_code/modules/assimilation/assim_model_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o ensemble_manager_mod.o model_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/assimilation/assim_model_mod.f90
options_mod.o: ../../../../assimilation_code/modules/utilities/options_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/options_mod.f90
default_location_mod.o: ../../../../assimilation_code/location/utilities/default_location_mod.f90 types_mod.o utilities_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/location/utilities/default_location_mod.f90
mpi_readobs.o: ../../../../observations/obs_converters/jedi2DART/enkf/mpi_readobs.f90 kinds.o radinfo.o readconvobs.o readsatobs.o readozobs.o mpisetup.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/mpi_readobs.f90
obs_kind_mod.o: ../../../../assimilation_code/modules/observations/obs_kind_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/observations/obs_kind_mod.f90
netcdf_mod.o: ../../../../observations/obs_converters/jedi2DART/enkf/netcdf_mod.f90 mpisetup.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/obs_converters/jedi2DART/enkf/netcdf_mod.f90
distributed_state_mod.o: ../../../../assimilation_code/modules/utilities/distributed_state_mod.f90 mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o no_cray_win_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/distributed_state_mod.f90
utilities_mod.o: ../../../../assimilation_code/modules/utilities/utilities_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../assimilation_code/modules/utilities/utilities_mod.f90
obs_def_mod.o: ../../../../observations/forward_operators/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o assim_model_mod.o obs_kind_mod.o ensemble_manager_mod.o obs_def_utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../observations/forward_operators/obs_def_mod.f90
./mpisetup.f90: ../../../../observations/obs_converters/jedi2DART/mpisetup.f90
	cp ../../../../observations/obs_converters/jedi2DART/mpisetup.f90 .
./sort_mod.f90: ../../../../assimilation_code/modules/utilities/sort_mod.f90
	cp ../../../../assimilation_code/modules/utilities/sort_mod.f90 .
./kinds.f90: ../../../../observations/obs_converters/jedi2DART/enkf/kinds.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/kinds.f90 .
./location_io_mod.f90: ../../../../assimilation_code/location/utilities/location_io_mod.f90
	cp ../../../../assimilation_code/location/utilities/location_io_mod.f90 .
./model_mod.f90: ../../../../models/wrf/model_mod.f90
	cp ../../../../models/wrf/model_mod.f90 .
./obs_sequence_mod.f90: ../../../../assimilation_code/modules/observations/obs_sequence_mod.f90
	cp ../../../../assimilation_code/modules/observations/obs_sequence_mod.f90 .
./random_seq_mod.f90: ../../../../assimilation_code/modules/utilities/random_seq_mod.f90
	cp ../../../../assimilation_code/modules/utilities/random_seq_mod.f90 .
./dart_time_io_mod.f90: ../../../../assimilation_code/modules/io/dart_time_io_mod.f90
	cp ../../../../assimilation_code/modules/io/dart_time_io_mod.f90 .
./obs_utilities_mod.f90: ../../../../observations/obs_converters/utilities/obs_utilities_mod.f90
	cp ../../../../observations/obs_converters/utilities/obs_utilities_mod.f90 .
./params.f90: ../../../../observations/obs_converters/jedi2DART/params.f90
	cp ../../../../observations/obs_converters/jedi2DART/params.f90 .
./state_structure_mod.f90: ../../../../assimilation_code/modules/io/state_structure_mod.f90
	cp ../../../../assimilation_code/modules/io/state_structure_mod.f90 .
./utilities_mod.f90: ../../../../assimilation_code/modules/utilities/utilities_mod.f90
	cp ../../../../assimilation_code/modules/utilities/utilities_mod.f90 .
./readozobs.f90: ../../../../observations/obs_converters/jedi2DART/enkf/readozobs.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/readozobs.f90 .
./default_location_mod.f90: ../../../../assimilation_code/location/utilities/default_location_mod.f90
	cp ../../../../assimilation_code/location/utilities/default_location_mod.f90 .
./netcdf_mod.f90: ../../../../observations/obs_converters/jedi2DART/enkf/netcdf_mod.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/netcdf_mod.f90 .
./radinfo.f90: ../../../../observations/obs_converters/jedi2DART/radinfo.f90
	cp ../../../../observations/obs_converters/jedi2DART/radinfo.f90 .
./mpi_readobs.f90: ../../../../observations/obs_converters/jedi2DART/enkf/mpi_readobs.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/mpi_readobs.f90 .
./readconvobs.f90: ../../../../observations/obs_converters/jedi2DART/enkf/readconvobs.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/readconvobs.f90 .
./assim_model_mod.f90: ../../../../assimilation_code/modules/assimilation/assim_model_mod.f90
	cp ../../../../assimilation_code/modules/assimilation/assim_model_mod.f90 .
./readsatobs.f90: ../../../../observations/obs_converters/jedi2DART/enkf/readsatobs.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/readsatobs.f90 .
./constants.f90: ../../../../observations/obs_converters/jedi2DART/enkf/constants.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/constants.f90 .
./types_mod.f90: ../../../../assimilation_code/modules/utilities/types_mod.f90
	cp ../../../../assimilation_code/modules/utilities/types_mod.f90 .
./jedi_to_dart.f90: ../../../../observations/obs_converters/jedi2DART/jedi_to_dart.f90
	cp ../../../../observations/obs_converters/jedi2DART/jedi_to_dart.f90 .
./netcdf_utilities_mod.f90: ../../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90
	cp ../../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 .
./time_manager_mod.f90: ../../../../assimilation_code/modules/utilities/time_manager_mod.f90
	cp ../../../../assimilation_code/modules/utilities/time_manager_mod.f90 .
./mpi_utilities_mod.f90: ../../../../assimilation_code/modules/utilities/mpi_utilities_mod.f90
	cp ../../../../assimilation_code/modules/utilities/mpi_utilities_mod.f90 .
./ensemble_manager_mod.f90: ../../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90
	cp ../../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 .
./location_mod.f90: ../../../../assimilation_code/location/threed_sphere/location_mod.f90
	cp ../../../../assimilation_code/location/threed_sphere/location_mod.f90 .
./dart_obs_seq_mod.f90: ../../../../observations/obs_converters/jedi2DART/dart_obs_seq_mod.f90
	cp ../../../../observations/obs_converters/jedi2DART/dart_obs_seq_mod.f90 .
./default_model_mod.f90: ../../../../models/utilities/default_model_mod.f90
	cp ../../../../models/utilities/default_model_mod.f90 .
./obs_def_mod.f90: ../../../../observations/forward_operators/obs_def_mod.f90
	cp ../../../../observations/forward_operators/obs_def_mod.f90 .
./module_map_utils.f90: ../../../../models/wrf/module_map_utils.f90
	cp ../../../../models/wrf/module_map_utils.f90 .
./jedi_time_routines_mod.f90: ../../../../observations/obs_converters/jedi2DART/enkf/jedi_time_routines_mod.f90
	cp ../../../../observations/obs_converters/jedi2DART/enkf/jedi_time_routines_mod.f90 .
./no_cray_win_mod.f90: ../../../../assimilation_code/modules/utilities/no_cray_win_mod.f90
	cp ../../../../assimilation_code/modules/utilities/no_cray_win_mod.f90 .
./obs_def_utilities_mod.f90: ../../../../observations/forward_operators/obs_def_utilities_mod.f90
	cp ../../../../observations/forward_operators/obs_def_utilities_mod.f90 .
./distributed_state_mod.f90: ../../../../assimilation_code/modules/utilities/distributed_state_mod.f90
	cp ../../../../assimilation_code/modules/utilities/distributed_state_mod.f90 .
./obs_kind_mod.f90: ../../../../assimilation_code/modules/observations/obs_kind_mod.f90
	cp ../../../../assimilation_code/modules/observations/obs_kind_mod.f90 .
./options_mod.f90: ../../../../assimilation_code/modules/utilities/options_mod.f90
	cp ../../../../assimilation_code/modules/utilities/options_mod.f90 .
SRC = ../../../../models/wrf/module_map_utils.f90 ../../../../assimilation_code/modules/utilities/mpi_utilities_mod.f90 ../../../../assimilation_code/modules/observations/obs_sequence_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/constants.f90 ../../../../observations/obs_converters/jedi2DART/params.f90 ../../../../assimilation_code/modules/utilities/time_manager_mod.f90 ../../../../assimilation_code/modules/io/dart_time_io_mod.f90 ../../../../observations/obs_converters/jedi2DART/radinfo.f90 ../../../../observations/obs_converters/jedi2DART/dart_obs_seq_mod.f90 ../../../../assimilation_code/location/threed_sphere/location_mod.f90 ../../../../assimilation_code/modules/utilities/sort_mod.f90 ../../../../observations/obs_converters/jedi2DART/jedi_to_dart.f90 ../../../../observations/obs_converters/jedi2DART/mpisetup.f90 ../../../../observations/obs_converters/utilities/obs_utilities_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/kinds.f90 ../../../../assimilation_code/modules/utilities/random_seq_mod.f90 ../../../../assimilation_code/modules/utilities/types_mod.f90 ../../../../models/utilities/default_model_mod.f90 ../../../../assimilation_code/modules/utilities/no_cray_win_mod.f90 ../../../../assimilation_code/modules/io/state_structure_mod.f90 ../../../../observations/forward_operators/obs_def_utilities_mod.f90 ../../../../models/wrf/model_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/readsatobs.f90 ../../../../observations/obs_converters/jedi2DART/enkf/jedi_time_routines_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/readconvobs.f90 ../../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 ../../../../assimilation_code/location/utilities/location_io_mod.f90 ../../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/readozobs.f90 ../../../../assimilation_code/modules/assimilation/assim_model_mod.f90 ../../../../assimilation_code/modules/utilities/options_mod.f90 ../../../../assimilation_code/location/utilities/default_location_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/mpi_readobs.f90 ../../../../assimilation_code/modules/observations/obs_kind_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/netcdf_mod.f90 ../../../../assimilation_code/modules/utilities/distributed_state_mod.f90 ../../../../assimilation_code/modules/utilities/utilities_mod.f90 ../../../../observations/forward_operators/obs_def_mod.f90
OBJ = module_map_utils.o mpi_utilities_mod.o obs_sequence_mod.o constants.o params.o time_manager_mod.o dart_time_io_mod.o radinfo.o dart_obs_seq_mod.o location_mod.o sort_mod.o jedi_to_dart.o mpisetup.o obs_utilities_mod.o kinds.o random_seq_mod.o types_mod.o default_model_mod.o no_cray_win_mod.o state_structure_mod.o obs_def_utilities_mod.o model_mod.o readsatobs.o jedi_time_routines_mod.o readconvobs.o netcdf_utilities_mod.o location_io_mod.o ensemble_manager_mod.o readozobs.o assim_model_mod.o options_mod.o default_location_mod.o mpi_readobs.o obs_kind_mod.o netcdf_mod.o distributed_state_mod.o utilities_mod.o obs_def_mod.o
OFF = ../../../../observations/obs_converters/jedi2DART/mpisetup.f90 ../../../../assimilation_code/modules/utilities/sort_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/kinds.f90 ../../../../assimilation_code/location/utilities/location_io_mod.f90 ../../../../models/wrf/model_mod.f90 ../../../../assimilation_code/modules/observations/obs_sequence_mod.f90 ../../../../assimilation_code/modules/utilities/random_seq_mod.f90 ../../../../assimilation_code/modules/io/dart_time_io_mod.f90 ../../../../observations/obs_converters/utilities/obs_utilities_mod.f90 ../../../../observations/obs_converters/jedi2DART/params.f90 ../../../../assimilation_code/modules/io/state_structure_mod.f90 ../../../../assimilation_code/modules/utilities/utilities_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/readozobs.f90 ../../../../assimilation_code/location/utilities/default_location_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/netcdf_mod.f90 ../../../../observations/obs_converters/jedi2DART/radinfo.f90 ../../../../observations/obs_converters/jedi2DART/enkf/mpi_readobs.f90 ../../../../observations/obs_converters/jedi2DART/enkf/readconvobs.f90 ../../../../assimilation_code/modules/assimilation/assim_model_mod.f90 ../../../../observations/obs_converters/jedi2DART/enkf/readsatobs.f90 ../../../../observations/obs_converters/jedi2DART/enkf/constants.f90 ../../../../assimilation_code/modules/utilities/types_mod.f90 ../../../../observations/obs_converters/jedi2DART/jedi_to_dart.f90 ../../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 ../../../../assimilation_code/modules/utilities/time_manager_mod.f90 ../../../../assimilation_code/modules/utilities/mpi_utilities_mod.f90 ../../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 ../../../../assimilation_code/location/threed_sphere/location_mod.f90 ../../../../observations/obs_converters/jedi2DART/dart_obs_seq_mod.f90 ../../../../models/utilities/default_model_mod.f90 ../../../../observations/forward_operators/obs_def_mod.f90 ../../../../models/wrf/module_map_utils.f90 ../../../../observations/obs_converters/jedi2DART/enkf/jedi_time_routines_mod.f90 ../../../../assimilation_code/modules/utilities/no_cray_win_mod.f90 ../../../../observations/forward_operators/obs_def_utilities_mod.f90 ../../../../assimilation_code/modules/utilities/distributed_state_mod.f90 ../../../../assimilation_code/modules/observations/obs_kind_mod.f90 ../../../../assimilation_code/modules/utilities/options_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) jedi_to_dart
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
jedi_to_dart: $(OBJ)
	$(MPILD) $(OBJ) -o jedi_to_dart  $(LDFLAGS)
