# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 

SRCROOT = /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/

CPPDEFS = -Duse_libMPI -Duse_netCDF -DSPMD


include /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../bin/mkmf.template.nyu


.DEFAULT:
	-echo $@ does not exist.
all: fms.x
a2b_edge.o: $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 fv_grid_utils.o fv_grid_tools.o fv_mp_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90
amip_interp.o: $(SRCROOT)shared/amip_interp/amip_interp.F90 time_interp.o time_manager.o horiz_interp.o fms.o fms_io.o constants.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/amip_interp/amip_interp.F90
atmos_model.o: $(SRCROOT)atmos_solo/atmos_model.F90 atmosphere.o time_manager.o fms.o fms_io.o mpp_domains.o mpp_io.o diag_manager.o field_manager.o tracer_manager.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_solo/atmos_model.F90
atmosphere.o: $(SRCROOT)atmos_cubed_sphere/driver/solo/atmosphere.F90 constants.o fms.o time_manager.o mpp_domains.o hswf.o fv_arrays.o fv_control.o fv_phys.o fv_diagnostics.o fv_timing.o fv_restart.o fv_dynamics.o fv_grid_tools.o lin_cloud_microphys.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/driver/solo/atmosphere.F90
axis_utils.o: $(SRCROOT)shared/axis_utils/axis_utils.F90 mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/axis_utils/axis_utils.F90
constants.o: $(SRCROOT)shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/constants/constants.F90
create_xgrid.o: $(SRCROOT)shared/mosaic/create_xgrid.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/create_xgrid.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic	$(SRCROOT)shared/mosaic/create_xgrid.c
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)shared/diag_manager/diag_grid.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/diag_manager/diag_grid.F90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h time_manager.o mpp_io.o fms.o mpp.o diag_axis.o diag_util.o diag_data.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.F90 mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_output.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_util.F90
dyn_core.o: $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 mpp_domains.o mpp_parameter.o fv_mp_mod.o fv_control.o sw_core.o a2b_edge.o nh_core.o fv_grid_tools.o fv_grid_utils.o fv_timing.o fv_diagnostics.o fv_nudge.o test_cases.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90
external_ic.o: $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 fms.o fms_io.o mpp.o mpp_parameter.o mpp_domains.o tracer_manager.o field_manager.o constants.o external_sst.o fv_arrays.o fv_diagnostics.o fv_grid_tools.o fv_grid_utils.o fv_io.o fv_mapz.o fv_mp_mod.o fv_surf_map.o fv_timing.o init_hydro.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90
external_sst.o: $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 amip_interp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)shared/field_manager/fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/field_manager/fm_util.F90
fms.o: $(SRCROOT)shared/fms/fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/fms/fms.F90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_2d.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_3d.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_4d.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms	$(SRCROOT)shared/fms/fms_io.F90
fv_arrays.o: $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90
fv_control.o: $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 constants.o field_manager.o mpp.o mpp_domains.o tracer_manager.o fv_io.o fv_restart.o fv_arrays.o fv_grid_utils.o fv_grid_tools.o fv_mp_mod.o test_cases.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_control.F90
fv_diagnostics.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 constants.o fms_io.o time_manager.o mpp_domains.o diag_manager.o fv_arrays.o fv_mapz.o fv_mp_mod.o fv_eta.o fv_grid_tools.o fv_grid_utils.o a2b_edge.o fv_surf_map.o fv_sg.o tracer_manager.o field_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90
fv_dynamics.o: $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 constants.o dyn_core.o fv_mapz.o fv_tracer2d.o fv_grid_tools.o fv_control.o fv_grid_utils.o fv_mp_mod.o fv_timing.o diag_manager.o fv_diagnostics.o mpp_domains.o field_manager.o tracer_manager.o fv_sg.o tp_core.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90
fv_eta.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 constants.o fv_mp_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90
fv_fill.o: $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
fv_grid_tools.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 constants.o fv_arrays.o fv_grid_utils.o fv_timing.o fv_mp_mod.o sorted_index.o mpp.o mpp_domains.o mpp_io.o mpp_parameter.o fms.o fms_io.o mosaic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90
fv_grid_utils.o: $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_parameter.o external_sst.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90
fv_io.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 fms.o fms_io.o mpp.o mpp_domains.o tracer_manager.o field_manager.o external_sst.o fv_arrays.o fv_mapz.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90
fv_mapz.o: $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 constants.o fv_grid_tools.o fv_grid_utils.o fv_fill.o fv_mp_mod.o mpp_domains.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90
fv_mp_mod.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 fms.o mpp.o mpp_domains.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90
fv_nudge.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 constants.o fms.o fms_io.o mpp.o mpp_domains.o time_manager.o external_sst.o fv_control.o fv_grid_utils.o fv_grid_tools.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90
fv_phys.o: $(SRCROOT)atmos_cubed_sphere/driver/solo/fv_phys.F90 time_manager.o fv_update_phys.o fv_timing.o hswf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/driver/solo/fv_phys.F90
fv_restart.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 constants.o fv_arrays.o fv_io.o fv_grid_tools.o fv_grid_utils.o fv_diagnostics.o init_hydro.o mpp_domains.o mpp.o test_cases.o fv_mp_mod.o fv_surf_map.o tracer_manager.o field_manager.o external_ic.o fv_eta.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90
fv_sg.o: $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 constants.o fv_mp_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90
fv_surf_map.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 fms.o mpp.o mpp_domains.o constants.o fms_io.o fv_grid_utils.o fv_mp_mod.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90
fv_timing.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 fv_mp_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90
fv_tracer2d.o: $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 tp_core.o fv_grid_tools.o fv_grid_utils.o fv_mp_mod.o mpp_domains.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90
fv_update_phys.o: $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 constants.o field_manager.o mpp_domains.o mpp_parameter.o mpp.o time_manager.o tracer_manager.o fv_arrays.o fv_control.o fv_mp_mod.o fv_eta.o fv_grid_utils.o fv_grid_tools.o fv_timing.o fv_diagnostics.o fv_nudge.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90
gradient.o: $(SRCROOT)shared/mosaic/gradient.F90 mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)shared/mosaic/gradient_c2l.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic	$(SRCROOT)shared/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)shared/mosaic/grid.F90 constants.o fms.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/grid.F90
horiz_interp.o: $(SRCROOT)shared/horiz_interp/horiz_interp.F90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 mpp.o fms.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
hswf.o: $(SRCROOT)atmos_cubed_sphere/driver/solo/hswf.F90 constants.o mpp_domains.o time_manager.o diag_manager.o lin_cloud_microphys.o mpp.o fv_grid_tools.o fv_grid_utils.o fv_mp_mod.o fv_diagnostics.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/driver/solo/hswf.F90
init_hydro.o: $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 constants.o fv_grid_utils.o fv_grid_tools.o fv_mp_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90
interp.o: $(SRCROOT)shared/mosaic/interp.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/interp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic	$(SRCROOT)shared/mosaic/interp.c
lin_cloud_microphys.o: $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 mpp.o diag_manager.o time_manager.o constants.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
mosaic.o: $(SRCROOT)shared/mosaic/mosaic.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)shared/mosaic/mosaic_util.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic	$(SRCROOT)shared/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/system_clock.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_a.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_s.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_a.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_s.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_int.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_scalar.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_comm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_define.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_misc.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_update_domains2D.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_update.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_updateV.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_check.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_checkV.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_redistribute.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_get_boundary.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_get_boundary.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_reduce.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_reduce.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum_tl.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_field.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_global_field.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_misc.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_connect.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_read.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_write.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)shared/mpp/mpp_utilities.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/mpp_utilities.F90
nh_core.o: $(SRCROOT)atmos_cubed_sphere/model/nh_core.F90 fms.o constants.o fv_control.o tp_core.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/nh_core.F90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
platform.o: $(SRCROOT)shared/platform/platform.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/platform/platform.F90
read_mosaic.o: $(SRCROOT)shared/mosaic/read_mosaic.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/read_mosaic.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic	$(SRCROOT)shared/mosaic/read_mosaic.c
sim_nc_mod.o: $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
sorted_index.o: $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
sw_core.o: $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 fv_mp_mod.o fv_grid_tools.o tp_core.o fv_grid_utils.o test_cases.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/sw_core.F90
test_cases.o: $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 constants.o init_hydro.o fv_mp_mod.o fv_grid_utils.o fv_surf_map.o fv_grid_tools.o fv_eta.o mpp.o mpp_domains.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90
test_fms_io.o: $(SRCROOT)shared/fms/test_fms_io.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/fms/test_fms_io.F90
test_mpp.o: $(SRCROOT)shared/mpp/test_mpp.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)shared/mpp/test_mpp_domains.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)shared/mpp/test_mpp_io.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)shared/mpp/test_mpp_pset.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include	$(SRCROOT)shared/mpp/test_mpp_pset.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
time_interp.o: $(SRCROOT)shared/time_interp/time_interp.F90 time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_interp/time_interp.F90
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 constants.o fms.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/time_manager/time_manager.F90
tp_core.o: $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 fv_mp_mod.o fv_grid_utils.o fv_grid_tools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/tp_core.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
./mpp_sum_nocomm.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_nocomm.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_nocomm.h .
./external_ic.F90: $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 .
./mpp_io_util.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_util.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_util.inc .
./mpp_do_checkV.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_checkV.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_checkV.h .
./tp_core.F90: $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 .
./mpp_domains_reduce.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_reduce.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_reduce.inc .
./fv_sg.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 .
./fv_grid_utils.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 .
./mpp_chksum_int.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_int.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_int.h .
./mpp_global_reduce.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_reduce.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_reduce.h .
./fv_update_phys.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 .
./fms_platform.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h .
./gradient_c2l.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/gradient_c2l.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/gradient_c2l.h .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./mpp_io_read.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_read.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_read.inc .
./test_mpp_pset.F90: $(SRCROOT)shared/mpp/test_mpp_pset.F90
	cp $(SRCROOT)shared/mpp/test_mpp_pset.F90 .
./mpp_io_write.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_write.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_write.inc .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./fv_grid_tools.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 .
./mpp_sum.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum.inc .
./fv_io.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 .
./mpp_update_domains2D.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_update_domains2D.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_update_domains2D.h .
./mpp_reduce_mpi.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_mpi.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_mpi.h .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./mosaic.F90: $(SRCROOT)shared/mosaic/mosaic.F90
	cp $(SRCROOT)shared/mosaic/mosaic.F90 .
./mpp_util_sma.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_sma.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_sma.inc .
./mpp_io_misc.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_misc.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_misc.inc .
./fv_dynamics.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 .
./mpp_transmit_mpi.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_mpi.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_mpi.h .
./mpp_read_2Ddecomp.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_read_2Ddecomp.h .
./mpp_reduce_sma.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_sma.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_sma.h .
./hswf.F90: $(SRCROOT)atmos_cubed_sphere/driver/solo/hswf.F90
	cp $(SRCROOT)atmos_cubed_sphere/driver/solo/hswf.F90 .
./test_cases.F90: $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 .
./fv_arrays.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 .
./mpp_do_redistribute.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_redistribute.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_redistribute.h .
./fv_timing.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./mpp_do_check.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_check.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_check.h .
./amip_interp.F90: $(SRCROOT)shared/amip_interp/amip_interp.F90
	cp $(SRCROOT)shared/amip_interp/amip_interp.F90 .
./read_data_4d.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_4d.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_4d.inc .
./sw_core.F90: $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 .
./gradient.F90: $(SRCROOT)shared/mosaic/gradient.F90
	cp $(SRCROOT)shared/mosaic/gradient.F90 .
./horiz_interp_type.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 .
./fv_tracer2d.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 .
./fv_surf_map.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 .
./mpp_data_sma.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_sma.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_sma.inc .
./constant.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h .
./mpp_utilities.F90: $(SRCROOT)shared/mpp/mpp_utilities.F90
	cp $(SRCROOT)shared/mpp/mpp_utilities.F90 .
./mpp_domains_comm.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_comm.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_comm.inc .
./axis_utils.F90: $(SRCROOT)shared/axis_utils/axis_utils.F90
	cp $(SRCROOT)shared/axis_utils/axis_utils.F90 .
./read_data_3d.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_3d.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_3d.inc .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./atmosphere.F90: $(SRCROOT)atmos_cubed_sphere/driver/solo/atmosphere.F90
	cp $(SRCROOT)atmos_cubed_sphere/driver/solo/atmosphere.F90 .
./test_mpp_domains.F90: $(SRCROOT)shared/mpp/test_mpp_domains.F90
	cp $(SRCROOT)shared/mpp/test_mpp_domains.F90 .
./mpp_chksum_scalar.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_scalar.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_scalar.h .
./mpp_util_nocomm.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_nocomm.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_nocomm.inc .
./interp.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/interp.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/interp.h .
./diag_axis.F90: $(SRCROOT)shared/diag_manager/diag_axis.F90
	cp $(SRCROOT)shared/diag_manager/diag_axis.F90 .
./mpp_reduce_nocomm.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_nocomm.h .
./mpp_util.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util.inc .
./horiz_interp_conserve.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 .
./dyn_core.F90: $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 .
./mpp_write.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write.h .
./init_hydro.F90: $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 .
./nh_core.F90: $(SRCROOT)atmos_cubed_sphere/model/nh_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/nh_core.F90 .
./mosaic_util.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h .
./mpp_error_s_s.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_s.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_s.h .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./diag_manager.F90: $(SRCROOT)shared/diag_manager/diag_manager.F90
	cp $(SRCROOT)shared/diag_manager/diag_manager.F90 .
./mpp_error_s_a.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_a.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_a.h .
./fv_fill.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 .
./horiz_interp_bicubic.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 .
./fms.F90: $(SRCROOT)shared/fms/fms.F90
	cp $(SRCROOT)shared/fms/fms.F90 .
./mpp_comm_nocomm.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_nocomm.inc .
./read_data_2d.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_2d.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_2d.inc .
./fv_restart.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 .
./mpp_do_update.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_update.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_update.h .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./mpp_sum_sma.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_sma.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_sma.h .
./fm_util.F90: $(SRCROOT)shared/field_manager/fm_util.F90
	cp $(SRCROOT)shared/field_manager/fm_util.F90 .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./sim_nc_mod.F90: $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 .
./mpp_comm_sma.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_sma.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_sma.inc .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./mpp_transmit_nocomm.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_nocomm.h .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./mpp_global_sum_tl.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum_tl.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum_tl.h .
./mpp_comm_mpi.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_mpi.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_mpi.inc .
./mpp_chksum.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum.h .
./mpp_do_global_field.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_global_field.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_global_field.h .
./gradient_c2l.c: $(SRCROOT)shared/mosaic/gradient_c2l.c
	cp $(SRCROOT)shared/mosaic/gradient_c2l.c .
./mpp_do_updateV.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_updateV.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_updateV.h .
./test_mpp_io.F90: $(SRCROOT)shared/mpp/test_mpp_io.F90
	cp $(SRCROOT)shared/mpp/test_mpp_io.F90 .
./diag_grid.F90: $(SRCROOT)shared/diag_manager/diag_grid.F90
	cp $(SRCROOT)shared/diag_manager/diag_grid.F90 .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./system_clock.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/system_clock.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/system_clock.h .
./mpp_global_sum.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum.h .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./a2b_edge.F90: $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 .
./mpp_memutils.F90: $(SRCROOT)shared/mpp/mpp_memutils.F90
	cp $(SRCROOT)shared/mpp/mpp_memutils.F90 .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./read_mosaic.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/read_mosaic.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/read_mosaic.h .
./lin_cloud_microphys.F90: $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
	cp $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 .
./test_fms_io.F90: $(SRCROOT)shared/fms/test_fms_io.F90
	cp $(SRCROOT)shared/fms/test_fms_io.F90 .
./fv_diagnostics.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 .
./horiz_interp_spherical.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./mpp_data_nocomm.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_nocomm.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_nocomm.inc .
./parse.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/field_manager/parse.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/field_manager/parse.inc .
./create_xgrid.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/create_xgrid.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/create_xgrid.h .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./mpp_io_connect.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_connect.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_connect.inc .
./test_mpp.F90: $(SRCROOT)shared/mpp/test_mpp.F90
	cp $(SRCROOT)shared/mpp/test_mpp.F90 .
./mpp_global_field.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_field.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_field.h .
./horiz_interp_bilinear.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 .
./mpp_transmit.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit.inc .
./mpp_error_a_s.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_s.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_s.h .
./mpp_domains_misc.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_misc.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_misc.inc .
./mosaic_util.c: $(SRCROOT)shared/mosaic/mosaic_util.c
	cp $(SRCROOT)shared/mosaic/mosaic_util.c .
./mpp_get_boundary.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_get_boundary.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_get_boundary.h .
./read_mosaic.c: $(SRCROOT)shared/mosaic/read_mosaic.c
	cp $(SRCROOT)shared/mosaic/read_mosaic.c .
./fv_mapz.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 .
./fv_nudge.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 .
./mpp_data_mpi.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_mpi.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_mpi.inc .
./interp.c: $(SRCROOT)shared/mosaic/interp.c
	cp $(SRCROOT)shared/mosaic/interp.c .
./mpp_error_a_a.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_a.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_a.h .
./sorted_index.F90: $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./mpp_write_2Ddecomp.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write_2Ddecomp.h .
./write_data.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/write_data.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/write_data.inc .
./fv_phys.F90: $(SRCROOT)atmos_cubed_sphere/driver/solo/fv_phys.F90
	cp $(SRCROOT)atmos_cubed_sphere/driver/solo/fv_phys.F90 .
./mpp_domains_util.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_util.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_util.inc .
./time_interp.F90: $(SRCROOT)shared/time_interp/time_interp.F90
	cp $(SRCROOT)shared/time_interp/time_interp.F90 .
./fv_eta.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 .
./horiz_interp.F90: $(SRCROOT)shared/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp.F90 .
./diag_output.F90: $(SRCROOT)shared/diag_manager/diag_output.F90
	cp $(SRCROOT)shared/diag_manager/diag_output.F90 .
./fv_control.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 .
./mpp_transmit_sma.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_sma.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_sma.h .
./mpp_comm.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm.inc .
./mpp_domains_define.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_define.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_define.inc .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./grid.F90: $(SRCROOT)shared/mosaic/grid.F90
	cp $(SRCROOT)shared/mosaic/grid.F90 .
./mpp_do_get_boundary.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_get_boundary.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_get_boundary.h .
./atmos_model.F90: $(SRCROOT)atmos_solo/atmos_model.F90
	cp $(SRCROOT)atmos_solo/atmos_model.F90 .
./external_sst.F90: $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 .
./mpp_sum_mpi.h: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_mpi.h
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_mpi.h .
./fv_mp_mod.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 .
./create_xgrid.c: $(SRCROOT)shared/mosaic/create_xgrid.c
	cp $(SRCROOT)shared/mosaic/create_xgrid.c .
./mpp_util_mpi.inc: /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_mpi.inc
	cp /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_mpi.inc .
SRC = $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)atmos_cubed_sphere/driver/solo/atmosphere.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)atmos_cubed_sphere/driver/solo/fv_phys.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/amip_interp/amip_interp.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)atmos_cubed_sphere/driver/solo/hswf.F90 $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)atmos_solo/atmos_model.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)atmos_cubed_sphere/model/nh_core.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_checkV.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_reduce.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_2d.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_update.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_int.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_reduce.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/gradient_c2l.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_read.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_write.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum_tl.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_update_domains2D.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_global_field.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_updateV.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_misc.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/system_clock.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/read_mosaic.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/field_manager/parse.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_redistribute.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/create_xgrid.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_connect.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_check.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_field.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_4d.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_s.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_misc.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_get_boundary.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_a.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_comm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_3d.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_scalar.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/interp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/write_data.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_define.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_get_boundary.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_s.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_a.h
OBJ = mpp_data.o tp_core.o sw_core.o gradient_c2l.o atmosphere.o fms_io.o tracer_manager.o mpp_memutils.o fv_phys.o diag_util.o mosaic_util.o a2b_edge.o external_ic.o time_interp.o mosaic.o memuse.o memutils.o threadloc.o horiz_interp_bicubic.o diag_axis.o diag_manager.o fv_grid_utils.o mpp_utilities.o mpp_pset.o fms.o test_mpp_io.o create_xgrid.o interp.o diag_output.o fv_io.o fv_nudge.o horiz_interp_conserve.o fv_update_phys.o mpp_parameter.o axis_utils.o external_sst.o platform.o amip_interp.o mpp_domains.o fv_dynamics.o diag_grid.o fv_mapz.o test_cases.o time_manager.o fv_sg.o hswf.o sorted_index.o fv_arrays.o atmos_model.o mpp_io.o fv_control.o sim_nc_mod.o fv_tracer2d.o horiz_interp_spherical.o init_hydro.o fv_restart.o fv_timing.o horiz_interp.o gradient.o dyn_core.o fv_diagnostics.o nh_core.o fv_mp_mod.o fm_util.o horiz_interp_bilinear.o diag_data.o grid.o read_mosaic.o constants.o fv_grid_tools.o fv_fill.o test_mpp_pset.o test_fms_io.o fv_eta.o field_manager.o nsclock.o mpp.o fv_surf_map.o horiz_interp_type.o test_mpp_domains.o test_mpp.o lin_cloud_microphys.o
OFF = /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_util.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_checkV.h $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_reduce.inc $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_int.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_reduce.h $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/include/fms_platform.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/gradient_c2l.h $(SRCROOT)shared/mpp/mpp_pset.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_read.inc $(SRCROOT)shared/mpp/test_mpp_pset.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_write.inc $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum.inc $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_update_domains2D.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mosaic/mosaic.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_misc.inc $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_mpi.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)atmos_cubed_sphere/driver/solo/hswf.F90 $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)shared/mpp/mpp_data.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/amip_interp/amip_interp.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_4d.inc $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_sma.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/constant.h $(SRCROOT)shared/mpp/mpp_utilities.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)shared/axis_utils/axis_utils.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_3d.inc $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)atmos_cubed_sphere/driver/solo/atmosphere.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum_scalar.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/interp.h $(SRCROOT)shared/diag_manager/diag_axis.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_reduce_nocomm.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util.inc $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write.h $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)atmos_cubed_sphere/model/nh_core.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/mosaic_util.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_s.h $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_s_a.h $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/fms/fms.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/read_data_2d.inc $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_update.h $(SRCROOT)shared/field_manager/field_manager.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_sma.h $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)shared/mpp/mpp.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)shared/memutils/memuse.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum_tl.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm_mpi.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_chksum.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_global_field.h $(SRCROOT)shared/mosaic/gradient_c2l.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_updateV.h $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/system_clock.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_sum.h $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/read_mosaic.h $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/memutils/memutils.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_nocomm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/field_manager/parse.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mosaic/create_xgrid.h $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_io_connect.inc $(SRCROOT)shared/mpp/test_mpp.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_global_field.h $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_s.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)shared/mosaic/mosaic_util.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_get_boundary.h $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_data_mpi.inc $(SRCROOT)shared/mosaic/interp.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_error_a_a.h $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)shared/mpp/threadloc.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/fms/write_data.inc $(SRCROOT)atmos_cubed_sphere/driver/solo/fv_phys.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_util.inc $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_transmit_sma.h /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_comm.inc /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_domains_define.inc $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)shared/mosaic/grid.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_do_get_boundary.h $(SRCROOT)atmos_solo/atmos_model.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)shared/mosaic/create_xgrid.c /home/ag4680/models/code/riga/cubed_sphere/exp/atmos_cubed_sphere/../../src/shared/mpp/include/mpp_util_mpi.inc
clean: neat
	-rm -f .fms.x.cppdefs $(OBJ) fms.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fms.x: $(OBJ) 
	$(LD) $(OBJ) -o fms.x  $(LDFLAGS)
