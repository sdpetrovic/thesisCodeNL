 #    Copyright (c) 2010-2016, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      120202    K. Kumar          Adapted to work with Earth orbiting satellite example.
 #      120523    A. Ronse          Updated terminology suiting application. Included script
 #                                  to find+load TudatCoreEnvironment.cmake.
 #      120205    K. Kumar          Removed obsolete TudatCoreEnvironment.cmake call; updated
 #                                  Tudat Core and Tudat minimum requirements; updated license;
 #                                  added static library so body.h shows up in project tree.
 #      130222    K. Kumar          Updated directory definitions; added unit test macro.
 #      130225    K. Kumar          Updated macros; updated set-function calls to check if
 #                                  variables have been set already; updated static library name
 #                                  and naming of applications.
 #     130326    K. Kumar           Added BUILD_DOCUMENTATION option for Doxygen html output
 #                                  generation; removed -O3 flag since it is set by "Release"
 #                                  build type if selected; updated install steps; re-organized
 #                                  file to isolate user-specified section.
 #      160401   S.D. Petrovic      Estimate of when file was created
 #      160428   S.D. Petrovic      Added own libraries to include self made header and source files without it being a class
 #
 #    References
 #
 #    Notes
 #

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(testMainFile_verification)

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(${CMAKE_MODULE_PATH}/compiler.cmake)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIB_ROOT)
  set(LIB_ROOT "${CODEROOT}/lib")
endif()
if(NOT BIN_ROOT)
  set(BIN_ROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BIN_ROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIB_ROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIB_ROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BIN_ROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Define the install targets to create a distribution.
if(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)
    set(TUDAT_BUNDLE_DISTRIBUTION_PATH "${CODEROOT}")
endif(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.45.0
             COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

# Set whether to use the SPICE library integration with Tudat or not. If it not supplied by the
# user (either directly as an argument or through the "UserSettings.txt" file, the default setting
# is "OFF". NB: It is important this section is placed AFTER the "UserSettings.txt" file is
# included, to ensure that that the correct indiciation is given for whether SPICE should be
# included or not has been selected.

# Add an option to toggle the inclusion of SPICE.
# Find CSPICE library on local system.
find_package(Spice)

# Include CSpice directories.
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
endif( )


############ Own additions #############

## Find PaGMO library on local system.
#find_package(PaGMO 1.1.5)

#if(PAGMO_VERSION_OK)
#message("PaGMO found")
#else(PAGMO_VERSION_OK)
#message("PaGMO not found!!")
#endif(PAGMO_VERSION_OK)

#message("Test")

## Include PaGMO directories.
#if(NOT APPLE)
#  include_directories(SYSTEM AFTER "${PaGMO_INCLUDE_DIR}")
#else( )
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PaGMO_INCLUDE_DIR}\"")
#endif( )

## PaGMO lib source files.
#SET(PAGMO_SOURCES
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/base.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/archipelago.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/base_island.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/island.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/population.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sea.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/de.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/jde.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/mde_pbx.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/de_1220.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pso.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pso_generational.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pso_generational_racing.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/bee_colony.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/firefly.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cmaes.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/inverover.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nn_tsp.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nsga2.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/moea_d.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sms_emoa.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cs.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ihs.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sga.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sga_gray.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/monte_carlo.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_co_evolution.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/null.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sa_corana.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/mbh.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ms.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/vega.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_self_adaptive.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pade.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_immune_system.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_core.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nspso.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/spea2.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/base.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/base_r_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/base_s_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/best_s_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/random_s_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/best_kill_s_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_greedy_s_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_best_s_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/fair_r_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/worst_r_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/random_r_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_greedy_r_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_fair_r_policy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_tsp.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_stochastic.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_dtlz.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_unc_mo.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/griewank.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/levy5.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/lennard_jones.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/lavor_maculan.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/ackley.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/branin.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/dtlz.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/golomb_ruler.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/identity.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/inventory.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/luksan_vlcek_1.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/luksan_vlcek_2.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/luksan_vlcek_3.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/himmelblau.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/bukin_f6.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rastrigin.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rosenbrock.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/schwefel.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/snopt_toyprob.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/string_match.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp_cs.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp_vrplc.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/michalewicz.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/dejong.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/sch.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/fon.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/pol.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/kur.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/zdt.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cec2006.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cec2009.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cec2013.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/con2mo.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/pressure_vessel.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/welded_beam.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tens_comp_string.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cstrs_co_evolution.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/death_penalty.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cstrs_self_adaptive.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/antibodies_problem.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/shifted.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/scaled.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rotated.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/normalized.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/decompose.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/noisy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/robust.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/con2uncon.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/barabasi_albert.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/clustered_ba.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/ageing_clustered_ba.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/base.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/custom.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/erdos_renyi.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/fully_connected.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/one_way_ring.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/pan.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/rim.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/ring.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/hypercube.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/unconnected.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/watts_strogatz.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/rng.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hypervolume.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/base.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv2d.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv3d.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/wfg.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/bf_approx.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/bf_fpras.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hoy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv4d_cpp_original/hv.c"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv4d.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/fpl_cpp_original/hv.c"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/fpl.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/racing.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/discrepancy.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/neighbourhood.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/race_pop.cpp"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/race_algo.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cassini_1.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cassini_2.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/gtoc_1.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/gtoc_2.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/sagas.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rosetta.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/messenger.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/messenger_full.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tandem.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp_ds.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/laplace.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/sample_return.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_1dsm_alpha.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_1dsm_tof.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_incipit.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_incipit_cstrs.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_part.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/mga.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/mga_dsm.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/misc4Tandem.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/Pl_Eph_An.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/PowSwingByInv.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/Lambert.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/Astro_Functions.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/time2distance.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/propagateKEP.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/base_gsl.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_bfgs.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_bfgs2.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_derivative_free.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_fr.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_pr.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_gradient.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_nm.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_nm2.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_nm2rand.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/spheres.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/gsl_init.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/base_nlopt.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_cobyla.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_bobyqa.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_sbplx.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_slsqp.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_mma.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_aug_lag.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_aug_lag_eq.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/snopt.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/snopt_cpp_wrapper/snoptProblem_PAGMO.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/snopt_cpp_wrapper/snfilewrapper_PAGMO.c"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ipopt.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ipopt_cpp_wrapper/ipopt_problem.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/worhp.cpp"
###                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/worhp_cpp_wrapper/set_get_params.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/mpi_environment.cpp"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/mpi_island.cpp"
#        )

## PaGMO lib header files.
#SET(PAGMO_HEADERS
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/base.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/archipelago.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/base_island.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/island.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/population.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sea.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/de.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/jde.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/mde_pbx.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/de_1220.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pso.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pso_generational.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pso_generational_racing.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/bee_colony.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/firefly.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cmaes.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/inverover.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nn_tsp.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nsga2.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/moea_d.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sms_emoa.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cs.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ihs.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sga.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sga_gray.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/monte_carlo.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_co_evolution.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/null.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/sa_corana.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/mbh.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ms.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/vega.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_self_adaptive.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/pade.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_immune_system.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/cstrs_core.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nspso.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/spea2.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/base.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/base_r_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/base_s_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/best_s_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/random_s_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/best_kill_s_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_greedy_s_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_best_s_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/fair_r_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/worst_r_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/random_r_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_greedy_r_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/migration/hv_fair_r_policy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_tsp.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_stochastic.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_dtlz.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base_unc_mo.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/griewank.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/levy5.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/lennard_jones.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/lavor_maculan.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/ackley.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/branin.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/dtlz.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/golomb_ruler.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/identity.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/inventory.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/luksan_vlcek_1.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/luksan_vlcek_2.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/luksan_vlcek_3.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/himmelblau.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/bukin_f6.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rastrigin.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rosenbrock.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/schwefel.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/snopt_toyprob.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/string_match.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp_cs.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp_vrplc.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/michalewicz.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/dejong.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/sch.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/fon.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/pol.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/kur.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/zdt.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cec2006.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cec2009.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cec2013.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/con2mo.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/pressure_vessel.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/welded_beam.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tens_comp_string.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cstrs_co_evolution.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/death_penalty.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cstrs_self_adaptive.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/antibodies_problem.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/shifted.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/scaled.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rotated.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/normalized.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/decompose.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/noisy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/robust.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/con2uncon.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/barabasi_albert.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/clustered_ba.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/ageing_clustered_ba.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/base.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/custom.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/erdos_renyi.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/fully_connected.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/one_way_ring.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/pan.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/rim.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/ring.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/hypercube.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/unconnected.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/topology/watts_strogatz.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/rng.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hypervolume.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/base.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv2d.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv3d.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/wfg.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/bf_approx.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/bf_fpras.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hoy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv4d_cpp_original/hv.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/hv4d.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/fpl_cpp_original/hv.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/hv_algorithm/fpl.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/racing.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/discrepancy.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/neighbourhood.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/race_pop.h"
#        "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/util/race_algo.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cassini_1.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/cassini_2.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/gtoc_1.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/gtoc_2.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/sagas.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/rosetta.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/messenger.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/messenger_full.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tandem.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/tsp_ds.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/laplace.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/sample_return.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_1dsm_alpha.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_1dsm_tof.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_incipit.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_incipit_cstrs.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/mga_part.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/mga.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/mga_dsm.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/misc4Tandem.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/Pl_Eph_An.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/PowSwingByInv.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/Lambert.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/Astro_Functions.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/time2distance.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/AstroToolbox/propagateKEP.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/base_gsl.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_bfgs.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_bfgs2.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_derivative_free.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_fr.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_pr.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_gradient.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_nm.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_nm2.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/gsl_nm2rand.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/spheres.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/gsl_init.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/base_nlopt.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_cobyla.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_bobyqa.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_sbplx.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_slsqp.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_mma.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_aug_lag.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/nlopt_aug_lag_eq.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/snopt.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/snopt_cpp_wrapper/snoptProblem_PAGMO.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/snopt_cpp_wrapper/snfilewrapper_PAGMO.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ipopt.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/ipopt_cpp_wrapper/ipopt_problem.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/worhp.h"
###                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/algorithm/worhp_cpp_wrapper/set_get_params.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/mpi_environment.h"
##                "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/mpi_island.h"
#        )


### Add PaGMO files
##set(PAGMO_SOURCES
##    "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/"
##)
##set(PAGMO_HEADERS
##    "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/"
##)


### Add static PaGMO libraries
#add_library(pagmo_static STATIC ${PAGMO_SOURCES} ${PAGMO_HEADERS})
##add_library(pagmo_static STATIC ${PAGMO_LIB_SRC_LIST})

# Add source files.
set(PROJECTLIBRARIES_SOURCES
#  "${SRCROOT}${PROJECTLIBRARIESDIR}/basicRecurrenceRelations.cpp"
    "${SRCROOT}/projectLibraries/basicRecurrenceRelations.cpp"
#    "${SRCROOT}/projectLibraries/allRecurrenceRelations.cpp"
    "${SRCROOT}/projectLibraries/allRecurrenceRelationsSpherical.cpp"
    "${SRCROOT}/projectLibraries/allRecurrenceRelationsCartesian.cpp"
#    "${SRCROOT}/projectLibraries/TaylorSeriesIntegration.cpp"
    "${SRCROOT}/projectLibraries/TaylorSeriesIntegrationSpherical.cpp"
    "${SRCROOT}/projectLibraries/TaylorSeriesIntegrationCartesian.cpp"
    "${SRCROOT}/projectLibraries/otherRequiredFunctions.cpp"
    "${SRCROOT}/projectLibraries/airTemperature.cpp"
    "${SRCROOT}/projectLibraries/airDensity.cpp"
    "${SRCROOT}/projectLibraries/ascentDragForce.cpp"
    "${SRCROOT}/projectLibraries/dragCoefficient.cpp"
    "${SRCROOT}/projectLibraries/trajectoryIntegration.cpp"
    "${SRCROOT}/projectLibraries/mavAscent.cpp"
#    "${SRCROOT}/projectLibraries/ascentStateDerivativeFunction.cpp"
#    "${SRCROOT}/projectLibraries/referenceFrameTransformationsUpdated.cpp"  # Does not quite work yet, so leaving it out for now
)

# Add header files.
set(PROJECTLIBRARIES_HEADERS
#  "${SRCROOT}${PROJECTLIBRARIESDIR}/basicRecurrenceRelations.h"
    "${SRCROOT}/projectLibraries/basicRecurrenceRelations.h"
#    "${SRCROOT}/projectLibraries/allRecurrenceRelations.h"
    "${SRCROOT}/projectLibraries/allRecurrenceRelationsSpherical.h"
    "${SRCROOT}/projectLibraries/allRecurrenceRelationsCartesian.h"
#    "${SRCROOT}/projectLibraries/TaylorSeriesIntegration.h"
    "${SRCROOT}/projectLibraries/TaylorSeriesIntegrationSpherical.h"
    "${SRCROOT}/projectLibraries/TaylorSeriesIntegrationCartesian.h"
    "${SRCROOT}/projectLibraries/otherRequiredFunctions.h"
    "${SRCROOT}/projectLibraries/airTemperature.h"
    "${SRCROOT}/projectLibraries/airDensity.h"
    "${SRCROOT}/projectLibraries/ascentDragForce.h"
    "${SRCROOT}/projectLibraries/dragCoefficient.h"
    "${SRCROOT}/projectLibraries/trajectoryIntegration.h"
    "${SRCROOT}/projectLibraries/mavAscent.h"
#    "${SRCROOT}/projectLibraries/ascentStateDerivativeFunction.h"
#    "${SRCROOT}/projectLibraries/referenceFrameTransformationsUpdated.h"   # Does not quite work yet, so leaving it out for now
)

## Include the directory to the binary tree
#include_directories("${SRCROOT}${PROJECTLIBRARIESDIR}")

# Add static libraries
add_library(project_libraries STATIC ${PROJECTLIBRARIES_SOURCES} ${PROJECTLIBRARIES_HEADERS})

# Add application.
#add_executable(simple_verification_testMainFileClasses "${SRCROOT}/testMainFileClasses.cpp")
#setup_executable_target(simple_verification_testMainFileClasses "${SRCROOT}")
#target_link_libraries(simple_verification_testMainFileClasses tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(simple_verification_testMainFileClassesAndHeaderAndSource "${SRCROOT}/testMainFileClassesAndHeaderAndSource.cpp")
#setup_executable_target(simple_verification_testMainFileClassesAndHeaderAndSource "${SRCROOT}")
#target_link_libraries(simple_verification_testMainFileClassesAndHeaderAndSource project_libraries tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(simple_verification_testMainFileTaylorSeriesIntegration "${SRCROOT}/testMainFileTaylorSeriesIntegration.cpp")
#setup_executable_target(simple_verification_testMainFileTaylorSeriesIntegration "${SRCROOT}")
#target_link_libraries(simple_verification_testMainFileTaylorSeriesIntegration project_libraries tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(simple_verification_testMainFileStateDerivativeFunction "${SRCROOT}/testMainFileStateDerivativeFunction.cpp")
#setup_executable_target(simple_verification_testMainFileStateDerivativeFunction "${SRCROOT}")
#target_link_libraries(simple_verification_testMainFileStateDerivativeFunction project_libraries tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(simple_verification_testMainFileRKF "${SRCROOT}/testMainFileRKF.cpp")
#setup_executable_target(simple_verification_testMainFileRKF "${SRCROOT}")
#target_link_libraries(simple_verification_testMainFileRKF project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

#add_executable(simple_verification_testMainFileFullIntegration "${SRCROOT}/testMainFileFullIntegration.cpp")
#setup_executable_target(simple_verification_testMainFileFullIntegration "${SRCROOT}")
#target_link_libraries(simple_verification_testMainFileFullIntegration project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(simple_verification_testMainFileFullIntegrationSpherical "${SRCROOT}/testMainFileFullIntegrationSpherical.cpp")
setup_executable_target(simple_verification_testMainFileFullIntegrationSpherical "${SRCROOT}")
target_link_libraries(simple_verification_testMainFileFullIntegrationSpherical project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(simple_verification_testMainFileFullIntegrationCartesian "${SRCROOT}/testMainFileFullIntegrationCartesian.cpp")
setup_executable_target(simple_verification_testMainFileFullIntegrationCartesian "${SRCROOT}")
target_link_libraries(simple_verification_testMainFileFullIntegrationCartesian project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )


add_executable(propagation_tool "${SRCROOT}/propagationTool.cpp")
setup_executable_target(propagation_tool "${SRCROOT}")
target_link_libraries(propagation_tool project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )


#add_executable(trajectory_optimisation "${SRCROOT}/trajectoryOptimisation.cpp")
#setup_executable_target(trajectory_optimisation "${SRCROOT}")
#target_link_libraries(trajectory_optimisation project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )
##${MANDATORY_LIBRARIES} pagmo_static


add_executable(validation_thrust_curve "${SRCROOT}/validationThrustCurve.cpp")
setup_executable_target(validation_thrust_curve "${SRCROOT}")
target_link_libraries(validation_thrust_curve project_libraries tudat_numerical_integrators tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

