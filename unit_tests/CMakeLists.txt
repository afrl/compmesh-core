include("$ENV{SHEFF_TOOLS}/cmake/IncludeSHEFF.cmake")

# Add all test as serial tests
file(GLOB_RECURSE _unit_tests "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")
foreach(_test ${_unit_tests})
  add_serial_test(${MAIN_COMPONENT} ${PUBLIC_MODULE_NAME}::${MAIN_COMPONENT} ${_test})
endforeach()

assign_source_groups_sheff()