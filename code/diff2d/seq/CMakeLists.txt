cmake_minimum_required( VERSION 3.20 )
project( seq VERSION 1.0 )

add_executable( ${PROJECT_NAME} )
target_compile_features( ${PROJECT_NAME} PRIVATE cxx_std_23 )
target_include_directories( ${PROJECT_NAME} PRIVATE ../lib )
target_sources( ${PROJECT_NAME} PRIVATE diff2d.cpp ../lib/options.cpp ../lib/seq.cpp )
target_sources( ${PROJECT_NAME} INTERFACE ../main.cpp )

# find_package( mdspan REQUIRED )
target_include_directories( ${PROJECT_NAME} PRIVATE ${MDSPAN_INC} )

install( TARGETS ${PROJECT_NAME} DESTINATION . )
