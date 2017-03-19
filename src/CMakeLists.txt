set(SRCs
	main.c
	fluid.c
	cg.c
	ui.c
	poissonfft.c
	fileio.c
)

find_package(SDL2 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(FFTW REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GNU_SOURCE -Wall -g -O3 -ffast-math ${OpenMP_C_FLAGS}")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g -O3 -ffast-math")
include_directories(${SDL2_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})

add_executable(nvst2 ${SRCs})


target_link_libraries(nvst2 -lm ${SDL2_LIBRARIES} ${FFTW_LIBRARIES} -lfftw3_omp)
#target_link_libraries(nvst2 -lm ${SDL2_LIBRARIES} ${FFTW_LIBRARIES} -lfftw3)
