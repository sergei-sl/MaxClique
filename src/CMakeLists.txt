set(SOURCES
	main.cpp
	graph.h
	graph.cpp
	utils.h
	simple_solver.h
	simple_solver.cpp
	cliquer_wrapper.h
	cliquer_wrapper.cpp
)

find_package(Threads REQUIRED)

add_executable(SlashchininMaxClique ${SOURCES})

target_link_libraries(SlashchininMaxClique cliquer ${CMAKE_THREAD_LIBS_INIT})
