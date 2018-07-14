function(glmCreateTestGTC NAME)
	set(SAMPLE_NAME test-${NAME})
	add_executable(${SAMPLE_NAME} ${NAME}.cpp)
	add_custom_command(
		TARGET ${SAMPLE_NAME} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data)
	add_test( 
		NAME ${SAMPLE_NAME}
		COMMAND $<TARGET_FILE:${SAMPLE_NAME}> )
	add_dependencies(${SAMPLE_NAME} ${TARGET_DATA_COPY})
endfunction()

add_subdirectory(bug)
add_subdirectory(core)
add_subdirectory(gtc)
add_subdirectory(gtx)
