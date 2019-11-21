CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(REDCMC_SRC utils.c arrays.c fdjac_ext.c fx_fmri.c gx_fmri.c reduce.c integrate.c rInterface.c)
SET(REDCMC_HDR utils.h arrays.h fx_fmri.h gx_fmri.h reduce.h integrate.h gsl_multifit_nlin_ext.h)

ADD_LIBRARY( redcmc SHARED ${REDCMC_SRC} ${REDCMC_HDR} )
TARGET_LINK_LIBRARIES(redcmc ${GSL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(redcmc PUBLIC "${PROJECT_SOURCE_DIR}/src" ${CMAKE_CURRENT_BINARY_DIR} ${GSL_INCLUDE_DIRS})

#ADD_EXECUTABLE( functest functest.c )
#TARGET_LINK_LIBRARIES(functest redcmc ${GSL_LIBRARIES})
#TARGET_INCLUDE_DIRECTORIES(functest PUBLIC "${PROJECT_SOURCE_DIR}/src" ${GSL_INCLUDE_DIRS})

INCLUDE(GenerateExportHeader)
generate_export_header(redcmc) # ${CMAKE_CURRENT_BINARY_DIR}/redcmc_export.h

install(TARGETS redcmc
        DESTINATION lib
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/redcmc_export.h
        DESTINATION lib
)


set (_resolvedFiles "")
foreach (_file ${GSL_LIBRARIES})
    get_filename_component(_resolvedFile "${_file}" REALPATH)
    list (APPEND _resolvedFiles "${_resolvedFile}")
endforeach()
install(FILES ${_resolvedFiles} DESTINATION lib)

#install(FILES ${GSL_LIBRARIES}
#        DESTINATION lib
#)
