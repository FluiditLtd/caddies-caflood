########### ADD CAFLOOD application ###############

# Application that performs a 2D flood modelling using cellular
# automata model(s) and it is implemented using the CAAPI.

### Set up the application configuration

# The exe name of the app.
set(CAAPI_APP_EXE  "caflood")

# The main source of the app.
set(CAAPI_APP_MAIN "main.cpp")

# The source of the app.
set(CAAPI_APP_SOURCES "")
file(GLOB CAAPI_APP_SOURCES "*.cpp")

# The headers of the app.
set(CAAPI_APP_HEADERS "")
file(GLOB CAAPI_APP_HEADERS "*.hpp")

# The data 
set(CAAPI_APP_DATA "NEWS;LICENSE;AUTHORS;Version.txt;CMakeLists.txt")

# Find the CA functions
set(CAAPI_APP_CAFUNCTIONS)
set(CAAPI_IMPL_HEADERS)

# CAFUNCTIONS
CAAPI_add_source_dir(CAfunctions "${CMAKE_CURRENT_SOURCE_DIR}/*.ca" CAAPI_APP_CAFUNCTIONS)

# Execute the module that manage ca for the specific implementation
if(EXISTS "${CAAPI_SPECIFIC_IMPL_DIR}/extra/CAAPI_ca.cmake")

  include("${CAAPI_SPECIFIC_IMPL_DIR}/extra/CAAPI_ca.cmake")

  CAAPI_ca("${CAAPI_SPECIFIC_IMPL_DIR}" CAAPI_IMPL_HEADERS CAAPI_APP_CAFUNCTIONS)
endif()

# Set up the custom target name
set(CAAPI_APP_CUSTOM_TARGET "${CAAPI_APP_EXE}_do_impl_headers")

# Create the target that make the eventual implementation headers
# necessary.
add_custom_target(${CAAPI_APP_CUSTOM_TARGET} DEPENDS 
  ${CAAPI_IMPL_HEADERS} ${CAAPI_HEADERS} ${CAAPI_APP_CAFUNCTIONS})

# Add the executable to the project using the specified source files.
add_executable("${CAAPI_APP_EXE}" "${CAAPI_APP_MAIN}"
  ${CAAPI_APP_SOURCES} ${CAAPI_SOURCES} 
  ${CAAPI_APP_CAFUNCTIONS} ${CAAPI_HEADERS} ${CAAPI_APP_HEADERS})
add_dependencies("${CAAPI_APP_EXE}" ${CAAPI_APP_CUSTOM_TARGET})
target_link_libraries("${CAAPI_APP_EXE}" ${CAAPI_IMPL_LIBRARIES} )

########################## VERSION  ##########################

file(STRINGS "Version.txt" CAAPI_APP_VERSION)

# Show the version number
message(STATUS "[${CAAPI_APP_NAME}]: Version ${CAAPI_APP_VERSION}")

if (WIN32)
  add_definitions( /DCAAPI_APP_VERSION=${CAAPI_APP_VERSION} )
  add_definitions( /DCAAPI_APP_NAME=${CAAPI_APP_NAME} )
else()
  add_definitions( -DCAAPI_APP_VERSION=${CAAPI_APP_VERSION} )
  add_definitions( -DCAAPI_APP_NAME=${CAAPI_APP_NAME} )
endif()


########################## INFO  ##########################

# Number of files
list(LENGTH  CAAPI_APP_CAFUNCTIONS   __num_ca)
list(LENGTH  CAAPI_APP_DATA    __num_data)
list(LENGTH  CAAPI_APP_HEADERS __num_headers)
list(LENGTH  CAAPI_APP_SOURCES __num_sources)
MATH(EXPR __num_sources "${__num_sources}+1")
message(STATUS "[${CAAPI_APP_NAME}]: File #: Source (${__num_sources}) Headers (${__num_headers}) CA (${__num_ca}) Data (${__num_data})")


########################## INSTALL  ##########################

# Set the component name which will influence the archive name.
set(CAAPI_APP_COMPONENT_NAME "app-${CAAPI_APP_NAME}-${CAAPI_APP_VERSION}")

# Add the file that are needed to compile this application into its
# specific component.
install(FILES   
  ${CAAPI_APP_MAIN}
  ${CAAPI_APP_SOURCES}
  ${CAAPI_APP_HEADERS}
  ${CAAPI_APP_DATA}
  ${CAAPI_APP_CAFUNCTIONS}
  DESTINATION apps/caflood
  COMPONENT ${CAAPI_APP_COMPONENT_NAME})


########################## PACKAGE  ##########################

# Add the newly created component
cpack_add_component(${CAAPI_APP_COMPONENT_NAME} DISPLAY_NAME "caflood version ${CAAPI_APP_VERSION}"
  DESCRIPTION "Application that performs a 2D flood modelling using cellular
    automata model(s) and it is implemented using the CAAPI"
  GROUP  apps)

