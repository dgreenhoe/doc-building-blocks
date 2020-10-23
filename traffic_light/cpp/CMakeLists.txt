#==============================================================================
# CMake Script file
# Daniel J. Greenhoe
#==============================================================================

#--------------------------------------
# Build Environment Setup
#--------------------------------------
 cmake_minimum_required( VERSION 3.17 )
 set( CMAKE_VERBOSE_MAKEFILE   $<BOOL:true>                )
 set( CMAKE_C_COMPILER_WORKS   $<BOOL:true>                )
 set( CMAKE_CXX_COMPILER_WORKS $<BOOL:true>                )
 set( CMAKE_SYSTEM_NAME        CMAKE_HOST_SYSTEM_NAME      )
 set( CMAKE_SYSTEM_PROCESSOR   CMAKE_HOST_SYSTEM_PROCESSOR )
 set( CMAKE_CROSSCOMPILING     $<BOOL:false>               )
 set( CMAKE_C_STANDARD         99                          )
 set( CMAKE_CXX_STANDARD       11                          )

#--------------------------------------
# Project
#--------------------------------------
 project( traffic_light
          VERSION     0
          DESCRIPTION "traffic light"
         #HOMEPAGE    https://github.com/dgreenhoe/doc-building-blocks/tree/master/traffic_light
          LANGUAGES   CXX
        )
 set( TARGET ${PROJECT_NAME} )

#--------------------------------------
# Set Directory locations
#--------------------------------------
 set( DIR_GNU C:/CVSSource/Libraries/GNU/mingw-w64/i686-8.1.0-posix-dwarf-rt_v6-rev0/mingw32 )
 set( DIR_BCC C:/CVSSource/Libraries/bcc102 )
 set( DIR_TOOLCHAIN ${DIR_GNU} )
 set( DIR_TOOLCHAIN_INCLUDES  ${DIR_TOOLCHAIN}/include/windows/crtl )
 set( CMAKE_CXX_COMPILER ${DIR_TOOLCHAIN}/bin/g++.exe )

#--------------------------------------
# Set Flags
#--------------------------------------
 set( FLAGS_COMPILE
      #-v
    )

 set( FLAGS_LINK
      #-v
    )

#--------------------------------------
# Source Files
#--------------------------------------
 set( SOURCES_TARGET
      main.cpp
    )

 set( SOURCES_CMD
      command.cpp
      help.cpp
    )

 set( SOURCES_DSP
      estimate_amplitude.cpp
      estimate_induction.cpp
      sequencer.cpp
      traffic.cpp
    )

 set( HEADERS_TARGET
      main.h
    )

 set( HEADERS_CMD
      command.h
      estimate_amplitude.h
      estimate_induction.h
      help.h
      sequencer.h
      traffic.h
    )

 set( HEADERS_DSP
      estimate_amplitude.h
      estimate_induction.h
      sequencer.h
      traffic.h
    )

#--------------------------------------
# Builds
#--------------------------------------
 add_library(                Lib_Cmd   STATIC                                         )
 set_target_properties(      Lib_Cmd   PROPERTIES  PREFIX ""                          )
 target_compile_options(     Lib_Cmd   PRIVATE     ${FLAGS_COMPILE}                   )
 target_sources(             Lib_Cmd   PRIVATE     ${SOURCES_CMD}                     )
 target_include_directories( Lib_Cmd   PRIVATE     ${DIR_TOOLCHAIN_INCLUDES} ..       )

 add_library(                Lib_DSP   STATIC                                         )
 set_target_properties(      Lib_DSP   PROPERTIES PREFIX ""                           )
 target_compile_options(     Lib_DSP   PRIVATE    ${FLAGS_COMPILE}                    )
 target_sources(             Lib_DSP   PRIVATE    ${SOURCES_DSP}                      )
 target_include_directories( Lib_DSP   PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ..        )

 add_executable(             ${TARGET}                                                )
 set_target_properties(      ${TARGET} PROPERTIES OUTPUT_NAME ${TARGET} SUFFIX ".exe" )
 target_compile_options(     ${TARGET} PRIVATE    ${FLAGS_COMPILE}                    )
 target_link_options(        ${TARGET} PRIVATE    ${FLAGS_LINK}                       )
 target_link_libraries(      ${TARGET} PRIVATE    Lib_Cmd Lib_DSP                     )
 target_sources(             ${TARGET} PRIVATE    ${SOURCES_TARGET}                   )
 target_include_directories( ${TARGET} PRIVATE    ${DIR_TOOLCHAIN_INCLUDES} ..        )
