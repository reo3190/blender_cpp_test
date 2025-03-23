# define our project name
project(AsPtrExample)

# Python のリリース版を明示的に検索
find_package(Python3 REQUIRED COMPONENTS Development)

# Python のヘッダーとライブラリを適切にセット
set(PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
set(PYTHON_LIBRARIES ${Python3_LIBRARIES})

# define the include directories our project may need
include_directories(
  ${PYTHON_INCLUDE_DIRS}  # ここを修正
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/source/blender/makesdna
  ${CMAKE_SOURCE_DIR}/source/blender/makesrna
  ${CMAKE_SOURCE_DIR}/source/blender/blenkernel
  ${CMAKE_SOURCE_DIR}/source/blender/blenlib
  ${CMAKE_SOURCE_DIR}/source/blender/bmesh
  ${CMAKE_SOURCE_DIR}/source/blender/imbuf
  ${CMAKE_SOURCE_DIR}/source/blender/depsgraph
  ${CMAKE_SOURCE_DIR}/source/blender/windowmanager
  ${CMAKE_SOURCE_DIR}/source/blender/editors
  ${CMAKE_SOURCE_DIR}/intern/guardedalloc
)

# definitions needed for buildinfo's..
add_definitions(-DBUILD_DATE)
add_definitions(-DWITH_BUILDINFO_HEADER)

# define our main cpp file to compile
set(SRC
  main.cpp
  buildinfo.c
)

add_library(makepyd_AsPtrExample SHARED ${SRC}) # needed to build a .pyd
set_target_properties(makepyd_AsPtrExample PROPERTIES PREFIX "") # needed to build a .pyd
set_target_properties(makepyd_AsPtrExample PROPERTIES OUTPUT_NAME "readmem" SUFFIX ".pyd") # .pyd拡張子を明示的に指定

# link libraries our project may need, may vary depending on the module we include
target_link_libraries(
  makepyd_AsPtrExample
  ${PYTHON_LIBRARIES}  # ここを修正
  bf_python
  # bf_collada
  bf_blenlib
  bf_dna
  bf_imbuf
  bf_editor_io
  bf_editor_scene
  bf_editor_object
  bf_editor_mesh
  bf_blenloader
  bf_blenkernel
  dbghelp
  Ws2_32
  Comctl32
  imm32
  Shcore
  vfw32
  version
  Dwmapi
)

# define destination folder, will be located in "buildRelease/bin/release/bakedpyd_AsPtrExample/"
install(TARGETS makepyd_AsPtrExample DESTINATION ${CMAKE_INSTALL_PREFIX}/bakedpyd_AsPtrExample)

# Pythonライブラリのパスを確認
message(STATUS "Using Python include dirs: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Using Python libraries: ${PYTHON_LIBRARIES}")
