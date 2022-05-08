# Remove quotes from input path
string(REPLACE "\"" "" VORO_CMAKE_FILE ${VORO_CMAKE_FILE})

# Actually apply the fix
file(READ ${VORO_CMAKE_FILE} filedata)
string(REPLACE "src/voro++.hh" "src/*.hh" filedata "${filedata}")
file(WRITE ${VORO_CMAKE_FILE} "${filedata}")