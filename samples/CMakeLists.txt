project(interior_df)
cmake_minimum_required(VERSION 3.0)

add_executable(diffMapTest diffMapsTest.cpp )
target_include_directories(diffMapTest PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(diffMapTest IDF)

add_executable(idfClassTest idfClassTest.cpp )
target_include_directories(idfClassTest PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(idfClassTest IDF)
