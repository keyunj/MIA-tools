CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(FFR)

FIND_PACKAGE(VTK REQUIRED)
FIND_PACKAGE(ITK REQUIRED)

INCLUDE(
   ${VTK_USE_FILE}
   ${ITK_USE_FILE}
)

INCLUDE_DIRECTORIES(
   ${CMAKE_CURRENT_BINARY_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/OpenCV/include
)

LINK_DIRECTORIES(
   ${CMAKE_CURRENT_SOURCE_DIR}/OpenCV/lib
)

ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

ADD_EXECUTABLE(FFR Project.cxx FFR.h FFR.cxx FFRCore.h FFRCore.cxx ExtendIdFilter.h ExtendIdFilter.cxx ExtendTubeFilter.h ExtendTubeFilter.cxx common.h common.cxx Learning.h Learning.cxx LearningImpl.h LearningImpl.cxx psimpl.h itkBinaryThinningImageFilter3D.h itkBinaryThinningImageFilter3D.txx ExtendSplineFilter.h ExtendSplineFilter.cxx)

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${VTK_LIBRARIES} ${ITK_LIBRARIES} optimized opencv_core300 optimized opencv_ml300 zlib)