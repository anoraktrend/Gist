# FFT Library Configuration - flags passed from parent cmake
# Supported flags: -DUSE_FFTW, -DUSE_KISS_FFT, -DUSE_ACCELERATE_FFT

# Add Kiss FFT subdirectory (always available as fallback)
add_subdirectory(../libs/kiss_fft130 ${CMAKE_CURRENT_BINARY_DIR}/kiss_fft130)

# Determine which FFT backend to use based on GIST_USE_FFT variable from parent
# Priority: ACCELERATE > FFTW > KISS (default)
if(GIST_USE_FFT STREQUAL "ACCELERATE")
    set(GIST_FFT_BACKEND "ACCELERATE_FFT")
elseif(GIST_USE_FFT STREQUAL "FFTW")
    set(GIST_FFT_BACKEND "FFTW")
else()
    set(GIST_FFT_BACKEND "KISS_FFT")
endif()

message(STATUS "Using FFT backend: ${GIST_FFT_BACKEND}")

add_library (
    Gist STATIC
    AccelerateFFT.cpp
    AccelerateFFT.h
    CoreFrequencyDomainFeatures.cpp
    CoreFrequencyDomainFeatures.h
    CoreTimeDomainFeatures.cpp
    CoreTimeDomainFeatures.h
    Gist.cpp
    Gist.h
    MFCC.cpp
    MFCC.h
    OnsetDetectionFunction.cpp
    OnsetDetectionFunction.h
    WindowFunctions.cpp
    WindowFunctions.h
    Yin.cpp
    Yin.h
)

target_include_directories(Gist PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../libs/kiss_fft130
)

source_group(Source src)

# Configure Gist based on selected FFT backend
if(GIST_FFT_BACKEND STREQUAL "ACCELERATE_FFT")
    target_link_libraries(Gist PUBLIC "-framework Accelerate")
    target_compile_definitions(Gist PUBLIC -DUSE_ACCELERATE_FFT)
    message(STATUS "Gist configured for Accelerate FFT")
elseif(GIST_FFT_BACKEND STREQUAL "FFTW")
    target_link_libraries(Gist PUBLIC FFTW3::FFTW)
    target_compile_definitions(Gist PUBLIC -DUSE_FFTW)
    message(STATUS "Gist configured for FFTW")
else()
    target_link_libraries(Gist PUBLIC KissFFT)
    target_compile_definitions(Gist PUBLIC -DUSE_KISS_FFT)
    message(STATUS "Gist configured for Kiss FFT")
endif()
