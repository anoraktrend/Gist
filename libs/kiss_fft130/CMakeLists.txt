add_library(KissFFT STATIC
    kiss_fft.h
    kiss_fft.c
    _kiss_fft_guts.h
    kissfft.hh
)

target_include_directories(KissFFT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
