set(SOURCES
  interval.cpp
  interval.hpp
  math.cpp
  math.hpp
  math.tpp
  morton.hpp
  sinc.cpp
  sinc.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
