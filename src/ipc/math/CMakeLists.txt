set(SOURCES
  heaviside_type.hpp
  interval.cpp
  interval.hpp
  math.cpp
  math.hpp
  math.tpp
  morton.hpp
  scalar_math.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})