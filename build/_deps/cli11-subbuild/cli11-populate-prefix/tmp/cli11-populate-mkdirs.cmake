# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-src"
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-build"
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix"
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix/tmp"
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix/src/cli11-populate-stamp"
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix/src"
  "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix/src/cli11-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix/src/cli11-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/maps/projects/fernandezguerra/apps/repos/genopack/build/_deps/cli11-subbuild/cli11-populate-prefix/src/cli11-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
