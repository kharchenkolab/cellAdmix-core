include(CMakeParseArguments)

function(celladmix_find_dependency package_name)
  set(options PREFER_PKG_CONFIG)
  set(one_value_args OUT_TARGET PKG_CONFIG_NAME INSTALL_HINT)
  set(multi_value_args TARGETS)
  cmake_parse_arguments(CELLADMIX_DEP "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  foreach(_target IN LISTS CELLADMIX_DEP_TARGETS)
    if(TARGET "${_target}")
      set(_selected_target "${_target}")
      break()
    endif()
  endforeach()

  if(NOT _selected_target AND CELLADMIX_DEP_PREFER_PKG_CONFIG AND CELLADMIX_DEP_PKG_CONFIG_NAME)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
      string(MAKE_C_IDENTIFIER "${package_name}" _pkg_id)
      set(_pc_prefix "CELLADMIX_PC_${_pkg_id}")
      pkg_check_modules(${_pc_prefix} QUIET IMPORTED_TARGET ${CELLADMIX_DEP_PKG_CONFIG_NAME})
      if(${_pc_prefix}_FOUND AND TARGET "PkgConfig::${_pc_prefix}")
        set(_selected_target "PkgConfig::${_pc_prefix}")
        get_target_property(_pc_options "${_selected_target}" INTERFACE_COMPILE_OPTIONS)
        if(_pc_options)
          list(FILTER _pc_options EXCLUDE REGEX "^-std=")
          set_target_properties("${_selected_target}" PROPERTIES INTERFACE_COMPILE_OPTIONS "${_pc_options}")
        endif()
      endif()
    endif()
  endif()

  if(NOT _selected_target)
    find_package(${package_name} QUIET)
    foreach(_target IN LISTS CELLADMIX_DEP_TARGETS)
      if(TARGET "${_target}")
        set(_selected_target "${_target}")
        break()
      endif()
    endforeach()
  endif()

  if(package_name STREQUAL "Eigen3" AND NOT TARGET Eigen3::Eigen)
    set(_eigen_hints)
    foreach(_prefix IN LISTS CMAKE_PREFIX_PATH)
      list(APPEND _eigen_hints
        "${_prefix}"
        "${_prefix}/include"
        "${_prefix}/include/eigen3")
    endforeach()
    if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
      list(APPEND _eigen_hints
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}"
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include"
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include/eigen3")
    endif()
    list(APPEND _eigen_hints
      /opt/homebrew
      /opt/homebrew/include
      /opt/homebrew/include/eigen3
      /usr/local
      /usr/local/include
      /usr/local/include/eigen3
      /usr
      /usr/include
      /usr/include/eigen3)

    find_path(EIGEN3_INCLUDE_DIR
      NAMES Eigen/Core
      HINTS ${_eigen_hints}
      PATH_SUFFIXES eigen3 include/eigen3)
    if(EIGEN3_INCLUDE_DIR)
      add_library(Eigen3::Eigen INTERFACE IMPORTED)
      target_include_directories(Eigen3::Eigen INTERFACE "${EIGEN3_INCLUDE_DIR}")
      set(Eigen3_FOUND TRUE)
      message(STATUS "Found Eigen3 headers: ${EIGEN3_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT _selected_target AND CELLADMIX_DEP_PKG_CONFIG_NAME)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
      string(MAKE_C_IDENTIFIER "${package_name}" _pkg_id)
      set(_pc_prefix "CELLADMIX_PC_${_pkg_id}")
      pkg_check_modules(${_pc_prefix} QUIET IMPORTED_TARGET ${CELLADMIX_DEP_PKG_CONFIG_NAME})
      if(${_pc_prefix}_FOUND AND TARGET "PkgConfig::${_pc_prefix}")
        set(_selected_target "PkgConfig::${_pc_prefix}")
        get_target_property(_pc_options "${_selected_target}" INTERFACE_COMPILE_OPTIONS)
        if(_pc_options)
          list(FILTER _pc_options EXCLUDE REGEX "^-std=")
          set_target_properties("${_selected_target}" PROPERTIES INTERFACE_COMPILE_OPTIONS "${_pc_options}")
        endif()
      endif()
    endif()
  endif()

  if(NOT _selected_target)
    message(FATAL_ERROR
      "Required dependency '${package_name}' was not found.\n"
      "${CELLADMIX_DEP_INSTALL_HINT}\n"
      "If the dependency is installed in a custom prefix, reconfigure with:\n"
      "  cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/prefix\n"
      "or set ${package_name}_DIR to the directory containing ${package_name}Config.cmake."
    )
  endif()

  if(CELLADMIX_DEP_OUT_TARGET)
    set(${CELLADMIX_DEP_OUT_TARGET} "${_selected_target}" PARENT_SCOPE)
  endif()
  message(STATUS "Using ${package_name} target: ${_selected_target}")
endfunction()
