# MegaMol
# Copyright (c) 2021, MegaMol Dev Team
# All rights reserved.
#

megamol_plugin(trailvolume
  BUILD_DEFAULT ON
  DEPENDS_PLUGINS
    geometry_calls
    mesh
    mmstd
  )

if (trailvolume_PLUGIN_ENABLED)
  # Do additional plugin specific stuff, not covered by the generic plugin creation,
  # i.e. find_package() dependencies, additional sources 'target_sources(megamolplugin PRIVATE ...)', ....
  find_package(voro CONFIG REQUIRED)
  find_package(Eigen3 CONFIG REQUIRED)

  target_link_libraries(trailvolume
    PRIVATE
      VORO::voro++
      Eigen3::Eigen)
endif ()
