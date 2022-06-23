# MegaMol
# Copyright (c) 2021, MegaMol Dev Team
# All rights reserved.
#

megamol_plugin(trailvolume
  BUILD_DEFAULT ON
  DEPENDS_PLUGINS
    geometry_calls
    mesh
    datatools
    # plugin2
  DEPENDS_EXTERNALS
    voro++
    Eigen
    # external2
  )

if (trailvolume_PLUGIN_ENABLED)
  # Do additional plugin specific stuff, not covered by the generic plugin creation,
  # i.e. find_package() dependencies, additional sources 'target_sources(megamolplugin PRIVATE ...)', ....
endif ()
