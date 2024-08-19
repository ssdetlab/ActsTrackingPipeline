file(REMOVE_RECURSE
  "../lib/libLUXESimEventDict.dylib"
  "../lib/libLUXESimEventDict.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/LUXESimEventDict.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
