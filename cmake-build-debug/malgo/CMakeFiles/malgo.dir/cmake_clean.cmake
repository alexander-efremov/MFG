file(REMOVE_RECURSE
  "libmalgo.pdb"
  "libmalgo.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/malgo.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
