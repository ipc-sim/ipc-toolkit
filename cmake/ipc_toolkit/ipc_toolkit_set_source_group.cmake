# Set source group for IDE like Visual Studio or XCode
function(ipc_toolkit_set_source_group)
  foreach(filepath IN ITEMS ${ARGN})
    get_filename_component(folderpath "${filepath}" DIRECTORY)
    get_filename_component(foldername "${folderpath}" NAME)
    source_group(foldername FILES "${filepath}")
  endforeach()
endfunction()
