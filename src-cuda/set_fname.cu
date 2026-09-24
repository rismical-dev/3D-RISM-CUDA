#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: set_fname (std::string control, std::string structure) {
  if (structure.empty()) {
    fname.append(control);
  } else {
    fname.append(structure);
  }

  // Strip only a trailing filename extension -- the last '.' after the
  // last '/' -- rather than the last '.' anywhere in the path. A directory
  // component that happens to contain a '.' (e.g.
  // ".../case.v2/input.ctrl") would otherwise get truncated at that dot
  // instead of the filename's own, discarding the rest of the path
  // (fname would end up as ".../case" instead of ".../case.v2/input").
  // If there's no '.' in the filename itself, fname is left as the full
  // path, same as the original rfind(".")-not-found behavior.
  size_t lastslash = fname.find_last_of('/');
  size_t searchfrom = (lastslash == std::string::npos) ? 0 : lastslash + 1;
  size_t dotpos = fname.rfind(".");
  if (dotpos != std::string::npos && dotpos >= searchfrom) {
    fname = fname.substr(0, dotpos);
  }
}
