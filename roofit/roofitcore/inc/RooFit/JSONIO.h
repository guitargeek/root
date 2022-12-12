#ifndef RooFit_JSONIO_h
#define RooFit_JSONIO_h

#include <fstream>
#include <memory>

class RooAbsPdf;
class RooWorkspace;

namespace RooFit {

void exportJSON(std::ostream &os, RooWorkspace const &ws);

void importJSON(std::istream &is, RooWorkspace &ws);

void exportJSON(std::ostream &os, RooAbsPdf const &pdf);

std::unique_ptr<RooAbsPdf> importJSON(std::istream &is, std::string const &pdfName);

} // namespace RooFit

#endif
