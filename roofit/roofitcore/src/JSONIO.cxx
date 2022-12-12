#include <RooFit/JSONIO.h>

#include <RooAbsPdf.h>
#include <RooArgProxy.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooFit/Detail/JSONInterface.h>

namespace RooFit {

namespace {

class RooFitJSONTree {
public:
   RooFitJSONTree() : _tree{RooFit::Detail::JSONTree::create()} { _tree->rootnode().set_map(); }

   RooFitJSONTree(std::istream &is) : _tree{RooFit::Detail::JSONTree::create(is)} {}

   auto &rootnode() { return _tree->rootnode(); }

   void writeJSON(std::ostream &os) { _tree->rootnode().writeJSON(os); }

private:
   std::unique_ptr<RooFit::Detail::JSONTree> _tree;
};

} // namespace

void exportJSON(std::ostream &os, RooWorkspace const &ws)
{
   RooFitJSONTree tree;
   auto &n = tree.rootnode();

   auto &distributions = n["distributions"];
   distributions.set_seq();

   auto &defaultDomainDict = n["domains"].append_child();
   defaultDomainDict.set_map();
   defaultDomainDict["name"] << "default";
   auto &defaultDomain = defaultDomainDict["values"];
   defaultDomain.set_seq();

   auto &defaultValuesDict = n["estimates"].append_child();
   defaultValuesDict.set_map();
   defaultValuesDict["name"] << "default";
   auto &defaultValues = defaultValuesDict["values"];
   defaultValues.set_seq();

   // actual import over workspace

   RooArgSet allVars{ws.allVars()};
   allVars.sort();
   for (auto *var : static_range_cast<RooRealVar *>(allVars)) {

      {
         auto &node = defaultValues.append_child();
         node.set_map();

         node["name"] << var->GetName();
         node["value"] << var->getVal();
      }

      {
         auto &node = defaultDomain.append_child();
         node.set_map();

         node["name"] << var->GetName();
         node["min"] << var->getMin();
         node["max"] << var->getMax();
      }
   }

   RooArgSet allPdfs{ws.allPdfs()};
   allPdfs.sort();
   for (auto *pdf : static_range_cast<RooAbsPdf *>(allPdfs)) {
      auto &node = distributions.append_child();
      node.set_map();
      node["name"] << pdf->GetName();
      node["title"] << pdf->GetTitle();
      node["type"] << pdf->ClassName();

      for (int iProxy = 0; iProxy < pdf->numProxies(); ++iProxy) {
         if (auto const *proxy = dynamic_cast<RooArgProxy const *>(pdf->getProxy(iProxy))) {
            node[proxy->name()] << proxy->absArg()->GetName();
         }
      }
   }

   tree.writeJSON(os);
}

void importJSON(std::istream &is, RooWorkspace &ws)
{

   std::vector<std::string> expressions;

   RooFitJSONTree tree(is);

   auto &n = tree.rootnode();

   RooFit::Detail::JSONNode const *defaultDomainVals = nullptr;
   RooFit::Detail::JSONNode const *defaultValuesVals = nullptr;

   for (auto const &domain : n["domains"].children()) {
      for (auto const &item : domain.children()) {
         if (item.key() == "name" && item.val() == "default") {
            defaultDomainVals = &domain["values"];
         }
      }
   }
   for (auto const &estimate : n["estimates"].children()) {
      for (auto const &item : estimate.children()) {
         if (item.key() == "name" && item.val() == "default") {
            defaultValuesVals = &estimate["values"];
         }
      }
   }

   assert(defaultDomainVals);
   assert(defaultValuesVals);

   struct VarInfo {
      double val = 0.0;
      double min = 0.0;
      double max = 0.0;
   };

   std::unordered_map<std::string, VarInfo> varInfos;

   for (auto const &node : defaultDomainVals->children()) {
      auto &info = varInfos[node["name"].val()];
      info.min = node["min"].val_float();
      info.max = node["max"].val_float();
   }

   for (auto const &node : defaultValuesVals->children()) {
      auto &info = varInfos[node["name"].val()];
      info.val = node["value"].val_float();
   }

   for (auto const &iter : varInfos) {
      std::string const &name = iter.first;
      VarInfo const &info = iter.second;
      std::stringstream expr;
      expr << name << "[" << info.val << "," << info.min << "," << info.max << "]";
      expressions.emplace_back(expr.str());
   }

   for (auto const &node : n["distributions"].children()) {
      if (node["type"].val() == "RooGaussian") {
         std::stringstream expr;
         expr << node["type"].val() << "::" << node["name"].val() << "(" << node["x"].val() << "," << node["mean"].val()
              << "," << node["sigma"].val() << ")";
         expressions.emplace_back(expr.str());
      } else {
         throw std::runtime_error("unknown PDF encountered!");
      }
   }

   for (std::string const &expr : expressions) {
      ws.factory(expr.c_str());
   }

   for (auto const &node : n["distributions"].children()) {
      ws.pdf(node["name"].val())->SetTitle(node["title"].val().c_str());
   }
}

void exportJSON(std::ostream &os, RooAbsPdf const &pdf)
{
   RooWorkspace ws;
   ws.import(pdf);
   exportJSON(os, ws);
}

std::unique_ptr<RooAbsPdf> importJSON(std::istream &is, std::string const &pdfName)
{
   RooWorkspace ws;
   importJSON(is, ws);
   return std::unique_ptr<RooAbsPdf>{static_cast<RooAbsPdf *>(ws.pdf(pdfName)->cloneTree())};
}

} // namespace RooFit
