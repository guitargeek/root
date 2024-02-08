void exampleMyMinuit2()
{
   std::string myMinuit2InstallDir = "/home/jonas/code/root/math/minuit2/install";
   gInterpreter->AddIncludePath((myMinuit2InstallDir + "/include").c_str());
   gSystem->AddDynamicPath((myMinuit2InstallDir + "/lib").c_str());
   gInterpreter->Declare("#include \"MyMinuit2/Minuit2Minimizer.h\"");

   gPluginMgr->AddHandler("ROOT::Math::Minimizer", "MyMinuit2", "ROOT::MyMinuit2::Minuit2Minimizer", "MyMinuit2",
                          "Minuit2Minimizer(const char *)");

   using namespace RooFit;

   RooRealVar x("x", "x", -10, 10);
   RooRealVar mean("mean", "mean of gaussian", 1, -10, 10);
   RooRealVar sigma("sigma", "width of gaussian", 3, 0.1, 10);
   RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
   std::unique_ptr<RooDataSet> data{gauss.generate(x, 10000)};
   std::unique_ptr<RooFitResult> result{gauss.fitTo(*data, PrintLevel(-1), Save(), Minimizer("MyMinuit2"))};
   result->Print();
}
