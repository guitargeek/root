// Tests for RooFit::setGlobalObservablesToAsimov()
// Authors: Jonas Rembser, CERN 2026

#include <RooFit/AsimovGlobalObservables.h>

#include <RooArgSet.h>
#include <RooBifurGauss.h>
#include <RooConstVar.h>
#include <RooFormulaVar.h>
#include <RooGamma.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooHelpers.h>
#include <RooLognormal.h>
#include <RooPoisson.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooRealVar.h>
#include <RooUniform.h>

#include "gtest/gtest.h"

#include <cmath>
#include <string>

namespace {

// A global observable: a constant variable with the value 0.
RooRealVar makeGlob(const char *name, double min = -100., double max = 100.)
{
   RooRealVar out{name, name, 0., min, max};
   out.setConstant(true);
   return out;
}

std::string messagesOf(RooFit::MsgLevel level, RooAbsCollection const &constraints, RooArgSet const &globs)
{
   RooHelpers::HijackMessageStream hijack{level, RooFit::Generation};
   RooFit::setGlobalObservablesToAsimov(constraints, globs);
   return hijack.str();
}

} // namespace

// The terms that relate the global observable to the parameter: the global observable takes its value, whichever
// of the two arguments it is.
TEST(AsimovGlobalObservables, Gaussian)
{
   RooRealVar alpha{"alpha", "alpha", 0.7, -5., 5.};
   RooRealVar beta{"beta", "beta", -1.3, -5., 5.};
   RooRealVar nomAlpha = makeGlob("nom_alpha");
   RooRealVar nomBeta = makeGlob("nom_beta");
   RooRealVar nomSum = makeGlob("nom_sum");
   // A constant width that is a variable, like the ones in old workspaces.
   RooRealVar sigma{"sigma", "sigma", 1.0, 0.1, 10.};
   sigma.setConstant(true);
   RooFormulaVar sum{"sum", "alpha + 2 * beta", {alpha, beta}};

   // Like HistFactory: the parameter is x, the global observable is the mean.
   RooGaussian asMean{"asMean", "asMean", alpha, nomAlpha, sigma};
   RooGaussian asX{"asX", "asX", nomBeta, beta, RooFit::RooConst(0.5)};
   // The mean is a function of several parameters.
   RooGaussian ofFunction{"ofFunction", "ofFunction", nomSum, sum, sigma};

   const RooArgSet globs{nomAlpha, nomBeta, nomSum};
   EXPECT_EQ(messagesOf(RooFit::WARNING, RooArgSet{asMean, asX, ofFunction}, globs), "");
   EXPECT_DOUBLE_EQ(nomAlpha.getVal(), 0.7);
   EXPECT_DOUBLE_EQ(nomBeta.getVal(), -1.3);
   EXPECT_DOUBLE_EQ(nomSum.getVal(), 0.7 - 2.6);
   // Each term is at its maximum now.
   RooArgSet nset{nomAlpha};
   const double atAsimov = asMean.getVal(&nset);
   alpha.setVal(0.8);
   EXPECT_GT(atAsimov, asMean.getVal(&nset));

   // Also for a parameter that is constant.
   alpha.setVal(-2.0);
   alpha.setConstant(true);
   RooFit::setGlobalObservablesToAsimov(RooArgSet{asMean}, globs);
   EXPECT_DOUBLE_EQ(nomAlpha.getVal(), -2.0);
   // Only the given constraint terms are used.
   EXPECT_DOUBLE_EQ(nomBeta.getVal(), -1.3);
}

TEST(AsimovGlobalObservables, BifurGauss)
{
   RooRealVar alpha{"alpha", "alpha", 0.4, -5., 5.};
   RooRealVar nom = makeGlob("nom");
   RooBifurGauss term{"term", "term", nom, alpha, RooFit::RooConst(1.0), RooFit::RooConst(2.0)};
   RooFit::setGlobalObservablesToAsimov(RooArgSet{term}, RooArgSet{nom});
   EXPECT_DOUBLE_EQ(nom.getVal(), 0.4);
}

// The expected number of events of the auxiliary measurement is no integer in general, so the rounding of the
// Poisson is switched off.
TEST(AsimovGlobalObservables, Poisson)
{
   RooRealVar gamma{"gamma", "gamma", 1.07, 0., 5.};
   RooRealVar nom = makeGlob("nom", 0., 1000.);
   RooProduct mean{"mean", "mean", {gamma, RooFit::RooConst(25.)}};
   RooPoisson term{"term", "term", nom, mean};
   ASSERT_FALSE(term.getNoRounding());

   RooFit::setGlobalObservablesToAsimov(RooArgSet{term}, RooArgSet{nom});
   EXPECT_DOUBLE_EQ(nom.getVal(), 1.07 * 25.);
   EXPECT_TRUE(term.getNoRounding());
}

// In the standard parametrization of the log-normal, the parameter is the logarithm of the median.
TEST(AsimovGlobalObservables, Lognormal)
{
   RooRealVar k{"k", "k", 1.3, 0.1, 10.};
   RooRealVar nom = makeGlob("nom", 0.01, 100.);
   const RooArgSet globs{nom};

   RooLognormal asX{"asX", "asX", nom, k, RooFit::RooConst(1.2)};
   RooFit::setGlobalObservablesToAsimov(RooArgSet{asX}, globs);
   EXPECT_DOUBLE_EQ(nom.getVal(), 1.3);

   RooLognormal asMedian{"asMedian", "asMedian", k, nom, RooFit::RooConst(1.2)};
   nom.setVal(1.0);
   RooFit::setGlobalObservablesToAsimov(RooArgSet{asMedian}, globs);
   EXPECT_DOUBLE_EQ(nom.getVal(), 1.3);

   RooLognormal standardAsX{"standardAsX", "standardAsX", nom, k, RooFit::RooConst(0.2), true};
   RooFit::setGlobalObservablesToAsimov(RooArgSet{standardAsX}, globs);
   EXPECT_DOUBLE_EQ(nom.getVal(), std::exp(1.3));

   RooLognormal standardAsMu{"standardAsMu", "standardAsMu", k, nom, RooFit::RooConst(0.2), true};
   RooFit::setGlobalObservablesToAsimov(RooArgSet{standardAsMu}, globs);
   EXPECT_DOUBLE_EQ(nom.getVal(), std::log(1.3));
}

// The Gamma constraint of HistFactory: Gamma(b | k = n + 1, theta), whose mode in b is n * theta.
TEST(AsimovGlobalObservables, Gamma)
{
   RooRealVar b{"b", "b", 1.2, 0., 10.};
   RooRealVar n = makeGlob("n", 0., 1000.);
   RooFormulaVar k{"k", "n + 1", {n}};
   RooRealVar theta{"theta", "theta", 0.05, 0.01, 1.};
   theta.setConstant(true);
   RooGamma term{"term", "term", b, k, theta, RooFit::RooConst(0.)};

   EXPECT_EQ(messagesOf(RooFit::WARNING, RooArgSet{term}, RooArgSet{n}), "");
   EXPECT_DOUBLE_EQ(n.getVal(), 1.2 / 0.05);
}

// Terms of other types: the global observable takes the value of the only argument that depends on parameters.
TEST(AsimovGlobalObservables, OtherTypes)
{
   RooRealVar alpha{"alpha", "alpha", 0.6, -5., 5.};
   RooRealVar nom = makeGlob("nom");
   RooGenericPdf term{"term", "exp(-0.5 * (nom - alpha) * (nom - alpha) / 4.)", {nom, alpha}};

   const std::string messages = messagesOf(RooFit::WARNING, RooArgSet{term}, RooArgSet{nom});
   EXPECT_NE(messages.find("constraint term term of type RooGenericPdf is a non-supported type"), std::string::npos)
      << messages;
   EXPECT_DOUBLE_EQ(nom.getVal(), 0.6);

   // With two such arguments, it is not known which one stands for the parameter.
   RooRealVar beta{"beta", "beta", 2.0, -5., 5.};
   RooGenericPdf twoParameters{"twoParameters", "exp(-0.5 * (nom - alpha) * (nom - alpha) / beta)", {nom, alpha, beta}};
   nom.setVal(0.);
   const std::string errors = messagesOf(RooFit::ERROR, RooArgSet{twoParameters}, RooArgSet{nom});
   EXPECT_NE(errors.find("has more than one server that depends on parameters"), std::string::npos) << errors;
   EXPECT_DOUBLE_EQ(nom.getVal(), 0.);
}

// Products that group constraint terms are set through their factors, and terms that can't be set are skipped.
TEST(AsimovGlobalObservables, ProductsAndSkippedTerms)
{
   RooRealVar alpha{"alpha", "alpha", 0.7, -5., 5.};
   RooRealVar beta{"beta", "beta", -1.3, -5., 5.};
   RooRealVar nomAlpha = makeGlob("nom_alpha");
   RooRealVar nomBeta = makeGlob("nom_beta");
   RooRealVar other = makeGlob("other", 0.1, 100.);
   RooGaussian cAlpha{"cAlpha", "cAlpha", alpha, nomAlpha, RooFit::RooConst(1.0)};
   RooGaussian cBeta{"cBeta", "cBeta", nomBeta, beta, RooFit::RooConst(1.0)};
   RooProdPdf group{"group", "group", {cAlpha, cBeta}};
   const RooArgSet globs{nomAlpha, nomBeta, other};

   EXPECT_EQ(messagesOf(RooFit::WARNING, RooArgSet{group}, globs), "");
   EXPECT_DOUBLE_EQ(nomAlpha.getVal(), 0.7);
   EXPECT_DOUBLE_EQ(nomBeta.getVal(), -1.3);

   // A uniform term has nothing to set.
   RooUniform uniform{"uniform", "uniform", alpha};
   EXPECT_EQ(messagesOf(RooFit::WARNING, RooArgSet{uniform}, globs), "");

   // A term without a global observable, e.g. a Gaussian that constrains the parameter to a number.
   RooGaussian toNumber{"toNumber", "toNumber", alpha, RooFit::RooConst(0.), RooFit::RooConst(1.0)};
   const std::string warnings = messagesOf(RooFit::WARNING, RooArgSet{toNumber}, globs);
   EXPECT_NE(warnings.find("constraint term toNumber has no global observables"), std::string::npos) << warnings;

   // A term in which the global observable is not one of the two sides.
   RooGaussian asSigma{"asSigma", "asSigma", alpha, beta, other};
   other.setVal(1.0);
   const std::string errors = messagesOf(RooFit::ERROR, RooArgSet{asSigma}, globs);
   EXPECT_NE(errors.find("constraint term asSigma has no direct dependence on global observable"), std::string::npos)
      << errors;
   EXPECT_DOUBLE_EQ(other.getVal(), 1.0);

   // A term with two global observables.
   RooGaussian twoGlobs{"twoGlobs", "twoGlobs", nomAlpha, nomBeta, RooFit::RooConst(1.0)};
   const std::string errors2 = messagesOf(RooFit::ERROR, RooArgSet{twoGlobs}, globs);
   EXPECT_NE(errors2.find("constraint term twoGlobs has multiple global observables"), std::string::npos) << errors2;
}
