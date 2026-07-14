#ifndef BETA_FITTER_MACRO_UTILS_H
#define BETA_FITTER_MACRO_UTILS_H

#include <TString.h>
#include <TSystem.h>

namespace BetaFitterMacro {

static TString result_dir(TString outputTag)
{
    while (outputTag.Length() > 1 && outputTag.EndsWith("/"))
        outputTag.Remove(outputTag.Length() - 1);
    if (outputTag.Length() == 0)
        outputTag = "plots";
    if (outputTag.BeginsWith("/") || outputTag.BeginsWith("macro/results/") ||
        outputTag.Contains("/"))
        return outputTag;
    return "macro/results/" + outputTag;
}

static TString output_path(const char *outputName,
                           const char *outputTag,
                           const TString &defaultName)
{
    if (outputName && outputName[0] != '\0')
        return outputName;
    TString directory = result_dir(outputTag ? outputTag : "");
    gSystem->mkdir(directory, true);
    return directory + "/" + defaultName;
}

} // namespace BetaFitterMacro

#endif
