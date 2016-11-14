#ifndef __FPANALYSISTHISTOGRAM_H
#define __FPANALYSISTHISTOGRAM_H

#include "FPBinaryBlob.h"
#include "FPCodeGen.h"
#include "FPConfig.h"
#include "FPAnalysis.h"

#include <gsl/gsl_histogram.h>
#include <cfloat>
#include <openssl/sha.h>
#include "count/hll.h"
using libcount::HLL;

#define HISTOGRAMS_PER_INST 2
#define HISTOGRAM_SIZE 77
#define HISTOGRAM_TIMESTEP 10000000 /* write out histogram whenever sum % X */

#define POI_COUNT 43
#define GAMMA 0.57721566490153286060651209008240243
#define LAPLACE_LIMIT 0.66274341934918158097474209710925290
#define ZETA3 1.2020569
#define PHI 1.61803398874989484820458683436563811
#define SQRT3 1.73205080756887729352744634150587236
#define SQRT5 2.23606797749978969640917366873127623

extern double hist_breaks[HISTOGRAM_SIZE];

namespace FPInst {

/**
 * Stores histogram data for a single instruction.
 */
struct FPAnalysisTHistogramInstData {
    FPSemantics *inst;
    long double min, max;
    unsigned long count;
    gsl_histogram *histogram_in;
    gsl_histogram *histogram_out;
};

/**
 * Performs histogram tracking analysis.
 */
class FPAnalysisTHistogram : public FPAnalysis {

    public:

        static bool existsInstance();
        static FPAnalysisTHistogram* getInstance(bool perInst);

        string getTag();
        string getDescription();

        void configure(FPConfig *config, FPDecoder *decoder,
                FPLog *log, FPContext *context);

        void configAddresses(FPConfig *config);
        bool isRestrictedByAddress();
        string listAddresses();

        bool shouldPreInstrument(FPSemantics *inst);
        bool shouldPostInstrument(FPSemantics *inst);
        bool shouldReplace(FPSemantics *inst);

        Snippet::Ptr buildPreInstrumentation(FPSemantics *inst,
                BPatch_addressSpace *app, bool &needsRegisters);
        Snippet::Ptr buildPostInstrumentation(FPSemantics *inst,
                BPatch_addressSpace *app, bool &needsRegisters);

        string finalInstReport();

        void registerInstruction(FPSemantics *inst);
        void handlePreInstruction(FPSemantics *inst);
        void handlePostInstruction(FPSemantics *inst);
        void handleReplacement(FPSemantics *inst);

        void checkHistogram(FPSemantics *inst, FPOperand *op, bool out);

        void fprintfHistogram(FILE * stream, const gsl_histogram * h, double min, double max);
        void fprintfHistogramTimeseries(FILE * stream, const gsl_histogram * h, double t);
        void finalOutput();

    private:

        FPAnalysisTHistogram(bool histPerInst);

        void expandInstData(size_t newSize);
        FPAnalysisTHistogramInstData *instData;
        bool perInst;
        size_t instCount;

        void* histogramAddresses[256];
        size_t numHistogramAddresses;

        size_t insnsInstrumented;

        gsl_histogram *histogram_in, *histogram_out, *hist_step_in, *hist_step_out;
        long double min_in, max_in, min_out, max_out;
	std::map<double, const char*> poi_names;
	std::map<double, size_t> poi_in, poi_out;

        // HyperLogLog++ objects used to track set cardinality.
        HLL *hll_in, *hll_out, *hll_all, *hll_step_in, *hll_step_out, *hll_step_all;

        bool isPOI(long double num);
        uint64_t hash(double i);
};

}

#endif

