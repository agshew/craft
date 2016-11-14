#include "FPAnalysisTHistogram.h"

double hist_breaks[HISTOGRAM_SIZE] =
{-DBL_MAX,
 -1000000000000000.0,
 -100000000000000.0,
 -10000000000000.0,
 -1000000000000.0,
 -100000000000.0,
 -10000000000.0,
 -1000000000.0,
 -100000000.0,
 -10000000.0,
 -1000000.0,
 -100000.0,
 -10000.0,
 -1000.0,
 -100.0,
 -10.0,
 -8.0,
 -4.0,
 -2.0,
 -1.0,
 -0.5,
 -0.25,
 -0.125,
 -0.0625,
 -0.03125,
 -0.015625,
 -0.0078125,
 -0.00390625,
 -0.001953125,
 -0.0009765625,
 -1.9073486328125e-06,
 -9.5367431640625e-07,
 -DBL_MIN,
 -2e-312,
 -2e-315,
 -2e-318,
 -2e-321,
 -5e-324,
 0.0,
 5e-324,		/* smallest double denormal */
 2e-321,
 2e-318,
 2e-315,
 2e-312,
 DBL_MIN,		/* smallest double normal, 2.22507e-308*/
 9.5367431640625e-07,	/* 1.0 / pow(2,19) */
 1.9073486328125e-06,	/* 1.0 / pow(2,20) */
 0.0009765625,		/* 1.0 / pow(2,10) */
 0.001953125,
 0.00390625,
 0.0078125,
 0.015625,
 0.03125,
 0.0625,
 0.125,
 0.25,
 0.5,
 1.0,
 2.0,
 4.0,
 8.0,
 10.0,
 100.0,
 1000.0,
 10000.0,
 100000.0,
 1000000.0,
 10000000.0,
 100000000.0,
 1000000000.0,
 10000000000.0,
 100000000000.0,
 1000000000000.0,
 10000000000000.0,
 100000000000000.0,
 1000000000000000.0,
 DBL_MAX};

/* Add float and half versions (i.e. *_F and *_H)? */
double poi_values[POI_COUNT] =
{-INFINITY,
 -0.0,
 0.0,
 0.1,
 0.2,
 0.3,
 M_1_PI,
 1.0 / 3.0,
 0.4,
 M_LOG10E,
 GAMMA,
 0.6,
 M_2_PI,
 LAPLACE_LIMIT,
 2.0 / 3.0,
 M_LN2,
 0.7,
 M_SQRT1_2,
 M_PI_4,
 0.8,
 0.9,
 1.0,
 1.1,
 M_2_SQRTPI,
 1.2,
 ZETA3,
 1.3,
 1.4,
 M_SQRT2,
 M_PI_2,
 1.6,
 PHI,
 1.7,
 SQRT3,
 1.8,
 1.9,
 SQRT5,
 M_LN10,
 M_E,
 M_PI,
 2 * M_PI,
 HUGE_VAL,
 INFINITY};

namespace FPInst {

extern "C" {
    FPAnalysisTHistogram *_INST_Main_THistogramAnalysis = NULL;
}

bool FPAnalysisTHistogram::existsInstance()
{
    return (_INST_Main_THistogramAnalysis != NULL);
}

FPAnalysisTHistogram* FPAnalysisTHistogram::getInstance(bool perInst)
{
    if (!_INST_Main_THistogramAnalysis) {
        _INST_Main_THistogramAnalysis = new FPAnalysisTHistogram(perInst);
    }
    return _INST_Main_THistogramAnalysis;
}

FPAnalysisTHistogram::FPAnalysisTHistogram(bool histPerInst)
    : FPAnalysis()
{
    numHistogramAddresses = 0;
    instCount = 0;
    insnsInstrumented = 0;
    instData = NULL;

    perInst = histPerInst;
    if (perInst)
        expandInstData(8192);

    histogram_in = gsl_histogram_alloc(HISTOGRAM_SIZE - 1);
    gsl_histogram_set_ranges(histogram_in, hist_breaks, HISTOGRAM_SIZE);

    histogram_out = gsl_histogram_alloc(HISTOGRAM_SIZE - 1);
    gsl_histogram_set_ranges(histogram_out, hist_breaks, HISTOGRAM_SIZE);

    hist_step_in = gsl_histogram_alloc(HISTOGRAM_SIZE - 1);
    gsl_histogram_set_ranges(hist_step_in, hist_breaks, HISTOGRAM_SIZE);

    hist_step_out = gsl_histogram_alloc(HISTOGRAM_SIZE - 1);
    gsl_histogram_set_ranges(hist_step_out, hist_breaks, HISTOGRAM_SIZE);

    poi_names[-INFINITY] = "-infinity";
    poi_names[-0.0] = "-0";
    poi_names[0.0] = "0";
    poi_names[0.1] = "0.1";
    poi_names[0.2] = "0.2";
    poi_names[0.3] = "0.3";
    poi_names[M_1_PI] = "1/pi";
    poi_names[1.0 / 3.0] = "1/3";
    poi_names[0.4] = "0.4";
    poi_names[M_LOG10E] = "log10(E)";
    poi_names[GAMMA] = "gamma";
    poi_names[0.6] = "0.6";
    poi_names[M_2_PI] = "2/pi";
    poi_names[LAPLACE_LIMIT] = "Laplace limit";
    poi_names[2.0 / 3.0] = "2/3";
    poi_names[M_LN2] = "ln(2)";
    poi_names[0.7] = "0.7";
    poi_names[M_SQRT1_2] = "1/sqrt(2)";
    poi_names[M_PI_4] = "pi/4";
    poi_names[0.8] = "0.8";
    poi_names[0.9] = "0.9";
    poi_names[1.0] = "1";
    poi_names[1.1] = "1.1";
    poi_names[M_2_SQRTPI] = "2/sqrt(pi)";
    poi_names[1.2] = "1.2";
    poi_names[ZETA3] = "zeta(3)";
    poi_names[1.3] = "1.3";
    poi_names[1.4] = "1.4";
    poi_names[M_SQRT2] = "sqrt(2)";
    poi_names[M_PI_2] = "pi/2";
    poi_names[1.6] = "1.6";
    poi_names[PHI] = "phi";
    poi_names[1.7] = "1.7";
    poi_names[SQRT3] = "sqrt(3)";
    poi_names[1.8] = "1.8";
    poi_names[1.9] = "1.9";
    poi_names[SQRT5] = "sqrt(5)";
    poi_names[M_LN10] = "ln(10)";
    poi_names[M_E] = "e";
    poi_names[M_PI] = "pi";
    poi_names[2 * M_PI] = "2*pi";
    poi_names[HUGE_VAL] = "HUGE_VAL";
    poi_names[INFINITY] = "infinity";

    hll_in = HLL::Create(12);
    hll_out = HLL::Create(12);
    hll_all = HLL::Create(12);
    hll_step_in = HLL::Create(12);
    hll_step_out = HLL::Create(12);
    hll_step_all = HLL::Create(12);
}

string FPAnalysisTHistogram::getTag()
{
    return "t_histogram";
}

string FPAnalysisTHistogram::getDescription()
{
    return "Histogram Tracking Analysis";
}

void FPAnalysisTHistogram::configure(FPConfig *config, FPDecoder *decoder,
        FPLog *log, FPContext *context)
{
    FPAnalysis::configure(config, decoder, log, context);
    configAddresses(config);
    if (isRestrictedByAddress()) {
        //status << "t_histogram: addresses=" << listAddresses();
    }
}

void FPAnalysisTHistogram::configAddresses(FPConfig *config)
{
    if (config->hasValue("thistogram_addresses")) {
        numHistogramAddresses = config->getAddressList(histogramAddresses, "thistogram_addresses");
        //printf("filtering by %d addresses\n", numHistogramAddresses);
    } /*else {
        printf("no filtering\n");
    }*/
}

bool FPAnalysisTHistogram::isRestrictedByAddress()
{
    return numHistogramAddresses > 0;
}

string FPAnalysisTHistogram::listAddresses() {
    stringstream ss;
    ss.clear();
    ss.str("");
    size_t i;
    ss << hex;
    for (i=0; i<numHistogramAddresses; i++) {
        ss << histogramAddresses[i] << " ";
    }
    return ss.str();
}

bool FPAnalysisTHistogram::shouldPreInstrument(FPSemantics * inst)
{
    bool isOperation = false, handle = false;
    FPOperation *op;
    FPOperationType opt;
    size_t i;

    // is this an operation?
    for (i=0; i<inst->numOps; i++) {
        op = (*inst)[i];
        opt = op->type;
        if ((opt != OP_INVALID && opt != OP_NONE &&
             opt != OP_MOV && opt != OP_CVT &&
             opt != OP_COM && opt != OP_COMI &&
             opt != OP_UCOM && opt != OP_UCOMI) &&
            (op->hasOperandOfType(IEEE_Single) ||
             op->hasOperandOfType(IEEE_Double) ||
             op->hasOperandOfType(C99_LongDouble))) {
            isOperation = true;
            break;
        }
    }

    // does it involve any floating-point numbers?

    // if so, is this on the list of addresses to handle?
    if (isOperation && numHistogramAddresses > 0) {
        void* addr = inst->getAddress();
        size_t i;
        bool found = false;
        for (i=0; i<numHistogramAddresses; i++) {
            if (histogramAddresses[i] == addr) {
                found = true;
                break;
            }
        }
        if (found) {
            handle = true;
        } else {
            handle = false;
        }
    } else {
        handle = isOperation;
    }

    return handle;
}

bool FPAnalysisTHistogram::shouldPostInstrument(FPSemantics * inst)
{
    return shouldPreInstrument(inst);;
}

bool FPAnalysisTHistogram::shouldReplace(FPSemantics * /*inst*/)
{
    return false;
}

Snippet::Ptr FPAnalysisTHistogram::buildPreInstrumentation(FPSemantics * /*inst*/,
        BPatch_addressSpace * /*app*/, bool & needsRegisters)
{
    insnsInstrumented++;
    needsRegisters = true;
    return Snippet::Ptr();
}

Snippet::Ptr FPAnalysisTHistogram::buildPostInstrumentation(FPSemantics * /*inst*/,
        BPatch_addressSpace * /*app*/, bool & needsRegisters)
{
    needsRegisters = true;
    return Snippet::Ptr();
}

string FPAnalysisTHistogram::finalInstReport()
{
    stringstream ss;
    ss << "THistogram: " << insnsInstrumented << " instrumented";
    return ss.str();
}

void FPAnalysisTHistogram::registerInstruction(FPSemantics *inst)
{
    if (!perInst)
        return;

    size_t idx = inst->getIndex();
    if (idx >= instCount) {
        expandInstData(idx+1);
    }
    instData[idx].inst = inst;
    instData[idx].min = INFINITY;
    instData[idx].max = -INFINITY;
    instData[idx].count = 0;
}

void FPAnalysisTHistogram::handlePreInstruction(FPSemantics *inst)
{
    FPOperation *op;
    FPOperandSet *ops;
    FPOperand *input;
    size_t i, j, k;

    // for each operation ...
    for (i=0; i<inst->numOps; i++) {
        op = (*inst)[i];
        ops = op->opSets;

        // for each input ...
        for (j=0; j<op->numOpSets; j++) {
            for (k=0; k<ops[j].nIn; k++) {
                input = ops[j].in[k];
                input->refresh(context);
                checkHistogram(inst, input, false);
            }
        }
    }
}

void FPAnalysisTHistogram::handlePostInstruction(FPSemantics *inst)
{
    FPOperation *op;
    FPOperandSet *ops;
    FPOperand *output;
    size_t i, j, k;

    // for each operation ...
    for (i=0; i<inst->numOps; i++) {
        op = (*inst)[i];
        ops = op->opSets;

        // for each output ...
        for (j=0; j<op->numOpSets; j++) {
            for (k=0; k<ops[j].nOut; k++) {
                output = ops[j].out[k];
                output->refresh(context);
                checkHistogram(inst, output, true);
            }
        }
    }
}

void FPAnalysisTHistogram::handleReplacement(FPSemantics * /*inst*/)
{ }


bool FPAnalysisTHistogram::isPOI(long double num)
{
    for (int i = 0; i < POI_COUNT; i++)
        if (num == poi_values[i])
            return true;

    if (isinf(num))
            return true;

    return false;
}

uint64_t FPAnalysisTHistogram::hash(double i) {
    // Structure that is 160 bits wide used to extract 64 bits from a SHA-1.
    struct hashval {
      uint64_t high64;
      char low96[12];
    } hash;
 
    // Calculate the SHA-1 hash of the integer.
    SHA_CTX ctx;
    SHA1_Init(&ctx);
    SHA1_Update(&ctx, (unsigned char*)&i, sizeof(i));
    SHA1_Final((unsigned char*)&hash, &ctx);
 
    // Return 64 bits of the hash.
    return hash.high64;
}

void FPAnalysisTHistogram::checkHistogram(FPSemantics *inst, FPOperand *op, bool out)
{
    FILE * outFile;
    char fileName[128];
    long double num;
    long unsigned int sum;
    double timestep;
    pid_t pid;
    uint64_t numhash, card, card_step;
    double bits, bits_step, density_step, density;

    // extract operand information
    num = op->getCurrentValueLD();

    numhash = hash(num);
    hll_all->Update(numhash);
    hll_step_all->Update(numhash);

    if (!out) {
        hll_in->Update(numhash);
        hll_step_in->Update(numhash);

        gsl_histogram_increment(histogram_in, (double)num);
        gsl_histogram_increment(hist_step_in, (double)num);

        if (num  < min_in) { min_in = num;  }
        if (num  > max_in) { max_in = num;  }
        if (isPOI(num)) { poi_in[num]++;  }

        sum = (long unsigned int) gsl_histogram_sum(histogram_in);
        timestep = (double) sum / HISTOGRAM_TIMESTEP;
        if (sum > 0 && sum % HISTOGRAM_TIMESTEP == 0) {
            pid = getpid();
            cout << "writing histogram " << pid << " (input), step " << timestep << endl;
            sprintf(fileName, "hist-%d-all-in-timeseries", pid);
            outFile = fopen(fileName, "a");
            fprintfHistogramTimeseries(outFile, hist_step_in, timestep);
            fclose(outFile);
            gsl_histogram_set_ranges(hist_step_in, hist_breaks, HISTOGRAM_SIZE);

            sprintf(fileName, "hist-%d-poi-in", getpid());
            printf("writing %s\n", fileName);
            outFile = fopen(fileName, "w");
            for(map<double, size_t>::const_iterator it = poi_in.begin(); it != poi_in.end(); ++it)
                fprintf(outFile, "%g \"%s\" %lu\n", it->first, poi_names[it->first], it->second);
            fclose(outFile);

            sprintf(fileName, "hll-%d-in-timeseries", getpid());
            printf("writing %s\n", fileName);
            outFile = fopen(fileName, "a");
            card_step = hll_step_in->Estimate();
            bits_step = log2(card_step);
            density_step = bits_step > 32 ? bits_step / 64 : bits_step / 32;
            card = hll_in->Estimate();
            bits = log2(card);
            density = bits > 32 ? bits / 64 : bits / 32;
            fprintf(outFile, "%g %lu %g %g %lu %g %g\n", timestep, card_step, bits_step, density_step, card, bits, density);
            fclose(outFile);
            delete hll_step_in;
            hll_step_in = HLL::Create(12);
        }
    } else {
        hll_out->Update(numhash);
        hll_step_out->Update(numhash);

        gsl_histogram_increment(histogram_out, (double)num);
        gsl_histogram_increment(hist_step_out, (double)num);

        if (num  < min_out) { min_out = num;  }
        if (num  > max_out) { max_out = num;  }
        if (isPOI(num)) { poi_out[num]++;  }

        sum = (long unsigned int) gsl_histogram_sum(histogram_out);
        timestep = (double) sum / HISTOGRAM_TIMESTEP;
        if (sum > 0 && sum % HISTOGRAM_TIMESTEP == 0) {
            pid = getpid();
            cout << "writing histogram " << pid << " (output), step " << timestep << endl;
            sprintf(fileName, "hist-%d-all-out-timeseries", pid);
            outFile = fopen(fileName, "a");
            fprintfHistogramTimeseries(outFile, hist_step_out, timestep);
            fclose(outFile);
            gsl_histogram_set_ranges(hist_step_out, hist_breaks, HISTOGRAM_SIZE);

            sprintf(fileName, "hist-%d-poi-out", getpid());
            printf("writing %s\n", fileName);
            outFile = fopen(fileName, "w");
            for(map<double, size_t>::const_iterator it = poi_out.begin(); it != poi_out.end(); ++it)
                fprintf(outFile, "%g \"%s\" %lu\n", it->first, poi_names[it->first], it->second);
            fclose(outFile);

            sprintf(fileName, "hll-%d-out-timeseries", getpid());
            printf("writing %s\n", fileName);
            outFile = fopen(fileName, "a");
            card_step = hll_step_out->Estimate();
            bits_step = log2(card_step);
            density_step = bits_step > 32 ? bits_step / 64 : bits_step / 32;
            card = hll_out->Estimate();
            bits = log2(card);
            density = bits > 32 ? bits / 64 : bits / 32;
            fprintf(outFile, "%g %lu %g %g %lu %g %g\n", timestep, card_step, bits_step, density_step, card, bits, density);
            fclose(outFile);
            delete hll_step_out;
            hll_step_out = HLL::Create(12);

            sprintf(fileName, "hll-%d-all-timeseries", getpid());
            printf("writing %s\n", fileName);
            outFile = fopen(fileName, "a");
            card_step = hll_step_all->Estimate();
            bits_step = log2(card_step);
            density_step = bits_step > 32 ? bits_step / 64 : bits_step / 32;
            card = hll_all->Estimate();
            bits = log2(card);
            density = bits > 32 ? bits / 64 : bits / 32;
            fprintf(outFile, "%g %lu %g %g %lu %g %g\n", timestep, card_step, bits_step, density_step, card, bits, density);
            fclose(outFile);
            delete hll_step_all;
            hll_step_all = HLL::Create(12);
        }
    }

    if (perInst) {
        size_t idx = inst->getIndex();
        if (num  < instData[idx].min) { instData[idx].min = num;  }
        if (num  > instData[idx].max) { instData[idx].max = num;  }

        if (!out)
            gsl_histogram_increment(instData[idx].histogram_in, num);
        else
            gsl_histogram_increment(instData[idx].histogram_out, num);
    }
}

void FPAnalysisTHistogram::expandInstData(size_t newSize)
{
    FPAnalysisTHistogramInstData *newInstData;
    size_t i = 0;
    newSize = newSize > instCount*2 ? newSize : instCount*2;
    //printf("expand_inst_data - old size: %u    new size: %u\n", (unsigned)instCount, (unsigned)newSize);
    newInstData = (FPAnalysisTHistogramInstData*)malloc(newSize * (sizeof(FPAnalysisTHistogramInstData) + (HISTOGRAMS_PER_INST)*2*(HISTOGRAM_SIZE + 1)*sizeof(double)));
    if (!newInstData) {
        fprintf(stderr, "OUT OF MEMORY!\n");
        exit(-1);
    }
    if (instData != NULL) {
        for (; i < instCount; i++) {
            newInstData[i].inst = instData[i].inst;
            newInstData[i].min = instData[i].min;
            newInstData[i].max = instData[i].max;
            newInstData[i].count = instData[i].count;
            newInstData[i].histogram_in = instData[i].histogram_in;
            newInstData[i].histogram_out = instData[i].histogram_out;

            gsl_histogram_free(instData[i].histogram_in);
            gsl_histogram_free(instData[i].histogram_out);
        }
        free(instData);
        instData = NULL;
    }
    for (; i < newSize; i++) {
        newInstData[i].inst = NULL;
        newInstData[i].min = INFINITY;
        newInstData[i].max = -INFINITY;
        newInstData[i].count = 0;

        newInstData[i].histogram_in = gsl_histogram_alloc(HISTOGRAM_SIZE - 1);
        gsl_histogram_set_ranges(newInstData[i].histogram_in, hist_breaks, HISTOGRAM_SIZE);

        newInstData[i].histogram_out = gsl_histogram_alloc(HISTOGRAM_SIZE - 1);
        gsl_histogram_set_ranges(newInstData[i].histogram_out, hist_breaks, HISTOGRAM_SIZE);
    }
    instData = newInstData;
    instCount = newSize;
}

void FPAnalysisTHistogram::fprintfHistogramTimeseries(FILE * stream, const gsl_histogram * h, double t)
{
    int status = 0;

    for (size_t i = 0; i < h->n; i++)
        status += fprintf(stream, "%g '[' %g %g ')' %g\n", t, h->range[i], h->range[i + 1], h->bin[i]);

    if (status < 0)
        printf("fprintf failed\n");
}

void FPAnalysisTHistogram::fprintfHistogram(FILE * stream, const gsl_histogram * h, double min, double max)
{
    int status = 0;
    size_t minbin, maxbin;

    if (!isfinite(min) || min < h->range[0])
        minbin = 0;

    if (!isfinite(max) || max > h->range[h->n])
        maxbin = h->n;

    if (gsl_histogram_sum(h) > 0) {
        gsl_histogram_find(h, min, &minbin);
        gsl_histogram_find(h, max, &maxbin);
    }

    status = fprintf(stream, "'[' %g %g ')' %g\n", min, h->range[minbin + 1], h->bin[minbin]);

    for (size_t i = minbin+1; i < maxbin-1; i++)
        status += fprintf(stream, "'[' %g %g ')' %g\n", h->range[i], h->range[i + 1], h->bin[i]);

    status += fprintf(stream, "'[' %g %g ']' %g\n", h->range[maxbin], max, h->bin[maxbin]);

    if (status < 0)
        printf("fprintf failed\n");
}

void FPAnalysisTHistogram::finalOutput()
{
    FILE * outFile;
    char fileName[128];
    size_t i;
    long unsigned int sum;
    double timestep;
    uint64_t card, card_step;
    double bits, bits_step, density_step, density;

    if (perInst) {
        for (i = 0; i < instCount; i++) {
            if (instData[i].inst) {
                sprintf(fileName, "hist-%d-%lu-in", getpid(), i);
                outFile = fopen(fileName, "w");
                fprintfHistogram(outFile, instData[i].histogram_in, instData[i].min, instData[i].max);
                fclose(outFile);
    
                sprintf(fileName, "hist-%d-%lu-out", getpid(), i);
                outFile = fopen(fileName, "w");
                fprintfHistogram(outFile, instData[i].histogram_out, instData[i].min, instData[i].max);
                fclose(outFile);
            }
        }
    }

    sprintf(fileName, "hist-%d-all-in", getpid());
    //printf("writing %s\n", fileName);
    outFile = fopen(fileName, "w");
    gsl_histogram_fprintf(outFile, histogram_in, "%g", "%g");
    fclose(outFile);

    sprintf(fileName, "hist-%d-all-out", getpid());
    //printf("writing %s\n", fileName);
    outFile = fopen(fileName, "w");
    gsl_histogram_fprintf(outFile, histogram_out, "%g", "%g");
    fclose(outFile);

    sprintf(fileName, "hist-%d-poi-in", getpid());
    //printf("writing %s\n", fileName);
    outFile = fopen(fileName, "w");
    for(map<double, size_t>::const_iterator it = poi_in.begin(); it != poi_in.end(); ++it)
        fprintf(outFile, "%g \"%s\" %lu\n", it->first, poi_names[it->first], it->second);
    fclose(outFile);

    sprintf(fileName, "hist-%d-poi-out", getpid());
    //printf("writing %s\n", fileName);
    outFile = fopen(fileName, "w");
    for(map<double, size_t>::const_iterator it = poi_out.begin(); it != poi_out.end(); ++it)
        fprintf(outFile, "%g \"%s\" %lu\n", it->first, poi_names[it->first], it->second);
    fclose(outFile);

    sum = (long unsigned int) gsl_histogram_sum(histogram_in);
    timestep = (double)sum / HISTOGRAM_TIMESTEP;
    sprintf(fileName, "hist-%d-all-in-timeseries", getpid());
    outFile = fopen(fileName, "a");
    fprintfHistogramTimeseries(outFile, hist_step_in, timestep);
    fclose(outFile);

    sprintf(fileName, "hll-%d-in-timeseries", getpid());
    outFile = fopen(fileName, "a");
    card_step = hll_step_in->Estimate();
    bits_step = log2(card_step);
    density_step = bits_step > 32 ? bits_step / 64 : bits_step / 32;
    card = hll_in->Estimate();
    bits = log2(card);
    density = bits > 32 ? bits / 64 : bits / 32;
    fprintf(outFile, "%g %lu %g %g %lu %g %g\n", timestep, card_step, bits_step, density_step, card, bits, density);
    fclose(outFile);

    sum = (long unsigned int) gsl_histogram_sum(histogram_out);
    timestep = (double)sum / HISTOGRAM_TIMESTEP;
    sprintf(fileName, "hist-%d-all-out-timeseries", getpid());
    outFile = fopen(fileName, "a");
    fprintfHistogramTimeseries(outFile, hist_step_out, timestep);
    fclose(outFile);

    sprintf(fileName, "hll-%d-out-timeseries", getpid());
    outFile = fopen(fileName, "a");
    card_step = hll_step_out->Estimate();
    bits_step = log2(card_step);
    density_step = bits_step > 32 ? bits_step / 64 : bits_step / 32;
    card = hll_out->Estimate();
    bits = log2(card);
    density = bits > 32 ? bits / 64 : bits / 32;
    fprintf(outFile, "%g %lu %g %g %lu %g %g\n", timestep, card_step, bits_step, density_step, card, bits, density);
    fclose(outFile);

    sprintf(fileName, "hist-%d-all-in", getpid());
    //printf("writing %s\n", fileName);
    outFile = fopen(fileName, "w");
    fprintfHistogram(outFile, histogram_in, min_in, max_in);
    fclose(outFile);

    sprintf(fileName, "hist-%d-all-out", getpid());
    //printf("writing %s\n", fileName);
    outFile = fopen(fileName, "w");
    fprintfHistogram(outFile, histogram_out, min_out, max_out);
    fclose(outFile);

    sprintf(fileName, "hll-%d-all-timeseries", getpid());
    printf("writing %s\n", fileName);
    outFile = fopen(fileName, "a");
    card_step = hll_step_all->Estimate();
    bits_step = log2(card_step);
    density_step = bits_step > 32 ? bits_step / 64 : bits_step / 32;
    card = hll_all->Estimate();
    bits = log2(card);
    density = bits > 32 ? bits / 64 : bits / 32;
    fprintf(outFile, "%g %lu %g %g %lu %g %g\n", timestep, card_step, bits_step, density_step, card, bits, density);
    fclose(outFile);
}

}

