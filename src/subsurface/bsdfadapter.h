/// Adapter to use BSDFs in the chi-square test
class BSDFAdapter {
public:
    BSDFAdapter(const BSDF *bsdf, Sampler *sampler, const Vector &wi, 
            int component, bool passSamplerToBSDF)
        : m_bsdf(bsdf), m_sampler(sampler), m_wi(wi), m_component(component),
          m_largestWeight(0), m_passSamplerToBSDF(passSamplerToBSDF) {
        m_fakeSampler = new FakeSampler(m_sampler);
    }

    std::pair<Vector, Float> generateSample() {
        Point2 sample(m_sampler->next2D());
        Intersection its;
        BSDFQueryRecord bRec(its);
        bRec.component = m_component;
        bRec.wi = m_wi;
        
        #if defined(MTS_DEBUG_FP)
            enableFPExceptions();
        #endif

        Float pdfVal;

        /* Only make the sampler available to the BSDF when requested
           by the testcase. This allows testing both sampling variants
           where applicable: those that can improve by having access to 
           an arbitrary random number stream vs. those that only use
           a single uniform 2D sample */

        if (m_passSamplerToBSDF)
            bRec.sampler = m_fakeSampler;

        /* Check the various sampling routines for agreement amongst each other */
        m_fakeSampler->clear();
        Spectrum f = m_bsdf->sample(bRec, pdfVal, sample);
        m_fakeSampler->rewind();
        Spectrum sampled = m_bsdf->sample(bRec, sample);

        if (f.isZero() || pdfVal == 0) {
            if (!sampled.isZero()) 
                Log(EWarn, "Inconsistency (1): f=%s, pdf=%f, sampled f/pdf=%s, bRec=%s",
                    f.toString().c_str(), pdfVal, sampled.toString().c_str(), bRec.toString().c_str());
            #if defined(MTS_DEBUG_FP)
                disableFPExceptions();
            #endif
            return std::make_pair(bRec.wo, 0.0f);
        } else if (sampled.isZero()) {
            if (!f.isZero() && pdfVal != 0)
                Log(EWarn, "Inconsistency (2): f=%s, pdf=%f, sampled f/pdf=%s, bRec=%s",
                    f.toString().c_str(), pdfVal, sampled.toString().c_str(), bRec.toString().c_str());
            #if defined(MTS_DEBUG_FP)
                disableFPExceptions();
            #endif
            return std::make_pair(bRec.wo, 0.0f);
        }

        Spectrum sampled2 = f/pdfVal;
        if (!sampled.isValid() || !sampled2.isValid()) {
            Log(EWarn, "Ooops: f=%s, pdf=%f, sampled f/pdf=%s, bRec=%s",
                f.toString().c_str(), pdfVal, sampled.toString().c_str(), bRec.toString().c_str());
            return std::make_pair(bRec.wo, 0.0f);
        }

        bool mismatch = false;
        for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
            Float a = sampled[i], b = sampled2[i];
            Float min = std::min(a, b);
            Float err = std::abs(a - b);
            m_largestWeight = std::max(m_largestWeight, a * std::abs(Frame::cosTheta(bRec.wo)));

            if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
                mismatch = true;
            else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
                mismatch = true;
        }

        if (mismatch)
            Log(EWarn, "Potential inconsistency (3): f/pdf=%s, sampled f/pdf=%s",
                sampled2.toString().c_str(), sampled.toString().c_str());
        
        #if defined(MTS_DEBUG_FP)
            disableFPExceptions();
        #endif

        return std::make_pair(bRec.wo, 1.0f);
    }

    Float pdf(const Vector &wo) {
        Intersection its;
        BSDFQueryRecord bRec(its);
        bRec.component = m_component;
        bRec.wi = m_wi;
        bRec.wo = wo;
        if (m_passSamplerToBSDF)
            bRec.sampler = m_sampler;

        #if defined(MTS_DEBUG_FP)
            enableFPExceptions();
        #endif

        if (m_bsdf->f(bRec).isZero())
            return 0.0f;
        Float result = m_bsdf->pdf(bRec);

        #if defined(MTS_DEBUG_FP)
            disableFPExceptions();
        #endif
        return result;
    }

    inline Float getLargestWeight() const { return m_largestWeight; }
private:
    ref<const BSDF> m_bsdf;
    ref<Sampler> m_sampler;
    ref<FakeSampler> m_fakeSampler;
    Vector m_wi;
    int m_component;
    Float m_largestWeight;
    bool m_passSamplerToBSDF;
};
