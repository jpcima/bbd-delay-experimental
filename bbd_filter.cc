#include "bbd_filter.h"
#include <vector>
#include <mutex>
#include <algorithm>
#include <cassert>

template <class Float>
static bool fp_equal(Float a, Float b, Float tol = 1e-5)
{
    return std::fabs(a - b) < tol;
}

///

template <class T>
static void interpolate_row(double d, unsigned rows, unsigned cols, const T *src, T *dst)
{
    assert(d >= 0);
    double row = d * (rows - 1);
    unsigned row1 = std::min((unsigned)row, rows - 1);
    unsigned row2 = std::min(row1 + 1, rows - 1);
    double mu = row - (unsigned)row;
    for (unsigned i = 0; i < cols; ++i)
        dst[i] = (1 - mu) * src[row1 * cols + i] + mu * src[row2 * cols + i];
}

void BBD_Filter_Coef::interpolate_G(double d, cdouble *g/*[M]*/) const noexcept
{
    interpolate_row(d, N, M, G.get(), g);
}

void BBD_Filter_Coef::interpolate_B0(double d, double *b0/*[M]*/) const noexcept
{
    interpolate_row(d, N, RealM, B0.get(), b0);
}

void BBD_Filter_Coef::interpolate_B1(double d, double *b1/*[M]*/) const noexcept
{
    interpolate_row(d, N, RealM, B1.get(), b1);
}

// static double ensure_real(cdouble x, double tol = 1e-5)
// {
//     if (!fp_equal(x.imag(), 0.0, tol))
//         throw std::runtime_error("The number must be real");
//     return x.real();
// }

static void ensure_complex_conjugates(const cdouble *p, unsigned m)
{
    if (m % 2 == 0)
        throw std::runtime_error("The number of coefficients must be odd");
    if (!fp_equal(p[0].imag(), 0.0))
        throw std::runtime_error("The [0] coefficient must be real-valued.");
    for (unsigned i = 1; i < m; i += 2) {
        if (!fp_equal(p[i].real(), p[i + 1].real()) || !fp_equal(p[i].imag(), -p[i + 1].imag()))
            throw std::runtime_error("The [1:M] coefficients must be pairwise complex conjugates.");
    }
}

BBD_Filter_Coef BBD::compute_filter(float fs, unsigned steps, const BBD_Filter_Spec &spec)
{
    BBD_Filter_Coef coef;
    double ts = 1 / fs;
    unsigned M = spec.M;

    coef.M = M;
    coef.N = steps;
    coef.G.reset(new cdouble[M * steps]);
    coef.P.reset(new cdouble[M]);

    cdouble *pm = coef.P.get();
    for (unsigned m = 0; m < M; ++m)
        pm[m] = std::exp(ts * spec.P[m]);

    for (unsigned step = 0; step < steps; ++step) {
        double d = (double)step / (steps - 1);
        cdouble *gm = &coef.G[step * M];
        for (unsigned m = 0; m < M; ++m)
            gm[m] = ts * spec.R[m] * std::pow(pm[m], d);
    }

    cdouble H = 0;
    for (unsigned m = 0; m < M; ++m)
        H -= spec.R[m] / spec.P[m];
    coef.H = H.real();

    /////////////////////////////////////////
    // convert to real-valued coefficients //
    /////////////////////////////////////////

    // sanity checking
    ensure_complex_conjugates(&coef.P[0], M);
    for (unsigned step = 0; step < steps; ++step)
        ensure_complex_conjugates(&coef.G[M * step], M);

    unsigned RealM = coef.RealM = (M + 1) / 2;
    coef.A1.reset(new double[RealM]);
    coef.A2.reset(new double[RealM]);
    coef.B0.reset(new double[RealM * steps]);
    coef.B1.reset(new double[RealM * steps]);

    // the first coefficient, real-only
    coef.A1[0] = coef.P[0].real();
    coef.A2[0] = 0;
    for (unsigned step = 0; step < steps; ++step) {
        double *b0 = &coef.B0[step * RealM];
        double *b1 = &coef.B1[step * RealM];
        cdouble *gm = &coef.G[step * M];
        b0[0] = gm[0].real();
        b1[0] = 0;
    }
    // the other coefficients, complex conjugates
    for (unsigned m = 1; m < RealM; ++m) {
        double p = coef.P[(m - 1) * 2 + 1].real();
        coef.A1[m] = 2 * std::cos(std::arg(p));
        coef.A2[m] = -(std::abs(p) * std::abs(p));
    }
    for (unsigned step = 0; step < steps; ++step) {
        double *b0 = &coef.B0[step * RealM];
        double *b1 = &coef.B1[step * RealM];
        for (unsigned m = 1; m < RealM; ++m) {
            cdouble p = coef.P[(m - 1) * 2 + 1];
            cdouble r = spec.R[m];
            double d = (double)step / (steps - 1);
            if (spec.Kind == BBD_Filter_Kind::Input) {
                double beta = 2 * ts * std::abs(r);
                b0[m] = beta * std::pow(std::abs(p), d) * std::cos(std::arg(r) + d * std::arg(p));
                b1[m] = -beta * std::pow(std::abs(p), d + 1) * std::cos(std::arg(r) + (d - 1) * std::arg(p));
            }
            else {
                double beta = 2 * std::abs(r / p);
                b0[m] = beta * std::pow(std::abs(p), 1 - d) * std::cos(std::arg(r) + (1 - d) * std::arg(p));
                b1[m] = -beta * std::pow(std::abs(p), 2 - d) * std::cos(std::arg(r) - d * std::arg(p));
            }
        }
    }

    return coef;
}

//------------------------------------------------------------------------------

struct Filter_Cache_Entry {
    float fs;
    unsigned steps;
    const BBD_Filter_Spec *spec;
    BBD_Filter_Coef coef;
};
static std::vector<std::unique_ptr<Filter_Cache_Entry>> filter_cache;
static std::mutex filter_cache_mutex;

const BBD_Filter_Coef &BBD::compute_filter_cached(float fs, unsigned steps, const BBD_Filter_Spec &spec)
{
    std::unique_lock<std::mutex> lock(filter_cache_mutex);

    for (const std::unique_ptr<Filter_Cache_Entry> &ent : filter_cache) {
        if (ent->fs == fs && ent->steps == steps && ent->spec == &spec)
            return ent->coef;
    }
    lock.unlock();

    std::unique_ptr<Filter_Cache_Entry> ent(new Filter_Cache_Entry);
    BBD_Filter_Coef &coef = ent->coef;

    ent->fs = fs;
    ent->steps = steps;
    ent->spec = &spec;
    coef = compute_filter(fs, steps, spec);

    lock.lock();
    filter_cache.emplace_back(std::move(ent));
    return coef;
}

void BBD::clear_filter_cache()
{
    std::lock_guard<std::mutex> lock(filter_cache_mutex);
    filter_cache.clear();
}

//------------------------------------------------------------------------------

namespace j60 {
static constexpr unsigned M_in = 5;
static constexpr cdouble R_in[M_in] = {{251589, 0}, {-130428, -4165}, {-130428, 4165}, {4634, -22873}, {4634, 22873}};
static constexpr cdouble P_in[M_in] = {{-46580, 0}, {-55482, 25082}, {-55482, -25082}, {-26929, -59437}, {-26929, 59437}};
static constexpr unsigned M_out = 5;
static constexpr cdouble R_out[M_out] = {{5092, 0}, {11256, -99566}, {11256, 99566}, {-13802, -24606}, {-13802, 24606}};
static constexpr cdouble P_out[M_out] = {{-176261, 0}, {-51468, 21437}, {-51468, -21437}, {-26276, -59699}, {-26276, 59699}};
} // namespace j60

const BBD_Filter_Spec bbd_fin_j60 = {BBD_Filter_Kind::Input, j60::M_in, j60::R_in, j60::P_in};
const BBD_Filter_Spec bbd_fout_j60 = {BBD_Filter_Kind::Output, j60::M_out, j60::R_out, j60::P_out};
