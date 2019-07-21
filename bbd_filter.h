#pragma once
#include <memory>
#include <complex>
typedef std::complex<double> cdouble;

enum class BBD_Filter_Kind : bool {
    Input,
    Output,
};

/*
  Analog specifications of BBD filters, input and output.
  M=order R=numerator P=denominator
  Analog transfer: H(s)=sum(m:1→M) (R[m]/(s-P[m]))
*/

struct BBD_Filter_Spec {
    BBD_Filter_Kind Kind;
    unsigned M;
    const cdouble *R;/*[M]*/
    const cdouble *P;/*[M]*/
};

/*
  Discretized matrix of filters coefficients.
  M=order, N=interpolation steps, H=feedback factor
*/
struct BBD_Filter_Coef {
    unsigned M;
    unsigned N;

    // complex-valued coefficients
    std::unique_ptr<cdouble[]> G;/*[M*N]*/
    std::unique_ptr<cdouble[]> P;/*[M]*/

    double H;

    //
    void interpolate_G(double d, cdouble *g/*[M]*/) const noexcept;

    // real-valued coefficients
    unsigned RealM;/*(M+1)/2*/
    std::unique_ptr<double[]> A1;/*[RealM]*/
    std::unique_ptr<double[]> A2;/*[RealM]*/
    std::unique_ptr<double[]> B0;/*[RealM*N]*/
    std::unique_ptr<double[]> B1;/*[RealM*N]*/

    //
    void interpolate_B0(double d, double *b0/*[M]*/) const noexcept;
    void interpolate_B1(double d, double *b1/*[M]*/) const noexcept;
};

namespace BBD {
BBD_Filter_Coef compute_filter(float fs, unsigned steps, const BBD_Filter_Spec &spec);
const BBD_Filter_Coef &compute_filter_cached(float fs, unsigned steps, const BBD_Filter_Spec &spec);
void clear_filter_cache();
} // namespace BBD

/*
  The model of BBD input and output filters from Juno 60.
*/
extern const BBD_Filter_Spec bbd_fin_j60;
extern const BBD_Filter_Spec bbd_fout_j60;
