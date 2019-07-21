#pragma once
#include "bbd_filter.h"
#include <algorithm>
#include <vector>
#include <memory>
#include <complex>

class BBD_Line {
public:
    explicit BBD_Line(double fs, unsigned ns, const BBD_Filter_Spec &fsin, const BBD_Filter_Spec &fsout);
    void set_delay_size(unsigned ns);
    void clear();
    void process(unsigned n, float *inout, const float *clock);

    const BBD_Filter_Coef &filter_in() const noexcept { return *fin_; }
    const BBD_Filter_Coef &filter_out() const noexcept { return *fout_; }

private:
    const double fs_; // sampling frequency
    unsigned ns_; // delay size
    std::vector<float> mem_; // delay memory
    unsigned imem_; // delay memory index
    double pclk_; // clock phase
    unsigned ptick_; // clock tick counter
    double ybbd_old_;
    const BBD_Filter_Coef *fin_;
    const BBD_Filter_Coef *fout_;
    std::unique_ptr<double[]> Xin_; // sample from input filter
    std::unique_ptr<double[]> Xout0_; // 1st sample to output filter
    std::unique_ptr<double[]> Xout1_; // 2nd sample to output filter
    std::unique_ptr<double[]> Z1in_; // 1st sample memory of input filter
    std::unique_ptr<double[]> Z2in_; // 2nd sample memory of input filter
    std::unique_ptr<double[]> Z1out_; // 1st sample memory of output filter
    std::unique_ptr<double[]> Z2out_; // 2nd sample memory of output filter
    std::unique_ptr<double[]> B0_; // temporary buffer in which b0 coefficients are interpolated
    std::unique_ptr<double[]> B1_; // temporary buffer in which b1 coefficients are interpolated
};
