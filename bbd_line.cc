#include "bbd_line.h"
#include <algorithm>
#include <cassert>

static unsigned interp_size = 128;

BBD_Line::BBD_Line(double fs, unsigned ns, const BBD_Filter_Spec &fsin, const BBD_Filter_Spec &fsout)
    : fs_(fs)
{
    mem_.reserve(8192);

    const BBD_Filter_Coef &fin = BBD::compute_filter_cached(fs, interp_size, fsin);
    const BBD_Filter_Coef &fout = BBD::compute_filter_cached(fs, interp_size, fsout);
    fin_ = &fin;
    fout_ = &fout;

    unsigned Min = fin.RealM;
    unsigned Mout = fout.RealM;
    Xin_.reset(new double[Min]);
    Xout0_.reset(new double[Mout]);
    Xout1_.reset(new double[Mout]);
    Z1in_.reset(new double[Min]);
    Z2in_.reset(new double[Min]);
    Z1out_.reset(new double[Mout]);
    Z2out_.reset(new double[Mout]);
    B0_.reset(new double[std::max(Min, Mout)]);
    B1_.reset(new double[std::max(Min, Mout)]);

    set_delay_size(ns);
    clear();
}

void BBD_Line::set_delay_size(unsigned ns)
{
    mem_.clear();
    mem_.resize(ns);
    imem_ = 0;
    ns_ = ns;
}

void BBD_Line::clear()
{
    std::fill(mem_.begin(), mem_.end(), 0);
    imem_ = 0;
    pclk_ = 0;
    ptick_ = 0;
    ybbd_old_ = 0;
    unsigned Min = fin_->RealM;
    unsigned Mout = fout_->RealM;
    std::fill(&Xin_[0], &Xin_[Min], 0);
    std::fill(&Xout0_[0], &Xout0_[Mout], 0);
    std::fill(&Xout1_[0], &Xout1_[Mout], 0);
    std::fill(&Z1in_[0], &Z1in_[Min], 0);
    std::fill(&Z2in_[0], &Z2in_[Min], 0);
    std::fill(&Z1out_[0], &Z1out_[Mout], 0);
    std::fill(&Z2out_[0], &Z2out_[Mout], 0);
}

void BBD_Line::process(unsigned n, float *inout, const float *clock)
{
    unsigned ns = ns_;
    float *mem = mem_.data();
    unsigned imem = imem_;
    double pclk = pclk_;
    unsigned ptick = ptick_;
    double ybbd_old = ybbd_old_;

    const BBD_Filter_Coef &fin = *fin_, &fout = *fout_;
    unsigned Min = fin.RealM, Mout = fout.RealM;
    double *Xin = Xin_.get(), *Xout0 = Xout0_.get(), *Xout1 = Xout1_.get();
    double *Z1in = Z1in_.get(), *Z2in = Z2in_.get();
    double *Z1out = Z1out_.get(), *Z2out = Z2out_.get();
    const double *A1in = fin.A1.get(), *A1out = fout.A1.get();
    const double *A2in = fin.A2.get(), *A2out = fout.A2.get();
    double *B0 = B0_.get(), *B1 = B1_.get();

    for (unsigned i = 0; i < n; ++i) {
        double fclk = clock[i];

        for (unsigned m = 0; m < Mout; ++m) {
            Xout0[m] = 0;
            Xout1[m] = 0;
        }

        if (fclk > 0) {
            double pclk_old = pclk;
            pclk += fclk;
            unsigned tick_count = (unsigned)pclk;
            pclk -= tick_count;
            for (unsigned tick = 0; tick < tick_count; ++tick) {
                double d = (1 - pclk_old + tick) * (1 / fclk);
                d -= (unsigned)d;
                if ((ptick & 1) == 0) {
                    fin.interpolate_B0(d, B0);
                    fin.interpolate_B1(d, B1);
                    double s = 0;
                    for (unsigned m = 0; m < Min; ++m) {
                        s += Xin[m] * B0[m];
                        s += Z1in[m] * B1[m];
                    }
                    mem[imem] = s;
                    imem = ((imem + 1) < ns) ? (imem + 1) : 0;
                }
                else {
                    fout.interpolate_B0(d, B0);
                    fout.interpolate_B1(d, B1);
                    double ybbd = mem[imem];
                    double delta = ybbd - ybbd_old;
                    ybbd_old = ybbd;
                    for (unsigned m = 0; m < Mout; ++m) {
                        Xout0[m] += B0[m] * delta;
                        Xout1[m] += B1[m] * delta;
                    }
                }
                ++ptick;
            }
        }

        for (unsigned m = 0; m < Min; ++m) {
            double t1 = A1in[m] * Z1in[m];
            double t2 = A2in[m] * Z2in[m];
            Z2in[m] = Z1in[m];
            Z1in[m] = Xin[m];
            Xin[m] = t1 + t2 + inout[i];
        }

        double y = fout.H * ybbd_old;
        for (unsigned m = 0; m < Mout; ++m) {
            double t = Z1out[m] + Xout0[m];
            Z1out[m] = Xout1[m] + A1out[m] * t;
            Z2out[m] = A2out[m] * t;
            y += t;
        }

        inout[i] = y;
    }

    imem_ = imem;
    pclk_ = pclk;
    ptick_ = ptick;
    ybbd_old_ = ybbd_old;
}
