#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
enum Constants { M_MAG_N = 32 };

std::array<std::complex<float>, M_MAG_N> m_fft_in = {
    {(0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f)}};
std::array<std::complex<float>, M_MAG_N> m_fft_out = {
    {(0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f)}};
std::array<float, M_MAG_N> m_fft_out_mag = {
    {(0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f), (0.0e+0f),
     (0.0e+0f), (0.0e+0f)}};

template <std::size_t N>
void ft(const std::array<std::complex<float>, N> &in,
        std::array<std::complex<float>, N> &out) {
  for (unsigned int k = 0; (k < N); k += 1) {
    out[k] = 0;
  }
  for (unsigned int k = 0; (k < N); k += 1) {
    for (unsigned int n = 0; (n < N); n += 1) {
      out[k] +=
          (std::exp(std::complex<float>((0.0e+0f), ((M_PI * -2 * k * n) / N))) *
           in[n]);
    }
  }
}
inline template <std::size_t N>
void bit_reverse_copy(const std::array<std::complex<float>, N> &in,
                      std::array<std::complex<float>, N> &out) {
  out[0] = in[0];
  out[16] = in[1];
  out[8] = in[2];
  out[24] = in[3];
  out[4] = in[4];
  out[20] = in[5];
  out[12] = in[6];
  out[28] = in[7];
  out[2] = in[8];
  out[18] = in[9];
  out[10] = in[10];
  out[26] = in[11];
  out[6] = in[12];
  out[22] = in[13];
  out[14] = in[14];
  out[30] = in[15];
  out[1] = in[16];
  out[17] = in[17];
  out[9] = in[18];
  out[25] = in[19];
  out[5] = in[20];
  out[21] = in[21];
  out[13] = in[22];
  out[29] = in[23];
  out[3] = in[24];
  out[19] = in[25];
  out[11] = in[26];
  out[27] = in[27];
  out[7] = in[28];
  out[23] = in[29];
  out[15] = in[30];
  out[31] = in[31];
}
template <std::size_t N>
void iterative_fft(const std::array<std::complex<float>, N> &in,
                   std::array<std::complex<float>, N> &out) {
  bit_reverse_copy(in, out);
  {
    auto w_m = const std::complex<float>((0.0e+0f), (-6.283185307179586e+0));
    for (auto k = 0; (k < N); k += 1) {
      {
        std::complex<float> w(1);
        for (unsigned int j = 0; (j <); j += 1) {
          {
            auto t((w * out[(k + j + nil)]));
            auto u(out[(k + j)]);
            out[(k + j)] = (u + t);
            out[(k + j + nil)] = (u - t);
            w = (w * w_m);
          }
        }
      }
    }
  }
  {
    auto w_m = const std::complex<float>((0.0e+0f), (-3.141592653589793e+0));
    for (auto k = 0; (k < N); k += 2) {
      {
        std::complex<float> w(1);
        for (unsigned int j = 0; (j < 1); j += 1) {
          {
            auto t((w * out[(k + j + 1)]));
            auto u(out[(k + j)]);
            out[(k + j)] = (u + t);
            out[(k + j + 1)] = (u - t);
            w = (w * w_m);
          }
        }
      }
    }
  }
  {
    auto w_m = const std::complex<float>((0.0e+0f), (-1.5707963267948966e+0));
    for (auto k = 0; (k < N); k += 4) {
      {
        std::complex<float> w(1);
        for (unsigned int j = 0; (j < 2); j += 1) {
          {
            auto t((w * out[(k + j + 2)]));
            auto u(out[(k + j)]);
            out[(k + j)] = (u + t);
            out[(k + j + 2)] = (u - t);
            w = (w * w_m);
          }
        }
      }
    }
  }
  {
    auto w_m = const std::complex<float>((0.0e+0f), (-7.853981633974483e-1));
    for (auto k = 0; (k < N); k += 8) {
      {
        std::complex<float> w(1);
        for (unsigned int j = 0; (j < 4); j += 1) {
          {
            auto t((w * out[(k + j + 4)]));
            auto u(out[(k + j)]);
            out[(k + j)] = (u + t);
            out[(k + j + 4)] = (u - t);
            w = (w * w_m);
          }
        }
      }
    }
  }
  {
    auto w_m = const std::complex<float>((0.0e+0f), (-3.9269908169872414e-1));
    for (auto k = 0; (k < N); k += 16) {
      {
        std::complex<float> w(1);
        for (unsigned int j = 0; (j < 8); j += 1) {
          {
            auto t((w * out[(k + j + 8)]));
            auto u(out[(k + j)]);
            out[(k + j)] = (u + t);
            out[(k + j + 8)] = (u - t);
            w = (w * w_m);
          }
        }
      }
    }
  }
}
template <std::size_t N> void fft(std::array<std::complex<float>, N> &zs) {
  {
    unsigned int j(0);
    // bit reversal;
    for (unsigned int i = 0; (i < (N - 1)); i += 1) {
      if ((i < j)) {
        std::swap(zs[i], zs[j]);
      }
      {
        int m((N / 2));
        j -= m;
        while ((0 == (j & m))) {
          m /= 2;
          j -= m;
        }
      }
    }
  }
  for (auto j = 1; (j < N); j *= 2) {
    for (unsigned int m = 0; (m < j); m += 1) {
      {
        const auto sign(1);
        float t(((M_PI * sign * m) / j));
        auto w(std::exp(std::complex<float>((0.0e+0f), t)));
        for (auto i = m; (i < N); i += (2 * j)) {
          {
            auto zi(zs[i]);
            auto tz((w * zs.at((i + j))));
            zs[i] = (zi + std::complex<float>(t));
            zs.at((i + j)) = (zi - t);
          }
        }
      }
    }
  }
}
int main() {
  for (unsigned int i = 0; (i < M_MAG_N); i += 1) {
    m_fft_in[i] = (0.0e+0f);
    m_fft_out[i] = (0.0e+0f);
    m_fft_out_mag[i] = (0.0e+0f);
  }
  m_fft_in[12] = (1.e+0f);
  ft(m_fft_in, m_fft_out);
  for (unsigned int i = 0; (i < M_MAG_N); i += 1) {
    m_fft_in[i] = (0.0e+0f);
    m_fft_out[i] = (0.0e+0f);
    m_fft_out_mag[i] = (0.0e+0f);
  }
  m_fft_in[12] = (1.e+0f);
  fft(m_fft_in);
  for (unsigned int i = 0; (i < M_MAG_N); i += 1) {
    (std::cout << i << " " << m_fft_out[i] << m_fft_in[i] << " "
               << std::exp(std::complex<float>((0.0e+0f),
                                               ((M_PI * -2 * i) / M_MAG_N)))
               << std::endl);
  }
}