#ifndef ZALCON_SAMP_H
#define ZALCON_SAMP_H

int64_t SampleLW(void);
int64_t SampleC1(int64_t c_in);
int32_t SampleArbitraryCenter128(int128_t num, uint64_t den);
int32_t SampleArbitraryCenter64(int64_t num, uint64_t den_ignored);
int64_t SampleLW(void);

#endif