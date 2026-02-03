#include <string.h>
#include "zint.h"
#include "modp.h"

uint32_t zint_sub(uint32_t *restrict a, const uint32_t *restrict b, size_t len,
	uint32_t ctl)
{
	size_t u;
	uint32_t cc, m;

	cc = 0;
	m = -ctl;
	for (u = 0; u < len; u ++) {
		uint32_t aw, w;

		aw = a[u];
		w = aw - b[u] - cc;
		cc = w >> 31;
		aw ^= ((w & 0x7FFFFFFF) ^ aw) & m;
		a[u] = aw;
	}
	return cc;
}


uint32_t zint_mul_small(uint32_t *m, size_t mlen, uint32_t x)
{
	size_t u;
	uint32_t cc;

	cc = 0;
	for (u = 0; u < mlen; u ++) {
		uint64_t z;

		z = (uint64_t)m[u] * (uint64_t)x + cc;
		m[u] = (uint32_t)z & 0x7FFFFFFF;
		cc = (uint32_t)(z >> 31);
	}
	return cc;
}


uint32_t zint_mod_small_unsigned(const uint32_t *d, size_t dlen,
	uint32_t p, uint32_t p0i, uint32_t R2)
{
	uint32_t x;
	size_t u;

	x = 0;
	u = dlen;
	while (u -- > 0) {
		uint32_t w;

		x = modp_montymul(x, R2, p, p0i);
		w = d[u] - p;
		w += p & -(w >> 31);
		x = modp_add(x, w, p);
	}
	return x;
}


uint32_t zint_mod_small_signed(const uint32_t *d, size_t dlen,
	uint32_t p, uint32_t p0i, uint32_t R2, uint32_t Rx)
{
	uint32_t z;

	if (dlen == 0) {
		return 0;
	}
	z = zint_mod_small_unsigned(d, dlen, p, p0i, R2);
	z = modp_sub(z, Rx & -(d[dlen - 1] >> 30), p);
	return z;
}


void zint_add_mul_small(uint32_t *restrict x,
	const uint32_t *restrict y, size_t len, uint32_t s)
{
	size_t u;
	uint32_t cc;

	cc = 0;
	for (u = 0; u < len; u ++) {
		uint32_t xw, yw;
		uint64_t z;

		xw = x[u];
		yw = y[u];
		z = (uint64_t)yw * (uint64_t)s + (uint64_t)xw + (uint64_t)cc;
		x[u] = (uint32_t)z & 0x7FFFFFFF;
		cc = (uint32_t)(z >> 31);
	}
	x[len] = cc;
}


void zint_norm_zero(uint32_t *restrict x, const uint32_t *restrict p, size_t len)
{
	size_t u;
	uint32_t r, bb;

	r = 0;
	bb = 0;
	u = len;
	while (u -- > 0) {
		uint32_t wx, wp, cc;

		wx = x[u];
		wp = (p[u] >> 1) | (bb << 30);
		bb = p[u] & 1;

		cc = wp - wx;
		cc = ((-cc) >> 31) | -(cc >> 31);

		r |= cc & ((r & 1) - 1);
	}

	zint_sub(x, p, len, r >> 31);
}


void zint_rebuild_CRT(uint32_t *restrict xx, size_t xlen, size_t xstride,
	size_t num, const small_prime *primes, int normalize_signed,
	uint32_t *restrict tmp)
{
	size_t u;
	uint32_t *x;

	tmp[0] = primes[0].p;
	for (u = 1; u < xlen; u ++) {

		uint32_t p, p0i, s, R2;
		size_t v;

		p = primes[u].p;
		s = primes[u].s;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);

		for (v = 0, x = xx; v < num; v ++, x += xstride) {
			uint32_t xp, xq, xr;
			xp = x[u];
			xq = zint_mod_small_unsigned(x, u, p, p0i, R2);

			xr = modp_montymul(s, modp_sub(xp, xq, p), p, p0i);
			zint_add_mul_small(x, tmp, u, xr);
		}

		tmp[u] = zint_mul_small(tmp, u, p);
	}

	if (normalize_signed) {
		for (u = 0, x = xx; u < num; u ++, x += xstride) {
			zint_norm_zero(x, tmp, xlen);
		}
	}
}


void zint_negate(uint32_t *a, size_t len, uint32_t ctl)
{
	size_t u;
	uint32_t cc, m;

	cc = ctl;
	m = -ctl >> 1;
	for (u = 0; u < len; u ++) {
		uint32_t aw;

		aw = a[u];
		aw = (aw ^ m) + cc;
		a[u] = aw & 0x7FFFFFFF;
		cc = aw >> 31;
	}
}

uint32_t zint_co_reduce(uint32_t *a, uint32_t *b, size_t len,
	int64_t xa, int64_t xb, int64_t ya, int64_t yb)
{
	size_t u;
	int64_t cca, ccb;
	uint32_t nega, negb;

	cca = 0;
	ccb = 0;
	for (u = 0; u < len; u ++) {
		uint32_t wa, wb;
		uint64_t za, zb;

		wa = a[u];
		wb = b[u];
		za = wa * (uint64_t)xa + wb * (uint64_t)xb + (uint64_t)cca;
		zb = wa * (uint64_t)ya + wb * (uint64_t)yb + (uint64_t)ccb;
		if (u > 0) {
			a[u - 1] = (uint32_t)za & 0x7FFFFFFF;
			b[u - 1] = (uint32_t)zb & 0x7FFFFFFF;
		}
		cca = *(int64_t *)&za >> 31;
		ccb = *(int64_t *)&zb >> 31;
	}
	a[len - 1] = (uint32_t)cca;
	b[len - 1] = (uint32_t)ccb;

	nega = (uint32_t)((uint64_t)cca >> 63);
	negb = (uint32_t)((uint64_t)ccb >> 63);
	zint_negate(a, len, nega);
	zint_negate(b, len, negb);
	return nega | (negb << 1);
}

void zint_finish_mod(uint32_t *a, size_t len, const uint32_t *m, uint32_t neg)
{
	size_t u;
	uint32_t cc, xm, ym;

	cc = 0;
	for (u = 0; u < len; u ++) {
		cc = (a[u] - m[u] - cc) >> 31;
	}

	xm = -neg >> 1;
	ym = -(neg | (1 - cc));
	cc = neg;
	for (u = 0; u < len; u ++) {
		uint32_t aw, mw;

		aw = a[u];
		mw = (m[u] ^ xm) & ym;
		aw = aw - mw - cc;
		a[u] = aw & 0x7FFFFFFF;
		cc = aw >> 31;
	}
}


void zint_co_reduce_mod(uint32_t *a, uint32_t *b, const uint32_t *m, size_t len,
	uint32_t m0i, int64_t xa, int64_t xb, int64_t ya, int64_t yb)
{
	size_t u;
	int64_t cca, ccb;
	uint32_t fa, fb;

	cca = 0;
	ccb = 0;
	fa = ((a[0] * (uint32_t)xa + b[0] * (uint32_t)xb) * m0i) & 0x7FFFFFFF;
	fb = ((a[0] * (uint32_t)ya + b[0] * (uint32_t)yb) * m0i) & 0x7FFFFFFF;
	for (u = 0; u < len; u ++) {
		uint32_t wa, wb;
		uint64_t za, zb;

		wa = a[u];
		wb = b[u];
		za = wa * (uint64_t)xa + wb * (uint64_t)xb
			+ m[u] * (uint64_t)fa + (uint64_t)cca;
		zb = wa * (uint64_t)ya + wb * (uint64_t)yb
			+ m[u] * (uint64_t)fb + (uint64_t)ccb;
		if (u > 0) {
			a[u - 1] = (uint32_t)za & 0x7FFFFFFF;
			b[u - 1] = (uint32_t)zb & 0x7FFFFFFF;
		}
		cca = *(int64_t *)&za >> 31;
		ccb = *(int64_t *)&zb >> 31;
	}
	a[len - 1] = (uint32_t)cca;
	b[len - 1] = (uint32_t)ccb;

	zint_finish_mod(a, len, m, (uint32_t)((uint64_t)cca >> 63));
	zint_finish_mod(b, len, m, (uint32_t)((uint64_t)ccb >> 63));
}


int zint_bezout(uint32_t *restrict u, uint32_t *restrict v,
	const uint32_t *restrict x, const uint32_t *restrict y,
	size_t len, uint32_t *restrict tmp)
{
	uint32_t *u0, *u1, *v0, *v1, *a, *b;
	uint32_t x0i, y0i;
	uint32_t num, rc;
	size_t j;

	if (len == 0) {
		return 0;
	}


	u0 = u;
	v0 = v;
	u1 = tmp;
	v1 = u1 + len;
	a = v1 + len;
	b = a + len;


	x0i = modp_ninv31(x[0]);
	y0i = modp_ninv31(y[0]);


	memcpy(a, x, len * sizeof *x);
	memcpy(b, y, len * sizeof *y);
	u0[0] = 1;
	memset(u0 + 1, 0, (len - 1) * sizeof *u0);
	memset(v0, 0, len * sizeof *v0);
	memcpy(u1, y, len * sizeof *u1);
	memcpy(v1, x, len * sizeof *v1);
	v1[0] --;

	for (num = 62 * (uint32_t)len + 30; num >= 30; num -= 30) {
		uint32_t c0, c1;
		uint32_t a0, a1, b0, b1;
		uint64_t a_hi, b_hi;
		uint32_t a_lo, b_lo;
		int64_t pa, pb, qa, qb;
		int i;
		uint32_t r;

		c0 = (uint32_t)-1;
		c1 = (uint32_t)-1;
		a0 = 0;
		a1 = 0;
		b0 = 0;
		b1 = 0;
		j = len;
		while (j -- > 0) {
			uint32_t aw, bw;

			aw = a[j];
			bw = b[j];
			a0 ^= (a0 ^ aw) & c0;
			a1 ^= (a1 ^ aw) & c1;
			b0 ^= (b0 ^ bw) & c0;
			b1 ^= (b1 ^ bw) & c1;
			c1 = c0;
			c0 &= (((aw | bw) + 0x7FFFFFFF) >> 31) - (uint32_t)1;
		}

		a1 |= a0 & c1;
		a0 &= ~c1;
		b1 |= b0 & c1;
		b0 &= ~c1;
		a_hi = ((uint64_t)a0 << 31) + a1;
		b_hi = ((uint64_t)b0 << 31) + b1;
		a_lo = a[0];
		b_lo = b[0];

		pa = 1;
		pb = 0;
		qa = 0;
		qb = 1;
		for (i = 0; i < 31; i ++) {

			uint32_t rt, oa, ob, cAB, cBA, cA;
			uint64_t rz;


			rz = b_hi - a_hi;
			rt = (uint32_t)((rz ^ ((a_hi ^ b_hi)
				& (a_hi ^ rz))) >> 63);

			oa = (a_lo >> i) & 1;
			ob = (b_lo >> i) & 1;
			cAB = oa & ob & rt;
			cBA = oa & ob & ~rt;
			cA = cAB | (oa ^ 1);

			a_lo -= b_lo & -cAB;
			a_hi -= b_hi & -(uint64_t)cAB;
			pa -= qa & -(int64_t)cAB;
			pb -= qb & -(int64_t)cAB;
			b_lo -= a_lo & -cBA;
			b_hi -= a_hi & -(uint64_t)cBA;
			qa -= pa & -(int64_t)cBA;
			qb -= pb & -(int64_t)cBA;

			a_lo += a_lo & (cA - 1);
			pa += pa & ((int64_t)cA - 1);
			pb += pb & ((int64_t)cA - 1);
			a_hi ^= (a_hi ^ (a_hi >> 1)) & -(uint64_t)cA;
			b_lo += b_lo & -cA;
			qa += qa & -(int64_t)cA;
			qb += qb & -(int64_t)cA;
			b_hi ^= (b_hi ^ (b_hi >> 1)) & ((uint64_t)cA - 1);
		}

		r = zint_co_reduce(a, b, len, pa, pb, qa, qb);
		pa -= (pa + pa) & -(int64_t)(r & 1);
		pb -= (pb + pb) & -(int64_t)(r & 1);
		qa -= (qa + qa) & -(int64_t)(r >> 1);
		qb -= (qb + qb) & -(int64_t)(r >> 1);
		zint_co_reduce_mod(u0, u1, y, len, y0i, pa, pb, qa, qb);
		zint_co_reduce_mod(v0, v1, x, len, x0i, pa, pb, qa, qb);
	}

	rc = a[0] ^ 1;
	for (j = 1; j < len; j ++) {
		rc |= a[j];
	}
	return (int)((1 - ((rc | -rc) >> 31)) & x[0] & y[0]);
}

void zint_add_scaled_mul_small(uint32_t *restrict x, size_t xlen,
	const uint32_t *restrict y, size_t ylen, int32_t k,
	uint32_t sch, uint32_t scl)
{
	size_t u;
	uint32_t ysign, tw;
	int32_t cc;

	if (ylen == 0) {
		return;
	}

	ysign = -(y[ylen - 1] >> 30) >> 1;
	tw = 0;
	cc = 0;
	for (u = sch; u < xlen; u ++) {
		size_t v;
		uint32_t wy, wys, ccu;
		uint64_t z;

		v = u - sch;
		wy = v < ylen ? y[v] : ysign;
		wys = ((wy << scl) & 0x7FFFFFFF) | tw;
		tw = wy >> (31 - scl);


		z = (uint64_t)((int64_t)wys * (int64_t)k + (int64_t)x[u] + cc);
		x[u] = (uint32_t)z & 0x7FFFFFFF;
		ccu = (uint32_t)(z >> 31);
		cc = *(int32_t *)&ccu;
	}
}

void zint_sub_scaled(uint32_t *restrict x, size_t xlen,
	const uint32_t *restrict y, size_t ylen, uint32_t sch, uint32_t scl)
{
	size_t u;
	uint32_t ysign, tw;
	uint32_t cc;

	if (ylen == 0) {
		return;
	}

	ysign = -(y[ylen - 1] >> 30) >> 1;
	tw = 0;
	cc = 0;
	for (u = sch; u < xlen; u ++) {
		size_t v;
		uint32_t w, wy, wys;

		/*
		 * Get the next word of y (scaled).
		 */
		v = u - sch;
		wy = v < ylen ? y[v] : ysign;
		wys = ((wy << scl) & 0x7FFFFFFF) | tw;
		tw = wy >> (31 - scl);

		w = x[u] - wys - cc;
		x[u] = w & 0x7FFFFFFF;
		cc = w >> 31;
	}
}


inline int32_t zint_one_to_plain(const uint32_t *x)
{
	uint32_t w;

	w = x[0];
	w |= (w & 0x40000000) << 1;
	return *(int32_t *)&w;
}


void poly_big_to_fp(double *d, const uint32_t *f, size_t flen, size_t fstride,
	unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	if (flen == 0) {
		for (u = 0; u < n; u ++) {
			d[u] = fpr_zero;
		}
		return;
	}
	for (u = 0; u < n; u ++, f += fstride) {
		size_t v;
		uint32_t neg, cc, xm;
		double x, fsc;

		/*
		 * Get sign of the integer; if it is negative, then we
		 * will load its absolute value instead, and negate the
		 * result.
		 */
		neg = -(f[flen - 1] >> 30);
		xm = neg >> 1;
		cc = neg & 1;
		x = fpr_zero;
		fsc = fpr_one;
		for (v = 0; v < flen; v ++, fsc = fpr_mul(fsc, fpr_ptwo31)) {
			uint32_t w;

			w = (f[v] ^ xm) + cc;
			cc = w >> 31;
			w &= 0x7FFFFFFF;
			w -= (w << 1) & neg;
			x = fpr_add(x, fpr_mul(fpr_of(*(int32_t *)&w), fsc));
		}
		d[u] = x;
	}
}


int poly_big_to_small(int8_t *d, const uint32_t *s, int lim, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		int32_t z;

		z = zint_one_to_plain(s + u);
		if (z < -lim || z > lim) {
			return 0;
		}
		d[u] = (int8_t)z;
	}
	return 1;
}