/* libmypaint - The MyPaint Brush Library
 * Copyright (C) 2007-2008 Martin Renold <martinxyz@gmx.ch>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#ifndef HELPERS_C
#define HELPERS_C

#include <config.h>

#include <assert.h>
#include <stdint.h>
#include <math.h>

#include "helpers.h"

float T_MATRIX[3][36] = {{5.47813E-05, 0.000184722, 0.000935514, 0.003096265, 0.009507714, 0.017351596, 0.022073595, 0.016353161, 0.002002407, -0.016177731, -0.033929391, -0.046158952, -0.06381706, -0.083911194, -0.091832385, -0.08258148, -0.052950086, -0.012727224, 0.037413037, 0.091701812, 0.147964686, 0.181542886, 0.210684154, 0.210058081, 0.181312094, 0.132064724, 0.093723787, 0.057159281, 0.033469657, 0.018235464, 0.009298756, 0.004023687, 0.002068643, 0.00109484, 0.000454231, 0.000255925},
{-4.65552E-05, -0.000157894, -0.000806935, -0.002707449, -0.008477628, -0.016058258, -0.02200529, -0.020027434, -0.011137726, 0.003784809, 0.022138944, 0.038965605, 0.063361718, 0.095981626, 0.126280277, 0.148575844, 0.149044804, 0.14239936, 0.122084916, 0.09544734, 0.067421931, 0.035691251, 0.01313278, -0.002384996, -0.009409573, -0.009888983, -0.008379513, -0.005606153, -0.003444663, -0.001921041, -0.000995333, -0.000435322, -0.000224537, -0.000118838, -4.93038E-05, -2.77789E-05},
{0.00032594, 0.001107914, 0.005677477, 0.01918448, 0.060978641, 0.121348231, 0.184875618, 0.208804428, 0.197318551, 0.147233899, 0.091819086, 0.046485543, 0.022982618, 0.00665036, -0.005816014, -0.012450334, -0.015524259, -0.016712927, -0.01570093, -0.013647887, -0.011317812, -0.008077223, -0.005863171, -0.003943485, -0.002490472, -0.001440876, -0.000852895, -0.000458929, -0.000248389, -0.000129773, -6.41985E-05, -2.71982E-05, -1.38913E-05, -7.35203E-06, -3.05024E-06, -1.71858E-06}};

float spectral_r[36] = {0.015603944435243, 0.015604290522261, 0.015605813688466, 0.015613374952138, 0.015641410348937, 0.015735182555386, 0.015963980154163, 
0.016413479597224, 0.017147856366573, 0.018216350851111, 0.019637333312242, 0.021433552239869, 0.023647340971203, 0.026423918541353, 
0.030050684039964, 0.035005703130505, 0.042098118739387, 0.052576115305368, 0.068652704521994, 0.094181329816583, 0.136134223719597, 
0.207308449035284, 0.328494607128972, 0.53005148318892, 0.838859915780766, 1.2465990882888, 1.69219313198829, 2.07675584506431, 
2.34910243251493, 2.51178891096814, 2.5974303446936, 2.63841433161008, 2.65888276907875, 2.6684963749706, 2.67229195424748, 2.67366156899577};

float spectral_g[36] = {0.047202088273281, 0.047206744553887, 0.047227230602721, 0.047328962072349, 0.04770743771169, 0.048991913602353, 0.052254925548439, 0.059258144233783, 0.072636771112054, 0.097300965154598, 0.142095959226619, 0.223445766596059, 0.369366506805273, 0.62287213146272, 1.00261672243663, 1.37101919554934, 1.44307825143871, 1.19020268861967, 0.850138320141449, 0.584247440487092, 0.411540207129038, 0.306351610576843, 0.242502999890411, 0.204181641595136, 0.181492529488003, 0.168348152754132, 0.160779551959378, 0.156644743701248, 0.154444959899402, 0.15332979076887, 0.152792337005596, 0.152546319127935, 0.152425840477336, 0.152369746608457, 0.152347684936107, 0.152339735430577};

float spectral_b[36] = {0.947081874376452, 0.947194900200616, 0.94767489139251, 0.949902270339965, 0.957380374473825, 0.978762576435505, 1.01726078899001, 1.05054089905291, 1.00787084940803, 0.821916356380445, 0.562943946673851, 0.345346773846951, 0.203906711045521, 0.120729683657514, 0.073588462825919, 0.047005376726015, 0.031749191014897, 0.022638662855379, 0.016998737102774, 0.013383881843011, 0.011007384216009, 0.009432817302141, 0.008375027656616, 0.0076802952958, 0.007239450229151, 0.006971693886012, 0.00681297038356, 0.006724786463696, 0.006677438813494, 0.006653323539191, 0.006641680210393, 0.006636347219181, 0.006633734150006, 0.006632517032454, 0.006632038255597, 0.006631865725982};

float rand_gauss (RngDouble * rng)
{
  double sum = 0.0;
  sum += rng_double_next(rng);
  sum += rng_double_next(rng);
  sum += rng_double_next(rng);
  sum += rng_double_next(rng);
  return sum * 1.73205080757 - 3.46410161514;
}

// C fmodf function is not "arithmetic modulo"; it doesn't handle negative dividends as you might expect
// if you expect 0 or a positive number when dealing with negatives, use
// this function instead.
float mod_arith(float a, float N)
{
    float ret = a - N * floor (a / N);
    return ret;
}

// Returns the smallest angular difference
float smallest_angular_difference(float angleA, float angleB)
{
    float a;
    a = angleB - angleA;
    a = mod_arith((a + 180), 360) - 180;
    a += (a>180) ? -360 : (a<-180) ? 360 : 0;
    return a;
}

// stolen from GIMP (gimpcolorspace.c)
// (from gimp_rgb_to_hsv)
void
rgb_to_hsv_float (float *r_ /*h*/, float *g_ /*s*/, float *b_ /*v*/)
{
  float max, min, delta;
  float h, s, v;
  float r, g, b;

  h = 0.0; // silence gcc warning

  r = *r_;
  g = *g_;
  b = *b_;

  r = CLAMP(r, 0.0, 1.0);
  g = CLAMP(g, 0.0, 1.0);
  b = CLAMP(b, 0.0, 1.0);

  max = MAX3(r, g, b);
  min = MIN3(r, g, b);

  v = max;
  delta = max - min;

  if (delta > 0.0001)
    {
      s = delta / max;

      if (r == max)
        {
          h = (g - b) / delta;
          if (h < 0.0)
            h += 6.0;
        }
      else if (g == max)
        {
          h = 2.0 + (b - r) / delta;
        }
      else if (b == max)
        {
          h = 4.0 + (r - g) / delta;
        }

      h /= 6.0;
    }
  else
    {
      s = 0.0;
      h = 0.0;
    }

  *r_ = h;
  *g_ = s;
  *b_ = v;
}

// (from gimp_hsv_to_rgb)
void
hsv_to_rgb_float (float *h_, float *s_, float *v_)
{
  int    i;
  double f, w, q, t;
  float h, s, v;
  float r, g, b;
  r = g = b = 0.0; // silence gcc warning

  h = *h_;
  s = *s_;
  v = *v_;

  h = h - floor(h);
  s = CLAMP(s, 0.0, 1.0);
  v = CLAMP(v, 0.0, 1.0);

  double hue;

  if (s == 0.0)
    {
      r = v;
      g = v;
      b = v;
    }
  else
    {
      hue = h;

      if (hue == 1.0)
        hue = 0.0;

      hue *= 6.0;

      i = (int) hue;
      f = hue - i;
      w = v * (1.0 - s);
      q = v * (1.0 - (s * f));
      t = v * (1.0 - (s * (1.0 - f)));

      switch (i)
        {
        case 0:
          r = v;
          g = t;
          b = w;
          break;
        case 1:
          r = q;
          g = v;
          b = w;
          break;
        case 2:
          r = w;
          g = v;
          b = t;
          break;
        case 3:
          r = w;
          g = q;
          b = v;
          break;
        case 4:
          r = t;
          g = w;
          b = v;
          break;
        case 5:
          r = v;
          g = w;
          b = q;
          break;
        }
    }

  *h_ = r;
  *s_ = g;
  *v_ = b;
}

// (from gimp_rgb_to_hsl)
void
rgb_to_hsl_float (float *r_, float *g_, float *b_)
{
  double max, min, delta;

  float h, s, l;
  float r, g, b;

  // silence gcc warnings
  h=0;

  r = *r_;
  g = *g_;
  b = *b_;

  r = CLAMP(r, 0.0, 1.0);
  g = CLAMP(g, 0.0, 1.0);
  b = CLAMP(b, 0.0, 1.0);

  max = MAX3(r, g, b);
  min = MIN3(r, g, b);

  l = (max + min) / 2.0;

  if (max == min)
    {
      s = 0.0;
      h = 0.0; //GIMP_HSL_UNDEFINED;
    }
  else
    {
      if (l <= 0.5)
        s = (max - min) / (max + min);
      else
        s = (max - min) / (2.0 - max - min);

      delta = max - min;

      if (delta == 0.0)
        delta = 1.0;

      if (r == max)
        {
          h = (g - b) / delta;
        }
      else if (g == max)
        {
          h = 2.0 + (b - r) / delta;
        }
      else if (b == max)
        {
          h = 4.0 + (r - g) / delta;
        }

      h /= 6.0;

      if (h < 0.0)
        h += 1.0;
    }

  *r_ = h;
  *g_ = s;
  *b_ = l;
}

static double
hsl_value (double n1,
           double n2,
           double hue)
{
  double val;

  if (hue > 6.0)
    hue -= 6.0;
  else if (hue < 0.0)
    hue += 6.0;

  if (hue < 1.0)
    val = n1 + (n2 - n1) * hue;
  else if (hue < 3.0)
    val = n2;
  else if (hue < 4.0)
    val = n1 + (n2 - n1) * (4.0 - hue);
  else
    val = n1;

  return val;
}


/**
 * gimp_hsl_to_rgb:
 * @hsl: A color value in the HSL colorspace
 * @rgb: The value converted to a value in the RGB colorspace
 *
 * Convert a HSL color value to an RGB color value.
 **/
void
hsl_to_rgb_float (float *h_, float *s_, float *l_)
{
  float h, s, l;
  float r, g, b;

  h = *h_;
  s = *s_;
  l = *l_;

  h = h - floor(h);
  s = CLAMP(s, 0.0, 1.0);
  l = CLAMP(l, 0.0, 1.0);

  if (s == 0)
    {
      /*  achromatic case  */
      r = l;
      g = l;
      b = l;
    }
  else
    {
      double m1, m2;

      if (l <= 0.5)
        m2 = l * (1.0 + s);
      else
        m2 = l + s - l * s;

      m1 = 2.0 * l - m2;

      r = hsl_value (m1, m2, h * 6.0 + 2.0);
      g = hsl_value (m1, m2, h * 6.0);
      b = hsl_value (m1, m2, h * 6.0 - 2.0);
    }

  *h_ = r;
  *s_ = g;
  *l_ = b;
}

void
rgb_to_hcy_float (float *r_, float *g_, float *b_) {
	
	float _HCY_RED_LUMA = 0.2162;
	float _HCY_GREEN_LUMA = 0.7152;
	float _HCY_BLUE_LUMA = 0.0722;
	float h, c, y;
	float r, g, b;
	float p, n, d;

	r = *r_;
	g = *g_;
	b = *b_;

	// Luma is just a weighted sum of the three components.
	y = _HCY_RED_LUMA*r + _HCY_GREEN_LUMA*g + _HCY_BLUE_LUMA*b;

	// Hue. First pick a sector based on the greatest RGB component, then add
	// the scaled difference of the other two RGB components.
	p = MAX3(r, g, b);
	n = MIN3(r, g, b);
	d = p - n; // An absolute measure of chroma: only used for scaling

	if (n == p){
		h = 0.0;
	} else if (p == r){
		h = (g - b)/d;
		if (h < 0){
			h += 6.0;
		}
	} else if (p == g){
		h = ((b - r)/d) + 2.0;
	} else {  // p==b
		h = ((r - g)/d) + 4.0;
	}
	h /= 6.0;
	h = fmod(h,1.0);

	// Chroma, relative to the RGB gamut envelope.
	if ((r == g) && (g == b)){
		// Avoid a division by zero for the achromatic case.
		c = 0.0;
	} else {
		// For the derivation, see the GLHS paper.
		c = MAX((y-n)/y, (p-y)/(1-y));
	}

	*r_ = h;
	*g_ = c;
	*b_ = y;
}

void
hcy_to_rgb_float (float *h_, float *c_, float *y_) {
	
	float _HCY_RED_LUMA = 0.2162;
	float _HCY_GREEN_LUMA = 0.7152;
	float _HCY_BLUE_LUMA = 0.0722;
	float h, c, y;
	float r, g, b;
	float th, tm;

	h = *h_;
	c = *c_;
	y = *y_;

	h = h - floor(h);
	c = CLAMP(c, 0.0f, 1.0f);
	y = CLAMP(y, 0.0f, 1.0f);

	if (c == 0)	{
	  /*  achromatic case  */
	  r = y;
	  g = y;
	  b = y;
	}

	h = fmod(h, 1.0);
	h *= 6.0;

	if (h < 1){
		// implies (p==r and h==(g-b)/d and g>=b)
		th = h;
		tm = _HCY_RED_LUMA + _HCY_GREEN_LUMA * th;
	} else if (h < 2) {
		// implies (p==g and h==((b-r)/d)+2.0 and b<r)
		th = 2.0 - h;
		tm = _HCY_GREEN_LUMA + _HCY_RED_LUMA * th;
	} else if (h < 3){
		// implies (p==g and h==((b-r)/d)+2.0 and b>=g)
		th = h - 2.0;
		tm = _HCY_GREEN_LUMA + _HCY_BLUE_LUMA * th;
	} else if (h < 4) {
		// implies (p==b and h==((r-g)/d)+4.0 and r<g)
		th = 4.0 - h;
		tm = _HCY_BLUE_LUMA + _HCY_GREEN_LUMA * th;
	} else if (h < 5){
		// implies (p==b and h==((r-g)/d)+4.0 and r>=g)
		th = h - 4.0;
		tm = _HCY_BLUE_LUMA + _HCY_RED_LUMA * th;
	} else {
		// implies (p==r and h==(g-b)/d and g<b)
		th = 6.0 - h;
		tm = _HCY_RED_LUMA + _HCY_BLUE_LUMA * th;
	}

	float n,p,o;
	// Calculate the RGB components in sorted order
	if (tm >= y){
		p = y + y*c*(1-tm)/tm;
		o = y + y*c*(th-tm)/tm;
		n = y - (y*c);
	}else{
		p = y + (1-y)*c;
		o = y + (1-y)*c*(th-tm)/(1-tm);
		n = y - (1-y)*c*tm/(1-tm);
	}

	// Back to RGB order
	if (h < 1){
		r = p;
		g = o;
		b = n;
	} else if (h < 2){
		r = o;
		g = p;
		b = n;
	} else if (h < 3){
		r = n;
		g = p;
		b = o;
	} else if (h < 4){
		r = n;
		g = o;
		b = p;
	} else if (h < 5){
		r = o;
		g = n;
		b = p;
	}else{ 
		r = p;
		g = n;
		b = o;
	}

	*h_ = r;
	*c_ = g;
	*y_ = b;
}


void
srgb_to_rgb_float (float *r_, float *g_, float *b_, float gamma) {
  
  float rgb[3] = {*r_, *g_, *b_};
  int i=0;
  for (i=0; i < 3; i++) { 
    rgb[i] = powf(rgb[i], gamma);
  }
  
  *r_ = rgb[0];
  *g_ = rgb[1];
  *b_ = rgb[2];
  
}

void
rgb_to_srgb_float (float *r_, float *g_, float *b_, float gamma) {

  float rgb[3] = {*r_, *g_, *b_};
  int i=0;
  for (i=0; i < 3; i++) {
    rgb[i] = powf(rgb[i],1/gamma);
  }
  *r_ = rgb[0];
  *g_ = rgb[1];
  *b_ = rgb[2];
}

//function to make it easy to blend normal and subtractive color blending modes w/ adjustable gamma
//a is the current smudge state, b is the get_color (get=1) or the brush color (get=0)
//mixing smudge_state+get_color is slightly different than mixing brush_color with smudge_color
//so I've used the bool'get' parameter to differentiate them.
float * mix_colors(float *a, float *b, float fac, float gamma, float normsub, gboolean get, float smudge_darken, float smudge_desat, float spectral)
{
  float rgbmixnorm[4] = {0};
  float rgbmixsub[4] = {0};
  float spectralmixnorm[4] = {0};
  float spectralmixsub[4] = {0};
  float normmix[4] = {0};
  float submix[4] = {0};
  static float result[4] = {0};
  //normsub is the ratio of normal to subtractive
  normsub = CLAMP(normsub, 0.0f, 1.0f);

  if (gamma < 1.0) {
    gamma = 1.0;
  }
  float ar=a[0];
  float ag=a[1];
  float ab=a[2];
  float aa=a[3];
  float br=b[0];
  float bg=b[1];
  float bb=b[2];
  float ba=b[3];

  //convert to linear rgb only if gamma is not 1.0 to avoid redundancy and rounding errors
  if (gamma != 1.0) {
    srgb_to_rgb_float(&ar, &ag, &ab, gamma);
    srgb_to_rgb_float(&br, &bg, &bb, gamma);
  }

  //do the mix
  //normal
  if (normsub < 1.0) {
    //premultiply alpha when getting color from canvas
    //the color in the smudge_state (a) is already premultiplied
    float bra =br;
    float bga =bg;
    float bba =bb;
    if (get == TRUE) {
      bra *=ba;
      bga *=ba;
      bba *=ba;
    }

    //RGB normal mixing
    if (spectral < 1.0) {

      rgbmixnorm[0] = fac * ar + (1-fac) * bra;
      rgbmixnorm[1] = fac * ag + (1-fac) * bga;
      rgbmixnorm[2] = fac * ab + (1-fac) * bba;
    }
    
    //do eraser_target alpha for smudge+brush color
    if (get == FALSE) {
      rgbmixnorm[0] /=ba;
      rgbmixnorm[1] /=ba;
      rgbmixnorm[2] /=ba;
    }    
    
    //spectral normal mixing
    if (spectral > 0.0) {
      
      //expand rgb to spectral primaries
      float spec_ar[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_ar[i] = spectral_r[i] * ar;
      }
      float spec_ag[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_ag[i] = spectral_g[i] * ag;
      }
      float spec_ab[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_ab[i] = spectral_b[i] * ab;
      }
      
      float spec_br[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_br[i] = spectral_r[i] * br;
      }
      float spec_bg[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_bg[i] = spectral_g[i] * bg;
      }
      float spec_bb[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_bb[i] = spectral_b[i] * bb;
      }
      
      //blend spectral primaries normal mode
      float spectralmix[3][36] = {0};
      for (int i=0; i < 36; i++) {
        spectralmix[0][i] = fac * spec_ar[i] + (1-fac) * spec_br[i];
        spectralmix[1][i] = fac * spec_ag[i] + (1-fac) * spec_bg[i];
        spectralmix[2][i] = fac * spec_ab[i] + (1-fac) * spec_bb[i];
      }
      //collapse into one SPD
      float spectral_combined[36] = {0};
      for (int i=0; i<36; i++) {
        spectral_combined[i] += spectralmix[0][i];
        spectral_combined[i] += spectralmix[1][i];
        spectral_combined[i] += spectralmix[2][i];
      }
      //convert to RGB
      for (int i=0; i<36; i++) {
        spectralmixnorm[0] += T_MATRIX[0][i] * spectral_combined[i];
        spectralmixnorm[1] += T_MATRIX[1][i] * spectral_combined[i];
        spectralmixnorm[2] += T_MATRIX[2][i] * spectral_combined[i];
      }
      spectralmixnorm[0] = CLAMP(spectralmixnorm[0], 0.0f, 1.0f);
      spectralmixnorm[1] = CLAMP(spectralmixnorm[1], 0.0f, 1.0f);
      spectralmixnorm[2] = CLAMP(spectralmixnorm[2], 0.0f, 1.0f);
    }
    
    //mix rgb and spectral normal modes:
    for (int i=0; i < 3; i++) {  
      normmix[i] = (((1-spectral)*rgbmixnorm[i]) + (spectral*spectralmixnorm[i]));
    }

  }

  //subtractive blend modes (spectral and rgb)
  if (normsub > 0.0) {
    float alpha_b;  //either brush_a or get_color_a
    if (get == FALSE) {
      alpha_b = 1;
    } else {
      alpha_b = ba;
    }
    
    //calculate a different smudge ratio for subtractive modes
    //may not be coherent but it looks good
    //this is a weighted geometric mean method inspired
    //by Scott Allen Burns (we're not using spectral though)
    //we want the sum of the ratios to be 1.0
    //but we're using alpha for that ratio so need to scale it
    float paint_a_ratio, paint_b_ratio;
    float subfac = fac;
    float alpha_sum = aa + alpha_b;
    if (alpha_sum >0.0) {
      paint_a_ratio = (aa/alpha_sum)*fac;
      paint_b_ratio = (ba/alpha_sum)*(1-fac);
      paint_a_ratio /=(paint_a_ratio+paint_b_ratio);
      paint_b_ratio /=(paint_a_ratio+paint_b_ratio);
    }
    
    float paint_ratio_sum = paint_a_ratio + paint_b_ratio;
    if (paint_ratio_sum > 0.0) {
      subfac = paint_a_ratio/(paint_ratio_sum);
    }
    
    if (spectral < 1.0) {
      //mix with Weighted Geometric Mean
      //don't allow absolute zero because it has infinite
      //darkening power
      rgbmixsub[0] = powf(MAX(ar, 0.0001), subfac) * powf(MAX(br, 0.0001),(1-subfac));
      rgbmixsub[1] = powf(MAX(ag, 0.0001), subfac) * powf(MAX(bg, 0.0001),(1-subfac));
      rgbmixsub[2] = powf(MAX(ab, 0.0001), subfac) * powf(MAX(bb, 0.0001),(1-subfac));
    }
    
    if (spectral > 0.0) {
      //expand rgb to spectral primaries
      float spec_ar[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_ar[i] = spectral_r[i] * ar;
      }
      float spec_ag[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_ag[i] = spectral_g[i] * ag;
      }
      float spec_ab[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_ab[i] = spectral_b[i] * ab;
      }
      
      float spec_br[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_br[i] = spectral_r[i] * br;
      }
      float spec_bg[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_bg[i] = spectral_g[i] * bg;
      }
      float spec_bb[36] = {0};
      for (int i=0; i < 36; i++) {
        spec_bb[i] = spectral_b[i] * bb;
      }
      
      //blend spectral primaries subtractive WGM
      float spectralmix[3][36] = {0};
      for (int i=0; i < 36; i++) {
        spectralmix[0][i] = powf(MAX(spec_ar[i], 0.0001), subfac) * powf(MAX(spec_br[i], 0.0001), (1-subfac));
        spectralmix[1][i] = powf(MAX(spec_ag[i], 0.0001), subfac) * powf(MAX(spec_bg[i], 0.0001), (1-subfac));
        spectralmix[2][i] = powf(MAX(spec_ab[i], 0.0001), subfac) * powf(MAX(spec_bb[i], 0.0001), (1-subfac));
      }
      //collapse into one SPD
      float spectral_combined[36] = {0};
      for (int i=0; i<36; i++) {
        spectral_combined[i] += spectralmix[0][i];
        spectral_combined[i] += spectralmix[1][i];
        spectral_combined[i] += spectralmix[2][i];
      }
      //convert to RGB
      for (int i=0; i<36; i++) {
        spectralmixsub[0] += T_MATRIX[0][i] * spectral_combined[i];
        spectralmixsub[1] += T_MATRIX[1][i] * spectral_combined[i];
        spectralmixsub[2] += T_MATRIX[2][i] * spectral_combined[i];
      }
      spectralmixsub[0] = CLAMP(spectralmixsub[0], 0.0f, 1.0f);
      spectralmixsub[1] = CLAMP(spectralmixsub[1], 0.0f, 1.0f);
      spectralmixsub[2] = CLAMP(spectralmixsub[2], 0.0f, 1.0f);
    }
    //mix rgb and spectral sub modes:
    for (int i=0; i < 3; i++) {  
      submix[i] = (((1-spectral)*rgbmixsub[i]) + (spectral*spectralmixsub[i]));
    }
  }



  //combine normal and subtractive RGB modes
  //combine in linear RGB if gamma was specified

  for (int i=0; i < 3; i++) {  
    result[i] = ((1-normsub)*normmix[i]) + (normsub*submix[i]);
  }
  //Now handle transform if necessary
  if (gamma != 1.0) {
    ar = result[0];
    ag = result[1];
    ab = result[2];

    //convert back to sRGB    
    rgb_to_srgb_float(&ar, &ag, &ab, gamma);

    result[0] = ar;
    result[1] = ag;
    result[2] = ab;
  }
  //alpha is simple
  result[0] = CLAMP(result[0], 0.0f, 1.0f);
  result[1] = CLAMP(result[1], 0.0f, 1.0f);
  result[2] = CLAMP(result[2], 0.0f, 1.0f);
  result[3] = CLAMP(fac * aa + (1-fac) * ba, 0.0f, 1.0f); 
  
  for (int i=0; i < 3; i++) {  
    if (isnan(result[i])) { result[i] = 0.0; }
  }
  
  //Chroma and Luminosity tweak
  //compare result to the smudge_state (a) and tweak the chroma and luma
  //based on hue angle difference in linear HCY
  if (smudge_desat != 0 || smudge_darken != 0) {
    float smudge_h = a[0];
    float smudge_c = a[1];
    float smudge_y = a[2];
    
    float result_h = result[0];
    float result_c = result[1];
    float result_y = result[2];
    
    //use specified gamma in case user wanted EOTF/linear
    srgb_to_rgb_float (&smudge_h, &smudge_c, &smudge_y, gamma);
    srgb_to_rgb_float (&result_h, &result_c, &result_y, gamma);
    rgb_to_hcy_float (&smudge_h, &smudge_c, &smudge_y);
    rgb_to_hcy_float (&result_h, &result_c, &result_y);

    //set our Brightness of the mix according to mode result. and also process saturation
    //Don't process achromatic colors (HCY has 3 achromatic states C=0 or Y=1 or Y=0)
    if (result_c != 0 && smudge_c != 0 && result_y != 0 && smudge_y != 0 && result_y != 1 && smudge_y != 1) {  

      //determine hue diff, proportional to smudge ratio
      //if fac is .5 the hueratio should be 1. When fac closer to 0 and closer to 1 should decrease towards zero
      //why- because if smudge is 0 or 1, only 100% of one of the brush or smudge color will be used so there is no comparison to make.
      //when fac is 0.5 the smudge and brush are mixed 50/50, so the huedifference should be respected 100%
      float huediff_sat;
      float huediff_bright;
      float hueratio = (0.5 - fabs(0.5 - fac)) / 0.5;
      float anglediff = fabs(smallest_angular_difference(result_h*360, smudge_h*360)/360);

      //calculate the adjusted hue difference and apply that to the saturation level and/or brightness
      //if smudge_desaturation setting is zero, the huediff will be zero.  Likewise when smudge (fac) is 0 or 1, the huediff will be zero.
      huediff_sat =  anglediff * smudge_desat * hueratio;
      //do the same for brightness
      huediff_bright = anglediff * smudge_darken * hueratio;

      //do the desaturation.  More strongly saturated colors will reduce S more than less saturated colors
      result_c = CLAMP(result_c*(1-huediff_sat), 0.0f, 1.0f);

      //attempt to simulate subtractive mode by darkening colors if they are different hues
      result_y = CLAMP(result_y*(1-huediff_bright), 0.0f, 1.0f);

      //convert back to companded sRGB
      hcy_to_rgb_float (&result_h, &result_c, &result_y);
      rgb_to_srgb_float (&result_h, &result_c, &result_y, gamma);

      result[0] = CLAMP(result_h, 0.0f, 1.0f);
      result[1] = CLAMP(result_c, 0.0f, 1.0f);
      result[2] = CLAMP(result_y, 0.0f, 1.0f);
    }
  }
  return result; 
}

#endif //HELPERS_C
