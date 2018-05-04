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
float * mix_colors(float *a, float *b, float fac, float gamma, float normsub, gboolean get, float smudge_darken, float smudge_desat)
{
  float rgbmixnorm[4] = {0};
  float rgbmixsub[4] = {0};
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

    rgbmixnorm[0] = fac * ar + (1-fac) * bra;
    rgbmixnorm[1] = fac * ag + (1-fac) * bga;
    rgbmixnorm[2] = fac * ab + (1-fac) * bba;
  }

  //subtractive
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
      //mix with Weighted Geometric Mean
      //don't allow absolute zero because it has infinite
      //darkening power
      rgbmixsub[0] = powf(MAX(ar, 0.0001), subfac) * powf(MAX(br, 0.0001),(1-subfac));
      rgbmixsub[1] = powf(MAX(ag, 0.0001), subfac) * powf(MAX(bg, 0.0001),(1-subfac));
      rgbmixsub[2] = powf(MAX(ab, 0.0001), subfac) * powf(MAX(bb, 0.0001),(1-subfac));
  }

  //do eraser_target alpha for smudge+brush color
  if (get == FALSE) {
    rgbmixnorm[0] /=ba;
    rgbmixnorm[1] /=ba;
    rgbmixnorm[2] /=ba;
  }    

  //combine normal and subtractive RGB modes
  //combine in linear RGB if gamma was specified

  for (int i=0; i < 3; i++) {  
    result[i] = ((1-normsub)*rgbmixnorm[i]) + (normsub*rgbmixsub[i]);
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