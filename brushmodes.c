/* libmypaint - The MyPaint Brush Library
 * Copyright (C) 2007-2014 Martin Renold <martinxyz@gmx.ch> et. al
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

#include <config.h>

#include <stdint.h>
#include <assert.h>
#include <math.h>
#include "fastapprox/fastpow.h"

#include "helpers.h"
#include "brushmodes.h"

// parameters to those methods:
//
// rgba: A pointer to 16bit rgba data with non-premultiplied alpha.
//       The range of each components is limited from 0 to 2^15.
//
// mask: Contains the dab shape, that is, the intensity of the dab at
//       each pixel. Usually rendering is done for one tile at a
//       time. The mask is LRE encoded to jump quickly over regions
//       that are not affected by the dab.
//
// opacity: overall strength of the blending mode. Has the same
//          influence on the dab as the values inside the mask.


// We are manipulating pixels with premultiplied alpha.
// This is an "over" operation (opa = topAlpha).
// In the formula below, topColor is assumed to be premultiplied.
//
//               opa_a      <   opa_b      >
// resultAlpha = topAlpha + (1.0 - topAlpha) * bottomAlpha
// resultColor = topColor + (1.0 - topAlpha) * bottomColor
//
void draw_dab_pixels_BlendMode_Normal (float *mask,
                                       float *rgba_buffer, DabBounds *b,
                                       float *brushcolor,
                                       float opacity,
                                       float volume) {

    for (int yp = b->y0; yp <= b->y1; yp++) {
      for (int xp = b->x0; xp <= b->x1; xp++) {
          const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
          const int a_chan = MYPAINT_NUM_CHANS-1;
          float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
          float opa_a = mask[offset]*opacity; // topAlpha
          float opa_b = 1.0-opa_a; // bottomAlpha
          for (int i=0; i < a_chan-1; i++) {
            rgba[i] = opa_a * brushcolor[i] + opa_b*rgba[i];
          }
          //rgba[a_chan-1] = opa_a * volume + rgba[a_chan-1]; //volume
          rgba[a_chan] = opa_a + opa_b * rgba[a_chan]; // alpha
          //printf("brushmode norm %f %f \n", rgba[a_chan-1], volume);

      }
    }
}

void draw_dab_pixels_BlendMode_Normal_Paint (float *mask,
                                       float *rgba_buffer, DabBounds *b,
                                       float *brushcolor,
                                       float opacity,
                                       float volume) {

    for (int yp = b->y0; yp <= b->y1; yp++) {
      for (int xp = b->x0; xp <= b->x1; xp++) {
          const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
          const int a_chan = MYPAINT_NUM_CHANS-1;
          float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
          float opa_a = mask[offset]*opacity; // topAlpha
          float opa_b = 1.0-opa_a; // bottomAlpha
          for (int i=0; i < a_chan-1; i++) {
            rgba[i] = opa_a * brushcolor[i] + opa_b*rgba[i];
          }
          //rgba[a_chan-1] = MAX((mask[offset] * volume) + brushcolor[a_chan-1] + rgba[a_chan-1], 0.0);
          rgba[a_chan] = opa_a + opa_b * rgba[a_chan];
          rgba[a_chan-1] = MAX(opa_a * brushcolor[a_chan-1], rgba[a_chan-1]); // use the larger vol, don't grow or shrink
          //rgba[a_chan-1] = opa_a * volume + rgba[a_chan-1];
          //printf("paint brushmode norm %f %f \n", rgba[a_chan-1], volume);

      }
    }
}

//Posterize.  Basically exactly like GIMP's posterize
//reduces colors by adjustable amount (posterize_num).
//posterize the canvas, then blend that via opacity
//does not affect alpha

void draw_dab_pixels_BlendMode_Posterize (float *mask,
                                          float *rgba_buffer, 
                                          DabBounds *b,
                                          float opacity,
                                          uint16_t posterize,
                                          uint16_t posterize_num) {

    for (int yp = b->y0; yp <= b->y1; yp++) {
      for (int xp = b->x0; xp <= b->x1; xp++) {
          const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
          const int a_chan = MYPAINT_NUM_CHANS-1;
          float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
          float opa_a = mask[offset]*opacity; // topAlpha
          float opa_b = 1.0-opa_a; // bottomAlpha
          for (int i=0; i < a_chan-1; i++) {
            float c = (float)rgba[i] * (1<<15);
            uint32_t post_c = (1<<15) * ROUND(c * posterize_num) / posterize_num;
            rgba[i] = ((float)opa_a*post_c + (float)opa_b*rgba[i])/(1<<30);
            }
          rgba[a_chan] = opa_a + opa_b * rgba[a_chan];
      }
  }
}

// Colorize: apply the source hue and saturation, retaining the target
// brightness. Same thing as in the PDF spec addendum, and upcoming SVG
// compositing drafts. Colorize should be used at either 1.0 or 0.0, values in
// between probably aren't very useful. This blend mode retains the target
// alpha, and any pure whites and blacks in the target layer.

#define MAX3(a, b, c) ((a)>(b)?MAX((a),(c)):MAX((b),(c)))
#define MIN3(a, b, c) ((a)<(b)?MIN((a),(c)):MIN((b),(c)))

// For consistency, these are the values used by MyPaint's Color and
// Luminosity layer blend modes, which in turn are defined by
// http://dvcs.w3.org/hg/FXTF/rawfile/tip/compositing/index.html.
// Same as ITU Rec. BT.601 (SDTV) rounded to 2 decimal places.

static const float LUMA_RED_COEFF   = 0.3;
static const float LUMA_GREEN_COEFF = 0.59;
static const float LUMA_BLUE_COEFF  = 0.11;


// See also http://en.wikipedia.org/wiki/YCbCr


/* Returns the sRGB luminance of an RGB triple, expressed as scaled ints. */

#define LUMA(r,g,b) \
   ((r)*LUMA_RED_COEFF + (g)*LUMA_GREEN_COEFF + (b)*LUMA_BLUE_COEFF)


/*
 * Sets the output RGB triple's luminance to that of the input, retaining its
 * colour. Inputs and outputs are scaled ints having factor 2**-15, and must
 * not store premultiplied alpha.
 */

inline static void
set_rgb16_lum_from_rgb16(const float topr,
                         const float topg,
                         const float topb,
                         float *botr,
                         float *botg,
                         float *botb)
{
    // Spec: SetLum()
    // Colours potentially can go out of band to both sides, hence the
    // temporary representation inflation.
    const float botlum = LUMA(*botr, *botg, *botb);
    const float toplum = LUMA(topr, topg, topb);
    const float diff = botlum - toplum;
    float r = topr + diff;
    float g = topg + diff;
    float b = topb + diff;

    // Spec: ClipColor()
    // Clip out of band values
    float lum = LUMA(r, g, b);
    float cmin = MIN3(r, g, b);
    float cmax = MAX3(r, g, b);
    if (cmin < 0) {
        r = lum + (((r - lum) * lum) / (lum - cmin));
        g = lum + (((g - lum) * lum) / (lum - cmin));
        b = lum + (((b - lum) * lum) / (lum - cmin));
    }
/*    if (cmax > (1<<15)) {*/
/*        r = lum + (((r - lum) * ((1<<15)-lum)) / (cmax - lum));*/
/*        g = lum + (((g - lum) * ((1<<15)-lum)) / (cmax - lum));*/
/*        b = lum + (((b - lum) * ((1<<15)-lum)) / (cmax - lum));*/
/*    }*/
#ifdef HEAVY_DEBUG
    assert(0 <= r);
    assert(0 <= g);
    assert(0 <= b);
#endif

    *botr = r;
    *botg = g;
    *botb = b;
}


// The method is an implementation of that described in the official Adobe "PDF
// Blend Modes: Addendum" document, dated January 23, 2006; specifically it's
// the "Color" nonseparable blend mode. We do however use different
// coefficients for the Luma value.

void
draw_dab_pixels_BlendMode_Color (float * mask,
                                 float * rgba_buffer, // b=bottom, premult
                                 DabBounds *b,
                                 float color_r,  // }
                                 float color_g,  // }-- a=top, !premult
                                 float color_b,  // }
                                 float opacity)
{

  for (int yp = b->y0; yp <= b->y1; yp++) {
    for (int xp = b->x0; xp <= b->x1; xp++) {
      const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
      float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);

    // De-premult
    float r, g, b;
    const float a = rgba[3];
    r = g = b = 0;
    if (rgba[3] != 0) {
      r = rgba[0] / a;
      g = rgba[1] / a;
      b = rgba[2] / a;
    }

    // Apply luminance
    set_rgb16_lum_from_rgb16(color_r, color_g, color_b, &r, &g, &b);

    // Re-premult
    r = r * a;
    g = g * a;
    b = b * a;

    // And combine as normal.
    float opa_a = mask[offset] * opacity; // topAlpha
    float opa_b = 1.0 - opa_a; // bottomAlpha
    rgba[0] = (opa_a*r + opa_b*rgba[0]);
    rgba[1] = (opa_a*g + opa_b*rgba[1]);
    rgba[2] = (opa_a*b + opa_b*rgba[2]);
    }
  }
}

// This blend mode is used for smudging and erasing.  Smudging
// allows to "drag" around transparency as if it was a color.  When
// smuding over a region that is 60% opaque the result will stay 60%
// opaque (color_a=0.6).  For normal erasing color_a is set to 0.0
// and color_r/g/b will be ignored. This function can also do normal
// blending (color_a=1.0).
//
void draw_dab_pixels_BlendMode_Normal_and_Eraser (float *mask,
                                                  float *rgba_buffer,
                                                  DabBounds *b,
                                                  float *brushcolor,
                                                  float color_a,
                                                  float opacity,
                                                  float volume) {

  //printf("brusmodes 227 colors are %f, %f, %f, %f, %f\n", color_r, color_g, color_b, opacity, color_a);
  for (int yp = b->y0; yp <= b->y1; yp++) {
    for (int xp = b->x0; xp <= b->x1; xp++) {
      const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
      float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
      const int a_chan = MYPAINT_NUM_CHANS-1;

      float opa_a = mask[offset]*opacity; // topAlpha
      float opa_b = 1.0-opa_a; // bottomAlpha
      opa_a = opa_a * color_a;
      rgba[a_chan] = opa_a + opa_b * rgba[a_chan];
      for (int i=0; i<a_chan-1; i++){
        rgba[i] = (opa_a*brushcolor[i] + opa_b*rgba[i]);
      }
    }
  }
}

void draw_dab_pixels_BlendMode_Normal_and_Eraser_Paint (float *mask,
                                                  float *rgba_buffer,
                                                  DabBounds *b,
                                                  float *brushcolor,
                                                  float color_a,
                                                  float opacity,
                                                  float volume) {

  for (int yp = b->y0; yp <= b->y1; yp++) {
    for (int xp = b->x0; xp <= b->x1; xp++) {
      const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
      float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
      const int a_chan = MYPAINT_NUM_CHANS-1;

      float opa_a = mask[offset]*opacity; // topAlpha
      float opa_b = 1.0-opa_a; // bottomAlpha
      opa_a = opa_a * color_a;
      rgba[a_chan] = opa_a + opa_b * rgba[a_chan];
      //rgba[a_chan-1] = opa_a * volume + rgba[a_chan-1];
      
      for (int i=0; i<a_chan; i++){
        rgba[i] = (opa_a*brushcolor[i] + opa_b*rgba[i]);
      }
      //printf("bm bv is %f, canvas is %f \n", brushcolor[a_chan-1], rgba[a_chan-1]);
      //rgba[a_chan-1] = opa_a * volume + rgba[a_chan-1];
      
    }
  }
}

// This is BlendMode_Normal with locked alpha channel.
//
void draw_dab_pixels_BlendMode_LockAlpha (float * mask,
                                          float * rgba_buffer,
                                          DabBounds *b,
                                          float *brushcolor,
                                          float opacity,
                                          float volume) {

  for (int yp = b->y0; yp <= b->y1; yp++) {
    for (int xp = b->x0; xp <= b->x1; xp++) {
        const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
        float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
        const int a_chan = MYPAINT_NUM_CHANS-1;

        float opa_a = mask[offset]*opacity; // topAlpha
        float opa_b = 1.0-opa_a; // bottomAlpha
        
        opa_a *= rgba[a_chan];
        for (int i=0; i<a_chan-1; i++){    
          rgba[i] = (opa_a*brushcolor[i] + opa_b*rgba[i]);
        }

    }
  }
}

void draw_dab_pixels_BlendMode_LockAlpha_Paint (float * mask,
                                          float * rgba_buffer,
                                          DabBounds *b,
                                          float *brushcolor,
                                          float opacity,
                                          float volume) {

  for (int yp = b->y0; yp <= b->y1; yp++) {
    for (int xp = b->x0; xp <= b->x1; xp++) {
        const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
        float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);
        const int a_chan = MYPAINT_NUM_CHANS-1;

        float opa_a = mask[offset]*opacity; // topAlpha
        float opa_b = 1.0-opa_a; // bottomAlpha
        
        opa_a *= rgba[a_chan];
        for (int i=0; i<a_chan; i++){
          rgba[i] = (opa_a*brushcolor[i] + opa_b*rgba[i]);
        }

    }
  }
}


// Sum up the color/alpha components inside the masked region.
// Called by get_color().
//
void get_color_pixels_accumulate (float *mask,
                                  float *rgba_buffer, DabBounds *bb,
                                  float * sum_weight,
                                  float * sum_color,
                                  float paint
                                  ) {

  for (int yp = bb->y0; yp <= bb->y1; yp++) {
    for (int xp = bb->x0; xp <= bb->x1; xp++) {
        const int offset = (yp*MYPAINT_TILE_SIZE)+xp;
        float *rgba = rgba_buffer + (offset*MYPAINT_NUM_CHANS);

      float opa = mask[offset];
      *sum_weight += opa;
      for (int i=0; i<MYPAINT_NUM_CHANS; i++){
        sum_color[i] += opa*rgba[i];
      }
      //printf("get vol is %f \n", rgba[MYPAINT_NUM_CHANS-2]);
    }
  }
}




