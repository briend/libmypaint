#ifndef BRUSHMODES_H
#define BRUSHMODES_H

void draw_dab_pixels_BlendMode_Normal (uint16_t * mask,
                                       uint16_t * rgba,
                                       uint16_t color_r,
                                       uint16_t color_g,
                                       uint16_t color_b,
                                       uint16_t opacity);
void
draw_dab_pixels_BlendMode_Color (uint16_t *mask,
                                 uint16_t *rgba, // b=bottom, premult
                                 uint16_t color_r,  // }
                                 uint16_t color_g,  // }-- a=top, !premult
                                 uint16_t color_b,  // }
                                 uint16_t opacity);
void
draw_dab_pixels_BlendMode_Posterize (uint16_t *mask,
                                 uint16_t *rgba, // b=bottom, premult
                                 uint16_t posterize,
                                 uint16_t posterize_num);
void draw_dab_pixels_BlendMode_Normal_and_Eraser (uint16_t * mask,
                                                  uint16_t * rgba,
                                                  uint16_t color_r,
                                                  uint16_t color_g,
                                                  uint16_t color_b,
                                                  uint16_t color_a,
                                                  uint16_t opacity);
void draw_dab_pixels_BlendMode_LockAlpha (uint16_t * mask,
                                          uint16_t * rgba,
                                          uint16_t color_r,
                                          uint16_t color_g,
                                          uint16_t color_b,
                                          uint16_t opacity);
void get_color_pixels_accumulate (uint16_t * mask,
                                  uint16_t * rgba,
                                  float * sum_weight,
                                  float * sum_r,
                                  float * sum_g,
                                  float * sum_b,
                                  float * sum_a
                                  );

void get_spectral_color_pixels_accumulate (uint16_t * mask,
                                  uint16_t * rgba,
                                  float * sum_weight,
                                  float * sum_spectral,
                                  float * sum_a,
                                  float gamma
                                  );


#endif // BRUSHMODES_H
