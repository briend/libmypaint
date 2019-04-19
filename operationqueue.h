#ifndef OPERATIONQUEUE_H
#define OPERATIONQUEUE_H

#include <stdint.h>
#include "tilemap.h"
#include "mypaint-config.h"

typedef struct {
    float x;
    float y;
    float radius;
    float color_r;
    float color_g;
    float color_b;
    float color_a;
    float brushcolor[MYPAINT_NUM_CHANS-1];
    float opaque;
    float hardness;
    float aspect_ratio;
    float angle;
    float normal;
    float lock_alpha;
    float colorize;
    float posterize;
    float posterize_num;
    float paint;
} OperationDataDrawDab;

typedef struct OperationQueue OperationQueue;

OperationQueue *operation_queue_new(void);
void operation_queue_free(OperationQueue *self);

int operation_queue_get_dirty_tiles(OperationQueue *self, TileIndex** tiles_out);
void operation_queue_clear_dirty_tiles(OperationQueue *self);

void operation_queue_add(OperationQueue *self, TileIndex index, OperationDataDrawDab *op);
OperationDataDrawDab *operation_queue_pop(OperationQueue *self, TileIndex index);

OperationDataDrawDab *operation_queue_peek_first(OperationQueue *self, TileIndex index);
OperationDataDrawDab *operation_queue_peek_last(OperationQueue *self, TileIndex index);

#endif // OPERATIONQUEUE_H
