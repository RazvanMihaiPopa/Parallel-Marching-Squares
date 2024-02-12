// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include "pthread_barrier_mac.h"

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

typedef struct {
    pthread_barrier_t *barrier;
    ppm_image *image, *new_image;
    ppm_image **contour_map;
    int step_x, step_y, thread_id, number_threads; 
    unsigned char sigma;
    unsigned char **grid;
} function_arguments;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

void *f(void *arg)
{
    function_arguments *args = (function_arguments *)arg;

    if (args->image->x > RESCALE_X || args->image->y > RESCALE_Y) {
        uint8_t sample[3];

        int start1 = args->thread_id * (double)args->new_image->x / args->number_threads;
        int end1 = (args->thread_id + 1) * (double)args->new_image->x / args->number_threads;
        if (end1 > args->new_image->x) {
            end1 = args->new_image->x;
        }

        // use bicubic interpolation for scaling
        for (int i = start1 ; i < end1; i++) {
            for (int j = 0; j < args->new_image->y; j++) {
                float u = (float)i / (float)(args->new_image->x - 1);
                float v = (float)j / (float)(args->new_image->y - 1);
                sample_bicubic(args->image, u, v, sample);

                args->new_image->data[i * args->new_image->y + j].red = sample[0];
                args->new_image->data[i * args->new_image->y + j].green = sample[1];
                args->new_image->data[i * args->new_image->y + j].blue = sample[2];
            }
        }
        pthread_barrier_wait(args->barrier);

        if (args->thread_id == 0) {
            free(args->image->data);
            free(args->image);
        }

        pthread_barrier_wait(args->barrier);
        args->image = args->new_image;
    }

    int p = args->image->x / args->step_x;
    int q = args->image->y / args->step_y;
    
    int start = args->thread_id * (double)p / args->number_threads;
    int end = (args->thread_id + 1) * (double)p / args->number_threads;
    if (end > p) {
        end = p;
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = args->image->data[i * args->step_x * args->image->y + j * args->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > args->sigma) {
                args->grid[i][j] = 0;
            } else {
                args->grid[i][j] = 1;
            }
        }
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = args->image->data[i * args->step_x * args->image->y + args->image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > args->sigma) {
            args->grid[i][q] = 0;
        } else {
            args->grid[i][q] = 1;
        }
    }

    if (args->thread_id == 0) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = args->image->data[(args->image->x - 1) * args->image->y + j * args->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > args->sigma) {
                args->grid[p][j] = 0;
            } else {
                args->grid[p][j] = 1;
            }
        }
    }

    pthread_barrier_wait(args->barrier);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * args->grid[i][j] + 4 * args->grid[i][j + 1] + 2 * args->grid[i + 1][j + 1] + 1 * args->grid[i + 1][j];
            update_image(args->image, args->contour_map[k], i * args->step_x, j * args->step_y);
        }
    }

	pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    int num_threads = atoi(argv[3]);
    int i, r;
	void *status;
	pthread_t threads[num_threads];
	function_arguments arguments[num_threads];
    
    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;

    // Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // Allocate memory for rescaled image
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // Allocate the grid
    unsigned char **grid;
    grid = (unsigned char **)malloc((RESCALE_X/step_x + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= RESCALE_X/step_x; i++) {
        grid[i] = (unsigned char *)malloc((RESCALE_Y/step_y + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    pthread_barrier_t barrier;
    r = pthread_barrier_init(&barrier, NULL, num_threads);
    if (r) {
        exit(-1);
    }

    for (i = 0; i < num_threads; i++) {
        arguments[i].grid = grid;
        arguments[i].thread_id = i;
        arguments[i].new_image = new_image;
        arguments[i].barrier = &barrier;
        arguments[i].number_threads = num_threads;
        arguments[i].step_x = step_x;
        arguments[i].step_y = step_y;
        arguments[i].sigma = SIGMA;
        arguments[i].image = image;
        arguments[i].contour_map = contour_map;

		r = pthread_create(&threads[i], NULL, f, &arguments[i]);

		if (r) {
			printf("Eroare la crearea thread-ului %d\n", i);
			exit(-1);
		}
	}

	for (i = 0; i < num_threads; i++) {
		r = pthread_join(threads[i], &status);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}

    pthread_barrier_destroy(&barrier);

    ppm_image *scaled_image = arguments[0].image;

    // Write output
    write_ppm(scaled_image, argv[2]);

    free_resources(scaled_image, contour_map, grid, step_x);

    return 0;
}
