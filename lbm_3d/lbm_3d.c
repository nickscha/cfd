#include "../cfd.h" /* Computational Fluid Dynamics API     */
#include "csr.h"    /* Software Renderer                    */
#include "vm.h"     /* Linear Algebra / Vetcor Math Library */
#include "../tests/cfd_renderer.h"

#include <stdlib.h>
#include <stdio.h>

static float cube_vertices[] = {
    -0.5f, -0.5f, 0.5f,
    0.5f, -0.5f, 0.5f,
    0.5f, 0.5f, 0.5f,
    -0.5f, 0.5f, 0.5f,
    -0.5f, -0.5f, -0.5f,
    0.5f, -0.5f, -0.5f,
    0.5f, 0.5f, -0.5f,
    -0.5f, 0.5f, -0.5f};

static int cube_indices[] = {
    0, 3, 2, 0, 2, 1, /* Front face (+z normal, facing camera)   */
    4, 5, 6, 4, 6, 7, /* Back face (-z normal, away from camera) */
    3, 7, 6, 3, 6, 2, /* Top face (+y normal)                    */
    0, 1, 5, 0, 5, 4, /* Bottom face (-y normal)                 */
    1, 2, 6, 1, 6, 5, /* Right face (+x normal)                  */
    0, 4, 7, 0, 7, 3  /* Left face (-x normal)                   */
};

static unsigned long cube_vertices_size = sizeof(cube_vertices) / sizeof(cube_vertices[0]);
static unsigned long cube_indices_size = sizeof(cube_indices) / sizeof(cube_indices[0]);

static void csr_save_ppm(char *filename_format, int frame, csr_context *context)
{
    FILE *fp;
    char filename[64];

    /* Format the filename with the frame number */
    sprintf(filename, filename_format, frame);

    fp = fopen(filename, "wb");

    if (!fp)
    {
        fprintf(stderr, "Error: Could not open file %s for writing.\n", filename);
        return;
    }

    /* PPM header */
    fprintf(fp, "P6\n%d %d\n255\n", context->width, context->height);

    /* Pixel data */
    fwrite(context->framebuffer, sizeof(csr_color), (size_t)(context->width * context->height), fp);

    fclose(fp);
}

static csr_context context = {0};
static csr_color clear_color = {40, 40, 40};
static m4x4 projection_view;

int init_graphics(int window_width, int window_height, v3 cam_position, v3 cam_look_at_pos)
{
    void *memory_graphics = malloc(csr_memory_size(window_width, window_height));

    if (!csr_init_model(&context, memory_graphics, csr_memory_size(window_width, window_height), window_width, window_height))
    {
        return 0;
    }

    {
        v3 world_up = vm_v3(0.0f, 1.0f, 0.0f);
        float cam_fov = 40.0f;

        m4x4 projection = vm_m4x4_perspective(vm_radf(cam_fov), (float)context.width / (float)context.height, 0.1f, 1000.0f);
        m4x4 view = vm_m4x4_lookAt(cam_position, cam_look_at_pos, world_up);
        projection_view = vm_m4x4_mul(projection, view);
    }

    return 1;
}

int main(void)
{
    /* Setup the lattice boltzmann parameters */
    int xdim = 50;
    int ydim = 50;
    int zdim = 50;

    float viscosity = 0.02f;
    float omega = 1.0f / (3.0f * viscosity + 0.5f);
    float speed = 0.1f;

    int frame_count = 250;
    int frame = 0;
    int frame_steps = 20;

    unsigned long memory_grid_size = cfd_lbm_3d_grid_memory_size(xdim, ydim, zdim);
    unsigned long memory_size = memory_grid_size;
    void *memory = malloc(memory_size);

    cfd_lbm_3d_grid grid = {0};

    {
        printf("[cfd][lbm][3d] mem. grid (MB): %10.4f\n", (double)memory_grid_size / 1024.0 / 1024.0);
        printf("[cfd][lbm][3d]      viscosity: %10.4f\n", (double)viscosity);
        printf("[cfd][lbm][3d]          omega: %10.4f\n", (double)omega);
        printf("[cfd][lbm][3d]          speed: %10.4f\n", (double)speed);
        printf("[cfd][lbm][3d]              x: %10i\n", xdim);
        printf("[cfd][lbm][3d]              y: %10i\n", ydim);
        printf("[cfd][lbm][3d]              z: %10i\n", zdim);
    }

    init_graphics(
        800, 600,                                                                                             /* Screen dimensions      */
        vm_v3(-((float)xdim * 0.5f), (float)ydim + ((float)ydim * 0.5f), (float)zdim + ((float)zdim * 0.5f)), /*Camera Position         */
        vm_v3(((float)xdim * 0.5f), ((float)ydim * 0.5f), ((float)zdim * 0.5f))                               /* Camera LookAt Position */
    );

    cfd_build_colormap();

    {
        /* Init lattice boltzmann */
        cfd_lbm_3d_init_grid(&grid, memory, xdim, ydim, zdim, omega);
        cfd_lbm_3d_init_fluid(&grid, speed);
        cfd_lbm_3d_init_barriers(&grid);
        cfd_lbm_3d_init_tracers(&grid);

        for (frame = 0; frame < frame_count; ++frame)
        {
            int step;
            int x, y, z, t;

            /* Run fluid simulation steps */
            for (step = 0; step < frame_steps; ++step)
            {
                cfd_lbm_3d_collide_and_stream(&grid);
                cfd_lbm_3d_move_tracers(&grid);
            }

            /* Clear Screen Frame and Depth Buffer */
            csr_render_clear_screen(&context, clear_color);

            /* Render Grid and Barriers */
            for (z = 0; z < grid.zdim; ++z)
            {
                for (y = 0; y < grid.ydim; ++y)
                {
                    for (x = 0; x < grid.xdim; ++x)
                    {
                        int xdim_ydim = grid.xdim * grid.ydim;
                        int idx = x + y * grid.xdim + z * xdim_ydim;

                        m4x4 model = vm_m4x4_translate(vm_m4x4_identity, vm_v3((float)x, (float)y, (float)z));

                        if (grid.barrier[idx])
                        {
                            m4x4 model_view_projection = vm_m4x4_mul(projection_view, model);

                            csr_color c0 = {200, 200, 200};

                            csr_render(
                                &context,
                                CSR_RENDER_SOLID,
                                CSR_CULLING_CCW_BACKFACE, 3,
                                cube_vertices, cube_vertices_size,
                                cube_indices, cube_indices_size,
                                model_view_projection.e, c0, c0, c0);
                        }
                        else
                        {
                            /*
                            m4x4 model_view_projection = vm_m4x4_mul(projection_view, vm_m4x4_scalef(model, 0.05f));
                            */
                            m4x4 model_view_projection = vm_m4x4_mul(
                                projection_view,
                                vm_m4x4_scale(model, vm_v3(0.05f, 0.05f, 0.05f)));

                            csr_color c0 = {120, 120, 120};

                            float value = cfd_lbm_3d_calculate_curl(&grid, x, y, z);

                            float contrast = 0.0f;
                            float contrastFactor = cfd_powf(1.2f, contrast);
                            int cIndex = (int)(400 * (value * contrastFactor + 0.5f));
                            cIndex = cIndex < 0 ? 0 : cIndex;
                            cIndex = cIndex > 400 ? 400 : cIndex;

                            c0.r = cfd_color_map[cIndex].r;
                            c0.g = cfd_color_map[cIndex].g;
                            c0.b = cfd_color_map[cIndex].b;

                            csr_render(
                                &context,
                                CSR_RENDER_SOLID,
                                CSR_CULLING_CCW_BACKFACE, 3,
                                cube_vertices, cube_vertices_size,
                                cube_indices, cube_indices_size,
                                model_view_projection.e, c0, c0, c0);
                        }
                    }
                }
            }

            /* Render Tracers */
            for (t = 0; t < CFD_LBM_3D_NUMBER_TRACERS; ++t)
            {
                csr_color tracer_color = {0, 255, 0};
                v3 tracer_position = vm_v3(grid.tracerX[t], grid.tracerY[t], grid.tracerZ[t]);

                if (tracer_position.x > xdim || tracer_position.y > ydim || tracer_position.z > zdim ||
                    tracer_position.x <= 0 || tracer_position.y <= 0 || tracer_position.z <= 0)
                {
                    continue;
                }

                {
                    m4x4 model = vm_m4x4_translate(vm_m4x4_identity, tracer_position);

                    m4x4 model_view_projection = vm_m4x4_mul(
                        projection_view,
                        vm_m4x4_scale(model, vm_v3(0.25f, 0.25f, 0.25f)));

                    csr_render(
                        &context,
                        CSR_RENDER_SOLID,
                        CSR_CULLING_CCW_BACKFACE, 3,
                        cube_vertices, cube_vertices_size,
                        cube_indices, cube_indices_size,
                        model_view_projection.e, tracer_color, tracer_color, tracer_color);
                }
            }

            csr_save_ppm("lbm_%05d.ppm", frame, &context);
        }
    }

    free(memory);

    return 0;
}