/**
 * MODELING SNOW CRYSTAL GROWTH II (normal/slow version)
 *
 *  Janko Gravner <gravner@math.ucdavis.edu>
 *  David Griffeath <griffeat@math.wisc.edu>
 *
 * Refactoring by Chengyu HAN, 2022/10/22
 * Revised version, September 2007
 */
#include <math.h>
#include <stdio.h>
#include <string.h> // strlen
#include <stdbool.h> // true; false
#include <stdlib.h> // srand48, drand48
#include <time.h>
#include <limits.h> // UCHAR_MAX, USHRT_MAX

// for dev
#include <assert.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>


// Use a faster, not more meaningful implementation
#define _USE_FAST_IMPL


#define NR_MAX 1002
#define NC_MAX 1002

#define KAPPA_MAX 64

// The cell is part of the snowflake.
#define HAS_SNOWFLAKE   true
// Cell is NOT yet part of the snowflake.
#define NO_SNOWFLAKE   false


/* ==== Input Parameters ==== */

/* --- [Initial state] 初始化参数 */

/** [initialize] 初始气相密度.
 * `init_gas_rho` the density of diffusing particles elsewhere
 */
double init_gas_rho;
/** [initialize] 初始晶种半径.
 * `h` is the radius of the initial hexagon of density `p` in the middle of the array
 * 
 * Note. h<0 is interpreted as special initialization with radius -h.
 */
int init_crystal_seed_radius;
int twelve_sided;
/** [initialize] 初始晶种结晶概率 [0, 1].
 * crystallization probability `p`
 */
double init_crystal_seed_probability;

/* --- [Dynamics] 结晶动力学参数 */

/** [freezing] 气相(直接)结晶率.
 * 
 * + kappa       凝华进入【固相】
 * + (1 - kappa) 液化进入【液相】
 */
double kappa;
/** [attachment] 1~2晶体邻居，【液相】结晶阈值下限.
 * 
 * 拥有 1~2 个已结晶的邻居时, 元胞的结晶阈值.
 * ( >= beta ) 即可结晶
 */
double beta;
/** [attachment] 3晶体邻居，【液相】结晶阈值下限.
 * 
 */
double alpha;
/** [attachment] 3晶体邻居，邻居【气相】结晶阈值上限.
 * 
 */
double theta;
/** [melting] 液相气化率 */
double mu;
/** [melting] 固相升华率 */
double gam;
/** [noise] 添加的噪声幅度 */
double sigma;


/** Size of the (square) LxL system.
 * 
 * nr = n_row;
 * nc = n_col;
 */
int nr, nc;
/** size of the pixel */
int sp;

/** Input/output files. */
#define MAX_IO_PATH_LEN  256
/** simulation status file */
FILE *g_state_file;
char g_in_file_path[MAX_IO_PATH_LEN];
char g_out_file_path[MAX_IO_PATH_LEN];
char g_graphics_file_path[MAX_IO_PATH_LEN];
char g_grahics_viewer_name[MAX_IO_PATH_LEN];
char g_comments[100];


/* ==== Global Variables ==== */

// ---- initialize()
/** Iteration Counter */
int g_pq;
/** Iteration stop flag.
 * Stop iterating when the snowflake grows to 2/3 of the grid.
 * 
 * Note: search for "g_stop = true" to check the condition.
 */
bool g_stop;
/** set and update, but [not use]. */
int g_par_update;
/** Counters, starting from 1, for updating `ash[]`. */
int g_par_ash;

/** Center point of the grid: (nr/2, nc/2).
 * 
 * Note: After initialization, is a Constant.
 */
int g_center_i, g_center_j;
/** The frontier of snowflake growing. */
int g_r_old, g_r_new;

/** diffusive mass (vapor).
 * 
 */
double  d_dif[NR_MAX][NC_MAX];
/** attachment flag, 
 * indicator of snowflake sites.
 * 
 * `true` if x belongs to the crystal at time t.
 */
bool    a_pic[NR_MAX][NC_MAX];
/** boundary mass (quasi-liquid).
 * 
 */
double  b__fr[NR_MAX][NC_MAX];
/** crystal mass (ice).
 * 
 */
double  c__lm[NR_MAX][NC_MAX];

/** rings pallette
 * 
 */
int     ash[NR_MAX][NC_MAX];

// ---- other global var
int g_noac;

// is `b__fr` changed?
bool g_is_fr_changed;


/* ==== X11 Window ==== */
Display *g_xDisplay;
// main window
Window g_xWindow;
GC g_xGC;
XEvent g_xEvent;
// windows size hint
XSizeHints g_xSizeHints;
// default screen number
int g_xScreen;
// white and black pixel value
unsigned long g_xBlack, g_xWhite;
// main while loop control flag.
bool g_exit_flag;

// x11 windows in args
const char gui_ICON_NAME_STR[] = "sn";
const char gui_WINDOW_NAME_STR[] = "digital snowflake";
const char gui_TIME_STR[] = "time:";
const char gui_ACTIVE_STR[] = "active area:";

// ---- color map
Colormap g_cmap;
XColor g_color[KAPPA_MAX];
XColor g_color_on[128];
XColor g_color_off[128];
XColor g_othp[20];
// temp color array
int g_red[125], g_green[125], g_blue[125];


/* ==== Helper function forward declaration ==== */
// --- colors
int gui_get_off_color_idx(int i, int j);
int gui_get_on_color_idx(int i, int j);
int gui_get_color_idx(int i, int j);
int gui_get_othp_color_idx(int i, int j);

/* ==== Speed helper functions ==== */
bool not_snowflake(bool pic) {
#ifdef _USE_FAST_IMPL
    return !pic;
#else
    return pic == NO_SNOWFLAKE;
#endif
} // not_snowflake

bool is_snowflake(bool pic) {
#ifdef _USE_FAST_IMPL
    return pic;
#else
    return pic == HAS_SNOWFLAKE;
#endif
} // is_snowflake

/* ==== Math helper functions ==== */

/**
 * nonnegative doubleprecision floating-point values uniformly distributed
 *  over the interval [0.0, 1.0).
 */
double uniform_01rand()
{
    return drand48();
} /* uniform_01rand() */

/** max( abs(i), abs(j) ) */
int norm_inf(int i, int j)
{
#ifdef _USE_FAST_IMPL
    i = abs(i);
    j = abs(j);
#else
    if (i < 0)
        i = -i;
    if (j < 0)
        j = -j;
#endif // _USE_FAST_IMPL

    if (i > j)
        return i;
    else
        return j;
} /* norm_inf(int i, int j) */

/** 
 * 计算当前网格相对于中心点(g_center_i, g_center_j)的距离.
 * 
 * ref: 
 *  - "max" form, https://www.redblobgames.com/grids/hexagons/#distances
 *  - Axial coordinates, https://www.redblobgames.com/grids/hexagons/#distances-axial
 */
int hex_coord_distance(int i, int j)
{
    int s = -i - j;
    int center_k = -g_center_i - g_center_j;
    int dq = i - g_center_i;
    int dr = j - g_center_j;
    int ds = s - center_k;

    return (abs(dq) + abs(dr) + abs(ds)) / 2;
} /* hex_coord_distance(int i, int j) */

/** abs(i + j) */
int semi_norm(int i, int j)
{
#ifdef _USE_FAST_IMPL
    return abs(i + j);
#else
    int k;

    k = i + j;
    if (k >= 0)
        return k;
    else
        return -k;
#endif // _USE_FAST_IMPL
} /* semi_norm(int i, int j) */


/* ==== Main simulation algorithm func ==== */

void initialize()
{
    time_t t1, t2;
    int i, j, k;
    double x;
    double x1, y1;

    printf(".initialize: creating init. state\n");
    g_pq = 0;
    g_stop = false;
    g_par_update = 0;

    // --- gen rand seed
    t1 = time(&t2);
    t1 = t1 % 1000;
    srand48(t1);
    printf("seed:%ld\n", t1);

    for (i = 1; i < t1; i++)
    {
        x = drand48();
    }

    g_center_i = nr / 2;
    g_center_j = nc / 2;

    g_r_old = 0;
    g_r_new = 0;

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            x = uniform_01rand();

            if (twelve_sided == 0)
            {
                if (
                    (norm_inf(i - g_center_i, j - g_center_j) <= init_crystal_seed_radius) 
                    && (semi_norm(i - g_center_i, j - g_center_j) <= init_crystal_seed_radius) 
                    && (x <= init_crystal_seed_probability))
                {
                    // seed for the flake
                    d_dif[i][j] = 0.0;
                    a_pic[i][j] = HAS_SNOWFLAKE;
                    b__fr[i][j] = 1.0;
                    ash[i][j] = 0;
                    c__lm[i][j] = 0.0;
                    k = norm_inf(i - g_center_i, j - g_center_j);
                    if (k > g_r_new)
                        g_r_new = k;
                }
                else
                {
                    // filled with water vapour
                    d_dif[i][j] = init_gas_rho;
                    a_pic[i][j] = NO_SNOWFLAKE;
                    b__fr[i][j] = 0.0;
                    ash[i][j] = 0;
                    c__lm[i][j] = 0.0;
                }
            }
            else
            {
                x1 = (double)(i - g_center_i) / init_crystal_seed_radius;
                y1 = (double)(j - g_center_j) / init_crystal_seed_radius;
                if (
                    (((i - g_center_i) == -(j - g_center_j)) && (i - g_center_i <= 0) && (i - g_center_i >= -init_crystal_seed_radius)) 
                    || ((i - g_center_i >= 0) && (i - g_center_i <= init_crystal_seed_radius) && (j - g_center_j == 0)) 
                    || ((j - g_center_j <= 0) && (j - g_center_j >= -init_crystal_seed_radius) && (i - g_center_i == 0)))
                {
                    d_dif[i][j] = 0.0;
                    a_pic[i][j] = HAS_SNOWFLAKE;
                    b__fr[i][j] = 0.0;
                    ash[i][j] = 0;
                    c__lm[i][j] = 1.0;
                    k = norm_inf(i - g_center_i, j - g_center_j);
                    if (k > g_r_new)
                        g_r_new = k;
                }
                else
                {
                    d_dif[i][j] = init_gas_rho;
                    a_pic[i][j] = NO_SNOWFLAKE;
                    b__fr[i][j] = 0.0;
                    ash[i][j] = 0;
                    c__lm[i][j] = 0.0;
                }
            }
        }
    }

    g_r_old = g_r_new;
    g_par_ash = 1;

    printf(".initialize: init. finished\n");
} /* initialize() */

/**
 * 对元胞及其邻居的气象质量(d)进行平均化.
 * 
 */
void dynamics_diffusion()
{
    // 更新 a_pic 用的临时数组 
    double b[NR_MAX][NC_MAX] = { {0.0} };
    int i, j;
    int id, iu, jl, jr;
    // 当前格子近邻未结晶格子计数.
    int not_flake_count;

    // --- b[] 置零
    // TODO: memset(b, 0, sizeof(b));
    for (i = 0; i < nr; i++)
        for (j = 0; (j < nc); j++)
        {
            b[i][j] = 0.0;
        }

    // --- 遍历所有非结晶格子, 计算新的气象质量, 结果保存到临时变量 b[] 中
    for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
        {
            if (not_snowflake(a_pic[i][j]))
            {
                // 使用循环边界
                id = (i + 1) % nr;
                iu = (i + nr - 1) % nr;
                jr = (j + 1) % nr;
                jl = (j + nr - 1) % nr;

                not_flake_count = 0;
                if (not_snowflake(a_pic[id][j]))
                    not_flake_count++;
                if (not_snowflake(a_pic[iu][j]))
                    not_flake_count++;
                if (not_snowflake(a_pic[i][jl]))
                    not_flake_count++;
                if (not_snowflake(a_pic[i][jr]))
                    not_flake_count++;
                if (not_snowflake(a_pic[iu][jr]))
                    not_flake_count++;
                if (not_snowflake(a_pic[id][jl]))
                    not_flake_count++;

                // Note: 实际上 if 分支可以合并. not_flake_count == 0 时, else 分支与 if 结果一致.
                if (not_flake_count == 0)
                    b[i][j] = d_dif[i][j];
                else
                {
                    /** 气相扩散后，剩余的气相质量.
                     * 
                     * 每个非雪花的格子，将均分当前格子 1/7 的气相质量.
                     */
                    double residual_gas_mass = (1.0 - (double)not_flake_count / 7.0) * d_dif[i][j];
                    /** 从其他格子扩散来的总气相质量.
                     * 
                     * - 已经变成雪花的格子，则不会有气相质量扩散过来.
                     * - 未结晶的格子, 会扩散 1/7 的质量过来.
                     */
                    double received_gas_mass =
                        (
                              d_dif[id][j] * not_snowflake(a_pic[id][j]) 
                            + d_dif[iu][j] * not_snowflake(a_pic[iu][j])
                            + d_dif[i][jl] * not_snowflake(a_pic[i][jl]) 
                            + d_dif[i][jr] * not_snowflake(a_pic[i][jr])
                            + d_dif[iu][jr] * not_snowflake(a_pic[iu][jr])
                            + d_dif[id][jl] * not_snowflake(a_pic[id][jl])
                        ) / 7.0;
                        
                    b[i][j] = residual_gas_mass + received_gas_mass;
                }
            } /* if (not_snowflake(a_pic[i][j])) */
        }

    // --- 更新非结晶格子的气相质量 d_dif
    for (i = 0; i < nr; i++)
        for (j = 0; (j < nc); j++)
        {
            if (not_snowflake(a_pic[i][j]))
                d_dif[i][j] = b[i][j];
        }
} /* dynamics_diffusion() */

/**
 * 处理气相液化与结晶过程.
 */
void dynamics_freezing()
{
    int i, j;
    int id, iu, jl, jr;
    int count;
    double offset;

    int ilo, iup, jlo, jup;

    ilo = g_center_i - g_r_new - 1;
    iup = g_center_i + g_r_new + 1;
    jlo = g_center_j - g_r_new - 1;
    jup = g_center_j + g_r_new + 1;
    g_is_fr_changed = false;

    // --- 遍历边界附近所有未结晶的格子
    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {
            if (not_snowflake(a_pic[i][j]))
            {
                id = (i + 1) % nr;
                iu = (i + nr - 1) % nr;
                jr = (j + 1) % nc;
                jl = (j + nc - 1) % nc;

                count = 0;
                if (is_snowflake(a_pic[id][j]))
                    count++;
                if (is_snowflake(a_pic[iu][j]))
                    count++;
                if (is_snowflake(a_pic[i][jl]))
                    count++;
                if (is_snowflake(a_pic[i][jr]))
                    count++;
                if (is_snowflake(a_pic[iu][jr]))
                    count++;
                if (is_snowflake(a_pic[id][jl]))
                    count++;

                // --- freezing
                if (count >= 1)
                {
                    // 占比 (1 - kappa) 的质量进入【液相】
                    offset = (1.0 - kappa) * d_dif[i][j];
                    b__fr[i][j] += offset;
                    offset = d_dif[i][j] - offset;
                    // 占比 (kappa) 的质量进入【固相】
                    c__lm[i][j] += offset;
                    // 【气相】无剩余质量
                    d_dif[i][j] = 0;
                }
            }
        }
    }
} /* dynamics_freezing() */

void dynamics_attachment()
{
    int bpic[NR_MAX][NC_MAX];

    int i, j, k;
    int id, iu, jl, jr;
    int count;
    double difmass;

    int ilo, iup, jlo, jup;

    ilo = g_center_i - g_r_new - 1;
    iup = g_center_i + g_r_new + 1;
    jlo = g_center_j - g_r_new - 1;
    jup = g_center_j + g_r_new + 1;
    g_is_fr_changed = false;

    // copy `a_pic => bpic`
    for (i = ilo; i <= iup; i++)
        for (j = jlo; j <= jup; j++)
        {
            bpic[i][j] = a_pic[i][j];
        }

    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {
            if (not_snowflake(a_pic[i][j]))
            {
                id = (i + 1) % nr;
                iu = (i + nr - 1) % nr;
                jr = (j + 1) % nc;
                jl = (j + nc - 1) % nc;

                count = 0;
                if (is_snowflake(a_pic[id][j]))
                    count++;
                if (is_snowflake(a_pic[iu][j]))
                    count++;
                if (is_snowflake(a_pic[i][jl]))
                    count++;
                if (is_snowflake(a_pic[i][jr]))
                    count++;
                if (is_snowflake(a_pic[iu][jr]))
                    count++;
                if (is_snowflake(a_pic[id][jl]))
                    count++;

                if (count >= 1)
                {
                    difmass = d_dif[i][j] + d_dif[id][j] * (1 - a_pic[id][j]) + d_dif[iu][j] * (1 - a_pic[iu][j]) +
                              d_dif[i][jl] * (1 - a_pic[i][jl]) + d_dif[i][jr] * (1 - a_pic[i][jr]) +
                              d_dif[iu][jr] * (1 - a_pic[iu][jr]) + d_dif[id][jl] * (1 - a_pic[id][jl]);

                    // [3a] count = [1, 2]
                    if (count <= 2)
                    {
                        if (b__fr[i][j] >= beta)
                        {
                            bpic[i][j] = 1;
                        }
                    }

                    // [3b]
                    if (count >= 3)
                    {
                        if ((b__fr[i][j] >= 1.0) || ((difmass <= theta) && (b__fr[i][j] >= alpha)))
                        {
                            bpic[i][j] = 1;
                        }
                    }
                    
                    // [3c]
                    if (count >= 4)
                        bpic[i][j] = 1;
                } /* if (count >= 1) */
            } /* if (not_snowflake(a_pic[i][j])) */
        }
    }

    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {
            // 使用 bpic 更新 a_pic
            if (a_pic[i][j] != bpic[i][j])
            {
                a_pic[i][j] = bpic[i][j];

                c__lm[i][j] += b__fr[i][j];
                b__fr[i][j] = 0.0;

                k = norm_inf(i - g_center_i, j - g_center_j);
                if (k > g_r_new)
                    g_r_new = k;
                if (g_r_new > 2 * nr / 3)
                    g_stop = true;
                ash[i][j] = g_par_ash;
                g_is_fr_changed = true;
            }
        }
    }

    g_par_update = 1 - g_par_update;
    if (g_r_new - g_r_old == 1)
    {
        g_par_ash = g_par_ash + 1;
        g_r_old = g_r_new;
    }
} /* dynamics_attachment() */

/** 边界上, 液相气化与固相升华.
 * 
 * - mu     液相气化率
 * - gamma  固相升华率
 */
void dynamics_melting()
{
    int i, j;
    int ilo, iup, jlo, jup;

    ilo = g_center_i - g_r_new - 1;
    iup = g_center_i + g_r_new + 1;
    jlo = g_center_j - g_r_new - 1;
    jup = g_center_j + g_r_new + 1;
    g_is_fr_changed = false;

    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {
            if (not_snowflake(a_pic[i][j]))
            {
                double b_ij, c_ij, dif_inc;

                // 占比 mu 的液相气化
                b_ij = b__fr[i][j];
                dif_inc = b_ij * mu;
                b__fr[i][j] -= dif_inc;
                d_dif[i][j] += dif_inc;

                // 占比 gamma 的固相升华
                c_ij = c__lm[i][j];
                if (c_ij > 0.0)
                {
                    dif_inc = c_ij * gam;
                    c__lm[i][j] -= dif_inc;
                    d_dif[i][j] += dif_inc;
                }
            }
        }
    }
} /* dynamics_melting() */

/**
 * Add Noise to `d_dif[]`
 */
void dynamics_add_noise()
{
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
        {
            double x = uniform_01rand();
            if (x < 0.5)
                d_dif[i][j] *= (1 + sigma);
            else
                d_dif[i][j] *= (1 - sigma);
        }
} /* dynamics_add_noise() */

/**
 * 
 * Note: 仅在从状态文件 (in file) 读取状态后调用.
 */
void dynamics_add_noise1()
{
    int i, j;
    double offset;

    if (sigma < 0)
    {
        for (i = 0; i < nr; i++)
            for (j = 0; (j < nc); j++)
            {
                if (not_snowflake(a_pic[i][j]))
                {
                    offset = sigma * d_dif[i][j];
                    d_dif[i][j] += offset;
                }
            }
    }
} /* dynamics_add_noise1() */

void dynamics()
{

    dynamics_diffusion();
    dynamics_freezing();
    dynamics_attachment();
    dynamics_melting();

    // add Noise
    if (sigma > 0.0)
        dynamics_add_noise();

    /* io_print_state(); */
} /* dynamics() */


/* ==== IO helper functions ==== */

/**
 * 跳过其他字符，直到遇见 `:`.
 * Skip any other characters until you meet a colon `:`.
 * 
 * 用于处理复制粘贴输入的情况。
 */
void io_skip()
{
    char dum;

    dum = getchar();
    while (dum != ':')
        dum = getchar();
} /* io_skip() */

/**
 * 从命令行读取参数.
 * Read parameters from the command line.
 */
void io_get_input_params() 
{
    printf("enter rho:");
    io_skip();
    scanf("%lf", &init_gas_rho);

    printf("enter h:");
    io_skip();
    scanf("%d", &init_crystal_seed_radius);

    twelve_sided = 0;
    if (init_crystal_seed_radius < 0)
    {
        init_crystal_seed_radius = -init_crystal_seed_radius;
        twelve_sided = 1;
    }

    printf("enter p:");
    io_skip();
    scanf("%lf", &init_crystal_seed_probability);

    printf("enter beta:");
    io_skip();
    scanf("%lf", &beta);

    printf("enter alpha:");
    io_skip();
    scanf("%lf", &alpha);

    printf("enter theta:");
    io_skip();
    scanf("%lf", &theta);

    printf("enter kappa:");
    io_skip();
    scanf("%lf", &kappa);

    printf("enter mu:");
    io_skip();
    scanf("%lf", &mu);

    printf("enter gamma:");
    io_skip();
    scanf("%lf", &gam);

    printf("enter sigma:");
    io_skip();
    scanf("%lf", &sigma);

    printf("enter no. of rows, L:");
    io_skip();
    scanf("%d", &nr);
    nc = nr;

    printf("size of the pixel, Zoom:");
    io_skip();
    scanf("%d", &sp);

    printf("input file:");
    io_skip();
    scanf("%s", g_in_file_path);

    printf("output file:");
    io_skip();
    scanf("%s", g_out_file_path);

    printf("graphics file:");
    io_skip();
    scanf("%s", g_graphics_file_path);

    printf("grahics viewer:");
    io_skip();
    scanf("%s", g_grahics_viewer_name);

    printf("comments (< 100 chars):");
    io_skip();
    scanf("%s", g_comments);

    printf("\n.io_get_input_params: Read params finished.\n");
} /* io_get_input_params() */

/** 
 * 读取 "input file" 里的模拟状态.
 * Read the simulation status from "input file".
 **/
void io_read_state()
{
    int i, j, k;
    double x;

    printf(".io_read_state: reading simulation state from file '%s'\n", g_in_file_path);
    g_state_file = fopen(g_in_file_path, "r");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            fscanf(g_state_file, "%lf", &x);
            d_dif[i][j] = x;
            fscanf(g_state_file, "%d", &k);
            a_pic[i][j] = k;
            fscanf(g_state_file, "%lf", &x);
            b__fr[i][j] = x;
            fscanf(g_state_file, "%d", &k);
            ash[i][j] = k;
            fscanf(g_state_file, "%lf", &x);
            c__lm[i][j] = x;
        }
    }
    fscanf(g_state_file, "%d", &k);
    g_r_old = k;
    printf(".io_read_state: read g_r_old = %d\n", g_r_old);
    fscanf(g_state_file, "%d", &k);
    g_r_new = k;
    printf(".io_read_state: read g_r_new = %d\n", g_r_new);
    fscanf(g_state_file, "%d", &k);
    g_pq = k;
    printf(".io_read_state: read g_pq = %d\n", g_pq);

    fclose(g_state_file);
    printf(".io_read_state: File read finished.\n");
} /* io_read_state() */

/** 
 * 保存模拟状态到 "output file"
 * Save the simulation state to "output file"
 **/
void io_save_state()
{
    int i, j;

    printf(".io_save_state: saving simulation state to file '%s'\n", g_out_file_path);
    g_state_file = fopen(g_out_file_path, "w");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            fprintf(g_state_file, "%.10lf %d %.10lf %d %.10lf ", d_dif[i][j], a_pic[i][j], b__fr[i][j], ash[i][j], c__lm[i][j]);
        }
    }
    fprintf(g_state_file, "%d %d ", g_r_old, g_r_new);
    fprintf(g_state_file, "%d ", g_pq);
    fclose(g_state_file);
    printf(".io_save_state: File written successfully.\n");
} /* io_save_state() */

/** 
 * 保存雪花图片到 g_graphics_file_path
 **/
void io_save_snowflake()
{
    int i, j, i1, j1, k;

    printf(".io_save_snowflake: saving snowflake image to file '%s'\n", g_graphics_file_path);
    g_state_file = fopen(g_graphics_file_path, "w");
    fprintf(g_state_file, "P3\n");

    fprintf(g_state_file, "#rho:%lf\n", init_gas_rho);
    fprintf(g_state_file, "#h:%d\n", init_crystal_seed_radius);
    fprintf(g_state_file, "#p:%lf\n", init_crystal_seed_probability);
    fprintf(g_state_file, "#beta:%lf\n", beta);
    fprintf(g_state_file, "#alpha:%lf\n", alpha);
    fprintf(g_state_file, "#theta:%lf\n", theta);
    fprintf(g_state_file, "#kappa:%lf\n", kappa);
    fprintf(g_state_file, "#mu:%lf\n", mu);
    fprintf(g_state_file, "#gamma:%lf\n", gam);
    fprintf(g_state_file, "#sigma:%lf\n", sigma);

    fprintf(g_state_file, "#L:%d\n", nr);
    fprintf(g_state_file, "#Z:%d\n", sp);

    fprintf(g_state_file, "#: no : no : no \n");

    fprintf(g_state_file, "#: %s\n", g_grahics_viewer_name);
    fprintf(g_state_file, "#: %s\n", g_comments);

    fprintf(g_state_file, "%d %d\n", nr, nr);
    fprintf(g_state_file, "255\n");
    printf("\n");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            i1 = i;
            j1 = j;
            if (g_pq % 2 == 1)
            {
                if (not_snowflake(a_pic[i1][j1]))
                {
                    k = gui_get_off_color_idx(i1, j1);
                    fprintf(g_state_file, "%d %d %d ", 
                        g_color_off[k].red * UCHAR_MAX / USHRT_MAX, 
                        g_color_off[k].green * UCHAR_MAX / USHRT_MAX,
                        g_color_off[k].blue * UCHAR_MAX / USHRT_MAX);
                }
                else
                {
                    k = gui_get_on_color_idx(i1, j1);
                    fprintf(g_state_file, "%d %d %d ", 
                        g_color_on[k].red * UCHAR_MAX / USHRT_MAX, 
                        g_color_on[k].green * UCHAR_MAX / USHRT_MAX,
                        g_color_on[k].blue * UCHAR_MAX / USHRT_MAX);
                }
            }
            else
            {
                if (not_snowflake(a_pic[i1][j1]))
                {
                    k = gui_get_off_color_idx(i1, j1);
                    fprintf(g_state_file, "%d %d %d ", 
                        g_color_off[k].red * UCHAR_MAX / USHRT_MAX, 
                        g_color_off[k].green * UCHAR_MAX / USHRT_MAX,
                        g_color_off[k].blue * UCHAR_MAX / USHRT_MAX);
                }
                else
                {
                    if (c__lm[i1][j1] > 1 + 0.5 * (beta - 1.0))
                    {
                        k = gui_get_othp_color_idx(i1, j1);
                        fprintf(g_state_file, "%d %d %d ", 
                            g_othp[k].red * UCHAR_MAX / USHRT_MAX, 
                            g_othp[k].green * UCHAR_MAX / USHRT_MAX,
                            g_othp[k].blue * UCHAR_MAX / USHRT_MAX);
                    }
                    else
                    {
                        k = gui_get_color_idx(i1, j1);
                        fprintf(g_state_file, "%d %d %d ", 
                            g_color[k].red * UCHAR_MAX / USHRT_MAX, 
                            g_color[k].green * UCHAR_MAX / USHRT_MAX,
                            g_color[k].blue * UCHAR_MAX / USHRT_MAX);
                    }
                }
            }
        }
        fprintf(g_state_file, "\n");
    }

    fclose(g_state_file);
    printf(".io_save_snowflake: File written successfully.\n");

    strcat(g_grahics_viewer_name, " ");
    strcat(g_grahics_viewer_name, g_graphics_file_path);
    popen(g_grahics_viewer_name, "r");
} /* io_save_snowflake() */

void io_print_state()
{
    int i, j;
    for (i = 0; i <= 9; i++)
    {
        for (j = 0; j <= 9; j++)
            printf("%.5lf|", d_dif[i][j]);

        printf("\n");
    }
    printf("\n");
} /* io_print_state() */

void io_check_state()
{
    int i, j, iup;

    iup = g_center_i + g_r_new + 1;
    for (i = 1; i <= iup; i++)
    {
        for (j = 1; ((j <= i) && (i + j <= nr - 1)); j++)
        {

            if ((is_snowflake(a_pic[i][j])) && (d_dif[i][j] > 0.0))
                printf("*%d %lf", a_pic[i][j], d_dif[i][j]);
        }
    }
} /* io_check_state() */


/* ==== X11 GUI helper functions ==== */

// ---- colors

void gui_blue_colors33()
{
    g_red[0]= 71; g_green[0]=204; g_blue[0]=231;
    g_red[1]= 70; g_green[1]=200; g_blue[1]=230;
    g_red[2]= 69; g_green[2]=196; g_blue[2]=229;
    g_red[3]= 68; g_green[3]=192; g_blue[3]=228;
    g_red[4]= 67; g_green[4]=188; g_blue[4]=227;
    g_red[5]= 66; g_green[5]=184; g_blue[5]=226;
    g_red[6]= 65; g_green[6]=180; g_blue[6]=225;
    g_red[7]= 64; g_green[7]=176; g_blue[7]=224;
    g_red[8]= 63; g_green[8]=172; g_blue[8]=223;
    g_red[9]= 62; g_green[9]=168; g_blue[9]=222;
    g_red[10]= 61; g_green[10]=164; g_blue[10]=221;
    g_red[11]= 60; g_green[11]=160; g_blue[11]=220;
    g_red[12]= 59; g_green[12]=156; g_blue[12]=219;
    g_red[13]= 58; g_green[13]=152; g_blue[13]=218;
    g_red[14]= 57; g_green[14]=148; g_blue[14]=217;
    g_red[15]= 56; g_green[15]=144; g_blue[15]=216;
    g_red[16]= 55; g_green[16]=140; g_blue[16]=215;
    g_red[17]= 54; g_green[17]=136; g_blue[17]=214;
    g_red[18]= 53; g_green[18]=132; g_blue[18]=213;
    g_red[19]= 52; g_green[19]=128; g_blue[19]=212;
    g_red[20]= 51; g_green[20]=124; g_blue[20]=211;
    g_red[21]= 50; g_green[21]=120; g_blue[21]=210;
    g_red[22]= 49; g_green[22]=116; g_blue[22]=209;
    g_red[23]= 48; g_green[23]=112; g_blue[23]=208;
    g_red[24]= 47; g_green[24]=108; g_blue[24]=207;
    g_red[25]= 46; g_green[25]=104; g_blue[25]=206;
    g_red[26]= 45; g_green[26]=100; g_blue[26]=205;
    g_red[27]= 44; g_green[27]=96; g_blue[27]=204;
    g_red[28]= 43; g_green[28]=92; g_blue[28]=203;
    g_red[29]= 42; g_green[29]=88; g_blue[29]=202;
    g_red[30]= 41; g_green[30]=84; g_blue[30]=201;
    g_red[31]= 40; g_green[31]=80; g_blue[31]=200;
    g_red[32]= 10; g_green[32]=20; g_blue[32]=100;
} /* gui_blue_colors33() */

void gui_braque_colors64()
{
    g_red[0]= 130; g_green[0]=166; g_blue[0]=167;
    g_red[1]= 140; g_green[1]=176; g_blue[1]=186;
    g_red[2]= 156; g_green[2]=193; g_blue[2]=200;
    g_red[3]= 163; g_green[3]=204; g_blue[3]=212;
    g_red[4]= 167; g_green[4]=213; g_blue[4]=212;
    g_red[5]= 169; g_green[5]=207; g_blue[5]=215;
    g_red[6]= 168; g_green[6]=207; g_blue[6]=211;
    g_red[7]= 157; g_green[7]=199; g_blue[7]=205;
    g_red[8]= 145; g_green[8]=162; g_blue[8]=155;
    g_red[9]= 122; g_green[9]=137; g_blue[9]=151;
    g_red[10]= 114; g_green[10]=128; g_blue[10]=136;
    g_red[11]= 101; g_green[11]=130; g_blue[11]=142;
    g_red[12]= 102; g_green[12]=130; g_blue[12]=157;
    g_red[13]= 96; g_green[13]=129; g_blue[13]=162;
    g_red[14]= 96; g_green[14]=130; g_blue[14]=165;
    g_red[15]= 98; g_green[15]=131; g_blue[15]=166;
    g_red[16]= 130; g_green[16]=166; g_blue[16]=167;
    g_red[17]= 140; g_green[17]=176; g_blue[17]=186;
    g_red[18]= 156; g_green[18]=193; g_blue[18]=200;
    g_red[19]= 163; g_green[19]=204; g_blue[19]=212;
    g_red[20]= 167; g_green[20]=213; g_blue[20]=212;
    g_red[21]= 169; g_green[21]=207; g_blue[21]=215;
    g_red[22]= 168; g_green[22]=207; g_blue[22]=211;
    g_red[23]= 157; g_green[23]=199; g_blue[23]=205;
    g_red[24]= 137; g_green[24]=175; g_blue[24]=189;
    g_red[25]= 130; g_green[25]=166; g_blue[25]=174;
    g_red[26]= 118; g_green[26]=152; g_blue[26]=164;
    g_red[27]= 118; g_green[27]=153; g_blue[27]=157;
    g_red[28]= 118; g_green[28]=158; g_blue[28]=160;
    g_red[29]= 123; g_green[29]=164; g_blue[29]=166;
    g_red[30]= 136; g_green[30]=177; g_blue[30]=177;
    g_red[31]= 146; g_green[31]=191; g_blue[31]=197;
    g_red[32]= 106; g_green[32]=140; g_blue[32]=143;
    g_red[33]= 116; g_green[33]=162; g_blue[33]=160;
    g_red[34]= 142; g_green[34]=183; g_blue[34]=185;
    g_red[35]= 184; g_green[35]=201; g_blue[35]=205;
    g_red[36]= 189; g_green[36]=224; g_blue[36]=229;
    g_red[37]= 214; g_green[37]=248; g_blue[37]=247;
    g_red[38]= 224; g_green[38]=232; g_blue[38]=240;
    g_red[39]= 201; g_green[39]=228; g_blue[39]=234;
    g_red[40]= 187; g_green[40]=218; g_blue[40]=222;
    g_red[41]= 170; g_green[41]=194; g_blue[41]=203;
    g_red[42]= 138; g_green[42]=175; g_blue[42]=176;
    g_red[43]= 115; g_green[43]=153; g_blue[43]=162;
    g_red[44]= 101; g_green[44]=137; g_blue[44]=151;
    g_red[45]= 83; g_green[45]=126; g_blue[45]=132;
    g_red[46]= 67; g_green[46]=106; g_blue[46]=111;
    g_red[47]= 70; g_green[47]=87; g_blue[47]=93;
    g_red[48]= 162; g_green[48]=197; g_blue[48]=198;
    g_red[49]= 166; g_green[49]=201; g_blue[49]=203;
    g_red[50]= 161; g_green[50]=200; g_blue[50]=203;
    g_red[51]= 153; g_green[51]=189; g_blue[51]=200;
    g_red[52]= 137; g_green[52]=175; g_blue[52]=185;
    g_red[53]= 121; g_green[53]=163; g_blue[53]=167;
    g_red[54]= 108; g_green[54]=147; g_blue[54]=151;
    g_red[55]= 106; g_green[55]=141; g_blue[55]=153;
    g_red[56]= 98; g_green[56]=131; g_blue[56]=161;
    g_red[57]= 96; g_green[57]=131; g_blue[57]=161;
    g_red[58]= 100; g_green[58]=136; g_blue[58]=155;
    g_red[59]= 111; g_green[59]=132; g_blue[59]=145;
    g_red[60]= 112; g_green[60]=134; g_blue[60]=144;
    g_red[61]= 105; g_green[61]=134; g_blue[61]=138;
    g_red[62]= 101; g_green[62]=134; g_blue[62]=138;
    g_red[63]= 105; g_green[63]=130; g_blue[63]=141;
} /* gui_braque_colors64() */

void gui_off_colors64()
{
    int i;

    for (i = 0; i < 64; i++)
    {
        g_red[i] = 129 + 2 * i;
        g_green[i] = 129 + 2 * i;
        g_blue[i] = 129 + 2 * i;
    }
} /* gui_off_colors64() */

// --- helper funcs to get color idx

/** idx for g_color_off[] */
int gui_get_off_color_idx(int i, int j)
{
    // 气相质量相对于初始浓度，[0, 64) 线性取值
    int k = floor(63.0 * (d_dif[i][j] / (init_gas_rho)));

    assert( 0 <= k && k <= 63 );
    return k;
} /* gui_get_off_color_idx(int i, int j) */

/** idx for g_color_on[] */
int gui_get_on_color_idx(int i, int j)
{
    // 液相、固相质量和
    double y = c__lm[i][j] + d_dif[i][j];
    int k = floor((33.0 * y - alpha) / (beta - alpha));

    // note: 修复初始化之后只有固相, y==0 的情况. 此时 k == -1
    if (k < 0)
        k = 0;

    if (k > 32)
        k = 32;

    assert( 0 <= k && k <= 32 );
    return k;
} /* gui_get_on_color_idx(int i, int j) */

/** idx for g_color[] */
int gui_get_color_idx(int i, int j)
{
    int k = ash[i][j];
    k = k % KAPPA_MAX;

    assert( 0 <= k && k <= 63 );
    return k;
} /* gui_get_color_idx(int i, int j) */

/** idx for g_othp[] */
int gui_get_othp_color_idx(int i, int j)
{
    int k = 0;
    double c_lm_ij = c__lm[i][j];
    assert( c_lm_ij > (1 + 0.5 * (beta - 1.0)) );
    
    if (c__lm[i][j] >= 1 + 0.2 * (beta - 1.0))
        k = 12;
    if (c__lm[i][j] >= 1 + 0.5 * (beta - 1.0))
        k = 13;
    if (c__lm[i][j] >= 1 + 0.7 * (beta - 1.0))
        k = 14;
    if (c__lm[i][j] >= beta)
        k = 15;

    assert( 0 <= k && k <= 18 );
    return k;
} /* gui_get_color_idx(int i, int j) */


void gui_X11init(int argc, char *argv[])
{
    int i;
    
    g_xDisplay = XOpenDisplay("");
    g_xScreen = DefaultScreen(g_xDisplay);
    g_xWhite = XWhitePixel(g_xDisplay, g_xScreen);
    g_xBlack = XBlackPixel(g_xDisplay, g_xScreen);

    g_xSizeHints.x = 0;
    g_xSizeHints.y = 0;

    g_xSizeHints.width = nc * sp + 100;
    g_xSizeHints.height = nr * sp + 60 + 40;

    g_xSizeHints.flags = PPosition | PSize;

    g_xWindow = XCreateSimpleWindow(g_xDisplay, DefaultRootWindow(g_xDisplay), g_xSizeHints.x, g_xSizeHints.y, g_xSizeHints.width, g_xSizeHints.height, 7, g_xBlack, g_xWhite);

    XSetWindowBorderWidth(g_xDisplay, g_xWindow, 100);

    XSetStandardProperties(g_xDisplay, g_xWindow, gui_WINDOW_NAME_STR, gui_ICON_NAME_STR, None, argv, argc, &g_xSizeHints);

    g_cmap = DefaultColormap(g_xDisplay, g_xScreen);

    gui_braque_colors64();
    for (i = 0; i < KAPPA_MAX; i++)
    {
        g_color[i].red = g_red[i] * USHRT_MAX / UCHAR_MAX;
        g_color[i].green = g_green[i] * USHRT_MAX / UCHAR_MAX;
        g_color[i].blue = g_blue[i] * USHRT_MAX / UCHAR_MAX;
        XAllocColor(g_xDisplay, g_cmap, &g_color[i]);
    }

    gui_blue_colors33();
    for (i = 0; i <= 32; i++)
    {

        g_color_on[i].red = g_red[i] * USHRT_MAX / UCHAR_MAX;
        g_color_on[i].green = g_green[i] * USHRT_MAX / UCHAR_MAX;
        g_color_on[i].blue = g_blue[i] * USHRT_MAX / UCHAR_MAX;
        XAllocColor(g_xDisplay, g_cmap, &g_color_on[i]);
    }

    gui_off_colors64();
    for (i = 0; i <= 63; i++)
    {

        g_color_off[63 - i].red = g_red[i] * USHRT_MAX / UCHAR_MAX;
        g_color_off[63 - i].green = g_green[i] * USHRT_MAX / UCHAR_MAX;
        g_color_off[63 - i].blue = g_blue[i] * USHRT_MAX / UCHAR_MAX;
        XAllocColor(g_xDisplay, g_cmap, &g_color_off[63 - i]);
    }

    XAllocNamedColor(g_xDisplay, g_cmap, "orange", &g_othp[0], &g_othp[0]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray90", &g_othp[1], &g_othp[1]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray80", &g_othp[2], &g_othp[2]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray70", &g_othp[3], &g_othp[3]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray60", &g_othp[4], &g_othp[4]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray50", &g_othp[5], &g_othp[5]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray40", &g_othp[6], &g_othp[6]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray30", &g_othp[7], &g_othp[7]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray25", &g_othp[8], &g_othp[8]);
    XAllocNamedColor(g_xDisplay, g_cmap, "gray20", &g_othp[9], &g_othp[9]);
    XAllocNamedColor(g_xDisplay, g_cmap, "black", &g_othp[10], &g_othp[10]);
    XAllocNamedColor(g_xDisplay, g_cmap, "azure", &g_othp[11], &g_othp[11]);
    XAllocNamedColor(g_xDisplay, g_cmap, "lightblue2", &g_othp[12], &g_othp[12]);
    XAllocNamedColor(g_xDisplay, g_cmap, "lightblue3", &g_othp[13], &g_othp[13]);
    XAllocNamedColor(g_xDisplay, g_cmap, "lightblue4", &g_othp[14], &g_othp[14]);
    XAllocNamedColor(g_xDisplay, g_cmap, "cornflowerblue", &g_othp[15], &g_othp[15]);
    XAllocNamedColor(g_xDisplay, g_cmap, "white", &g_othp[16], &g_othp[16]);
    XAllocNamedColor(g_xDisplay, g_cmap, "palegreen", &g_othp[17], &g_othp[17]);
    XAllocNamedColor(g_xDisplay, g_cmap, "red", &g_othp[18], &g_othp[18]);

    g_xGC = XCreateGC(g_xDisplay, g_xWindow, 0, 0);

    XSetBackground(g_xDisplay, g_xGC, g_xWhite);

    XSelectInput(g_xDisplay, g_xWindow, (ButtonPressMask | ExposureMask));

    XMapRaised(g_xDisplay, g_xWindow);

    g_exit_flag = false;

    XNextEvent(g_xDisplay, &g_xEvent);
} /* gui_X11init(int argc, char *argv[]) */

void gui_X11clean()
{
    XFreeGC(g_xDisplay, g_xGC);
    XDestroyWindow(g_xDisplay, g_xWindow);
    XCloseDisplay(g_xDisplay);
} /* gui_X11clean() */

void gui_draw_buttons()
{
    const char quitstring[] = "QUIT";
    const char pausestring[] = "pause";
    const char playstring[] = "play";
    const char savestring[] = "save";
    const char readstring[] = "read";

    XSetForeground(g_xDisplay, g_xGC, g_xBlack);
    XSetBackground(g_xDisplay, g_xGC, g_xWhite);

    XDrawRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 20, 50, nc * sp + 20, nr * sp + 20);

    XDrawRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 10, 10, 50, 20);
    XDrawRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 65, 10, 50, 20);
    XDrawRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 120, 10, 50, 20);
    XDrawRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 175, 10, 50, 20);
    XDrawRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 230, 10, 50, 20);

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 20, 25, quitstring, strlen(quitstring));

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 75, 25, pausestring, strlen(pausestring));

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 130, 25, playstring, strlen(playstring));

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 185, 25, savestring, strlen(savestring));

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 240, 25, readstring, strlen(readstring));
} /* gui_draw_buttons() */

/**
 * 实时绘制雪花到 GUI.
 * Draw snowflakes to the GUI in real time.
 */
void gui_picture_big()
{
    int i, j, k, pqn, kf;
    char pqc[10];

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            if (not_snowflake(a_pic[i][j]))
            {
                k = gui_get_off_color_idx(i, j);
                XSetForeground(g_xDisplay, g_xGC, g_color_off[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
            }
            else
            {
                k = gui_get_on_color_idx(i, j);
                XSetForeground(g_xDisplay, g_xGC, g_color_on[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
            }
        }
    }

    /* g_pq to string: sprintf(pqc, "%d", g_pq); */
    if (g_pq == 0)
    {
        pqc[0] = '0';
        pqc[1] = '\0';
    }
    else
    {
        pqn = g_pq;
        for (k = 9; k >= 0; k--)
        {

            pqc[k] = pqn % 10 + '0';
            pqn = pqn / 10;
        }
        for (kf = 0; pqc[kf] == '0'; kf++)
            ;

        for (k = 0; k <= 9 - kf; k++)
            pqc[k] = pqc[k + kf];
        pqc[10 - kf] = '\0';
    }

    XSetForeground(g_xDisplay, g_xGC, g_xBlack);
    XSetBackground(g_xDisplay, g_xGC, g_xWhite);

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 10, 45, gui_TIME_STR, strlen(gui_TIME_STR));
    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 40, 45, pqc, strlen(pqc));
} /* gui_picture_big() */

void gui_picture_rings()
{
    int i, j, k, pqn, kf;
    char pqc[10];

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            if (not_snowflake(a_pic[i][j]))
            {
                k = gui_get_off_color_idx(i, j);
                XSetForeground(g_xDisplay, g_xGC, g_color_off[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
            }
            else
            {
                k = gui_get_color_idx(i, j);
                XSetForeground(g_xDisplay, g_xGC, g_color[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
                if (c__lm[i][j] > 1 + 0.5 * (beta - 1.0))
                {
                    k = gui_get_othp_color_idx(i, j);
                    XSetForeground(g_xDisplay, g_xGC, g_othp[k].pixel);
                    XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
                }
            }
        }
    }

    /* g_pq to string: sprintf(pqc, "%d", g_pq); */
    if (g_pq == 0)
    {
        pqc[0] = '0';
        pqc[1] = '\0';
    }
    else
    {
        pqn = g_pq;
        for (k = 9; k >= 0; k--)
        {

            pqc[k] = pqn % 10 + '0';
            pqn = pqn / 10;
        }
        for (kf = 0; pqc[kf] == '0'; kf++)
            ;

        for (k = 0; k <= 9 - kf; k++)
            pqc[k] = pqc[k + kf];
        pqc[10 - kf] = '\0';
    }

    XSetForeground(g_xDisplay, g_xGC, g_xBlack);
    XSetBackground(g_xDisplay, g_xGC, g_xWhite);

    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 10, 45, gui_TIME_STR, strlen(gui_TIME_STR));
    XDrawImageString(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, 40, 45, pqc, strlen(pqc));
} /* gui_picture_rings() */


int main(int argc, char *argv[])
{
    int posx, posy;
    Window rw, cw;
    int rootx, rooty;
    unsigned int kgb;


    // ---- enter data
    io_get_input_params();

    // ---- X11 GUI
    gui_X11init(argc, argv);
    gui_draw_buttons();

    // ---- init state
    initialize();
    gui_picture_big();
    /*io_print_state(); */

    // ---- main loop
    g_pq = 0;
    while (g_exit_flag == false)
    {
        // 获取事件 | Get Events
        XNextEvent(g_xDisplay, &g_xEvent);
        /** 根据事件类型进行处理
         *  Processing based on event type
         */
        switch (g_xEvent.type)
        {

        case ButtonPress:
            // 查询鼠标位置 | Query mouse position
            XQueryPointer(g_xDisplay, g_xWindow, &rw, &cw, &rootx, &rooty, &posx, &posy, &kgb);
            if ((posx >= 10) && (posx <= 60) && (posy >= 10) && (posy <= 30))
            {
                printf("[QUIT]\n");
                g_exit_flag = true;
            }
            else if ((posx >= 65) && (posx <= 115) && (posy >= 10) && (posy <= 30))
            {
                printf("[pause]\n");
                if (g_pq % 2 == 0)
                    gui_picture_rings();
                else
                    gui_picture_big();
            }
            else if ((posx >= 120) && (posx <= 170) && (posy >= 10) && (posy <= 30))
            {
                printf("[play]\n");
                while ((XEventsQueued(g_xDisplay, QueuedAfterReading) == 0) && (g_pq != -1) && (g_stop == false))
                {
                    g_noac = 0;
                    g_pq++;
                    dynamics();
                    if (g_pq % 10 == 0)
                    {
                        gui_picture_big();
                    }
                }
            }
            else if ((posx >= 175) && (posx <= 225) && (posy >= 10) && (posy <= 30))
            {
                printf("[save] to file\n");
                io_save_state();
                io_save_snowflake();
            }
            else if ((posx >= 230) && (posx <= 280) && (posy >= 10) && (posy <= 30))
            {
                printf("[read] from file\n");
                io_read_state();
                dynamics_add_noise1();
                gui_picture_big();
            }
            else
            {
                printf("[step]\n");
                g_noac = 0;
                g_pq++;
                dynamics();
                gui_picture_big();
            }
            break;

        case Expose:
            if (g_xEvent.xexpose.count == 0)
            {
                gui_picture_big();
                gui_draw_buttons();
            }
            break;

        default:
            /* do nothing */
            break;
        } /* switch (g_xEvent.type) */
    } /* while (g_exit_flag == false) */

    // do clean befor exit.
    gui_X11clean();
    
    return 0;
} // main
