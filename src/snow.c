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

#include <X11/Xlib.h>
#include <X11/Xutil.h>


#define NR_MAX 1002
#define NC_MAX 1002

#define KAPPA_MAX 64


/* ==== Input Parameters ==== */
/* --- [Initial state] */
/**
 * `rho` the density of diffusing particles elsewhere
 */
double rho;
/**
 * `h` is the radius of the initial hexagon of density `p` in the middle of the array
 * 
 * Note. h<0 is interpreted as special initialization with radius -h.
 */
int r_init;
int twelve_sided;
/**
 * density `p`
 */
double rhor_init;

/* --- [Dynamics] */
/** Dynamics params
 * 
 */
double beta, kappa, mu, theta, alpha, gam, sigma;

/** Size of the (square) LxL system.
 * 
 * nr = n_row;
 * nc = n_col;
 */
int nr, nc;
/** size of the pixel */
int sp;

/** Input/output files. */
FILE *g_picture_file;
char g_in_file_path[30];
char g_out_file_path[30];
char g_graphics_file_path[30];
char g_grahics_viewer_name[30];
char g_comments[100];


/* ==== Global Variables ==== */
// ---- initialize()
int g_pq;
int g_stop;
int g_par_update;
int g_par_ash;

// center (1, 1)
int g_center_i, g_center_j;
int g_r_old, g_r_new;

/** diffusion field */
double  adif[NR_MAX][NC_MAX];
/** indicator of snowflake sites */
int     apic[NR_MAX][NC_MAX];
/** boundary mass */
double  afr[NR_MAX][NC_MAX];
/** crystal mass */
double  alm[NR_MAX][NC_MAX];

/** rings pallette */
int     ash[NR_MAX][NC_MAX];

// ---- other global var
int g_noac;
int g_frchange;


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
int g_exit_flag;

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


/* ==== Math helper functions ==== */

/**
 * nonnegative doubleprecision floating-point values uniformly distributed
 *  over the interval [0.0, 1.0).
 */
double uniform_01rand()
{
    double drand48();

    return drand48();
} /* uniform_01rand() */

int norm_inf(int i, int j)
{
    if (i < 0)
        i = -i;
    if (j < 0)
        j = -j;
    if (i > j)
        return i;
    else
        return j;
} /* norm_inf(int i, int j) */

int semi_norm(int i, int j)
{
    int k;

    k = i + j;
    if (k >= 0)
        return k;
    else
        return -k;
} /* semi_norm(int i, int j) */


/* ==== Main simulation algorithm func ==== */

void initialize()
{
    int t1, t2;
    int i, j, k;
    double x;
    double drand48();
    double x1, y1;

    g_pq = 0;

    g_stop = false;
    g_par_update = 0;

    srand48();
    t1 = time(&t2);
    srand48();
    t1 = t1 % 1000;
    printf("seed:%d\n", t1);

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
                if ((norm_inf(i - g_center_i, j - g_center_j) <= r_init) && (semi_norm(i - g_center_i, j - g_center_j) <= r_init) &&
                    (x <= rhor_init))
                {
                    adif[i][j] = 0.0;
                    apic[i][j] = 1;
                    afr[i][j] = 1;
                    ash[i][j] = 0;
                    alm[i][j] = 0.0;
                    k = norm_inf(i - g_center_i, j - g_center_j);
                    if (k > g_r_new)
                        g_r_new = k;
                }
                else
                {
                    adif[i][j] = rho;
                    apic[i][j] = 0;
                    afr[i][j] = 0.0;
                    ash[i][j] = 0;
                    alm[i][j] = 0.0;
                }
            }
            else
            {

                x1 = (double)(i - g_center_i) / r_init;
                y1 = (double)(j - g_center_j) / r_init;
                if ((((i - g_center_i) == -(j - g_center_j)) && (i - g_center_i <= 0) && (i - g_center_i >= -r_init)) ||
                    ((i - g_center_i >= 0) && (i - g_center_i <= r_init) && (j - g_center_j == 0)) ||
                    ((j - g_center_j <= 0) && (j - g_center_j >= -r_init) && (i - g_center_i == 0)))
                {
                    adif[i][j] = 0.0;
                    apic[i][j] = 1;
                    afr[i][j] = 0.0;
                    ash[i][j] = 0;
                    alm[i][j] = 1.0;
                    k = norm_inf(i - g_center_i, j - g_center_j);
                    if (k > g_r_new)
                        g_r_new = k;
                }

                else
                {
                    adif[i][j] = rho;
                    apic[i][j] = 0;
                    afr[i][j] = 0.0;
                    ash[i][j] = 0;
                    alm[i][j] = 0.0;
                }
            }
        }
    }
    g_r_old = g_r_new;
    g_par_ash = 1;
} /* initialize() */

void dynamics_dif()
{
    double b[NR_MAX][NC_MAX];
    double x, y;
    int i, j, k;
    int id, iu, jl, jr;
    int jend;
    int part;
    int count;
    double b1, b2;
    double masscorrection;
    int nrhalf;

    for (i = 0; i < nr; i++)
        for (j = 0; (j < nc); j++)
        {
            b[i][j] = 0.0;
        }

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            if (apic[i][j] == 0)
            {
                id = (i + 1) % nr;
                iu = (i + nr - 1) % nr;
                jr = (j + 1) % nr;
                jl = (j + nr - 1) % nr;
                count = 0;
                if (apic[id][j] == 0)
                    count++;
                if (apic[iu][j] == 0)
                    count++;
                if (apic[i][jl] == 0)
                    count++;
                if (apic[i][jr] == 0)
                    count++;
                if (apic[iu][jr] == 0)
                    count++;
                if (apic[id][jl] == 0)
                    count++;

                if (count == 0)
                    b[i][j] = adif[i][j];
                else
                {

                    b[i][j] = (1.0 - (double)count / 7.0) * adif[i][j] +
                              (adif[id][j] * (1.0 - apic[id][j]) + adif[iu][j] * (1.0 - apic[iu][j]) +
                               adif[i][jl] * (1.0 - apic[i][jl]) + adif[i][jr] * (1.0 - apic[i][jr]) +
                               adif[iu][jr] * (1.0 - apic[iu][jr]) + adif[id][jl] * (1.0 - apic[id][jl])) /
                                  7.0;
                }
            }
        }
    }

    for (i = 0; i < nr; i++)
    {
        for (j = 0; (j < nc); j++)
        {
            if (apic[i][j] == 0)
                adif[i][j] = b[i][j];
        }
    }
} /* dynamics_dif() */

/**
 * Add Noise to `adif[]`
 */
void dynamics_pop()
{
    double b[NR_MAX][NC_MAX];
    double x, y;
    int i, j, k;
    int id, iu, jl, jr;
    int part;
    int count;
    double offset;

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            x = uniform_01rand();
            if (x < 0.5)
            {
                adif[i][j] = adif[i][j] * (1 + sigma);
            }
            else
                adif[i][j] = adif[i][j] * (1 - sigma);
        }
    }
} /* dynamics_pop() */

void dynamics_pop1()
{
    int i, j;
    double offset;

    if (sigma < 0)
    {
        for (i = 0; i < nr; i++)
            for (j = 0; (j < nc); j++)
            {
                if (apic[i][j] == 0)
                {
                    offset = sigma * adif[i][j];
                    adif[i][j] += offset;
                }
            }
    }
} /* dynamics_pop1() */

void dynamics_unfre()
{
    double x, y, afrij;
    int i, j, k;
    int id, iu, jl, jr;
    int part;
    int count;
    double offset;
    double difmass;

    int ilo, iup, jlo, jup;

    ilo = g_center_i - g_r_new - 1;
    iup = g_center_i + g_r_new + 1;
    jlo = g_center_j - g_r_new - 1;
    jup = g_center_j + g_r_new + 1;
    g_frchange = false;

    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {
            if (apic[i][j] == 0)
            {

                afrij = afr[i][j];
                y = afrij * mu;
                afr[i][j] = afr[i][j] - y;
                adif[i][j] = adif[i][j] + y;

                afrij = alm[i][j];
                if (afrij > 0.0)
                {
                    y = afrij * gam;
                    alm[i][j] = alm[i][j] - y;
                    adif[i][j] = adif[i][j] + y;
                }
            }
        }
    }
} /* dynamics_unfre() */

void dynamics_fre()
{
    int bpic[NR_MAX][NC_MAX];

    double x, y, afrij;
    int i, j, k;
    int id, iu, jl, jr;
    int part;
    int count;
    double offset;
    double nfrsum;
    double frmass, difmass;

    int ilo, iup, jlo, jup;

    ilo = g_center_i - g_r_new - 1;
    iup = g_center_i + g_r_new + 1;
    jlo = g_center_j - g_r_new - 1;
    jup = g_center_j + g_r_new + 1;
    g_frchange = false;

    for (i = ilo; i <= iup; i++)
        for (j = jlo; j <= jup; j++)
        {

            bpic[i][j] = apic[i][j];
        }
    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {

            if (apic[i][j] == 0)
            {

                id = (i + 1) % nr;
                iu = (i + nr - 1) % nr;
                jr = (j + 1) % nc;
                jl = (j + nc - 1) % nc;
                count = 0;
                if (apic[id][j] == 1)
                    count++;
                if (apic[iu][j] == 1)
                    count++;
                if (apic[i][jl] == 1)
                    count++;
                if (apic[i][jr] == 1)
                    count++;
                if (apic[iu][jr] == 1)
                    count++;
                if (apic[id][jl] == 1)
                    count++;

                if (count >= 1)
                {

                    difmass = adif[i][j] + adif[id][j] * (1 - apic[id][j]) + adif[iu][j] * (1 - apic[iu][j]) +
                              adif[i][jl] * (1 - apic[i][jl]) + adif[i][jr] * (1 - apic[i][jr]) +
                              adif[iu][jr] * (1 - apic[iu][jr]) + adif[id][jl] * (1 - apic[id][jl]);

                    if (count <= 2)
                    {

                        if (afr[i][j] >= beta)
                        {
                            bpic[i][j] = 1;
                        }
                    }

                    if (count >= 3)
                    {

                        if ((afr[i][j] >= 1.0) || ((difmass <= theta) && (afr[i][j] >= alpha)))
                        {
                            bpic[i][j] = 1;
                        }
                    }
                    if (count >= 4)
                        bpic[i][j] = 1;
                }
            }
        }
    }
    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {

            if (apic[i][j] != bpic[i][j])
            {
                apic[i][j] = bpic[i][j];

                alm[i][j] += afr[i][j];
                afr[i][j] = 0.0;
                k = norm_inf(i - g_center_i, j - g_center_j);
                if (k > g_r_new)
                    g_r_new = k;
                if (g_r_new > 2 * nr / 3)
                    g_stop = true;
                ash[i][j] = g_par_ash;
                g_frchange = true;
            }
        }
    }
    g_par_update = 1 - g_par_update;
    if (g_r_new - g_r_old == 1)
    {
        g_par_ash = g_par_ash + 1;
        g_r_old = g_r_new;
    }
} /* dynamics_fre() */

void dynamics_fre1()
{
    double x, y;
    int i, j, k;
    int id, iu, jl, jr;
    int part;
    int count;
    double offset;
    double frmass;
    double blockmass;

    int ilo, iup, jlo, jup;
    double epsilon;

    ilo = g_center_i - g_r_new - 1;
    iup = g_center_i + g_r_new + 1;
    jlo = g_center_j - g_r_new - 1;
    jup = g_center_j + g_r_new + 1;
    g_frchange = false;

    for (i = ilo; i <= iup; i++)
    {
        for (j = jlo; j <= jup; j++)
        {

            if (apic[i][j] == 0)
            {

                id = (i + 1) % nr;
                iu = (i + nr - 1) % nr;
                jr = (j + 1) % nc;
                jl = (j + nc - 1) % nc;
                count = 0;
                if (apic[id][j] == 1)
                    count++;
                if (apic[iu][j] == 1)
                    count++;
                if (apic[i][jl] == 1)
                    count++;
                if (apic[i][jr] == 1)
                    count++;
                if (apic[iu][jr] == 1)
                    count++;
                if (apic[id][jl] == 1)
                    count++;

                if (count >= 1)
                {
                    offset = (1.0 - kappa) * adif[i][j];
                    afr[i][j] = afr[i][j] + offset;
                    offset = adif[i][j] - offset;
                    adif[i][j] = 0;
                    alm[i][j] += offset;
                }
            }
        }
    }
} /* dynamics_fre1() */

void dynamics()
{
    int i;

    dynamics_dif();
    dynamics_fre1();
    dynamics_fre();
    dynamics_unfre();

    // add Noise
    if (sigma > 0.0)
        dynamics_pop();

    /* io_plot_state(); */
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
    scanf("%lf", &rho);

    printf("enter h:");
    io_skip();
    scanf("%d", &r_init);

    twelve_sided = 0;
    if (r_init < 0)
    {
        r_init = -r_init;
        twelve_sided = 1;
    }

    printf("enter p:");
    io_skip();
    scanf("%lf", &rhor_init);

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

    printf("\n %d\n", nr);

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
} /* io_get_input_params() */

/** 
 * 读取 g_in_file_path 里的模拟状态
 **/
void io_read_picture()
{
    int i, j, k;
    double x;

    g_picture_file = fopen(g_in_file_path, "r");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            fscanf(g_picture_file, "%lf", &x);
            adif[i][j] = x;
            fscanf(g_picture_file, "%d", &k);
            apic[i][j] = k;
            fscanf(g_picture_file, "%lf", &x);
            afr[i][j] = x;
            fscanf(g_picture_file, "%d", &k);
            ash[i][j] = k;
            fscanf(g_picture_file, "%lf", &x);
            alm[i][j] = x;
        }
    }
    fscanf(g_picture_file, "%d", &k);
    g_r_old = k;
    fscanf(g_picture_file, "%d", &k);
    g_r_new = k;
    fscanf(g_picture_file, "%d", &k);
    g_pq = k;
    fclose(g_picture_file);
} /* io_read_picture() */

/** 
 * 保存模拟状态到 g_out_file_path
 **/
void io_save_picture()
{
    int i, j;

    g_picture_file = fopen(g_out_file_path, "w");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            fprintf(g_picture_file, "%.10lf %d %.10lf %d %.10lf ", adif[i][j], apic[i][j], afr[i][j], ash[i][j], alm[i][j]);
        }
    }
    fprintf(g_picture_file, "%d %d ", g_r_old, g_r_new);
    fprintf(g_picture_file, "%d ", g_pq);
    fclose(g_picture_file);
} /* io_save_picture() */

/** 
 * 保存雪花图片到 g_graphics_file_path
 **/
void io_save_snowflake()
{

    int i, j, i1, j1, k, pqn, kf;

    char pqc[10];

    double totalmass;

    double y;

    FILE *dum;

    g_picture_file = fopen(g_graphics_file_path, "w");
    fprintf(g_picture_file, "P3\n");

    fprintf(g_picture_file, "#rho:%lf\n", rho);
    fprintf(g_picture_file, "#h:%d\n", r_init);
    fprintf(g_picture_file, "#p:%lf\n", rhor_init);
    fprintf(g_picture_file, "#beta:%lf\n", beta);
    fprintf(g_picture_file, "#alpha:%lf\n", alpha);
    fprintf(g_picture_file, "#theta:%lf\n", theta);
    fprintf(g_picture_file, "#kappa:%lf\n", kappa);
    fprintf(g_picture_file, "#mu:%lf\n", mu);
    fprintf(g_picture_file, "#gam:%lf\n", gam);
    fprintf(g_picture_file, "#sigma:%lf\n", sigma);

    fprintf(g_picture_file, "#L:%d\n", nr);
    fprintf(g_picture_file, "#Z:%d\n", sp);

    fprintf(g_picture_file, "#: no : no : no \n");

    fprintf(g_picture_file, "#: %s\n", g_grahics_viewer_name);
    fprintf(g_picture_file, "#: %s\n", g_comments);

    fprintf(g_picture_file, "%d %d\n", nr, nr);
    fprintf(g_picture_file, "255\n");
    printf("\n");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            i1 = i;
            j1 = j;
            if (g_pq % 2 == 1)
            {
                if (apic[i1][j1] == 0)
                {

                    k = floor(63.0 * (adif[i1][j1] / (rho)));
                    fprintf(g_picture_file, "%d %d %d ", g_color_off[k].red * 255 / 65535, g_color_off[k].green * 255 / 65535,
                            g_color_off[k].blue * 255 / 65535);
                }
                else
                {

                    y = alm[i1][j1] + adif[i1][j1];

                    k = floor((33.0 * y - alpha) / (beta - alpha));
                    if (k > 32)
                        k = 32;

                    fprintf(g_picture_file, "%d %d %d ", g_color_on[k].red * 255 / 65535, g_color_on[k].green * 255 / 65535,
                            g_color_on[k].blue * 255 / 65535);
                }
            }
            else
            {
                if (apic[i1][j1] == 0)
                {
                    k = floor(63.0 * (adif[i1][j1] / (rho)));
                    fprintf(g_picture_file, "%d %d %d ", g_color_off[k].red * 255 / 65535, g_color_off[k].green * 255 / 65535,
                            g_color_off[k].blue * 255 / 65535);
                }
                else
                {
                    if (alm[i1][j1] > 1 + 0.5 * (beta - 1.0))
                    {
                        if (alm[i1][j1] >= 1 + 0.2 * (beta - 1.0))
                            k = 12;
                        if (alm[i1][j1] >= 1 + 0.5 * (beta - 1.0))
                            k = 13;
                        if (alm[i1][j1] >= 1 + 0.7 * (beta - 1.0))
                            k = 14;
                        if (alm[i1][j1] >= beta)
                            k = 15;
                        fprintf(g_picture_file, "%d %d %d ", g_othp[k].red * 255 / 65535, g_othp[k].green * 255 / 65535,
                                g_othp[k].blue * 255 / 65535);
                    }
                    else
                    {
                        k = ash[i1][j1];
                        k = k % KAPPA_MAX;
                        fprintf(g_picture_file, "%d %d %d ", g_color[k].red * 255 / 65535, g_color[k].green * 255 / 65535,
                                g_color[k].blue * 255 / 65535);
                    }
                }
            }
        }
        fprintf(g_picture_file, "\n");
    }

    fclose(g_picture_file);

    strcat(g_grahics_viewer_name, " ");
    strcat(g_grahics_viewer_name, g_graphics_file_path);
    dum = popen(g_grahics_viewer_name, "r");
} /* io_save_snowflake() */

void io_plot_state()
{
    int i, j;
    for (i = 0; i <= 9; i++)
    {
        for (j = 0; j <= 9; j++)
            printf("%.5lf|", adif[i][j]);

        printf("\n");
    }
    printf("\n");
} /* io_plot_state() */

void io_check_state()
{
    int i, j, iup;

    iup = g_center_i + g_r_new + 1;
    for (i = 1; i <= iup; i++)
    {
        for (j = 1; ((j <= i) && (i + j <= nr - 1)); j++)
        {

            if ((apic[i][j] == 1) && (adif[i][j] > 0.0))
                printf("*%d %lf", apic[i][j], adif[i][j]);
        }
    }
} /* io_check_state() */


/* ==== X11 GUI helper functions ==== */

// ---- colors

void gui_blue_colors33()
{
    int i;
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
    int i;

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
        g_color[i].red = g_red[i] * 65535 / 255;
        g_color[i].green = g_green[i] * 65535 / 255;
        g_color[i].blue = g_blue[i] * 65535 / 255;
        XAllocColor(g_xDisplay, g_cmap, &g_color[i]);
    }

    gui_blue_colors33();

    for (i = 0; i <= 32; i++)
    {

        g_color_on[i].red = g_red[i] * 65535 / 255;
        g_color_on[i].green = g_green[i] * 65535 / 255;
        g_color_on[i].blue = g_blue[i] * 65535 / 255;
        XAllocColor(g_xDisplay, g_cmap, &g_color_on[i]);
    }

    gui_off_colors64();

    for (i = 0; i <= 63; i++)
    {

        g_color_off[63 - i].red = g_red[i] * 65535 / 255;
        g_color_off[63 - i].green = g_green[i] * 65535 / 255;
        g_color_off[63 - i].blue = g_blue[i] * 65535 / 255;
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

void gui_picture_big()
{
    int i, j, k, pqn, kf;

    double y;

    char pqc[10];

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {

            if (apic[i][j] == 0)
            {

                k = floor(63.0 * (adif[i][j] / (rho)));
                XSetForeground(g_xDisplay, g_xGC, g_color_off[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
            }
            else
            {

                y = alm[i][j] + adif[i][j];

                k = floor((33.0 * y - alpha) / (beta - alpha));
                if (k > 32)
                    k = 32;

                XSetForeground(g_xDisplay, g_xGC, g_color_on[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
            }
        }
    }

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

            if (apic[i][j] == 0)
            {

                k = floor(63.0 * (adif[i][j] / (rho)));
                XSetForeground(g_xDisplay, g_xGC, g_color_off[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
            }
            else
            {
                k = ash[i][j];
                k = k % KAPPA_MAX;
                XSetForeground(g_xDisplay, g_xGC, g_color[k].pixel);
                XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
                if (alm[i][j] > 1 + 0.5 * (beta - 1.0))
                {
                    if (alm[i][j] >= 1 + 0.2 * (beta - 1.0))
                        k = 12;
                    if (alm[i][j] >= 1 + 0.5 * (beta - 1.0))
                        k = 13;
                    if (alm[i][j] >= 1 + 0.7 * (beta - 1.0))
                        k = 14;
                    if (alm[i][j] >= beta)
                        k = 15;

                    XSetForeground(g_xDisplay, g_xGC, g_othp[k].pixel);
                    XFillRectangle(g_xEvent.xexpose.display, g_xEvent.xexpose.window, g_xGC, j * sp + 30, i * sp + 60, sp, sp);
                }
            }
        }
    }

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

    printf("creating init. st.\n");
    initialize();
    gui_picture_big();
    /*io_plot_state(); */

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
                printf("QUIT\n");
                g_exit_flag = true;
            }
            else if ((posx >= 65) && (posx <= 115) && (posy >= 10) && (posy <= 30))
            {
                printf("pause\n");
                if (g_pq % 2 == 0)
                    gui_picture_rings();
                else
                    gui_picture_big();
            }
            else if ((posx >= 120) && (posx <= 170) && (posy >= 10) && (posy <= 30))
            {
                printf("play\n");
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
                printf("save to file %s\n", g_out_file_path);
                io_save_picture();
                io_save_snowflake();
            }
            else if ((posx >= 230) && (posx <= 280) && (posy >= 10) && (posy <= 30))
            {
                printf("read from file %s\n", g_in_file_path);
                io_read_picture();
                dynamics_pop1();
                gui_picture_big();
            }
            else
            {
                printf("step\n");
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
