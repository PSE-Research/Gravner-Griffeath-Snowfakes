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

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#define T 1
#define F 0

#define nrmax 1002
#define ncmax 1002

#define kappamax 64

int sp;

double adif[nrmax][ncmax]; /*diffusion field*/
int apic[nrmax][ncmax];    /*indicator of snowflake sites*/
double afr[nrmax][ncmax];  /* boundary mass*/
double alm[nrmax][ncmax];  /*crystal mass*/

int ash[nrmax][ncmax]; /*rings pallette*/

/* parameters*/

double beta, kappa, mu, theta, alpha, gam, sigma;

int nr, nc;

int change, centeri, centerj;

int parupdate;

int parash;

int rold, rnew;

int rinit;
double rhorinit;

int frchange;
int twelvesided;

double rho;

int stop;

char po[30];

Display *td;
Window tw;
GC tgc;
XEvent te;
KeySym tk;
XSizeHints th;
int ts;
unsigned long tf, tb;
char tbuf[8];
int kc;
int fin;

char timestring[] = "time:";
char activestring[] = "active area:";

int noac, pq;

Colormap cmap;
XColor clr[kappamax];

XColor clron[128];
XColor clroff[128];
XColor othp[20];

int flags = {DoRed | DoGreen | DoBlue};

char in[] = "sn";
char wn[] = "digital snowflake";

int red[125], green[125], blue[125];

FILE *picf;

FILE *prof;

char infile[30];

char outfile[30];

char graphicsfile[30];

char comments[100];

void bluecolors33()
{
    int i;
    red[0]= 71; green[0]=204; blue[0]=231;
    red[1]= 70; green[1]=200; blue[1]=230;
    red[2]= 69; green[2]=196; blue[2]=229;
    red[3]= 68; green[3]=192; blue[3]=228;
    red[4]= 67; green[4]=188; blue[4]=227;
    red[5]= 66; green[5]=184; blue[5]=226;
    red[6]= 65; green[6]=180; blue[6]=225;
    red[7]= 64; green[7]=176; blue[7]=224;
    red[8]= 63; green[8]=172; blue[8]=223;
    red[9]= 62; green[9]=168; blue[9]=222;
    red[10]= 61; green[10]=164; blue[10]=221;
    red[11]= 60; green[11]=160; blue[11]=220;
    red[12]= 59; green[12]=156; blue[12]=219;
    red[13]= 58; green[13]=152; blue[13]=218;
    red[14]= 57; green[14]=148; blue[14]=217;
    red[15]= 56; green[15]=144; blue[15]=216;
    red[16]= 55; green[16]=140; blue[16]=215;
    red[17]= 54; green[17]=136; blue[17]=214;
    red[18]= 53; green[18]=132; blue[18]=213;
    red[19]= 52; green[19]=128; blue[19]=212;
    red[20]= 51; green[20]=124; blue[20]=211;
    red[21]= 50; green[21]=120; blue[21]=210;
    red[22]= 49; green[22]=116; blue[22]=209;
    red[23]= 48; green[23]=112; blue[23]=208;
    red[24]= 47; green[24]=108; blue[24]=207;
    red[25]= 46; green[25]=104; blue[25]=206;
    red[26]= 45; green[26]=100; blue[26]=205;
    red[27]= 44; green[27]=96; blue[27]=204;
    red[28]= 43; green[28]=92; blue[28]=203;
    red[29]= 42; green[29]=88; blue[29]=202;
    red[30]= 41; green[30]=84; blue[30]=201;
    red[31]= 40; green[31]=80; blue[31]=200;
    red[32]= 10; green[32]=20; blue[32]=100;
}

void braquecolors64()
{
    int i;

    red[0]= 130; green[0]=166; blue[0]=167;
    red[1]= 140; green[1]=176; blue[1]=186;
    red[2]= 156; green[2]=193; blue[2]=200;
    red[3]= 163; green[3]=204; blue[3]=212;
    red[4]= 167; green[4]=213; blue[4]=212;
    red[5]= 169; green[5]=207; blue[5]=215;
    red[6]= 168; green[6]=207; blue[6]=211;
    red[7]= 157; green[7]=199; blue[7]=205;
    red[8]= 145; green[8]=162; blue[8]=155;
    red[9]= 122; green[9]=137; blue[9]=151;
    red[10]= 114; green[10]=128; blue[10]=136;
    red[11]= 101; green[11]=130; blue[11]=142;
    red[12]= 102; green[12]=130; blue[12]=157;
    red[13]= 96; green[13]=129; blue[13]=162;
    red[14]= 96; green[14]=130; blue[14]=165;
    red[15]= 98; green[15]=131; blue[15]=166;
    red[16]= 130; green[16]=166; blue[16]=167;
    red[17]= 140; green[17]=176; blue[17]=186;
    red[18]= 156; green[18]=193; blue[18]=200;
    red[19]= 163; green[19]=204; blue[19]=212;
    red[20]= 167; green[20]=213; blue[20]=212;
    red[21]= 169; green[21]=207; blue[21]=215;
    red[22]= 168; green[22]=207; blue[22]=211;
    red[23]= 157; green[23]=199; blue[23]=205;
    red[24]= 137; green[24]=175; blue[24]=189;
    red[25]= 130; green[25]=166; blue[25]=174;
    red[26]= 118; green[26]=152; blue[26]=164;
    red[27]= 118; green[27]=153; blue[27]=157;
    red[28]= 118; green[28]=158; blue[28]=160;
    red[29]= 123; green[29]=164; blue[29]=166;
    red[30]= 136; green[30]=177; blue[30]=177;
    red[31]= 146; green[31]=191; blue[31]=197;
    red[32]= 106; green[32]=140; blue[32]=143;
    red[33]= 116; green[33]=162; blue[33]=160;
    red[34]= 142; green[34]=183; blue[34]=185;
    red[35]= 184; green[35]=201; blue[35]=205;
    red[36]= 189; green[36]=224; blue[36]=229;
    red[37]= 214; green[37]=248; blue[37]=247;
    red[38]= 224; green[38]=232; blue[38]=240;
    red[39]= 201; green[39]=228; blue[39]=234;
    red[40]= 187; green[40]=218; blue[40]=222;
    red[41]= 170; green[41]=194; blue[41]=203;
    red[42]= 138; green[42]=175; blue[42]=176;
    red[43]= 115; green[43]=153; blue[43]=162;
    red[44]= 101; green[44]=137; blue[44]=151;
    red[45]= 83; green[45]=126; blue[45]=132;
    red[46]= 67; green[46]=106; blue[46]=111;
    red[47]= 70; green[47]=87; blue[47]=93;
    red[48]= 162; green[48]=197; blue[48]=198;
    red[49]= 166; green[49]=201; blue[49]=203;
    red[50]= 161; green[50]=200; blue[50]=203;
    red[51]= 153; green[51]=189; blue[51]=200;
    red[52]= 137; green[52]=175; blue[52]=185;
    red[53]= 121; green[53]=163; blue[53]=167;
    red[54]= 108; green[54]=147; blue[54]=151;
    red[55]= 106; green[55]=141; blue[55]=153;
    red[56]= 98; green[56]=131; blue[56]=161;
    red[57]= 96; green[57]=131; blue[57]=161;
    red[58]= 100; green[58]=136; blue[58]=155;
    red[59]= 111; green[59]=132; blue[59]=145;
    red[60]= 112; green[60]=134; blue[60]=144;
    red[61]= 105; green[61]=134; blue[61]=138;
    red[62]= 101; green[62]=134; blue[62]=138;
    red[63]= 105; green[63]=130; blue[63]=141;
}

void offcolors()
{
    int i;

    for (i = 0; i < 64; i++)
    {
        red[i] = 129 + 2 * i;
        green[i] = 129 + 2 * i;
        blue[i] = 129 + 2 * i;
    }
}

double myrand()

{
    double drand48();

    return drand48();
}

int norminf(i, j)

int i, j;

{
    if (i < 0)
        i = -i;
    if (j < 0)
        j = -j;
    if (i > j)
        return i;
    else
        return j;
}

int seminorm(i, j)

int i, j;

{
    int k;

    k = i + j;
    if (k >= 0)
        return k;
    else
        return -k;
}

int shape12(x, y)

double x, y;

{
    if (x < 0)
        x = -x;
    if (y < 0)
        y = -y;
    if ((sqrt(2.0) * y < 1.0) && (y + sqrt(3.0) * x <= sqrt(2.0)) && (sqrt(2.0) * x < 1.0) &&
        (x + sqrt(3.0) * y <= sqrt(2.0)))
        return 1;
    else
        return 0;
}

int shapecircle(x, y)

double x, y;

{

    if (x * x + y * y <= 1)
        return 1;
    else
        return 0;
}

int chi(i)

int i;

{
    if (i == 0)
        return 0;
    else
        return 1;
}

double sqr(x)

double x;

{
    return x * x;
}

void plotstate()

{
    int i, j;
    for (i = 0; i <= 9; i++)
    {
        for (j = 0; j <= 9; j++)
            printf("%.5lf|", adif[i][j]);

        printf("\n");
    }
    printf("\n");
}

void check()

{

    int i, j, iup;

    iup = centeri + rnew + 1;
    for (i = 1; i <= iup; i++)
    {
        for (j = 1; ((j <= i) && (i + j <= nr - 1)); j++)
        {

            if ((apic[i][j] == 1) && (adif[i][j] > 0.0))
                printf("*%d %lf", apic[i][j], adif[i][j]);
        }
    }
}

void initialize()

{
    int t1, t2;
    int i, j, k;
    double x;
    double drand48();
    double x1, y1;

    pq = 0;

    stop = F;
    parupdate = 0;

    srand48();
    t1 = time(&t2);
    srand48();
    t1 = t1 % 1000;
    printf("seed:%d\n", t1);

    for (i = 1; i < t1; i++)
    {
        x = drand48();
    }

    centeri = nr / 2;
    centerj = nc / 2;

    rold = 0;
    rnew = 0;

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            x = myrand();

            if (twelvesided == 0)
            {
                if ((norminf(i - centeri, j - centerj) <= rinit) && (seminorm(i - centeri, j - centerj) <= rinit) &&
                    (x <= rhorinit))
                {
                    adif[i][j] = 0.0;
                    apic[i][j] = 1;
                    afr[i][j] = 1;
                    ash[i][j] = 0;
                    alm[i][j] = 0.0;
                    k = norminf(i - centeri, j - centerj);
                    if (k > rnew)
                        rnew = k;
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

                x1 = (double)(i - centeri) / rinit;
                y1 = (double)(j - centerj) / rinit;
                if ((((i - centeri) == -(j - centerj)) && (i - centeri <= 0) && (i - centeri >= -rinit)) ||
                    ((i - centeri >= 0) && (i - centeri <= rinit) && (j - centerj == 0)) ||
                    ((j - centerj <= 0) && (j - centerj >= -rinit) && (i - centeri == 0)))
                {
                    adif[i][j] = 0.0;
                    apic[i][j] = 1;
                    afr[i][j] = 0.0;
                    ash[i][j] = 0;
                    alm[i][j] = 1.0;
                    k = norminf(i - centeri, j - centerj);
                    if (k > rnew)
                        rnew = k;
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
    rold = rnew;
    parash = 1;
}

void dynamicsdif()

{

    double b[nrmax][ncmax];
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
}

void dynamicspop()

{

    double b[nrmax][ncmax];
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

            x = myrand();
            if (x < 0.5)
            {
                adif[i][j] = adif[i][j] * (1 + sigma);
            }
            else
                adif[i][j] = adif[i][j] * (1 - sigma);
        }
    }
}

void dynamicspop1()

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
}

void dynamicsunfre()

{

    double x, y, afrij;
    int i, j, k;
    int id, iu, jl, jr;
    int part;
    int count;
    double offset;
    double difmass;

    int ilo, iup, jlo, jup;

    ilo = centeri - rnew - 1;
    iup = centeri + rnew + 1;
    jlo = centerj - rnew - 1;
    jup = centerj + rnew + 1;
    frchange = F;

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
}

void dynamicsfre()

{

    int bpic[nrmax][ncmax];

    double x, y, afrij;
    int i, j, k;
    int id, iu, jl, jr;
    int part;
    int count;
    double offset;
    double nfrsum;
    double frmass, difmass;

    int ilo, iup, jlo, jup;

    ilo = centeri - rnew - 1;
    iup = centeri + rnew + 1;
    jlo = centerj - rnew - 1;
    jup = centerj + rnew + 1;
    frchange = F;

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
                k = norminf(i - centeri, j - centerj);
                if (k > rnew)
                    rnew = k;
                if (rnew > 2 * nr / 3)
                    stop = T;
                ash[i][j] = parash;
                frchange = T;
            }
        }
    }
    parupdate = 1 - parupdate;
    if (rnew - rold == 1)
    {
        parash = parash + 1;
        rold = rnew;
    }
}

void dynamicsfre1()

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

    ilo = centeri - rnew - 1;
    iup = centeri + rnew + 1;
    jlo = centerj - rnew - 1;
    jup = centerj + rnew + 1;
    frchange = F;

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
}

void dynamics()

{
    int i;

    dynamicsdif();
    dynamicsfre1();
    dynamicsfre();
    dynamicsunfre();

    if (sigma > 0.0)
        dynamicspop();

    /*plotstate(); */
}

void picturebig()

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
                XSetForeground(td, tgc, clroff[k].pixel);
                XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j * sp + 30, i * sp + 60, sp, sp);
            }
            else
            {

                y = alm[i][j] + adif[i][j];

                k = floor((33.0 * y - alpha) / (beta - alpha));
                if (k > 32)
                    k = 32;

                XSetForeground(td, tgc, clron[k].pixel);
                XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j * sp + 30, i * sp + 60, sp, sp);
            }
        }
    }

    if (pq == 0)
    {
        pqc[0] = '0';
        pqc[1] = '\0';
    }
    else
    {
        pqn = pq;
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

    XSetForeground(td, tgc, tf);
    XSetBackground(td, tgc, tb);

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 10, 45, timestring, strlen(timestring));
    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 40, 45, pqc, strlen(pqc));
}

void picturerings()

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
                XSetForeground(td, tgc, clroff[k].pixel);
                XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j * sp + 30, i * sp + 60, sp, sp);
            }
            else
            {
                k = ash[i][j];
                k = k % kappamax;
                XSetForeground(td, tgc, clr[k].pixel);
                XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j * sp + 30, i * sp + 60, sp, sp);
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

                    XSetForeground(td, tgc, othp[k].pixel);
                    XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j * sp + 30, i * sp + 60, sp, sp);
                }
            }
        }
    }

    if (pq == 0)
    {
        pqc[0] = '0';
        pqc[1] = '\0';
    }
    else
    {
        pqn = pq;
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

    XSetForeground(td, tgc, tf);
    XSetBackground(td, tgc, tb);

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 10, 45, timestring, strlen(timestring));
    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 40, 45, pqc, strlen(pqc));
}
void drawbuttons()

{
    char quitstring[] = "QUIT";
    char pausestring[] = "pause";
    char playstring[] = "play";
    char savestring[] = "save";
    char readstring[] = "read";

    XSetForeground(td, tgc, tf);
    XSetBackground(td, tgc, tb);

    XDrawRectangle(te.xexpose.display, te.xexpose.window, tgc, 20, 50, nc * sp + 20, nr * sp + 20);

    XDrawRectangle(te.xexpose.display, te.xexpose.window, tgc, 10, 10, 50, 20);
    XDrawRectangle(te.xexpose.display, te.xexpose.window, tgc, 65, 10, 50, 20);
    XDrawRectangle(te.xexpose.display, te.xexpose.window, tgc, 120, 10, 50, 20);
    XDrawRectangle(te.xexpose.display, te.xexpose.window, tgc, 175, 10, 50, 20);
    XDrawRectangle(te.xexpose.display, te.xexpose.window, tgc, 230, 10, 50, 20);

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 20, 25, quitstring, strlen(quitstring));

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 75, 25, pausestring, strlen(pausestring));

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 130, 25, playstring, strlen(playstring));

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 185, 25, savestring, strlen(savestring));

    XDrawImageString(te.xexpose.display, te.xexpose.window, tgc, 240, 25, readstring, strlen(readstring));
}

void skip()

{
    char dum;

    dum = getchar();
    while (dum != ':')
        dum = getchar();
}

void readpicture()

{
    int i, j, k;
    double x;

    picf = fopen(infile, "r");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            fscanf(picf, "%lf", &x);
            adif[i][j] = x;
            fscanf(picf, "%d", &k);
            apic[i][j] = k;
            fscanf(picf, "%lf", &x);
            afr[i][j] = x;
            fscanf(picf, "%d", &k);
            ash[i][j] = k;
            fscanf(picf, "%lf", &x);
            alm[i][j] = x;
        }
    }
    fscanf(picf, "%d", &k);
    rold = k;
    fscanf(picf, "%d", &k);
    rnew = k;
    fscanf(picf, "%d", &k);
    pq = k;
    fclose(picf);
}

void savepicture()

{
    int i, j;

    picf = fopen(outfile, "w");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            fprintf(picf, "%.10lf %d %.10lf %d %.10lf ", adif[i][j], apic[i][j], afr[i][j], ash[i][j], alm[i][j]);
        }
    }
    fprintf(picf, "%d %d ", rold, rnew);
    fprintf(picf, "%d ", pq);
    fclose(picf);
}

void savesnowflake()

{

    int i, j, i1, j1, k, pqn, kf;

    char pqc[10];

    double totalmass;

    double y;

    FILE *dum;

    picf = fopen(graphicsfile, "w");
    fprintf(picf, "P3\n");

    fprintf(picf, "#rho:%lf\n", rho);
    fprintf(picf, "#h:%d\n", rinit);
    fprintf(picf, "#p:%lf\n", rhorinit);
    fprintf(picf, "#beta:%lf\n", beta);
    fprintf(picf, "#alpha:%lf\n", alpha);
    fprintf(picf, "#theta:%lf\n", theta);
    fprintf(picf, "#kappa:%lf\n", kappa);
    fprintf(picf, "#mu:%lf\n", mu);
    fprintf(picf, "#gam:%lf\n", gam);
    fprintf(picf, "#sigma:%lf\n", sigma);

    fprintf(picf, "#L:%d\n", nr);
    fprintf(picf, "#Z:%d\n", sp);

    fprintf(picf, "#: no : no : no \n");

    fprintf(picf, "#: %s\n", po);
    fprintf(picf, "#: %s\n", comments);

    fprintf(picf, "%d %d\n", nr, nr);
    fprintf(picf, "255\n");
    printf("\n");

    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            i1 = i;
            j1 = j;
            if (pq % 2 == 1)
            {
                if (apic[i1][j1] == 0)
                {

                    k = floor(63.0 * (adif[i1][j1] / (rho)));
                    fprintf(picf, "%d %d %d ", clroff[k].red * 255 / 65535, clroff[k].green * 255 / 65535,
                            clroff[k].blue * 255 / 65535);
                }
                else
                {

                    y = alm[i1][j1] + adif[i1][j1];

                    k = floor((33.0 * y - alpha) / (beta - alpha));
                    if (k > 32)
                        k = 32;

                    fprintf(picf, "%d %d %d ", clron[k].red * 255 / 65535, clron[k].green * 255 / 65535,
                            clron[k].blue * 255 / 65535);
                }
            }
            else
            {
                if (apic[i1][j1] == 0)
                {
                    k = floor(63.0 * (adif[i1][j1] / (rho)));
                    fprintf(picf, "%d %d %d ", clroff[k].red * 255 / 65535, clroff[k].green * 255 / 65535,
                            clroff[k].blue * 255 / 65535);
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
                        fprintf(picf, "%d %d %d ", othp[k].red * 255 / 65535, othp[k].green * 255 / 65535,
                                othp[k].blue * 255 / 65535);
                    }
                    else
                    {
                        k = ash[i1][j1];
                        k = k % kappamax;
                        fprintf(picf, "%d %d %d ", clr[k].red * 255 / 65535, clr[k].green * 255 / 65535,
                                clr[k].blue * 255 / 65535);
                    }
                }
            }
        }
        fprintf(picf, "\n");
    }

    fclose(picf);

    strcat(po, " ");
    strcat(po, graphicsfile);
    dum = popen(po, "r");
}

main(argc, argv)

    int argc;
char *argv[];

{

    int i, j, k, nop;
    char ch;
    int posx, posy;

    Window rw, cw;
    int rootx, rooty;
    unsigned int kgb;

    /* enter data */

    printf("enter rho:");
    skip();
    scanf("%lf", &rho);

    printf("enter rinit:");
    skip();
    scanf("%d", &rinit);

    twelvesided = 0;
    if (rinit < 0)
    {
        rinit = -rinit;
        twelvesided = 1;
    }

    printf("enter rhorinit:");
    skip();
    scanf("%lf", &rhorinit);

    printf("enter beta:");
    skip();
    scanf("%lf", &beta);

    printf("enter alpha:");
    skip();
    scanf("%lf", &alpha);

    printf("enter theta:");
    skip();
    scanf("%lf", &theta);

    printf("enter kappa:");
    skip();
    scanf("%lf", &kappa);

    printf("enter mu:");
    skip();
    scanf("%lf", &mu);

    printf("enter gam:");
    skip();
    scanf("%lf", &gam);

    printf("enter sigma:");
    skip();
    scanf("%lf", &sigma);

    printf("enter no. of rows:");
    skip();
    scanf("%d", &nr);

    printf("\n %d\n", nr);

    nc = nr;

    printf("size of the pixel:");
    skip();
    scanf("%d", &sp);

    printf("infile:");
    skip();
    scanf("%s", infile);

    printf("outfile:");
    skip();
    scanf("%s", outfile);

    printf("graphicsfile:");
    skip();
    scanf("%s", graphicsfile);

    printf("po:");
    skip();
    scanf("%s", po);

    printf("comments:");
    skip();
    scanf("%s", comments);

    /* end data*/

    td = XOpenDisplay("");
    ts = DefaultScreen(td);
    tb = XWhitePixel(td, ts);
    tf = XBlackPixel(td, ts);

    th.x = 0;
    th.y = 0;

    th.width = nc * sp + 100;
    th.height = nr * sp + 60 + 40;

    th.flags = PPosition | PSize;

    tw = XCreateSimpleWindow(td, DefaultRootWindow(td), th.x, th.y, th.width, th.height, 7, tf, tb);

    XSetWindowBorderWidth(td, tw, 100);

    XSetStandardProperties(td, tw, wn, in, None, argv, argc, &th);

    cmap = DefaultColormap(td, ts);

    braquecolors64();
    for (i = 0; i < kappamax; i++)
    {
        clr[i].red = red[i] * 65535 / 255;
        clr[i].green = green[i] * 65535 / 255;
        clr[i].blue = blue[i] * 65535 / 255;
        XAllocColor(td, cmap, &clr[i]);
    }

    bluecolors33();

    for (i = 0; i <= 32; i++)
    {

        clron[i].red = red[i] * 65535 / 255;
        clron[i].green = green[i] * 65535 / 255;
        clron[i].blue = blue[i] * 65535 / 255;
        XAllocColor(td, cmap, &clron[i]);
    }

    offcolors();

    for (i = 0; i <= 63; i++)
    {

        clroff[63 - i].red = red[i] * 65535 / 255;
        clroff[63 - i].green = green[i] * 65535 / 255;
        clroff[63 - i].blue = blue[i] * 65535 / 255;
        XAllocColor(td, cmap, &clroff[63 - i]);
    }

    XAllocNamedColor(td, cmap, "orange", &othp[0], &othp[0]);
    XAllocNamedColor(td, cmap, "gray90", &othp[1], &othp[1]);
    XAllocNamedColor(td, cmap, "gray80", &othp[2], &othp[2]);
    XAllocNamedColor(td, cmap, "gray70", &othp[3], &othp[3]);
    XAllocNamedColor(td, cmap, "gray60", &othp[4], &othp[4]);
    XAllocNamedColor(td, cmap, "gray50", &othp[5], &othp[5]);
    XAllocNamedColor(td, cmap, "gray40", &othp[6], &othp[6]);
    XAllocNamedColor(td, cmap, "gray30", &othp[7], &othp[7]);
    XAllocNamedColor(td, cmap, "gray25", &othp[8], &othp[8]);
    XAllocNamedColor(td, cmap, "gray20", &othp[9], &othp[9]);
    XAllocNamedColor(td, cmap, "black", &othp[10], &othp[10]);
    XAllocNamedColor(td, cmap, "azure", &othp[11], &othp[11]);
    XAllocNamedColor(td, cmap, "lightblue2", &othp[12], &othp[12]);
    XAllocNamedColor(td, cmap, "lightblue3", &othp[13], &othp[13]);
    XAllocNamedColor(td, cmap, "lightblue4", &othp[14], &othp[14]);
    XAllocNamedColor(td, cmap, "cornflowerblue", &othp[15], &othp[15]);
    XAllocNamedColor(td, cmap, "white", &othp[16], &othp[16]);
    XAllocNamedColor(td, cmap, "palegreen", &othp[17], &othp[17]);
    XAllocNamedColor(td, cmap, "red", &othp[18], &othp[18]);

    tgc = XCreateGC(td, tw, 0, 0);

    XSetBackground(td, tgc, tb);

    XSelectInput(td, tw, (ButtonPressMask | ExposureMask));

    XMapRaised(td, tw);

    fin = F;

    XNextEvent(td, &te);

    drawbuttons();

    printf("creating init. st.\n");
    initialize();
    picturebig();
    /*plotstate(); */

    pq = 0;

    while (fin == F)
    {
        XNextEvent(td, &te);
        switch (te.type)
        {

        case ButtonPress:

            XQueryPointer(td, tw, &rw, &cw, &rootx, &rooty, &posx, &posy, &kgb);

            if ((posx >= 10) && (posx <= 60) && (posy >= 10) && (posy <= 30))
            {

                printf("QUIT\n");
                fin = T;
            }
            else if ((posx >= 65) && (posx <= 115) && (posy >= 10) && (posy <= 30))
            {
                printf("pause\n");
                if (pq % 2 == 0)
                    picturerings();
                else
                    picturebig();
            }
            else if ((posx >= 120) && (posx <= 170) && (posy >= 10) && (posy <= 30))
            {

                printf("play\n");

                while ((XEventsQueued(td, QueuedAfterReading) == 0) && (pq != -1) && (stop == F))
                {
                    noac = 0;
                    pq++;
                    dynamics();
                    if (pq % 10 == 0)
                    {
                        picturebig();
                    }
                }
            }
            else if ((posx >= 175) && (posx <= 225) && (posy >= 10) && (posy <= 30))
            {
                printf("save to file %s\n", outfile);
                savepicture();
                savesnowflake();
            }

            else if ((posx >= 230) && (posx <= 280) && (posy >= 10) && (posy <= 30))
            {

                printf("read from file %s\n", infile);
                readpicture();
                dynamicspop1();
                picturebig();
            }

            else
            {
                printf("step\n");

                noac = 0;
                pq++;
                dynamics();
                picturebig();
            }

            break;

        case Expose:

            if (te.xexpose.count == 0)
            {
                picturebig();
                drawbuttons();
            }
            break;
        }
    }

    XFreeGC(td, tgc);
    XDestroyWindow(td, tw);
    XCloseDisplay(td);
}
