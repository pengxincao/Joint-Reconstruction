/*根据MATLAB radon C代码改写
可变像素和探测器大小


*/


#include <math.h>  
#include "mex.h"  
#include <omp.h>
#define MAXX(x,y) ((x) > (y) ? (x) : (y))  
#define MIXX(x,y) ((x) < (y) ? (x) : (y))  


double
interW(double* LTp, double R, double FirstD, double deltaD)
{
    int idnum;
    double delta,K;
    double RR;
    RR = (R - FirstD) / deltaD;
    idnum = (int)RR;
    delta = RR - idnum;
    K = LTp[idnum] * (1 - delta) + LTp[idnum + 1] * delta;
    return K;
}

static void
radon(double* pPtr, double* image, double* thetaPtr, mwSize N,
    double d_pixel, double d_dec, mwSignedIndex num_dec, mwSize numAngles, double* mask, double* LT, mwSize LT_N, double FirstD, double deltaD, double* LTD)
{
    int i, j, k, posPr;             
    mwSize m, n;                 /* loop counters */
    double angle;             /* radian angle value */
    double cosine, sine;      /* cosine and sine of current angle */
    double* pr;               /* points inside output array */
    double* pixelPtr, * pixelPtrmask, * LTp;         /* points inside input array */
    double pixel,maskv;             /* current pixel value */
    double d_pos, ini_posx, ini_posy;   /*子像素的初始位置和间隔*/
    double R, R1, R2;
    //mwSignedIndex idxnum;
    double Firstdecnorm, contrb;
    double *Xpos, *Ypos;   /* 子像素位置   */
    double *Xsin, *Ycos;  /* xcos ysin 值   */
    mwSignedIndex xsubpos, ysubpos;
    //double leftdec, rightdec;  /*对探测器有贡献的左右边界*/
    double LimitD; //当前最大投影距离
    int beginnum, endnum;  //投影起始坐标
    //double W;
    Xpos = (double*)mxCalloc(N, sizeof(double));
    Ypos = (double*)mxCalloc(N, sizeof(double));
    Xsin= (double*)mxCalloc(numAngles *N, sizeof(double));
    Ycos = (double*)mxCalloc(numAngles *N, sizeof(double));



    Firstdecnorm = (-(double)num_dec / 2 + 0.5) ;  /*归一化探测器首单元中心位置*/

    //提前对探测器位置归一化
    d_pos = d_pixel / d_dec;
    ini_posx = -d_pixel * ((double)N / 2 - 0.5) / d_dec;
    ini_posy = d_pixel * ((double)N / 2 - 0.5) / d_dec;
    


    for (i = 0; i < N; i++)
    {   
        *(Xpos+i) = ini_posx + i * d_pos;
        *(Ypos+i) = ini_posy - i * d_pos;
    }
    /*d对Xsin,Ycos初始化，位置优先，角度循环在外*/
   #pragma omp parallel for private(k,angle,cosine,sine,posPr,i)
    for (k = 0; k < numAngles; k++)
    {   
        angle = thetaPtr[k];
        cosine = cos(angle);
        sine = sin(angle);
        posPr = k * N;
        for (i = 0; i <N; i++)
        {  
            *(Xsin+posPr + i) = *(Xpos + i) * sine;
            *(Ycos+posPr + i) = *(Ypos + i) * cosine;
        }
    }
    

    /*并行操作*/
    #pragma omp parallel for   private(k,pr,posPr,LTp,pixelPtr,pixelPtrmask,LimitD,n,m,pixel,maskv,R,R1,R2,i,beginnum,endnum)
    for (k = 0; k < numAngles; k++) 
    {
        pr = pPtr + k * num_dec;  /* pointer to the top of the output column */      
        posPr = k * N;
        LTp = LT + k * LT_N;
        /* MATLAB矩阵为列优先存储，即循环时X（列坐标）在外，Y（行坐标）在内  */
        pixelPtr = image;      //指向图像
        pixelPtrmask = mask;  //指向mask
        LimitD = *(LTD+k);    //此处最大投影距离
        /*     对每个像素循环     */
        for (n = 0; n < N; n++)  /*  X坐标 */
        {
            for (m = 0; m < N; m++) /*  Y坐标 */
            {
                pixel = *pixelPtr++;
                maskv = *pixelPtrmask++;
                if (maskv != 0)  //此处存在物质
                {    /*   子像素循环*/
 
                  R = -*(Xsin+ posPr + n) + *(Ycos + posPr + m) - Firstdecnorm;   /* 归一化投影位置  */
                  R1 = R - LimitD;
                  R2 = R + LimitD;
                  if (R2 > 0 && R1 < num_dec - 1)
                  {
                      if (R1 < 0)
                      { 
                          beginnum = 0;
                      }
                      else
                      {
                          beginnum = (int)R1+1;
                      }

                      endnum = MIXX(num_dec - 1,(int)R2);

                      for (i = beginnum; i <= endnum; i++)
                      {   
                          //W = interW(LTp, i - R, FirstD, deltaD);
                          *(pr + i) += pixel* interW(LTp, i - R, FirstD, deltaD);
                      }


                  }
                }
            }
        }
    }

    /*释放内存*/
    mxFree((void*)Xpos);
    mxFree((void*)Ypos);
    mxFree((void*)Xsin);
    mxFree((void*)Ycos);

}

/* Input Arguments */
#define image  (prhs[0])   /* 图像    N*N大小*/
#define THETA  (prhs[1])  /* 投影角度 弧度值 序列 */
/*     d_pixel  prhs[2]   像素宽度  */
/*     d_dec    prhs[3]   探测单元宽度  */
/*     num_dec  prhs[4]   探测器单元数目  */
//     mask     prhs[5]   图像掩膜  
//     LT       prhs[6]   查询表
//     FirstD   prhs[7]   查询表起始距离
//     deltaD   prhs[8]   查询表间隔距离
//     LTD      prhs[9]   最大投影距离


/* Output Arguments */
#define P      (plhs[0])  /* 投影 */


void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    mwSize numAngles;       /* number of theta values */
    /*double* thetaPtr;       /* pointer to theta values in radians */
    double* pr1, * pr2;      /* double pointers used in loop */
  /*  double deg2rad;         /* conversion factor */
    mwSize k;               /* loop counter */
    mwSize  N, numpixel,LT_N;            /* input image size */
    double d_pixel;       /* 像素宽度 */
    double d_dec;       /* 探测器宽度 */
    mwSignedIndex num_dec;  /* 探测器数目 */
    double* mask, * LT, * LTD;
    double FirstD, deltaD;
    /* 不对输入输出进行检查 */
 

    /* Get THETA values */
    numAngles = mxGetM(THETA) * mxGetN(THETA);

    N = mxGetN(image); /* 图像大小  */
    d_pixel = *(mxGetPr(prhs[2]));
    d_dec = *(mxGetPr(prhs[3]));
    num_dec = *(mxGetPr(prhs[4]));
    //subpixel = *(mxGetPr(prhs[5]));
    mask = mxGetPr(prhs[5]);
    LT = mxGetPr(prhs[6]);
    LT_N= mxGetM(prhs[6]);
    FirstD = *(mxGetPr(prhs[7]));
    deltaD = *(mxGetPr(prhs[8]));
    LTD = mxGetPr(prhs[9]);




   
    
    P = mxCreateDoubleMatrix(num_dec, numAngles, mxREAL);
    radon(mxGetPr(P), mxGetPr(image), mxGetPr(THETA), N, d_pixel, d_dec, num_dec, numAngles, mask, LT, LT_N, FirstD, deltaD, LTD);
 
}
