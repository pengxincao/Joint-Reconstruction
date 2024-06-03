/*
对照相应编写的前向投影文件 MyforwardProjectC.c

*/


#include <math.h>  
#include "mex.h"  
#include <omp.h>
#define MAXX(x,y) ((x) > (y) ? (x) : (y))  
#define MIXX(x,y) ((x) < (y) ? (x) : (y))  

static void
FBPiradon(double* pPtr, double* image, double* thetaPtr, mwSize N,
    double d_pixel, double d_dec, mwSignedIndex num_dec, mwSignedIndex subpixel, mwSize numAngles)
{
    int i, j;
    int k;                 /* unsigned loop counter */
    int m, n;                 /* loop counters */
    double angle;             /* radian angle value */
    double* pr;               /* points inside output array */
    double* pixelPtr, * pixelCol;        /* 指向当前的像素 */
    double cosine, sine;                       
    /* tables for x*cos(angle) and y*sin(angle) */
    //double dr;
    //double con_k, con_b;
    mwSize numAngSubpixel = numAngles * subpixel;  /* 角度数目*子像素数目  */
    mwSize Lsubpixel = N * subpixel;  /* 像素分割后总长度  */
    double d_pos, ini_posx, ini_posy;   /*子像素的初始位置和间隔*/
    double R, delta;
    mwSignedIndex idxnum;
    double Firstdecnorm, contrb;
    double Xpos,Ypos;   /* 子像素位置   */
    double* Listsin, * Listcos; //角度的sin和cos值
    double* Xsin, * Ycos;  /* xcos ysin 值   提前对探测器投影位置归一化*/
    mwSignedIndex xsubpos, ysubpos, xsubposAng, ysubposAng;
    double leftdec, rightdec;  /*对探测器有贡献的左右边界*/
    double sub_pixel_coef; /* 子像素的贡献系数  */
    mwSignedIndex posPrk, posPr;

    sub_pixel_coef = 1/ ((double)subpixel* (double)subpixel);
    //dr = d_pixel / (double)(2 * subpixel * d_dec); /*投影分割边界*/
    //con_k = -1 / (2 * dr);                          /*投影贡献系数k*/
    //con_b = 1 / (4 * dr) + 0.5;                     /*投影贡献系数b*/
    leftdec = -1;                            //投影到探测器的左边界
    rightdec = (double)num_dec;          //投影到探测器的右边界
    d_pos = d_pixel / (double)subpixel;              
    ini_posx = -d_pixel * ((double)N / 2) + 0.5 * d_pos;
    ini_posy = d_pixel * ((double)N / 2) - 0.5 * d_pos;
    Firstdecnorm = (-(double)num_dec / 2 + 0.5);  //归一化探测器首单元中心位置

    Listsin = (double*)mxCalloc(numAngles, sizeof(double));
    Listcos = (double*)mxCalloc(numAngles, sizeof(double));

    Xsin = (double*)mxCalloc((Lsubpixel*numAngles), sizeof(double));
    Ycos = (double*)mxCalloc((Lsubpixel*numAngles), sizeof(double));


    //*Listsin,Listcos初始化
    for (k = 0; k < numAngles; k++)
    {    
        angle = thetaPtr[k];
        //合并计算，提前归一化投影位置
        *(Listsin + k) = sin(angle) / d_dec;
        *(Listcos + k) = cos(angle) / d_dec;
        //mexprintf("%f \n", *(Listsin + k));
    }

    //Xsin,Ypos初始化，顺序与大循环一致
   #pragma omp parallel for private(i,posPr,posPrk,Xpos,Ypos,k,j,sine,cosine)
    for (i = 0; i < N; i++)  //像素 
    {   
        posPr = i * numAngSubpixel; 
        Xpos = ini_posx + i * subpixel * d_pos;
        Ypos = ini_posy - i * subpixel * d_pos;
        for (k = 0; k < numAngles; k++)  //角度           
        {  
            posPrk = posPr + k * subpixel;
            sine = *(Listsin + k);
            cosine = *(Listcos + k);
            for (j = 0; j < subpixel; j++)  //子像素
            {   
               *(Xsin + posPrk + j) = (Xpos+j* d_pos) * sine;
               *(Ycos + posPrk + j) = (Ypos-j* d_pos) * cosine;
            }
        }
    }


    //对图像进行循环
   #pragma omp parallel for private(n,pixelCol,xsubpos,m,pixelPtr,ysubpos,k,xsubposAng,ysubposAng,pr,i,j,R,idxnum,delta,contrb)
    for (n = 0; n < N; n++)  /*  X坐标 */
    {
        pixelCol = image + n * N;  //指向图像当前列
        xsubpos = n * numAngSubpixel;
        for (m = 0; m < N; m++) /*  Y坐标 */
        {
            pixelPtr = pixelCol + m;//当前的像素
            ysubpos= m * numAngSubpixel;

                 //对角度进行循环
            for (k = 0; k < numAngles; k++)
            {
                xsubposAng = xsubpos + k * subpixel;
                ysubposAng = ysubpos + k * subpixel;
                pr = pPtr + k * num_dec;//指向当前角度的探测器
                    /*   子像素循环*/
                for (i = 0; i < subpixel; i++)
                {
                    for (j = 0; j < subpixel; j++)
                    {
                        R = -Xsin[xsubposAng + i] + Ycos[ysubposAng + j] - Firstdecnorm;   /* 归一化投影位置  */
                        /*根据投影位置细分情况*/
                        if (R > 0 && R < num_dec - 1)/*全部投影到探测器*/
                        {
                            idxnum = (mwSignedIndex)R;
                            delta = R - idxnum;                  
                            (*pixelPtr) += (1- delta) * pr[idxnum]+ delta * pr[idxnum + 1];
                        }
                        else if (R > leftdec && R <= 0)/*仅投影到右侧*/
                        {
                            (*pixelPtr) += (*pr);
                        }
                        else if (R >= num_dec - 1 && R < rightdec)/*仅投影到左侧*/
                        {
                            (*pixelPtr) += pr[num_dec - 1];
                        }
                    }
                }
            }

            (*pixelPtr) *= sub_pixel_coef; //投影系数
        }
    }

    /*释放内存*/
    mxFree((void*)Listsin);
    mxFree((void*)Listcos);
    mxFree((void*)Xsin);
    mxFree((void*)Ycos);
}

/* Input Arguments */
#define  P     (prhs[0])   /* 投影数据  大小*/
#define THETA  (prhs[1])  /* 投影角度 序列 */
/*     d_pixel  prhs[2]   像素宽度  */
/*      N       prhs[3]    图像大小   */
/*     d_dec    prhs[4]   探测单元宽度  */
/*     num_dec  prhs[5]   探测器单元数目  */
/*     subpixel prhs[6]   分割像素数目  subpixel>2*d_pixel/d_dec  */


/* Output Arguments */
#define image      (plhs[0])  /* 图像 */


void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    mwSize numAngles;       /* number of theta values */
    mwSize  N, numpixel;            /* input image size */
    double d_pixel;       /* 像素宽度 */
    double d_dec;       /* 探测器宽度 */
    mwSignedIndex num_dec;  /* 探测器数目 */
    mwSignedIndex subpixel;  /* 子像素分割数 */



    /* 不对输入输出进行检查 */


    /* Get THETA values  直接以弧度输入*/
    numAngles = mxGetM(THETA) * mxGetN(THETA);

    d_pixel = *(mxGetPr(prhs[2])); 
    N = *(mxGetPr(prhs[3]));/* 图像大小*/
    d_dec = *(mxGetPr(prhs[4]));
    num_dec = *(mxGetPr(prhs[5]));
    subpixel = *(mxGetPr(prhs[6]));

    /*  输入检查
      mexPrintf("N=  %d \n", N);
      mexPrintf("d_pixel=  %f \n", d_pixel);
      mexPrintf("d_dec=  %f \n", d_dec);
      mexPrintf("num_dec=  %d \n", num_dec);
      mexPrintf("subpixel=  %d \n", subpixel);
      mexPrintf("IOrigin=  %f \n", IOrigin);
   */

    image = mxCreateDoubleMatrix(N, N, mxREAL);
    FBPiradon(mxGetPr(P), mxGetPr(image), mxGetPr(THETA), N, d_pixel, d_dec, num_dec, subpixel, numAngles);

}
