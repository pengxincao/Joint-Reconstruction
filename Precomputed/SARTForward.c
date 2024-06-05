/*����MATLAB radon C�����д
�ɱ����غ�̽������С


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
    double d_pos, ini_posx, ini_posy;   /*�����صĳ�ʼλ�úͼ��*/
    double R, R1, R2;
    //mwSignedIndex idxnum;
    double Firstdecnorm, contrb;
    double *Xpos, *Ypos;   /* ������λ��   */
    double *Xsin, *Ycos;  /* xcos ysin ֵ   */
    mwSignedIndex xsubpos, ysubpos;
    //double leftdec, rightdec;  /*��̽�����й��׵����ұ߽�*/
    double LimitD; //��ǰ���ͶӰ����
    int beginnum, endnum;  //ͶӰ��ʼ����
    //double W;
    Xpos = (double*)mxCalloc(N, sizeof(double));
    Ypos = (double*)mxCalloc(N, sizeof(double));
    Xsin= (double*)mxCalloc(numAngles *N, sizeof(double));
    Ycos = (double*)mxCalloc(numAngles *N, sizeof(double));



    Firstdecnorm = (-(double)num_dec / 2 + 0.5) ;  /*��һ��̽�����׵�Ԫ����λ��*/

    //��ǰ��̽����λ�ù�һ��
    d_pos = d_pixel / d_dec;
    ini_posx = -d_pixel * ((double)N / 2 - 0.5) / d_dec;
    ini_posy = d_pixel * ((double)N / 2 - 0.5) / d_dec;
    


    for (i = 0; i < N; i++)
    {   
        *(Xpos+i) = ini_posx + i * d_pos;
        *(Ypos+i) = ini_posy - i * d_pos;
    }
    /*d��Xsin,Ycos��ʼ����λ�����ȣ��Ƕ�ѭ������*/
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
    

    /*���в���*/
    #pragma omp parallel for   private(k,pr,posPr,LTp,pixelPtr,pixelPtrmask,LimitD,n,m,pixel,maskv,R,R1,R2,i,beginnum,endnum)
    for (k = 0; k < numAngles; k++) 
    {
        pr = pPtr + k * num_dec;  /* pointer to the top of the output column */      
        posPr = k * N;
        LTp = LT + k * LT_N;
        /* MATLAB����Ϊ�����ȴ洢����ѭ��ʱX�������꣩���⣬Y�������꣩����  */
        pixelPtr = image;      //ָ��ͼ��
        pixelPtrmask = mask;  //ָ��mask
        LimitD = *(LTD+k);    //�˴����ͶӰ����
        /*     ��ÿ������ѭ��     */
        for (n = 0; n < N; n++)  /*  X���� */
        {
            for (m = 0; m < N; m++) /*  Y���� */
            {
                pixel = *pixelPtr++;
                maskv = *pixelPtrmask++;
                if (maskv != 0)  //�˴���������
                {    /*   ������ѭ��*/
 
                  R = -*(Xsin+ posPr + n) + *(Ycos + posPr + m) - Firstdecnorm;   /* ��һ��ͶӰλ��  */
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

    /*�ͷ��ڴ�*/
    mxFree((void*)Xpos);
    mxFree((void*)Ypos);
    mxFree((void*)Xsin);
    mxFree((void*)Ycos);

}

/* Input Arguments */
#define image  (prhs[0])   /* ͼ��    N*N��С*/
#define THETA  (prhs[1])  /* ͶӰ�Ƕ� ����ֵ ���� */
/*     d_pixel  prhs[2]   ���ؿ��  */
/*     d_dec    prhs[3]   ̽�ⵥԪ���  */
/*     num_dec  prhs[4]   ̽������Ԫ��Ŀ  */
//     mask     prhs[5]   ͼ����Ĥ  
//     LT       prhs[6]   ��ѯ��
//     FirstD   prhs[7]   ��ѯ����ʼ����
//     deltaD   prhs[8]   ��ѯ��������
//     LTD      prhs[9]   ���ͶӰ����


/* Output Arguments */
#define P      (plhs[0])  /* ͶӰ */


void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    mwSize numAngles;       /* number of theta values */
    /*double* thetaPtr;       /* pointer to theta values in radians */
    double* pr1, * pr2;      /* double pointers used in loop */
  /*  double deg2rad;         /* conversion factor */
    mwSize k;               /* loop counter */
    mwSize  N, numpixel,LT_N;            /* input image size */
    double d_pixel;       /* ���ؿ�� */
    double d_dec;       /* ̽������� */
    mwSignedIndex num_dec;  /* ̽������Ŀ */
    double* mask, * LT, * LTD;
    double FirstD, deltaD;
    /* ��������������м�� */
 

    /* Get THETA values */
    numAngles = mxGetM(THETA) * mxGetN(THETA);

    N = mxGetN(image); /* ͼ���С  */
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
