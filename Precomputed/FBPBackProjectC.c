/*
������Ӧ��д��ǰ��ͶӰ�ļ� MyforwardProjectC.c

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
    double* pixelPtr, * pixelCol;        /* ָ��ǰ������ */
    double cosine, sine;                       
    /* tables for x*cos(angle) and y*sin(angle) */
    //double dr;
    //double con_k, con_b;
    mwSize numAngSubpixel = numAngles * subpixel;  /* �Ƕ���Ŀ*��������Ŀ  */
    mwSize Lsubpixel = N * subpixel;  /* ���طָ���ܳ���  */
    double d_pos, ini_posx, ini_posy;   /*�����صĳ�ʼλ�úͼ��*/
    double R, delta;
    mwSignedIndex idxnum;
    double Firstdecnorm, contrb;
    double Xpos,Ypos;   /* ������λ��   */
    double* Listsin, * Listcos; //�Ƕȵ�sin��cosֵ
    double* Xsin, * Ycos;  /* xcos ysin ֵ   ��ǰ��̽����ͶӰλ�ù�һ��*/
    mwSignedIndex xsubpos, ysubpos, xsubposAng, ysubposAng;
    double leftdec, rightdec;  /*��̽�����й��׵����ұ߽�*/
    double sub_pixel_coef; /* �����صĹ���ϵ��  */
    mwSignedIndex posPrk, posPr;

    sub_pixel_coef = 1/ ((double)subpixel* (double)subpixel);
    //dr = d_pixel / (double)(2 * subpixel * d_dec); /*ͶӰ�ָ�߽�*/
    //con_k = -1 / (2 * dr);                          /*ͶӰ����ϵ��k*/
    //con_b = 1 / (4 * dr) + 0.5;                     /*ͶӰ����ϵ��b*/
    leftdec = -1;                            //ͶӰ��̽��������߽�
    rightdec = (double)num_dec;          //ͶӰ��̽�������ұ߽�
    d_pos = d_pixel / (double)subpixel;              
    ini_posx = -d_pixel * ((double)N / 2) + 0.5 * d_pos;
    ini_posy = d_pixel * ((double)N / 2) - 0.5 * d_pos;
    Firstdecnorm = (-(double)num_dec / 2 + 0.5);  //��һ��̽�����׵�Ԫ����λ��

    Listsin = (double*)mxCalloc(numAngles, sizeof(double));
    Listcos = (double*)mxCalloc(numAngles, sizeof(double));

    Xsin = (double*)mxCalloc((Lsubpixel*numAngles), sizeof(double));
    Ycos = (double*)mxCalloc((Lsubpixel*numAngles), sizeof(double));


    //*Listsin,Listcos��ʼ��
    for (k = 0; k < numAngles; k++)
    {    
        angle = thetaPtr[k];
        //�ϲ����㣬��ǰ��һ��ͶӰλ��
        *(Listsin + k) = sin(angle) / d_dec;
        *(Listcos + k) = cos(angle) / d_dec;
        //mexprintf("%f \n", *(Listsin + k));
    }

    //Xsin,Ypos��ʼ����˳�����ѭ��һ��
   #pragma omp parallel for private(i,posPr,posPrk,Xpos,Ypos,k,j,sine,cosine)
    for (i = 0; i < N; i++)  //���� 
    {   
        posPr = i * numAngSubpixel; 
        Xpos = ini_posx + i * subpixel * d_pos;
        Ypos = ini_posy - i * subpixel * d_pos;
        for (k = 0; k < numAngles; k++)  //�Ƕ�           
        {  
            posPrk = posPr + k * subpixel;
            sine = *(Listsin + k);
            cosine = *(Listcos + k);
            for (j = 0; j < subpixel; j++)  //������
            {   
               *(Xsin + posPrk + j) = (Xpos+j* d_pos) * sine;
               *(Ycos + posPrk + j) = (Ypos-j* d_pos) * cosine;
            }
        }
    }


    //��ͼ�����ѭ��
   #pragma omp parallel for private(n,pixelCol,xsubpos,m,pixelPtr,ysubpos,k,xsubposAng,ysubposAng,pr,i,j,R,idxnum,delta,contrb)
    for (n = 0; n < N; n++)  /*  X���� */
    {
        pixelCol = image + n * N;  //ָ��ͼ��ǰ��
        xsubpos = n * numAngSubpixel;
        for (m = 0; m < N; m++) /*  Y���� */
        {
            pixelPtr = pixelCol + m;//��ǰ������
            ysubpos= m * numAngSubpixel;

                 //�ԽǶȽ���ѭ��
            for (k = 0; k < numAngles; k++)
            {
                xsubposAng = xsubpos + k * subpixel;
                ysubposAng = ysubpos + k * subpixel;
                pr = pPtr + k * num_dec;//ָ��ǰ�Ƕȵ�̽����
                    /*   ������ѭ��*/
                for (i = 0; i < subpixel; i++)
                {
                    for (j = 0; j < subpixel; j++)
                    {
                        R = -Xsin[xsubposAng + i] + Ycos[ysubposAng + j] - Firstdecnorm;   /* ��һ��ͶӰλ��  */
                        /*����ͶӰλ��ϸ�����*/
                        if (R > 0 && R < num_dec - 1)/*ȫ��ͶӰ��̽����*/
                        {
                            idxnum = (mwSignedIndex)R;
                            delta = R - idxnum;                  
                            (*pixelPtr) += (1- delta) * pr[idxnum]+ delta * pr[idxnum + 1];
                        }
                        else if (R > leftdec && R <= 0)/*��ͶӰ���Ҳ�*/
                        {
                            (*pixelPtr) += (*pr);
                        }
                        else if (R >= num_dec - 1 && R < rightdec)/*��ͶӰ�����*/
                        {
                            (*pixelPtr) += pr[num_dec - 1];
                        }
                    }
                }
            }

            (*pixelPtr) *= sub_pixel_coef; //ͶӰϵ��
        }
    }

    /*�ͷ��ڴ�*/
    mxFree((void*)Listsin);
    mxFree((void*)Listcos);
    mxFree((void*)Xsin);
    mxFree((void*)Ycos);
}

/* Input Arguments */
#define  P     (prhs[0])   /* ͶӰ����  ��С*/
#define THETA  (prhs[1])  /* ͶӰ�Ƕ� ���� */
/*     d_pixel  prhs[2]   ���ؿ��  */
/*      N       prhs[3]    ͼ���С   */
/*     d_dec    prhs[4]   ̽�ⵥԪ���  */
/*     num_dec  prhs[5]   ̽������Ԫ��Ŀ  */
/*     subpixel prhs[6]   �ָ�������Ŀ  subpixel>2*d_pixel/d_dec  */


/* Output Arguments */
#define image      (plhs[0])  /* ͼ�� */


void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    mwSize numAngles;       /* number of theta values */
    mwSize  N, numpixel;            /* input image size */
    double d_pixel;       /* ���ؿ�� */
    double d_dec;       /* ̽������� */
    mwSignedIndex num_dec;  /* ̽������Ŀ */
    mwSignedIndex subpixel;  /* �����طָ��� */



    /* ��������������м�� */


    /* Get THETA values  ֱ���Ի�������*/
    numAngles = mxGetM(THETA) * mxGetN(THETA);

    d_pixel = *(mxGetPr(prhs[2])); 
    N = *(mxGetPr(prhs[3]));/* ͼ���С*/
    d_dec = *(mxGetPr(prhs[4]));
    num_dec = *(mxGetPr(prhs[5]));
    subpixel = *(mxGetPr(prhs[6]));

    /*  ������
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
