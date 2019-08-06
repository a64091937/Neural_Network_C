/*	
 Copyright © 2019-2019, Zhang Hongda, a64091937. All Rights Reserved.
 Description:A Standard Build C Project
 Author:Zhang Hongda
 Creation time:2019-7-28
 */
 
#ifndef BP_H
#define BP_H
#include <stdio.h> 
#include <time.h> 
#include <math.h> 
#include <stdlib.h>
#include "def.h"

typedef struct {
    INT32U ulDataNum;    //样本数据量
    INT32U ulTestDataNum;//测试数据                                
    INT32U ulAllData; //1数据库
    INT16U NN_In;		// 输入
    INT16U NN_Out;     //输出，正确的输出那一位为1，测试的输出为概率
    INT16U NN_Neuron;	//神经元 #隐层0
    INT16U NN_Neuron1;	//神经元 #隐层1
    INT32U NN_TrainC;	//训练次数
    FL64S NN_A;        //动量值
    FL64S NN_B;        //学习速率    
    
} BP_NN_Input_Struct;

typedef struct {
    INT32U ulDataNum;        //动量值
    INT32U ulTestDataNum;        //学习速率    
    INT32U ulAllData;        //动量值  
    INT16U uwDataLen;
    INT8S scDataFile[20];
} BP_In_Data_File_Struct;


extern BP_NN_Input_Struct BP_In_Str;

#define Data	    25000	//样本数据量
#define testData	35000	//测试数据                                
#define allData	    60000  //1数据库

#define In 39		// 输入
#define Out 41		//输出，正确的输出那一位为1，测试的输出为概率
#define Neuron	 80	//神经元 #隐层0
#define Neuron1	 80	//神经元 #隐层1
#define TrainC	20	//训练次数
#define A	0.01	//动量值
#define B	0.01	//学习速率


extern double d_all[allData][40];	//从dd[Data][40]分别得到输出和输出
extern double d_allgai[allData];
extern double d_paixu[allData];

extern double max_data;

extern double d_in[Data][In],d_out[Data][Out];	//从dd[Data][40]分别得到输出和输出
extern double dd[Data][40],dd_out[Data][Out];	//将数据文件中的输出和输出提取出来储存到这个数组中
extern double test_dd[testData][40];	//将测试文件中的输出和输出提取出来储存到这个数组中

extern double Max[testData];//每一条输出数据的最大值是多少
extern int    MaxNum[testData];//概率最大的那一位是几

//正向过程：d_in->w->o->u->o1->v->d_out
//                    #隐层0                       #隐层1
extern double w[Neuron][In],o[Neuron],u[Neuron1][Neuron],o1[Neuron1],v[Out][Neuron1];	//神经元的权值

extern double Maxin[In],Minin[In],Maxout[Out],Minout[Out];//最大值最小值，用于输入输出归一化

extern double OutputData[Out];	//神经元的输出

extern double resultOut[Out];	//测试的结果

extern double dw[Neuron][In],du[Neuron1][Neuron],dv[Out][Neuron1];	//神经元的权值的修正量

extern double e;//每一次训练的平均误差

extern double dd_e;




typedef struct {
    
    
} BP_NN_Internal_Struct;

typedef struct {
    
    
} BP_NN_Output_Struct;

void initBPNework();
void computO(int var);
void backUpdate(int var);
void result(double var1,double var2,double var3,double var4,double var5,double var6,double var7,double var8,
			  double var9,double var10,
			  double var11,double var12,double var13,double var14,double var15,double var16,double var17,double var18,
			  double var19,double var20,
			  double var21,double var22,double var23,double var24,double var25,double var26,double var27,double var28,
			  double var29,double var30,
			  double var31,double var32,double var33,double var34,double var35,double var36,double var37,
			  double var38,double var39);
void trainNetwork();              


#endif