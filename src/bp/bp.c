/*	
 Copyright © 2019-2019, Zhang Hongda, a64091937. All Rights Reserved.
 Description:A Standard Build C Project
 Author:Zhang Hongda
 Creation time:2019-7-28
 */
#include <stdio.h>
#include "bp.h"

BP_NN_Input_Struct BP_In_Str = {25000,35000,60000,39,41,80,80,20,0.01,0.01};
double d_alltest[Data][40]={0};
double d_all[allData][40];	//从dd[Data][40]分别得到输出和输出
double d_allgai[allData];
double d_paixu[allData];

double max_data=0;

double d_in[Data][In],d_out[Data][Out];	//从dd[Data][40]分别得到输出和输出
double dd[Data][40],dd_out[Data][Out];	//将数据文件中的输出和输出提取出来储存到这个数组中
double test_dd[testData][40];	//将测试文件中的输出和输出提取出来储存到这个数组中

double Max[testData];//每一条输出数据的最大值是多少
int    MaxNum[testData];//概率最大的那一位是几

//正向过程：d_in->w->o->u->o1->v->d_out
//                    #隐层0                       #隐层1
double w[Neuron][In],o[Neuron],u[Neuron1][Neuron],o1[Neuron1],v[Out][Neuron1];	//神经元的权值

double Maxin[In],Minin[In],Maxout[Out],Minout[Out];//最大值最小值，用于输入输出归一化

double OutputData[Out];	//神经元的输出

double resultOut[Out];	//测试的结果

double dw[Neuron][In],du[Neuron1][Neuron],dv[Out][Neuron1];	//神经元的权值的修正量

double e;//每一次训练的平均误差

double dd_e;

//初始化神经网络 
void initBPNework()
{
	FILE *fp1,*fp2,*fp3; 
	int i,j;

	//输入数据最大值，最小值 
	for(i=0; i<In; i++){

		Minin[i]=Maxin[i]=d_in[0][i];

		for(j=0; j<Data; j++){	
			Maxin[i]=Maxin[i]>d_in[j][i]?Maxin[i]:d_in[j][i]; 
			Minin[i]=Minin[i]<d_in[j][i]?Minin[i]:d_in[j][i];

		}

	}

	//输出数据的最大最小值 
	for(i=0; i<Out; i++){

		Minout[i]=Maxout[i]=d_out[0][i]; 
		for(j=0; j<Data; j++){

			Maxout[i]=Maxout[i]>d_out[j][i]?Maxout[i]:d_out[j][i]; 
			Minout[i]=Minout[i]<d_out[j][i]?Minout[i]:d_out[j][i];
		
		}

	}

	//输入数据归一化处理 ，在0~1之间
	for (i = 0; i < In; i++)

		for(j = 0; j < Data; j++) 
			d_in[j][i]=(d_in[j][i]-Minin[i]+1)/(Maxin[i]-Minin[i]+1); //归一化处理就是统一转换成分数的形式，0<x<1，便于统一处理

/*	//输出数据归一化处理 
	for (i = 0; i < Out; i++)

		for(j = 0; j < Data; j++) 
			d_out[j][i]=(d_out[j][i]-Minout[i]+1)/(Maxout[i]-Minout[i]+1);

*/
	srand((unsigned)time(NULL)); 
	//读取储存权值w的文件，没有则随机生成权值w
	if((fp1=fopen("neuronw.txt","r"))==NULL){
	
		//隐层0输入权值随机
		for (i = 0; i < Neuron; ++i)
			for (j = 0; j < In; ++j){
				w[i][j]=(rand()*2.0/RAND_MAX-1)/2; //神经元的权值,一般为负，其绝对值在0和1之间
				dw[i][j]=0;					       //神经元的权值的修正量
			}
	}
	else{
		for (i = 0; i < Neuron; ++i){ 
			for (j = 0; j < In; ++j){
				fscanf(fp1,"%lf",&w[i][j]);  //从neuronw.txt中读取隐层0的输入权值
				dw[i][j]=0;					//神经元的权值的修正量
			}	
		}
		fclose(fp1);
	}

	//读取储存权值u的文件，没有则随机生成权值u
	if((fp2=fopen("neuronu.txt","r"))==NULL){
	
		//隐层1输入权值随机  
		for (i = 0; i < Neuron; ++i)
			for (j = 0; j < Neuron1; ++j){
				u[j][i]=(rand()*2.0/RAND_MAX-1)/2; 	//神经元的权值,一般为负，其绝对值在0和1之间
				du[j][i]=0;					        //神经元的权值的修正量
			}		
	}
	else{
		for (i = 0; i < Neuron; ++i)
			for (j = 0; j < Neuron1; ++j){
				fscanf(fp2,"%lf",&u[j][i]);  //从neuronu.txt中读取隐层1的输入权值
				du[j][i]=0;					//神经元的权值的修正量
			}
		fclose(fp2);
	}

	//读取储存权值u的文件，没有则随机生成权值u	
	if((fp3=fopen("neuronv.txt","r"))==NULL){
	
     	//隐层1输出权值随机  
		for (i = 0; i < Neuron1; ++i)
			for (j = 0; j < Out; ++j){
				v[j][i]=(rand()*2.0/RAND_MAX-1)/2; 	//神经元的权值,一般为负，其绝对值在0和1之间
				dv[j][i]=0;					//神经元的权值的修正量
			}		
	}
	else{
		for (i = 0; i < Neuron1; ++i)
			for (j = 0; j < Out; ++j){
				fscanf(fp3,"%lf ",&v[j][i]);   //从neuronv.txt中读取隐层1的输出权值
				dv[j][i]=0;					//神经元的权值的修正量
			}
		fclose(fp3);
	}
}
//comput0是正向传递
void computO(int var)
{	//var 为某一个样本数据
	int i,j; 
	double sum,y=0,ei[Out];

	for (i = 0; i < Neuron; ++i){ //Neuron 为#隐层0-神经元个数

		sum=0;

		for (j = 0; j < In; ++j) 
			sum+=w[i][j]*d_in[var][j];

		o[i]=1/(1+exp(-1*sum));//sigmoid函数----是神经元输出的激活函数
								
	}

	for (i = 0; i < Neuron1; ++i){ //Neuron1 为#隐层1-神经元个数

		sum=0;

		for (j = 0; j < Neuron; ++j) 
			sum+=u[i][j]*o[j];

		o1[i]=1/(1+exp(-1*sum));	
								
	}

	for (i = 0; i < Out; ++i){

		sum=0;

		for (j = 0; j < Neuron1; ++j) 
			sum+=v[i][j]*o1[j]; 

		ei[i]=exp(sum);	//输出值，为全局变量，有Out个
		y += ei[i];
	}
	for (i = 0; i < Out; ++i){

		OutputData[i] = ei[i]/y;	//输出值，为全局变量，有Out个

	}
}


//反向传播算法部分
void backUpdate(int var)
{   //var 为某一个样本数据
	int i,j,k,m;
    double t[Neuron1],t1[Neuron];
    for (k = 0; k < Neuron1; ++k)   
    {
        t[k]=0;
        for (m = 0; m < Out; ++m){	
            t[k]+=(OutputData[m]-d_out[var][m])*v[m][k];

            dv[m][k]=A*dv[m][k]+B*(OutputData[m]-d_out[var][m])*o1[k];
            v[m][k]-=dv[m][k];
        }

        for (j = 0; j < Neuron; ++j){	
            du[k][j]=A*du[k][j]+B*t[k]*o1[k]*(1-o1[k])*o[j];
            u[k][j]-=du[k][j];
        }
		
    }

    for (j = 0; j < Neuron; ++j)   
    {
        t1[j]=0;
        for (k = 0; k < Neuron1; ++k){	
				t1[j]+=t[k]*o1[k]*(1-o1[k])*u[k][j];	
        }

        for (i = 0; i < In; ++i){	
            dw[j][i]=A*dw[j][i]+B*t1[j]*o[j]*(1-o[j])*d_in[var][i];
            w[j][i]-=dw[j][i];
        }	
    }
}
//用于测试准确率的函数，属于正向传递
void result(double var1,double var2,double var3,double var4,double var5,double var6,double var7,double var8,
			  double var9,double var10,
			  double var11,double var12,double var13,double var14,double var15,double var16,double var17,double var18,
			  double var19,double var20,
			  double var21,double var22,double var23,double var24,double var25,double var26,double var27,double var28,
			  double var29,double var30,
			  double var31,double var32,double var33,double var34,double var35,double var36,double var37,
			  double var38,double var39)
{	//如果输入或者输出变量变化，此处需要修改函数。这里是*最后测试准确度*的时候要用到的函数

	int i,j; 
	double sum,res=0;

	var1=(var1-Minin[0]+1)/(Maxin[0]-Minin[0]+1); //测试数据归一化 
	var2=(var2-Minin[1]+1)/(Maxin[1]-Minin[1]+1); //测试数据归一化
	var3=(var3-Minin[2]+1)/(Maxin[2]-Minin[2]+1); //测试数据归一化 
	var4=(var4-Minin[3]+1)/(Maxin[3]-Minin[3]+1); //测试数据归一化
	var5=(var5-Minin[4]+1)/(Maxin[4]-Minin[4]+1); //测试数据归一化 
	var6=(var6-Minin[5]+1)/(Maxin[5]-Minin[5]+1); //测试数据归一化
	var7=(var7-Minin[6]+1)/(Maxin[6]-Minin[6]+1); //测试数据归一化 
	var8=(var8-Minin[7]+1)/(Maxin[7]-Minin[7]+1); //测试数据归一化
	var9=(var9-Minin[8]+1)/(Maxin[8]-Minin[8]+1); //测试数据归一化 
	var10=(var10-Minin[9]+1)/(Maxin[9]-Minin[9]+1); //测试数据归一化

	var11=(var11-Minin[10]+1)/(Maxin[10]-Minin[10]+1); //测试数据归一化 
	var12=(var12-Minin[11]+1)/(Maxin[11]-Minin[11]+1); //测试数据归一化
	var13=(var13-Minin[12]+1)/(Maxin[12]-Minin[12]+1); //测试数据归一化 
	var14=(var14-Minin[13]+1)/(Maxin[13]-Minin[13]+1); //测试数据归一化
	var15=(var15-Minin[14]+1)/(Maxin[14]-Minin[14]+1); //测试数据归一化 
	var16=(var16-Minin[15]+1)/(Maxin[15]-Minin[15]+1); //测试数据归一化
	var17=(var17-Minin[16]+1)/(Maxin[16]-Minin[16]+1); //测试数据归一化 
	var18=(var18-Minin[17]+1)/(Maxin[17]-Minin[17]+1); //测试数据归一化
	var19=(var19-Minin[18]+1)/(Maxin[18]-Minin[18]+1); //测试数据归一化 
	var20=(var20-Minin[19]+1)/(Maxin[19]-Minin[19]+1); //测试数据归一化

	var21=(var21-Minin[20]+1)/(Maxin[20]-Minin[20]+1); //测试数据归一化 
	var22=(var22-Minin[21]+1)/(Maxin[21]-Minin[21]+1); //测试数据归一化
	var23=(var23-Minin[22]+1)/(Maxin[22]-Minin[22]+1); //测试数据归一化 
	var24=(var24-Minin[23]+1)/(Maxin[23]-Minin[23]+1); //测试数据归一化
	var25=(var25-Minin[24]+1)/(Maxin[24]-Minin[24]+1); //测试数据归一化 
	var26=(var26-Minin[25]+1)/(Maxin[25]-Minin[25]+1); //测试数据归一化
	var27=(var27-Minin[26]+1)/(Maxin[26]-Minin[26]+1); //测试数据归一化 
	var28=(var28-Minin[27]+1)/(Maxin[27]-Minin[27]+1); //测试数据归一化
	var29=(var29-Minin[28]+1)/(Maxin[28]-Minin[28]+1); //测试数据归一化 
	var30=(var30-Minin[29]+1)/(Maxin[29]-Minin[29]+1); //测试数据归一化

	var31=(var31-Minin[30]+1)/(Maxin[30]-Minin[30]+1); //测试数据归一化 
	var32=(var32-Minin[31]+1)/(Maxin[31]-Minin[31]+1); //测试数据归一化
	var33=(var33-Minin[32]+1)/(Maxin[32]-Minin[32]+1); //测试数据归一化 
	var34=(var34-Minin[33]+1)/(Maxin[33]-Minin[33]+1); //测试数据归一化
	var35=(var35-Minin[34]+1)/(Maxin[34]-Minin[34]+1); //测试数据归一化 
	var36=(var36-Minin[35]+1)/(Maxin[35]-Minin[35]+1); //测试数据归一化
	var37=(var37-Minin[36]+1)/(Maxin[36]-Minin[36]+1); //测试数据归一化 
	var38=(var38-Minin[37]+1)/(Maxin[37]-Minin[37]+1); //测试数据归一化
	var39=(var39-Minin[38]+1)/(Maxin[38]-Minin[38]+1); //测试数据归一化 




	for (i = 0; i < Neuron; ++i){

		sum=0;

		sum=w[i][0]*var1+w[i][1]*var2+w[i][2]*var3+w[i][3]*var4+w[i][4]*var5+w[i][5]*var6+w[i][6]*var7+w[i][7]*var8+w[i][8]*var9+w[i][9]*var10
			+w[i][10]*var11+w[i][11]*var12+w[i][12]*var13+w[i][13]*var14+w[i][14]*var15+w[i][15]*var16+w[i][16]*var17+w[i][17]*var18+w[i][18]*var19+w[i][19]*var20
			+w[i][20]*var21+w[i][21]*var22+w[i][22]*var23+w[i][23]*var24+w[i][24]*var25+w[i][25]*var26+w[i][26]*var27+w[i][27]*var28+w[i][28]*var29+w[i][29]*var30
			+w[i][30]*var31+w[i][31]*var32+w[i][32]*var33+w[i][33]*var34+w[i][34]*var35+w[i][35]*var36+w[i][36]*var37+w[i][37]*var38+w[i][38]*var39; //出入变量变化，增加输入 
		o[i]=1/(1+exp(-1*sum));

	}
	for (i = 0; i < Neuron1; ++i){ 

		sum=0;

		for (j = 0; j < Neuron; ++j) 
			sum+=u[i][j]*o[j];

		o1[i]=1/(1+exp(-1*sum));	
								
	}
	for (i = 0; i < Out; ++i){ 
		sum=0;
		for (j = 0; j < Neuron1; ++j)
			sum+=v[i][j]*o1[j];
		resultOut[i]=(sum*(Maxout[i]-Minout[i]+1)+Minout[i]-1);

	}

	//下面是softmax函数，将结果概率化
	for (i = 0; i < Out; ++i){ 
		
		resultOut[i]=exp(resultOut[i]);//softmax函数：Pi=ei/(所有ei之和)
		res+=resultOut[i];

	}
	for (i = 0; i < Out; ++i){ 		
		resultOut[i]=resultOut[i] /res;
	}

}



//训练神经网络的函数
void trainNetwork()
{
	int i,cc=0,yzc=0,j;
	double error=0;
	FILE *fp1; 
	if((fp1=fopen("zhunquelv_out.txt","a"))==NULL){	//读取inout.txt，不可以就弹出错误警告
		printf("can not read the  fanhua_out file\n"); 
		exit(0);
	}
	do{

		e=0;
		yzc=0;
		for (i = 0; i < Data; ++i){			//每运行一整轮Data次，就得到一个偏差百分比之和

			computO(i);						//正向传播获取输出值,输出存在 output[out]数组里面
			error =0;
			for(j = 0;j<Out;++j)
				error+=fabs( (OutputData[j]-d_out[i][j]) );  //计算所有输出的平均误差
			e+=error;	//偏差的百分比之和

			backUpdate(i);					//反向传播误差计算


		}


		for (i = 0; i < Data; ++i){			//每运行一整轮Data次，就得到一个偏差百分比之和

			result(dd[i][0], dd[i][1], dd[i][2], dd[i][3], dd[i][4], dd[i][5], dd[i][6], dd[i][7], dd[i][8], dd[i][9],
				   dd[i][10],dd[i][11],dd[i][12],dd[i][13],dd[i][14],dd[i][15],dd[i][16],dd[i][17],dd[i][18],dd[i][19],
			       dd[i][20],dd[i][21],dd[i][22],dd[i][23],dd[i][24],dd[i][25],dd[i][26],dd[i][27],dd[i][28],dd[i][29],
			       dd[i][30],dd[i][31],dd[i][32],dd[i][33],dd[i][34],dd[i][35],dd[i][36],dd[i][37],dd[i][38]);
			Max[i]=resultOut[0];//找出最大的概率是多少
			MaxNum[i]=0;
			for(j=0;j<Out;j++){
		
				if(Max[i]<resultOut[j]){
					Max[i]=resultOut[j];
					MaxNum[i]=j;
				}
			}


			if(MaxNum[i] != dd[i][39]){	//把每一组误差的每一个输出的概率写入special文件中
				yzc++;
			} 


		}


		dd_e = 	((double)yzc/(double)Data) ;
		fprintf(fp1,"%d	%lf %lf\n",cc,e/Data,((double)yzc/(double)Data) );
		printf("%d	%lf %lf\n",cc,e/Data,((double)yzc/(double)Data) );		//cc是当前训练次数，e/Data是平均偏差百分比

		cc++;
		
	}while(cc<TrainC && e/Data>0.0005);//当训练次数超过TrainC或者平均偏差百分比小于0.3%就退出循环，不再训练了
	fclose(fp1);
}

