/*	
 Copyright © 2019-2019, Zhang Hongda, a64091937. All Rights Reserved.
 Description:A Standard Build C Project
 Author:Zhang Hongda
 Creation time:2019-7-28
 */
#include "file_manager.h"


//读取数据
void readData()
{
	FILE *fp1,*fp2,*fp3,*fp4; 
	int i,j,k,cx,outi;
	double temp;

	if((fp1=fopen("60000.txt","r"))==NULL){	//读取inout.txt，不可以就弹出错误警告
		printf("can not read the  inout file\n"); 
		exit(0);
	}

	for(i=0;i<allData;i++){
		for(j=0; j<40; j++){
			fscanf(fp1,"%lf",&d_all[i][j]);    //将读取到的数据先存储于dd[i][j]数组中
		}
	}	
	fclose(fp1);

//读取完所有数据了，下面要随机抽取了
	srand((unsigned)time(NULL)); 

	for (i = 0; i < allData; ++i) {
		d_allgai[i]=(rand()*2)/2+100; //
		d_paixu[i]=d_allgai[i];
	}

    for (j = 0; j < allData - 1; j++) {
        for (i = 0; i < allData - 1 - j; i++) {
            if(d_paixu[i] > d_paixu[i+1]) {
                temp = d_paixu[i];
                d_paixu[i] = d_paixu[i+1];
                d_paixu[i+1] = temp;
            }
        }

	}

	max_data=d_paixu[Data];


	cx=0;
    for (i = 0; i < Data; i++) {
		
        for (j = cx; j < allData; j++) {
            if(d_allgai[j] <= max_data) {
				for (outi = 0; outi < 40; outi++) {
					dd[i][outi]=d_all[j][outi];
				
				}
				cx=j+1;
				d_allgai[j] = 1;
				break;
            }
        }

	}
	for(i=0;i<Data;i++){
		for(j=0; j<In; j++)
			d_in[i][j]  = dd[i][j];	         //将dd[i][j]中的输入部分先存储于d_in[i][j] 中
		
		k = dd[i][39];
		
		for(j=0; j<Out; j++){
			d_out[i][j] = 0;         //将dd[i][j]中的输出部分先存储于d_out[i][j]中
			dd_out[i][j] = 0; 
		}
		d_out[i][k] = 1; 
		dd_out[i][k] = 1; 

	}

	cx=0;
    for (i = 0; i < testData; i++){
		
        for (j = cx; j < allData; j++)
        {
            if(d_allgai[j] != 1)
            {
				for (outi = 0; outi < 40; outi++){
					test_dd[i][outi]=d_all[j][outi];
				
				}
				cx=j+1;
				break;
            }
        }

	}


}

//将神经元的权值写入各个文件中
void writeNeuron_2() 
{  
	FILE *fp1,*fp2,*fp3; 
	int i,j;

	if((fp1=fopen("neuronw.txt","w"))==NULL){

		printf("can not open the neuron file\n"); 
		exit(0);

	}

	for (i = 0; i < Neuron; ++i){ 
		for (j = 0; j < In; ++j){

			fprintf(fp1,"%lf ",w[i][j]);  //把隐层0的输入权值写入neuronw.txt中

		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	if((fp2=fopen("neuronu.txt","w"))==NULL){

		printf("can not open the neuron file\n"); 
		exit(0);

	}


	for (i = 0; i < Neuron; ++i){

		for (j = 0; j < Neuron1; ++j){

			fprintf(fp2,"%lf ",u[j][i]);  //把隐层1的输入权值写入neuronu.txt中

		}
		fprintf(fp2,"\n");
	}
	fclose(fp2);

	if((fp3=fopen("neuronv.txt","w"))==NULL){

		printf("can not open the neuron file\n"); 
		exit(0);

	}
	for (i = 0; i < Neuron1; ++i){

		for (j = 0; j < Out; ++j){

			fprintf(fp3,"%lf ",v[j][i]);  //把隐层1的输出权值写入neuronv.txt中

		}
		fprintf(fp3,"\n");
	}
	fclose(fp3);

}
