/*
 Copyright ? 2019-2019, Zhang Hongda, a64091937. All Rights Reserved.
 Description:A Standard Build C Project
 Author:Zhang Hongda
 Creation time:2019-7-28
 */

#ifndef FILE_MANEGER_H
#define FILE_MANEGER_H

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "def.h"
#include "bp.h"

typedef struct {
    INT32U ulDataNum;        //训练数据数
    INT32U ulTestDataNum;    //测试数据数
    INT32U ulAllData;        //总数据数
    INT16U uwDataLen;        //单条数据长度
    INT8U ucInNum;
    INT8U ucOutNum;
    INT8S scDataFile[20];    //储存数据的文件名 
} In_Data_File_Struct;       //输入数据结构体

// typedef struct {
    // FP64* fpTrainData;   //训练数据
    // FP64* fpTestData;    //测试数据
    // FP64* fpData;        //总数据
// } *PIn_Data_Struct;       //处理后的输入数据结构体

void readData();
void writeNeuron_2();
FP64* readData_func(In_Data_File_Struct File_Size, FP64* fpTrainData, FP64* fpTestData, FP64* fpInTrainData,FP64* fpOutTrainData);
void saveData_func(INT32U line,INT32U column, FP64* fpTrainData);

#endif
