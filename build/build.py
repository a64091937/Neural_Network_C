#!/usr/bin/python
# -*- coding: UTF-8 -*- 
#Copyright © 2019-2019, Zhang Hongda, a64091937. All Rights Reserved.
#Description:A Build GCC Python Project
#Author:Zhang Hongda
#Creation time:2019-7-28

import os
import platform
import subprocess

root_Path = os.path.join(os.getcwd(), "..")
Path = [   #所有的需要编译的文件列表
    'src',
    #'src\\func',
]
file_Paths = []               #每个文件的绝对路径
object_file = 'output\\main'  #输出路径
gcc_Path = 'C:\\MinGW\\bin'   #gcc目录
source_C_file = ''            #将所有.c文件的路径写到一起
include_Path = ''             #将所有.h文件的路径写到一起

for i in range(0,len(Path)):
    this_Path = [name for name in os.listdir(root_Path+'\\'+Path[i]) if name.endswith('.c')]  #当前目录里面的文件赋值给列表
    #this_Path = os.listdir(root_Path+'\\'+Path[i])
    for j in range(0,len(this_Path)):
        file_Paths.append(root_Path + '\\' + Path[i] + '\\' + this_Path[j])


for i in range(0,len(file_Paths)):
    source_C_file = source_C_file + ' ' + file_Paths[i]  #Path[0] + '\\' + 'test.c'
#print(source_C_file)

for i in range(0,len(Path)):
    include_Path = include_Path + ' -I ' + root_Path + '\\' +Path[i]
#print(include_Path)

os.chdir(gcc_Path)

build_cmd = gcc_Path + '\\' + "gcc.exe " + source_C_file + include_Path + " -o " + root_Path + '\\' + object_file
#print('build_cmd',build_cmd)
#print("\n")
if 0 == subprocess.call(build_cmd, shell = True):
    print('Build successfully!')

os.chdir(root_Path)
