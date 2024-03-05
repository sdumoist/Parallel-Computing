# 并行计算程序使用说明



通用版本的编译脚本是***pivot***

使用方法：

1.命令行输入 ./pivot 回车

2.出现“请输入线程数：”的提示后，***\*输入线程数\****，回车

 

仅能适用于 k = 2时的版本的编译脚本是***pivotMini***

使用方法：

1.命令行输入 ./pivotMini 回车

2.出现“请输入线程数：”的提示后，***\*输入线程数\****，回车

 

**PS：如果脚本出现问题，可手动输入编译命令：**

gcc -fopenmp -std=c++11 pivot.c -o pivot -O3 -lstdc++ -lm

或

gcc -fopenmp -std=c++11 pivotMini.c -o pivotMini -O3 -lstdc++ -lm



 

 