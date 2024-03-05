# MPI安装过程详解（含可能遇到的错误及解决方法）

## 环境：VMWare + Ubuntu16
## 安装步骤

### 第一步：官网下载MPI 安装包（http://www.mpich.org/downloads/ ）

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\1)

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\2)

### 第二步：移动到Ubuntu指定文件夹并解压

- 解压命令：tar -zxvf mpich-3.3.2.tar.gz
- 比如这里解压到/home/michael/Downloads/mpi	

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\3)

### 第三步：进入mpich-3.3.2文件夹

- cd mpich-3.3.2

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\4)

#### 第四步：配置编译环境

- 在此之前，新建一个文件夹作为安装路径，如这里在/home/michael/mpi新建文件夹mpich3，作为安装路径
- 在刚才解压的文件目录下，即：/home/michael/Downloads/mpi/mpich-3.3.2，键入：
  ./configure --prefix=/home/michael/mpi/mpich3
- ![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\5)

- 配置期间可能会出现如下情况

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\6)

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\7)

- 依次键入以下两条命令，即可解决
  sudo apt-get install fort77
  sudo apt-get install gfortran

### 第五步：编译
- make

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\8)

### 第六步：安装

- make install

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\9)

### 第七步：设置环境变量

- vim ~/.bashrc 或者 vi ~/.bashrc
  在最后一行添加export PATH=/home/michael/mpi/mpich3/bin:$PATH
  注意：$PATH一定要添加！！！
- 保存后退出（i 为插入模式，添加后键入ESC 和 :wq，注意冒号）

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\10)

### 第八步：更新环境变量

- source ~/.bashrc
  PS：若在添加路径时，忘记输入$PATH，则更新后ls、vi等指令都不能使用，这时键入
  export PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/root/bin
  即可恢复，然后再vim ~/.bashrc将$PATH输入再更新即可

### 第九步：至此，mpi安装完成

### 第十步：检验效果
- 进入刚解压的路径：/home/michael/Downloads/mpi/mpich-3.3.2/examples

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\11)

- 编译得到目标文件：/home/michael/mpi/mpich3/bin/mpicc hellow.c -o hellow
- 若直接mpicc无法编译（mpirun也一样），则键入mpicc的绝对路径再编译即可（如上面的命令）

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\12)

- 运行hellow文件：/home/michael/mpi/mpich3/bin/mpirun -np 4 ./hellow

![在这里插入图片描述](C:\Users\Dell\Desktop\MarkDown学习\python\csdmtomd\图片\13)

### 第十一步：END



## 错误

### 1.make

# 

## 参考资料：

- https://blog.csdn.net/qq_39709535/article/details/82858793
- https://blog.csdn.net/ibless/article/details/80383760
  ————————————————
  版权声明：本文为CSDN博主「望~」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
  原文链接：https://blog.csdn.net/jiacong_wang/article/details/105593209