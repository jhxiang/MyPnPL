# MyPnPL
## 安装要求
- g2o https://github.com/RainerKuemmerle/g2o.git git clone后编译
- Eigen库 下载后编译，添加到系统环境变量（linux：/usr/local/include）
- OpenCV 下载后编译，添加到系统环境变量（linux：/usr/local/include）

## 说明
1、这是一个通过图优化(g2o库)求解PnP问题和PnL问题的cpp源文件，可以直接使用，也可以作为头文件使用。使用时需要将源文件中的头文件路径改为自己的路径。

2、cpp文件中flags：
 - 1：仅通过线求解
 - 2：仅通过点求解
 - 3：通过线和点共同求解
